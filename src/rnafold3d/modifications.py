"""Modified nucleotides: a catalogue, the sequence notation, their effect on 2D folding and their atoms in 3D.

Notation (one sequence line): a modified nucleotide is written in square brackets, by its MODOMICS short name, a
common alias or its PDB chemical-component code, e.g. ``GGC[m6A]CU[Y]AG[8oxoG]`` or ``[A>I]`` for an A-to-I edited
site. ``*`` between two nucleotides is a phosphorothioate linkage (``A*G``: the phosphate of G carries a sulfur).

How a modification acts (``energy`` of a catalogue entry):

* ``('vienna', name)`` -- measured nearest-neighbour parameters shipped with ViennaRNA >= 2.6 (m6A, pseudouridine,
  inosine, dihydrouridine, purine, 7-deaza-adenosine; Lorenz et al., Varenyk et al. 2023). Applied as soft
  constraints (``fold_compound.sc_mod_*``).
* ``('pairs', {partner: (reference pair, ddG per stack)})`` -- parameters generated here: the Turner 2004 stacking
  energies of a reference canonical pair plus a correction per stack (kcal/mol). This also lets a modification pair
  with a partner its parent cannot pair with (8-oxo-G with A in the syn conformation, lysidine with A). The
  corrections are *estimates* from the literature trend (melting studies), marked as such in the catalogue.
* ``('blocked',)`` -- the Watson-Crick edge carries a substituent (m1A, m1G, m3C, m3U, wybutosine ...): the
  nucleotide is kept unpaired in the 2D structures (hard constraint).
* ``None`` -- no measurable effect on canonical pairing is known; the parent base is folded.

In 3D the engine samples the parent (``build``) bases; afterwards the modification's atoms are placed on the model
(``decorate``): substituents in the base plane or on the sugar, the alternative with the fewest contacts first
(e.g. the N6-methyl of m6A anti when paired, syn when free), the C-glycosidic ring of pseudouridine renamed as in
the PDB chemical component, 8-oxo purines turned syn when they are not paired. Pairs that only a modification allows
(non-canonical for the parent bases) are kept as distance restraints of an A-form helix instead of rigid pairs.
Sugar-pucker preferences are reported but not sampled (the engine samples C3'-endo sugars).
"""
from __future__ import annotations

import json
import math
import re
from functools import lru_cache

import numpy as np

CANON = ('AU', 'UA', 'GC', 'CG', 'GU', 'UG')
PSI_MAP = {'N1': 'C5', 'C2': 'C4', 'O2': 'O4', 'N3': 'N3', 'C4': 'C2', 'O4': 'O2', 'C5': 'N1', 'C6': 'C6'}
BACKBONE = ('P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'")

CATEGORIES = {
    'methyl': 'Methylation',
    'isomer': 'Isomerisation and reduction',
    'editing': 'Editing',
    'damage': 'Damage and lesions',
    'thio': 'Thiolation and halogenation',
    'hyper': 'Hypermodified (tRNA)',
    'sugar': 'Sugar and backbone (2\'-O-methyl, 2\'-F, LNA, DNA, phosphorothioate)',
    'analog': 'Base analogues (synthetic)',
}

# ----------------------------------------------------------------------------------------------- 3D building blocks
# ('ring', name, R, P, Q, bond): substituent on ring atom R, in the plane, on the bisector away from neighbours P, Q
# ('z', name, bonded, angle_ref, dihedral_ref, bond, angle, dihedrals[, flip_when_paired])
#      a stereocentre takes its first substituent with free rotamers and the second relative to the first
#      (dihedral_ref = the sibling, one fixed dihedral near +-120), so that no rotamer choice can invert it
# ('swap', old, new, anchor, bond): atom replaced (O -> S, O2' -> F2'), moved along the anchor bond
# ('del', name) / ('ren', old, new) / ('psi',) C-glycoside permutation / ('bridge', name, A, B, dA, dB) LNA methylene
OME = (('z', 'CM2', "O2'", "C2'", "C1'", 1.42, 116.0, (80.0, 180.0, -60.0)),)
DEOXY = (('del', "O2'"),)
FLUORO = (('swap', "O2'", "F2'", "C2'", 1.39),)
LNA = (('bridge', "C2''", "O2'", "C4'", 1.44, 1.56),)
METHYL5 = (('ring', 'C5M', 'C5', 'C4', 'C6', 1.50),)


def _c5chain(*atoms):
    """C5-CH2-... side chains of the uridine wobble modifications: C5A on the ring, then a z-matrix chain."""
    return (('ring', 'C5A', 'C5', 'C4', 'C6', 1.50),) + atoms


CM5 = _c5chain(('z', 'C5B', 'C5A', 'C5', 'C4', 1.52, 112.0, (90.0, -90.0)),
               ('z', 'O5B', 'C5B', 'C5A', 'C5', 1.25, 118.0, (0.0,)),
               ('z', 'O5C', 'C5B', 'C5A', 'C5', 1.25, 118.0, (180.0,)))
MCM5 = CM5 + (('z', 'C5D', 'O5C', 'C5B', 'C5A', 1.44, 116.0, (180.0,)),)
NCM5 = _c5chain(('z', 'C5B', 'C5A', 'C5', 'C4', 1.52, 112.0, (90.0, -90.0)),
                ('z', 'O5B', 'C5B', 'C5A', 'C5', 1.23, 121.0, (0.0,)),
                ('z', 'N5C', 'C5B', 'C5A', 'C5', 1.33, 116.0, (180.0,)))
MNM5 = _c5chain(('z', 'N5B', 'C5A', 'C5', 'C4', 1.47, 112.0, (90.0, -90.0)),
                ('z', 'C5C', 'N5B', 'C5A', 'C5', 1.47, 112.0, (180.0, 60.0)))
CMNM5 = MNM5 + (('z', 'C5D', 'C5C', 'N5B', 'C5A', 1.52, 112.0, (180.0,)),
                ('z', 'O5E', 'C5D', 'C5C', 'N5B', 1.25, 118.0, (0.0,)),
                ('z', 'O5F', 'C5D', 'C5C', 'N5B', 1.25, 118.0, (180.0,)))
TM5 = MNM5 + (('z', 'C5D', 'C5C', 'N5B', 'C5A', 1.52, 112.0, (180.0,)),
              ('z', 'S5E', 'C5D', 'C5C', 'N5B', 1.78, 112.0, (180.0,)),
              ('z', 'O1S', 'S5E', 'C5D', 'C5C', 1.45, 107.0, (60.0,)),
              ('z', 'O2S', 'S5E', 'C5D', 'C5C', 1.45, 107.0, (180.0,)),
              ('z', 'O3S', 'S5E', 'C5D', 'C5C', 1.45, 107.0, (-60.0,)))
HO5 = (('ring', 'O5A', 'C5', 'C4', 'C6', 1.36),)
MO5 = HO5 + (('z', 'C5B', 'O5A', 'C5', 'C4', 1.43, 117.0, (180.0, 0.0)),)
CMO5 = HO5 + (('z', 'C5B', 'O5A', 'C5', 'C4', 1.43, 117.0, (180.0, 0.0)),
              ('z', 'C5C', 'C5B', 'O5A', 'C5', 1.52, 108.0, (180.0, 90.0)),
              ('z', 'O5D', 'C5C', 'C5B', 'O5A', 1.25, 118.0, (0.0,)),
              ('z', 'O5E', 'C5C', 'C5B', 'O5A', 1.25, 118.0, (180.0,)))
MCMO5 = CMO5 + (('z', 'C5F', 'O5E', 'C5C', 'C5B', 1.44, 116.0, (180.0,)),)
S2U = (('swap', 'O2', 'S2', 'C2', 1.67),)
S4U = (('swap', 'O4', 'S4', 'C4', 1.67),)
MS2 = (('ring', 'S2', 'C2', 'N1', 'N3', 1.75), ('z', 'CS2', 'S2', 'C2', 'N1', 1.80, 101.0, (0.0, 180.0)))
M6A = (('z', 'CM6', 'N6', 'C6', 'N1', 1.46, 123.0, (0.0, 180.0), True),)        # syn when free, anti when paired
I6A = (('z', 'C10', 'N6', 'C6', 'N1', 1.46, 123.0, (0.0, 180.0), True),
       ('z', 'C11', 'C10', 'N6', 'C6', 1.50, 112.0, (90.0, -90.0, 180.0)),
       ('z', 'C12', 'C11', 'C10', 'N6', 1.34, 125.0, (180.0, 100.0, -100.0)),
       ('z', 'C13', 'C12', 'C11', 'C10', 1.50, 123.0, (0.0,)),
       ('z', 'C14', 'C12', 'C11', 'C10', 1.50, 123.0, (180.0,)))
# N6-threonylcarbamoyl: the ureido group lies over the Watson-Crick edge (N11-H...N1), which blocks pairing
T6A = (('z', 'C10', 'N6', 'C6', 'N1', 1.38, 127.0, (0.0, 180.0)),
       ('z', 'O10', 'C10', 'N6', 'C6', 1.23, 121.0, (180.0,)),
       ('z', 'N11', 'C10', 'N6', 'C6', 1.36, 116.0, (0.0,)),
       ('z', 'C12', 'N11', 'C10', 'N6', 1.45, 121.0, (180.0,)),
       ('z', 'C13', 'C12', 'N11', 'C10', 1.53, 111.0, (-60.0, 180.0, 60.0)),
       ('z', 'O13', 'C13', 'C12', 'N11', 1.43, 109.5, (60.0, 180.0, -60.0)),
       ('z', 'C14', 'C13', 'C12', 'O13', 1.52, 109.5, (-120.0,)),            # C-beta of L-threonine (3R)
       ('z', 'C15', 'C12', 'N11', 'C13', 1.53, 109.5, (120.0,)),             # C-alpha of L-threonine (2S)
       ('z', 'O16', 'C15', 'C12', 'N11', 1.25, 118.0, (0.0,)),
       ('z', 'O17', 'C15', 'C12', 'N11', 1.25, 118.0, (180.0,)))
M2G = (('z', 'CN2', 'N2', 'C2', 'N1', 1.46, 121.0, (180.0, 0.0)),)
M22G = (('z', 'CN2', 'N2', 'C2', 'N1', 1.46, 121.0, (180.0,)), ('z', 'CN3', 'N2', 'C2', 'N1', 1.46, 121.0, (0.0,)))
DEAZA7 = (('ren', 'N7', 'C7'),)
Q_SIDE = (('ring', 'C10', 'C7', 'C8', 'C5', 1.50),
          ('z', 'N11', 'C10', 'C7', 'C8', 1.47, 112.0, (90.0, -90.0)),
          ('z', 'C12', 'N11', 'C10', 'C7', 1.47, 112.0, (180.0, 70.0, -70.0)),
          ('z', 'C13', 'C12', 'N11', 'C10', 1.51, 110.4, (90.0, 180.0, -60.0)),
          ('z', 'C16', 'C12', 'N11', 'C13', 1.55, 110.4, (115.3,)),             # C1 of the ring is S
          ('z', 'C14', 'C13', 'C12', 'C16', 1.31, 111.1, (16.4,)),              # C13=C14
          ('z', 'C15', 'C14', 'C13', 'C12', 1.51, 111.1, (0.0,)),
          ('z', 'O15', 'C15', 'C14', 'C13', 1.43, 110.4, (102.4,)),             # cis-diol, trans to the amine
          ('z', 'O16', 'C16', 'C12', 'C13', 1.43, 111.0, (-142.7,)))
PREQ1 = (('ring', 'C10', 'C7', 'C8', 'C5', 1.50), ('z', 'N11', 'C10', 'C7', 'C8', 1.47, 112.0, (90.0, -90.0)))
PREQ0 = (('ring', 'C10', 'C7', 'C8', 'C5', 1.43), ('z', 'N11', 'C10', 'C7', 'C8', 1.15, 178.0, (0.0,)))
ARCH = (('ring', 'C10', 'C7', 'C8', 'C5', 1.47), ('z', 'N11', 'C10', 'C7', 'C8', 1.32, 120.0, (0.0,)),
        ('z', 'N12', 'C10', 'C7', 'C8', 1.32, 120.0, (180.0,)))
ACP3 = (('ring', 'C3A', 'N3', 'C2', 'C4', 1.47),
        ('z', 'C3B', 'C3A', 'N3', 'C2', 1.52, 112.0, (90.0, -90.0)),
        ('z', 'C3C', 'C3B', 'C3A', 'N3', 1.53, 112.0, (180.0, 60.0)),
        ('z', 'N3D', 'C3C', 'C3B', 'C3A', 1.47, 110.0, (60.0,)),
        ('z', 'C3E', 'C3C', 'C3B', 'C3A', 1.53, 110.0, (-60.0,)),
        ('z', 'O3F', 'C3E', 'C3C', 'C3B', 1.25, 118.0, (0.0,)),
        ('z', 'O3G', 'C3E', 'C3C', 'C3B', 1.25, 118.0, (180.0,)))
AC4C = (('z', 'C7', 'N4', 'C4', 'N3', 1.37, 127.0, (180.0, 0.0)),
        ('z', 'O7', 'C7', 'N4', 'C4', 1.22, 122.0, (0.0,)),
        ('z', 'C8', 'C7', 'N4', 'C4', 1.50, 115.0, (180.0,)))
K2C = (('del', 'O2'), ('ring', 'N2', 'C2', 'N1', 'N3', 1.29),
       ('z', 'CEL', 'N2', 'C2', 'N3', 1.46, 120.0, (180.0, 0.0)),
       ('z', 'CDL', 'CEL', 'N2', 'C2', 1.52, 111.0, (180.0, 90.0)),
       ('z', 'CGL', 'CDL', 'CEL', 'N2', 1.52, 112.0, (180.0,)),
       ('z', 'CBL', 'CGL', 'CDL', 'CEL', 1.52, 112.0, (180.0,)),
       ('z', 'CAL', 'CBL', 'CGL', 'CDL', 1.53, 112.0, (180.0,)),
       ('z', 'NAL', 'CAL', 'CBL', 'CGL', 1.47, 110.0, (60.0,)),
       ('z', 'CXL', 'CAL', 'CBL', 'CGL', 1.53, 110.0, (-60.0,)),
       ('z', 'OL1', 'CXL', 'CAL', 'CBL', 1.25, 118.0, (0.0,)),
       ('z', 'OL2', 'CXL', 'CAL', 'CBL', 1.25, 118.0, (180.0,)))
PS = (('swap', 'OP1', 'S1P', 'P', 1.96),)


def _e(code, name, parent, category, energy=None, *, build=None, aliases=(), atoms=(), syn=False,
       pucker=None, evidence='', note='', core_only=False, kind='base'):
    # 'ccd' (PDB chemical component), 'names' and 'also' are filled from COMPONENTS below
    return {'code': code, 'name': name, 'parent': parent, 'build': build or parent, 'category': category,
            'energy': energy, 'ccd': '', 'aliases': tuple(aliases), 'atoms': tuple(atoms), 'syn': syn,
            'pucker': pucker, 'evidence': evidence, 'note': note, 'core_only': core_only, 'kind': kind}


V = lambda name: ('vienna', name)                       # measured parameters shipped with ViennaRNA
BLOCK = ('blocked',)
MEAS = 'measured'
EST = 'estimated'


def _pairs(**kw):
    return ('pairs', {k: v for k, v in kw.items()})


def _nm(parent, d=-0.1):
    """2'-O-methyl, 2'-F ... : the parent's canonical partners with a correction d per stack."""
    return _pairs(**{p: (parent + p, d) for p in 'ACGU' if parent + p in CANON})


ENTRIES = [
    # ---------------------------------------------------------------- adenosine
    _e('m6A', 'N6-methyladenosine', 'A', 'methyl', V('m6A'), aliases=('6mA', 'N6mA'), atoms=M6A,
       evidence=MEAS, note='The most abundant internal mRNA modification. Pairs with U, but the methyl must turn anti '
       'to pair, which costs about 0.5-1.5 kcal/mol per pair; it stacks better when unpaired (Kierzek et al. 2022).'),
    _e('m1A', 'N1-methyladenosine', 'A', 'methyl', BLOCK, atoms=(('ring', 'CM1', 'N1', 'C6', 'C2', 1.47),),
       evidence=MEAS, note='Positively charged methyl on the Watson-Crick edge: blocks A-U pairing (Hoogsteen pairs '
       'remain possible, e.g. m1A58 in tRNA T-loops). Kept unpaired.'),
    _e('m2A', '2-methyladenosine', 'A', 'methyl', _pairs(U=('AU', 0.0)),
       atoms=(('ring', 'C2M', 'C2', 'N1', 'N3', 1.50),), evidence=EST, note='Methyl in the minor groove; pairing with U is '
       'essentially unchanged.'),
    _e('m6,6A', 'N6,N6-dimethyladenosine', 'A', 'methyl', BLOCK, aliases=('m66A', 'm62A'),
       atoms=(('z', 'CM6', 'N6', 'C6', 'N1', 1.46, 123.0, (0.0,)), ('z', 'CM7', 'N6', 'C6', 'N1', 1.46, 123.0, (180.0,))),
       evidence=EST, note='No N6 hydrogen is left to donate, so the Watson-Crick pair with U cannot form (rRNA 3\' end). Kept unpaired.'),
    _e('m8A', '8-methyladenosine', 'A', 'methyl', _pairs(U=('AU', 0.4)), atoms=(('ring', 'C8M', 'C8', 'N7', 'N9', 1.50),),
       syn=True, evidence=EST, note='8-substituted purines prefer the syn conformation; pairing in a helix is weakened.'),
    _e('Am', "2'-O-methyladenosine", 'A', 'sugar', _nm('A'), atoms=OME, pucker="C3'-endo", evidence=EST,
       note="2'-O-methylation stabilises the C3'-endo sugar and helices by about 0.1-0.2 kcal/mol per nucleotide."),
    _e('m6Am', "N6,2'-O-dimethyladenosine", 'A', 'methyl', V('m6A'), atoms=M6A + OME, pucker="C3'-endo", evidence=MEAS,
       note="Cap-adjacent mRNA nucleotide: m6A pairing parameters plus a 2'-O-methyl group."),
    _e('i6A', 'N6-isopentenyladenosine', 'A', 'hyper', _pairs(U=('AU', 0.6)), atoms=I6A, evidence=EST,
       note='tRNA position 37: the hydrophobic isopentenyl group stacks on the anticodon and weakens pairing.'),
    _e('ms2i6A', '2-methylthio-N6-isopentenyladenosine', 'A', 'hyper', _pairs(U=('AU', 0.8)), atoms=I6A + MS2,
       evidence=EST, note='tRNA position 37 (codons starting with U); stacks strongly, pairs weakly.'),
    _e('t6A', 'N6-threonylcarbamoyladenosine', 'A', 'hyper', BLOCK, atoms=T6A, evidence=EST,
       note='tRNA position 37 (ANN codons): the ureido group lies over the Watson-Crick edge (N11-H...N1), which blocks '
       'pairing and pre-orders the anticodon loop. Kept unpaired.'),
    _e('ms2t6A', '2-methylthio-N6-threonylcarbamoyladenosine', 'A', 'hyper', BLOCK, atoms=T6A + MS2, evidence=EST,
       note='As t6A, with a 2-methylthio group. Kept unpaired.'),
    _e('m6t6A', 'N6-methyl-N6-threonylcarbamoyladenosine', 'A', 'hyper', BLOCK, atoms=T6A, evidence=EST,
       note='As t6A; the N6-methyl removes the last N6 hydrogen. Kept unpaired (methyl not placed).'),
    _e('8oxoA', '8-oxo-7,8-dihydroadenosine', 'A', 'damage', _pairs(U=('AU', 0.4)), aliases=('8-oxoA', 'oxoA', 'o8A'),
       atoms=(('ring', 'O8', 'C8', 'N7', 'N9', 1.23),), syn=True, evidence=EST,
       note='Oxidative lesion; prefers syn when unpaired and pairs more weakly with U.'),
    _e('c7A', '7-deazaadenosine (tubercidin)', 'A', 'analog', V('7DA'), aliases=('7DA', '7-deazaA', 'tubercidin'),
       atoms=DEAZA7, evidence=MEAS, note='N7 replaced by CH: no Hoogsteen acceptor; A-U pairing slightly weaker '
       '(Richardson & Znosko 2016).'),
    _e('neb', 'Purine riboside (nebularine)', 'A', 'analog', V('purine'), aliases=('purine', 'P9', 'nebularine'),
       atoms=(('del', 'N6'),), evidence=MEAS, note='Adenine without the 6-amino group: one hydrogen bond less with U '
       '(Jolley & Znosko 2017).'),
    _e('2AP', '2-aminopurine riboside', 'A', 'analog', _pairs(U=('AU', 0.4)), aliases=('2-aminopurine', '2-AP'),
       atoms=(('del', 'N6'), ('ring', 'N2', 'C2', 'N1', 'N3', 1.34)), evidence=EST,
       note='Fluorescent adenine analogue; pairs with U through N1 and the 2-amino group, more weakly than A.'),
    _e('DAP', '2,6-diaminopurine riboside', 'A', 'analog', _pairs(U=('AU', -0.4)), aliases=('2,6-diaminopurine', 'n2A'),
       atoms=(('ring', 'N2', 'C2', 'N1', 'N3', 1.34),), evidence=EST,
       note='Adds a third hydrogen bond to U in the minor groove: DAP-U is more stable than A-U.'),
    # ---------------------------------------------------------------- editing
    _e('I', 'Inosine (A-to-I editing)', 'A', 'editing', V('inosine'), build='G',
       aliases=('A>I', 'A-to-I', 'AtoI', 'Ino', 'inosine', 'ADAR'), atoms=(('del', 'N2'),), evidence=MEAS,
       note='Product of ADAR editing. Inosine pairs with C (Watson-Crick-like, two hydrogen bonds) and with U (wobble, '
       'weaker than A-U), so an edited A-U pair becomes an I-U wobble and an A-C mismatch can become an I-C pair '
       '(Wright et al. 2007, 2018). Built on a guanosine without N2.'),
    _e('m1I', '1-methylinosine', 'A', 'editing', BLOCK, build='G',
       atoms=(('del', 'N2'), ('ring', 'CM1', 'N1', 'C6', 'C2', 1.47)), evidence=EST,
       note='tRNA position 37 (from A-to-I then N1 methylation): the Watson-Crick edge is blocked. Kept unpaired.'),
    _e('Im', "2'-O-methylinosine", 'A', 'editing', V('inosine'), build='G', atoms=(('del', 'N2'),) + OME,
       pucker="C3'-endo", evidence=MEAS, note="Inosine pairing parameters plus a 2'-O-methyl group."),
    _e('C>U', 'C-to-U editing (uridine)', 'C', 'editing', None, build='U', aliases=('C-to-U', 'CtoU', 'APOBEC'),
       evidence=MEAS, note='APOBEC editing turns C into a normal uridine: it then pairs with A (and G as a wobble). '
       'Folded and built as U; the edit is recorded.'),
    # ---------------------------------------------------------------- guanosine
    _e('m1G', '1-methylguanosine', 'G', 'methyl', BLOCK, atoms=(('ring', 'CM1', 'N1', 'C6', 'C2', 1.47),),
       evidence=MEAS, note='Methyl on N1 of the Watson-Crick edge (tRNA position 37, 9): blocks G-C pairing. Kept unpaired.'),
    _e('m2G', 'N2-methylguanosine', 'G', 'methyl', _pairs(C=('GC', 0.0), U=('GU', 0.0)), atoms=M2G,
       evidence=EST, note='The methyl turns towards the minor groove (cis to N3) and leaves G-C pairing unchanged.'),
    _e('m2,2G', 'N2,N2-dimethylguanosine', 'G', 'methyl', _pairs(U=('GU', 0.0), A=('GU', 0.6)),
       aliases=('m22G', 'm2,2G'), atoms=M22G, evidence=EST,
       note='No N2 hydrogen is left and one methyl sits where a C would pair: the Watson-Crick pair with C is lost; the '
       'G-U wobble and the imino G-A pair (tRNA m2,2G26-A44) remain (Pallan et al. 2008). G-A pairs are kept as helix '
       'distance restraints in 3D.'),
    _e('m7G', '7-methylguanosine', 'G', 'methyl', _pairs(C=('GC', 0.3), U=('GU', 0.3)),
       atoms=(('ring', 'CM7', 'N7', 'C8', 'C5', 1.47),), evidence=EST,
       note='Positively charged (internal m7G, e.g. tRNA position 46; also the mRNA cap). The Watson-Crick edge is free; '
       'pairing is slightly weakened.'),
    _e('Gm', "2'-O-methylguanosine", 'G', 'sugar', _nm('G'), atoms=OME, pucker="C3'-endo", evidence=EST,
       note="2'-O-methylation stabilises the C3'-endo sugar and helices by about 0.1-0.2 kcal/mol per nucleotide."),
    _e('m2Gm', "N2,2'-O-dimethylguanosine", 'G', 'methyl', _nm('G'), atoms=M2G + OME, pucker="C3'-endo", evidence=EST,
       note="N2-methyl (minor groove) and 2'-O-methyl."),
    _e('m8G', '8-methylguanosine', 'G', 'methyl', _pairs(C=('GC', 0.6), U=('GU', 0.6)),
       atoms=(('ring', 'C8M', 'C8', 'N7', 'N9', 1.50),), syn=True, evidence=EST,
       note='The 8-methyl group forces the guanosine syn: m8G stabilises left-handed Z-RNA, which r(CGC[m8G]CG) forms '
            'at physiological salt (Balasubramaniyam et al. 2018). Weakens A-form pairing.'),
    _e('8BrG', '8-bromoguanosine', 'G', 'damage', _pairs(C=('GC', 0.6), U=('GU', 0.6)), aliases=('br8G',),
       atoms=(('ring', 'BR8', 'C8', 'N7', 'N9', 1.87),), syn=True, evidence=EST,
       note='Bulky 8-bromo group: syn guanosine, used to lock Z-DNA and Z-RNA; weakens A-form pairing.'),
    _e('8oxoG', '8-oxo-7,8-dihydroguanosine', 'G', 'damage', _pairs(C=('GC', 0.75), A=('GU', 0.8)),
       aliases=('8-oxoG', 'oxoG', 'o8G', '8OG', '8-oxo-G'), atoms=(('ring', 'O8', 'C8', 'N7', 'N9', 1.23),), syn=True,
       evidence=EST, note='The main oxidative RNA lesion. In anti it pairs with C (weaker than G-C); in syn its Hoogsteen '
       'edge pairs with A (oxoG(syn)-A mismatch), which is why it miscodes. Both pairs are allowed in the 2D; oxoG-A pairs '
       'are kept as helix distance restraints in 3D, and an unpaired 8-oxoG is turned syn.'),
    _e('O6mG', 'O6-methylguanosine', 'G', 'damage', _pairs(U=('GU', 0.2), C=('GC', 1.2)), aliases=('O6-meG', 'm6G', 'O6meG'),
       atoms=(('z', 'CM6', 'O6', 'C6', 'N1', 1.44, 117.0, (0.0, 180.0)),), evidence=EST,
       note='Alkylation lesion: the methyl on O6 favours pairing with U (Watson-Crick-like) over C, so it miscodes.'),
    _e('s6G', '6-thioguanosine', 'G', 'thio', _pairs(C=('GC', 0.5), U=('GU', 0.3)), aliases=('6sG', 'thioG'),
       atoms=(('swap', 'O6', 'S6', 'C6', 1.67),), evidence=EST,
       note='Thiopurine drug metabolite and photo-crosslinker; the larger sulfur weakens pairing with C.'),
    _e('c7G', '7-deazaguanosine', 'G', 'analog', _pairs(C=('GC', 0.2), U=('GU', 0.2)), aliases=('7-deazaG',),
       atoms=DEAZA7, evidence=EST, note='N7 replaced by CH: removes the Hoogsteen acceptor and a cation site.'),
    _e('Q', 'Queuosine', 'G', 'hyper', _pairs(C=('GC', 0.0), U=('GU', 0.0)), aliases=('queuosine',),
       atoms=DEAZA7 + Q_SIDE, evidence=EST, note='tRNA wobble position 34 (GUN anticodons): 7-deazaguanosine with a '
       'cyclopentenediol side chain in the major groove; pairs with C and U.'),
    _e('preQ1', '7-aminomethyl-7-deazaguanosine', 'G', 'hyper', _pairs(C=('GC', 0.0), U=('GU', 0.0)), atoms=DEAZA7 + PREQ1,
       evidence=EST, note='Queuosine precursor.'),
    _e('preQ0', '7-cyano-7-deazaguanosine', 'G', 'hyper', _pairs(C=('GC', 0.0), U=('GU', 0.0)), atoms=DEAZA7 + PREQ0,
       evidence=EST, note='Queuosine and archaeosine precursor.'),
    _e('G+', 'Archaeosine', 'G', 'hyper', _pairs(C=('GC', 0.0), U=('GU', 0.0)), aliases=('archaeosine',),
       atoms=DEAZA7 + ARCH, evidence=EST, note='Archaeal tRNA position 15: a positively charged formamidino group in '
       'the major groove stabilises the tertiary fold.'),
    _e('yW', 'Wybutosine', 'G', 'hyper', BLOCK, aliases=('wybutosine', 'Y-base'), core_only=True, evidence=MEAS,
       note='tRNA-Phe position 37: tricyclic base whose extra ring covers the Watson-Crick edge; stacks strongly. Kept '
       'unpaired. Only the guanine core is built.'),
    _e('imG', 'Wyosine', 'G', 'hyper', BLOCK, aliases=('wyosine',), core_only=True, evidence=MEAS,
       note='Tricyclic G of archaeal tRNA position 37; Watson-Crick edge blocked. Kept unpaired; guanine core only.'),
    # ---------------------------------------------------------------- cytidine
    _e('m5C', '5-methylcytidine', 'C', 'methyl', _pairs(G=('CG', -0.1)), atoms=METHYL5, evidence=EST,
       note='Major-groove methyl: slightly stabilises C-G pairs and stacking (Hopfinger et al. 2020).'),
    _e('m3C', '3-methylcytidine', 'C', 'methyl', BLOCK, atoms=(('ring', 'CM3', 'N3', 'C2', 'C4', 1.47),), evidence=MEAS,
       note='Methyl on N3 of the Watson-Crick edge (tRNA position 32): blocks pairing. Kept unpaired.'),
    _e('m4C', 'N4-methylcytidine', 'C', 'methyl', _pairs(G=('CG', 0.2)),
       atoms=(('z', 'CN4', 'N4', 'C4', 'N3', 1.46, 122.0, (0.0, 180.0), True),), evidence=EST,
       note='The methyl turns away from the Watson-Crick edge to pair with G, at a small cost.'),
    _e('m4,4C', 'N4,N4-dimethylcytidine', 'C', 'methyl', BLOCK, aliases=('m44C',),
       atoms=(('z', 'CN4', 'N4', 'C4', 'N3', 1.46, 122.0, (0.0,)), ('z', 'CN5', 'N4', 'C4', 'N3', 1.46, 122.0, (180.0,))),
       evidence=EST, note='No N4 hydrogen left and a methyl on the pairing edge: no Watson-Crick pair with G. Kept unpaired.'),
    _e('ac4C', 'N4-acetylcytidine', 'C', 'methyl', _pairs(G=('CG', -0.2)), aliases=('ac4C',), atoms=AC4C,
       pucker="C3'-endo", evidence=EST, note='The acetyl lies in the major groove next to C5-H and stabilises C-G pairs '
       'and the C3\'-endo sugar (Bartee et al. 2022).'),
    _e('hm5C', '5-hydroxymethylcytidine', 'C', 'methyl', _pairs(G=('CG', 0.0)),
       atoms=(('ring', 'C5A', 'C5', 'C4', 'C6', 1.50), ('z', 'O5B', 'C5A', 'C5', 'C4', 1.43, 111.0, (90.0, -90.0))),
       evidence=EST, note='Oxidation product of m5C; pairing with G unchanged.'),
    _e('f5C', '5-formylcytidine', 'C', 'methyl', _pairs(G=('CG', 0.1)),
       atoms=(('ring', 'C5A', 'C5', 'C4', 'C6', 1.46), ('z', 'O5B', 'C5A', 'C5', 'C4', 1.22, 124.0, (0.0, 180.0))),
       evidence=EST, note='Mitochondrial tRNA-Met wobble base 34: lets the anticodon read AUA as well as AUG.'),
    _e('ho5C', '5-hydroxycytidine', 'C', 'methyl', _pairs(G=('CG', 0.0)), atoms=HO5, evidence=EST,
       note='Oxidation product (bacterial rRNA).'),
    _e('Cm', "2'-O-methylcytidine", 'C', 'sugar', _nm('C'), atoms=OME, pucker="C3'-endo", evidence=EST,
       note="2'-O-methylation stabilises the C3'-endo sugar and helices by about 0.1-0.2 kcal/mol per nucleotide."),
    _e('m5Cm', "5,2'-O-dimethylcytidine", 'C', 'methyl', _nm('C', -0.15), atoms=METHYL5 + OME, pucker="C3'-endo",
       evidence=EST, note="5-methyl and 2'-O-methyl."),
    _e('ac4Cm', "N4-acetyl-2'-O-methylcytidine", 'C', 'methyl', _nm('C', -0.25), atoms=AC4C + OME, pucker="C3'-endo",
       evidence=EST, note="Thermophile tRNA: acetyl and 2'-O-methyl both stiffen the nucleotide."),
    _e('s2C', '2-thiocytidine', 'C', 'thio', _pairs(G=('CG', 0.3)), atoms=S2U, evidence=EST,
       note='tRNA position 32; the 2-thio group weakens pairing with G slightly and restricts the sugar.'),
    _e('k2C', 'Lysidine', 'C', 'hyper', _pairs(A=('CG', 0.6)), aliases=('lysidine', 'L'), atoms=K2C, evidence=EST,
       note='Bacterial tRNA-Ile wobble base 34: the lysine on C2 switches pairing from G to A (reads AUA). The C-A pair '
       'is allowed in the 2D and kept as a helix distance restraint in 3D.'),
    _e('5BrC', '5-bromocytidine', 'C', 'thio', _pairs(G=('CG', -0.1)), aliases=('br5C',),
       atoms=(('ring', 'BR5', 'C5', 'C4', 'C6', 1.89),), evidence=EST, note='Heavy-atom and crosslinking analogue.'),
    # ---------------------------------------------------------------- uridine
    _e('Y', 'Pseudouridine (Ψ)', 'U', 'isomer', V('pseudouridine'), aliases=('Ψ', 'psi', 'pseudoU', 'pU', 'PSU'),
       atoms=(('psi',),), pucker="C3'-endo", evidence=MEAS, note='C5-glycoside isomer of U with an extra N1-H (water-mediated contact to '
       'the backbone): stabilises Ψ-A pairs and stacking (Hudson et al. 2013) and the C3\'-endo sugar.'),
    _e('m1Y', 'N1-methylpseudouridine (m1Ψ)', 'U', 'isomer', V('pseudouridine'), aliases=('m1Ψ', 'm1psi', 'N1mY'),
       atoms=(('psi',), ('ring', 'CM1', 'N1', 'C2', 'C6', 1.47)), pucker="C3'-endo", evidence=EST,
       note='Used in place of every U in mRNA vaccines. Pairs with A like Ψ (approximated with the Ψ parameters).'),
    _e('m3Y', 'N3-methylpseudouridine', 'U', 'isomer', BLOCK, aliases=('m3Ψ',),
       atoms=(('psi',), ('ring', 'CM3', 'N3', 'C2', 'C4', 1.47)), evidence=EST,
       note='Methyl on N3 of the Watson-Crick edge (rRNA): blocks pairing. Kept unpaired.'),
    _e('Ym', "2'-O-methylpseudouridine", 'U', 'isomer', V('pseudouridine'), aliases=('Ψm',), atoms=(('psi',),) + OME,
       pucker="C3'-endo", evidence=MEAS, note="Ψ pairing parameters plus a 2'-O-methyl group."),
    _e('D', 'Dihydrouridine', 'U', 'isomer', V('dihydrouridine'), aliases=('DHU', 'dihydroU'), pucker="C2'-endo",
       evidence=MEAS, note='Saturated C5-C6 bond: a non-planar base that stacks poorly and favours the flexible C2\'-endo '
       'sugar (D-loop of tRNA); destabilises helices (Varenyk et al. 2023).'),
    _e('m3U', '3-methyluridine', 'U', 'methyl', BLOCK, atoms=(('ring', 'CM3', 'N3', 'C2', 'C4', 1.47),),
       evidence=MEAS, note='Methyl on N3 of the Watson-Crick edge: blocks pairing. Kept unpaired.'),
    _e('m5U', '5-methyluridine (ribothymidine, T)', 'U', 'methyl', _pairs(A=('UA', -0.1), G=('UG', 0.0)),
       aliases=('rT', 'T', 'ribothymidine'), atoms=METHYL5, evidence=EST,
       note='tRNA position 54 (T-loop): the methyl stacks on its neighbour and slightly stabilises U-A pairs.'),
    _e('Um', "2'-O-methyluridine", 'U', 'sugar', _nm('U'), atoms=OME, pucker="C3'-endo", evidence=EST,
       note="2'-O-methylation stabilises the C3'-endo sugar and helices by about 0.1-0.2 kcal/mol per nucleotide."),
    _e('m5Um', "5,2'-O-dimethyluridine", 'U', 'methyl', _nm('U', -0.15), atoms=METHYL5 + OME, pucker="C3'-endo",
       evidence=EST, note="5-methyl and 2'-O-methyl."),
    _e('s2U', '2-thiouridine', 'U', 'thio', _pairs(A=('UA', -0.4), G=('UG', 0.6)), aliases=('2sU',), atoms=S2U,
       pucker="C3'-endo", evidence=EST, note='Locks the C3\'-endo sugar: stabilises pairing with A and discriminates '
       'against the G wobble (tRNA wobble position 34; Kumar & Davis 1997).'),
    _e('s4U', '4-thiouridine', 'U', 'thio', _pairs(A=('UA', 0.2), G=('UG', 0.0)), aliases=('4sU',), atoms=S4U,
       evidence=EST, note='tRNA position 8 and metabolic labelling (photo-crosslinker); pairs slightly more weakly with A.'),
    _e('s2Um', "2-thio-2'-O-methyluridine", 'U', 'thio', _pairs(A=('UA', -0.5), G=('UG', 0.6)), atoms=S2U + OME,
       pucker="C3'-endo", evidence=EST, note='Thermophile tRNA: a very rigid C3\'-endo nucleotide.'),
    _e('m5s2U', '5-methyl-2-thiouridine', 'U', 'thio', _pairs(A=('UA', -0.5), G=('UG', 0.6)), aliases=('s2T',),
       atoms=METHYL5 + S2U, pucker="C3'-endo", evidence=EST, note='Thermophile tRNA position 54.'),
    _e('ho5U', '5-hydroxyuridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), atoms=HO5, evidence=EST,
       note='Precursor of mo5U and cmo5U.'),
    _e('mo5U', '5-methoxyuridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), atoms=MO5, evidence=EST,
       note='Gram-positive tRNA wobble base 34 (extended wobble); also used in therapeutic mRNA.'),
    _e('cmo5U', 'Uridine 5-oxyacetic acid', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), aliases=('V',), atoms=CMO5,
       evidence=EST, note='Gram-negative tRNA wobble base 34: reads A, G and U.'),
    _e('mcmo5U', 'Uridine 5-oxyacetic acid methyl ester', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), atoms=MCMO5,
       evidence=EST, note='tRNA wobble base 34.'),
    _e('cm5U', '5-carboxymethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.1)), atoms=CM5, evidence=EST,
       note='tRNA wobble base 34.'),
    _e('mcm5U', '5-methoxycarbonylmethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.1)), atoms=MCM5,
       evidence=EST, note='Eukaryotic tRNA wobble base 34.'),
    _e('mcm5s2U', '5-methoxycarbonylmethyl-2-thiouridine', 'U', 'hyper', _pairs(A=('UA', -0.4), G=('UG', 0.6)),
       atoms=MCM5 + S2U, pucker="C3'-endo", evidence=EST, note='Eukaryotic tRNA-Lys/Glu/Gln wobble base 34: reads A '
       'strongly, G poorly.'),
    _e('ncm5U', '5-carbamoylmethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.1)), atoms=NCM5, evidence=EST,
       note='Eukaryotic tRNA wobble base 34.'),
    _e('mnm5U', '5-methylaminomethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), atoms=MNM5, evidence=EST,
       note='Bacterial tRNA wobble base 34.'),
    _e('mnm5s2U', '5-methylaminomethyl-2-thiouridine', 'U', 'hyper', _pairs(A=('UA', -0.4), G=('UG', 0.3)),
       atoms=MNM5 + S2U, pucker="C3'-endo", evidence=EST, note='Bacterial tRNA-Lys/Glu wobble base 34.'),
    _e('cmnm5U', '5-carboxymethylaminomethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.0)), atoms=CMNM5,
       evidence=EST, note='Mitochondrial and bacterial tRNA wobble base 34.'),
    _e('cmnm5s2U', '5-carboxymethylaminomethyl-2-thiouridine', 'U', 'hyper', _pairs(A=('UA', -0.4), G=('UG', 0.3)),
       atoms=CMNM5 + S2U, pucker="C3'-endo", evidence=EST, note='Bacterial tRNA wobble base 34.'),
    _e('tm5U', '5-taurinomethyluridine', 'U', 'hyper', _pairs(A=('UA', 0.0), G=('UG', 0.3)), atoms=TM5, evidence=EST,
       note='Human mitochondrial tRNA-Leu(UUR) wobble base 34 (lost in MELAS).'),
    _e('tm5s2U', '5-taurinomethyl-2-thiouridine', 'U', 'hyper', _pairs(A=('UA', -0.4), G=('UG', 0.6)), atoms=TM5 + S2U,
       pucker="C3'-endo", evidence=EST, note='Human mitochondrial tRNA-Lys wobble base 34 (lost in MERRF).'),
    _e('acp3U', '3-(3-amino-3-carboxypropyl)uridine', 'U', 'hyper', BLOCK, atoms=ACP3, evidence=EST,
       note='tRNA position 47 (variable loop): the substituent on N3 blocks pairing. Kept unpaired.'),
    _e('5BrU', '5-bromouridine', 'U', 'thio', _pairs(A=('UA', -0.1), G=('UG', 0.0)), aliases=('br5U', 'BrU'),
       atoms=(('ring', 'BR5', 'C5', 'C4', 'C6', 1.89),), evidence=EST, note='Heavy-atom, photo-crosslinking analogue.'),
    _e('5IU', '5-iodouridine', 'U', 'thio', _pairs(A=('UA', -0.1), G=('UG', 0.0)), aliases=('i5U', 'IU'),
       atoms=(('ring', 'I5', 'C5', 'C4', 'C6', 2.09),), evidence=EST, note='Heavy-atom, photo-crosslinking analogue.'),
    _e('5FU', '5-fluorouridine', 'U', 'thio', _pairs(A=('UA', 0.0), G=('UG', -0.2)), aliases=('f5U', 'FU'),
       atoms=(('ring', 'F5', 'C5', 'C4', 'C6', 1.35),), evidence=EST,
       note='Anticancer drug metabolite incorporated into RNA; favours the G wobble a little.'),
    _e('5EU', '5-ethynyluridine', 'U', 'analog', _pairs(A=('UA', -0.1), G=('UG', 0.0)), aliases=('EU',),
       atoms=(('ring', 'C5A', 'C5', 'C4', 'C6', 1.43), ('z', 'C5B', 'C5A', 'C5', 'C4', 1.20, 178.0, (0.0,))), evidence=EST,
       note='Metabolic labelling of new RNA (click chemistry).'),
    # ---------------------------------------------------------------- sugar and backbone (any base)
    *[_e(f'f{b}', f"2'-fluoro-2'-deoxy{nm}", b, 'sugar', _nm(b, -0.15), aliases=(f"2'F{b}", f'2F{b}'), atoms=FLUORO,
         pucker="C3'-endo", evidence=EST, note="siRNA/aptamer chemistry: 2'-F favours C3'-endo and raises helix stability "
         "(about 0.2-0.5 kcal/mol per nucleotide).") for b, nm in (('A', 'adenosine'), ('C', 'cytidine'), ('G', 'guanosine'), ('U', 'uridine'))],
    *[_e(f'+{b}', f'LNA {nm}', 'U' if b == 'T' else b, 'sugar', _nm('U' if b == 'T' else b, -0.4), aliases=(f'ln{b}', f'LNA-{b}', f'L{b}'),
         atoms=LNA + (METHYL5 if b in 'TC' else ()), pucker="C3'-endo (locked)", evidence=EST,
         note="Locked nucleic acid: a 2'-O,4'-C methylene bridge locks the C3'-endo sugar and stabilises helices strongly "
         "(about 0.5-1 kcal/mol per nucleotide)." + (' LNA-C and LNA-T carry a 5-methyl.' if b in 'TC' else ''))
      for b, nm in (('A', 'adenosine'), ('C', '5-methylcytidine'), ('G', 'guanosine'), ('T', 'thymidine'))],
    *[_e(f'd{b}', f'2\'-deoxy{nm} (DNA)', 'U' if b == 'T' else b, 'sugar', _nm('U' if b == 'T' else b, 0.2), aliases=(f'D{b}',),
         atoms=DEOXY + (METHYL5 if b == 'T' else ()), pucker="C2'-endo", evidence=EST,
         note="A DNA nucleotide in RNA (chimera, hybrid): no 2'-OH, prefers C2'-endo and destabilises an RNA helix a little.")
      for b, nm in (('A', 'adenosine'), ('C', 'cytidine'), ('G', 'guanosine'), ('T', 'thymidine'), ('U', 'uridine'))],
    _e('ps', 'Phosphorothioate linkage', 'N', 'sugar', None, aliases=('*', 'PS'), atoms=PS, kind='backbone', evidence=EST,
       note='Written as * between two nucleotides (A*G): one non-bridging oxygen of the phosphate of the second '
       'nucleotide is replaced by sulfur (nuclease-resistant antisense and siRNA backbones). No effect on the 2D.'),
]
BY_CODE = {e['code']: e for e in ENTRIES}

# ----------------------------------------------------------------------------- PDB chemical components
# code: (component written to PDB/mmCIF files, {RCV atom name: atom name in that component}
#        [, {other component read as this modification: its atom names}])
# Every row is checked against the wwPDB chemical component dictionary (tests/data/pdb_components.json, read by
# tests/test_components.py): the residue RCV builds has the atoms, the bonds and the stereocentres of the component,
# and the files carry the dictionary's atom names. Entries that are not listed have no nucleotide component of
# exactly that chemistry (m6Am, 8-oxo-G, the xm5U family without sulfur ...) and are written under the parent
# base's name with RCV's atom names. Components read but not written: 7MG and QUO are reduced forms (8-hydro-m7G,
# a saturated cyclopentyl ring) that older tRNA entries use for m7G and queuosine, B8Q a reduced m3C, YG the second
# wybutosine component, XSX the (3R) epimer of acp3U, A1I9V a five-character code that a PDB-format file cannot
# hold. DNA nucleotides (dA ...) keep the parent name on purpose: a DA/DT entry here would make every DNA chain of
# an opened file read as RNA.
COMPONENTS = {
    'm6A': ('6MZ', {'CM6': 'C9'}),
    'm1A': ('1MA', {}),
    'm2A': ('2MA', {'C2M': 'CM2'}),
    'm6,6A': ('MA6', {'CM6': 'C9', 'CM7': 'C10'}),
    'Am': ('A2M', {'CM2': "CM'"}),
    'i6A': ('6IA', {'C10': 'C12', 'C11': 'C13', 'C12': 'C14', 'C13': 'C15', 'C14': 'C16'}),
    'ms2i6A': ('MIA', {'C10': 'C12', 'C11': 'C13', 'C12': 'C14', 'C13': 'C15', 'C14': 'C16', 'S2': 'S10',
        'CS2': 'C11'}),
    't6A': ('T6A', {'C13': 'C14', 'O13': 'O14', 'C14': 'C15', 'C15': 'C13', 'O16': 'ODA', 'O17': 'ODB'}),
    'ms2t6A': ('12A', {'C10': 'CC', 'O10': 'OO', 'N11': 'N', 'C12': 'CA', 'C13': 'CB', 'O13': 'OG1', 'C14': 'CG2',
        'C15': 'C', 'O16': 'O', 'O17': 'OXT', 'CS2': 'C2M'}),
    'neb': ('P5P', {}),
    'DAP': ('N6G', {}),
    'I': ('I', {}),
    'm1I': ('7S3', {}),
    'm1G': ('1MG', {}),
    'm2G': ('2MG', {'CN2': 'CM2'}),
    'm2,2G': ('M2G', {'CN2': 'CM1', 'CN3': 'CM2'}),
    'm7G': ('G7M', {'CM7': 'CN7'}, {'7MG': {}}),
    'Gm': ('OMG', {}),
    '8BrG': ('GRB', {'BR8': 'BR'}),
    'O6mG': ('B8W', {'CM6': 'C61'}),
    'Q': ('56B', {'C10': 'C9', 'N11': 'N10', 'C12': 'C10', 'C13': 'C14', 'C16': 'C11', 'C14': 'C13', 'C15': 'C12',
        'O15': 'O12', 'O16': 'O11'}, {'QUO': {'C13': 'C16', 'C16': 'C13', 'C14': 'C15', 'C15': 'C14',
        'O15': 'O14', 'O16': 'O13'}}),
    'preQ1': ('PQ1', {}),
    'yW': ('YYG', {}, {'YG': {}}),
    'm5C': ('5MC', {'C5M': 'CM5'}),
    'm3C': ('JMH', {'CM3': 'C31'}, {'B8Q': {'CM3': 'C31'}}),
    'm4C': ('B8T', {'CN4': 'C41'}),
    'm4,4C': ('LV2', {'CN4': 'CM4', 'CN5': 'C1'}),
    'ac4C': ('4AC', {'C8': 'CM7'}),
    'hm5C': ('5HM', {'C5A': 'CM5', 'O5B': 'OM5'}),
    'f5C': ('RSQ', {'C5A': 'C10', 'O5B': 'O30'}),
    'ho5C': ('TKW', {'O5A': 'OAC'}),
    'Cm': ('OMC', {}),
    'ac4Cm': ('LHH', {'C8': 'CM7', 'CM2': 'C1'}),
    's2C': ('RSP', {}),
    'k2C': ('M3X', {'CEL': 'CE', 'CDL': 'CD', 'CGL': 'CG', 'CBL': 'CB', 'CAL': 'CZ', 'NAL': 'N', 'CXL': 'C',
        'OL1': 'O', 'OL2': 'OXT'}),
    'Y': ('PSU', {}),
    'm1Y': ('B8H', {'CM1': 'CN1'}),
    'm3Y': ('3TD', {'CM3': 'C10'}),
    'Ym': ('UY1', {}),
    'D': ('H2U', {}),
    'm3U': ('UR3', {'CM3': 'C3U'}),
    'm5U': ('5MU', {}),
    'Um': ('OMU', {}),
    'm5Um': ('2MU', {'CM2': "C6'"}),
    's2U': ('SUR', {}),
    's4U': ('4SU', {}),
    'm5s2U': ('', {}, {'A1I9V': {}}),
    'cmo5U': ('CM0', {'O5A': 'O5', 'C5B': 'C7', 'C5C': 'C8', 'O5D': 'O8', 'O5E': 'O9'}),
    'mcm5s2U': ('70U', {'C5A': 'C5M', 'C5B': 'C8', 'O5B': 'O8', 'O5C': 'O9', 'C5D': 'C9'}),
    'mnm5U': ('MNU', {'C5A': 'C7', 'N5B': 'N8', 'C5C': 'C9'}),
    'mnm5s2U': ('U8U', {'C5A': 'C', 'N5B': 'N', 'C5C': 'CA'}),
    'acp3U': ('', {}, {'XSX': {'C3A': 'C1', 'C3B': 'C3', 'C3C': 'C7', 'N3D': 'N4', 'C3E': 'C9', 'O3F': 'O10',
        'O3G': 'O11'}}),
    '5BrU': ('5BU', {'BR5': 'BR'}),
    '5IU': ('IU', {}),
    '5FU': ('5FU', {}),
    'fA': ('AF2', {"F2'": 'F'}),
    'fC': ('CFZ', {}),
    'fG': ('GF2', {"F2'": 'F'}),
    'fU': ('UFT', {}),
    '+A': ('LCA', {"C2''": "C6'"}),
    '+C': ('LCC', {"C2''": "C6'"}),
    '+G': ('LCG', {"C2''": "C6'"}),
    '+T': ('TLN', {"C2''": "C6'"}),
}
for _ent in ENTRIES:
    _ccd, _names, *_also = COMPONENTS.get(_ent['code'], ('', {}))
    _ent['ccd'], _ent['names'], _ent['also'] = _ccd, dict(_names), dict(_also[0]) if _also else {}
_INDEX = {}
for _ent in ENTRIES:
    for _k in (_ent['code'], *_ent['aliases'], *([_ent['ccd']] if _ent['ccd'] else []), *_ent['also']):
        _INDEX.setdefault(_k, _ent['code'])
_LOWER = {}
for _k, _c in _INDEX.items():
    _LOWER.setdefault(_k.lower(), set()).add(_c)


def lookup(token):
    """Catalogue entry for a code, alias or PDB component id (exact first, then case-insensitive if unique)."""
    t = token.strip().replace('Ψ', 'Ψ').replace('ψ', 'Ψ')
    if t in _INDEX:
        return BY_CODE[_INDEX[t]]
    hits = _LOWER.get(t.lower(), set())
    if len(hits) == 1:
        return BY_CODE[next(iter(hits))]
    t2 = t.replace(' ', '').replace('_', '')
    if t2 != t:
        return lookup(t2)
    return None


def compatible(entry, base):
    """Can this catalogue entry be written at a nucleotide whose parent base is `base`?"""
    return entry['kind'] == 'backbone' or entry['parent'] == base


# ------------------------------------------------------------------------------------------------- notation
TOKEN = re.compile(r'\[([^\[\]]{1,32})\]')


def parse(text):
    """Sequence line with modifications -> (sequence of parent/build bases, modifications, problems).

    modifications: [{'pos': 1-based, 'code'}] (one base/sugar modification per position, plus 'ps' for a
    phosphorothioate). problems: messages (unknown names, wrong letters)."""
    s = re.sub(r'\s+', '', str(text or ''))
    seq, mods, problems = [], [], []
    pending_ps = False
    k = 0
    while k < len(s):
        ch = s[k]
        if ch == '[':
            m = TOKEN.match(s, k)
            if not m:
                problems.append(f'unclosed "[" at character {k + 1}')
                break
            tok = m.group(1)
            ent = lookup(tok)
            if ent is None:
                problems.append(f'unknown modification [{tok}]')
                seq.append('N')
            elif ent['kind'] == 'backbone':
                problems.append(f'[{tok}] is a linkage: write * between two nucleotides')
            else:
                seq.append(ent['build'])
                mods.append({'pos': len(seq), 'code': ent['code']})
                if pending_ps:
                    mods.append({'pos': len(seq), 'code': 'ps'})
            pending_ps = False
            k = m.end()
            continue
        if ch == '*':
            if not seq:
                problems.append('"*" (phosphorothioate) must stand between two nucleotides')
            pending_ps = True
            k += 1
            continue
        if ch.isdigit():               # position numbers of a copied sequence ("1 ggcagatctg 11 agcc...")
            k += 1
            continue
        up = ch.upper()
        if up == 'T':
            up = 'U'
        if up in 'ACGU':
            seq.append(up)
            if pending_ps:
                mods.append({'pos': len(seq), 'code': 'ps'})
        else:
            seq.append('N')
            problems.append(f'letter {ch}')
        pending_ps = False
        k += 1
    if pending_ps:
        problems.append('"*" at the end of the sequence')
    return ''.join(seq), mods, problems


def has_notation(text):
    return '[' in str(text or '') or '*' in str(text or '')


def format_sequence(seq, mods):
    """Inverse of parse: sequence with [code] tokens and * linkages."""
    by = {m['pos']: m['code'] for m in mods if m['code'] != 'ps'}
    ps = {m['pos'] for m in mods if m['code'] == 'ps'}
    out = []
    for i, b in enumerate(seq, 1):
        if i in ps and i > 1:
            out.append('*')
        out.append(f'[{by[i]}]' if i in by else b)
    return ''.join(out)


def normalize(mods, seq):
    """Validate a modification list against the sequence: [{'pos', 'code'}] sorted, codes canonical."""
    out, seen = [], set()
    for m in mods or ():
        ent = lookup(str(m.get('code', '')))
        pos = int(m.get('pos', 0))
        if ent is None:
            raise ValueError(f'unknown modification {m.get("code")!r}')
        if not 1 <= pos <= len(seq):
            raise ValueError(f'modification {ent["code"]} at position {pos} is outside the sequence (1-{len(seq)})')
        if ent['kind'] != 'backbone' and seq[pos - 1] != ent['build']:
            raise ValueError(f'{ent["code"]} at position {pos} needs {ent["build"]} in the sequence, found {seq[pos - 1]}')
        key = (pos, ent['kind'])
        if key in seen:
            raise ValueError(f'two modifications at position {pos}')
        seen.add(key)
        out.append({'pos': pos, 'code': ent['code']})
    return sorted(out, key=lambda m: (m['pos'], m['code'] == 'ps'))


def summary(mods):
    """Short text: 'm6A ×1, Ψ ×2'."""
    count = {}
    for m in mods:
        count[m['code']] = count.get(m['code'], 0) + 1
    return ', '.join(f'{c} ×{n}' for c, n in count.items())


# ------------------------------------------------------------------------------------------------- 2D (ViennaRNA)
VIENNA_UNMOD = {'m6A': 'A', 'pseudouridine': 'U', 'inosine': 'A', 'dihydrouridine': 'U', 'purine': 'A', '7DA': 'A'}
_PT = {'CG': 1, 'GC': 2, 'GU': 3, 'UG': 4, 'AU': 5, 'UA': 6}
_LETTERS = '0123456789'


@lru_cache(maxsize=1)
def _stack_table():
    import RNA
    p = RNA.param(RNA.md())
    return [list(row) for row in p.stack], p.TerminalAU / 100.0


def _stack(outer, inner):
    """Turner stacking energy (kcal/mol) of 5'-outer[0] inner[0]-3' / 3'-outer[1] inner[1]-5'."""
    table, _ = _stack_table()
    return table[_PT[outer]][_PT[inner[1] + inner[0]]] / 100.0


def generated_json(entry, code_letter, unmodified):
    """ViennaRNA modified-base JSON for an entry with ('pairs', {partner: (ref pair, ddG per stack)})."""
    spec = entry['energy'][1]
    _, term_au = _stack_table()
    M = code_letter
    stacks, terms = {}, {}
    for partner, (ref, ddg) in spec.items():
        r, rp = ref[0], ref[1]
        for yz in CANON:
            y, z = yz
            stacks[f'{M}{y}{partner}{z}'] = round(_stack(r + rp, y + z) + ddg, 2)     # M at the 5' end of the top strand
            stacks[f'{y}{M}{z}{partner}'] = round(_stack(y + z, r + rp) + ddg, 2)     # M at the 3' end of the top strand
            stacks[f'{partner}{z}{M}{y}'] = round(_stack(rp + r, z + y) + ddg, 2)     # M in the bottom strand
            stacks[f'{z}{partner}{y}{M}'] = round(_stack(z + y, rp + r) + ddg, 2)
        t = 0.0 if ref in ('GC', 'CG') else term_au
        terms[f'{M}{partner}'] = t
        terms[f'{partner}{M}'] = t
    return json.dumps({'modified_base': {
        'name': entry['name'], 'unmodified': unmodified, 'one_letter_code': M, 'fallback': entry['build'],
        'pairing_partners': sorted(spec), 'stacking_energies': stacks, 'terminal_energies': terms,
        'sources': [{'authors': 'RCV', 'title': 'estimated from Turner 2004 stacks of a reference pair plus a correction',
                     'journal': '-', 'year': 2026, 'doi': '-'}]}})


def vienna_sequence(seq, mods):
    """The sequence ViennaRNA folds: a site with built-in parameters carries their 'unmodified' letter."""
    s = list(seq)
    for m in mods:
        e = BY_CODE[m['code']]
        if e['energy'] and e['energy'][0] == 'vienna':
            s[m['pos'] - 1] = VIENNA_UNMOD[e['energy'][1]]
    return ''.join(s)


def apply_to_fold_compound(fc, mods):
    """Soft and hard constraints of the modifications on a ViennaRNA fold compound (sequence from vienna_sequence).
    Returns {'measured': [...], 'estimated': [...], 'blocked': [...], 'none': [...]} (1-based positions)."""
    import RNA
    used = {'measured': [], 'estimated': [], 'blocked': [], 'none': []}
    groups = {}
    for m in mods:
        e = BY_CODE[m['code']]
        en = e['energy']
        if e['kind'] == 'backbone' and en is None:
            continue
        if en is None:
            used['none'].append(m['pos'])
        elif en[0] == 'blocked':
            fc.hc_add_up(m['pos'], RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
            used['blocked'].append(m['pos'])
        else:
            groups.setdefault((en[0], en[1] if en[0] == 'vienna' else e['code']), []).append(m['pos'])
    gen_k = 0
    for (kind, key), sites in groups.items():
        if kind == 'vienna':
            fn = getattr(fc, f'sc_mod_{key}', None)
            if fn is not None:
                fn(list(sites))
                used['measured'] += sites
                continue
            continue                                     # old ViennaRNA without these parameters: parent base
        e = BY_CODE[key]
        letter = _LETTERS[gen_k % len(_LETTERS)]
        gen_k += 1
        fc.sc_mod_json(generated_json(e, letter, e['build']), list(sites))
        used['estimated'] += sites
    for k in used:
        used[k] = sorted(used[k])
    return used


def fold_compound(seq, mods, md=None, options=None):
    """ViennaRNA fold compound of a (build) sequence with its modifications applied."""
    import RNA
    vseq = vienna_sequence(seq, mods)
    fc = RNA.fold_compound(vseq, md) if options is None else RNA.fold_compound(vseq, md, options)
    info = apply_to_fold_compound(fc, mods) if mods else {}
    return fc, info


def blocked_positions(mods):
    return {m['pos'] - 1 for m in mods if (BY_CODE[m['code']]['energy'] or ('',))[0] == 'blocked'}


def split_pairs(seq, pairs, mods):
    """Pairs of a 2D hypothesis (0-based) -> (pairs the engine builds as helices, pairs only a modification allows
    (non-canonical for the parent bases), pairs dropped because a nucleotide is blocked)."""
    blocked = blocked_positions(mods)
    allowed = {}
    for m in mods:
        e = BY_CODE[m['code']]
        if e['energy'] and e['energy'][0] == 'pairs':
            allowed[m['pos'] - 1] = set(e['energy'][1])
        elif e['energy'] and e['energy'][0] == 'vienna' and e['code'] in ('I', 'Im'):
            allowed[m['pos'] - 1] = {'C', 'U'}
    rigid, special, dropped = [], [], []
    for i, j in pairs:
        if i in blocked or j in blocked:
            dropped.append((i, j))
        elif seq[i] + seq[j] in CANON:
            rigid.append((i, j))
        elif seq[j] in allowed.get(i, ()) or seq[i] in allowed.get(j, ()):
            special.append((i, j))
        else:
            dropped.append((i, j))
    return rigid, special, dropped


def special_pair_restraints(seq, pairs, special):
    """C4'-C4' distance restraints that hold a modification-only pair (8-oxoG(syn)-A, lysidine-A ...) in its helix:
    the distances of an ideal A-form stem in which that pair is a G-U / G-C stand-in."""
    if not special:
        return None
    from .pseudoknots import helix_restraints
    from .fragments import merge_restraints
    allp = sorted(set(pairs) | set(special))
    stand = list(seq)
    for i, j in special:
        stand[i], stand[j] = ('G', 'U') if seq[j] in 'AU' else ('G', 'C')
    stand = ''.join(stand)
    out = []
    for i, j in special:
        run = [(a, b) for a, b in allp if abs(a - i) <= 2 and abs(b - j) <= 2 and (a - i) == -(b - j)]
        run = sorted(run)
        if len(run) < 2:
            run = [(i, j), (i + 1, j - 1)] if j - i > 4 else [(i, j)]
        if len(run) >= 2:
            out.append(helix_restraints(stand, run, weight=1.0))
    return merge_restraints(*out) if out else None


# ------------------------------------------------------------------------------------------------- 3D atoms
def _ring(at, r, p, q, bond):
    v = (at[r] - at[p]) + (at[r] - at[q])
    return at[r] + bond * v / np.linalg.norm(v)


# Distances in the bicyclic sugar of a locked nucleotide (PDB component LCA; O2' and the methylene to the ring atoms)
LNA_DIST = {("O2'", "C1'"): 2.375, ("O2'", "C2'"): 1.455, ("O2'", "C3'"): 2.34, ("O2'", "O4'"): 2.82, ("O2'", "O3'"): 2.78,
            ("X", "C1'"): 2.73, ("X", "C2'"): 2.33, ("X", "C3'"): 2.38, ("X", "O4'"): 2.375, ("X", "C5'"): 2.63,
            ("X", "O3'"): 2.83}


def _bridge(at, a, b, da, db):
    """LNA: positions of O2' (a) and of the methylene between O2' and C4' (b). The 2'-O,4'-C bridge is shorter than
    the O2'...C4' distance of a ribose, so O2' moves with it (about 0.6 Å on a C3'-endo sugar): both atoms are
    fitted to the distances of the bicyclic sugar, the three bonds of the bridge weighted most, the ring atoms
    kept. On a sugar far from C3'-endo no such bridge exists; O2' then stays and the methylene takes the point at
    da from it that is nearest to db from C4' (strained, not overlapping)."""
    ring = ("C1'", "C2'", "C3'", "O4'", "C5'", "O3'")
    if all(k in at for k in ring) and "C2'" in at:
        fixed = {k: np.asarray(at[k], float) for k in ring + (b,)}
        o0 = np.asarray(at[a], float)
        out_dir = o0 + fixed[b] - fixed["C2'"] - fixed["C3'"]
        x = np.concatenate([o0, 0.5 * (o0 + fixed[b]) + 0.6 * out_dir / np.linalg.norm(out_dir)])
        terms = [(0, k2, d, 1.0 if k2 != "C2'" else 3.0) for (k1, k2), d in LNA_DIST.items() if k1 == a and k2 in fixed]
        terms += [(3, k2, d, 1.0) for (k1, k2), d in LNA_DIST.items() if k1 == 'X' and k2 in fixed]
        terms.append((3, b, db, 3.0))

        def residual(x):
            r = [w * (np.linalg.norm(x[o:o + 3] - fixed[k]) - d) for o, k, d, w in terms]
            r.append(3.0 * (np.linalg.norm(x[:3] - x[3:]) - da))
            return np.array(r)
        for _ in range(40):
            r = residual(x)
            jac = np.array([(residual(x + h) - r) / 1e-5 for h in 1e-5 * np.eye(6)]).T
            step = np.linalg.lstsq(jac.T @ jac + 1e-6 * np.eye(6), -jac.T @ r, rcond=None)[0]
            x = x + np.clip(step, -0.3, 0.3)
            if np.linalg.norm(step) < 1e-5:
                break
        bonds = (np.linalg.norm(x[:3] - fixed["C2'"]) - 1.455, np.linalg.norm(x[:3] - x[3:]) - da,
                 np.linalg.norm(x[3:] - fixed[b]) - db)
        if max(abs(v) for v in bonds) < 0.12:
            return x[:3], x[3:]
    A, B = at[a], at[b]
    others = np.array([at[k] for k in ("C1'", "C2'", "C3'", "O3'", "O4'", "C5'") if k in at])
    rng = np.random.default_rng(0)
    v = rng.normal(size=(600, 3))
    pts = A + da * v / np.linalg.norm(v, axis=1)[:, None]
    score = (np.linalg.norm(pts - B, axis=1) - db) ** 2
    if len(others):
        d = np.linalg.norm(pts[:, None, :] - others[None, :, :], axis=2)
        score = score + np.sum(np.clip(2.35 - d, 0, None) ** 2, axis=1) * 4.0
    return np.asarray(A, float), pts[int(np.argmin(score))]


def _min_dist(x, env):
    if env is None or not len(env):
        return 99.0
    return float(np.sqrt(np.min(np.sum((env - x) ** 2, axis=1))))


def _apply_ops(at, ops, choice, flip):
    """Apply the 3D operations to a residue's atoms (in place). choice: dihedral per rotatable z-operation (None:
    the first). Returns [(new atom, its reference atoms)]."""
    from .geomutil import place
    new, k = [], 0
    for op in ops:
        kind = op[0]
        if kind == 'psi':
            old = {n: at.pop(n) for n in list(at) if n in PSI_MAP}
            for n, x in old.items():
                at[PSI_MAP[n]] = x
        elif kind == 'del':
            at.pop(op[1], None)
        elif kind == 'ren':
            if op[1] in at:
                at[op[2]] = at.pop(op[1])
        elif kind == 'swap':
            _, old, name, anchor, bond = op
            if old in at and anchor in at:
                v = at.pop(old) - at[anchor]
                at[name] = at[anchor] + bond * v / np.linalg.norm(v)
                new.append((name, (anchor,)))
        elif kind == 'ring':
            _, name, r, p, q, bond = op
            if all(n in at for n in (r, p, q)):
                at[name] = _ring(at, r, p, q, bond)
                new.append((name, (r, p, q)))
        elif kind == 'bridge':
            _, name, a, b, da, db = op
            if a in at and b in at:
                at[a], at[name] = _bridge(at, a, b, da, db)
                new.append((name, (a, b, "C2'", "C3'")))
        elif kind == 'z':
            name, c, b, a, bond, ang, dihs = op[1:8]
            rotatable = len(dihs) > 1
            if not all(n in at for n in (a, b, c)):
                k += rotatable
                continue
            order = list(dihs)[::-1] if (len(op) > 8 and op[8] and flip) else list(dihs)
            d = order[0]
            if rotatable:
                if choice is not None and choice[k] is not None:
                    d = choice[k]
                k += 1
            at[name] = np.array(place(at[a], at[b], at[c], bond, ang, d))
            new.append((name, (a, b, c)))
    return new


def _rotamer_options(ops, flip):
    out = []
    for op in ops:
        if op[0] == 'z' and len(op[7]) > 1:
            out.append(list(op[7])[::-1] if (len(op) > 8 and op[8] and flip) else list(op[7]))
    return out


def _placement_penalty(at, new, env):
    """Overlap of the new atoms with the other residues and with the atoms of their own residue that are not their
    bonded neighbours: sum of (3.0 Å - d)^2 over the closer contacts (own atoms with a 0.6 Å allowance for 1-3 and
    1-4 neighbours). 0 when the new atoms meet nothing."""
    pen = 0.0
    names = list(at)
    X = np.array([at[n] for n in names])
    E = env if env is not None and len(env) else np.zeros((0, 3))
    for name, refs in new:
        x = at[name]
        children = {n for n, r in new if name in r}
        mask = np.array([n not in refs and n != name and n not in children for n in names])
        d_env = np.linalg.norm(E - x, axis=1) if len(E) else np.zeros(0)
        d_own = np.linalg.norm(X[mask] - x, axis=1) + 0.6 if mask.any() else np.zeros(0)
        d = np.concatenate([d_env, d_own])
        pen += float(np.sum(np.clip(3.0 - d, 0.0, None) ** 2))
    return pen


def _place(at, ops, env, flip):
    """Apply ops choosing the side-chain rotamers: the preferred combination when its atoms meet nothing, else the
    combination (of the listed ones, then of a 30-degree scan for one- and two-rotor groups) with the least overlap."""
    import itertools
    opts = _rotamer_options(ops, flip)
    combos = list(itertools.product(*opts)) if opts else [None]
    if len(combos) > 400:
        combos = combos[:400]
    best, best_p, best_new = None, None, []
    for combo in combos:
        trial = dict(at)
        new = _apply_ops(trial, ops, list(combo) if combo is not None else None, flip)
        pen = _placement_penalty(trial, new, env)
        if best_p is None or pen < best_p - 1e-9:
            best, best_p, best_new = trial, pen, new
        if pen == 0.0:
            break
    if best_p and best_p > 0.5 and 1 <= len(opts) <= 2:
        grid = [float(v) for v in range(-180, 180, 30)]
        for combo in itertools.product(*([grid] * len(opts))):
            trial = dict(at)
            new = _apply_ops(trial, ops, list(combo), flip)
            pen = _placement_penalty(trial, new, env)
            if pen < 0.7 * best_p:
                best, best_p, best_new = trial, pen, new
    at.clear()
    at.update(best)
    return best_new


def _relieve_chi(at, new_atoms, env, rotation_about_axis):
    """An unpaired base whose new atoms meet other atoms is turned about its glycosidic bond (up to 40 degrees) to the
    position with the least overlap (loop bases are flexible; paired bases are never moved)."""
    base = [k for k in at if k not in BACKBONE and k not in ('CM2', "C2''", "F2'", 'S1P')]
    if not base:
        return
    glyc = min(base, key=lambda k: float(np.linalg.norm(at[k] - at["C1'"])))
    sub = [k for k in new_atoms if k in base]
    sugar = np.array([at[k] for k in at if k in BACKBONE and k != "C1'"])
    E = env if env is not None and len(env) else np.zeros((0, 3))

    def penalty(pos):
        pen = 0.0
        for k, x in pos.items():
            if k == glyc:
                continue
            if len(E):
                pen += float(np.sum(np.clip(3.0 - np.linalg.norm(E - x, axis=1), 0, None) ** 2))
            if k in sub:
                pen += float(np.sum(np.clip(3.0 - (np.linalg.norm(sugar - x, axis=1) + 0.6), 0, None) ** 2))
        return pen
    cur = {k: at[k] for k in base}
    p0 = penalty(cur)
    if p0 < 0.05:
        return
    best, best_p = None, p0
    axis = at[glyc] - at["C1'"]
    for ang in (-40, -30, -20, -10, 10, 20, 30, 40):
        rot = rotation_about_axis(axis, math.radians(ang))
        pos = {k: (x - at[glyc]) @ rot.T + at[glyc] for k, x in cur.items()}
        pp = penalty(pos)
        if pp < best_p - 0.05 and pp < 0.7 * p0:
            best, best_p = pos, pp
    if best is not None:
        at.update(best)


def decorate(seq, residues, mods, pairs=(), sampled=False):
    """Place the atoms of the modifications on an all-atom model of the parent bases.

    residues: list of {atom: xyz} per position (copied, not changed). pairs: 0-based pairs of the model (to choose
    the substituent orientation that a pair needs). Returns (residues, comps, report) where comps are the residue
    names to write (PDB chemical component, or the parent letter when there is none; a residue with a component
    carries that component's atom names, see COMPONENTS) and report lists the atoms placed, the residues turned syn and close contacts (< 2.2 Å) of the new atoms. sampled=True: the model comes from
    a search that already had the modification in it (pucker, syn/anti, contacts), so the glycosidic conformation it
    chose is kept (no syn flip afterwards)."""
    from .geomutil import rotation_about_axis
    out = [{a: np.asarray(x, float).copy() for a, x in r.items()} for r in residues]
    comps = list(seq)
    paired = {i for p in pairs for i in p}
    by_pos = {}
    for m in mods:
        by_pos.setdefault(m['pos'] - 1, []).append(BY_CODE[m['code']])
    report = {'placed_atoms': 0, 'syn': [], 'contacts': [], 'core_only': [], 'pucker_preferences': []}
    for i in sorted(by_pos):
        ents = sorted(by_pos[i], key=lambda e: e['kind'] == 'backbone')
        at = out[i]
        pts = [x for k, r in enumerate(out) if k != i for x in r.values()]
        env = np.array(pts) if pts else None
        ops = []
        for e in ents:
            if e['kind'] != 'backbone':
                comps[i] = e['ccd'] or e['build']
            if e['core_only']:
                report['core_only'].append(i + 1)
            if e['pucker'] and not sampled:
                report['pucker_preferences'].append({'residue': i + 1, 'code': e['code'], 'pucker': e['pucker']})
            ops += list(e['atoms'])
        try:
            new_atoms = [n for n, _r in _place(at, ops, env, i in paired)]
        except (KeyError, ValueError, ZeroDivisionError, FloatingPointError):
            new_atoms = []
        if not sampled and any(e['syn'] for e in ents) and i not in paired and 'N9' in at and "C1'" in at:
            # an 8-substituted purine meets its own sugar-phosphate in anti: turn it syn when nothing is in the way
            base = [k for k in at if k not in BACKBONE and k not in ('CM2', "C2''", "F2'", 'S1P')]
            sub8 = [k for k in base if k in ('O8', 'C8M')]
            sugar = np.array([at[k] for k in at if k in BACKBONE and k not in ("C1'", "C2'", "O4'")])
            rot = rotation_about_axis(at['N9'] - at["C1'"], math.pi)
            moved = {k: (at[k] - at['N9']) @ rot.T + at['N9'] for k in base}

            def score(pos):
                far = [x for k, x in pos.items() if k not in ('N9', 'C8', 'C4')]
                d_env = min(_min_dist(x, env) for x in far)
                d_own = min((_min_dist(pos[k], sugar) for k in sub8), default=9.0)
                return min(d_env, d_own + 0.6)
            if score(moved) >= min(2.8, score({k: at[k] for k in base}) + 0.2):
                at.update(moved)
                report['syn'].append(i + 1)
        if i not in paired and new_atoms and "C1'" in at:
            _relieve_chi(at, new_atoms, env, rotation_about_axis)
        backbone_own = np.array([at[k] for k in at if k in BACKBONE and k not in ("C1'", "C2'", "O2'", "C3'", "C4'", 'P')])
        for name in new_atoms:
            if name in at:
                d = _min_dist(at[name], env)
                if name not in ("C2''", "F2'", 'CM2', 'S1P'):
                    d = min(d, _min_dist(at[name], backbone_own))
                if d < 2.2:
                    report['contacts'].append({'residue': i + 1, 'atom': name, 'distance': round(d, 2)})
        report['placed_atoms'] += len(new_atoms)
    for i, names in written_names(mods).items():        # the dictionary's atom names, last: all of the above
        out[i - 1] = {names.get(a, a): x for a, x in out[i - 1].items()}      # works on RCV's names
        for c in report['contacts']:
            if c['residue'] == i:
                c['atom'] = names.get(c['atom'], c['atom'])
    return out, comps, report


def written_names(mods):
    """{position (1-based): {RCV atom name: name written to files}} for the modified residues whose PDB chemical
    component names its atoms differently."""
    out = {}
    for m in mods or ():
        e = BY_CODE.get(m['code'])
        if e is not None and e['kind'] != 'backbone' and e['ccd'] and e['names']:
            out[int(m['pos'])] = e['names']
    return out


# ---------------------------------------------------------------------------- the modification in the 3D search
# How strongly the sampler follows a modification's sugar pucker: the probability that a library conformer drawn for a
# non-helical residue is taken from the preferred pucker (the unmodified library is 94 % C3'-endo, 6 % C2'-endo).
# Helical residues keep the A-form C3'-endo of their stem (DNA nucleotides in an RNA duplex adopt it too).
PUCKER_STRENGTH = {'locked': 1.0, 'fluoro': 0.9, 'thio2': 0.85, 'D': 0.7, 'deoxy': 0.6, 'default': 0.75, 'psi': 0.6}
# Preferred glycosidic conformation: syn for 8-substituted purines (the 8-substituent meets the sugar in anti), anti
# reinforced for locked and 2'-fluoro sugars (C3'-endo couples to anti). (probability that a drawn conformer takes the
# preferred mode, preference in kT in the chi prior of the search).
SYN_STRENGTH = {'m8A': (0.85, 3.0), 'm8G': (0.85, 3.0), '8BrG': (0.85, 3.0)}
SYN_DEFAULT = (0.6, 1.5)
ANTI_STRENGTH = (0.5, 1.0)
SYN_CHI, ANTI_CHI = 62.0, -160.0


def sampling_preferences(mods):
    """Per residue (0-based): {'pucker': {r: ('C3endo'|'C2endo', p)}, 'chi': {r: ('syn'|'anti', p, kT)}}."""
    puck, chi = {}, {}
    for m in mods or ():
        e = BY_CODE.get(m['code'])
        if e is None or e['kind'] == 'backbone':
            continue
        r = m['pos'] - 1
        pk = e['pucker'] or ''
        if pk:
            if 'locked' in pk:
                p = PUCKER_STRENGTH['locked']
            elif e['code'] == 'D':
                p = PUCKER_STRENGTH['D']
            elif e['code'].startswith('d'):
                p = PUCKER_STRENGTH['deoxy']
            elif e['code'].startswith('f'):
                p = PUCKER_STRENGTH['fluoro']
            elif 's2' in e['code']:
                p = PUCKER_STRENGTH['thio2']
            elif ('psi',) in e['atoms'] and not e['code'].endswith('m'):
                p = PUCKER_STRENGTH['psi']
            else:
                p = PUCKER_STRENGTH['default']
            puck[r] = ('C2endo' if pk.startswith("C2'") else 'C3endo', p)
        if e['syn']:
            chi[r] = ('syn', *SYN_STRENGTH.get(e['code'], SYN_DEFAULT))
        elif 'locked' in pk or e['code'].startswith('f'):
            chi[r] = ('anti', *ANTI_STRENGTH)
    return {'pucker': puck, 'chi': chi}


def _frame(a, b, c):
    e1 = b - a
    e1 = e1 / np.linalg.norm(e1)
    v = c - a
    e2 = v - np.dot(v, e1) * e1
    e2 = e2 / np.linalg.norm(e2)
    return np.array([e1, e2, np.cross(e1, e2)])


BASE_FRAME = {'A': ('N9', 'C4', 'C8'), 'G': ('N9', 'C4', 'C8'), 'C': ('N1', 'C2', 'C6'), 'U': ('N1', 'C2', 'C6')}
SUGAR_FRAME = ("C1'", "C2'", "C3'")
PHOS_FRAME = ('P', 'OP1', 'OP2')


def rider_templates(rb, seq, mods, paired=()):
    """The modification atoms as riders of the sampled model: each new atom in the local frame of the rigid group it
    hangs on (base ring, sugar ring for each pucker, or phosphate), so that the search can place it from the parent
    atoms at any moment. Rotamers are chosen on the ideal nucleotide (the substituent orientation a pair needs when the
    residue is paired). Returns (riders, ghosts): riders [{'res', 'name', 'group', 'frame', 'local': {pucker: xyz}}],
    ghosts [(res, atom)]: parent atoms the modification removes or replaces (O2' of a deoxynucleotide, N2 of inosine,
    O2 of 2-thio-U), left out of the contact scoring."""
    from .build import BASE_ATOMS, SUGAR_ATOMS
    paired = set(paired)
    by_pos = {}
    for m in mods or ():
        by_pos.setdefault(m['pos'] - 1, []).append(BY_CODE[m['code']])
    riders, ghosts = [], []
    chi0 = rb.g.aform['chi']
    for r in sorted(by_pos):
        if not 0 <= r < len(seq):
            continue
        b = seq[r]
        ops, phos = [], []
        for e in by_pos[r]:
            for op in e['atoms']:
                (phos if op[0] == 'swap' and op[3] == 'P' else ops).append(op)
        for op in phos:                          # phosphorothioate: S along P-OP1
            if r > 0:
                riders.append({'res': r, 'name': op[2], 'group': 'phosphate', 'frame': PHOS_FRAME,
                               'local': {pk: np.array([op[4], 0.0, 0.0]) for pk in ('C3endo', 'C2endo')}})
                ghosts.append((r, op[1]))
        if not ops:
            continue
        base_set, sugar_set = set(BASE_ATOMS[b]), set(SUGAR_ATOMS)
        per_pucker, new_ref, ghost_names = {}, {}, set()
        try:
            for pk in ('C3endo', 'C2endo'):
                parent = {k: np.array(v, float) for k, v in rb.unit(b, pk, chi0).items()}
                trial = dict(parent)
                new = _place(trial, ops, None, r in paired)
                per_pucker[pk] = (parent, trial)
                for name, refs in new:
                    if name not in parent:
                        new_ref[name] = refs
                ghost_names |= set(parent) - set(trial)
        except (KeyError, ValueError, ZeroDivisionError, FloatingPointError):
            continue
        for name, refs in new_ref.items():
            root, seen = refs[0], 0
            while root in new_ref and seen < 10:
                root, seen = new_ref[root][0], seen + 1
            if root in base_set:
                group, frame = 'base', BASE_FRAME[b]
            elif root in sugar_set:
                group, frame = 'sugar', SUGAR_FRAME
            else:
                continue
            local = {}
            for pk, (parent, trial) in per_pucker.items():
                if name not in trial:
                    break
                o = parent[frame[0]]
                local[pk] = _frame(o, parent[frame[1]], parent[frame[2]]) @ (trial[name] - o)
            if len(local) == 2:
                riders.append({'res': r, 'name': name, 'group': group, 'frame': frame, 'local': local})
        ghosts += [(r, nm) for nm in sorted(ghost_names)]
    return riders, ghosts


def conformation_report(seq, residues, mods):
    """Sugar pucker and glycosidic conformation of each modified residue in a model, against its preference:
    [{'residue', 'code', 'pucker', 'pucker_preferred', 'chi', 'glycosidic', 'glycosidic_preferred'}]."""
    from .geomutil import dihedral, pseudorotation
    prefs = sampling_preferences(mods)
    out = []
    for m in mods or ():
        e = BY_CODE.get(m['code'])
        r = m['pos'] - 1
        if e is None or e['kind'] == 'backbone' or not 0 <= r < len(residues):
            continue
        at = residues[r]
        row = {'residue': r + 1, 'code': m['code']}
        ring = ("C4'", "O4'", "C1'", "C2'", "C3'")
        if all(k in at for k in ring):
            ph, _ = pseudorotation([np.asarray(at[k], float) for k in ring])
            row['pucker'] = "C3'-endo" if (ph < 50 or ph > 340) else ("C2'-endo" if 130 <= ph <= 200 else 'other')
        if r in prefs['pucker']:
            row['pucker_preferred'] = "C2'-endo" if prefs['pucker'][r][0] == 'C2endo' else "C3'-endo"
        b = seq[r]
        g, ref = ('N9', 'C4') if b in 'AG' else ('N1', 'C2')
        if all(k in at for k in ("O4'", "C1'", g, ref)):
            chi = float(dihedral(*[np.asarray(at[k], float) for k in ("O4'", "C1'", g, ref)]))
            row['chi'] = round(chi, 1)
            row['glycosidic'] = 'syn' if -20.0 <= chi <= 120.0 else 'anti'
        if r in prefs['chi']:
            row['glycosidic_preferred'] = prefs['chi'][r][0]
        out.append(row)
    return out


# -------------------------------------------------------------------------------- reading modified residues back
CCD_TO_CODE, CCD_PARENTS, CCD_NAMES = {}, {}, {}
for _ent in ENTRIES:
    for _k, _names in ((_ent['ccd'], _ent['names']), *_ent['also'].items()):
        if _k and _ent['kind'] != 'backbone':
            CCD_TO_CODE.setdefault(_k, _ent['code'])
            CCD_PARENTS.setdefault(_k, _ent['build'])
            CCD_NAMES.setdefault(_k, {v: k for k, v in _names.items()})       # component's atom name -> RCV's
PSI_CODES = {c for c, e in BY_CODE.items() if ('psi',) in e['atoms']}
_BACK = {'S2': 'O2', 'S4': 'O4', 'S6': 'O6', "F2'": "O2'", 'C7': 'N7', 'S1P': 'OP1'}


def canonical_atoms(atoms, comp=None, code=None):
    """Atom names of a modified residue mapped back to its parent base (for pair annotation and validation):
    the component's own atom names first become RCV's (the fluorine of AF2/GF2 is F, RCV's F2'), then
    pseudouridine-type C-glycosides are permuted back to uridine names, S/F replacements back to O, C7 to N7."""
    back = CCD_NAMES.get(str(comp or '').upper()) or CCD_NAMES.get(str(comp or ''))
    if back and any(k in atoms for k in back):
        atoms = {back.get(k, k): x for k, x in atoms.items()}
    code = code or CCD_TO_CODE.get(str(comp or '').upper()) or CCD_TO_CODE.get(str(comp or ''))
    if code is None and str(comp or '').upper() not in ('PSU',):
        if not any(k in atoms for k in _BACK):
            return atoms
    out = dict(atoms)
    if code in PSI_CODES or str(comp or '').upper() == 'PSU':
        inv = {v: k for k, v in PSI_MAP.items()}
        ring = {k: out.pop(k) for k in list(out) if k in inv}
        for k, x in ring.items():
            out[inv[k]] = x
    for k, v in _BACK.items():
        if k in out and v not in out:
            out[v] = out.pop(k)
    return out


def catalogue():
    """The catalogue for the interface (JSON-serialisable)."""
    rows = []
    for e in ENTRIES:
        en = e['energy']
        if en is None:
            pairing, partners = 'unchanged', None
        elif en[0] == 'blocked':
            pairing, partners = 'blocked', []
        elif en[0] == 'vienna':
            pairing = 'measured'
            partners = {'m6A': ['U'], 'pseudouridine': ['A'], 'inosine': ['C', 'U'], 'dihydrouridine': ['A'],
                        'purine': ['U'], '7DA': ['U']}[en[1]]
        else:
            pairing = 'estimated'
            partners = sorted(en[1])
        rows.append({'code': e['code'], 'name': e['name'], 'parent': e['parent'], 'build': e['build'],
                     'category': e['category'], 'ccd': e['ccd'], 'ccd_also': sorted(e['also']),
                     'aliases': list(e['aliases']), 'pairing': pairing,
                     'partners': partners, 'ddg': ({p: v[1] for p, v in en[1].items()} if en and en[0] == 'pairs' else None),
                     'syn': e['syn'], 'pucker': e['pucker'], 'evidence': e['evidence'], 'note': e['note'],
                     'core_only': e['core_only'], 'kind': e['kind'], 'atoms3d': bool(e['atoms'])})
    return {'categories': CATEGORIES, 'entries': rows}
