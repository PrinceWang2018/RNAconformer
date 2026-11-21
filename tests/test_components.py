"""Modified nucleotides against the wwPDB chemical component dictionary (0.7.6).

Each catalogue entry that is written to files as a PDB chemical component (modifications.COMPONENTS) must be that
component: the residue the engine builds has the component's atoms under the component's names, the same bonds
and the same configuration at every stereocentre. tests/data/pdb_components.json holds the dictionary's
definitions (atoms, bonds, ideal coordinates) of the components concerned."""
import itertools
import json
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from rnafold3d import modifications as M                                         # noqa: E402
from rnafold3d.structure_io import element_of                                    # noqa: E402

REF = json.loads((Path(__file__).parent / 'data' / 'pdb_components.json').read_text())
SEQ = 'GACUGGAAACAGUC'
PAIRS = [(0, 13), (1, 12), (2, 11), (3, 10), (4, 9)]
POS = {'A': 2, 'C': 3, 'U': 4, 'G': 5}                   # helical residues (C3'-endo) of the 5' strand, 1-based
BOND = {'S': 2.15, 'P': 2.15, 'BR': 2.15, 'I': 2.35}
OLD = {'O1P': 'OP1', 'O2P': 'OP2'}                       # phosphate oxygens of components defined before 2007
# components that are read as a modification although they are another stereoisomer or a reduced form of it
NOT_THE_SAME_ISOMER = {('Q', 'QUO'), ('acp3U', 'XSX')}


def _volume(x, centre, neighbours):
    v = [(x[n] - x[centre]) / np.linalg.norm(x[n] - x[centre]) for n in neighbours[:3]]
    return float(np.linalg.det(np.array(v)))


class TestPdbComponents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from rnafold3d.fold import FoldModel
        m = FoldModel(SEQ, PAIRS)
        m.sample_torsions(np.random.default_rng(1))
        m.set_loops_active(True)
        m.build()
        cls.base = m.residues_atoms()

    def built(self, code):
        """(atoms under RCV's names, bonds, residue name written, atoms as written) of one modified residue."""
        e = M.BY_CODE[code]
        pos = POS[e['build']]
        out, comps, _ = M.decorate(SEQ, self.base, M.normalize([{'pos': pos, 'code': code}], SEQ), PAIRS)
        back = {v: k for k, v in e['names'].items()} if e['ccd'] else {}
        at = {back.get(k, k): np.asarray(x, float) for k, x in out[pos - 1].items()}
        made = {}
        for op in e['atoms']:
            if op[0] in ('z', 'ring'):
                made[op[1]] = [op[2]]
            elif op[0] == 'swap':
                made[op[2]] = [op[3]]
            elif op[0] == 'bridge':
                made[op[1]] = [op[2], op[3]]
        el = {a: element_of(a) for a in at}
        parent = sorted(a for a in at if a not in made)
        bonds = {(a, b) for i, a in enumerate(parent) for b in parent[i + 1:]
                 if np.linalg.norm(at[a] - at[b]) < max(BOND.get(el[a], 1.75), BOND.get(el[b], 1.75))}
        bonds |= {tuple(sorted((a, b))) for a, bs in made.items() for b in bs if a in at and b in at}
        new = sorted(a for a in made if a in at)                 # ring closures between built atoms (queuosine)
        bonds |= {(a, b) for i, a in enumerate(new) for b in new[i + 1:] if np.linalg.norm(at[a] - at[b]) < 1.7}
        return at, bonds, comps[pos - 1], out[pos - 1]

    def test_every_component_is_the_residue_built(self):
        checked = 0
        for code, (ccd, names, *also) in M.COMPONENTS.items():
            e = M.BY_CODE[code]
            if e['core_only']:
                continue                                         # wybutosine: only the guanosine core is built
            at, bonds, comp, written = self.built(code)
            self.assertEqual(comp, ccd or e['build'], code)
            for cid, nm in ([(ccd, names)] if ccd else []) + list((also[0] if also else {}).items()):
                with self.subTest(code=code, component=cid):
                    ref = REF[cid]
                    ren = lambda a: nm.get(a, a)                                         # noqa: E731
                    ref_atoms = {OLD.get(a, a): x for a, x in ref['atoms'].items()}
                    self.assertEqual({ren(a) for a in at}, set(ref_atoms))
                    self.assertEqual({tuple(sorted((ren(a), ren(b)))) for a, b in bonds},
                                     {tuple(sorted((OLD.get(a, a), OLD.get(b, b)))) for a, b in ref['bonds']})
                    for a in at:                                 # element of the name as RCV writes and reads it
                        self.assertEqual(element_of(ren(a)).upper(), ref_atoms[ren(a)].upper(), ren(a))
                    if cid == ccd:
                        self.assertEqual(set(written), set(ref_atoms))
                    if (code, cid) in NOT_THE_SAME_ISOMER:
                        continue
                    x = {OLD.get(a, a): np.array(v) for a, v in ref['xyz'].items()}
                    adj = {}
                    for a, b in bonds:
                        adj.setdefault(a, []).append(b)
                        adj.setdefault(b, []).append(a)
                    for a, nb in adj.items():
                        nb = sorted(nb)
                        if a == 'P' or len(nb) < 3:
                            continue                             # the free phosphate of a component has no hand
                        v1, v2 = _volume(at, a, nb), _volume(x, ren(a), [ren(n) for n in nb])
                        if abs(v1) > 0.3 and abs(v2) > 0.3:
                            self.assertGreater(v1 * v2, 0, f'configuration at {a} ({ren(a)} of {cid})')
                    checked += 1
        self.assertGreaterEqual(checked, 63)

    def test_side_chain_rotamers_keep_the_configuration(self):
        """No rotamer combination may turn L-threonine (t6A) or the queuosine ring into another stereoisomer."""
        for code, centres in (('t6A', {'C12': ('N11', 'C13', 'C15'), 'C13': ('C12', 'O13', 'C14')}),
                              ('Q', {'C12': ('N11', 'C13', 'C16'), 'C15': ('C14', 'C16', 'O15'),
                                     'C16': ('C12', 'C15', 'O16')})):
            e = M.BY_CODE[code]
            start = {k: np.asarray(v, float) for k, v in self.base[POS[e['build']] - 1].items()}
            signs = set()
            for combo in itertools.product(*M._rotamer_options(e['atoms'], False)):
                at = dict(start)
                M._apply_ops(at, e['atoms'], list(combo), False)
                signs.add(tuple(_volume(at, c, list(nb)) > 0 for c, nb in centres.items()))
            self.assertEqual(len(signs), 1, code)
        at, _, _, _ = self.built('Q')
        self.assertAlmostEqual(float(np.linalg.norm(at['C15'] - at['C16'])), 1.54, delta=0.04)     # the ring closes

    def test_locked_sugar(self):
        """LNA: the 2'-O,4'-C methylene makes three bonds of normal length on a C3'-endo sugar."""
        for code in ('+A', '+G', '+C', '+T'):
            at, _, comp, _ = self.built(code)
            self.assertIn(comp, ('LCA', 'LCG', 'LCC', 'TLN'))
            for a, b, d in (("O2'", "C2'", 1.455), ("O2'", "C2''", 1.44), ("C2''", "C4'", 1.56)):
                self.assertAlmostEqual(float(np.linalg.norm(at[a] - at[b])), d, delta=0.07, msg=f'{code} {a}-{b}')

    def test_reading_back(self):
        for code, (ccd, names, *also) in M.COMPONENTS.items():
            for cid in ([ccd] if ccd else []) + list(also[0] if also else ()):
                self.assertEqual(M.CCD_TO_CODE[cid], code)
                self.assertEqual(M.CCD_PARENTS[cid], M.BY_CODE[code]['build'])
                self.assertEqual(M.lookup(cid)['code'], code)
                if len(cid) > 3:                                 # a PDB-format file has three columns for the name
                    self.assertNotEqual(cid, ccd)
        for code in ('dA', 'dC', 'dG', 'dT', 'dU'):              # DNA chains of an opened file must not turn into RNA
            self.assertEqual(M.BY_CODE[code]['ccd'], '')
        for cid in ('DA', 'DC', 'DG', 'DT'):
            self.assertNotIn(cid, M.CCD_PARENTS)
        # 2'-fluoro-A as the PDB writes it (fluorine named F): the parent atoms come back for pair annotation
        at, _, comp, written = self.built('fA')
        self.assertEqual(comp, 'AF2')
        self.assertIn('F', written)
        self.assertIn("O2'", M.canonical_atoms(written, 'AF2'))
        _, _, comp, written = self.built('m1Y')                  # a pseudouridine: C-glycoside permuted back
        back = M.canonical_atoms(written, comp)
        self.assertAlmostEqual(float(np.linalg.norm(back['N1'] - back["C1'"])), 1.5, delta=0.12)

    def test_alphafold3_input_takes_the_components(self):
        """AlphaFold 3 takes modified nucleotides as PDB components: every entry with one is passed on."""
        from rcv.alphafold import input_json
        seq, mods, problems = M.parse('GG[m1Y]C[+A]C[Q]U[k2C]G[fC][m6Am]C')
        self.assertFalse(problems)
        rna = input_json('x', seq, modifications=mods)['sequences'][0]['rna']
        self.assertEqual(rna['sequence'], 'GGUCACGUCGCAC')
        self.assertEqual([(m['modificationType'], m['basePosition']) for m in rna['modifications']],
                         [('B8H', 3), ('LCA', 5), ('56B', 7), ('M3X', 9), ('CFZ', 11)])      # m6Am has no component

    def test_lysidine_carboxyl_is_carbon(self):
        names = {op[1] for op in M.K2C if op[0] in ('z', 'ring')}
        self.assertTrue(all(element_of(n) in ('C', 'N', 'O') for n in names), names)


if __name__ == '__main__':
    unittest.main()
