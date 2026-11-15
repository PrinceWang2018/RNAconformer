"""The RNA as the scoring function sees it: typed atoms, charges, and the planes of the bases.

A docking score for RNA needs three things that a protein score does not bother with, and they all come
from here:

* **The bases are planes.** Stacking on a base face is the dominant interaction of RNA-binding small
  molecules, so every nucleotide contributes a ring centre and a ring normal, and the stacking term is
  written against those rather than against atoms.
* **The backbone is charged.** One −1 per phosphate, every 6 Å along the chain. A compound's own charge
  decides much of what it does here, so the charges are explicit rather than folded into a generic polar
  term.
* **Types come from the names, not from perception.** A nucleotide is a known object: which atom donates a
  hydrogen bond, which accepts one, which carbon is apolar (C5 of a pyrimidine, and nothing else) is a
  table, not a guess. Modified and unknown residues fall back to perception by element and connectivity.

The charge model is coarse on purpose: a formal −1 per phosphate, split over the two free oxygens, plus
dipole charges on the polar atoms of each base, balanced so that everything except the phosphate is
neutral. It is for ranking compounds against each other, not for computing an energy.
"""
from __future__ import annotations

import numpy as np

from .chem import XS_RADIUS, DEFAULT_RADIUS

# Per nucleotide: donors (N/O carrying a hydrogen), acceptors, apolar carbons, and the ring atoms that make
# the base plane. Everything not listed is a polar-neighbour carbon: neither apolar nor a hydrogen bonder.
PURINE_RING = ('N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4')
PYRIMIDINE_RING = ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')
BASES = {
    'A': {'donor': ('N6',), 'acceptor': ('N1', 'N3', 'N7'), 'apolar': (), 'ring': PURINE_RING},
    'G': {'donor': ('N1', 'N2'), 'acceptor': ('N3', 'N7', 'O6'), 'apolar': (), 'ring': PURINE_RING},
    'C': {'donor': ('N4',), 'acceptor': ('O2', 'N3'), 'apolar': ('C5',), 'ring': PYRIMIDINE_RING},
    'U': {'donor': ('N3',), 'acceptor': ('O2', 'O4'), 'apolar': ('C5',), 'ring': PYRIMIDINE_RING},
    'T': {'donor': ('N3',), 'acceptor': ('O2', 'O4'), 'apolar': ('C5', 'C7'), 'ring': PYRIMIDINE_RING},
}
# the residue names that are a chain rather than a compound, used only when a file has no ATOM records
POLYMER_RESIDUES = {'A', 'G', 'C', 'U', 'T', 'DA', 'DG', 'DC', 'DT', 'PSU', '1MA', '2MG', 'M2G', '7MG',
                    'OMG', 'OMC', '5MC', '5MU', 'H2U', '4SU', 'YG'}
PARENT = {'A': 'A', 'ADE': 'A', 'G': 'G', 'GUA': 'G', 'C': 'C', 'CYT': 'C', 'U': 'U', 'URA': 'U',
          'T': 'T', 'THY': 'T', 'DA': 'A', 'DG': 'G', 'DC': 'C', 'DT': 'T', 'PSU': 'U', '1MA': 'A',
          '2MG': 'G', 'M2G': 'G', '7MG': 'G', 'OMG': 'G', 'OMC': 'C', '5MC': 'C', '5MU': 'U', 'H2U': 'U',
          '4SU': 'U', 'YG': 'G', 'I': 'A', 'N': 'A'}
# backbone, the same in every nucleotide
BACKBONE_ACCEPTOR = ('OP1', 'OP2', 'O1P', 'O2P', "O5'", "O3'", "O4'", "O2'")
BACKBONE_DONOR = ("O2'",)                       # the 2'-hydroxyl, the one donor of the backbone
PHOSPHATE_FREE = ('OP1', 'OP2', 'O1P', 'O2P')
# charges: a phosphate carries −1, split over its free oxygens; the rest is a dipole model, neutral per base
CHARGE = {'P': 0.6, 'OP1': -0.8, 'OP2': -0.8, 'O1P': -0.8, 'O2P': -0.8, "O5'": 0.0, "O3'": 0.0,
          "O4'": -0.25, "O2'": -0.35, "C1'": 0.25, "C2'": 0.35, "C3'": 0.0, "C4'": 0.0, "C5'": 0.0}
# the sugar charges above cancel and the base is balanced in _type(), so a nucleotide carries exactly the
# −1 of its phosphate: the number that matters when a charged compound comes near the backbone
BASE_CHARGE = {'N': -0.35, 'O': -0.45, 'C': 0.0}        # polar base atoms; carbons take the balance


class Receptor:
    """Typed receptor atoms and the base planes, ready for the scoring function.

    xyz         (n, 3) heavy-atom coordinates
    element     element per atom
    radius      van der Waals radius per atom
    donor/acceptor/hydrophobic/aromatic   flags per atom
    charge      partial charge per atom (one −1 per phosphate plus a neutral dipole model per base)
    residue     residue index per atom
    ring_centre, ring_normal, ring_residue   one per nucleotide base found: the plane that stacking uses
    """

    def __init__(self, xyz, element, residue, name, resname):
        self.xyz = np.asarray(xyz, float)
        self.element = [str(e).upper() for e in element]
        self.residue = np.asarray(residue, int)
        self.name = [str(a).strip() for a in name]
        self.resname = [str(r).strip().upper() for r in resname]
        n = len(self.xyz)
        self.radius = np.array([XS_RADIUS.get(e, DEFAULT_RADIUS) for e in self.element])
        self.donor = np.zeros(n, bool)
        self.acceptor = np.zeros(n, bool)
        self.hydrophobic = np.zeros(n, bool)
        self.aromatic = np.zeros(n, bool)
        self.charge = np.zeros(n)
        self.metal = np.array([e in ('MG', 'MN', 'ZN', 'CA', 'K', 'NA') for e in self.element])
        self._type()
        self._planes()

    def __len__(self):
        return len(self.xyz)

    def _type(self):
        by_res = {}
        for a in range(len(self.xyz)):
            by_res.setdefault(int(self.residue[a]), []).append(a)
        for r, atoms in by_res.items():
            base = PARENT.get(self.resname[atoms[0]], None)
            spec = BASES.get(base)
            names = {self.name[a]: a for a in atoms}
            free = [names[x] for x in PHOSPHATE_FREE if x in names]
            for a in atoms:
                nm, el = self.name[a], self.element[a]
                if self.metal[a]:
                    self.charge[a] = 2.0 if el in ('MG', 'MN', 'ZN', 'CA') else 1.0
                    continue
                if spec is None:                       # an unknown residue: by element, as a ligand would be
                    self.acceptor[a] = el in ('O', 'N')
                    self.donor[a] = el in ('O', 'N')
                    self.hydrophobic[a] = el == 'C'
                    self.charge[a] = {'O': -0.3, 'N': -0.2, 'P': 0.6}.get(el, 0.0)
                    continue
                self.acceptor[a] = nm in BACKBONE_ACCEPTOR or nm in spec['acceptor']
                self.donor[a] = nm in BACKBONE_DONOR or nm in spec['donor']
                self.hydrophobic[a] = nm in spec['apolar']
                self.aromatic[a] = nm in spec['ring']          # the base ring: what a compound stacks on
                if nm in PHOSPHATE_FREE and free:
                    self.charge[a] = -1.6 / len(free)              # with P (+0.6): one −1 per phosphate
                elif nm in CHARGE:
                    self.charge[a] = CHARGE[nm]
                elif not nm.endswith("'"):
                    self.charge[a] = BASE_CHARGE.get(el, 0.0)
            # the base, phosphate apart, is neutral: the carbons take the balance of the dipole charges
            base_atoms = [a for a in atoms if not self.name[a].endswith("'") and self.name[a] not in
                          ('P',) + PHOSPHATE_FREE]
            carbons = [a for a in base_atoms if self.element[a] == 'C']
            if base_atoms and carbons:
                self.charge[carbons] -= self.charge[base_atoms].sum() / len(carbons)

    def _planes(self):
        """One plane per nucleotide base: the centre of its ring atoms and the normal of that ring."""
        centres, normals, res = [], [], []
        by_res = {}
        for a in range(len(self.xyz)):
            by_res.setdefault(int(self.residue[a]), []).append(a)
        for r, atoms in sorted(by_res.items()):
            spec = BASES.get(PARENT.get(self.resname[atoms[0]], None))
            if spec is None:
                continue
            ring = [a for a in atoms if self.name[a] in spec['ring']]
            if len(ring) < 5:
                continue
            p = self.xyz[ring]
            c = p.mean(0)
            normal = np.linalg.svd(p - c)[2][2]
            centres.append(c)
            normals.append(normal / (np.linalg.norm(normal) or 1.0))
            res.append(r)
        self.ring_centre = np.array(centres).reshape(-1, 3)
        self.ring_normal = np.array(normals).reshape(-1, 3)
        self.ring_residue = np.array(res, int)

    # ------------------------------------------------------------------------------- where a receptor comes from
    @classmethod
    def from_document(cls, doc, coords, k=0):
        """From an RCV document and its coordinates (conformation k)."""
        xyz = np.asarray(coords)[k]
        res = np.asarray(doc['atoms']['residue'], int)
        resname = [doc['residues'][int(r)].get('comp') or doc['residues'][int(r)].get('name') or 'N' for r in res]
        return cls(xyz, doc['atoms']['elements'], res, doc['atoms']['names'], resname)

    @classmethod
    def from_pdb(cls, path_or_text, chain=None, models=1):
        """From the ATOM records of a PDB file (HETATM is left out: a ligand is not part of the receptor,
        except the metals, which are)."""
        from .chem import _text, _element
        xyz, el, res, name, resname = [], [], [], [], []
        seen, model = {}, 0
        for ln in _text(path_or_text).splitlines():
            if ln.startswith('MODEL'):
                model += 1
                if model > models:
                    break
            if not ln.startswith('ATOM') and not (ln.startswith('HETATM') and
                                                  ln[17:20].strip().upper() in ('MG', 'MN', 'ZN', 'K', 'NA', 'CA')):
                continue
            if ln[16] not in (' ', 'A') or (chain and ln[21] != chain):
                continue
            e = _element(ln[76:78], ln[12:16])
            if e == 'H':
                continue
            key = (ln[21], ln[22:27])
            if key not in seen:
                seen[key] = len(seen)
            xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            el.append(e)
            res.append(seen[key])
            name.append(ln[12:16].strip())
            resname.append(ln[17:20].strip())
        if not xyz:
            raise ValueError('no receptor atoms found in the file')
        return cls(xyz, el, res, name, resname)

    def within(self, centre, radius):
        """Indices of the atoms within `radius` of a point: the receptor a pocket actually sees."""
        d = np.linalg.norm(self.xyz - np.asarray(centre, float), axis=1)
        return np.nonzero(d <= radius)[0]

    def subset(self, idx):
        """A receptor of the given atoms only (the base planes are recomputed from them)."""
        idx = np.asarray(idx, int)
        return Receptor(self.xyz[idx], [self.element[a] for a in idx], self.residue[idx],
                        [self.name[a] for a in idx], [self.resname[a] for a in idx])


def read_complex(path_or_text, chain=None, min_atoms=6):
    """A PDB entry split into what binds and what is bound: `(receptor, [compounds])`.

    Getting this split right is most of what a redocking benchmark does, and it is not the file's HETATM
    flag: a modified nucleotide (a 5-bromo-uridine, a GTP at the 5′ end) is written as HETATM and is part of
    the RNA, while the ligand of interest is a separate molecule. So a HETATM residue belongs to the
    receptor when one of its atoms is within covalent distance of the chain; otherwise it is a compound,
    together with any other such residue it is bonded to (an aminoglycoside is deposited as its sugars).
    The name alone will not do it: a free adenine (the ligand of an adenine riboswitch) is `ADE` and an
    argininamide is `ARG`, both the names of residues. Only when a file has no ATOM records at all are the
    residue names used, to find the chain in the first place. Waters, ions and the usual crystallisation
    additives are dropped, and only the first model of an NMR entry is read.
    """
    from scipy.spatial import cKDTree

    from .chem import _element, _text, prepare
    AMINO = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
             'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MSE'}
    SOLVENT = {'HOH', 'DOD', 'WAT', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'PG4', 'MPD', 'ACT', 'FMT', 'DMS',
               'TRS', 'EPE', 'CIT', 'MES', 'NO3', 'ACY', 'IPA', 'BME', 'NH4', 'UNX', 'UNL'}
    IONS = {'MG', 'MN', 'ZN', 'CA', 'K', 'NA', 'CL', 'BR', 'IOD', 'CS', 'SR', 'CD', 'CO', 'NI', 'CU', 'BA',
            'RB', 'LI', 'F', 'OS', 'IR', 'PT', 'AU', 'HG', 'SM', 'EU', 'TB', 'YB', 'PB'}
    rows = []
    for ln in _text(path_or_text).splitlines():
        if ln.startswith('ENDMDL'):
            break                                     # an NMR entry: the first model is the structure
        if not ln.startswith(('ATOM  ', 'HETATM')):
            continue
        if ln[16] not in (' ', 'A') or (chain and ln[21] != chain):
            continue
        el = _element(ln[76:78], ln[12:16])
        res = ln[17:20].strip().upper()
        if el == 'H' or res in SOLVENT:
            continue
        rows.append({'het': ln.startswith('HETATM'), 'res': res, 'chain': ln[21], 'num': ln[22:27].strip(),
                     'name': ln[12:16].strip(), 'el': el,
                     'xyz': [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]})
    if not rows:
        raise ValueError('no atoms found in the file')
    xyz = np.array([r['xyz'] for r in rows])
    polymer = [k for k, r in enumerate(rows) if not r['het']]
    if not polymer:                                   # a file written entirely as HETATM: fall back to names
        polymer = [k for k, r in enumerate(rows) if r['res'] in POLYMER_RESIDUES or r['res'] in AMINO]
    groups = {}
    for k, r in enumerate(rows):
        groups.setdefault((r['res'], r['chain'], r['num']), []).append(k)
    tree = cKDTree(xyz[polymer]) if polymer else None
    rec_atoms, candidates, known = list(polymer), [], set(polymer)
    for key, idx in groups.items():
        if set(idx) <= known:
            continue
        linked = tree is not None and float(np.min(tree.query(xyz[idx])[0])) < 2.0
        if key[0] in IONS or linked:                  # an ion, or a modified residue of the chain: receptor
            rec_atoms += idx
            continue
        candidates.append((key, idx))
    # An aminoglycoside is deposited as its sugars, one residue each: residues bonded to each other but not
    # to the chain are one compound, or paromomycin would be docked as three unrelated rings.
    merged, seen_c = [], set()
    for c in range(len(candidates)):
        if c in seen_c:
            continue
        group, stack = [c], [c]
        seen_c.add(c)
        while stack:
            x = stack.pop()
            for y in range(len(candidates)):
                if y in seen_c:
                    continue
                d = cKDTree(xyz[candidates[x][1]]).query(xyz[candidates[y][1]])[0]
                if float(np.min(d)) < 2.0:
                    seen_c.add(y)
                    group.append(y)
                    stack.append(y)
        merged.append(group)
    ligands = []
    for group in merged:
        idx = sorted(k for c in group for k in candidates[c][1])
        if len(idx) < min_atoms:
            continue
        names = [candidates[c][0][0] for c in group]
        key = candidates[group[0]][0]
        m = prepare([rows[k]['el'] for k in idx], xyz[idx], None,
                    name='-'.join(names) + f'_{key[1]}{key[2]}')
        m.properties.update({'residue': names[0], 'residues': names, 'chain': key[1], 'number': key[2]})
        ligands.append(m)
    rec_atoms = sorted(set(rec_atoms))
    seen, res_index = {}, []
    for k in rec_atoms:
        res_index.append(seen.setdefault((rows[k]['chain'], rows[k]['num']), len(seen)))
    rec = Receptor(xyz[rec_atoms], [rows[k]['el'] for k in rec_atoms], res_index,
                   [rows[k]['name'] for k in rec_atoms], [rows[k]['res'] for k in rec_atoms])
    ligands.sort(key=lambda m: -len(m))
    return rec, ligands
