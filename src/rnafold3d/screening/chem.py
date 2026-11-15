"""Small molecules for the screen: read a compound, work out what the scoring function needs, nothing else.

The scoring core must not depend on RDKit (`docs/NEXT_STEPS.md` §1.6): it takes prepared atoms, types and
charges, so a prepared library works on a machine without it. This module is that preparation, written for
the few things an RNA docking score actually uses:

* **Heavy atoms only.** Hydrogens are read if present and then dropped, their count kept on the atom they
  hang from. A compound from the PDB has no hydrogens at all, so donors have to be inferred from valence in
  either case, and a united-atom model gives the same answer at a fraction of the cost.
* **Atom radii and types** in the X-Score / AutoDock Vina spirit (the Vina scoring function is Apache-2.0 and
  may be used as a reference; no GPL docking code or parameters are copied): a radius per element, and the
  three flags the terms need — hydrophobic, hydrogen-bond donor, hydrogen-bond acceptor — plus aromatic,
  which an RNA score needs for stacking and a protein score does not.
* **Formal charges** of the groups that are ionised at physiological pH (carboxylate, phosphate, sulfonate,
  protonated amine, amidinium, guanidinium) spread over the atoms that share them, and small partial charges
  on polar atoms. This is a charge model for ranking, not a force field: it exists because the RNA backbone
  is a line of −1 charges and a compound's own charge decides most of what it does with that.
* **Rotatable bonds and the torsion tree**, so a pose can be sampled by moving branches instead of atoms.

Readers: SDF/MOL (V2000, the format libraries ship in), MOL2, and PDB/mmCIF HETATM records (how a ligand
arrives when it is taken out of a complex, with no bond orders — those are inferred from the geometry).

Nothing here is chemistry for its own sake: tautomers, protonation states and ring conformers are listed in
the plan as RDKit's job, and a compound is read as it is given.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

import numpy as np

# van der Waals radii used by the scoring function (X-Score radii, as in AutoDock Vina)
XS_RADIUS = {'C': 1.9, 'N': 1.8, 'O': 1.7, 'S': 2.0, 'P': 2.1, 'F': 1.5, 'CL': 1.8, 'BR': 2.0, 'I': 2.2,
             'SE': 2.0, 'B': 1.9, 'SI': 2.1, 'MG': 1.2, 'MN': 1.2, 'ZN': 1.2, 'CA': 1.2, 'FE': 1.2, 'NA': 1.2,
             'K': 1.2, 'H': 1.0}
DEFAULT_RADIUS = 1.9
COVALENT = {'H': 0.37, 'C': 0.77, 'N': 0.75, 'O': 0.73, 'F': 0.71, 'P': 1.06, 'S': 1.02, 'CL': 0.99,
            'BR': 1.14, 'I': 1.33, 'SE': 1.16, 'B': 0.82, 'SI': 1.11, 'MG': 1.30, 'MN': 1.39, 'ZN': 1.31,
            'CA': 1.74, 'FE': 1.25, 'NA': 1.54, 'K': 1.96}
HALOGEN = ('F', 'CL', 'BR', 'I')
METALS = ('MG', 'MN', 'ZN', 'CA', 'FE', 'NA', 'K')
# how many bonds an element normally makes, to tell a donor (free valence on N/O) from a full one
VALENCE = {'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 5, 'F': 1, 'CL': 1, 'BR': 1, 'I': 1, 'SE': 2, 'B': 3, 'SI': 4}


@dataclass
class Molecule:
    """A compound as the scoring function sees it: heavy atoms, bonds, and the flags of each atom.

    elements    element symbol per heavy atom, upper case
    xyz         (n, 3) coordinates, Å
    bonds       [(i, j, order)], order 1, 2, 3 or 4 for aromatic
    hydrogens   how many hydrogens hang from each heavy atom
    charge      partial charge per atom (formal charges spread over the group, plus small polar charges)
    radius      van der Waals radius per atom
    hydrophobic carbon, halogen or sulfur with no polar neighbour: the atoms a buried-surface term likes
    donor       can give a hydrogen bond (N, O or S with a hydrogen, or with a free valence)
    acceptor    can take one (O, and N or S without a positive charge)
    aromatic    in a planar conjugated ring: what stacks against a nucleobase
    rings       [[atom indices]] of the smallest rings found
    rotatable   [(i, j)] single acyclic bonds the pose may turn about
    name        whatever the file called it
    """
    elements: list
    xyz: np.ndarray
    bonds: list = field(default_factory=list)
    hydrogens: np.ndarray = None
    charge: np.ndarray = None
    radius: np.ndarray = None
    hydrophobic: np.ndarray = None
    donor: np.ndarray = None
    acceptor: np.ndarray = None
    aromatic: np.ndarray = None
    rings: list = field(default_factory=list)
    rotatable: list = field(default_factory=list)
    name: str = ''
    properties: dict = field(default_factory=dict)

    def __len__(self):
        return len(self.elements)

    @property
    def n_rotatable(self):
        return len(self.rotatable)

    @property
    def formula(self):
        out, seen = [], {}
        for e in self.elements:
            seen[e] = seen.get(e, 0) + 1
        for e in sorted(seen, key=lambda x: (x != 'C', x != 'N', x)):
            out.append(e.capitalize() + (str(seen[e]) if seen[e] > 1 else ''))
        h = int(self.hydrogens.sum()) if self.hydrogens is not None else 0
        return ''.join(out) + (f'H{h}' if h else '')

    @property
    def weight(self):
        """Molecular weight including the hydrogens that were dropped (Da)."""
        w = {'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'P': 30.974, 'F': 18.998, 'CL': 35.45,
             'BR': 79.904, 'I': 126.90, 'SE': 78.97, 'B': 10.81, 'SI': 28.085}
        h = float(self.hydrogens.sum()) if self.hydrogens is not None else 0.0
        return round(sum(w.get(e, 12.0) for e in self.elements) + 1.008 * h, 2)

    @property
    def neighbours(self):
        if getattr(self, '_nb', None) is None or len(self._nb) != len(self.elements):
            nb = [[] for _ in self.elements]
            for i, j, _o in self.bonds:
                nb[i].append(j)
                nb[j].append(i)
            self._nb = nb
        return self._nb

    def copy_with(self, xyz):
        """The same molecule at new coordinates (a pose)."""
        m = Molecule(**{k: getattr(self, k) for k in
                        ('elements', 'bonds', 'hydrogens', 'charge', 'radius', 'hydrophobic', 'donor',
                         'acceptor', 'aromatic', 'rings', 'rotatable', 'name', 'properties')}, xyz=np.asarray(xyz, float))
        return m


# ------------------------------------------------------------------------------------------- readers
def _element(symbol, name=''):
    s = re.sub(r'[^A-Za-z]', '', str(symbol)).upper()
    if s in XS_RADIUS or s in METALS:
        return s
    s2 = re.sub(r'[^A-Za-z]', '', str(name)).upper()[:2]
    return s2 if s2 in XS_RADIUS else (s[:1] or 'C')


def read_sdf(text_or_path, limit=None):
    """Compounds in an SDF or MOL file (V2000). Records that cannot be read are skipped, with a reason."""
    text = _text(text_or_path)
    out = []
    for block in text.split('$$$$'):
        if not block.strip():
            continue
        try:
            out.append(_one_mol(block))
        except (ValueError, IndexError) as err:
            out.append(_broken(block.strip().splitlines()[0][:60] if block.strip() else '', str(err)))
        if limit and len(out) >= limit:
            break
    return out


def _one_mol(block):
    lines = block.lstrip('\n').splitlines()
    if len(lines) < 4:
        raise ValueError('not a molfile record')
    name = lines[0].strip()
    counts = lines[3]
    if 'V3000' in counts:
        raise ValueError('V3000 molfiles are not read; convert to V2000 or SDF V2000')
    na, nb = int(counts[0:3]), int(counts[3:6])
    elements, xyz, charges = [], [], []
    for k in range(na):
        ln = lines[4 + k]
        xyz.append([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])])
        elements.append(_element(ln[31:34]))
        ch = int(ln[36:39]) if len(ln) >= 39 and ln[36:39].strip() else 0
        charges.append({0: 0, 1: 3, 2: 2, 3: 1, 5: -1, 6: -2, 7: -3}.get(ch, 0))
    bonds = []
    for k in range(nb):
        ln = lines[4 + na + k]
        i, j, o = int(ln[0:3]) - 1, int(ln[3:6]) - 1, int(ln[6:9])
        bonds.append((i, j, o))
    formal = np.array(charges, float)
    for ln in lines[4 + na + nb:]:                      # M  CHG lines override the (obsolete) charge column
        if ln.startswith('M  CHG'):
            f = ln.split()[3:]
            formal[:] = 0.0
            for a, c in zip(f[0::2], f[1::2]):
                formal[int(a) - 1] = float(c)
    props = {}
    for k, ln in enumerate(lines):
        if ln.startswith('> '):
            key = re.findall(r'<([^>]+)>', ln)
            if key and k + 1 < len(lines):
                props[key[0]] = lines[k + 1].strip()
    return prepare(elements, xyz, bonds, formal=formal, name=name, properties=props)


def read_mol2(text_or_path, limit=None):
    """Compounds in a MOL2 file. The SYBYL atom types are read for the bond orders only; the flags are
    worked out here, so a file typed by another program is treated like any other."""
    text = _text(text_or_path)
    out = []
    for block in text.split('@<TRIPOS>MOLECULE')[1:]:
        lines = block.splitlines()
        name = lines[0].strip() if lines else ''
        sec, atoms, bonds = None, [], []
        for ln in lines:
            if ln.startswith('@<TRIPOS>'):
                sec = ln.strip()[9:]
                continue
            f = ln.split()
            if sec == 'ATOM' and len(f) >= 6:
                atoms.append((_element(f[5].split('.')[0], f[1]), [float(f[2]), float(f[3]), float(f[4])],
                              float(f[8]) if len(f) > 8 else 0.0))
            elif sec == 'BOND' and len(f) >= 4:
                o = {'1': 1, '2': 2, '3': 3, 'ar': 4, 'am': 1, 'du': 1, 'un': 1}.get(f[3], 1)
                bonds.append((int(f[1]) - 1, int(f[2]) - 1, o))
        if not atoms:
            continue
        out.append(prepare([a[0] for a in atoms], [a[1] for a in atoms], bonds, name=name))
        if limit and len(out) >= limit:
            break
    return out


def read_pdb_ligand(text_or_path, resname=None, min_atoms=6):
    """Ligands in the HETATM records of a PDB file (or the ATOM records of an mmCIF converted to one).

    This is how a compound arrives when it is taken out of a complex, so: water, ions and the common
    crystallisation additives are left out, each residue is one compound, and the bonds are worked out from
    the distances, since a PDB file has none."""
    text = _text(text_or_path)
    skip = {'HOH', 'DOD', 'WAT', 'SO4', 'PO4', 'GOL', 'EDO', 'PEG', 'MPD', 'ACT', 'FMT', 'DMS', 'TRS', 'EPE',
            'CIT', 'MES', 'NO3', 'CL', 'NA', 'K', 'MG', 'MN', 'ZN', 'CA', 'SR', 'CD', 'CO', 'NI', 'CU', 'IOD',
            'BR', 'F', 'SPM', 'SPD', 'PUT'}
    groups = {}
    for ln in text.splitlines():
        if not ln.startswith('HETATM'):
            continue
        res = ln[17:20].strip().upper()
        if (resname and res != str(resname).upper()) or (not resname and res in skip):
            continue
        alt = ln[16]
        if alt not in (' ', 'A'):
            continue
        key = (res, ln[21], ln[22:27].strip())
        el = _element(ln[76:78], ln[12:16])
        if el == 'H':
            continue
        groups.setdefault(key, []).append((el, [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    out = []
    for (res, chain, num), atoms in groups.items():
        if len(atoms) < min_atoms:
            continue
        m = prepare([a[0] for a in atoms], [a[1] for a in atoms], None, name=f'{res}_{chain}{num}')
        m.properties.update({'residue': res, 'chain': chain, 'number': num})
        out.append(m)
    return out


def read(path, **kw):
    """Compounds from a file, by its extension (.sdf, .mol, .mol2, .pdb, .ent)."""
    s = str(path).lower()
    if s.endswith(('.sdf', '.mol', '.sd')):
        return read_sdf(path, **kw)
    if s.endswith('.mol2'):
        return read_mol2(path, **kw)
    if s.endswith(('.pdb', '.ent')):
        return read_pdb_ligand(path, **kw)
    raise ValueError(f'unknown compound file type: {path}')


def _text(x):
    """The content of a file, or the text it already is (a path, a Path, or the file's text)."""
    from pathlib import Path
    if hasattr(x, 'read_text'):
        return x.read_text(encoding='utf-8', errors='replace')
    if isinstance(x, bytes):
        return x.decode('utf-8', 'replace')
    s = str(x)
    if '\n' not in s and len(s) < 4096:
        p = Path(s)
        if p.exists():
            return p.read_text(encoding='utf-8', errors='replace')
        raise FileNotFoundError(s)
    return s


def _broken(name, why):
    m = Molecule(elements=[], xyz=np.zeros((0, 3)), name=name)
    m.properties['error'] = why
    return m


# ------------------------------------------------------------------------------------- what the score needs
def bonds_from_distance(elements, xyz, tolerance=0.45):
    """Bonds of a molecule that came without any: two atoms are bonded when they are closer than the sum of
    their covalent radii plus a tolerance. Bond orders are guessed from the length afterwards."""
    from scipy.spatial import cKDTree
    xyz = np.asarray(xyz, float)
    rad = np.array([COVALENT.get(e, 0.77) for e in elements])
    pairs = cKDTree(xyz).query_pairs(float(2 * rad.max() + tolerance), output_type='ndarray')
    out = []
    for i, j in pairs:
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        if d <= rad[i] + rad[j] + tolerance and not (elements[i] in METALS or elements[j] in METALS):
            out.append((int(i), int(j), _order_from_length(elements[i], elements[j], d)))
    return out


def _order_from_length(a, b, d):
    """1, 2 or 3 from the bond length. Aromatic bonds come out as 2 or 1 and are fixed by the ring search."""
    key = tuple(sorted((a, b)))
    table = {('C', 'C'): (1.41, 1.25), ('C', 'N'): (1.38, 1.22), ('C', 'O'): (1.34, 1.19),
             ('N', 'N'): (1.33, 1.20), ('N', 'O'): (1.30, 1.18), ('C', 'S'): (1.71, 1.60),
             ('O', 'P'): (1.55, 1.45), ('O', 'S'): (1.55, 1.44)}
    lim = table.get(key)
    if lim is None:
        return 1
    return 1 if d > lim[0] else (2 if d > lim[1] else 3)


def find_rings(n, bonds, max_size=8):
    """The small rings of a molecule: for every bond in a cycle, the shortest cycle through it."""
    nb = [[] for _ in range(n)]
    for i, j, _o in bonds:
        nb[i].append(j)
        nb[j].append(i)
    rings, seen = [], set()
    for i, j, _o in bonds:
        path = _shortest_path(nb, i, j, max_size)           # i to j without using the bond itself
        if path is None:
            continue
        key = frozenset(path)
        if key not in seen:
            seen.add(key)
            rings.append(path)
    return sorted(rings, key=len)


def _shortest_path(nb, start, goal, max_size):
    """Shortest path from `start` to `goal` that does not use the direct bond between them."""
    from collections import deque
    q = deque([(start, [start])])
    best = None
    while q:
        at, path = q.popleft()
        if len(path) > max_size:
            continue
        for k in nb[at]:
            if k == goal and len(path) > 2:
                return path + [k]
            if k in path or (at == start and k == goal):
                continue
            q.append((k, path + [k]))
    return best


def _planar(xyz, idx, tol=0.2):
    p = np.asarray(xyz, float)[list(idx)]
    if len(p) < 4:
        return True
    c = p.mean(0)
    normal = np.linalg.svd(p - c)[2][2]
    return float(np.abs((p - c) @ normal).max()) <= tol


def prepare(elements, xyz, bonds=None, formal=None, hydrogens=None, name='', properties=None, ph=7.0):
    """Everything the scoring function asks of a compound: heavy atoms, radii, flags, charges, rings and
    rotatable bonds. Hydrogens in the input are counted onto their heavy atom and then dropped."""
    elements = [_element(e) for e in elements]
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    if len(elements) != len(xyz):
        raise ValueError('one element per atom is needed')
    stated = bonds is not None               # orders from the file, or guessed from the bond lengths below
    if bonds is None:
        bonds = bonds_from_distance(elements, xyz)
    bonds = [(int(i), int(j), int(o)) for i, j, o in bonds if i != j]
    formal = np.zeros(len(elements)) if formal is None else np.asarray(formal, float)
    h = np.zeros(len(elements), int) if hydrogens is None else np.asarray(hydrogens, int)
    # drop the hydrogens, keeping them as a count on the heavy atom they hang from
    heavy = [k for k, e in enumerate(elements) if e != 'H']
    if len(heavy) != len(elements):
        at = {old: new for new, old in enumerate(heavy)}
        for i, j, _o in bonds:
            if elements[i] == 'H' and elements[j] != 'H':
                h[j] += 1
            elif elements[j] == 'H' and elements[i] != 'H':
                h[i] += 1
        bonds = [(at[i], at[j], o) for i, j, o in bonds if i in at and j in at]
        elements = [elements[k] for k in heavy]
        xyz, formal, h = xyz[heavy], formal[heavy], h[heavy]
    n = len(elements)
    m = Molecule(elements=elements, xyz=xyz, bonds=bonds, hydrogens=h, name=name,
                 properties=dict(properties or {}))
    nb = m.neighbours
    order = {}
    for i, j, o in bonds:
        order[(i, j)] = order[(j, i)] = max(o, order.get((i, j), 0))
    m.rings = find_rings(n, bonds)
    in_ring = np.zeros(n, bool)
    for r in m.rings:
        in_ring[list(r)] = True
    # aromatic: a planar ring of five or six atoms in which no atom is saturated
    m.aromatic = np.zeros(n, bool)
    for r in m.rings:
        ring = sorted(set(r))
        if len(ring) not in (5, 6) or not _planar(xyz, ring):
            continue
        if all(elements[a] in ('C', 'N', 'O', 'S') and len(nb[a]) + h[a] <= 3 for a in ring):
            m.aromatic[ring] = True
    # flags
    m.radius = np.array([XS_RADIUS.get(e, DEFAULT_RADIUS) for e in elements])
    polar_neighbour = np.array([any(elements[k] in ('N', 'O', 'P', 'S') for k in nb[a]) for a in range(n)])
    m.hydrophobic = np.array([(e == 'C' and not polar_neighbour[a]) or e in HALOGEN
                              for a, e in enumerate(elements)])
    free = np.array([max(0, VALENCE.get(e, 4) - sum(order.get((a, k), 1) for k in nb[a]) - h[a])
                     for a, e in enumerate(elements)])
    polar = [e in ('N', 'O', 'S') for e in elements]
    if h.sum() == 0:
        # The file brought no hydrogens — the usual case for a compound taken out of a PDB entry. Which
        # tautomer it means cannot be known, so every polar atom with room for a hydrogen is treated as both
        # a donor and an acceptor. A library file with hydrogens (SDF, MOL2) is read exactly instead.
        m.donor = np.array([polar[a] and len(nb[a]) < VALENCE.get(elements[a], 3) for a in range(n)])
        m.acceptor = np.array([polar[a] and formal[a] <= 0 for a in range(n)])
    else:
        m.donor = np.array([polar[a] and (h[a] > 0 or (free[a] > 0 and not m.aromatic[a]))
                            for a in range(n)])
        m.acceptor = np.array([(elements[a] == 'O') or (polar[a] and formal[a] <= 0 and
                               sum(order.get((a, k), 1) for k in nb[a]) + h[a] <= VALENCE.get(elements[a], 3))
                               for a in range(n)])
    m.charge = _charges(m, formal, order, ph, stated)
    m.rotatable = _rotatable(m, order)
    return m


def ionised(m, order, ph=7.0, implicit_h=None, stated_orders=True):
    """Formal charge per atom at `ph` (only pH 7 is modelled), for the groups that are charged in water.

    An RNA backbone is a line of −1 charges, so whether a compound carries +2 or −1 decides most of what it
    does with it, and a compound drawn as a neutral acid in a file is not neutral at the bench. The groups:

    * carboxylate, sulfonate and phosphate: −1 per group (−2 for a phosphate with three free oxygens),
      shared over the oxygens that carry it, whether the file drew them protonated or not;
    * aliphatic amine (not an amide, nitrile or aromatic nitrogen), amidine and guanidine: +1, shared over
      the nitrogens of the group — also when the file brought no hydrogens, since the amine has one in
      water whether the file drew it or not;
    * a nitrogen with four bonds, a sulfur with three, and a nitrogen with three bonds in an aromatic ring
      that *also* holds a double bond inside that ring: a pyridinium or a thiazolium is charged, while the
      N9 of a guanine, the N10 of a flavin and the two methylated nitrogens of a theophylline are not —
      they spend their lone pair on the ring and keep their three single bonds. Where the file brought no
      bond orders at all (a compound taken from a PDB entry, orders guessed from the bond lengths) the test
      falls back to the pyrrole position: in a five-ring, a sulfur or an oxygen takes it if the ring has
      one, so a thiazolium nitrogen is charged and a guanine N9 is not; in a six-ring no atom has to, so
      nothing is charged, which is the right answer far more often than the opposite.

    An amine on an aromatic ring (an aniline, an aminopyrimidine) is left neutral: its pair is in the ring
    and it is not protonated at pH 7, which is what separates it from an aliphatic amine.

    Anything else is neutral. Tautomers, unusual pKa's and anything between pH 4 and 9 are not modelled;
    that is the preparation step RDKit does (`docs/NEXT_STEPS.md` §1.6).
    """
    n, el, nb, h = len(m), m.elements, m.neighbours, m.hydrogens
    q = np.zeros(n)
    if ph is None:
        return q
    # a file without any hydrogens (a compound from a PDB entry) still has them in water: an amine with
    # room for one is protonated, and the rules below are applied as if it carried it
    implicit = bool(h.sum() == 0) if implicit_h is None else bool(implicit_h)
    for a in range(n):
        if el[a] in ('C', 'S', 'P'):
            terminal = [k for k in nb[a] if el[k] == 'O' and len(nb[k]) == 1]
            if el[a] == 'C' and len(terminal) >= 2:                      # carboxylate
                for k in terminal[:2]:
                    q[k] = -0.5
            elif el[a] == 'S' and len(terminal) >= 3:                    # sulfonate
                for k in terminal:
                    q[k] = -1.0 / len(terminal)
            elif el[a] == 'P' and len(terminal) >= 2:                    # phosphate / phosphonate
                total = -2.0 if len(terminal) >= 3 else -1.0
                for k in terminal:
                    q[k] = total / len(terminal)
        if el[a] == 'S' and len([k for k in nb[a] if el[k] != 'H']) >= 3:
            q[a] = 1.0                                                   # sulfonium
        if el[a] == 'N' and len([k for k in nb[a] if el[k] != 'H']) + h[a] >= 4:
            q[a] = 1.0                                                   # a quaternary nitrogen
    # Does the file state bond orders at all? With Kekulé orders an aromatic ring alternates, so at least
    # one of its bonds is a double one; a set of orders guessed from bond lengths does not reach that.
    kekule = stated_orders and any(
        order.get(tuple(sorted((ring[i], ring[(i + 1) % len(ring)]))), 1) >= 2
        for ring in m.rings if m.aromatic[sorted(set(ring))].all()
        for i in range(len(ring)))
    for ring in m.rings:                                                 # pyridinium, thiazolium, imidazolium
        ring = sorted(set(ring))
        if not m.aromatic[ring].all():
            continue
        three = [a for a in ring if el[a] == 'N' and len([k for k in nb[a] if el[k] != 'H']) >= 3]
        if not three:
            continue
        if kekule:
            # the nitrogen is charged when it holds a double bond in the ring: it has spent no lone pair on
            # the ring's π system and still carries three σ-bonds
            for a in three:
                if any(order.get(tuple(sorted((a, k))), 1) >= 2 for k in nb[a] if k in ring):
                    q[a] = 1.0
            continue
        if len(ring) == 6:                                               # no pyrrole position is needed
            continue
        # In a five-ring one atom gives the ring its second pair of π electrons — the pyrrole position. A
        # sulfur or an oxygen takes it if there is one (thiazolium), otherwise one of the nitrogens does and
        # is neutral: that is why the N9 of a guanine, with its three bonds, carries no charge.
        donor = [a for a in ring if el[a] in ('S', 'O')] or three[:1]
        for a in three:
            if a not in donor:
                q[a] = 1.0
    for a in range(n):                                                   # amidinium / guanidinium
        if el[a] == 'C' and not m.aromatic[a]:
            ns = [k for k in nb[a] if el[k] == 'N' and not m.aromatic[k]]
            if len(ns) >= 2 and any(order.get((a, k), 1) >= 2 for k in ns) and all(q[k] == 0 for k in ns):
                for k in ns:
                    q[k] = 1.0 / len(ns)
    for a in range(n):                                                   # aliphatic amine
        if el[a] != 'N' or q[a] != 0 or m.aromatic[a]:
            continue
        heavy = [k for k in nb[a] if el[k] != 'H']
        if len(heavy) > 3 or (h[a] == 0 and not implicit) or any(order.get((a, k), 1) > 1 for k in heavy):
            continue
        conjugated = any(el[k] in ('C', 'S', 'P') and any(order.get((k, x), 1) >= 2 and el[x] in ('O', 'S', 'N')
                                                          for x in nb[k] if x != a) for k in heavy)
        aniline = any(m.aromatic[k] for k in heavy)                      # its pair is in the ring, pKa ≈ 5
        if not conjugated and not aniline:                               # an amide or a nitro nitrogen is neutral
            q[a] = 1.0
    return q


def _charges(m, formal, order, ph=7.0, stated=True):
    """Charge per atom: the ionised groups (from the file when it gives them, else perceived at pH 7) plus a
    small polarity charge on every other polar atom, shifted so that the polarity charges cancel and the net
    charge of the molecule is the sum of the formal ones. A charge the file states and `ionised` does not
    assign is kept; `ph=None` takes the file's charges as they are, for a library prepared elsewhere. A ranking model, not a force field: it is here because the sign and size of a
    compound's charge is most of its electrostatics against a polyanion."""
    n, el, nb, h = len(m), m.elements, m.neighbours, m.hydrogens
    q = ionised(m, order, ph, stated_orders=stated)
    formal = np.asarray(formal, float)
    keep = (formal != 0) & (q == 0)                   # a charge the file states and the rules above do not
    q[keep] = formal[keep]
    if ph is None:
        q = formal.copy()
    base = {'O': -0.38, 'N': -0.30, 'S': -0.15, 'F': -0.20, 'CL': -0.12, 'BR': -0.10, 'I': -0.08, 'P': 0.35}
    small = np.zeros(n)
    for a in range(n):
        if abs(q[a]) > 1e-6:
            continue
        if el[a] == 'C':
            small[a] = 0.14 * sum(1 for k in nb[a] if el[k] in ('N', 'O'))
        else:
            small[a] = base.get(el[a], 0.0) + (0.20 * h[a] if el[a] in ('N', 'O', 'S') else 0.0)
    free = small != 0
    if free.any():                                    # the polarity charges must not change the net charge
        small[free] -= small.sum() / free.sum()
    return q + small


def _rotatable(m, order):
    """Single bonds the pose may turn about: acyclic, with a heavy atom on each side to move, and not one of
    the bonds that do not really turn (amide C–N, a bond to a terminal group of equal atoms, a triple bond)."""
    el, nb = m.elements, m.neighbours
    ring_bond = {tuple(sorted((r[x], r[(x + 1) % len(r)]))) for r in m.rings for x in range(len(r))}
    out = []
    for i, j, o in m.bonds:
        if o != 1 or tuple(sorted((i, j))) in ring_bond:     # a bond of a ring, not merely touching one
            continue
        hi = [k for k in nb[i] if k != j]
        hj = [k for k in nb[j] if k != i]
        if not hi or not hj:
            continue                                      # a terminal bond moves nothing
        if {el[i], el[j]} == {'C', 'N'}:                   # amide / ester bond: conjugated, treated as rigid
            c, nn = (i, j) if el[i] == 'C' else (j, i)
            if any(order.get((c, k), 1) >= 2 and el[k] in ('O', 'S', 'N') for k in nb[c] if k != nn):
                continue
        if all(el[k] == el[hj[0]] and len(nb[k]) == 1 for k in hj) and len(hj) >= 2:
            continue                                      # turning a symmetric terminal group changes nothing
        if all(el[k] == el[hi[0]] and len(nb[k]) == 1 for k in hi) and len(hi) >= 2:
            continue
        out.append((int(i), int(j)))
    return out


def torsion_tree(m, root=None):
    """The compound as a rigid root with branches that turn: [(i, j, [atoms that move when i-j turns])].

    The branches are ordered so that a branch's atoms are moved after its parent's, which is what lets a
    pose be built by applying the torsions one after another (the standard flexible-ligand representation).
    """
    n = len(m)
    nb = m.neighbours
    rot = {tuple(sorted(b)) for b in m.rotatable}
    # the rigid pieces: connected components when the rotatable bonds are cut
    comp = [-1] * n
    pieces = []
    for a in range(n):
        if comp[a] >= 0:
            continue
        stack, members = [a], []
        comp[a] = len(pieces)
        while stack:
            x = stack.pop()
            members.append(x)
            for k in nb[x]:
                if comp[k] < 0 and tuple(sorted((x, k))) not in rot:
                    comp[k] = comp[a]
                    stack.append(k)
        pieces.append(sorted(members))
    if root is None:                                      # the largest piece stays still: fewest atoms moved
        root = int(np.argmax([len(p) for p in pieces]))
    # walk outwards from the root piece
    out, seen, queue = [], {root}, [root]
    while queue:
        p = queue.pop(0)
        for a in pieces[p]:
            for k in nb[a]:
                if comp[k] in seen or tuple(sorted((a, k))) not in rot:
                    continue
                child = comp[k]
                seen.add(child)
                queue.append(child)
                moved = []
                stack, done = [child], {p, child}
                while stack:
                    c = stack.pop()
                    moved += pieces[c]
                    for b in pieces[c]:
                        for x in nb[b]:
                            if comp[x] not in done:
                                done.add(comp[x])
                                stack.append(comp[x])
                out.append((int(a), int(k), sorted(moved)))
    return pieces[root], out


def apply_torsions(m, angles, tree=None, xyz=None):
    """The compound with its rotatable bonds turned by `angles` (radians), in the order of `torsion_tree`."""
    _root, branches = tree if tree is not None else torsion_tree(m)
    p = np.array(m.xyz if xyz is None else xyz, float)
    for (i, j, moved), ang in zip(branches, np.asarray(angles, float)):
        axis = p[j] - p[i]
        norm = np.linalg.norm(axis)
        if norm < 1e-9 or abs(ang) < 1e-12:
            continue
        axis = axis / norm
        c, s = math.cos(ang), math.sin(ang)
        v = p[moved] - p[j]
        p[moved] = p[j] + v * c + np.cross(axis, v) * s + np.outer(v @ axis, axis) * (1 - c)
    return p


def write_sdf(molecules, path=None, xyz=None, properties=None):
    """Write compounds (or poses of one compound) as an SDF: the format every other program reads.

    `xyz` replaces the coordinates of a single molecule — how a list of poses is written out — and
    `properties` is a list of dicts written as SDF data fields (the score of each pose, for instance).
    """
    molecules = [molecules] if isinstance(molecules, Molecule) else list(molecules)
    coords = [m.xyz for m in molecules] if xyz is None else [np.asarray(x, float) for x in xyz]
    if xyz is not None and len(molecules) == 1:
        molecules = molecules * len(coords)
    out = []
    for k, (m, c) in enumerate(zip(molecules, coords)):
        props = (properties or [{}] * len(molecules))[k]
        lines = [m.name or 'compound', '  rnafold3d screening', '',
                 f'{len(m):3d}{len(m.bonds):3d}  0  0  0  0  0  0  0  0999 V2000']
        for a in range(len(m)):
            x, y, z = c[a]
            lines.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {m.elements[a].capitalize():<3s} 0  0  0  0  0  0  0  0  0  0  0  0')
        for i, j, o in m.bonds:
            lines.append(f'{i + 1:3d}{j + 1:3d}{min(o, 4):3d}  0  0  0  0')
        charged = [(a, int(round(q))) for a, q in enumerate(m.charge) if abs(round(q)) >= 1]
        for a, q in charged:
            lines.append(f'M  CHG  1 {a + 1:3d} {q:3d}')
        lines.append('M  END')
        for key, value in (props or {}).items():
            lines += [f'> <{key}>', str(value), '']
        lines.append('$$$$')
        out.append('\n'.join(lines))
    text = '\n'.join(out) + '\n'
    if path is not None:
        from pathlib import Path
        Path(path).write_text(text, encoding='utf-8')
    return text
