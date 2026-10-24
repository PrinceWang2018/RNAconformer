"""Independent all-atom validation of RNA models (v2).

Independent of the builder: works from atom names and coordinates only, with
reference values from the calibration file. Applicable to predicted models and to
experimental structures (used as a sanity check of the thresholds).

Hard gate (``passed``):
  * every residue has its full heavy-atom set (5'-terminal phosphate optional)
  * all bond lengths within max(6 sd, 0.08 A) of calibrated means
  * all bond angles within max(6 sd, 8 deg) of calibrated means
  * O3'(i)-P(i+1) present and 1.50-1.72 A for every consecutive pair (no chain breaks)
  * correct D-ribose chirality at C1', C2', C3', C4'
  * no severe clashes (non-bonded heavy atoms, >3 bonds apart, closer than 2.2 A)
Diagnostics: 4-sd outlier counts, heavy-atom clash proxy, pucker distribution, input-pair recovery.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .build import BASE_ATOMS, GLY, load_geometry, atom_order
from .annotate import canonical_pairs
from .geomutil import angle as _angle, pseudorotation

BASE_BONDS = {
    'A': [('N9', 'C8'), ('C8', 'N7'), ('N7', 'C5'), ('C5', 'C6'), ('C6', 'N6'), ('C6', 'N1'), ('N1', 'C2'),
          ('C2', 'N3'), ('N3', 'C4'), ('C4', 'C5'), ('C4', 'N9')],
    'G': [('N9', 'C8'), ('C8', 'N7'), ('N7', 'C5'), ('C5', 'C6'), ('C6', 'O6'), ('C6', 'N1'), ('N1', 'C2'),
          ('C2', 'N2'), ('C2', 'N3'), ('N3', 'C4'), ('C4', 'C5'), ('C4', 'N9')],
    'C': [('N1', 'C2'), ('C2', 'O2'), ('C2', 'N3'), ('N3', 'C4'), ('C4', 'N4'), ('C4', 'C5'), ('C5', 'C6'), ('C6', 'N1')],
    'U': [('N1', 'C2'), ('C2', 'O2'), ('C2', 'N3'), ('N3', 'C4'), ('C4', 'O4'), ('C4', 'C5'), ('C5', 'C6'), ('C6', 'N1')],
}
BACKBONE_BONDS = [('P', 'OP1'), ('P', 'OP2'), ('P', "O5'"), ("O5'", "C5'"), ("C5'", "C4'"), ("C4'", "O4'"),
                  ("O4'", "C1'"), ("C1'", "C2'"), ("C2'", "C3'"), ("C3'", "C4'"), ("C3'", "O3'"), ("C2'", "O2'")]
ANGLES = ["P-O5'-C5'", "O5'-C5'-C4'", "C5'-C4'-C3'", "C5'-C4'-O4'", "C4'-C3'-O3'", "C2'-C3'-O3'", "C4'-O4'-C1'",
          "O4'-C1'-C2'", "C1'-C2'-C3'", "C2'-C3'-C4'", "C3'-C4'-O4'", "C1'-C2'-O2'", "C3'-C2'-O2'",
          "OP1-P-OP2", "OP1-P-O5'", "OP2-P-O5'"]
CHIRAL = {"C1'": ("O4'", "C2'", 'N'), "C2'": ("C1'", "C3'", "O2'"), "C3'": ("C2'", "C4'", "O3'"),
          "C4'": ("C3'", "O4'", "C5'")}
VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80}


def _bond_graph(sequence, residues):
    """Adjacency over (res, atom) keys for present atoms."""
    adj = {}

    def add(a, b):
        if a[1] in residues[a[0]] and b[1] in residues[b[0]]:
            adj.setdefault(a, set()).add(b)
            adj.setdefault(b, set()).add(a)
    for r, base in enumerate(sequence):
        for x, y in BACKBONE_BONDS + BASE_BONDS[base]:
            add((r, x), (r, y))
        add((r, "C1'"), (r, GLY[base]))
        if r + 1 < len(sequence):
            add((r, "O3'"), (r + 1, 'P'))
    return adj


def _within_bonds(adj, a, b, depth=3):
    frontier = {a}
    seen = {a}
    for _ in range(depth):
        nxt = set()
        for u in frontier:
            for v in adj.get(u, ()):
                if v == b:
                    return True
                if v not in seen:
                    seen.add(v)
                    nxt.add(v)
        frontier = nxt
    return False


def _expected_chirality(geom):
    s = geom.sugar('C3endo')
    out = {}
    for c, (a, b, d) in CHIRAL.items():
        out[c] = float(np.sign(np.dot(s[a] - s[c], np.cross(s[b] - s[c], s[d] - s[c]))))
    return out


def validate_model(sequence, residues, pairs=None, geom=None, max_report=60):
    """residues: list of {atom: xyz} (0-based list, one per nucleotide); pairs: 0-based input pairs."""
    g = geom or load_geometry()
    bonds_ref = g.raw['bonds']
    angles_ref = g.raw['angles']
    issues = []
    stats = {'bond_outliers_4sd': 0, 'angle_outliers_4sd': 0, 'n_bonds': 0, 'n_angles': 0}

    def issue(kind, r, detail, **kw):
        issues.append({'kind': kind, 'residue': r + 1, 'detail': detail, **kw})

    n = len(sequence)
    if len(residues) != n:
        raise ValueError('residue count differs from sequence length')
    # completeness
    for r, base in enumerate(sequence):
        need = [a for a in atom_order(base) if not (r == 0 and a in ('P', 'OP1', 'OP2'))]
        missing = [a for a in need if a not in residues[r]]
        if missing:
            issue('missing_atoms', r, ','.join(missing))
        for a, xyz in residues[r].items():
            if not np.all(np.isfinite(xyz)):
                issue('nonfinite', r, a)
    base_ref = {b: {k: np.array(v) for k, v in g.raw['base_templates'][b]['atoms'].items()} for b in 'ACGU'}
    for r, base in enumerate(sequence):
        at = residues[r]
        for x, y in BACKBONE_BONDS:
            if x in at and y in at:
                ref = bonds_ref[f'{x}-{y}']
                d = float(np.linalg.norm(at[x] - at[y]))
                stats['n_bonds'] += 1
                dev = abs(d - ref['mean'])
                if dev > 4 * ref['sd']:
                    stats['bond_outliers_4sd'] += 1
                if dev > max(6 * ref['sd'], 0.08):
                    issue('bond_length', r, f'{x}-{y}', observed=round(d, 3), expected=ref['mean'])
        for x, y in BASE_BONDS[base] + [("C1'", GLY[base])]:
            if x in at and y in at:
                ref_d = float(np.linalg.norm(base_ref[base][x] - base_ref[base][y])) if x != "C1'" else bonds_ref["C1'-N"]['mean']
                d = float(np.linalg.norm(at[x] - at[y]))
                stats['n_bonds'] += 1
                if abs(d - ref_d) > 0.08:
                    issue('bond_length', r, f'{x}-{y}', observed=round(d, 3), expected=round(ref_d, 3))
        for name in ANGLES:
            x, y, z = name.split('-')
            if x in at and y in at and z in at:
                ref = angles_ref[name]
                v = _angle(at[x], at[y], at[z])
                stats['n_angles'] += 1
                dev = abs(v - ref['mean'])
                if dev > 4 * ref['sd']:
                    stats['angle_outliers_4sd'] += 1
                if dev > max(6 * ref['sd'], 8.0):
                    issue('bond_angle', r, name, observed=round(v, 1), expected=ref['mean'])
    # linkage
    lk = bonds_ref["O3'-P"]
    for r in range(n - 1):
        a, b = residues[r], residues[r + 1]
        if "O3'" not in a or 'P' not in b:
            issue('chain_break', r, "missing O3'-P link")
            continue
        d = float(np.linalg.norm(a["O3'"] - b['P']))
        stats['n_bonds'] += 1
        if abs(d - lk['mean']) > 4 * lk['sd']:
            stats['bond_outliers_4sd'] += 1
        if not 1.50 <= d <= 1.72:
            issue('chain_break', r, "O3'-P(next)", observed=round(d, 3))
        for name, (x, y, z) in (("C3'-O3'-P", (a.get("C3'"), a.get("O3'"), b.get('P'))),
                                ("O3'-P-O5'", (a.get("O3'"), b.get('P'), b.get("O5'")))):
            if x is None or y is None or z is None:
                continue
            ref = angles_ref[name]
            v = _angle(x, y, z)
            stats['n_angles'] += 1
            if abs(v - ref['mean']) > 4 * ref['sd']:
                stats['angle_outliers_4sd'] += 1
            if abs(v - ref['mean']) > max(6 * ref['sd'], 8.0):
                issue('bond_angle', r, name, observed=round(v, 1), expected=ref['mean'])
    # chirality
    exp = _expected_chirality(g)
    for r, base in enumerate(sequence):
        at = dict(residues[r])
        if GLY[base] in at:
            at['N'] = at[GLY[base]]
        for c, (a1, b1, d1) in CHIRAL.items():
            if all(k in at for k in (c, a1, b1, d1)):
                v = float(np.dot(at[a1] - at[c], np.cross(at[b1] - at[c], at[d1] - at[c])))
                if np.sign(v) != exp[c] or abs(v) < 0.3:
                    issue('chirality', r, c, signed_volume=round(v, 3))
    # clashes
    keys, pts, elem = [], [], []
    for r, at in enumerate(residues):
        for a, xyz in at.items():
            keys.append((r, a))
            pts.append(xyz)
            elem.append('O' if a.startswith('OP') else a[0])
    pts = np.array(pts)
    adj = _bond_graph(sequence, residues)
    severe = 0
    overlaps = 0
    min_nb = np.inf
    if len(pts) > 1:
        for i, j in cKDTree(pts).query_pairs(3.4, output_type='ndarray'):
            ki, kj = keys[i], keys[j]
            if abs(ki[0] - kj[0]) <= 1 and _within_bonds(adj, ki, kj, 3):
                continue
            if ki[0] == kj[0] and _within_bonds(adj, ki, kj, 3):
                continue
            d = float(np.linalg.norm(pts[i] - pts[j]))
            min_nb = min(min_nb, d)
            polar = elem[i] in 'NO' and elem[j] in 'NO'
            if d < 2.2:
                severe += 1
                if severe <= 25:
                    issue('severe_clash', ki[0], f'{ki[1]} vs {kj[0] + 1}:{kj[1]}', distance=round(d, 2))
            lim = 2.5 if polar else VDW.get(elem[i], 1.7) + VDW.get(elem[j], 1.7) - 0.4 - 0.5
            if d < lim:
                overlaps += 1
    # sugar puckers
    puck = {'C3endo': 0, 'C2endo': 0, 'other': 0}
    for at in residues:
        if all(k in at for k in ("C4'", "O4'", "C1'", "C2'", "C3'")):
            ph, _ = pseudorotation([at[k] for k in ("C4'", "O4'", "C1'", "C2'", "C3'")])
            if ph < 50 or ph > 340:
                puck['C3endo'] += 1
            elif 130 <= ph <= 200:
                puck['C2endo'] += 1
            else:
                puck['other'] += 1
    # input-pair recovery
    pair_report = None
    if pairs is not None:
        class _R:
            def __init__(self, i, b, a):
                self.index, self.base, self.atoms = i, b, a
        found = {(i, j) for i, j, _ in canonical_pairs([_R(k, sequence[k], residues[k]) for k in range(n)])}
        want = {(min(i, j), max(i, j)) for i, j in pairs}
        canon = {p for p in want if ''.join(sorted(sequence[p[0]] + sequence[p[1]])) in ('AU', 'CG', 'GU')}
        rec = len(want & found)
        pair_report = {'input_pairs': len(want), 'canonical_input_pairs': len(canon), 'recovered': rec,
                       'recovery': round(rec / len(canon), 3) if canon else None,
                       'extra_canonical_pairs': len(found - want)}
    hard = [i for i in issues]
    kinds = {}
    for i in hard:
        kinds[i['kind']] = kinds.get(i['kind'], 0) + 1
    passed = not hard
    n_atoms = len(pts)
    return {
        'passed': passed,
        'issue_counts': kinds,
        'issues': hard[:max_report],
        'n_atoms': n_atoms,
        'severe_clashes': severe,
        'heavy_atom_overlaps_per_1000_atoms': round(1000.0 * overlaps / max(1, n_atoms), 2),
        'min_nonbonded_distance': None if not np.isfinite(min_nb) else round(min_nb, 3),
        **stats,
        'sugar_puckers': puck,
        'pair_recovery': pair_report,
        'criteria': {'bond': 'max(6 sd, 0.08 A)', 'angle': 'max(6 sd, 8 deg)', 'link': "O3'-P 1.50-1.72 A",
                     'severe_clash': '<2.2 A, >3 bonds apart', 'chirality': 'D-ribose C1-C4 signs'},
        'note': 'Heavy-atom checks only (no hydrogens); overlap proxy is not a MolProbity clashscore.',
    }
