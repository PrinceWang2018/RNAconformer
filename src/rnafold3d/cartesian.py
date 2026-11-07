"""Restrained Cartesian refinement of a full-atom model.

The torsion-space fold model builds stems as ideal A-form blocks and loops from ideal residue geometry,
so it reproduces an experimental template only to about 1.5-2.5 A. This step moves every atom freely:
target positions (template atoms superposed on the model) pull the model, while the covalent geometry
of the input model is kept by restraining

  - all distances within the sugar (C1' C2' C3' C4' O4' O2' O3' C5' + glycosidic N): ring, pucker, chirality;
  - all distances within the base plus C1' (planarity; the glycosidic angle chi stays free);
  - the phosphate group with the O3' of the previous residue, and the 1-3 distances P-C5', O5'-C4',
    C3'(i-1)-P(i) of the backbone (bond angles; the backbone torsions stay free),

and atom overlaps are pushed apart. All terms are harmonic in interatomic distances, minimised with L-BFGS.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree

from .build import BASE_ATOMS, GLY
from .checks import BACKBONE_BONDS, BASE_BONDS

SUGAR = ("C1'", "C2'", "C3'", "C4'", "O4'", "O2'", "O3'", "C5'")
PHOSPHATE = ('P', 'OP1', 'OP2', "O5'")


def _flatten(residues):
    keys, pts = [], []
    for r, at in enumerate(residues):
        for a, xyz in at.items():
            keys.append((r, a))
            pts.append(xyz)
    return keys, np.asarray(pts, float)


def _topology(seq, keys):
    """Restrained pairs (covalent geometry) and the 1-2/1-3/1-4 exclusion set for overlaps."""
    idx = {k: i for i, k in enumerate(keys)}
    pairs = set()

    def group(atoms):
        ids = [idx[a] for a in atoms if a in idx]
        for x in range(len(ids)):
            for y in range(x + 1, len(ids)):
                pairs.add((min(ids[x], ids[y]), max(ids[x], ids[y])))
    bonds = []
    for r, b in enumerate(seq):
        g = GLY[b]
        group([(r, a) for a in SUGAR + (g,)])
        group([(r, a) for a in BASE_ATOMS[b]] + [(r, "C1'")])
        group([(r, a) for a in PHOSPHATE] + ([(r - 1, "O3'")] if r else []))
        group([(r, 'P'), (r, "C5'")])
        group([(r, "O5'"), (r, "C4'")])
        group([(r, "O5'"), (r, "C5'")])
        if r:
            group([(r - 1, "C3'"), (r, 'P')])
        for x, y in BACKBONE_BONDS + BASE_BONDS[b] + [("C1'", g)]:
            if (r, x) in idx and (r, y) in idx:
                bonds.append((idx[(r, x)], idx[(r, y)]))
        if r + 1 < len(seq) and (r, "O3'") in idx and (r + 1, 'P') in idx:
            bonds.append((idx[(r, "O3'")], idx[(r + 1, 'P')]))
    n = len(keys)
    b = np.array(bonds, int).reshape(-1, 2)
    A = coo_matrix((np.ones(2 * len(b)), (np.r_[b[:, 0], b[:, 1]], np.r_[b[:, 1], b[:, 0]])), shape=(n, n)).tocsr()
    reach = A + A @ A
    reach = reach + reach @ A
    reach = reach.tocoo()
    excl = set(zip(np.minimum(reach.row, reach.col).tolist(), np.maximum(reach.row, reach.col).tolist()))
    p = np.array(sorted(pairs), int).reshape(-1, 2)
    return p, excl


def _ideal_links(seq, keys, p, d0, x0, geom):
    """Replace the restrained distances across the O3'(i-1)-P(i) linkage by ideal values (bond length and
    the C3'-O3'-P, O3'-P-O5', O3'-P-OP1/OP2 angles of the geometry file): closure strain is repaired."""
    import math
    b, a = geom.raw['bonds'], geom.raw['angles']

    def third(d1, d2, ang):
        return math.sqrt(d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * math.cos(math.radians(ang)))
    o3p = b["O3'-P"]['mean']
    want = {}
    idx = {k: i for i, k in enumerate(keys)}
    for r in range(1, len(seq)):
        o3, c3 = idx.get((r - 1, "O3'")), idx.get((r - 1, "C3'"))
        P = idx.get((r, 'P'))
        if o3 is None or P is None:
            continue
        want[(o3, P)] = o3p
        if c3 is not None:
            want[(c3, P)] = third(float(np.linalg.norm(x0[c3] - x0[o3])), o3p, a["C3'-O3'-P"]['mean'])
        for nm, bond, ang in (("O5'", "P-O5'", "O3'-P-O5'"), ('OP1', 'P-OP1', "O3'-P-OP1"), ('OP2', 'P-OP2', "O3'-P-OP2")):
            j = idx.get((r, nm))
            if j is not None:
                want[(o3, j)] = third(o3p, b[bond]['mean'], a[ang]['mean'])
    d0 = d0.copy()
    for q, (i, j) in enumerate(p):
        v = want.get((int(i), int(j))) or want.get((int(j), int(i)))
        if v is not None:
            d0[q] = v
    return d0


def _pair_terms(x, p, d0, k):
    diff = x[p[:, 0]] - x[p[:, 1]]
    d = np.sqrt(np.sum(diff * diff, axis=1)) + 1e-12
    dev = d - d0
    e = float(np.sum(k * dev * dev))
    f = (2.0 * k * dev / d)[:, None] * diff
    g = np.zeros_like(x)
    np.add.at(g, p[:, 0], f)
    np.add.at(g, p[:, 1], -f)
    return e, g


def refine_cartesian(seq, residues, targets, k_cov=60.0, k_pos=2.0, pos_tol=0.25, k_clash=20.0,
                     clash=2.9, polar_clash=2.6, rounds=(0.2, 0.5, 1.0, 1.0, 1.0), iters=200, geom=None,
                     reference=None):
    """Pull `residues` (list of {atom: xyz}) onto `targets` {(res, atom): xyz} keeping covalent geometry.

    rounds: target-weight multipliers for successive minimisations (gradual pull; the overlap pair list is
    rebuilt for every round). geom: use ideal O3'-P linkage geometry instead of the input's (repairs closure
    strain). reference: model whose covalent geometry is kept (default: the input; same atoms).
    Returns (new residues, stats)."""
    keys, x0 = _flatten(residues)
    n = len(keys)
    p, excl = _topology(seq, keys)
    xr = x0 if reference is None else np.array([reference[r][a] for r, a in keys], float)
    d0 = np.linalg.norm(xr[p[:, 0]] - xr[p[:, 1]], axis=1)
    if geom is not None:
        d0 = _ideal_links(seq, keys, p, d0, xr, geom)
    kc = np.full(len(p), k_cov)
    idx = {k: i for i, k in enumerate(keys)}
    tk = [(idx[k], v) for k, v in targets.items() if k in idx and v is not None and np.all(np.isfinite(v))]
    if not tk:
        return residues, None
    ti = np.array([a for a, _ in tk], int)
    tx = np.array([v for _, v in tk], float)
    elem = np.array([('O' if a.startswith('OP') else a[0]) for _, a in keys])
    polar = np.isin(elem, ['N', 'O'])
    res_of = np.array([r for r, _ in keys])

    def clash_pairs(x):
        q = cKDTree(x).query_pairs(max(clash, polar_clash) + 1.5, output_type='ndarray')
        if not len(q):
            return np.zeros((0, 2), int), np.zeros(0)
        keep = np.array([(int(a), int(b)) not in excl for a, b in q], bool)
        q = q[keep]
        lim = np.where(polar[q[:, 0]] & polar[q[:, 1]], polar_clash, clash)
        return q, lim

    def rms_to_targets(x):
        return float(np.sqrt(np.mean(np.sum((x[ti] - tx) ** 2, axis=1))))

    x = x0.copy()
    start = rms_to_targets(x)
    for scale in rounds:
        cp, lim = clash_pairs(x)

        def fun(v):
            X = v.reshape(n, 3)
            e, g = _pair_terms(X, p, d0, kc)
            dv = X[ti] - tx
            dd = np.sqrt(np.sum(dv * dv, axis=1)) + 1e-12
            over = np.maximum(0.0, dd - pos_tol)
            e += float(np.sum(k_pos * scale * over * over))
            np.add.at(g, ti, (2.0 * k_pos * scale * over / dd)[:, None] * dv)
            if len(cp):
                diff = X[cp[:, 0]] - X[cp[:, 1]]
                d = np.sqrt(np.sum(diff * diff, axis=1)) + 1e-12
                ov = np.maximum(0.0, lim - d)
                e += float(np.sum(k_clash * ov * ov))
                f = (-2.0 * k_clash * ov / d)[:, None] * diff
                np.add.at(g, cp[:, 0], f)
                np.add.at(g, cp[:, 1], -f)
            return e, g.ravel()
        r = minimize(fun, x.ravel(), jac=True, method='L-BFGS-B', options={'maxiter': iters, 'gtol': 1e-4})
        x = r.x.reshape(n, 3)
    out = [dict() for _ in residues]
    for (r, a), v in zip(keys, x):
        out[r][a] = v
    cov_dev = np.abs(np.linalg.norm(x[p[:, 0]] - x[p[:, 1]], axis=1) - d0)
    stats = {'targets': len(ti), 'rms_to_targets_start_A': round(start, 3), 'rms_to_targets_A': round(rms_to_targets(x), 3),
             'max_covalent_deviation_A': round(float(cov_dev.max()), 3),
             'atoms_moved_rms_A': round(float(np.sqrt(np.mean(np.sum((x - x0) ** 2, axis=1)))), 3)}
    return out, stats


def template_targets(seq, residues, template_atoms, atom_names, anchor="C4'", group=None):
    """Superpose template atoms {res: [[x, y, z] | None per atom name]} on the model (on the anchor atom of
    every mapped residue) and return {(res, atom): xyz} targets. 'NG' is the glycosidic N of the query base.
    group: {res: g} when the atoms come from several templates; each is superposed on its own."""
    k = atom_names.index(anchor)
    groups = {}
    for r in template_atoms:
        groups.setdefault((group or {}).get(r, 0), []).append(r)
    out = {}
    for members in groups.values():
        q = [r for r in members if template_atoms[r][k] is not None and anchor in residues[r]]
        if len(q) < 3:
            continue
        A = np.array([template_atoms[r][k] for r in q], float)
        B = np.array([residues[r][anchor] for r in q], float)
        a0, b0 = A.mean(0), B.mean(0)
        u, _, vt = np.linalg.svd((A - a0).T @ (B - b0))
        d = np.sign(np.linalg.det(u @ vt))
        R = u @ np.diag([1.0, 1.0, d]) @ vt
        for r in members:
            for nm, v in zip(atom_names, template_atoms[r]):
                if v is not None:
                    out[(r, GLY[seq[r]] if nm == 'NG' else nm)] = (np.asarray(v, float) - a0) @ R + b0
    return out


def regularize(seq, residues, geom, k_pos=1.0, pos_tol=0.3, reference=None, k_cov=200.0):
    """Repair local strain (linkage geometry, overlaps) of a model while keeping its shape: every atom is
    restrained loosely to its own position; covalent geometry from `reference` (default: the model)."""
    targets = {(r, a): np.asarray(v, float) for r, at in enumerate(residues) for a, v in at.items()}
    return refine_cartesian(seq, residues, targets, k_cov=k_cov, k_pos=k_pos, pos_tol=pos_tol,
                            rounds=(1.0, 1.0, 1.0), geom=geom, reference=reference)
