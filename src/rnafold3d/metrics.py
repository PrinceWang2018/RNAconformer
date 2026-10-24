"""Structure-comparison metrics for RNA models against a reference.

* C4' RMSD after optimal superposition over residues present in both (matched by number)
* eRMSD (Bottaro, Di Palma & Bussi, NAR 2014): base-centred G-vectors from C2/C4/C6 frames,
  ellipsoidal scaling (5, 5, 3 A), cutoff 2.4; eRMSD < ~1 indicates the same base arrangement
* INF-WC (Parisien et al., RNA 2009): sqrt(PPV * sensitivity) over canonical base pairs
  annotated with rnafold3d.annotate (canonical cWW + GU only; not INF-all)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .annotate import canonical_pairs
from .structure_io import Residue, FALLBACK_PARENTS, CANONICAL


def read_numbered_pdb(path, chain=None):
    """{residue number: Residue} using the PDB residue number as the sequence index."""
    out = {}
    for line in Path(path).read_text(errors='replace').splitlines():
        if line.startswith('ENDMDL'):
            break
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        if line[16].strip() not in ('', 'A', '1'):
            continue
        if chain is not None and line[21].strip() != chain:
            continue
        comp = line[17:20].strip().upper()
        base = comp if comp in CANONICAL else FALLBACK_PARENTS.get(comp)
        if base is None:
            continue
        name = line[12:16].strip().replace('*', "'")
        name = {'O1P': 'OP1', 'O2P': 'OP2'}.get(name, name)
        if name.startswith('H'):
            continue
        num = int(line[22:26])
        res = out.setdefault(num, Residue(num, base, comp, {}))
        res.atoms.setdefault(name, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    return out


def kabsch_rmsd(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((a @ r - b) ** 2, axis=1))))


def c4_rmsd(ref, pred):
    common = sorted(k for k in ref if k in pred and "C4'" in ref[k].atoms and "C4'" in pred[k].atoms
                    and ref[k].base == pred[k].base)
    if len(common) < 3:
        return None, len(common)
    return kabsch_rmsd([ref[k].atoms["C4'"] for k in common], [pred[k].atoms["C4'"] for k in common]), len(common)


def _gmat(coords, cutoff=2.4):
    """coords: (N, 3, 3) array of C2, C4, C6 per residue."""
    origo = coords.mean(axis=1)
    x = coords[:, 0] - origo
    x /= np.linalg.norm(x, axis=1)[:, None]
    c = coords[:, 1] - coords[:, 2]
    z = np.cross(x, c)
    z /= np.linalg.norm(z, axis=1)[:, None]
    y = np.cross(z, x)
    lcs = np.stack([x, y, z], axis=2)                   # columns are the axes
    dvec = origo[None, :, :] - origo[:, None, :]        # [i, j] = o_j - o_i
    proj = np.einsum('ijk,ikl->ijl', dvec, lcs) * np.array([1 / 5.0, 1 / 5.0, 1 / 3.0])
    dd = np.linalg.norm(proj, axis=2)
    gamma = np.pi / cutoff
    g = np.zeros(proj.shape[:2] + (4,))
    m = (dd < cutoff) & (dd > 0)
    g[m, :3] = (np.sin(gamma * dd[m]) / dd[m])[:, None] * proj[m]
    g[m, 3] = 1.0 + np.cos(gamma * dd[m])
    return g / gamma


def ermsd(ref, pred, cutoff=2.4):
    keys = sorted(k for k in ref if k in pred and all(a in ref[k].atoms and a in pred[k].atoms for a in ('C2', 'C4', 'C6')))
    if len(keys) < 2:
        return None, len(keys)
    a = np.array([[ref[k].atoms[n] for n in ('C2', 'C4', 'C6')] for k in keys])
    b = np.array([[pred[k].atoms[n] for n in ('C2', 'C4', 'C6')] for k in keys])
    ga, gb = _gmat(a, cutoff), _gmat(b, cutoff)
    return float(np.sqrt(np.sum((ga - gb) ** 2) / len(keys))), len(keys)


def inf_wc(ref, pred):
    """INF over canonical pairs; residues absent from the reference are ignored in the model."""
    rp = {(i, j) for i, j, _ in canonical_pairs(list(ref.values()))}
    pp = {(i, j) for i, j, _ in canonical_pairs([pred[k] for k in pred if k in ref])}
    tp = len(rp & pp)
    if not rp and not pp:
        return 1.0, {'tp': 0, 'fp': 0, 'fn': 0}
    ppv = tp / len(pp) if pp else 0.0
    sty = tp / len(rp) if rp else 0.0
    return float(np.sqrt(ppv * sty)), {'tp': tp, 'fp': len(pp - rp), 'fn': len(rp - pp)}


def compare(reference_pdb, model_pdb):
    ref = read_numbered_pdb(reference_pdb)
    pred = read_numbered_pdb(model_pdb)
    r, n = c4_rmsd(ref, pred)
    e, ne = ermsd(ref, pred)
    inf, counts = inf_wc(ref, pred)
    return {'c4_rmsd': None if r is None else round(r, 3), 'n_c4': n, 'ermsd': None if e is None else round(e, 3),
            'inf_wc': round(inf, 3), 'pairs': counts, 'reference_residues': len(ref)}
