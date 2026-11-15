"""Knowledge-based base-base potential for ranking models (statistics over the structure library).

For every two nucleotides that are not in the same helix and at least MIN_SEP apart in sequence, the distance
between their base centroids and how far each base lies above the other's plane (0 = side by side, as in a pair
or a base triple; 1 = on top, as in a stack) are binned, separately for the ten nucleotide-pair types and two
sequence-separation classes. The potential is -ln(P_type(d, h) / P_all(d, h)) (a reference state of all types
pooled, as in RASP and 3dRNAscore), counted over the library chains released before 2024-05-01 that are not
homologous to the development benchmark (tools/build_base_potential.py, data/basepot_v1.npz).

Pairs inside one helix are left out both in the statistics and in the score: the helices are given by the 2D,
which the ranking already scores through its free energy, so the potential measures only the tertiary
arrangement: how helices, loops and junction strands pack against each other. That is the part the hand-set
terms of the search (compaction, coaxial stacking, phosphate repulsion) only approximate.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / 'data' / 'basepot_v1.npz'
MIN_SEP = 3
D_MIN, D_MAX, D_STEP = 2.0, 16.0, 0.5
H_BINS = 4
H_CUT = 10.0                       # orientation resolved below this distance only
N_D = int(round((D_MAX - D_MIN) / D_STEP))
TYPES = {}
_k = 0
for _a, _x in enumerate('ACGU'):
    for _y in 'ACGU'[_a:]:
        TYPES[(_x, _y)] = TYPES[(_y, _x)] = _k
        _k += 1
N_TYPES = _k                        # 10


def base_frames(seq, atoms):
    """Base centroids and unit normals from the six-membered ring atoms C2, C4, C6 (NaN where missing).

    atoms: {'C2': (n, 3), 'C4': (n, 3), 'C6': (n, 3)}."""
    c2, c4, c6 = (np.asarray(atoms[k], float) for k in ('C2', 'C4', 'C6'))
    cen = (c2 + c4 + c6) / 3.0
    nrm = np.cross(c4 - c2, c6 - c2)
    ln = np.linalg.norm(nrm, axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        nrm = nrm / ln[:, None]
    return cen, nrm


def helix_id(n, pairs):
    """Helix index per residue (-1 for unpaired) from nested canonical pairs (stacked runs)."""
    from .fold import stems_from_pairs
    hid = np.full(n, -1)
    for k, st in enumerate(stems_from_pairs(sorted((min(i, j), max(i, j)) for i, j in pairs))):
        for i, j in st:
            hid[i] = hid[j] = k
    return hid


def pair_features(seq, cen, nrm, hid, d_max=D_MAX):
    """Index arrays (type, separation class, distance bin, orientation bin) of the counted pairs."""
    from scipy.spatial import cKDTree
    ok = np.where(~np.isnan(cen[:, 0]) & ~np.isnan(nrm[:, 0]))[0]
    if len(ok) < 2:
        return (np.zeros(0, int),) * 4
    pr = cKDTree(cen[ok]).query_pairs(d_max - 1e-6, output_type='ndarray')
    if not len(pr):
        return (np.zeros(0, int),) * 4
    i, j = ok[pr[:, 0]], ok[pr[:, 1]]
    keep = (np.abs(i - j) >= MIN_SEP) & ~((hid[i] == hid[j]) & (hid[i] >= 0))
    i, j = i[keep], j[keep]
    v = cen[j] - cen[i]
    d = np.linalg.norm(v, axis=1)
    keep = d >= D_MIN
    i, j, v, d = i[keep], j[keep], v[keep], d[keep]
    u = v / d[:, None]
    h = 0.5 * (np.abs(np.sum(nrm[i] * u, axis=1)) + np.abs(np.sum(nrm[j] * u, axis=1)))
    hb = np.minimum((h * H_BINS).astype(int), H_BINS - 1)
    hb = np.where(d < H_CUT, hb, 0)
    db = np.minimum(((d - D_MIN) / D_STEP).astype(int), N_D - 1)
    codes = np.array([TYPES.get((a, b), -1) for a, b in zip((seq[x] for x in i), (seq[x] for x in j))], int)
    sc = np.where(np.abs(i - j) < 6, 0, 1)
    keep = codes >= 0
    return codes[keep], sc[keep], db[keep], hb[keep]


def potential_from_counts(counts, pseudo=50.0, cap=3.0):
    """counts[type, sep, d, h] -> energy table of the same shape (kT units, capped).

    P_ref(d, h | sep): all types pooled; P_type(d, h | sep) with `pseudo` pseudo-observations drawn from P_ref, so
    that sparse bins fall back to zero energy instead of to noise."""
    c = counts.astype(float)
    ref = c.sum(axis=0, keepdims=True)
    p_r = ref / ref.sum(axis=(2, 3), keepdims=True).clip(1)
    n_t = c.sum(axis=(2, 3), keepdims=True)
    p_t = (c + pseudo * p_r) / (n_t + pseudo)
    with np.errstate(divide='ignore', invalid='ignore'):
        e = -np.log(p_t / p_r)
    e[~np.isfinite(e)] = 0.0
    return np.clip(e, -cap, cap)


@lru_cache(maxsize=1)
def table(path=None):
    p = Path(path) if path else DATA
    if not p.exists():
        return None
    z = np.load(p)
    return z['energy']


def score(seq, residues, pairs, tab=None):
    """Sum of the potential over the counted pairs of an all-atom model (residues: list of {atom: xyz})."""
    tab = table() if tab is None else tab
    if tab is None:
        return 0.0
    n = len(seq)
    atoms = {k: np.array([r.get(k, (np.nan,) * 3) for r in residues], float) for k in ('C2', 'C4', 'C6')}
    cen, nrm = base_frames(seq, atoms)
    t, s, d, h = pair_features(seq, cen, nrm, helix_id(n, pairs))
    return float(tab[t, s, d, h].sum()) if len(t) else 0.0
