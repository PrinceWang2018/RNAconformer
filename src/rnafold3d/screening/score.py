"""The scoring function: how well a compound in a given pose fits an RNA.

Eight terms. The first five are the shape-and-contact terms every docking score has, in the form AutoDock
Vina uses (its scoring function is Apache-2.0 and may be used as a reference; no GPL docking code or
parameters are copied here). The last three are the ones an RNA score cannot do without, and they are why
this file exists rather than a call to a protein docking program:

    gauss1, gauss2   contact: two Gaussians of the surface distance, one tight, one loose
    repulsion        overlap: the square of the surface distance where it is negative
    hydrophobic      apolar atom against apolar atom at contact
    hbond            donor against acceptor at contact (receptor geometry only, see `stack` for why)
    stacking         an aromatic compound ring against the face of a base, with the angle between the two
                     planes: the dominant interaction of RNA binders and the term protein scores get worst
    electrostatic    Coulomb with a distance-dependent dielectric (ε = 4r), explicit phosphate charges
    desolvation      the cost of burying a charged compound atom against the RNA

and a final division by (1 + 0.0585 · rotatable bonds), the entropy a compound pays for its own flexibility.

All terms are functions of the *surface* distance d = r − R_i − R_j (so atom size is already in them) and
vanish beyond 8 Å, which is what lets them be put on a grid (`grid.py`) and evaluated as a lookup.

**The weights of the three RNA terms are not fitted.** They are set to put a good stack, a salt bridge and
a buried charge on the scale of a hydrogen bond, and the plan (`docs/NEXT_STEPS.md` §1.4) says they may only
ever be tuned against the enrichment benchmark, never against the selectivity cases. Until that benchmark
exists, what this function produces is a plausible ranking and not a measurement. It is not an affinity and
is not reported in kcal/mol anywhere in the interface.
"""
from __future__ import annotations

import numpy as np

CUTOFF = 8.0                        # Å, beyond which every term is zero
EPS_MIN = 3.0                       # Å, the closest two charges are allowed to count: an ion pair in water
                                    # does not come closer than its two hydration shells allow
COULOMB = 332.0                     # kcal·Å/(mol·e²)

WEIGHTS = {
    'gauss1': -0.0356, 'gauss2': -0.00516, 'repulsion': 0.840,      # Vina's shape terms
    'hydrophobic': -0.0351, 'hbond': -0.587,
    'stacking': -0.250,             # per aromatic compound atom: a face-to-face ring pair is about −1.4
    'electrostatic': 0.050,         # scales a screened Coulomb energy into the same units
    'desolvation': 0.110,           # the cost of burying one unit charge against the RNA
    'metal': -0.700,                # an acceptor or anion at coordinating distance of Mg²⁺ or Mn²⁺
    'rotatable': 0.0585,            # the entropy term, applied as a division at the end
}
TERMS = ('gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hbond', 'stacking', 'electrostatic',
         'desolvation', 'metal')


def _slope(x, lo, hi):
    """1 below `lo`, 0 above `hi`, linear in between (the ramp Vina uses for contact terms)."""
    return np.clip((hi - x) / (hi - lo), 0.0, 1.0)


def pair_terms(lig_xyz, lig, rec, cutoff=CUTOFF):
    """Every term, summed over the compound-receptor atom pairs, for one pose.

    Returns a dict of raw term values (before the weights), which is what `score` weights and what the grid
    in `grid.py` has to reproduce atom by atom.
    """
    from scipy.spatial import cKDTree
    lig_xyz = np.asarray(lig_xyz, float)
    out = {t: 0.0 for t in TERMS}
    if not len(lig_xyz) or not len(rec):
        return out
    pairs = cKDTree(lig_xyz).query_ball_tree(cKDTree(rec.xyz), cutoff)
    i = np.concatenate([np.full(len(p), a, int) for a, p in enumerate(pairs)]) if any(pairs) else np.zeros(0, int)
    j = np.concatenate([np.asarray(p, int) for p in pairs]) if any(pairs) else np.zeros(0, int)
    if not len(i):
        return out
    r = np.linalg.norm(lig_xyz[i] - rec.xyz[j], axis=1)
    d = r - lig.radius[i] - rec.radius[j]
    w = _slope(r, cutoff - 2.0, cutoff)              # the same smooth cut the grid applies, so both agree
    out['gauss1'] = float(np.sum(np.exp(-(d / 0.5) ** 2) * w))
    out['gauss2'] = float(np.sum(np.exp(-((d - 3.0) / 2.0) ** 2) * w))
    out['repulsion'] = float(np.sum(np.where(d < 0, d * d, 0.0)))
    phob = lig.hydrophobic[i] & rec.hydrophobic[j]
    out['hydrophobic'] = float(np.sum(_slope(d[phob], 0.5, 1.5)))
    # a pair counts once per direction: a compound hydroxyl facing a 2′-hydroxyl can bond either way
    hb = (lig.donor[i] & rec.acceptor[j]).astype(float) + (lig.acceptor[i] & rec.donor[j]).astype(float)
    out['hbond'] = float(np.sum(_slope(d, -0.7, 0.0) * hb))
    q = lig.charge[i] * rec.charge[j]
    out['electrostatic'] = float(np.sum(COULOMB * q / (4.0 * np.maximum(r, EPS_MIN) ** 2) * w))
    vol = 4 / 3 * np.pi * rec.radius[j] ** 3
    out['desolvation'] = float(np.sum(np.abs(lig.charge[i]) * vol * np.exp(-(r / 3.5) ** 2) / 100.0))
    if rec.metal.any():
        met = (rec.metal[j] & (lig.acceptor[i] | (lig.charge[i] < -0.3)))
        out['metal'] = float(np.sum(np.exp(-((r[met] - 2.1) / 0.5) ** 2) * w[met]))
    out['stacking'] = stacking(lig_xyz, lig, rec)
    return out


def stacking(lig_xyz, lig, rec):
    """Aromatic rings of the compound against the faces of the bases.

    Per aromatic compound atom and base plane: a Gaussian on the distance from the base plane (3.6 Å, the
    stacking distance), a Gaussian on how far the atom sits off the ring axis, and the square of the cosine
    between the two ring normals, so face-to-face counts and edge-on does not. Summed over the aromatic
    atoms of the compound, a well-stacked six-membered ring reaches about 5 to 6 units.

    This is the term that has to be in the geometry rather than in a pairwise distance table: the base is a
    plane, not a ball, and a score that cannot tell a stack from a side-on contact ranks RNA binders badly.
    """
    normals = ring_normals(lig_xyz, lig)
    if normals is None or not len(rec.ring_centre):
        return 0.0
    idx = np.nonzero(lig.aromatic)[0]
    if not len(idx):
        return 0.0
    v = lig_xyz[idx][:, None, :] - rec.ring_centre[None, :, :]           # (atoms, bases, 3)
    perp = np.abs(np.einsum('abk,bk->ab', v, rec.ring_normal))
    lateral = np.sqrt(np.maximum(0.0, np.sum(v * v, axis=2) - perp ** 2))
    cos2 = np.einsum('ak,bk->ab', normals[idx], rec.ring_normal) ** 2
    e = np.exp(-((perp - 3.6) / 0.9) ** 2) * np.exp(-(lateral / 2.2) ** 2) * cos2
    return float(e.sum())


def ring_normals(lig_xyz, lig):
    """The normal of its aromatic ring for every aromatic atom of the compound (zeros elsewhere)."""
    out = np.zeros((len(lig_xyz), 3))
    xyz = np.asarray(lig_xyz, float)
    for ring in lig.rings:
        ring = sorted(set(ring))
        if not lig.aromatic[ring].all():
            continue
        p = xyz[ring]
        n = np.linalg.svd(p - p.mean(0))[2][2]
        out[ring] = n / (np.linalg.norm(n) or 1.0)
    return out


def score(lig_xyz, lig, rec, weights=None, cutoff=CUTOFF, breakdown=False):
    """The score of one pose: lower is better. Not an affinity, and not in kcal/mol — see the module text.

    With `breakdown`, also the weighted contribution of each term, which is what a report shows and what
    makes a wrong ranking debuggable."""
    w = dict(WEIGHTS, **(weights or {}))
    raw = pair_terms(lig_xyz, lig, rec, cutoff)
    parts = {t: w[t] * raw[t] for t in TERMS}
    total = sum(parts.values()) / (1.0 + w['rotatable'] * lig.n_rotatable)
    if breakdown:
        return total, {'raw': raw, 'weighted': parts, 'rotatable': lig.n_rotatable}
    return total


def strain(lig_xyz, lig, k=10.0):
    """How much the pose distorts the compound itself: the change of its bond lengths from the input, as a
    penalty. The search only turns rotatable bonds, so this stays near zero; it exists to catch a pose built
    from a bad conformer, and to keep the ranking honest when ring conformers are enumerated later."""
    xyz = np.asarray(lig_xyz, float)
    if not len(lig.bonds):
        return 0.0
    i = np.array([b[0] for b in lig.bonds])
    j = np.array([b[1] for b in lig.bonds])
    d0 = np.linalg.norm(lig.xyz[i] - lig.xyz[j], axis=1)
    d = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    return float(k * np.sum((d - d0) ** 2))
