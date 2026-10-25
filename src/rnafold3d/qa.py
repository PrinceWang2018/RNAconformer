"""RCV confidence model (optional): predicted local accuracy of an RNA 3D model from its coordinates alone.

For every nucleotide the model predicts its lDDT against the (unknown) true structure - the fraction of
interatomic distances to the neighbouring nucleotides (15 A) kept within 0.5, 1, 2 and 4 A, as used in CASP. The
mean over the nucleotides is the model's global confidence, which can rank models of the same sequence.

The network is small on purpose: hand-made structural features of each nucleotide (identity, pairing, stacking,
backbone pseudo-torsions, sugar pucker, glycosidic torsion, packing, clashes, chain integrity) plus their averages
over the spatial neighbourhood, two dense layers (tanh) and a sigmoid - about 3,200 parameters, numpy only,
milliseconds per model. Weights: data/qa_v1.npz, trained by tools/qa_model on structures that share no family with
the benchmark. It never changes a structure; it only scores.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from .annotate import canonical_pairs

FEATURE_VERSION = 1
WEIGHTS = Path(__file__).resolve().parent / 'data' / 'qa_v1.npz'
GLYCO = {'A': 'N9', 'G': 'N9', 'C': 'N1', 'U': 'N1'}
CHI_REF = {'A': 'C4', 'G': 'C4', 'C': 'C2', 'U': 'C2'}
RING = ("C4'", "O4'", "C1'", "C2'", "C3'")
PAIR_TYPES = ({'G', 'C'}, {'A', 'U'}, {'G', 'U'})


def _dih(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    n1 = np.linalg.norm(b1)
    if n1 < 1e-9:
        return None
    b1 = b1 / n1
    v = b0 - (b0 @ b1) * b1
    w = b2 - (b2 @ b1) * b1
    return math.atan2(float(np.cross(b1, v) @ w), float(v @ w))


def _pucker(atoms):
    """Pseudorotation phase (radians) of the ribose, or None."""
    if not all(a in atoms for a in RING):
        return None
    ring = [atoms[a] for a in ("C3'", "C4'", "O4'", "C1'", "C2'")]       # nu0..nu4 order of Altona-Sundaralingam
    nu = []
    for k in range(5):
        t = _dih(ring[(k + 4) % 5], ring[k], ring[(k + 1) % 5], ring[(k + 2) % 5])
        if t is None:
            return None
        nu.append(t)
    a = sum(nu[k] * math.cos(4 * math.pi * k / 5) for k in range(5))
    b = -sum(nu[k] * math.sin(4 * math.pi * k / 5) for k in range(5))
    return math.atan2(b, a)


def feature_names():
    base = ['A', 'C', 'G', 'U', 'rel_pos', 'log_len', 'paired', 'pair_gc', 'pair_au', 'pair_gu', 'stacked_pairs',
            'eta_sin', 'eta_cos', 'theta_sin', 'theta_cos', 'torsion_missing', 'chi_sin', 'chi_cos', 'pucker_sin',
            'pucker_cos', 'c4_8', 'c4_12', 'c4_16', 'burial', 'clashes', 'centre_dist', 'rg_ratio', 'link_strain',
            'chain_break', 'frac_paired', 'clash_density']
    nbr = ['nb_paired', 'nb_stacked', 'nb_c4_8', 'nb_burial', 'nb_clashes', 'nb_link_strain', 'nb_chain_break',
           'seq_paired', 'seq_pucker_cos', 'seq_burial']
    return base + nbr


def residue_features(seq, residues):
    """Feature matrix (n x F) for residues (list in chain order, one per nucleotide; atoms: name -> xyz)."""
    n = len(residues)
    F = np.zeros((n, len(feature_names())))
    if n == 0:
        return F
    pairs = canonical_pairs(residues)                 # by residue .index
    pos = {r.index: k for k, r in enumerate(residues)}
    partner = [-1] * n
    for i, j, _ in pairs:
        if i in pos and j in pos:
            partner[pos[i]], partner[pos[j]] = pos[j], pos[i]
    c4 = np.array([r.atoms.get("C4'", r.atoms.get("C1'", next(iter(r.atoms.values())) if r.atoms else np.zeros(3))) for r in residues], float)
    P = [r.atoms.get('P') for r in residues]
    O3 = [r.atoms.get("O3'") for r in residues]
    # all heavy atoms with their residue
    xyz, owner = [], []
    for k, r in enumerate(residues):
        for v in r.atoms.values():
            xyz.append(v)
            owner.append(k)
    xyz = np.asarray(xyz, float)
    owner = np.asarray(owner)
    centre = c4.mean(axis=0)
    rg = float(np.sqrt(np.mean(np.sum((c4 - centre) ** 2, axis=1)))) or 1.0
    d_c4 = np.linalg.norm(c4[:, None] - c4[None], axis=2)
    # burial and clashes: other-residue heavy atoms near each residue's atoms (cell grid via sorting is overkill here)
    burial = np.zeros(n)
    clashes = np.zeros(n)
    for k, r in enumerate(residues):
        if not r.atoms:
            continue
        near = np.where(d_c4[k] < 20.0)[0]
        mask = np.isin(owner, near) & (np.abs(owner - k) > 0)
        if not mask.any():
            continue
        own = np.asarray(list(r.atoms.values()), float)
        other = xyz[mask]
        d = np.linalg.norm(own[:, None] - other[None], axis=2)
        burial[k] = float((d < 4.5).any(axis=0).sum())
        seq_far = np.abs(owner[mask] - k) > 1
        clashes[k] = float(((d < 2.6) & seq_far[None, :]).any(axis=0).sum())
    frac_paired = sum(p >= 0 for p in partner) / n
    clash_density = clashes.sum() / n
    for k, r in enumerate(residues):
        b = seq[k] if k < len(seq) else r.base
        f = F[k]
        if b in 'ACGU':
            f['ACGU'.index(b)] = 1
        f[4] = k / max(1, n - 1)
        f[5] = math.log(n) / 5
        j = partner[k]
        if j >= 0:
            f[6] = 1
            pb = {b, seq[j] if j < len(seq) else residues[j].base}
            for t, pt in enumerate(PAIR_TYPES):
                if pb == pt:
                    f[7 + t] = 1
            f[10] = sum(1 for d in (-1, 1) if 0 <= k + d < n and 0 <= j - d < n and partner[k + d] == j - d)
        # pseudo-torsions eta (C4'(k-1)-P(k)-C4'(k)-P(k+1)) and theta (P(k)-C4'(k)-P(k+1)-C4'(k+1))
        eta = theta = None
        if 0 < k < n - 1 and P[k] is not None and P[k + 1] is not None:
            eta = _dih(c4[k - 1], P[k], c4[k], P[k + 1])
            theta = _dih(P[k], c4[k], P[k + 1], c4[k + 1])
        if eta is not None and theta is not None:
            f[11:15] = [math.sin(eta), math.cos(eta), math.sin(theta), math.cos(theta)]
        else:
            f[15] = 1
        at = r.atoms
        g, ref = GLYCO.get(b), CHI_REF.get(b)
        if g in at and ref in at and "O4'" in at and "C1'" in at:
            chi = _dih(at["O4'"], at["C1'"], at[g], at[ref])
            if chi is not None:
                f[16:18] = [math.sin(chi), math.cos(chi)]
        ph = _pucker(at)
        if ph is not None:
            f[18:20] = [math.sin(ph), math.cos(ph)]
        f[20] = (np.sum(d_c4[k] < 8) - 1) / 10
        f[21] = (np.sum(d_c4[k] < 12) - 1) / 20
        f[22] = (np.sum(d_c4[k] < 16) - 1) / 40
        f[23] = burial[k] / 60
        f[24] = min(clashes[k], 5) / 5
        f[25] = float(np.linalg.norm(c4[k] - centre)) / rg
        f[26] = rg / (5.5 * n ** (1 / 3))
        strain = 0.0
        for d in (-1, 1):
            if 0 <= k + d < n:
                strain = max(strain, abs(d_c4[k, k + d] - 6.0) / 6.0)
        f[27] = min(strain, 1.5)
        if k > 0 and P[k] is not None and O3[k - 1] is not None:
            f[28] = float(np.linalg.norm(P[k] - O3[k - 1]) > 2.0)
        f[29] = frac_paired
        f[30] = min(clash_density, 2) / 2
    # neighbourhood averages: spatial (C4' within 16 A, weight exp(-d/8)) and along the chain (+-2)
    W = np.exp(-d_c4 / 8.0) * (d_c4 < 16.0)
    np.fill_diagonal(W, 0.0)
    ws = W.sum(axis=1, keepdims=True)
    ws[ws == 0] = 1.0
    nb = (W @ F[:, [6, 10, 20, 23, 24, 27, 28]]) / ws
    F[:, 31:38] = nb
    for k in range(n):
        lo, hi = max(0, k - 2), min(n, k + 3)
        F[k, 38] = F[lo:hi, 6].mean()
        F[k, 39] = F[lo:hi, 19].mean()
        F[k, 40] = F[lo:hi, 23].mean()
    return F


# ----------------------------------------------------------------------------- network
@lru_cache(maxsize=2)
def load_weights(path=None):
    p = Path(path) if path else WEIGHTS
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    return {k: z[k] for k in z.files}


def available():
    return load_weights() is not None


def forward(X, w):
    Z = (X - w['mean']) / w['std']
    H1 = np.tanh(Z @ w['W1'] + w['b1'])
    H2 = np.tanh(H1 @ w['W2'] + w['b2'])
    return 1.0 / (1.0 + np.exp(-(H2 @ w['W3'] + w['b3']).ravel()))


def predict(seq, residues, weights=None):
    """{'per_residue': predicted lDDT per nucleotide (0-1), 'global': mean} or None without weights."""
    w = weights or load_weights()
    if w is None:
        return None
    if int(w.get('feature_version', np.array(FEATURE_VERSION))) != FEATURE_VERSION:
        return None
    X = residue_features(seq, residues)
    y = forward(X, w)
    return {'per_residue': y, 'global': float(y.mean()) if len(y) else None}


# ----------------------------------------------------------------------------- labels (training only)
def lddt_per_residue(ref, model, radius=15.0, thresholds=(0.5, 1.0, 2.0, 4.0)):
    """All-atom lDDT of every residue of model against ref ({number: Residue}); residues/atoms present in both."""
    keys = [k for k in sorted(ref) if k in model]
    atoms, owner = [], []
    for k in keys:
        for name, v in ref[k].atoms.items():
            if name in model[k].atoms:
                atoms.append((v, model[k].atoms[name]))
                owner.append(k)
    if not atoms:
        return {}
    R = np.array([a for a, _ in atoms], float)
    M = np.array([b for _, b in atoms], float)
    owner = np.array(owner)
    dR = np.linalg.norm(R[:, None] - R[None], axis=2)
    dM = np.linalg.norm(M[:, None] - M[None], axis=2)
    incl = (dR < radius) & (owner[:, None] != owner[None, :])
    diff = np.abs(dR - dM)
    score = np.zeros(len(R))
    for t in thresholds:
        score += ((diff < t) & incl).sum(axis=1)
    cnt = incl.sum(axis=1) * len(thresholds)
    out = {}
    for k in keys:
        sel = (owner == k) & (cnt > 0)
        if sel.any():
            out[k] = float(score[sel].sum() / cnt[sel].sum())
    return out
