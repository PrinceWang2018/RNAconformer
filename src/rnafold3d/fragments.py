"""Motif-level 3D fragments from experimental structures, grafted onto the fold tree.

Between whole-chain templates (templates.py) and junction geometry (motifs.py) sits the fragment: a hairpin
loop, a bulge, an internal loop such as a k-turn or a loop E, taken from any library chain whose sequence
matches the motif. This is what FARFAR2 and Vfold do with their fragment libraries, and it matters most where
there is no template.

For every hairpin loop of a 2D hypothesis, a library hairpin with the same loop length is drawn per trajectory,
strongly preferring the same loop sequence (a GNRA tetraloop gets a GNRA tetraloop, a UUCG a UUCG, a T-loop a
T-loop). Its backbone torsions become the starting torsions of the loop, which is then closed from them, and
C4'-C4' distances around it become soft restraints. Internal loops and multi-way junctions get the same
treatment through motifs.junction_restraints (library junction of the same topology, sequence-weighted), which
here also passes on the torsions of the strands and of the residues where the strands leave and enter the
helices: the fragment is grafted, not only restrained.

As for junction geometry, library chains homologous to the query are excluded, so a fragment never copies a
related structure, and the restraint energy guides the sampling but is not part of the ranking score.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

MIN_LOOP, MAX_LOOP = 3, 12
CANON = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}


def hairpin_loops(pairs):
    """Closing pairs (i, j) of the hairpin loops of a nested 2D (no pair inside)."""
    ps = sorted((min(i, j), max(i, j)) for i, j in pairs)
    inner = {}
    for i, j in ps:
        inner[(i, j)] = False
    for i, j in ps:
        for k, l in ps:
            if i < k and l < j:
                inner[(i, j)] = True
                break
    return [(i, j) for i, j in ps if not inner[(i, j)] and MIN_LOOP <= j - i - 1 <= MAX_LOOP]


@lru_cache(maxsize=2)
def hairpin_index(lib_path=None):
    """{loop length: [(chain, i, j)]} for library hairpins closed by a canonical pair with a stacked pair below
    it and complete C4' atoms and torsions over the loop."""
    from .templates import load_library
    from .motifs import library_nested
    lib = load_library(lib_path)
    if lib is None:
        return None
    index = {}
    for k in range(lib.n_chains):
        if not len(lib.chain_pairs(k)):
            continue
        nested = library_nested(lib, k)
        pset = set(nested)
        seq = lib.chain_seq(k)
        x = tor = None
        for i, j in hairpin_loops(nested):
            if (i - 1, j + 1) not in pset or seq[i] + seq[j] not in CANON:
                continue
            if x is None:
                x, tor = lib.chain_xyz(k), lib.chain_tor(k)
            if np.isnan(x[i - 1:j + 2, 0]).any() or np.isnan(tor[i + 1:j, :6]).any():
                continue
            index.setdefault(j - i - 1, []).append((k, i, j))
    return index


def hairpin_candidates(seq, i, j, lib, allowed, max_n=40):
    """Library hairpins for the query hairpin closed by (i, j): (similarity, chain, i, j), best first.

    Similarity: identity over the loop residues and the closing pair, the loop counting double; library
    hairpins with a pyrimidine/purine pattern different at more than half of the loop are not used."""
    index = hairpin_index()
    if not index:
        return []
    n = j - i - 1
    q = seq[i:j + 1]
    out = []
    for k, a, b in index.get(n, ()):
        if not allowed[k]:
            continue
        t = lib.chain_seq(k)[a:b + 1]
        loop_same = sum(x == y for x, y in zip(q[1:-1], t[1:-1]))
        cls_same = sum((x in 'AG') == (y in 'AG') for x, y in zip(q[1:-1], t[1:-1]))
        if cls_same * 2 < n:
            continue
        sim = (2.0 * loop_same + (q[0] == t[0]) + (q[-1] == t[-1])) / (2.0 * n + 2.0)
        out.append((sim, k, a, b))
    out.sort(key=lambda c: (-c[0], c[1]))
    return out[:max_n]


def _torsion_row(vals, keep):
    """Library torsion row -> list for FoldModel.apply_torsions (None = keep the sampled value)."""
    return [None if (nm not in keep or np.isnan(v)) else round(float(v), 1)
            for nm, v in zip(('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi'), vals)]


ALL = {'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi'}


def hairpin_fragments(seq, pairs, lib, allowed, rng, beta=10.0, weight=0.5, cache=None):
    """One library hairpin per hairpin loop of the 2D: restraints, torsions and a record of what was used."""
    I, J, D, T, W = [], [], [], [], []
    torsions, used, loops = {}, [], []
    for i, j in hairpin_loops(pairs):
        key = ('hairpin', i, j)
        cands = cache.get(key) if cache is not None else None
        if cands is None:
            cands = hairpin_candidates(seq, i, j, lib, allowed)
            if cache is not None:
                cache[key] = cands
        if not cands:
            continue
        sims = np.array([c[0] for c in cands])
        p = np.exp(beta * (sims - sims.max()))
        sim, k, a, b = cands[int(rng.choice(len(cands), p=p / p.sum()))]
        x, tor = lib.chain_xyz(k), lib.chain_tor(k)
        q = list(range(i - 1, j + 2))
        t = list(range(a - 1, b + 2))
        X = x[t]
        for u in range(len(q)):
            for v in range(u + 2, len(q)):
                if not (i < q[u] < j or i < q[v] < j):
                    continue                          # both in the stem: rigid anyway
                d = float(np.linalg.norm(X[u] - X[v]))
                I.append(q[u])
                J.append(q[v])
                D.append(round(d, 2))
                T.append(round(1.0 + 0.05 * d, 2))
                W.append(weight)
        for r in range(i + 1, j):
            torsions[r] = _torsion_row(tor[a + (r - i)], ALL)
        torsions[i] = _torsion_row(tor[a], {'epsilon', 'zeta'})
        torsions[j] = _torsion_row(tor[b], {'alpha'})
        loops.append((i, j))
        m = lib.meta[k]
        used.append({'hairpin': (i + 1, j + 1), 'loop': seq[i + 1:j], 'source': f"{m['pdb']}_{m['chain']}",
                     'source_loop': lib.chain_seq(k)[a + 1:b], 'source_residues': (a + 1, b + 1),
                     'identity': round(float(sim), 2)})
    return {'i': I, 'j': J, 'd': D, 'tol': T, 'w': W}, torsions, loops, used


def merge_restraints(*drs):
    out = {'i': [], 'j': [], 'd': [], 'tol': [], 'w': []}
    for dr in drs:
        if not dr:
            continue
        n = len(dr.get('i', ()))
        for key in ('i', 'j', 'd', 'tol'):
            out[key] += list(dr[key])
        out['w'] += list(dr.get('w') or [1.0] * n)
    return out
