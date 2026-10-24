"""Geometric annotation of canonical base pairs in 3D RNA coordinates.

Detects cis Watson-Crick pairs (A-U, G-C) and G-U wobble pairs from heavy-atom
hydrogen-bond distances, base coplanarity and C1'-C1' distance. Non-canonical
pairs and stacking are not annotated; INF computed from these annotations is
therefore an *INF over canonical pairs* (INF-WC), not the full INF-all.
"""
from __future__ import annotations

import numpy as np

RING_ATOMS = {
    'A': ('N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'),
    'G': ('N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3', 'C4'),
    'C': ('N1', 'C2', 'N3', 'C4', 'C5', 'C6'),
    'U': ('N1', 'C2', 'N3', 'C4', 'C5', 'C6'),
}
# Donor/acceptor atom pairs (atom on first base, atom on second base).
HBONDS = {
    ('G', 'C'): (('N1', 'N3'), ('O6', 'N4'), ('N2', 'O2')),
    ('A', 'U'): (('N1', 'N3'), ('N6', 'O4')),
    ('G', 'U'): (('N1', 'O2'), ('O6', 'N3')),
}


def base_plane(atoms, base):
    """Return (centroid, unit normal) of the base ring atoms or None if incomplete."""
    names = RING_ATOMS[base]
    if not all(n in atoms for n in names):
        return None
    pts = np.array([atoms[n] for n in names])
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c)
    return c, vt[2]


def _pair_score(ri, rj):
    """Return (is_pair, pair_type, mean H-bond distance) for residues ri, rj (Residue-like)."""
    key = (ri.base, rj.base)
    swap = False
    if key not in HBONDS:
        key = (rj.base, ri.base)
        swap = True
    if key not in HBONDS:
        return False, None, None
    a, b = (rj, ri) if swap else (ri, rj)
    dists = []
    for x, y in HBONDS[key]:
        if x not in a.atoms or y not in b.atoms:
            return False, None, None
        dists.append(float(np.linalg.norm(a.atoms[x] - b.atoms[y])))
    if dists[0] > 3.5:                       # central N1-N3 / N1-O2 bond is mandatory
        return False, None, None
    good = sum(d <= 3.5 for d in dists)
    if good < min(2, len(dists)):
        return False, None, None
    pa, pb = base_plane(a.atoms, a.base), base_plane(b.atoms, b.base)
    if pa is None or pb is None:
        return False, None, None
    cos = abs(float(np.dot(pa[1], pb[1])))
    if cos < np.cos(np.radians(40)):
        return False, None, None
    if "C1'" in a.atoms and "C1'" in b.atoms:
        c1 = float(np.linalg.norm(a.atoms["C1'"] - b.atoms["C1'"]))
        if not 9.0 <= c1 <= 12.0:
            return False, None, None
    return True, ''.join(sorted((ri.base, rj.base))), float(np.mean(dists))


def canonical_pairs(residues):
    """Annotate canonical pairs; residues are objects with .index, .base, .atoms.

    Returns sorted list of (index_i, index_j, type) with index_i < index_j and each
    residue in at most one pair (conflicts resolved by shortest mean H-bond).
    """
    res = [r for r in residues if r.base in 'ACGU']
    # Pre-filter with C1' distances when available (fast for 200 nt).
    cands = []
    for a in range(len(res)):
        for b in range(a + 1, len(res)):
            ri, rj = res[a], res[b]
            if abs(ri.index - rj.index) < 4:
                continue
            if "C1'" in ri.atoms and "C1'" in rj.atoms:
                if np.linalg.norm(ri.atoms["C1'"] - rj.atoms["C1'"]) > 12.5:
                    continue
            ok, typ, score = _pair_score(ri, rj)
            if ok:
                cands.append((score, ri.index, rj.index, typ))
    cands.sort()
    used = set()
    pairs = []
    for score, i, j, typ in cands:
        if i in used or j in used:
            continue
        used.update((i, j))
        pairs.append((min(i, j), max(i, j), typ))
    return sorted(pairs)


def nested_subset(pairs):
    """Maximum-cardinality non-crossing subset of (i, j) pairs (O(N^2) dynamic programming).

    Returns (nested_pairs, removed_pairs). Each residue must appear in at most one pair.
    """
    pairs = sorted((min(p[0], p[1]), max(p[0], p[1])) for p in pairs)
    if not pairs:
        return [], []
    partner = {}
    for i, j in pairs:
        partner[i] = j
        partner[j] = i
    idx = sorted(partner)
    lo, hi = idx[0], idx[-1]
    pos = {v: k for k, v in enumerate(range(lo, hi + 1))}
    n = hi - lo + 1
    best = [[0] * (n + 1) for _ in range(n + 1)]
    choice = [[0] * (n + 1) for _ in range(n + 1)]
    for length in range(1, n + 1):
        for a in range(0, n - length + 1):
            b = a + length - 1
            i = a + lo
            val = best[a + 1][b] if a + 1 <= b else 0
            ch = 0
            k = partner.get(i)
            if k is not None and i < k <= b + lo:
                kk = k - lo
                inner = best[a + 1][kk - 1] if a + 1 <= kk - 1 else 0
                outer = best[kk + 1][b] if kk + 1 <= b else 0
                if 1 + inner + outer > val:
                    val = 1 + inner + outer
                    ch = 1
            best[a][b] = val
            choice[a][b] = ch
    kept = []
    stack = [(0, n - 1)]
    while stack:
        a, b = stack.pop()
        if a > b:
            continue
        if choice[a][b]:
            i = a + lo
            kk = partner[i] - lo
            kept.append((i, partner[i]))
            stack.append((a + 1, kk - 1))
            stack.append((kk + 1, b))
        else:
            stack.append((a + 1, b))
    kept = sorted(kept)
    removed = sorted(set(pairs) - set(kept))
    return kept, removed


def pairs_to_dotbracket(n, pairs, offset=1):
    """Pairs use 1-based indices by default; crossing pairs get [] {} <> levels."""
    chars = ['.'] * n
    levels = []
    brackets = ('()', '[]', '{}', '<>')
    for i, j in sorted((min(p[0], p[1]), max(p[0], p[1])) for p in pairs):
        lvl = next((k for k, assigned in enumerate(levels)
                    if all(not (a < i < b < j or i < a < j < b) for a, b in assigned)), None)
        if lvl is None:
            if len(levels) == len(brackets):
                raise ValueError('too many pseudoknot levels')
            levels.append([])
            lvl = len(levels) - 1
        levels[lvl].append((i, j))
        chars[i - offset] = brackets[lvl][0]
        chars[j - offset] = brackets[lvl][1]
    return ''.join(chars)


def dotbracket_to_pairs(structure, offset=1):
    stacks = {o: [] for o in '([{<'}
    close = {')': '(', ']': '[', '}': '{', '>': '<'}
    pairs = []
    for k, ch in enumerate(structure):
        if ch in stacks:
            stacks[ch].append(k + offset)
        elif ch in close:
            if not stacks[close[ch]]:
                raise ValueError(f'unbalanced bracket at {k + 1}')
            pairs.append((stacks[close[ch]].pop(), k + offset))
        elif ch not in '.-_,:':
            raise ValueError(f'unsupported character {ch!r}')
    if any(stacks.values()):
        raise ValueError('unbalanced brackets')
    return sorted(pairs)
