"""Topology-aware ensemble map (Energy tab map, Pathway tab, plot studio, interactive HTML page).

Conformations are placed so that map distances follow their pairwise C4' RMSD under a topology constraint: every
secondary-structure topology - the abstract shape of its 2D structure, helix nesting only (level 5 of Giegerich,
Voss & Rehmsmeier 2004) - occupies its own disk that no other topology enters. The disks are arranged by how close
the topologies are in 3D; inside a disk the conformations are placed by RMSD and grouped by 2D state.

  1. shape and 2D state of every conformation;
  2. inside each topology: stress majorisation (SMACOF) of the RMSD with a separation term between different 2D
     states; the disk radius is the extent of that layout;
  3. disk centres: SMACOF on the median RMSD between topologies, every target distance at least the sum of the two
     radii plus a margin, then overlap removal (no disk touches another);
  4. orientation: the arrangement and every local layout are rotated or reflected (Procrustes) to agree with
     classical MDS of all conformations;
  5. joint refinement: SMACOF over all conformations, each point kept inside its disk.
Map units are Angstrom of RMSD. Quality is reported for this map and for plain classical MDS: normalised stress,
Spearman correlation of distances, trustworthiness and continuity of neighbourhoods (Venna & Kaski 2001) and the
topology purity of each point's map neighbours.
"""
from __future__ import annotations

import math

import numpy as np

from .switch import ensemble_states, pairs_of

VERSION = 1


# ----------------------------------------------------------------------------- shapes
def shape_of(dbn):
    """Abstract shape (helix nesting only; stacks, bulges and interior loops merged): '[[][]]', '_' if unpaired."""
    pairs = pairs_of(dbn)
    if not pairs:
        return '_'
    pairs.sort()
    children = {p: [] for p in pairs}
    roots, stack = [], []
    for i, j in sorted(pairs, key=lambda p: p[0]):
        while stack and stack[-1][1] < i:
            stack.pop()
        (children[stack[-1]] if stack else roots).append((i, j))
        stack.append((i, j))

    def sh(p):
        while len(children[p]) == 1:
            p = children[p][0]
        return '[' + ''.join(sh(c) for c in children[p]) + ']'
    return ''.join(sh(r) for r in roots)


def levenshtein(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]


# ----------------------------------------------------------------------------- layout tools
def classical_mds(D, dim=2):
    n = len(D)
    if n == 1:
        return np.zeros((1, dim))
    J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (np.asarray(D, float) ** 2) @ J
    w, v = np.linalg.eigh(B)
    order = np.argsort(w)[::-1][:dim]
    w = np.clip(w[order], 0, None)
    X = v[:, order] * np.sqrt(w)
    if X.shape[1] < dim:
        X = np.hstack([X, np.zeros((n, dim - X.shape[1]))])
    return X


def _pdist(X):
    return np.sqrt(np.maximum(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1), 0.0))


def smacof(D, X, iters=300, project=None, tol=1e-7):
    """Unweighted stress majorisation (Guttman transform) from X; project(X) after every step if given."""
    D = np.asarray(D, float)
    n = len(D)
    if n < 2:
        return X
    prev = None
    for _ in range(iters):
        d = _pdist(X)
        with np.errstate(divide='ignore', invalid='ignore'):
            B = np.where(d > 1e-12, -D / d, 0.0)
        np.fill_diagonal(B, 0.0)
        np.fill_diagonal(B, -B.sum(axis=1))
        X = B @ X / n
        if project:
            X = project(X)
        s = float(((_pdist(X) - D) ** 2).sum())
        if prev is not None and abs(prev - s) <= tol * max(prev, 1e-12):
            break
        prev = s
    return X


def procrustes_rotation(L, Y):
    """Orthogonal matrix Q (rotation or reflection) minimising |L Q - Y| for centred L, Y (n x 2)."""
    if len(L) < 2 or np.allclose(L, 0) or np.allclose(Y, 0):
        return np.eye(2)
    U, _, Vt = np.linalg.svd(L.T @ Y)
    return U @ Vt


def _ranks(x):
    order = np.argsort(x, kind='mergesort')
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    # average ties
    xs = np.asarray(x)[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = (i + j) / 2
        i = j + 1
    return r


def quality(D, X, cls, k=None):
    """Stress-1, Spearman rho of distances, trustworthiness, continuity and topology purity of a map X."""
    n = len(D)
    if n < 3:
        return None
    D = np.asarray(D, float)
    d = _pdist(X)
    iu = np.triu_indices(n, 1)
    dv, Dv = d[iu], D[iu]
    s = float((dv * Dv).sum() / max((dv * dv).sum(), 1e-12))     # best scale of the map
    stress = float(math.sqrt(((s * dv - Dv) ** 2).sum() / max((Dv ** 2).sum(), 1e-12)))
    rho = float(np.corrcoef(_ranks(Dv), _ranks(dv))[0, 1]) if len(Dv) > 2 and Dv.std() > 0 and dv.std() > 0 else None
    k = k or max(2, min(7, (n - 1) // 4))
    k = min(k, max(1, (n - 1) // 2))

    def ranks_rows(M):
        R = np.zeros((n, n), int)
        for i in range(n):
            o = np.argsort(np.where(np.arange(n) == i, -1, M[i]), kind='mergesort')
            R[i, o] = np.arange(n)          # self gets rank 0, neighbours 1..n-1
        return R
    RD, Rd = ranks_rows(D), ranks_rows(d)
    nn_D = [set(np.where((RD[i] >= 1) & (RD[i] <= k))[0]) for i in range(n)]
    nn_d = [set(np.where((Rd[i] >= 1) & (Rd[i] <= k))[0]) for i in range(n)]
    norm = 2.0 / (n * k * (2 * n - 3 * k - 1)) if (2 * n - 3 * k - 1) > 0 else 0.0
    trust = 1 - norm * sum(RD[i, j] - k for i in range(n) for j in nn_d[i] - nn_D[i])
    cont = 1 - norm * sum(Rd[i, j] - k for i in range(n) for j in nn_D[i] - nn_d[i])
    purity = float(np.mean([np.mean([cls[j] == cls[i] for j in nn_d[i]]) for i in range(n)])) if n > 1 else 1.0
    return {'stress': round(stress, 4), 'spearman': round(rho, 4) if rho is not None else None,
            'trustworthiness': round(float(trust), 4), 'continuity': round(float(cont), 4),
            'topology_purity': round(purity, 4), 'k': int(k)}


# ----------------------------------------------------------------------------- the map
def layout(rmsd, dbns, sep_factor=0.1, margin_factor=0.4):
    """Topology-aware positions for conformations with RMSD matrix rmsd and nested 2D structures dbns."""
    D = np.asarray(rmsd, float)
    n = len(D)
    shapes = [shape_of(s) for s in dbns]
    shape_ids = list(dict.fromkeys(shapes))
    cls = [shape_ids.index(s) for s in shapes]
    state_ids = list(dict.fromkeys(dbns))
    state = [state_ids.index(s) for s in dbns]
    nz = D[np.triu_indices(n, 1)] if n > 1 else np.array([1.0])
    G = float(np.median(nz[nz > 0])) if (nz > 0).any() else 1.0
    Y = classical_mds(D) if n > 1 else np.zeros((1, 2))
    members = [[i for i in range(n) if cls[i] == c] for c in range(len(shape_ids))]
    # 2. local layouts
    local, radius = [], []
    for c, mem in enumerate(members):
        Dl = D[np.ix_(mem, mem)].copy()
        if len(mem) > 1:
            Dw = Dl[np.triu_indices(len(mem), 1)]
            lam = max(0.5, sep_factor * float(np.median(Dw))) if len(Dw) else 0.5
            st = np.array([state[i] for i in mem])
            Dl = Dl + lam * (st[:, None] != st[None, :])
            np.fill_diagonal(Dl, 0.0)
            L = smacof(Dl, classical_mds(Dl))
            L -= L.mean(axis=0)
            Yc = Y[mem] - Y[mem].mean(axis=0)
            L = L @ procrustes_rotation(L, Yc)
        else:
            L = np.zeros((1, 2))
        local.append(L)
        r = float(np.max(np.linalg.norm(L, axis=1))) if len(mem) > 1 else 0.0
        radius.append(r + max(0.15 * G, 0.4))
    # 3. disk centres
    m = len(members)
    margin = margin_factor * G
    if m == 1:
        C = np.zeros((1, 2))
    else:
        Dc = np.zeros((m, m))
        for a in range(m):
            for b in range(a + 1, m):
                med = float(np.median(D[np.ix_(members[a], members[b])]))
                Dc[a, b] = Dc[b, a] = max(med, radius[a] + radius[b] + margin)
        cent = np.array([Y[mem].mean(axis=0) for mem in members])
        C = smacof(Dc, cent.copy() if not np.allclose(cent, cent[0]) else classical_mds(Dc))
        for _ in range(500):                      # overlap removal
            moved = False
            for a in range(m):
                for b in range(a + 1, m):
                    v = C[b] - C[a]
                    dist = float(np.linalg.norm(v))
                    need = radius[a] + radius[b] + margin
                    if dist < need - 1e-9:
                        u = v / dist if dist > 1e-9 else np.array([math.cos(a + b), math.sin(a + b)])
                        push = (need - dist) / 2 + 1e-6
                        C[a] -= u * push
                        C[b] += u * push
                        moved = True
            if not moved:
                break
        C -= C.mean(axis=0)
        cc = cent - cent.mean(axis=0)
        C = C @ procrustes_rotation(C, cc)
    X = np.zeros((n, 2))
    for c, mem in enumerate(members):
        X[mem] = C[c] + local[c]
    # 5. joint refinement inside the disks
    cen = C[cls]
    rad = np.array([radius[c] for c in cls]) - 0.1 * margin

    anchor = cen.mean(axis=0)

    def project(Z):
        Z = Z - Z.mean(axis=0) + anchor            # the Guttman transform centres the cloud; keep the disk frame
        v = Z - cen
        dist = np.linalg.norm(v, axis=1)
        over = dist > rad
        Z[over] = cen[over] + v[over] / dist[over, None] * rad[over, None]
        return Z
    if m > 1 and n > 2:
        X = smacof(D, X, iters=150, project=project)
    q_map = quality(D, X, cls)
    q_mds = quality(D, Y, cls)
    edges = []
    for a in range(m):
        for b in range(a + 1, m):
            if levenshtein(shape_ids[a], shape_ids[b]) <= 2:
                edges.append({'a': a, 'b': b, 'kind': 'one helix', 'rmsd': round(float(np.median(D[np.ix_(members[a], members[b])])), 2)})
    if m > 1:                                     # keep the topology graph connected (minimum spanning tree on RMSD)
        linked = {0}
        while len(linked) < m:
            best = None
            for a in linked:
                for b in range(m):
                    if b in linked:
                        continue
                    v = float(np.median(D[np.ix_(members[a], members[b])]))
                    if best is None or v < best[0]:
                        best = (v, a, b)
            linked.add(best[2])
            if not any({e['a'], e['b']} == {best[1], best[2]} for e in edges):
                edges.append({'a': min(best[1], best[2]), 'b': max(best[1], best[2]), 'kind': 'nearest', 'rmsd': round(best[0], 2)})
    return {
        'version': VERSION,
        'positions': [[round(float(x), 3), round(float(y), 3)] for x, y in X],
        'mds': [[round(float(x), 3), round(float(y), 3)] for x, y in Y],
        'class_of': cls, 'state_of': state,
        'classes': [{'shape': shape_ids[c], 'label': f'T{c + 1}', 'members': members[c],
                     'center': [round(float(C[c][0]), 3), round(float(C[c][1]), 3)], 'radius': round(radius[c], 3),
                     'states': sorted({state[i] for i in members[c]})} for c in range(m)],
        'states': [{'dbn': s, 'members': [i for i in range(n) if state[i] == k], 'class': cls[state.index(k)]} for k, s in enumerate(state_ids)],
        'edges': edges,
        'quality': {'map': q_map, 'mds': q_mds},
        'units': 'Å (C4′ RMSD)',
        'method': 'topology-aware stress layout: one disk per abstract shape (level 5), SMACOF inside and between '
                  'disks, Procrustes orientation to classical MDS, joint refinement inside the disks',
    }


def layout_doc(doc):
    """Map of an ensemble document (2D states as in the Switch tab)."""
    n = len(doc['conformations'])
    dbn = [''] * n
    for st in ensemble_states(doc):
        for k in st['conformations']:
            dbn[k] = st['dbn']
    rmsd = doc.get('ensemble', {}).get('rmsd_c4') or [[0.0] * n for _ in range(n)]
    return layout(rmsd, dbn)
