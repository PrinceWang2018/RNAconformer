"""Junction geometry from experimental structures: template-free 3D interactions between helices.

Where helices meet (internal loops, bulges and multi-way junctions) their relative orientation is set by
coaxial stacking, the lengths of the connecting strands and loop-helix contacts. Instead of rules, RCV takes
this geometry from the same library of experimental structures that serves as templates: every junction of
every library chain is indexed by its number of helices and strand lengths. For each junction of a query
2D, a library junction with the same topology (any cyclic rotation; preferring similar loop sequence) is
drawn per trajectory, and C4'-C4' distances between the residues around it become soft restraints.

Junction representation: a cycle of ports; port m is the base pair at which the chain enters (a) and leaves
(b) helix m, so the helix continues with pairs (a + t, b - t); strands[m] is the number of unpaired residues
between b_m and a_{m+1}. Multiloops and internal loops/bulges closed by a pair are junctions; the exterior
loop is not (its helices are not held together by a closing pair).

Library chains homologous to the query (significant sequence hits) are excluded, so these restraints
never copy a related structure; that is the job of templates.py.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

MAX_DEPTH = 3


def nested_greedy(pairs):
    """Nested subset by removing the pair with most crossings until none cross (fast for large chains)."""
    p = np.array(sorted(pairs), int).reshape(-1, 2)
    if len(p) < 2:
        return [tuple(map(int, x)) for x in p]
    i, j = p[:, 0], p[:, 1]
    keep = np.ones(len(p), bool)
    while True:
        ii, jj = i[keep], j[keep]
        cr = ((ii[:, None] < ii[None, :]) & (ii[None, :] < jj[:, None]) & (jj[:, None] < jj[None, :]))
        cr = cr | cr.T
        cnt = cr.sum(axis=1)
        if not cnt.any():
            break
        idx = np.where(keep)[0]
        keep[idx[int(np.argmax(cnt))]] = False
    return [tuple(map(int, x)) for x in p[keep]]


def junctions(pairs):
    """Closed junctions of a nested 2D: [{'ports': [(a, b)...], 'strands': [...], 'depth': [...]}]."""
    partner = {}
    for i, j in pairs:
        partner[i] = j
        partner[j] = i
    out = []
    for i, j in sorted(pairs):
        if i > j:
            i, j = j, i
        kids = []
        k = i + 1
        while k < j:
            l = partner.get(k)
            if l is not None and l > k:
                kids.append((k, l))
                k = l + 1
            else:
                k += 1
        if not kids or (len(kids) == 1 and kids[0] == (i + 1, j - 1)):
            continue
        ports = [(j, i)] + kids
        strands = []
        for m in range(len(ports)):
            leave = ports[m][1]
            arrive = ports[(m + 1) % len(ports)][0]
            strands.append(arrive - leave - 1)
        depth = []
        for a, b in ports:
            t = 0
            while t < MAX_DEPTH and partner.get(a + t) == b - t:
                t += 1
            depth.append(t)
        out.append({'ports': ports, 'strands': strands, 'depth': depth})
    return out


def strand_residues(junction, m):
    leave = junction['ports'][m][1]
    return list(range(leave + 1, leave + 1 + junction['strands'][m]))


# ------------------------------------------------------------------------------------------- index
@lru_cache(maxsize=2)
def library_index(lib_path=None):
    """{(h, strands rotated): [(chain, junction, rotation)]} for every closed junction of the library."""
    from .templates import load_library
    lib = load_library(lib_path)
    if lib is None:
        return None, None
    entries = []
    index = {}
    for k in range(lib.n_chains):
        pairs = lib.chain_pairs(k)
        if len(pairs) < 2:
            continue
        nested = nested_greedy(pairs)
        for J in junctions(nested):
            if min(J['depth']) < 1:
                continue
            jid = len(entries)
            entries.append((k, J))
            h = len(J['ports'])
            for r in range(h):
                key = (h, tuple(J['strands'][r:] + J['strands'][:r]))
                index.setdefault(key, []).append((jid, r))
    return index, entries


def _mapping(Jq, Jl, rot):
    """Query residue -> library residue for a junction match; units label ports/strands."""
    h = len(Jq['ports'])
    mp, unit = {}, {}
    for m in range(h):
        ml = (m + rot) % h
        (a, b), (al, bl) = Jq['ports'][m], Jl['ports'][ml]
        for t in range(min(Jq['depth'][m], Jl['depth'][ml])):
            mp[a + t], unit[a + t] = al + t, ('p', m)
            mp[b - t], unit[b - t] = bl - t, ('p', m)
        for x, y in zip(strand_residues(Jq, m), strand_residues(Jl, ml)):
            mp[x], unit[x] = y, ('s', m)
    return mp, unit


def candidates(seq, Jq, lib, allowed, max_n=40):
    """Library junctions with the same helix count and strand lengths, scored by loop-sequence identity."""
    index, entries = library_index()
    if index is None:
        return []
    h = len(Jq['ports'])
    key = (h, tuple(Jq['strands']))
    out = []
    for jid, r in index.get(key, ()):
        k, Jl = entries[jid]
        if not allowed[k]:
            continue
        # rotation r of the library junction aligns its strand r with query strand 0
        mp, unit = _mapping(Jq, Jl, r)
        lseq = lib.chain_seq(k)
        same = sum(seq[q] == lseq[t] for q, t in mp.items())
        out.append((same / max(1, len(mp)), k, Jl, r))
    out.sort(key=lambda x: -x[0])
    return out[:max_n]


def junction_restraints(seq, pairs, lib, allowed, rng, beta=4.0, weight=0.5, cache=None):
    """Distance restraints for one trajectory: one library junction drawn for every junction of the 2D."""
    I, J, D, T, W = [], [], [], [], []
    used = []
    for Jq in junctions(pairs):
        key = (tuple(Jq['ports']), tuple(Jq['strands']))
        cands = cache.get(key) if cache is not None else None
        if cands is None:
            cands = candidates(seq, Jq, lib, allowed)
            if cache is not None:
                cache[key] = cands
        if not cands:
            continue
        sims = np.array([c[0] for c in cands])
        p = np.exp(beta * (sims - sims.max()))
        sim, k, Jl, r = cands[int(rng.choice(len(cands), p=p / p.sum()))]
        mp, unit = _mapping(Jq, Jl, r)
        xyz = lib.chain_xyz(k)
        qs = [q for q in sorted(mp) if not np.isnan(xyz[mp[q], 0])]
        X = np.array([xyz[mp[q]] for q in qs])
        for a in range(len(qs)):
            for b in range(a + 1, len(qs)):
                qa, qb = qs[a], qs[b]
                if unit[qa] == unit[qb] and (unit[qa][0] == 'p' or abs(qa - qb) < 2):
                    continue
                d = float(np.linalg.norm(X[a] - X[b]))
                I.append(qa)
                J.append(qb)
                D.append(round(d, 2))
                T.append(round(1.0 + 0.05 * d, 2))
                W.append(weight)
        m = lib.meta[k]
        used.append({'junction': [(a + 1, b + 1) for a, b in Jq['ports']], 'strands': Jq['strands'],
                     'source': f"{m['pdb']}_{m['chain']}", 'identity': round(float(sim), 2)})
    return {'i': I, 'j': J, 'd': D, 'tol': T, 'w': W}, used
