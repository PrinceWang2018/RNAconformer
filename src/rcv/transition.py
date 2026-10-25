"""Transition analysis between two conformations of an ensemble (Pathway tab).

Secondary-structure level (ViennaRNA Turner energies under the chosen conditions, see conditions.py):
  * the direct refolding path from A to B (findpath heuristic, Flamm et al. 2001: only pairs of A are opened and
    pairs of B closed) with its barrier, saddle and the order of base-pair events;
  * indirect routes through intermediates - the ensemble's own 2D states and local minima reached by gradient
    descent from Boltzmann-sampled structures. The route whose highest saddle is lowest (minimax route) and the
    flux-ranked routes of transition path theory (Metzner, Schütte & Vanden-Eijnden 2009) on the network of states,
    with Arrhenius rates k = k0 exp(-(S - E)/RT) over the findpath saddles S;
  * the cooperative unit - the base pairs that must open together on the climb to the highest saddle - and the
    nucleus, the pairs that form right after it;
  * an in-silico phi-value scan: for every point mutation the change of the barrier (ddG_TS) and of the A/B balance
    (ddG_eq), phi = ddG_TS / ddG_eq where |ddG_eq| >= 0.6 kcal/mol (the route is kept, its energies re-evaluated).
3D level (coordinates of the two conformations):
  * rigid modules: groups of residues that move as one body, seeded by helices present in both conformations and
    grown by merging neighbouring groups while their joint rigid fit stays within 2 A; hinge residues; rotation,
    translation and screw axis of every module relative to the largest one;
  * when each module moves, from the 2D events of its residues (for the mechanism morph).

Barriers are those of the paths found (heuristics give upper bounds of the lowest possible barrier). Rates and
half-lives are order-of-magnitude estimates that scale with the prefactor k0.
"""
from __future__ import annotations

import heapq
import math
import threading

import numpy as np

from .conditions import Conditions
from .switch import compatible, ensemble_states, pairs_of

K0_DEFAULT = 1.0e6            # s^-1, attempt frequency for the Arrhenius rates (adjustable)
PHI_MIN_DDG = 0.6             # kcal/mol; below this |ddG_eq| phi values are not interpretable
RIGID_TAU = 2.0               # A; a group is rigid if its A->B fit RMSD is below this
SAMPLE_SEED = 20261025        # ViennaRNA random seed of the Boltzmann samples (reproducible results)
_RNG_LOCK = threading.Lock()  # the ViennaRNA generator is global; seed and sample together


# ----------------------------------------------------------------------------- structures
def pair_set(dbn):
    return set(pairs_of(dbn))


def bp_distance(a, b):
    return len(pair_set(a) ^ pair_set(b))


def helix_segments(pairs):
    """Group pairs into stacked runs (i, j), (i+1, j-1), ...; returns lists of pairs, 5' first."""
    left = sorted(pairs)
    pset = set(left)
    seen, out = set(), []
    for p in left:
        if p in seen:
            continue
        i, j = p
        while (i - 1, j + 1) in pset and (i - 1, j + 1) not in seen:
            i, j = i - 1, j + 1
        run = []
        while (i, j) in pset:
            run.append((i, j))
            seen.add((i, j))
            i, j = i + 1, j - 1
        out.append(run)
    return sorted(out, key=lambda r: r[0])


class Landscape:
    """Energies, refolding paths and local minima of one sequence under given conditions (cached)."""

    def __init__(self, seq, cond: Conditions, width=None):
        self.seq = seq
        self.cond = cond
        self.fc = cond.fold_compound(seq)
        self.sites = cond.motif_sites(seq)
        n = len(seq)
        self.width = width or (16 if n <= 40 else 10 if n <= 100 else 6)
        self._e, self._paths = {}, {}

    def E(self, s):
        if s not in self._e:
            self._e[s] = round(self.cond.energy(self.fc, self.seq, s, self.sites), 2)
        return self._e[s]

    def path(self, a, b):
        """Lower-saddle findpath of the two directions, oriented a -> b: [(structure, energy)]."""
        if a == b:
            return [(a, self.E(a))]
        key = (a, b) if a < b else (b, a)
        if key not in self._paths:
            best = None
            for x, y in ((key[0], key[1]), (key[1], key[0])):
                p = self.fc.path_findpath(x, y, width=self.width)
                steps = [(st.s, self.E(st.s)) for st in p]
                if not steps or steps[0][0] != x:
                    steps = [(x, self.E(x))] + steps
                if steps[-1][0] != y:
                    steps.append((y, self.E(y)))
                if x != key[0]:
                    steps = steps[::-1]
                if best is None or max(e for _, e in steps) < max(e for _, e in best):
                    best = steps
            self._paths[key] = best
        p = self._paths[key]
        return p if key == (a, b) else p[::-1]

    def saddle(self, a, b):
        return max(e for _, e in self.path(a, b))

    def descend(self, s, path=False):
        """Steepest descent with the insertion/deletion move set to a local minimum (Turner energies); with path=True
        the visited structures (start ... minimum)."""
        import RNA
        fc = self.fc
        key = ('descent', s)
        if key in self._paths:
            visited = self._paths[key]
            return visited if path else visited[-1]
        s = list(s)
        visited = [''.join(s)]
        for _ in range(4 * len(self.seq)):
            cur = ''.join(s)
            best = None
            for mv in fc.neighbors(RNA.ptable(cur)):
                d = fc.eval_move(cur, mv.pos_5, mv.pos_3)
                if d < -1e-6 and (best is None or d < best[0] or (d == best[0] and (mv.pos_5, mv.pos_3) < best[1:])):
                    best = (d, mv.pos_5, mv.pos_3)
            if best is None:
                break
            _, i, j = best
            if i > 0:
                s[i - 1], s[j - 1] = '(', ')'
            else:
                s[-i - 1], s[-j - 1] = '.', '.'
            visited.append(''.join(s))
        self._paths[key] = visited
        return visited if path else visited[-1]

    def sample_minima(self, n_samples):
        fc = self.cond.fold_compound(self.seq, uniq_ml=True)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        import RNA
        with _RNG_LOCK:
            RNA.init_rand(SAMPLE_SEED)
            samples = list(fc.pbacktrack(int(n_samples)))
        seen = {}
        for s in samples:
            if s not in seen:
                seen[s] = self.descend(s)
        return sorted(set(seen.values()), key=self.E)


# ----------------------------------------------------------------------------- network, routes, flux
def _minimax_route(nodes, edges, a, b):
    """Route a -> b minimising the highest saddle, then the number of legs. edges: {(i, j): saddle}."""
    adj = {i: [] for i in range(len(nodes))}
    for (i, j), s in edges.items():
        adj[i].append((j, s))
        adj[j].append((i, s))
    best = {a: (-math.inf, 0)}
    prev = {}
    heap = [(-math.inf, 0, a)]
    while heap:
        m, legs, u = heapq.heappop(heap)
        if (m, legs) > best.get(u, (math.inf, 0)):
            continue
        if u == b:
            break
        for v, s in adj[u]:
            cand = (max(m, s), legs + 1)
            if cand < best.get(v, (math.inf, math.inf)):
                best[v] = cand
                prev[v] = u
                heapq.heappush(heap, (cand[0], cand[1], v))
    if b not in prev and a != b:
        return None
    route = [b]
    while route[-1] != a:
        route.append(prev[route[-1]])
    return route[::-1]


def transition_path_theory(E, edges, a, b, kT, k0=K0_DEFAULT, max_paths=6):
    """Committors, reactive flux, rate and flux-ranked pathways on a network with Arrhenius rates over saddles.

    E: node energies (kcal/mol); edges {(i, j): saddle}. Detailed balance holds by construction (the saddle is
    shared by both directions)."""
    n = len(E)
    E = np.asarray(E, float)
    beta = 1.0 / kT
    nbr = {i: {} for i in range(n)}
    for (i, j), s in edges.items():
        s = max(s, E[i], E[j])
        nbr[i][j] = nbr[j][i] = s
    # only the part of the network connected to a takes part
    comp, todo = {a}, [a]
    while todo:
        u = todo.pop()
        for v in nbr[u]:
            if v not in comp:
                comp.add(v)
                todo.append(v)
    nbr = {i: ({j: s for j, s in nbr[i].items() if j in comp} if i in comp else {}) for i in range(n)}
    # committor q (0 at a, 1 at b): sum_j k_ij (q_j - q_i) = 0, each row scaled by its fastest rate
    inner = [i for i in range(n) if i not in (a, b) and nbr[i]]
    q = np.zeros(n)
    q[b] = 1.0
    if inner:
        idx = {i: k for k, i in enumerate(inner)}
        M = np.zeros((len(inner), len(inner)))
        rhs = np.zeros(len(inner))
        for i in inner:
            smin = min(nbr[i].values())
            r = idx[i]
            for j, s in nbr[i].items():
                w = math.exp(-beta * (s - smin))
                M[r, r] -= w
                if j == b:
                    rhs[r] -= w
                elif j != a:
                    M[r, idx[j]] += w
        sol = np.linalg.lstsq(M, rhs, rcond=None)[0]
        for i in inner:
            q[i] = min(1.0, max(0.0, sol[idx[i]]))
    # reactive flux f_ij = pi_i k_ij (1 - q_i) q_j, in units of k0/Z * exp(-beta*S_ref)
    s_ref = min(nbr[i][j] for i in range(n) for j in nbr[i]) if any(nbr.values()) else 0.0
    net = {}
    for i in range(n):
        for j, s in nbr[i].items():
            f = math.exp(-beta * (s - s_ref)) * (1 - q[i]) * q[j]
            g = math.exp(-beta * (s - s_ref)) * (1 - q[j]) * q[i]
            if f > g:
                net[(i, j)] = f - g
    total = sum(v for (i, _), v in net.items() if i == a)
    # rate from a: k_AB = F / sum_i pi_i (1 - q_i)
    emin = float(E.min())
    denom = sum(math.exp(-beta * (E[i] - emin)) * (1 - q[i]) for i in range(n) if nbr[i] or i in (a, b))
    k_ab = k0 * total * math.exp(-beta * (s_ref - emin)) / denom if denom > 0 and total > 0 else 0.0
    # pathway decomposition: repeatedly take the widest path in the net flux graph and subtract its bottleneck
    paths = []
    rem = dict(net)
    for _ in range(max_paths):
        width = {a: math.inf}
        prev = {}
        heap = [(-math.inf, a)]
        while heap:
            w, u = heapq.heappop(heap)
            w = -w
            if w < width.get(u, 0):
                continue
            for (x, y), f in rem.items():
                if x != u or f <= 0:
                    continue
                c = min(w, f)
                if c > width.get(y, 0):
                    width[y] = c
                    prev[y] = u
                    heapq.heappush(heap, (-c, y))
        if b not in width or total <= 0 or width[b] < 0.01 * total:
            break
        route = [b]
        while route[-1] != a:
            route.append(prev[route[-1]])
        route = route[::-1]
        c = width[b]
        for x, y in zip(route, route[1:]):
            rem[(x, y)] -= c
        paths.append({'nodes': [int(x) for x in route], 'flux_share': round(float(c / total), 4)})
    # mean first-passage time from a to b: sum_j k_ij (tau_j - tau_i) = -1 for i != b, tau_b = 0 (rows scaled)
    mfpt = occupancy = None
    others = [i for i in range(n) if i != b and (nbr[i] or i == a)]
    if a in others and nbr[a]:
        idx = {i: k for k, i in enumerate(others)}
        M = np.zeros((len(others), len(others)))
        rhs = np.zeros(len(others))
        for i in others:
            if not nbr[i]:
                M[idx[i], idx[i]] = 1.0
                continue
            lograte = {j: -beta * (sv - E[i]) for j, sv in nbr[i].items()}
            top = max(lograte.values())
            r = idx[i]
            for j, lr in lograte.items():
                wgt = math.exp(lr - top)
                M[r, r] -= wgt
                if j != b:
                    M[r, idx[j]] += wgt
            rhs[r] = -math.exp(-top) / k0            # -1 / (k0 * exp(top))
        tau = np.linalg.lstsq(M, rhs, rcond=None)[0]
        if np.isfinite(tau[idx[a]]) and tau[idx[a]] > 0:
            mfpt = float(tau[idx[a]])
            # where that time is spent: expected time in each node before b is reached, x = e_a^T (-Q)^-1. With
            # Q = S M (S = fastest rate out of each node, the row scaling above): M^T y = -e_a and x = y / S.
            e = np.zeros(len(others))
            e[idx[a]] = -1.0
            y = np.linalg.lstsq(M.T, e, rcond=None)[0]
            occ = np.zeros(n)
            for i in others:
                if nbr[i]:
                    top = max(-beta * (sv - E[i]) for sv in nbr[i].values())
                    occ[i] = max(0.0, y[idx[i]]) * math.exp(-top) / k0
            if occ.sum() > 0:
                occupancy = [round(float(v / occ.sum()), 4) for v in occ]
    return {'committor': [round(float(x), 4) for x in q], 'rate': float(k_ab), 'mfpt': mfpt, 'occupancy': occupancy, 'flux_total_relative': float(total),
            'pathways': paths, 'net_flux': {f'{i}-{j}': round(float(v / total), 5) for (i, j), v in net.items()} if total > 0 else {}}


# ----------------------------------------------------------------------------- the 2D analysis
def profile_of(land, route_structs):
    """Concatenated refolding path along a route: list of {'s', 'E', 'leg'}."""
    out = []
    for leg, (x, y) in enumerate(zip(route_structs, route_structs[1:])):
        steps = land.path(x, y)
        for k, (s, e) in enumerate(steps):
            if out and k == 0:
                continue
            out.append({'s': s, 'E': e, 'leg': leg})
    if not out:
        out.append({'s': route_structs[0], 'E': land.E(route_structs[0]), 'leg': 0})
    return out


def events_of(profile):
    """Base-pair events between consecutive structures of a profile."""
    ev = []
    m = max(1, len(profile) - 1)
    for k in range(1, len(profile)):
        p0, p1 = pair_set(profile[k - 1]['s']), pair_set(profile[k]['s'])
        for i, j in sorted(p0 - p1):
            ev.append({'step': k, 't': k / m, 'type': 'open', 'i': i, 'j': j, 'E': profile[k]['E']})
        for i, j in sorted(p1 - p0):
            ev.append({'step': k, 't': k / m, 'type': 'close', 'i': i, 'j': j, 'E': profile[k]['E']})
    return ev


def bottleneck_of(profile):
    """Highest saddle, the cooperative unit opened on the climb to it and the nucleus formed after it."""
    E = [p['E'] for p in profile]
    top = int(np.argmax(E))
    leg = profile[top]['leg']
    lo = min(k for k, p in enumerate(profile) if p['leg'] == leg) - (1 if leg > 0 else 0)
    lo = max(0, lo)
    hi = max(k for k, p in enumerate(profile) if p['leg'] == leg)
    before = range(lo, top + 1)
    m = min(before, key=lambda k: (E[k], -k))          # lowest point before the saddle (the last of equals)
    after = range(top, hi + 1)
    nxt = min(after, key=lambda k: (E[k], k))
    p_m, p_top, p_n = pair_set(profile[m]['s']), pair_set(profile[top]['s']), pair_set(profile[nxt]['s'])
    unit = sorted(p_m - p_top)
    nucleus = sorted(p_n - p_top)
    return {'saddle_step': top, 'saddle_structure': profile[top]['s'], 'saddle_energy': E[top],
            'climb_from_step': m, 'climb_energy': round(E[top] - E[m], 2), 'unit': [list(p) for p in unit],
            'unit_helices': [[list(p) for p in h] for h in helix_segments(unit)],
            'nucleus': [list(p) for p in nucleus], 'nucleus_helices': [[list(p) for p in h] for h in helix_segments(nucleus)],
            'descend_to_step': nxt}


def residue_timing(n, events, unit, nucleus):
    """Per residue: first and last time its pairing changes on the route (None if never), and its role."""
    first, last = [None] * n, [None] * n
    for e in events:
        for r in (e['i'], e['j']):
            first[r] = e['t'] if first[r] is None else min(first[r], e['t'])
            last[r] = e['t'] if last[r] is None else max(last[r], e['t'])
    role = [''] * n
    for i, j in unit:
        role[i] = role[j] = 'unit'
    for i, j in nucleus:
        if not role[i]:
            role[i] = 'nucleus'
        if not role[j]:
            role[j] = 'nucleus'
    return [{'first': first[r], 'last': last[r], 'role': role[r]} for r in range(n)]


def find_trap(energies, committor, occupancy=None, a=0, b=1, margin=0.5, share=0.5):
    """A trap: a minimum deeper than A's basin that has to be left again before B is reached, which slows the switch.

    With the occupancy (share of the mean first-passage time spent in each node): when the nodes below A together
    take at least `share` of the time, the trap is the one that holds the RNA longest (near-degenerate wells split
    the time between them; 'nodes' lists all that take 5 % or more). Its kind: 'a-side' (committor < 0.5), 'deep'
    (deeper than B too) or 'intermediate' (on the way to B). Without the occupancy: the deepest node below A on the
    A side, or deeper than both end basins with 0.02 < committor < 0.98. None when there is no such node."""
    eMa, eMb = energies[a], energies[b]

    def kind_of(i):
        return 'a-side' if committor[i] < 0.5 else 'deep' if energies[i] < eMb - margin else 'intermediate'
    if occupancy is not None:
        below = [i for i in range(len(energies)) if i not in (a, b) and energies[i] < eMa - margin]
        total = sum(occupancy[i] for i in below)
        if not below or total < share:
            return None
        i = max(below, key=lambda k: (occupancy[k], -energies[k]))
        return {'node': i, 'energy': energies[i], 'depth': round(eMa - energies[i], 2), 'kind': kind_of(i),
                'time_share': round(total, 4), 'nodes': sorted((k for k in below if occupancy[k] >= 0.05), key=lambda k: -occupancy[k])}
    cand = []
    for i, q in enumerate(committor):
        if i in (a, b) or q is None:
            continue
        if q < 0.5 and energies[i] < eMa - margin:
            cand.append((energies[i], i, 'a-side'))
        elif energies[i] < min(eMa, eMb) - margin and 0.02 < q < 0.98:
            cand.append((energies[i], i, 'deep'))
    if not cand:
        return None
    e, i, kind = min(cand)
    return {'node': i, 'energy': e, 'depth': round(eMa - e, 2), 'kind': kind}


def analyse_2d(seq, sa, sb, cond: Conditions, states=(), k0=K0_DEFAULT, intermediates=True, n_samples=None,
               max_nodes=36, max_edges=360):
    """Refolding from 2D structure sa to sb. Kinetics are computed between basins - the local minima that sa and sb
    reach by steepest descent - as in barrier-tree analyses: a structure that is not itself a minimum relaxes into
    its basin first, and the switch is the crossing from A's basin to B's. The direct path between the two
    structures themselves is reported for comparison."""
    land = Landscape(seq, cond)
    n = len(seq)
    sa, sb = compatible(seq, sa), compatible(seq, sb)
    eA, eB = land.E(sa), land.E(sb)
    direct = profile_of(land, [sa, sb])
    s_direct = max(p['E'] for p in direct)
    down_a, down_b = land.descend(sa, path=True), land.descend(sb, path=True)
    ma, mb = down_a[-1], down_b[-1]
    eMa, eMb = land.E(ma), land.E(mb)
    same_basin = ma == mb
    dab = bp_distance(sa, sb)
    # --- network: the two basins, the basins of the ensemble's other states and sampled local minima
    nodes, kinds, labels, members = [ma, mb], ['A', 'B'], ['A', 'B'], [[sa], [sb]]
    s_ab = land.saddle(ma, mb) if not same_basin else eMa
    if intermediates and not same_basin:
        cand = []
        for st in states:
            s0 = compatible(seq, st['dbn'])
            cand.append((land.descend(s0), 'ensemble', st.get('label'), s0))
        if n_samples is None:
            n_samples = int(min(150, max(20, 4000 // max(1, n))))
        for m in land.sample_minima(n_samples):
            cand.append((m, 'minimum', None, None))
        E = [p['E'] for p in direct]
        for k in range(1, len(direct) - 1):          # minima along the direct path
            if E[k] < E[k - 1] and E[k] <= E[k + 1]:
                cand.append((land.descend(direct[k]['s']), 'minimum', None, None))
        found = {}
        for m, kind, lab, src in cand:
            if m in (ma, mb):
                if src is not None and src not in members[nodes.index(m)]:
                    members[nodes.index(m)].append(src)
                continue
            f = found.setdefault(m, {'kind': kind, 'label': lab, 'members': []})
            if kind == 'ensemble':
                f['kind'], f['label'] = 'ensemble', f['label'] if f['kind'] == 'ensemble' else lab
            if src is not None:
                f['members'].append(src)
        keep = []
        for m, f in found.items():
            e = land.E(m)
            if e >= s_ab - 0.3:
                continue                              # cannot lower the barrier between the basins
            if bp_distance(ma, m) + bp_distance(m, mb) > 1.6 * bp_distance(ma, mb) + 6:
                continue                              # far off the way between A and B
            keep.append((e if f['kind'] != 'ensemble' else e - 100, m, f))
        keep.sort(key=lambda x: (x[0], x[1]))
        for _, m, f in keep[:max(0, max_nodes - 2)]:
            nodes.append(m)
            kinds.append(f['kind'])
            labels.append(f['label'])
            members.append(f['members'])
    energies = [land.E(m) for m in nodes]
    # --- edges: findpath saddles between nearby nodes (always A-B)
    edges = {(0, 1): s_ab} if not same_basin else {}
    pairs = sorted((bp_distance(nodes[i], nodes[j]), i, j) for i in range(len(nodes)) for j in range(i + 1, len(nodes)) if (i, j) != (0, 1))
    limit = max(8, int(0.75 * bp_distance(ma, mb)) + 4)
    for d, i, j in pairs[:max_edges]:
        if d <= limit:
            edges[(i, j)] = land.saddle(nodes[i], nodes[j])
    # --- best route between the basins, and the full profile A -> basin(A) -> ... -> basin(B) -> B
    route = (_minimax_route(nodes, edges, 0, 1) or [0, 1]) if not same_basin else [0]
    core = profile_of(land, [nodes[k] for k in route])
    relax_a = [{'s': x, 'E': land.E(x), 'leg': -1, 'phase': 'relax_a'} for x in down_a[:-1]]
    relax_b = [{'s': x, 'E': land.E(x), 'leg': len(route) - 1, 'phase': 'relax_b'} for x in down_b[-2::-1]]
    for p_ in core:
        p_['phase'] = 'route'
    full = relax_a + core + relax_b
    kT = cond.kT
    tpt = transition_path_theory(energies, edges, 0, 1, kT, k0) if not same_basin else None
    s_best = max(p_['E'] for p_ in core)
    bott = None
    if not same_basin:
        bott = bottleneck_of(core)
        bott['saddle_step'] += len(relax_a)
        bott['climb_from_step'] += len(relax_a)
        bott['descend_to_step'] += len(relax_a)
    ev = events_of(full)
    timing = residue_timing(n, ev, bott['unit'] if bott else [], bott['nucleus'] if bott else [])
    # time to switch: mean first-passage time from A's basin to B's on the network (falls back to the TPT rate, then
    # to the Arrhenius rate of the best route)
    if same_basin:
        rate = None
    elif tpt and tpt.get('mfpt'):
        rate = 1.0 / tpt['mfpt']
    elif tpt and tpt['rate'] > 0:
        rate = tpt['rate']
    else:
        rate = k0 * math.exp(-(s_best - eMa) / kT)
    k_direct = k0 * math.exp(-(s_direct - eA) / kT) if sa != sb else None
    trap = find_trap(energies, tpt['committor'], tpt.get('occupancy')) if tpt else None
    if trap and rate and rate > 0.5 * k0 * math.exp(-(s_best - eMa) / kT):
        trap = None                                   # the RNA waits there, but the switch is not slower for it

    def summary(profile, start):
        e = [p_['E'] for p_ in profile]
        return {'barrier': round(max(e) - start, 2), 'saddle': round(max(e), 2), 'steps': len(profile) - 1,
                'profile': [{'s': p_['s'], 'E': p_['E'], 'leg': p_['leg'], **({'phase': p_['phase']} if 'phase' in p_ else {})} for p_ in profile]}
    best = {**summary(full, eMa), 'route': route, 'route_barrier': round(s_best - eMa, 2)}
    out = {
        'sequence': seq, 'conditions': cond.to_json(), 'conditions_text': cond.describe(), 'kT': kT, 'k0': k0,
        'A': {'dbn': sa, 'energy': eA, 'basin': ma, 'basin_energy': eMa, 'above_basin': round(eA - eMa, 2)},
        'B': {'dbn': sb, 'energy': eB, 'basin': mb, 'basin_energy': eMb, 'above_basin': round(eB - eMb, 2)},
        'same_basin': same_basin, 'dG': round(eB - eA, 2), 'dG_basins': round(eMb - eMa, 2), 'bp_distance': dab,
        'direct': summary(direct, eA), 'best': best,
        'nodes': [{'dbn': m, 'energy': e, 'kind': k, 'label': lab, 'members': mem} for m, e, k, lab, mem in zip(nodes, energies, kinds, labels, members)],
        'edges': [{'a': i, 'b': j, 'saddle': round(sv, 2)} for (i, j), sv in edges.items()],
        'tpt': tpt, 'trap': trap, 'bottleneck': bott, 'events': ev, 'timing': timing,
        'rate': rate, 'rate_direct': k_direct,
        'half_life_s': (math.log(2) / rate) if rate else None,
        'effective_barrier': (round(-kT * math.log(rate / k0), 2) if rate else None),
        'method': {'path': f'findpath (Flamm et al. 2001), search width {land.width}', 'width': land.width, 'intermediates': bool(intermediates),
                   'samples': n_samples, 'seed': SAMPLE_SEED, 'vienna': _vienna(), 'kinetics': 'between basins (steepest-descent minima)'},
    }
    return out


def _vienna():
    try:
        import RNA
        return RNA.__version__
    except Exception:
        return None


def phi_scan(seq, route_structs, cond: Conditions, width=None, progress=None):
    """In-silico phi-value scan over every point mutation, with the route kept and its energies re-evaluated."""
    width = width or (6 if len(seq) <= 60 else 4 if len(seq) <= 120 else 2)
    def barrier(land, route):
        e0 = land.E(route[0])
        s = max(land.saddle(x, y) for x, y in zip(route, route[1:])) if len(route) > 1 else e0
        return s - e0, land.E(route[-1]) - e0
    wt = Landscape(seq, cond, width=width)
    b_wt, g_wt = barrier(wt, [compatible(seq, s) for s in route_structs])
    rows = []
    for i, base0 in enumerate(seq):
        row = {}
        for base in 'ACGU':
            if base == base0:
                continue
            m = seq[:i] + base + seq[i + 1:]
            land = Landscape(m, cond, width=width)
            route = [compatible(m, s) for s in route_structs]
            b, g = barrier(land, route)
            dts, deq = b - b_wt, g - g_wt
            row[base] = {'ddg_ts': round(dts, 2), 'ddg_eq': round(deq, 2),
                         'phi': round(dts / deq, 2) if abs(deq) >= PHI_MIN_DDG else None,
                         'breaks': any(r != s for r, s in zip(route, [compatible(seq, x) for x in route_structs]))}
        rows.append(row)
        if progress:
            progress(i + 1, len(seq))
    per_pos = []
    for row in rows:
        vals = [v['ddg_ts'] for v in row.values()]
        phis = [v['phi'] for v in row.values() if v['phi'] is not None]
        per_pos.append({'max_abs_ddg_ts': round(max(abs(x) for x in vals), 2), 'min_ddg_ts': min(vals), 'max_ddg_ts': max(vals),
                        'phi': round(float(np.median(phis)), 2) if phis else None, 'n_phi': len(phis)})
    flat = [(i, b, v) for i, row in enumerate(rows) for b, v in row.items()]
    faster = sorted((x for x in flat if x[2]['ddg_ts'] < -0.05), key=lambda x: x[2]['ddg_ts'])[:8]
    slower = sorted((x for x in flat if x[2]['ddg_ts'] > 0.05), key=lambda x: -x[2]['ddg_ts'])[:8]
    fmt = lambda x: {'pos': x[0], 'from': seq[x[0]], 'to': x[1], **x[2]}
    return {'barrier_wt': round(b_wt, 2), 'dg_wt': round(g_wt, 2), 'matrix': rows, 'positions': per_pos,
            'faster': [fmt(x) for x in faster], 'slower': [fmt(x) for x in slower], 'width': width,
            'phi_min_ddg': PHI_MIN_DDG}


# ----------------------------------------------------------------------------- 3D: rigid modules
def _kabsch(P, Q):
    """R, t minimising |P R + t - Q|; returns (R, t, rmsd)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1.0, 1.0, d])
    R = U @ D @ Vt
    t = qc - pc @ R
    rmsd = float(np.sqrt(np.mean(np.sum((P @ R + t - Q) ** 2, axis=1))))
    return R, t, rmsd


def _rotation_angle(R):
    return math.degrees(math.acos(max(-1.0, min(1.0, (np.trace(R) - 1) / 2))))


def _quaternion(R):
    """Unit quaternion (w, x, y, z) of the rotation x -> x R (row vectors), i.e. of the matrix R.T acting on columns."""
    M = R.T
    w = math.sqrt(max(0.0, 1 + M[0, 0] + M[1, 1] + M[2, 2])) / 2
    x = math.sqrt(max(0.0, 1 + M[0, 0] - M[1, 1] - M[2, 2])) / 2
    y = math.sqrt(max(0.0, 1 - M[0, 0] + M[1, 1] - M[2, 2])) / 2
    z = math.sqrt(max(0.0, 1 - M[0, 0] - M[1, 1] + M[2, 2])) / 2
    x = math.copysign(x, M[2, 1] - M[1, 2])
    y = math.copysign(y, M[0, 2] - M[2, 0])
    z = math.copysign(z, M[1, 0] - M[0, 1])
    nrm = math.sqrt(w * w + x * x + y * y + z * z) or 1.0
    return [w / nrm, x / nrm, y / nrm, z / nrm]


def _screw(R, t, centre):
    """Screw axis of x -> x R + t: axis direction, point on the axis near centre, angle (deg), shift along axis (A)."""
    M = R.T
    ang = _rotation_angle(R)
    axis = np.array([M[2, 1] - M[1, 2], M[0, 2] - M[2, 0], M[1, 0] - M[0, 1]])
    if np.linalg.norm(axis) < 1e-8:
        return {'angle': round(ang, 1), 'axis': None, 'point': None, 'shift': round(float(np.linalg.norm(t)), 2)}
    axis /= np.linalg.norm(axis)
    # column form x' = M x + t; points p on the axis satisfy M p + t = p + s u with s = t.u (Chasles)
    shift = float(t @ axis)
    p = np.linalg.lstsq(np.eye(3) - M, t - shift * axis, rcond=None)[0]
    point = p + axis * float((centre - p) @ axis)          # the point of the axis nearest to the module centre
    return {'angle': round(ang, 1), 'axis': [round(float(v), 4) for v in axis],
            'point': [round(float(v), 2) for v in point], 'shift': round(shift, 2)}


def residue_atoms(doc):
    """Per residue the indices of the backbone/sugar atoms used for rigid fits (C4', C1', P; O3' if few)."""
    names = doc['atoms']['names']
    out = []
    for r in doc['residues']:
        idx = [r['first_atom'] + q for q in range(r['n_atoms']) if names[r['first_atom'] + q] in ("C4'", "C1'", 'P', "O3'")]
        out.append(idx)
    return out


def _pseudo_torsions(X, trace, phos):
    """eta (C4'(i-1)-P(i)-C4'(i)-P(i+1)) and theta (P(i)-C4'(i)-P(i+1)-C4'(i+1)) per residue, degrees or None."""
    def dih(a, b, c, d):
        b0, b1, b2 = a - b, c - b, d - c
        b1 /= np.linalg.norm(b1) or 1.0
        v = b0 - (b0 @ b1) * b1
        w = b2 - (b2 @ b1) * b1
        return math.degrees(math.atan2(float(np.cross(b1, v) @ w), float(v @ w)))
    n = len(trace)
    eta, theta = [None] * n, [None] * n
    for i in range(1, n - 1):
        if min(trace[i - 1], phos[i], trace[i], phos[i + 1]) >= 0:
            eta[i] = dih(X[trace[i - 1]], X[phos[i]], X[trace[i]], X[phos[i + 1]])
        if min(phos[i], trace[i], phos[i + 1], trace[i + 1]) >= 0:
            theta[i] = dih(X[phos[i]], X[trace[i]], X[phos[i + 1]], X[trace[i + 1]])
    return eta, theta


def rigid_modules(doc, XA, XB, pairs_a, pairs_b, timing=None, tau=RIGID_TAU, min_size=4):
    """Rigid modules between two conformations (coordinates XA, XB: atoms x 3, same atom table)."""
    n = len(doc['residues'])
    ratoms = residue_atoms(doc)
    names = doc['atoms']['names']
    trace = [next((r['first_atom'] + q for q in range(r['n_atoms']) if names[r['first_atom'] + q] == "C4'"), -1) for r in doc['residues']]
    phos = [next((r['first_atom'] + q for q in range(r['n_atoms']) if names[r['first_atom'] + q] == 'P'), -1) for r in doc['residues']]

    def atoms_of(group):
        return [a for r in group for a in ratoms[r]]

    def fit(group):
        idx = atoms_of(group)
        if len(idx) < 3:
            return None
        return _kabsch(XA[idx], XB[idx])

    # seeds: helices kept in both conformations (both strands), then short single-stranded chunks
    shared = set(map(tuple, pairs_a)) & set(map(tuple, pairs_b))
    groups, owner = [], [-1] * n
    for h in helix_segments(shared):
        if len(h) < 2:
            continue
        g = sorted({r for p in h for r in p})
        if any(owner[r] >= 0 for r in g):
            continue
        for r in g:
            owner[r] = len(groups)
        groups.append(g)
    run = []
    for r in list(range(n)) + [None]:
        if r is not None and owner[r] < 0:
            run.append(r)
            if len(run) == 3:
                for x in run:
                    owner[x] = len(groups)
                groups.append(run)
                run = []
            continue
        if run:
            for x in run:
                owner[x] = len(groups)
            groups.append(run)
            run = []
    # agglomerative merging of neighbouring groups while the union stays rigid
    c4 = np.array([XA[trace[r]] if trace[r] >= 0 else XA[ratoms[r][0]] for r in range(n)])
    c4b = np.array([XB[trace[r]] if trace[r] >= 0 else XB[ratoms[r][0]] for r in range(n)])
    near = (np.linalg.norm(c4[:, None] - c4[None], axis=2) < 7.5) | (np.linalg.norm(c4b[:, None] - c4b[None], axis=2) < 7.5)
    for r in range(n - 1):
        near[r, r + 1] = near[r + 1, r] = True
    alive = {k: g for k, g in enumerate(groups)}
    cost = {}

    def touching(g1, g2):
        return bool(near[np.ix_(g1, g2)].any())

    def merge_cost(k1, k2):
        f = fit(alive[k1] + alive[k2])
        return f[2] if f else math.inf

    keys = sorted(alive)
    for x in range(len(keys)):
        for y in range(x + 1, len(keys)):
            if touching(alive[keys[x]], alive[keys[y]]):
                cost[(keys[x], keys[y])] = merge_cost(keys[x], keys[y])
    while cost:
        (k1, k2), c = min(cost.items(), key=lambda kv: kv[1])
        if c > tau:
            break
        alive[k1] = sorted(alive[k1] + alive.pop(k2))
        cost = {k: v for k, v in cost.items() if k1 not in k and k2 not in k}
        for k in alive:
            if k != k1 and touching(alive[k1], alive[k]):
                cost[(min(k, k1), max(k, k1))] = merge_cost(k, k1)
    mods = sorted((g for g in alive.values() if len(g) >= min_size and (fit(g) or (0, 0, math.inf))[2] <= tau), key=lambda g: (-len(g), g[0]))
    module_of = [-1] * n
    for m, g in enumerate(mods):
        for r in g:
            module_of[r] = m
    # motions relative to the largest module (B superposed on A by that module)
    out_mods = []
    if mods:
        ref = atoms_of(mods[0])
        Rr, tr, _ = _kabsch(XB[ref], XA[ref])
        XBr = XB @ Rr + tr
    for m, g in enumerate(mods):
        idx = atoms_of(g)
        R, t, rmsd = _kabsch(XA[idx], XB[idx])                     # display frame (for the morph)
        R2, t2, _ = _kabsch(XA[idx], XBr[idx])                     # relative to the largest module
        ca, cb = XA[idx].mean(axis=0), XB[idx].mean(axis=0)
        ev = [timing[r] for r in g] if timing else []
        firsts = [e['first'] for e in ev if e['first'] is not None]
        lasts = [e['last'] for e in ev if e['last'] is not None]
        out_mods.append({'id': m, 'residues': g, 'size': len(g), 'fit_rmsd': round(rmsd, 2),
                         'rotation': round(_rotation_angle(R2), 1), 'translation': round(float(np.linalg.norm(XBr[idx].mean(axis=0) - ca)), 2),
                         'screw': _screw(R2, t2, ca), 'quaternion': [round(v, 6) for v in _quaternion(R)],
                         'centre_a': [round(float(v), 3) for v in ca], 'centre_b': [round(float(v), 3) for v in cb],
                         'window': [min(firsts), max(lasts)] if firsts else None})
    # hinges: flexible residues next to modules, ranked by backbone pseudo-torsion change
    etaA, thA = _pseudo_torsions(XA, trace, phos)
    etaB, thB = _pseudo_torsions(XB, trace, phos)

    def dang(x, y):
        if x is None or y is None:
            return 0.0
        d = abs(x - y) % 360
        return min(d, 360 - d)
    change = [dang(etaA[r], etaB[r]) + dang(thA[r], thB[r]) for r in range(n)]
    flexible = [r for r in range(n) if module_of[r] < 0]
    hinges = sorted((r for r in flexible if any(0 <= r + d < n and module_of[r + d] >= 0 for d in (-2, -1, 1, 2))),
                    key=lambda r: -change[r])
    # blend weights for the morph: flexible residue -> neighbouring modules along the chain
    blend = []
    for r in range(n):
        if module_of[r] >= 0:
            blend.append([[module_of[r], 1.0]])
            continue
        up = next(((r - k, module_of[r - k]) for k in range(1, n) if r - k >= 0 and module_of[r - k] >= 0), None)
        dn = next(((r + k, module_of[r + k]) for k in range(1, n) if r + k < n and module_of[r + k] >= 0), None)
        if up and dn and up[1] != dn[1]:
            u = (r - up[0]) / (dn[0] - up[0])
            blend.append([[up[1], round(1 - u, 4)], [dn[1], round(u, 4)]])
        elif up or dn:
            blend.append([[(up or dn)[1], 1.0]])
        else:
            blend.append([])
    return {'modules': out_mods, 'module_of': module_of, 'hinges': hinges[:12],
            'torsion_change': [round(c, 1) for c in change], 'blend': blend, 'tau': tau,
            'method': "rigid groups seeded by helices shared by A and B, merged while the joint fit of C4', C1', P "
                      f"and O3' atoms stays below {tau:g} Å; motions relative to the largest module"}


def couple(modules, bottleneck, timing):
    """Which module motion depends on the cooperative unit (the rate-limiting coordinated motion)."""
    if not bottleneck or not modules['modules']:
        return None
    unit = {r for p in bottleneck['unit'] for r in p}
    nuc = {r for p in bottleneck['nucleus'] for r in p}
    best = None
    for m in modules['modules']:
        res = set(m['residues'])
        # residues of the module and the flexible residues next to it
        ext = set(res)
        for r in res:
            for d in (-2, -1, 1, 2):
                if 0 <= r + d < len(modules['module_of']) and modules['module_of'][r + d] < 0:
                    ext.add(r + d)
        overlap = len(ext & unit) + 0.5 * len(ext & nuc)
        score = overlap * (1 + m['rotation'] / 30)
        if overlap and (best is None or score > best[0]):
            best = (score, m['id'], sorted(ext & unit), sorted(ext & nuc))
    if not best:
        return None
    return {'module': best[1], 'unit_residues': best[2], 'nucleus_residues': best[3]}


# ----------------------------------------------------------------------------- entry point
def state_dbn(doc, k):
    for st in ensemble_states(doc):
        if k in st['conformations']:
            return st['dbn']
    raise ValueError(f'conformation {k} not found')


def analyse(doc, a, b, cond: Conditions, k0=K0_DEFAULT, intermediates=True, coords=None):
    seq = doc['sequence']
    n_conf = len(doc['conformations'])
    if not (0 <= a < n_conf and 0 <= b < n_conf):
        raise ValueError('conformation index out of range')
    states = ensemble_states(doc)
    for i, st in enumerate(states):
        st['label'] = f'S{i + 1}'
    sa, sb = state_dbn(doc, a), state_dbn(doc, b)
    res = analyse_2d(seq, sa, sb, cond, states=states, k0=k0, intermediates=intermediates)
    res['a'], res['b'] = a, b
    # network nodes (basins) -> the conformations of the ensemble states that relax into them (closest to A first)
    rmsd = doc.get('ensemble', {}).get('rmsd_c4')
    for idx, node in enumerate(res['nodes']):
        mem = set(node.get('members') or [])
        confs = sorted({k for st in states if compatible(seq, st['dbn']) in mem for k in st['conformations']})
        if idx == 0:
            confs = [a] + [k for k in confs if k != a]
        elif idx == 1:
            confs = [b] + [k for k in confs if k != b]
        elif rmsd and confs:
            confs = sorted(confs, key=lambda k: rmsd[a][k])
        node['conformations'] = confs
    if coords is not None:
        XA, XB = np.asarray(coords[a], float), np.asarray(coords[b], float)
        pa = [p[:2] for p in doc['conformations'][a]['pairs']]
        pb = [p[:2] for p in doc['conformations'][b]['pairs']]
        mods = rigid_modules(doc, XA, XB, pa, pb, res['timing'])
        res['motion'] = mods
        res['coupling'] = couple(mods, res['bottleneck'], res['timing'])
        res['rmsd_ab'] = round(float(rmsd[a][b]), 2) if rmsd else None
    return res
