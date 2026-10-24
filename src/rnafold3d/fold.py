"""v2 fold-tree model: rigid all-atom A-form stems connected by torsion-space linkers.

Topology
--------
* Every stem (run of stacked pairs) is a rigid body built by build.build_stem.
* Linkers are the backbone segments between consecutive stem strands. A minimum
  spanning tree over inter-stem linkers (weight = number of unpaired residues) is
  chosen; tree linkers are built kinematically (their torsions position the child
  stem), the remaining linkers and all hairpin loops are *closure* linkers that are
  closed afterwards by cyclic coordinate descent (loops.py).
* 5'/3' tails hang from the first/last stem.

Moves are rotations about backbone bonds (alpha..zeta) or glycosidic bonds (chi), so
covalent geometry stays exact during the search; torsion values are tracked and the
final model is rebuilt from them with NeRF.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np

from .build import (ResidueBuilder, load_geometry, build_stem, atom_order, BASE_ATOMS, GLY, SUGAR_ATOMS)
from .geomutil import kabsch, rotation_about_axis, dihedral

TOR = {'alpha': 0, 'beta': 1, 'gamma': 2, 'epsilon': 3, 'zeta': 4, 'chi': 5}
LIB_COLS = ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi')


@dataclass
class Body:
    id: int
    kind: str                 # 'stem' or 'residue'
    residues: list
    pairs: list = field(default_factory=list)
    atoms: np.ndarray = None  # global atom indices owned by the body
    local: np.ndarray = None  # local coordinates of those atoms


@dataclass
class Linker:
    id: int
    a: int | None             # upstream residue (strand end) or None (5' tail)
    b: int | None             # downstream residue (strand start) or None (3' tail)
    loop: tuple
    kind: str                 # tree | closure | hairpin | tail5 | tail3
    body_a: int | None = None
    body_b: int | None = None
    direction: str = 'forward'
    parent: int | None = None
    child: int | None = None
    atoms: np.ndarray = None  # owned atom indices

    @property
    def n(self):
        return len(self.loop)


@dataclass
class Dof:
    res: int
    name: str
    u: int                    # axis start atom index
    v: int                    # axis end atom index (rotation centre)
    moved: np.ndarray
    linker: int


def stems_from_pairs(pairs):
    pairset = set(pairs)
    starts = sorted(p for p in pairs if (p[0] - 1, p[1] + 1) not in pairset)
    stems = []
    for i, j in starts:
        run = []
        while (i, j) in pairset:
            run.append((i, j))
            i += 1
            j -= 1
        stems.append(run)
    return stems


def check_pairs(n, pairs, min_hairpin=3):
    seen = set()
    for i, j in pairs:
        if not 0 <= i < j < n:
            raise ValueError(f'pair out of range: {(i + 1, j + 1)}')
        if i in seen or j in seen:
            raise ValueError(f'residue paired twice: {(i + 1, j + 1)}')
        seen.update((i, j))
        if j - i - 1 < min_hairpin:
            # only an error if (i, j) closes a hairpin (no pair nested inside)
            if not any(i < k < l < j for k, l in pairs):
                raise ValueError(f'hairpin loop shorter than {min_hairpin} nt at {(i + 1, j + 1)}')
    for i, j in pairs:
        for k, l in pairs:
            if i < k < j < l:
                raise ValueError('crossing (pseudoknotted) pairs are not supported by the v2 fold tree')


class FoldModel:
    """Kinematic all-atom model for one sequence and one nested secondary structure."""

    def __init__(self, sequence, pairs, rb: ResidueBuilder | None = None, explicit_loops: bool = True):
        self.explicit_loops = explicit_loops
        self.seq = sequence
        self.n = len(sequence)
        self.pairs = sorted((min(i, j), max(i, j)) for i, j in pairs)
        check_pairs(self.n, self.pairs)
        self.rb = rb or ResidueBuilder(load_geometry())
        self.geom = self.rb.g
        self._index_atoms()
        self._make_bodies()
        self._make_linkers()
        self._make_tree()
        self._assign_atoms()
        self._make_dofs()
        self._make_exclusions()
        self.tor = np.zeros((self.n, 6))
        aform = self.geom.aform
        for k, name in TOR.items():
            self.tor[:, name] = aform[k]
        self.puck = np.array(['C3endo'] * self.n, dtype=object)
        self.xyz = np.zeros((self.n_atoms, 3))
        self.present = np.zeros(self.n_atoms, bool)
        self.set_loops_active(explicit_loops)

    def set_loops_active(self, flag):
        """All loops (closure/hairpin linkers) explicit in search (True) or absent (False)."""
        self.set_active_loops({L.id for L in self.loop_linkers} if flag else set())

    def set_active_loops(self, ids):
        """Make the given closure/hairpin linkers explicit (built, sampled, chain-break restrained)."""
        self.active_loops = set(ids) if self.explicit_loops else set()
        self.loops_active = bool(self.active_loops)
        loop_ids = {L.id for L in self.loop_linkers}
        self.active_dofs = [k for k, d in enumerate(self.dofs) if d.linker not in loop_ids or d.linker in self.active_loops]
        # loops whose closure is perturbed by each DOF (anchors on opposite sides of the moved set,
        # or the DOF lies inside the loop itself)
        if not hasattr(self, '_moved_sets'):
            self._moved_sets = [set(d.moved.tolist()) for d in self.dofs]
        self.dof_loops = {}
        active = [L for L in self.loop_linkers if L.id in self.active_loops]
        for k in self.active_dofs:
            d = self.dofs[k]
            ms = self._moved_sets[k]
            hit = []
            for L in active:
                if d.linker == L.id and d.name != 'chi':
                    hit.append(L)
                elif d.name != 'chi' and ((self.idx[(L.a, "C4'")] in ms) != (self.idx[(L.b, "C4'")] in ms)):
                    hit.append(L)
            if hit:
                self.dof_loops[k] = hit
        res = {}
        for k in self.active_dofs:
            res.setdefault(self.dofs[k].res, []).append(k)
        self.active_dof_by_res = res

    # ------------------------------------------------------------------ non-bonded exclusions
    def _make_exclusions(self):
        """Lookup of atom-name pairs within 3 covalent bonds (same residue or across one O3'-P link)."""
        from .checks import BACKBONE_BONDS, BASE_BONDS
        names = sorted({nm for b in 'ACGU' for nm in atom_order(b)})
        code = {nm: k for k, nm in enumerate(names)}
        K = len(names)
        self._name_code = np.array([code[nm] for nm in self.atom_name])
        self._K = K
        same = np.zeros((4, K, K), bool)      # per base type of the residue
        link = np.zeros((4, 4, K, K), bool)   # (base i, base i+1): atom of i, atom of i+1
        bidx = {b: k for k, b in enumerate('ACGU')}

        def graph(bases):
            adj = {}
            for r, b in enumerate(bases):
                for x, y in BACKBONE_BONDS + BASE_BONDS[b] + [("C1'", GLY[b])]:
                    adj.setdefault((r, x), set()).add((r, y))
                    adj.setdefault((r, y), set()).add((r, x))
            if len(bases) == 2:
                adj.setdefault((0, "O3'"), set()).add((1, 'P'))
                adj.setdefault((1, 'P'), set()).add((0, "O3'"))
            return adj

        def within3(adj, a):
            seen = {a: 0}
            frontier = [a]
            for depth in range(1, 4):
                nxt = []
                for u in frontier:
                    for v in adj.get(u, ()):
                        if v not in seen:
                            seen[v] = depth
                            nxt.append(v)
                frontier = nxt
            return seen
        for b1 in 'ACGU':
            for b2 in 'ACGU':
                adj = graph([b1, b2])
                for (r, nm) in list(adj):
                    for (r2, nm2) in within3(adj, (r, nm)):
                        if r == 0 and r2 == 0:
                            same[bidx[b1], code[nm], code[nm2]] = True
                        elif r == 0 and r2 == 1:
                            link[bidx[b1], bidx[b2], code[nm], code[nm2]] = True
        self._excl_same = same
        self._excl_link = link
        self._res_bidx = np.array([bidx[b] for b in self.seq])

    def nonbonded_mask(self, i, j):
        """For atom index arrays i, j: True where the pair is NOT within 3 bonds."""
        ri, rj = self.atom_res[i], self.atom_res[j]
        swap = ri > rj
        a = np.where(swap, j, i)
        b = np.where(swap, i, j)
        ra, rb_ = self.atom_res[a], self.atom_res[b]
        ca, cb = self._name_code[a], self._name_code[b]
        excl = np.zeros(len(i), bool)
        s0 = ra == rb_
        if s0.any():
            excl[s0] = self._excl_same[self._res_bidx[ra[s0]], ca[s0], cb[s0]]
        s1 = rb_ - ra == 1
        if s1.any():
            excl[s1] = self._excl_link[self._res_bidx[ra[s1]], self._res_bidx[rb_[s1]], ca[s1], cb[s1]]
        return ~excl

    def atom_overlap_energy(self, cutoff=2.8, polar_cutoff=2.5, skin=2.0):
        """Sum of squared non-bonded heavy-atom overlaps, using a Verlet neighbour list.

        Exact (identical to atom_clashes) as long as no atom moved more than skin/2 since the
        list was built; the list is rebuilt automatically otherwise or when atoms appear/disappear."""
        key = (int(self.present.sum()), cutoff, polar_cutoff)
        vl = getattr(self, '_verlet', None)
        if vl is None or vl['key'] != key or \
                np.max(np.sum((self.xyz[vl['ids']] - vl['ref']) ** 2, axis=1)) > (skin / 2) ** 2:
            from scipy.spatial import cKDTree
            ids = np.where(self.present)[0]
            pairs = cKDTree(self.xyz[ids]).query_pairs(cutoff + skin, output_type='ndarray')
            i, j = ids[pairs[:, 0]], ids[pairs[:, 1]]
            own_i, own_j = self.atom_owner[i], self.atom_owner[j]
            keep = ~((own_i == own_j) & (own_i >= 0) & (own_i < 1000))
            i, j = i[keep], j[keep]
            keep = self.nonbonded_mask(i, j)
            i, j = i[keep], j[keep]
            polar = np.isin(self.element[i], ('N', 'O')) & np.isin(self.element[j], ('N', 'O'))
            vl = {'key': key, 'ids': ids, 'ref': self.xyz[ids].copy(), 'i': i, 'j': j,
                  'lim': np.where(polar, polar_cutoff, cutoff)}
            self._verlet = vl
        if not len(vl['i']):
            return 0.0
        d = np.linalg.norm(self.xyz[vl['i']] - self.xyz[vl['j']], axis=1)
        ov = vl['lim'] - d
        ov = ov[ov > 0]
        return float(np.sum(ov * ov))

    def atom_clashes(self, cutoff=2.8, polar_cutoff=2.5, subset=None):
        """Non-bonded heavy-atom overlaps among present atoms (excluding pairs inside one stem body).

        Returns (i, j, overlap) arrays. Polar N/O pairs use polar_cutoff (hydrogen bonds)."""
        from scipy.spatial import cKDTree
        ids = np.where(self.present)[0]
        if len(ids) < 2:
            return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
        pts = self.xyz[ids]
        tree = cKDTree(pts)
        pairs = tree.query_pairs(cutoff, output_type='ndarray')
        if not len(pairs):
            return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
        i, j = ids[pairs[:, 0]], ids[pairs[:, 1]]
        if subset is not None:
            keep = subset[i] | subset[j]
            i, j = i[keep], j[keep]
        own_i, own_j = self.atom_owner[i], self.atom_owner[j]
        keep = ~((own_i == own_j) & (own_i >= 0) & (own_i < 1000))
        i, j = i[keep], j[keep]
        keep = self.nonbonded_mask(i, j)
        i, j = i[keep], j[keep]
        d = np.linalg.norm(self.xyz[i] - self.xyz[j], axis=1)
        polar = np.isin(self.element[i], ('N', 'O')) & np.isin(self.element[j], ('N', 'O'))
        lim = np.where(polar, polar_cutoff, cutoff)
        ov = lim - d
        m = ov > 0
        return i[m], j[m], ov[m]

    # ------------------------------------------------------------------ setup
    def _index_atoms(self):
        self.idx = {}
        names, res = [], []
        for r, b in enumerate(self.seq):
            for nm in atom_order(b):
                self.idx[(r, nm)] = len(names)
                names.append(nm)
                res.append(r)
        self.atom_name = names
        self.atom_res = np.array(res)
        self.n_atoms = len(names)
        self.element = np.array([('O' if nm.startswith('OP') else nm[0]) for nm in names])

    def _make_bodies(self):
        self.bodies = []
        self.res_body = np.full(self.n, -1)
        self.stem_reports = []
        for k, run in enumerate(stems_from_pairs(self.pairs)):
            res = sorted([p[0] for p in run] + [p[1] for p in run])
            body = Body(len(self.bodies), 'stem', res, run)
            local, rep = build_stem(self.rb, self.seq, run)
            body._local_atoms = local
            self.stem_reports.append(rep)
            for r in res:
                self.res_body[r] = body.id
            self.bodies.append(body)
        if not self.bodies:
            unit = dict(self.rb.unit(self.seq[0], 'C3endo', self.geom.aform['chi']))
            body = Body(0, 'residue', [0])
            body._local_atoms = {0: unit}
            self.res_body[0] = 0
            self.bodies.append(body)

    def _segments(self):
        segs = []
        r = 0
        while r < self.n:
            b = self.res_body[r]
            if b < 0:
                r += 1
                continue
            s = r
            while r + 1 < self.n and self.res_body[r + 1] == b and r + 1 in self.bodies[b].residues \
                    and self._same_strand(b, r, r + 1):
                r += 1
            segs.append((b, s, r))
            r += 1
        return segs

    def _same_strand(self, b, r1, r2):
        body = self.bodies[b]
        if body.kind != 'stem':
            return False
        s5 = {p[0] for p in body.pairs}
        return (r1 in s5) == (r2 in s5)

    def _make_linkers(self):
        segs = self._segments()
        self.segments = segs
        self.linkers = []
        first, last = segs[0], segs[-1]
        if first[1] > 0:
            self.linkers.append(Linker(len(self.linkers), None, first[1], tuple(range(0, first[1])), 'tail5',
                                       None, first[0], 'backward', first[0]))
        for (b1, s1, e1), (b2, s2, e2) in zip(segs, segs[1:]):
            loop = tuple(range(e1 + 1, s2))
            kind = 'hairpin' if b1 == b2 else 'inter'
            self.linkers.append(Linker(len(self.linkers), e1, s2, loop, kind, b1, b2))
        if last[2] < self.n - 1:
            self.linkers.append(Linker(len(self.linkers), last[2], None, tuple(range(last[2] + 1, self.n)), 'tail3',
                                       last[0], None, 'forward', last[0]))

    def _make_tree(self):
        nb = len(self.bodies)
        parent = list(range(nb))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        inter = [L for L in self.linkers if L.kind == 'inter']
        # Kruskal: short linkers first; ties -> prefer linkers nearer the 5' end
        tree = []
        for L in sorted(inter, key=lambda L: (L.n, L.a)):
            ra, rb_ = find(L.body_a), find(L.body_b)
            if ra != rb_:
                parent[ra] = rb_
                tree.append(L)
                L.kind = 'tree'
            else:
                L.kind = 'closure'
        # root: tree centre (minimise eccentricity), tie -> larger body
        adj = {k: [] for k in range(nb)}
        for L in tree:
            adj[L.body_a].append(L.body_b)
            adj[L.body_b].append(L.body_a)

        def ecc(s):
            dist = {s: 0}
            q = [s]
            while q:
                u = q.pop(0)
                for v in adj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        q.append(v)
            return max(dist.values())
        self.root = min(range(nb), key=lambda k: (ecc(k), -len(self.bodies[k].residues)))
        # orient
        self.tree_order = []
        placed = {self.root}
        frontier = [self.root]
        while frontier:
            u = frontier.pop(0)
            for L in tree:
                if L.parent is not None:
                    continue
                if L.body_a == u and L.body_b not in placed:
                    L.parent, L.child, L.direction = u, L.body_b, 'forward'
                elif L.body_b == u and L.body_a not in placed:
                    L.parent, L.child, L.direction = u, L.body_a, 'backward'
                else:
                    continue
                placed.add(L.child)
                frontier.append(L.child)
                self.tree_order.append(L)
        for L in self.linkers:
            if L.kind in ('closure', 'hairpin'):
                L.direction = 'forward'
                # where CCD closes the loop: at O3'-P with b's beta/gamma kept canonical ('p'),
                # or with b's beta/gamma free ('c5'; needed for triloops and 0-nt junctions)
                L.close_at = 'c5' if (L.n == 0 or (L.kind == 'hairpin' and L.n <= 3)) else 'p'
        # tails are built once their parent body is placed
        self.tail_linkers = [L for L in self.linkers if L.kind in ('tail5', 'tail3')]

    def _assign_atoms(self):
        owner = np.full(self.n_atoms, -1)
        strand_starts = {s for (_, s, _) in self.segments}
        for body in self.bodies:
            ids, loc = [], []
            for r in body.residues:
                for nm, xyz in body._local_atoms[r].items():
                    if nm in ('P', 'OP1', 'OP2'):
                        if r in strand_starts:
                            continue
                    if nm == "O5'" and r in strand_starts and r != 0:
                        continue
                    ids.append(self.idx[(r, nm)])
                    loc.append(xyz)
                if r == 0 and "O5'" not in body._local_atoms[r]:
                    o5 = self.rb.o5_from_sugar(body._local_atoms[r], self.geom.aform['gamma'])
                    ids.append(self.idx[(r, "O5'")])
                    loc.append(o5)
            body.atoms = np.array(ids)
            body.local = np.array(loc)
            owner[body.atoms] = body.id
        for L in self.linkers:
            ids = []
            for r in L.loop:
                for nm in atom_order(self.seq[r]):
                    if r == 0 and nm in ('P', 'OP1', 'OP2'):
                        continue
                    ids.append(self.idx[(r, nm)])
            if L.b is not None:
                ids += [self.idx[(L.b, nm)] for nm in ('P', 'OP1', 'OP2', "O5'")]
            L.atoms = np.array(ids, dtype=int)
            owner[L.atoms] = 1000 + L.id
        # residue 0 phosphate never exists
        self.never = np.array([self.idx[(0, nm)] for nm in ('P', 'OP1', 'OP2')])
        self.atom_owner = owner

    # ------------------------------------------------------------------ build plans
    def _plan(self, L, dangling=False):
        """Ordered steps: ('dof', res, name, u, v) / ('atoms', [ids]) / ('chi', res, u, v, [ids]).

        dangling=True (explicit closure/hairpin loops): the chain stops after alpha of b; beta and
        gamma of b stay fixed and the free end is pulled onto b by the chain-break restraint."""
        I = self.idx
        steps = []
        if L.direction == 'forward':
            prev = L.a
            seq = list(L.loop) + ([L.b] if L.b is not None else [])
            for m in seq:
                steps.append(('dof', prev, 'epsilon', I[(prev, "C3'")], I[(prev, "O3'")]))
                steps.append(('atoms', [I[(m, 'P')]]))
                steps.append(('dof', prev, 'zeta', I[(prev, "O3'")], I[(m, 'P')]))
                steps.append(('atoms', [I[(m, "O5'")], I[(m, 'OP1')], I[(m, 'OP2')]]))
                steps.append(('dof', m, 'alpha', I[(m, 'P')], I[(m, "O5'")]))
                if m == L.b and dangling:
                    break
                if m == L.b:
                    steps.append(('dof', m, 'beta', I[(m, "O5'")], I[(m, "C5'")]))
                    steps.append(('dof', m, 'gamma', I[(m, "C5'")], I[(m, "C4'")]))
                    break
                steps.append(('atoms', [I[(m, "C5'")]]))
                steps.append(('dof', m, 'beta', I[(m, "O5'")], I[(m, "C5'")]))
                steps.append(('atoms', [I[(m, "C4'")]]))
                steps.append(('dof', m, 'gamma', I[(m, "C5'")], I[(m, "C4'")]))
                base = [I[(m, nm)] for nm in BASE_ATOMS[self.seq[m]]]
                steps.append(('atoms', [I[(m, nm)] for nm in ("C3'", "O3'", "O4'", "C2'", "O2'", "C1'")] + base))
                steps.append(('chi', m, I[(m, "C1'")], I[(m, GLY[self.seq[m]])], base))
                prev = m
        else:
            nxt = L.b
            seq = list(reversed(L.loop)) + ([L.a] if L.a is not None else [])
            for m in seq:
                steps.append(('dof', nxt, 'gamma', I[(nxt, "C4'")], I[(nxt, "C5'")]))
                steps.append(('atoms', [I[(nxt, "O5'")]]))
                steps.append(('dof', nxt, 'beta', I[(nxt, "C5'")], I[(nxt, "O5'")]))
                steps.append(('atoms', [I[(nxt, 'P')]]))
                steps.append(('dof', nxt, 'alpha', I[(nxt, "O5'")], I[(nxt, 'P')]))
                if m == L.a:
                    steps.append(('atoms', [I[(nxt, 'OP1')], I[(nxt, 'OP2')]]))
                    steps.append(('dof', m, 'zeta', I[(nxt, 'P')], I[(m, "O3'")]))
                    steps.append(('dof', m, 'epsilon', I[(m, "O3'")], I[(m, "C3'")]))
                    break
                steps.append(('atoms', [I[(m, "O3'")], I[(nxt, 'OP1')], I[(nxt, 'OP2')]]))
                steps.append(('dof', m, 'zeta', I[(nxt, 'P')], I[(m, "O3'")]))
                steps.append(('atoms', [I[(m, "C3'")]]))
                steps.append(('dof', m, 'epsilon', I[(m, "O3'")], I[(m, "C3'")]))
                base = [I[(m, nm)] for nm in BASE_ATOMS[self.seq[m]]]
                steps.append(('atoms', [I[(m, nm)] for nm in ("C4'", "C5'", "O4'", "C2'", "O2'", "C1'")] + base))
                steps.append(('chi', m, I[(m, "C1'")], I[(m, GLY[self.seq[m]])], base))
                nxt = m
            if L.a is None:          # 5' tail: O5' of residue 0 hangs from its sugar
                steps.append(('dof', 0, 'gamma', I[(0, "C4'")], I[(0, "C5'")]))
                steps.append(('atoms', [I[(0, "O5'")]]))
        return steps

    def _subtree_atoms(self, body_id):
        out = [self.bodies[body_id].atoms]
        for L in self.tree_order:
            if L.parent == body_id:
                out.append(L.atoms)
                out.append(self._subtree_atoms(L.child))
        for L in self.tail_linkers:
            if L.parent == body_id:
                out.append(L.atoms)
        # explicit closure/hairpin loops are built forward from their upstream anchor a:
        # they move rigidly with the body that contains a
        for L in self.linkers:
            if L.kind in ('closure', 'hairpin') and L.body_a == body_id:
                out.append(L.atoms)
        return np.concatenate(out) if out else np.zeros(0, int)

    def _make_dofs(self):
        self.dofs = []
        self.plans = {}
        self.loop_linkers = [L for L in self.linkers if L.kind in ('closure', 'hairpin')]
        movable = list(self.tree_order) + self.tail_linkers + (self.loop_linkers if self.explicit_loops else [])
        for L in movable:
            steps = self._plan(L, dangling=L.kind in ('closure', 'hairpin'))
            self.plans[L.id] = steps
            child_atoms = self._subtree_atoms(L.child) if L.kind == 'tree' else np.zeros(0, int)
            for k, st in enumerate(steps):
                if st[0] == 'dof':
                    after = [np.asarray(s[1], int) for s in steps[k + 1:] if s[0] == 'atoms']
                    moved = np.unique(np.concatenate(after + [child_atoms])) if (after or len(child_atoms)) else np.zeros(0, int)
                    self.dofs.append(Dof(st[1], st[2], st[3], st[4], moved, L.id))
                elif st[0] == 'chi':
                    self.dofs.append(Dof(st[1], 'chi', st[2], st[3], np.asarray(st[4], int), L.id))
        if not self.explicit_loops:
            for L in self.loop_linkers:
                self.plans[L.id] = self._plan(L)
        self.dof_by_res = {}
        for k, d in enumerate(self.dofs):
            self.dof_by_res.setdefault(d.res, []).append(k)

    # ------------------------------------------------------------------ torsion initialisation
    def sample_torsions(self, rng, p_aform_junction=0.5):
        """Draw torsions for all degrees of freedom (loop residues: whole library rows)."""
        lib = self.geom.loop_library
        aform = self.geom.aform
        free_res = set()
        for L in self.linkers:
            free_res.update(L.loop)
        for r in free_res:
            row = lib[rng.integers(len(lib))]
            for k, nm in enumerate(LIB_COLS):
                if nm in TOR:
                    self.tor[r, TOR[nm]] = row[k]
            self.puck[r] = 'C2endo' if row[3] > 115 else 'C3endo'
        for L in self.linkers:
            junction = []
            if L.a is not None:
                junction += [(L.a, 'epsilon'), (L.a, 'zeta')]
            if L.b is not None:
                junction += [(L.b, 'alpha')] if L.kind in ('closure', 'hairpin') else \
                    [(L.b, 'alpha'), (L.b, 'beta'), (L.b, 'gamma')]
            if L.kind in ('closure', 'hairpin'):
                self.tor[L.b, TOR['beta']] = aform['beta']
                self.tor[L.b, TOR['gamma']] = aform['gamma']
            use_aform = L.n == 0 and rng.random() < p_aform_junction
            row = lib[rng.integers(len(lib))]
            for r, nm in junction:
                self.tor[r, TOR[nm]] = aform[nm] if use_aform else row[LIB_COLS.index(nm)]
        if 0 not in free_res:
            self.tor[0, TOR['gamma']] = aform['gamma']

    # ------------------------------------------------------------------ building
    def _res_atoms(self, r, names):
        return {nm: self.xyz[self.idx[(r, nm)]] for nm in names}

    def _write(self, r, atoms):
        for nm, v in atoms.items():
            k = self.idx.get((r, nm))
            if k is not None:
                self.xyz[k] = v
                self.present[k] = True

    def _place_body(self, body, src_names_res, built):
        """Place body so that its atoms (res, names) coincide with `built` coordinates."""
        loc_index = {int(a): k for k, a in enumerate(body.atoms)}
        src = np.array([body.local[loc_index[self.idx[key]]] for key in src_names_res])
        R, t = kabsch(src, np.asarray(built))
        self.xyz[body.atoms] = body.local @ R.T + t
        self.present[body.atoms] = True

    def _build_linker(self, L):
        rb = self.rb
        T = self.tor
        if L.direction == 'forward':
            prev = L.a
            prev_atoms = self._res_atoms(prev, ("C4'", "C3'", "O3'"))
            seq = list(L.loop) + ([L.b] if L.b is not None else [])
            for m in seq:
                P, O5, C5, C4, C3 = rb.forward_link(prev_atoms, T[prev, 3], T[prev, 4], T[m, 0], T[m, 1], T[m, 2], self.puck[m])
                op = rb.phosphate_oxygens(prev_atoms["O3'"], P, O5)
                self._write(m, {'P': P, "O5'": O5, 'OP1': op[0], 'OP2': op[1]})
                if m == L.b:
                    return np.array([C5, C4, C3])
                unit = rb.place_unit_forward(self.seq[m], self.puck[m], T[m, 5], C5, C4, C3)
                self._write(m, unit)
                prev, prev_atoms = m, unit
            return None
        nxt = L.b
        nxt_atoms = self._res_atoms(nxt, ("C3'", "C4'", "C5'"))
        seq = list(reversed(L.loop)) + ([L.a] if L.a is not None else [])
        for m in seq:
            O5, P, O3, C3, C4 = rb.backward_link(nxt_atoms, T[nxt, 2], T[nxt, 1], T[nxt, 0], T[m, 4], T[m, 3], self.puck[m])
            op = rb.phosphate_oxygens(O3, P, O5)
            self._write(nxt, {'P': P, "O5'": O5, 'OP1': op[0], 'OP2': op[1]})
            if m == L.a:
                return np.array([O3, C3, C4])
            unit = rb.place_unit_backward(self.seq[m], self.puck[m], T[m, 5], O3, C3, C4)
            self._write(m, unit)
            nxt, nxt_atoms = m, unit
        if L.a is None:
            self._write(0, {"O5'": rb.o5_from_sugar(nxt_atoms, T[0, 2])})
        return None

    def build(self):
        """Place every body and tree/tail linker from the current torsions (closure loops absent)."""
        self.xyz[:] = 0.0
        self.present[:] = False
        root = self.bodies[self.root]
        self.xyz[root.atoms] = root.local
        self.present[root.atoms] = True
        for L in self.tree_order:
            built = self._build_linker(L)
            child = self.bodies[L.child]
            if L.direction == 'forward':
                self._place_body(child, [(L.b, "C5'"), (L.b, "C4'"), (L.b, "C3'")], built)
            else:
                self._place_body(child, [(L.a, "O3'"), (L.a, "C3'"), (L.a, "C4'")], built)
        for L in self.tail_linkers:
            self._build_linker(L)
        for L in self.loop_linkers:
            if L.id in self.active_loops:
                self._build_linker(L)
        self.present[self.never] = False

    def chainbreak(self, L):
        """Squared distance (A^2, summed over P, O5', C5') between the dangling end of loop
        linker L and the positions required by residue b's fixed sugar (beta/gamma of b fixed)."""
        I = self.idx
        rb = self.rb
        resb = {n: self.xyz[I[(L.b, n)]] for n in ("C3'", "C4'", "C5'")}
        O5t, Pt = rb.phosphate_from_sugar(resb, self.tor[L.b, 2], self.tor[L.b, 1])
        prev_o3 = self.xyz[I[((L.b - 1), "O3'")]]
        P = self.xyz[I[(L.b, 'P')]]
        O5 = self.xyz[I[(L.b, "O5'")]]
        from .geomutil import place as _place
        C5g = _place(prev_o3, P, O5, rb.b_O5C5, rb.a_PO5C5, self.tor[L.b, 0])
        return float(np.sum((P - Pt) ** 2) + np.sum((O5 - O5t) ** 2) + np.sum((C5g - resb["C5'"]) ** 2))

    # ------------------------------------------------------------------ moves
    def rotate(self, k, delta_deg):
        d = self.dofs[k]
        if len(d.moved) == 0 or delta_deg == 0.0:
            self.tor[d.res, TOR[d.name]] += delta_deg
            return
        u = self.xyz[d.u]
        v = self.xyz[d.v]
        R = rotation_about_axis(v - u, math.radians(delta_deg))
        self.xyz[d.moved] = (self.xyz[d.moved] - v) @ R.T + v
        val = self.tor[d.res, TOR[d.name]] + delta_deg
        self.tor[d.res, TOR[d.name]] = (val + 180.0) % 360.0 - 180.0

    def measured_torsion(self, k):
        """Recompute a DOF torsion from coordinates (for consistency tests)."""
        d = self.dofs[k]
        r, nm = d.res, d.name
        I = self.idx
        x = self.xyz
        if nm == 'alpha':
            return dihedral(x[I[(r - 1, "O3'")]], x[I[(r, 'P')]], x[I[(r, "O5'")]], x[I[(r, "C5'")]])
        if nm == 'beta':
            return dihedral(x[I[(r, 'P')]], x[I[(r, "O5'")]], x[I[(r, "C5'")]], x[I[(r, "C4'")]])
        if nm == 'gamma':
            return dihedral(x[I[(r, "O5'")]], x[I[(r, "C5'")]], x[I[(r, "C4'")]], x[I[(r, "C3'")]])
        if nm == 'epsilon':
            return dihedral(x[I[(r, "C4'")]], x[I[(r, "C3'")]], x[I[(r, "O3'")]], x[I[(r + 1, 'P')]])
        if nm == 'zeta':
            return dihedral(x[I[(r, "C3'")]], x[I[(r, "O3'")]], x[I[(r + 1, 'P')]], x[I[(r + 1, "O5'")]])
        b = self.seq[r]
        ref = 'C4' if b in 'AG' else 'C2'
        return dihedral(x[I[(r, "O4'")]], x[I[(r, "C1'")]], x[I[(r, GLY[b])]], x[I[(r, ref)]])

    # ------------------------------------------------------------------ output
    def residues_atoms(self):
        """List (per residue) of {atom: xyz} for present atoms, in canonical atom order."""
        out = []
        for r, b in enumerate(self.seq):
            atoms = {}
            for nm in atom_order(b):
                k = self.idx[(r, nm)]
                if self.present[k]:
                    atoms[nm] = self.xyz[k].copy()
            out.append(atoms)
        return out

    def summary(self):
        kinds = {}
        for L in self.linkers:
            kinds[L.kind] = kinds.get(L.kind, 0) + 1
        return {'n': self.n, 'n_pairs': len(self.pairs), 'n_stems': sum(b.kind == 'stem' for b in self.bodies),
                'linkers': kinds, 'n_dofs': len(self.dofs), 'root_body': self.root,
                'closure_linkers': [{'a': L.a + 1, 'b': L.b + 1, 'n': L.n, 'kind': L.kind}
                                    for L in self.linkers if L.kind in ('closure', 'hairpin')],
                'max_stem_pair_fit_rmsd': max((r['pair_fit_rmsd'] for r in self.stem_reports), default=0.0)}
