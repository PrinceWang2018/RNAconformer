"""Coarse scoring and simulated-annealing search over fold-tree torsions (v2 engine).

Score terms (dimensionless, kT-like model units; NOT free energies):
  clash        soft overlap between coarse sites (P, C4', C1', base centroid) of residues
               in different rigid groups, |i-j| > 1
  feasibility  closure linkers whose anchors are farther apart than the calibrated
               99.5th-percentile reach of an n-nt loop
  compaction   one-sided harmonic restraint Rg > Rg0 = 5.5 N^(1/3) A (empirical scaling of
               folded RNAs, Hyeon, Dima & Thirumalai 2006) - a mean-field stand-in for
               unknown tertiary contacts
  stacking     reward for coaxially stacked helix ends (distance and facing geometry)
  hairpin      excluded volume of not-yet-built hairpin loops (pseudo-site)
  electro      screened Debye-Hueckel phosphate repulsion (lambda = 10 A, reduced charge)
  torsion      -ln p(torsion) of calibrated loop residues for every free torsion
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math

import numpy as np
from scipy.spatial import cKDTree

from .build import BASE_ATOMS
from .fold import FoldModel, TOR, LIB_COLS
from .loops import TorsionPrior

SITE_P, SITE_C4, SITE_C1, SITE_B = range(4)
THR = np.array([[4.5, 3.6, 3.6, 3.6],
                [3.6, 3.6, 3.6, 3.4],
                [3.6, 3.6, 3.6, 3.4],
                [3.6, 3.4, 3.4, 3.0]])


@dataclass(frozen=True)
class ScoreWeights:
    clash: float = 3.0
    atom_clash: float = 3.0
    feasibility: float = 3.0
    chainbreak: float = 3.0
    compaction: float = 0.3
    rg_coefficient: float = 5.5
    compaction_branched_only: bool = True   # apply Rg prior only to branched 2D (multiloop or >1 root stem)
    stacking: float = 4.0
    hairpin: float = 1.0
    electro: float = 1.5
    debye_length: float = 10.0
    torsion: float = 1.0
    restraint: float = 0.5                  # crossing (pseudoknot) / tertiary base-pair restraints
    pull: float = 2.0                       # interactive pull restraints (bias; removed before relaxation)


def is_branched(pairs):
    """True if the nested 2D has a multiloop (a pair enclosing >= 2 helices) or >= 2 root helices."""
    from .fold import stems_from_pairs
    stems = stems_from_pairs(sorted(pairs))
    if not stems:
        return False
    outer = [(s[0][0], s[0][1]) for s in stems]
    inner = [(s[-1][0], s[-1][1]) for s in stems]

    def children(lo, hi):
        kids = [o for o in outer if lo < o[0] and o[1] < hi]
        return [k for k in kids if not any(o != k and o[0] < k[0] and k[1] < o[1] and lo < o[0] and o[1] < hi
                                               for o in outer)]
    if len(children(-1, 10 ** 9)) >= 2:
        return True
    return any(len(children(i, j)) >= 2 for i, j in inner)


CENTRAL_HBOND = {('G', 'C'): ('N1', 'N3'), ('C', 'G'): ('N3', 'N1'), ('A', 'U'): ('N1', 'N3'),
                 ('U', 'A'): ('N3', 'N1'), ('G', 'U'): ('N1', 'O2'), ('U', 'G'): ('O2', 'N1')}


class Scorer:
    def __init__(self, model: FoldModel, weights: ScoreWeights = ScoreWeights(), prior: TorsionPrior | None = None,
                 restraints=(), pull_restraints=()):
        self.m = model
        self.w = weights
        self.prior = prior or TorsionPrior(model.geom)
        I = model.idx
        n = model.n
        self.p_idx = np.array([I[(r, 'P')] for r in range(n)])
        self.c4_idx = np.array([I[(r, "C4'")] for r in range(n)])
        self.c1_idx = np.array([I[(r, "C1'")] for r in range(n)])
        flat, offs = [], []
        for r, b in enumerate(model.seq):
            offs.append(len(flat))
            flat += [I[(r, nm)] for nm in BASE_ATOMS[b] if nm not in ('N6', 'O6', 'N2', 'N4', 'O4', 'O2')]
        self.base_flat = np.array(flat)
        self.base_offs = np.array(offs)
        self.base_cnt = np.diff(np.append(self.base_offs, len(flat)))
        # residue group: body for stem residues, linker id otherwise
        group = np.array(model.res_body, dtype=int)
        for L in model.linkers:
            for r in L.loop:
                group[r] = 1000 + L.id
        self.group = group
        self.stem_body = np.array([model.bodies[g].kind == 'stem' if g >= 0 and g < 1000 else False for g in group])
        # closure feasibility
        reach = model.geom.raw['loop_reach']
        self.closures = []
        for L in model.linkers:
            if L.kind == 'closure':
                nn = min(L.n, len(reach['p995']) - 1)
                hi = reach['p995'][nn] + 1.0 + max(0, L.n - nn) * 5.5
                self.closures.append((I[(L.a, "O3'")], I[(L.b, "C5'")], hi, 2.5 if L.n == 0 else 0.0, L.id))
        # helix ends
        self.ends = []      # (body, (i, j) end pair, (i2, j2) neighbour pair or None)
        for body in model.bodies:
            if body.kind != 'stem':
                continue
            p = body.pairs
            self.ends.append((body.id, p[0], p[1] if len(p) > 1 else None, -1))
            self.ends.append((body.id, p[-1], p[-2] if len(p) > 1 else None, +1))
        self.hairpins = []
        self.loop_linkers = [L for L in model.linkers if L.kind in ('closure', 'hairpin')]
        for L in model.linkers:
            if L.kind == 'hairpin':
                end = next(k for k, e in enumerate(self.ends) if e[1] == (L.a, L.b))
                self.hairpins.append((end, 4.5 + 0.35 * max(0, L.n - 4), 4.0 + 0.5 * max(0, L.n - 4), 1000 + L.id, L.id))
        self.rg0 = weights.rg_coefficient * n ** (1.0 / 3.0)
        # base-pair restraints (pairs that cannot be part of the nested fold tree)
        loop_of = {}
        for L in model.linkers:
            if L.kind in ('closure', 'hairpin'):
                for r in L.loop:
                    loop_of[r] = L
        self.restraints = []
        for i, j in restraints:
            i, j = min(i, j), max(i, j)
            hb = CENTRAL_HBOND.get((model.seq[i], model.seq[j]))
            item = {'i': i, 'j': j, 'c1': (I[(i, "C1'")], I[(j, "C1'")]),
                    'hb': None if hb is None else (I[(i, hb[0])], I[(j, hb[1])])}
            for key, r in (('pi', i), ('pj', j)):
                L = loop_of.get(r)
                item[key] = None if L is None else (I[(L.a, "C4'")], I[(L.b, "C4'")], L.id,
                                                      min(6.0 * min(r - L.a, L.b - r), 6.0 + 1.0 * L.n))
            self.restraints.append(item)
        # pull restraints: target distances between atoms taken from an interactively deformed structure.
        # Atoms in loops that are not built yet are represented by the midpoint of the loop's flanking C4'
        # atoms, with the loop's reach added to the tolerance.
        self.pull = []
        self.pull_scale = 1.0
        for p in pull_restraints:
            ends = []
            for r, atom in ((p['i'], p.get('a', "C4'")), (p['j'], p.get('b', "C4'"))):
                L = loop_of.get(r)
                proxy = None if L is None else (I[(L.a, "C4'")], I[(L.b, "C4'")],
                                                min(6.0 * min(r - L.a, L.b - r), 6.0 + 1.0 * L.n))
                ends.append((I[(r, atom)], proxy))
            self.pull.append((ends[0], ends[1], float(p['target']), float(p.get('tol', 1.5))))
        self.branched = is_branched(model.pairs)
        self.compaction_on = (not weights.compaction_branched_only) or self.branched
        self.clash_scale = 1.0
        self.atom_clash_from = 1.0      # all-atom sterics only once the steric ramp reaches this scale
        self.dof_names = [d.name for d in model.dofs]
        self.dof_res = np.array([d.res for d in model.dofs], dtype=int)
        self.dof_col = np.array([TOR[d.name] for d in model.dofs], dtype=int)

    # ------------------------------------------------------------------ sites
    def base_centroids(self):
        x = self.m.xyz[self.base_flat]
        return np.add.reduceat(x, self.base_offs, axis=0) / self.base_cnt[:, None]

    def _end_geometry(self, bc):
        """centres, outward unit vectors and validity of helix ends."""
        m = self.m
        x = m.xyz
        cen = np.zeros((len(self.ends), 3))
        out = np.zeros((len(self.ends), 3))
        for k, (body, (i, j), nb, sign) in enumerate(self.ends):
            c = 0.5 * (bc[i] + bc[j])
            cen[k] = c
            if nb is not None:
                v = c - 0.5 * (bc[nb[0]] + bc[nb[1]])
            else:
                a = x[self.c1_idx[j]] - x[self.c1_idx[i]]
                g = c - 0.5 * (x[self.c1_idx[i]] + x[self.c1_idx[j]])
                v = np.cross(a, g) * sign
            nv = np.linalg.norm(v)
            out[k] = v / nv if nv > 1e-6 else 0.0
        return cen, out

    # ------------------------------------------------------------------ energy
    def __call__(self, detail=False):
        m = self.m
        w = self.w
        x = m.xyz
        pres = m.present
        bc = self.base_centroids()
        res_ok = pres[self.c4_idx]
        rids = np.where(res_ok)[0]
        # coarse sites
        pts = [x[self.c4_idx[rids]], x[self.c1_idx[rids]], bc[rids]]
        typ = [np.full(len(rids), SITE_C4), np.full(len(rids), SITE_C1), np.full(len(rids), SITE_B)]
        rr = [rids, rids, rids]
        p_ok = rids[pres[self.p_idx[rids]]]
        pts.append(x[self.p_idx[p_ok]])
        typ.append(np.full(len(p_ok), SITE_P))
        rr.append(p_ok)
        S = np.vstack(pts)
        T = np.concatenate(typ)
        R = np.concatenate(rr)
        G = self.group[R]
        tree = cKDTree(S)
        pairs = tree.query_pairs(4.5, output_type='ndarray')
        e_clash = 0.0
        if len(pairs):
            i, j = pairs[:, 0], pairs[:, 1]
            keep = (np.abs(R[i] - R[j]) > 1) & ~((G[i] == G[j]) & (G[i] < 1000))
            i, j = i[keep], j[keep]
            if len(i):
                d = np.linalg.norm(S[i] - S[j], axis=1)
                ov = THR[T[i], T[j]] - d
                ov = ov[ov > 0]
                e_clash = float(np.sum(ov ** 2))
        # electrostatics on phosphates
        e_el = 0.0
        if w.electro and len(p_ok) > 1:
            P = x[self.p_idx[p_ok]]
            gp = self.group[p_ok]
            pp = cKDTree(P).query_pairs(2.0 * w.debye_length, output_type='ndarray')
            if len(pp):
                a, b = pp[:, 0], pp[:, 1]
                keep = ~((gp[a] == gp[b]) & (gp[a] < 1000)) & (np.abs(p_ok[a] - p_ok[b]) > 1)
                a, b = a[keep], b[keep]
                d = np.maximum(np.linalg.norm(P[a] - P[b], axis=1), 3.0)
                e_el = float(np.sum(1.8 * np.exp(-d / w.debye_length) / d))
        # compaction (C4' sites of present residues)
        c4 = x[self.c4_idx[rids]]
        rg = float(np.sqrt(np.mean(np.sum((c4 - c4.mean(axis=0)) ** 2, axis=1))))
        e_rg = max(0.0, rg - self.rg0) ** 2 if self.compaction_on else 0.0
        # feasibility of closure linkers
        e_f = 0.0
        for ia, ib, hi, lo, lid in self.closures:
            if lid in m.active_loops:
                continue
            d = float(np.linalg.norm(x[ia] - x[ib]))
            if d > hi:
                e_f += (d - hi) ** 2
            elif d < lo:
                e_f += (lo - d) ** 2
        # helix-end stacking
        e_st = 0.0
        e_hp = 0.0
        n_stacks = 0
        if self.ends:
            cen, out = self._end_geometry(bc)
            if len(self.ends) > 1:
                ends_body = np.array([e[0] for e in self.ends])
                diff = cen[None, :, :] - cen[:, None, :]
                dist = np.linalg.norm(diff, axis=2)
                iu, ju = np.where(np.triu(dist < 7.5, 1))
                for a, b in zip(iu, ju):
                    if ends_body[a] == ends_body[b]:
                        continue
                    d = dist[a, b]
                    v = diff[a, b] / d
                    f1 = float(np.dot(v, out[a]))
                    f2 = float(-np.dot(v, out[b]))
                    if f1 <= 0.3 or f2 <= 0.3:
                        continue
                    s = math.exp(-0.5 * ((d - 3.9) / 1.0) ** 2) * min(1.0, f1 / 0.8) * min(1.0, f2 / 0.8)
                    e_st -= s
                    n_stacks += s > 0.5
            # hairpin loop pseudo-sites
            hps = [h for h in self.hairpins if h[4] not in m.active_loops]
            if hps:
                hp = np.array([cen[k] + out[k] * off for k, off, rad, _, _ in hps])
                rad = np.array([r for _, _, r, _, _ in hps])
                hp_body = [self.ends[k][0] for k, _, _, _, _ in hps]
                hp_loop = [g for _, _, _, g, _ in hps]
                hits = cKDTree(S).query_ball_point(hp, rad.max() + 1.0)
                for h, (pt, rr_, body, own) in enumerate(zip(hp, rad, hp_body, hp_loop)):
                    if not hits[h]:
                        continue
                    idx = np.array(hits[h])
                    idx = idx[(G[idx] != body) & (G[idx] != own)]
                    if len(idx):
                        d = np.linalg.norm(S[idx] - pt, axis=1)
                        ov = rr_ - d
                        e_hp += float(np.sum(ov[ov > 0] ** 2))
        # torsion prior on free torsions
        e_t = 0.0
        act = m.active_dofs
        if act:
            key = len(act)
            if getattr(self, '_tor_key', None) != key:
                self._tor_key = key
                self._tor_sel = {name: np.array([k for k in act if self.dof_names[k] == name], int)
                                 for name in ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi')}
            for name, sel in self._tor_sel.items():
                if len(sel):
                    e_t += float(np.sum(self.prior(name, m.tor[self.dof_res[sel], self.dof_col[sel]])))
        e_ac = 0.0
        if w.atom_clash and self.clash_scale >= self.atom_clash_from:
            e_ac = m.atom_overlap_energy()
        e_rs = 0.0
        for rs in self.restraints:
            ia, ib = rs['c1']
            if pres[ia] and pres[ib]:
                d = float(np.linalg.norm(x[ia] - x[ib]))
                e_rs += max(0.0, abs(d - 10.5) - 0.7) ** 2
                if rs['hb'] is not None:
                    d = float(np.linalg.norm(x[rs['hb'][0]] - x[rs['hb'][1]]))
                    e_rs += max(0.0, abs(d - 2.95) - 0.3) ** 2
            else:
                # proxies for residues in loops not yet built: closing-pair midpoint + reach allowance
                pts, allow = [], 0.0
                for key, k in (('pi', ia), ('pj', ib)):
                    if pres[k]:
                        pts.append(x[k])
                    else:
                        pa, pb, _, al = rs[key]
                        pts.append(0.5 * (x[pa] + x[pb]))
                        allow += al
                d = float(np.linalg.norm(pts[0] - pts[1]))
                e_rs += 0.25 * max(0.0, d - 10.5 - allow) ** 2
        e_pull = 0.0
        if self.pull and self.pull_scale:
            e_pull = self.pull_energy()
        e_cb = 0.0
        for L in self.loop_linkers:
            if L.id in m.active_loops:
                e_cb += m.chainbreak(L)
        cs = self.clash_scale
        total = (cs * w.clash * e_clash + cs * w.atom_clash * e_ac + cs * w.chainbreak * e_cb
                 + w.restraint * e_rs + w.feasibility * e_f + w.compaction * e_rg + w.stacking * e_st
                 + w.hairpin * e_hp + w.electro * e_el + w.torsion * e_t)
        pull_bias = self.pull_scale * w.pull * e_pull
        total += pull_bias
        if not detail:
            return total
        return {'total': total, 'clash': cs * w.clash * e_clash, 'atom_clash': cs * w.atom_clash * e_ac,
                'clash_scale': cs, 'restraint': w.restraint * e_rs, 'chainbreak': cs * w.chainbreak * e_cb, 'chainbreak_sq_A2': e_cb, 'feasibility': w.feasibility * e_f,
                'compaction': w.compaction * e_rg, 'stacking': w.stacking * e_st, 'hairpin': w.hairpin * e_hp,
                'electro': w.electro * e_el, 'torsion': w.torsion * e_t, 'rg': rg, 'rg0': self.rg0,
                'compaction_applied': self.compaction_on,
                'n_stacked_helix_ends': int(n_stacks), 'pull_bias': pull_bias}

    def pull_energy(self):
        """Flat-bottom harmonic energy of the pull restraints (A^2)."""
        x, pres = self.m.xyz, self.m.present
        e = 0.0
        for (ia, pa), (ib, pb), target, tol in self.pull:
            allow = 0.0
            pts = []
            for idx, proxy in ((ia, pa), (ib, pb)):
                if pres[idx]:
                    pts.append(x[idx])
                elif proxy is not None and pres[proxy[0]] and pres[proxy[1]]:
                    pts.append(0.5 * (x[proxy[0]] + x[proxy[1]]))
                    allow += proxy[2]
                else:
                    pts = None
                    break
            if pts is None:
                continue
            d = float(np.linalg.norm(pts[0] - pts[1]))
            e += max(0.0, abs(d - target) - tol - allow) ** 2
        return e

    def pull_deviation(self):
        """RMS deviation (A) of built pull-restraint distances from their targets."""
        x, pres = self.m.xyz, self.m.present
        dev = [float(np.linalg.norm(x[ia] - x[ib])) - t for (ia, _), (ib, _), t, _ in self.pull if pres[ia] and pres[ib]]
        return float(np.sqrt(np.mean(np.square(dev)))) if dev else None


@dataclass(frozen=True)
class AnnealConfig:
    steps: int = 0                 # 0 -> automatic (3000 + 40 * n_dofs, capped)
    t_start: float = 6.0
    t_end: float = 0.15
    sigma_start: float = 35.0
    sigma_end: float = 6.0
    p_residue_jump: float = 0.2
    p_junction_aform: float = 0.05
    max_steps: int = 10000
    reclose_sweeps: int = 10         # CCD sweeps to re-close explicit loops after every move (0 = off)
    clash_scale_start: float = 0.3   # soft sterics early (chains may pass through each other)
    clash_scale_end: float = 25.0    # hard sterics at the end
    quench_fraction: float = 0.2     # final fraction at t_end and full clash weight


def reclose(model: FoldModel, ks, sweeps=10, tol=0.05):
    """Re-close the explicit loops perturbed by DOFs ks with a few CCD sweeps (in place)."""
    if not sweeps:
        return
    from .loops import ccd_close
    seen = set()
    for k in ks:
        for L in model.dof_loops.get(k, ()):
            if L.id in seen:
                continue
            seen.add(L.id)
            close_at = getattr(L, 'close_at', 'p')
            frozen = {(L.b, 'beta'), (L.b, 'gamma')} if close_at == 'p' else set()
            ok, rmsd, tor = ccd_close(model, L, model.tor, model.puck, max_sweeps=sweeps, tol=tol,
                                      close_at=close_at, frozen=frozen)
            model.tor[:] = tor
            model._build_linker(L)
    model.present[model.never] = False


def anneal(model: FoldModel, scorer: Scorer, rng, cfg: AnnealConfig = AnnealConfig()):
    """Simulated annealing with a soft-to-hard steric ramp and a final quench.

    The returned state is the lowest-scoring state seen during the quench (full clash weight).
    """
    act = model.active_dofs
    nd = len(act)
    steps = cfg.steps or min(cfg.max_steps, 3000 + 40 * nd)
    ramp = max(1, int(steps * (1.0 - cfg.quench_fraction)))
    scorer.clash_scale = cfg.clash_scale_start
    e = scorer()
    stats = {'steps': steps, 'accepted': 0, 'initial_score': e}
    if nd == 0:
        scorer.clash_scale = cfg.clash_scale_end
        e = scorer()
        stats.update(final_score=e, best_score=e)
        return stats
    lib = model.geom.loop_library
    free_res = sorted(model.active_dof_by_res)
    junction_linkers = [L for L in model.tree_order if L.n == 0]
    aform = model.geom.aform
    best = None
    for step in range(steps):
        if step < ramp:
            f = step / ramp
            T = cfg.t_start * (cfg.t_end / cfg.t_start) ** f
            sigma = cfg.sigma_start * (cfg.sigma_end / cfg.sigma_start) ** f
            cs = cfg.clash_scale_start * (cfg.clash_scale_end / cfg.clash_scale_start) ** f
        else:
            T, sigma, cs = cfg.t_end, cfg.sigma_end, cfg.clash_scale_end
        if step % 25 == 0 or step == ramp:
            scorer.clash_scale = cs
            e = scorer()
            if step == ramp:
                best = (e, model.tor.copy(), model.xyz.copy())
        old_xyz = model.xyz.copy()
        old_tor = model.tor.copy()
        u = rng.random()
        touched = []
        if u < cfg.p_residue_jump:
            r = free_res[rng.integers(len(free_res))]
            row = lib[rng.integers(len(lib))]
            for k in model.active_dof_by_res[r]:
                nm = model.dofs[k].name
                target = row[LIB_COLS.index(nm)]
                model.rotate(k, ((target - model.tor[r, TOR[nm]]) + 180) % 360 - 180)
                touched.append(k)
        elif u < cfg.p_residue_jump + cfg.p_junction_aform and junction_linkers:
            L = junction_linkers[rng.integers(len(junction_linkers))]
            for k, d in enumerate(model.dofs):
                if d.linker == L.id and d.name != 'chi':
                    model.rotate(k, ((aform[d.name] - model.tor[d.res, TOR[d.name]]) + 180) % 360 - 180)
                    touched.append(k)
        else:
            k = act[int(rng.integers(nd))]
            model.rotate(k, float(rng.normal(0.0, sigma)))
            touched.append(k)
        if model.dof_loops:
            reclose(model, touched, cfg.reclose_sweeps)
        e_new = scorer()
        if e_new <= e or rng.random() < math.exp(-(e_new - e) / T):
            e = e_new
            stats['accepted'] += 1
            if best is not None and e < best[0]:
                best = (e, model.tor.copy(), model.xyz.copy())
        else:
            model.xyz[:] = old_xyz
            model.tor[:] = old_tor
    if best is not None:
        model.tor[:] = best[1]
        model.xyz[:] = best[2]
    scorer.clash_scale = cfg.clash_scale_end
    stats.update(final_score=e, best_score=best[0] if best else e)
    return stats

