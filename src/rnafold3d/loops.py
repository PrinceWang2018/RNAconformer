"""All-atom loop closure by cyclic coordinate descent (Canutescu & Dunbrack 2003) in
RNA backbone torsion space.

A closure linker runs from residue a (fixed, in one rigid body) through the unpaired
residues to residue b (fixed, possibly in another body). The free torsions are
epsilon/zeta of a, alpha/beta/gamma/epsilon/zeta of every loop residue and
alpha/beta/gamma of b. CCD moves a *copy* of b's C5'/C4'/C3' onto the fixed atoms;
bond lengths and angles stay exact everywhere, and the remaining gap at the
O5'(b)-C5'(b) junction is reported.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.spatial import cKDTree

from .fold import TOR, LIB_COLS
from .geomutil import place, dihedral


class TorsionPrior:
    """-ln p(torsion) from a von Mises kernel density of calibrated loop residues (kT units)."""

    def __init__(self, geom, kappa=12.0, cap=6.0):
        lib = geom.loop_library
        grid = np.radians(np.arange(-180, 180))
        self.table = {}
        for name in ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi'):
            vals = np.radians(lib[:, LIB_COLS.index(name)])
            dens = np.exp(kappa * np.cos(grid[:, None] - vals[None, :])).sum(axis=1)
            e = -np.log(dens / dens.max() + 1e-12)
            self.table[name] = np.minimum(e, cap)

    def __call__(self, name, deg):
        k = (np.asarray(deg, float) + 180.0).astype(int) % 360
        return self.table[name][k]


def _nerf(a, b, c, bond, ca, sa, ct, st):
    """Tuple NeRF: d from a, b, c with |cd| = bond, cos/sin of angle(b,c,d) and of the dihedral."""
    ux, uy, uz = c[0] - b[0], c[1] - b[1], c[2] - b[2]
    lu = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux, uy, uz = ux / lu, uy / lu, uz / lu
    vx, vy, vz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    nx, ny, nz = vy * uz - vz * uy, vz * ux - vx * uz, vx * uy - vy * ux
    ln = math.sqrt(nx * nx + ny * ny + nz * nz)
    nx, ny, nz = nx / ln, ny / ln, nz / ln
    mx, my, mz = ny * uz - nz * uy, nz * ux - nx * uz, nx * uy - ny * ux
    d0 = -bond * ca
    d1 = bond * sa * ct
    d2 = bond * sa * st
    return (c[0] + d0 * ux + d1 * mx + d2 * nx, c[1] + d0 * uy + d1 * my + d2 * ny, c[2] + d0 * uz + d1 * mz + d2 * nz)


def _consts(rb):
    k = getattr(rb, '_ccd_consts', None)
    if k is None:
        r = math.radians
        k = {'O3P': (rb.b_O3P, math.cos(r(rb.a_C3O3P)), math.sin(r(rb.a_C3O3P))),
             'PO5': (rb.b_PO5, math.cos(r(rb.a_O3PO5)), math.sin(r(rb.a_O3PO5))),
             'O5C5': (rb.b_O5C5, math.cos(r(rb.a_PO5C5)), math.sin(r(rb.a_PO5C5))),
             'puck': {}}
        for p, sp in rb.sug.items():
            k['puck'][p] = {'C5C4': (sp['b_C5C4'], math.cos(r(rb.a_O5C5C4)), math.sin(r(rb.a_O5C5C4))),
                            'C4C3': (sp['b_C4C3'], math.cos(r(sp['a_C5C4C3'])), math.sin(r(sp['a_C5C4C3']))),
                            'C3O3': (sp['b_C3O3'], math.cos(r(sp['a_C4C3O3'])), math.sin(r(sp['a_C4C3O3'])),
                                     math.cos(r(sp['delta'])), math.sin(r(sp['delta'])))}
        rb._ccd_consts = k
    return k


def _mainchain(model, L, tor, puck, close_at='p'):
    """Main-chain coordinates (list of tuples) and DOF axes for a forward closure linker.

    close_at='p' : targets are P, O5', C5' of residue b (b's beta/gamma fixed); the chain
                   is closed at the O3'(b-1)-P(b) bond (DOFs: eps/zeta of a, loop torsions, alpha of b).
    close_at='c5': targets are C5', C4', C3' of b (alpha/beta/gamma of b free) - used for 0-nt linkers.
    """
    K = _consts(model.rb)
    I = model.idx
    x = model.xyz
    a = L.a
    c4p, c3p, o3p = (tuple(x[I[(a, n)]].tolist()) for n in ("C4'", "C3'", "O3'"))
    chain = [c4p, c3p, o3p]
    dofs = [(a, 'epsilon', 1, 2), (a, 'zeta', 2, 3)]
    prev_r = a
    rad = math.radians
    for m in list(L.loop) + [L.b]:
        e, z = rad(tor[prev_r, 3]), rad(tor[prev_r, 4])
        al, be, ga = rad(tor[m, 0]), rad(tor[m, 1]), rad(tor[m, 2])
        kp = K['puck'][puck[m]]
        P = _nerf(c4p, c3p, o3p, *K['O3P'], math.cos(e), math.sin(e))
        O5 = _nerf(c3p, o3p, P, *K['PO5'], math.cos(z), math.sin(z))
        C5 = _nerf(o3p, P, O5, *K['O5C5'], math.cos(al), math.sin(al))
        o = len(chain)
        if m == L.b and close_at == 'p':
            chain += [P, O5, C5]
            dofs.append((m, 'alpha', o, o + 1))
            break
        C4 = _nerf(P, O5, C5, *kp['C5C4'], math.cos(be), math.sin(be))
        C3 = _nerf(O5, C5, C4, *kp['C4C3'], math.cos(ga), math.sin(ga))
        chain += [P, O5, C5, C4, C3]
        dofs += [(m, 'alpha', o, o + 1), (m, 'beta', o + 1, o + 2), (m, 'gamma', o + 2, o + 3)]
        if m == L.b:
            break
        b3, ca3, sa3, cd, sd = kp['C3O3']
        O3 = _nerf(C5, C4, C3, b3, ca3, sa3, cd, sd)
        chain.append(O3)
        dofs += [(m, 'epsilon', o + 4, o + 5), (m, 'zeta', o + 5, o + 6)]
        c4p, c3p, o3p = C4, C3, O3
        prev_r = m
    return chain, dofs


def _targets(model, L, tor, close_at):
    I = model.idx
    rb = model.rb
    resb = {n: model.xyz[I[(L.b, n)]] for n in ("C3'", "C4'", "C5'")}
    if close_at == 'p':
        O5, P = rb.phosphate_from_sugar(resb, tor[L.b, 2], tor[L.b, 1])
        return np.array([P, O5, resb["C5'"]])
    return np.array([resb["C5'"], resb["C4'"], resb["C3'"]])


def ccd_close(model, L, tor, puck, max_sweeps=150, tol=0.03, close_at='p', frozen=()):
    """CCD on a copy of torsions; returns (converged, rmsd, tor).

    Backward sweeps (end effector -> start): rotations applied earlier in a sweep only move
    atoms downstream of later axes, so each step needs just the three end-effector atoms;
    the chain is rebuilt from torsions (NeRF) once per sweep.
    """
    tor = tor.copy()
    target = _targets(model, L, tor, close_at)
    T = target.tolist()
    rmsd = float('inf')
    for sweep in range(max_sweeps + 1):
        chain, dofs = _mainchain(model, L, tor, puck, close_at)
        if frozen:
            dofs = [d for d in dofs if (d[0], d[1]) not in frozen]
        E = [list(p) for p in chain[-3:]]
        rmsd = math.sqrt(sum((E[q][0] - T[q][0]) ** 2 + (E[q][1] - T[q][1]) ** 2 + (E[q][2] - T[q][2]) ** 2
                             for q in range(3)) / 3.0)
        if rmsd <= tol or sweep == max_sweeps:
            break
        pts = chain
        for res, name, u, v in reversed(dofs):
            ox, oy, oz = pts[v]
            ax, ay, az = ox - pts[u][0], oy - pts[u][1], oz - pts[u][2]
            la = math.sqrt(ax * ax + ay * ay + az * az)
            ax, ay, az = ax / la, ay / la, az / la
            num = den = 0.0
            for q in range(3):
                rx, ry, rz = E[q][0] - ox, E[q][1] - oy, E[q][2] - oz
                pr = rx * ax + ry * ay + rz * az
                px, py, pz = rx - pr * ax, ry - pr * ay, rz - pr * az
                fx, fy, fz = T[q][0] - ox, T[q][1] - oy, T[q][2] - oz
                num += fx * (ay * pz - az * py) + fy * (az * px - ax * pz) + fz * (ax * py - ay * px)
                den += fx * px + fy * py + fz * pz
            theta = math.atan2(num, den)
            if abs(theta) < 1e-7:
                continue
            c, sn = math.cos(theta), math.sin(theta)
            for q in range(3):
                rx, ry, rz = E[q][0] - ox, E[q][1] - oy, E[q][2] - oz
                dot = rx * ax + ry * ay + rz * az
                cx_, cy_, cz_ = ay * rz - az * ry, az * rx - ax * rz, ax * ry - ay * rx
                E[q] = [ox + rx * c + cx_ * sn + ax * dot * (1 - c),
                        oy + ry * c + cy_ * sn + ay * dot * (1 - c),
                        oz + rz * c + cz_ * sn + az * dot * (1 - c)]
            k = TOR[name]
            tor[res, k] = (tor[res, k] + math.degrees(theta) + 180.0) % 360.0 - 180.0
    return rmsd <= tol, rmsd, tor


def loop_torsion_residues(L):
    out = [(L.a, 'epsilon'), (L.a, 'zeta')]
    for m in L.loop:
        out += [(m, nm) for nm in ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi')]
    out += [(L.b, 'alpha'), (L.b, 'beta'), (L.b, 'gamma')]
    return out


def _loop_clashes(model, L, cutoff=2.8):
    """Non-bonded heavy-atom overlaps involving the linker's atoms (bond-graph aware)."""
    subset = np.zeros(model.n_atoms, bool)
    subset[L.atoms] = True
    i, j, ov = model.atom_clashes(cutoff=cutoff, subset=subset)
    return int(len(ov)), float(np.sum(ov ** 2))


def _residue_ids(model, m):
    from .build import atom_order
    return np.array([model.idx[(m, nm)] for nm in atom_order(model.seq[m]) if not (m == 0 and nm in ('P', 'OP1', 'OP2'))])


def grow_loop(model, L, rng, tor, puck, reach, tries=25, cutoff=2.8):
    """Sequentially place loop residues from library conformers, avoiding clashes and keeping
    the far anchor reachable. Modifies tor/puck copies in place; leaves model atoms of the
    linker absent. Returns number of residues placed without a clash."""
    from scipy.spatial import cKDTree
    rb = model.rb
    I = model.idx
    lib = model.geom.loop_library
    model.present[L.atoms] = False
    stat = np.where(model.present)[0]
    tree = cKDTree(model.xyz[stat]) if len(stat) else None
    aform = model.geom.aform
    tor[L.b, 1] = aform['beta']
    tor[L.b, 2] = aform['gamma']
    resb = {n: model.xyz[I[(L.b, n)]] for n in ("C3'", "C4'", "C5'")}
    _, P_b = rb.phosphate_from_sugar(resb, tor[L.b, 2], tor[L.b, 1])
    prev_r = L.a
    prev = {n: model.xyz[I[(L.a, n)]] for n in ("C4'", "C3'", "O3'")}
    grown = []
    clean = 0
    nloop = len(L.loop)
    for k, m in enumerate(L.loop):
        remaining = nloop - k - 1          # residues still to place after m
        lim = reach['p995'][min(remaining, len(reach['p995']) - 1)] + 1.5
        best = None
        for t in range(tries):
            row = lib[rng.integers(len(lib))]
            e_prev, z_prev = (row[4], row[5]) if t else (tor[prev_r, 3], tor[prev_r, 4])
            prow = lib[rng.integers(len(lib))]
            e_prev, z_prev = prow[4], prow[5]
            pk = 'C2endo' if row[3] > 115 else 'C3endo'
            P, O5, C5, C4, C3 = rb.forward_link(prev, e_prev, z_prev, row[0], row[1], row[2], pk)
            unit = rb.place_unit_forward(model.seq[m], pk, row[6], C5, C4, C3)
            op = rb.phosphate_oxygens(prev["O3'"], P, O5)
            atoms = dict(unit)
            atoms.update({'P': P, "O5'": O5, 'OP1': op[0], 'OP2': op[1]})
            ids = np.array([I[(m, nm)] for nm in atoms])
            pts = np.array(list(atoms.values()))
            # distance from this residue's O3' to the next P target region
            d_reach = float(np.linalg.norm(atoms["O3'"] - P_b))
            bad = max(0.0, d_reach - lim)
            n_cl = 0
            if tree is not None:
                hits = tree.query_ball_point(pts, cutoff)
                ii, jj = [], []
                for q, h in enumerate(hits):
                    for hh in h:
                        ii.append(ids[q])
                        jj.append(stat[hh])
                if ii:
                    ii, jj = np.array(ii), np.array(jj)
                    keep = model.nonbonded_mask(ii, jj)
                    n_cl += int(keep.sum())
            if grown:
                gids = np.concatenate(grown)
                d = np.linalg.norm(pts[:, None, :] - model.xyz[gids][None, :, :], axis=2)
                qi, qj = np.where(d < cutoff)
                if len(qi):
                    keep = model.nonbonded_mask(ids[qi], gids[qj])
                    n_cl += int(keep.sum())
            score = n_cl + 3.0 * bad
            if best is None or score < best[0]:
                best = (score, row, e_prev, z_prev, pk, atoms, ids)
            if score == 0:
                break
        score, row, e_prev, z_prev, pk, atoms, ids = best
        clean += score == 0
        tor[prev_r, 3], tor[prev_r, 4] = e_prev, z_prev
        tor[m, 0], tor[m, 1], tor[m, 2], tor[m, 3], tor[m, 4], tor[m, 5] = row[0], row[1], row[2], row[4], row[5], row[6]
        puck[m] = pk
        for nm, v in atoms.items():
            model.xyz[I[(m, nm)]] = v
        grown.append(ids)
        prev, prev_r = atoms, m
    model.present[L.atoms] = False
    return clean


def _apply(model, L, tor, puck):
    """Build the linker from torsions; return the junction error (A): RMS deviation of the
    O5'(b)-C5'(b) bond length and of the two junction angles expressed as distances."""
    model.tor[:] = tor
    model.puck[:] = puck
    model._build_linker(L)
    model.present[model.never] = False
    I = model.idx
    rb = model.rb
    o5 = model.xyz[I[(L.b, "O5'")]]
    p = model.xyz[I[(L.b, 'P')]]
    c5 = model.xyz[I[(L.b, "C5'")]]
    c4 = model.xyz[I[(L.b, "C4'")]]
    d_bond = np.linalg.norm(o5 - c5) - rb.b_O5C5
    # angle deviations converted to displacement at bond length ~1.5 A
    from .geomutil import angle as _ang
    da1 = math.radians(_ang(p, o5, c5) - rb.a_PO5C5) * 1.5
    da2 = math.radians(_ang(o5, c5, c4) - rb.a_O5C5C4) * 1.5
    return float(math.sqrt((d_bond ** 2 + da1 ** 2 + da2 ** 2) / 3.0))


def _loop_score(model, L, prior):
    subset = np.zeros(model.n_atoms, bool)
    subset[L.atoms] = True
    i, j, ov = model.atom_clashes(cutoff=2.8, subset=subset)
    d = np.linalg.norm(model.xyz[i] - model.xyz[j], axis=1) if len(i) else np.zeros(0)
    severe = int(np.sum(d < 2.3))
    e_tor = float(sum(prior(nm, model.tor[r, TOR[nm]]) for r, nm in loop_torsion_residues(L)))
    return 30.0 * float(np.sum(ov ** 2)) + 2.0 * len(ov) + 40.0 * severe, e_tor, int(len(ov))


def relax_loop(model, L, rng, prior, steps=150, sigma=20.0, tol=0.03, close_at='p'):
    """Monte Carlo in loop torsion space with CCD re-closure after every move."""
    frozen = {(L.b, 'beta'), (L.b, 'gamma')} if close_at == 'p' else set()
    free = [(r, nm) for r, nm in loop_torsion_residues(L) if (r, nm) not in frozen and nm != 'chi'
            and not (close_at == 'p' and r == L.b and nm != 'alpha')]
    tor0, puck0 = model.tor.copy(), model.puck.copy()
    cl, et, n_cl = _loop_score(model, L, prior)
    cur = cl + et
    best = (cur, tor0.copy(), n_cl)
    for step in range(steps):
        T = 3.0 * (0.2 / 3.0) ** (step / max(1, steps - 1))
        tor = model.tor.copy()
        r, nm = free[rng.integers(len(free))]
        tor[r, TOR[nm]] = (tor[r, TOR[nm]] + rng.normal(0, sigma) + 180) % 360 - 180
        ok, rmsd, tor = ccd_close(model, L, tor, model.puck, max_sweeps=40, tol=tol, close_at=close_at, frozen=frozen)
        if not ok:
            continue
        saved = model.tor.copy()
        gap = _apply(model, L, tor, model.puck)
        cl, et, n_cl = _loop_score(model, L, prior)
        new = cl + et + 100.0 * gap
        if new <= cur or rng.random() < math.exp(-(new - cur) / T):
            cur = new
            if new < best[0]:
                best = (new, model.tor.copy(), n_cl)
            if n_cl == 0 and step >= 30:
                break
        else:
            _apply(model, L, saved, model.puck)
    _apply(model, L, best[1], model.puck)
    return best[2]


def close_linker(model, L, rng, prior: TorsionPrior, attempts=24, keep=4, tol=0.03, relax_steps=150):
    """Close linker L in place (writes model.tor/puck/xyz/present). Returns a report dict."""
    reach = model.geom.raw['loop_reach']
    close_at = getattr(L, 'close_at', 'p' if L.n >= 1 else 'c5')
    frozen = {(L.b, 'beta'), (L.b, 'gamma')} if close_at == 'p' else set()
    best = None
    tried = converged = 0
    for attempt in range(attempts):
        tried += 1
        tor = model.tor.copy()
        puck = model.puck.copy()
        if L.n >= 1:
            grow_loop(model, L, rng, tor, puck, reach)
        else:
            lib = model.geom.loop_library
            row = lib[rng.integers(len(lib))]
            for r, nm in ((L.a, 'epsilon'), (L.a, 'zeta'), (L.b, 'alpha'), (L.b, 'beta'), (L.b, 'gamma')):
                tor[r, TOR[nm]] = row[LIB_COLS.index(nm)]
        ok, rmsd, tor = ccd_close(model, L, tor, puck, tol=tol, close_at=close_at, frozen=frozen)
        if not ok:
            continue
        converged += 1
        saved = (model.tor.copy(), model.puck.copy())
        gap = _apply(model, L, tor, puck)
        clash_score, e_tor, n_cl = _loop_score(model, L, prior)
        score = clash_score + e_tor + 100.0 * gap
        cand = {'score': score, 'clashes': n_cl, 'torsion_energy': e_tor, 'gap_rmsd': gap,
                'tor': model.tor.copy(), 'puck': model.puck.copy()}
        model.tor[:], model.puck[:] = saved
        if best is None or score < best['score']:
            best = cand
        if converged >= keep and best['clashes'] == 0:
            break
    report = {'a': L.a + 1, 'b': L.b + 1, 'n': L.n, 'kind': L.kind, 'attempts': tried,
              'converged_attempts': converged, 'closure_point': "O3'-P" if close_at == 'p' else "O5'-C5'"}
    if best is None:
        model.present[L.atoms] = False
        report['closed'] = False
        return report
    gap = _apply(model, L, best['tor'], best['puck'])
    relaxed = False
    if best['clashes'] and relax_steps:
        relax_loop(model, L, rng, prior, steps=relax_steps, tol=tol, close_at=close_at)
        relaxed = True
    gap = _apply(model, L, model.tor, model.puck)
    _, e_tor, n_cl = _loop_score(model, L, prior)
    report.update({'closed': True, 'gap_rmsd': round(gap, 4), 'clashes': n_cl,
                   'torsion_energy': round(e_tor, 2), 'relaxed': relaxed})
    return report


def refine_chi(model, L, prior, choices=(-160.0, -120.0, -75.0, 60.0)):
    """Greedy glycosidic-angle choice for loop residues to reduce atom overlaps."""
    from .build import BASE_ATOMS, GLY
    changed = 0
    for m in L.loop:
        k = next((i for i, d in enumerate(model.dofs) if d.res == m and d.name == 'chi'), None)
        b = model.seq[m]
        base_ids = np.array([model.idx[(m, nm)] for nm in BASE_ATOMS[b]])
        subset = np.zeros(model.n_atoms, bool)
        subset[base_ids] = True

        def score():
            _, _, ov = model.atom_clashes(subset=subset)
            return 30.0 * float(np.sum(ov ** 2)) + float(prior('chi', model.tor[m, TOR['chi']]))
        best_val = model.tor[m, TOR['chi']]
        best_score = score()
        for target in choices:
            delta = ((target - model.tor[m, TOR['chi']]) + 180) % 360 - 180
            _rotate_base(model, m, base_ids, k, delta)
            sc = score()
            if sc < best_score - 1e-9:
                best_score, best_val = sc, model.tor[m, TOR['chi']]
        back = ((best_val - model.tor[m, TOR['chi']]) + 180) % 360 - 180
        _rotate_base(model, m, base_ids, k, back)
        changed += 1
    return changed


def _rotate_base(model, m, base_ids, k, delta):
    from .fold import TOR as _T
    from .geomutil import rotation_about_axis
    if k is not None:
        model.rotate(k, delta)
        return
    from .build import GLY
    b = model.seq[m]
    u = model.xyz[model.idx[(m, "C1'")]]
    v = model.xyz[model.idx[(m, GLY[b])]]
    R = rotation_about_axis(v - u, math.radians(delta))
    model.xyz[base_ids] = (model.xyz[base_ids] - v) @ R.T + v
    model.tor[m, _T['chi']] = (model.tor[m, _T['chi']] + delta + 180) % 360 - 180


def final_relief(model, rng, prior, relax_steps=120):
    """Post-closure clash relief with every atom present.

    Only moves that cannot break other loops are used: closure loops are relaxed with CCD
    re-closure, tail torsions are rotated (tails are chain ends), and loop glycosidic angles
    are re-chosen. Tree-linker backbone torsions are not touched here (they carry stems).
    """
    report = {'relaxed_loops': 0, 'tail_steps': 0}
    for L in model.linkers:
        if L.kind not in ('closure', 'hairpin'):
            continue
        n_cl, _ = _loop_clashes(model, L)
        if n_cl:
            close_at = getattr(L, 'close_at', 'p' if L.n >= 1 else 'c5')
            relax_loop(model, L, rng, prior, steps=relax_steps, close_at=close_at)
            report['relaxed_loops'] += 1
    tails = [L for L in model.linkers if L.kind in ('tail5', 'tail3')]
    for L in tails:
        ks = [k for k, d in enumerate(model.dofs) if d.linker == L.id]
        if not ks:
            continue
        sub = np.zeros(model.n_atoms, bool)
        sub[L.atoms] = True

        def tail_score():
            _, _, ov = model.atom_clashes(subset=sub)
            e = sum(float(prior(model.dofs[k].name, model.tor[model.dofs[k].res, TOR[model.dofs[k].name]])) for k in ks)
            return 30.0 * float(np.sum(ov ** 2)) + 2.0 * len(ov) + e
        cur = tail_score()
        for step in range(200):
            k = ks[rng.integers(len(ks))]
            delta = float(rng.normal(0, 25.0))
            old_xyz, old_tor = model.xyz.copy(), model.tor.copy()
            model.rotate(k, delta)
            new = tail_score()
            if new <= cur:
                cur = new
            else:
                model.xyz[:], model.tor[:] = old_xyz, old_tor
            report['tail_steps'] += 1
    for L in model.linkers:
        if L.loop:
            refine_chi(model, L, prior)
    return report


def finish_closure(model, L, rng, prior, tol=0.03):
    """Close an explicit (already sampled) loop exactly with CCD from its current torsions;
    fall back to regrowth (close_linker) if CCD does not converge or leaves severe clashes."""
    close_at = getattr(L, 'close_at', 'p')
    frozen = {(L.b, 'beta'), (L.b, 'gamma')} if close_at == 'p' else set()
    start_cb = model.chainbreak(L)
    ok, rmsd, tor = ccd_close(model, L, model.tor, model.puck, tol=tol, close_at=close_at, frozen=frozen)
    report = {'a': L.a + 1, 'b': L.b + 1, 'n': L.n, 'kind': L.kind,
              'closure_point': "O3'-P" if close_at == 'p' else "O5'-C5'",
              'chainbreak_before_A': round(math.sqrt(start_cb / 3.0), 3)}
    if ok:
        saved = (model.tor.copy(), model.puck.copy())
        gap = _apply(model, L, tor, model.puck)
        _, e_tor, n_cl = _loop_score(model, L, prior)
        subset = np.zeros(model.n_atoms, bool)
        subset[L.atoms] = True
        i, j, ov = model.atom_clashes(cutoff=2.2, polar_cutoff=2.2, subset=subset)
        if len(ov) == 0:
            report.update({'closed': True, 'method': 'ccd_from_search', 'gap_rmsd': round(gap, 4), 'clashes': n_cl,
                           'torsion_energy': round(e_tor, 2), 'relaxed': False})
            return report
        relax_loop(model, L, rng, prior, steps=120, tol=tol, close_at=close_at)
        gap = _apply(model, L, model.tor, model.puck)
        i, j, ov = model.atom_clashes(cutoff=2.2, polar_cutoff=2.2, subset=subset)
        if len(ov) == 0:
            _, e_tor, n_cl = _loop_score(model, L, prior)
            report.update({'closed': True, 'method': 'ccd_from_search+relax', 'gap_rmsd': round(gap, 4),
                           'clashes': n_cl, 'torsion_energy': round(e_tor, 2), 'relaxed': True})
            return report
        model.tor[:], model.puck[:] = saved
        _apply(model, L, model.tor, model.puck)
    rep = close_linker(model, L, rng, prior)
    rep['method'] = 'regrow'
    rep['chainbreak_before_A'] = report['chainbreak_before_A']
    return rep
