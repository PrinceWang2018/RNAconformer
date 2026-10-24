"""Internal-coordinate all-atom RNA builder (v2 engine).

Every nucleotide is assembled from calibrated rigid pieces (sugar unit per pucker,
planar base, phosphate group) joined by NeRF placement along the backbone
P-O5'-C5'-C4'-C3'-O3'-P with the torsions alpha..zeta and chi. Bond lengths, bond
angles, ring geometry and chirality are therefore correct by construction; only
loop-closure points (see loops.py) can carry a residual, which is measured.

Coordinates are Angstrom. Residues are 0-based internally; output numbering is 1-based.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources

import numpy as np

from .geomutil import place, kabsch, dihedral, rotation_about_axis, to_frame, from_frame

SUGAR_ATOMS = ("C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'")
BASE_ATOMS = {
    'A': ('N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'),
    'G': ('N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'),
    'C': ('N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'),
    'U': ('N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'),
}
GLY = {'A': 'N9', 'G': 'N9', 'C': 'N1', 'U': 'N1'}
CHI_REF = {'A': 'C4', 'G': 'C4', 'C': 'C2', 'U': 'C2'}
PHOSPHATE_ATOMS = ('P', 'OP1', 'OP2', "O5'")
TORSIONS = ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi')
PUCKERS = ('C3endo', 'C2endo')


def atom_order(base):
    return PHOSPHATE_ATOMS + SUGAR_ATOMS + BASE_ATOMS[base]


@dataclass(frozen=True)
class Geometry:
    """Calibrated local geometry loaded from data/geometry_v1.json."""
    raw: dict

    @property
    def version(self):
        return self.raw['version']

    def bond(self, name):
        return self.raw['bonds'][name]['mean']

    def angle(self, name):
        return self.raw['angles'][name]['mean']

    @property
    def aform(self):
        return dict(self.raw.get('aform_torsions_refined') or self.raw['aform_torsions'])

    def sugar(self, pucker):
        return {k: np.array(v) for k, v in self.raw['sugar_templates'][pucker]['atoms'].items()}

    def base(self, b):
        return {k: np.array(v) for k, v in self.raw['base_templates'][b]['atoms'].items()}

    def phosphate(self):
        return {k: np.array(v) for k, v in self.raw['phosphate_template']['atoms'].items()}

    def pair(self, bi, bj):
        key = bi + bj
        t = self.raw['pair_templates'].get(key)
        if t is None:
            return None
        return np.array(t['R']), np.array(t['t'])

    @property
    def loop_library(self):
        lib = self.raw['loop_conformer_library']
        return np.array(lib['rows'])


@lru_cache(maxsize=4)
def load_geometry(name='geometry_v1.json'):
    text = resources.files('rnafold3d').joinpath('data', name).read_text()
    return Geometry(json.loads(text))


class ResidueBuilder:
    """Pre-computes rigid sugar+base units and places them from backbone atoms."""

    def __init__(self, geom: Geometry):
        self.g = geom
        self._sugar = {p: geom.sugar(p) for p in PUCKERS}
        self._base = {b: geom.base(b) for b in 'ACGU'}
        ph = geom.phosphate()
        self._ph_frame = np.array([ph["O3'prev"], ph['P'], ph["O5'"]])
        self._ph_op = np.array([ph['OP1'], ph['OP2']])
        self._ph_local = to_frame(self._ph_op, ph["O3'prev"], ph['P'], ph["O5'"])
        self.b_O3P = geom.bond("O3'-P")
        self.b_PO5 = geom.bond("P-O5'")
        self.b_O5C5 = geom.bond("O5'-C5'")
        self.a_C3O3P = geom.angle("C3'-O3'-P")
        self.a_O3PO5 = geom.angle("O3'-P-O5'")
        self.a_PO5C5 = geom.angle("P-O5'-C5'")
        self.a_O5C5C4 = geom.angle("O5'-C5'-C4'")
        # pucker-dependent values are read from the sugar templates themselves so that
        # NeRF placement and rigid-unit fitting are exactly consistent.
        self.sug = {}
        for p, s in self._sugar.items():
            from .geomutil import angle as ang
            self.sug[p] = {
                'b_C5C4': float(np.linalg.norm(s["C5'"] - s["C4'"])),
                'b_C4C3': float(np.linalg.norm(s["C4'"] - s["C3'"])),
                'b_C3O3': float(np.linalg.norm(s["C3'"] - s["O3'"])),
                'a_C5C4C3': ang(s["C5'"], s["C4'"], s["C3'"]),
                'a_C4C3O3': ang(s["C4'"], s["C3'"], s["O3'"]),
                'delta': dihedral(s["C5'"], s["C4'"], s["C3'"], s["O3'"]),
            }
        self._unit_cache = {}

    # ------------------------------------------------------------------ rigid units
    def unit(self, base, pucker, chi):
        """Sugar + base atoms (dict) in the sugar-template frame for a given chi (deg)."""
        key = (base, pucker, round(float(chi), 1))
        hit = self._unit_cache.get(key)
        if hit is not None:
            return hit
        s = self._sugar[pucker]
        bt = self._base[base]
        gly = GLY[base]
        # align base C1'->N direction with sugar C1'->N, then rotate to target chi
        names = BASE_ATOMS[base]
        bxyz = np.array([bt[n] for n in names])            # base template, N at origin
        c1_b = bt["C1'"]
        v_tpl = -c1_b / np.linalg.norm(c1_b)                # direction C1'->N in template
        v_sug = s['N'] - s["C1'"]
        v_sug /= np.linalg.norm(v_sug)
        r = _rotation_between(v_tpl, v_sug)
        placed = bxyz @ r.T + s['N']
        cur = dihedral(s["O4'"], s["C1'"], s['N'], placed[names.index(CHI_REF[base])])
        rot = rotation_about_axis(v_sug, math.radians(chi - cur))
        placed = (placed - s['N']) @ rot.T + s['N']
        out = {n: s[n] for n in SUGAR_ATOMS}
        for n, xyz in zip(names, placed):
            out[n] = xyz
        if len(self._unit_cache) > 20000:
            self._unit_cache.clear()
        self._unit_cache[key] = out
        return out

    def _unit_local(self, base, pucker, chi):
        key = ('L', base, pucker, round(float(chi), 1))
        hit = self._unit_cache.get(key)
        if hit is None:
            u = self.unit(base, pucker, chi)
            names = list(u)
            xyz = np.array([u[k] for k in names])
            fwd = to_frame(xyz, u["C5'"], u["C4'"], u["C3'"])
            bwd = to_frame(xyz, u["O3'"], u["C3'"], u["C4'"])
            hit = (names, fwd, bwd)
            self._unit_cache[key] = hit
        return hit

    def place_unit_forward(self, base, pucker, chi, c5, c4, c3):
        names, fwd, _ = self._unit_local(base, pucker, chi)
        xyz = from_frame(fwd, c5, c4, c3)
        return dict(zip(names, xyz))

    def place_unit_backward(self, base, pucker, chi, o3, c3, c4):
        names, _, bwd = self._unit_local(base, pucker, chi)
        xyz = from_frame(bwd, o3, c3, c4)
        return dict(zip(names, xyz))

    def phosphate_oxygens(self, o3_prev, p, o5):
        return from_frame(self._ph_local, o3_prev, p, o5)

    # ------------------------------------------------------------------ NeRF steps
    def forward_link(self, prev, eps, zeta, alpha, beta, gamma, pucker_next):
        """From residue atoms `prev` (needs C4', C3', O3') return P, O5', C5', C4', C3' of next."""
        sp = self.sug[pucker_next]
        P = place(prev["C4'"], prev["C3'"], prev["O3'"], self.b_O3P, self.a_C3O3P, eps)
        O5 = place(prev["C3'"], prev["O3'"], P, self.b_PO5, self.a_O3PO5, zeta)
        C5 = place(prev["O3'"], P, O5, self.b_O5C5, self.a_PO5C5, alpha)
        C4 = place(P, O5, C5, sp['b_C5C4'], self.a_O5C5C4, beta)
        C3 = place(O5, C5, C4, sp['b_C4C3'], sp['a_C5C4C3'], gamma)
        return P, O5, C5, C4, C3

    def backward_link(self, nxt, gamma, beta, alpha, zeta, eps, pucker_prev):
        """From residue atoms `nxt` (needs C3', C4', C5') return O5', P (of nxt) and O3', C3', C4' of previous."""
        O5 = place(nxt["C3'"], nxt["C4'"], nxt["C5'"], self.b_O5C5, self.a_O5C5C4, gamma)
        P = place(nxt["C4'"], nxt["C5'"], O5, self.b_PO5, self.a_PO5C5, beta)
        O3 = place(nxt["C5'"], O5, P, self.b_O3P, self.a_O3PO5, alpha)
        sp = self.sug[pucker_prev]
        C3 = place(O5, P, O3, sp['b_C3O3'], self.a_C3O3P, zeta)
        C4 = place(P, O3, C3, sp['b_C4C3'], sp['a_C4C3O3'], eps)
        return O5, P, O3, C3, C4

    def o5_from_sugar(self, res, gamma):
        """5'-terminal O5' placed from the residue's own sugar (no phosphate)."""
        return place(res["C3'"], res["C4'"], res["C5'"], self.b_O5C5, self.a_O5C5C4, gamma)

    def phosphate_from_sugar(self, res, gamma, beta):
        """O5' and P of a residue placed from its own sugar (used for 5' phosphate of strand starts)."""
        O5 = place(res["C3'"], res["C4'"], res["C5'"], self.b_O5C5, self.a_O5C5C4, gamma)
        P = place(res["C4'"], res["C5'"], O5, self.b_PO5, self.a_PO5C5, beta)
        return O5, P


def _rotation_between(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = float(np.dot(a, b))
    if s < 1e-9:
        if c > 0:
            return np.eye(3)
        perp = np.cross(a, [1.0, 0, 0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(a, [0, 1.0, 0])
        return rotation_about_axis(perp, math.pi)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1 - c) / s ** 2)


# ---------------------------------------------------------------------- strands and stems
def build_strand(rb: ResidueBuilder, bases, torsions=None, puckers=None, origin_unit=None):
    """Build a continuous single strand 5'->3' with given per-residue torsions.

    torsions: list of dicts with alpha, beta, gamma, epsilon, zeta, chi (defaults: A-form).
    The first residue's sugar unit is placed in the template frame (or `origin_unit`).
    Returns list of residue atom dicts (O5'/P/OP of the first residue included only as O5').
    """
    n = len(bases)
    aform = rb.g.aform
    torsions = torsions or [aform] * n
    puckers = puckers or ['C3endo'] * n
    res = []
    first = origin_unit or rb.unit(bases[0], puckers[0], torsions[0]['chi'])
    first = dict(first)
    first["O5'"] = rb.o5_from_sugar(first, torsions[0]['gamma'])
    res.append(first)
    for k in range(1, n):
        prev = res[-1]
        tp, tc = torsions[k - 1], torsions[k]
        P, O5, C5, C4, C3 = rb.forward_link(prev, tp['epsilon'], tp['zeta'], tc['alpha'], tc['beta'], tc['gamma'], puckers[k])
        atoms = rb.place_unit_forward(bases[k], puckers[k], tc['chi'], C5, C4, C3)
        op = rb.phosphate_oxygens(prev["O3'"], P, O5)
        atoms.update({'P': P, "O5'": O5, 'OP1': op[0], 'OP2': op[1]})
        res.append(atoms)
    return res


def base_frame(rb: ResidueBuilder, atoms, base):
    """(R, t) mapping the calibrated base template onto the residue's base atoms."""
    tpl = rb._base[base]
    names = BASE_ATOMS[base]
    src = np.array([tpl[n] for n in names])
    dst = np.array([atoms[n] for n in names])
    return kabsch(src, dst)


CENTRAL_HB = {('G', 'C'): ('N1', 'N3'), ('C', 'G'): ('N3', 'N1'), ('A', 'U'): ('N1', 'N3'),
              ('U', 'A'): ('N3', 'N1'), ('G', 'U'): ('N1', 'O2'), ('U', 'G'): ('O2', 'N1')}
# every Watson-Crick / wobble hydrogen bond (donor/acceptor heavy atoms), first base listed first
PAIR_HBONDS = {('G', 'C'): (('N1', 'N3'), ('O6', 'N4'), ('N2', 'O2')),
               ('A', 'U'): (('N1', 'N3'), ('N6', 'O4')),
               ('G', 'U'): (('N1', 'O2'), ('O6', 'N3'))}
_STEM_CACHE = {}


def _pair_targets(rb, sequence, pairs, strand1):
    """Target coordinates (base atoms + C1') of every strand-2 residue from its partner's base frame."""
    targets, generic, names_all = [], [], []
    for k, (i, j) in enumerate(pairs):
        bi, bj = sequence[i], sequence[j]
        tf = rb.g.pair(bi, bj)
        if tf is None:
            tf = rb.g.pair('G', 'C') if bj in 'CU' else rb.g.pair('C', 'G')
            generic.append((i, j))
        R_ij, t_ij = tf
        Ri, ti = base_frame(rb, strand1[k], bi)
        Rj = Ri @ R_ij
        tj = Ri @ t_ij + ti
        tpl = rb._base[bj]
        names = BASE_ATOMS[bj] + ("C1'",)
        targets.append(np.array([tpl[n] for n in names]) @ Rj.T + tj)
        names_all.append(names)
    return targets, names_all, generic


def _stem_quality(model, sequence, pairs):
    """(largest compression of any pair H-bond below 2.9 A, minimum non-bonded cross-residue distance).

    Stretched H-bonds (e.g. in wobble pairs) are not penalised: they cannot cause clashes."""
    worst = 0.0
    for i, j in pairs:
        a, b = sequence[i], sequence[j]
        if (a, b) in PAIR_HBONDS:
            hbs = [(model[i][x], model[j][y]) for x, y in PAIR_HBONDS[(a, b)]]
        elif (b, a) in PAIR_HBONDS:
            hbs = [(model[j][x], model[i][y]) for x, y in PAIR_HBONDS[(b, a)]]
        else:
            hbs = []
        for p, q in hbs:
            worst = max(worst, 2.9 - float(np.linalg.norm(p - q)))
    from scipy.spatial import cKDTree
    keys = [(r, a) for r in model for a in model[r]]
    pts = np.array([model[r][a] for r, a in keys])
    dmin = np.inf
    for a, b in cKDTree(pts).query_pairs(2.75):
        ra, rb_ = keys[a][0], keys[b][0]
        if abs(ra - rb_) <= 1:
            continue
        pol = keys[a][1][0] in 'NO' and keys[b][1][0] in 'NO' and not keys[a][1].startswith('OP')
        d = float(np.linalg.norm(pts[a] - pts[b]))
        if pol and (ra, rb_) in {(i, j) for i, j in pairs} | {(j, i) for i, j in pairs}:
            continue                      # Watson-Crick hydrogen bonds are judged by `worst`
        dmin = min(dmin, d if not pol else d + 0.3)
    return worst, dmin


def build_stem(rb: ResidueBuilder, sequence, pairs, refine=True):
    """Ideal A-form duplex for consecutive stacked pairs [(i, j), (i+1, j-1), ...] (0-based).

    Strand 1 (i..) and strand 2 (..j) are each built from calibrated A-form torsions, so both
    strands have exact covalent geometry. Strand 2 is rigidly fitted so every base matches the
    calibrated Watson-Crick/wobble pair geometry relative to its partner. If a pair is distorted
    (any pair H-bond compressed below 2.4 A) or two residues come closer than 2.5 A (e.g. next to G-U
    wobbles or at stem ends), strand-2 torsions are refined by regularised least squares
    (pair-geometry residuals + deviation from A-form in units of the calibrated spread).
    Results are cached per (stem sequence) since they do not depend on the rest of the RNA.
    """
    key = (''.join(sequence[i] for i, _ in pairs), ''.join(sequence[j] for _, j in pairs), refine)
    hit = _STEM_CACHE.get(key)
    if hit is not None:
        base_i, base_j = pairs[0]
        loc, rep, s1_off, s2_off = hit
        out = {}
        for d, atoms in loc.items():
            r = (base_i + d[1]) if d[0] == 1 else (base_j - d[1])
            out[r] = {k: v.copy() for k, v in atoms.items()}
        return out, dict(rep)
    L = len(pairs)
    s1 = [p[0] for p in pairs]
    s2 = [p[1] for p in reversed(pairs)]           # 5'->3' order of strand 2
    aform = rb.g.aform
    strand1 = build_strand(rb, [sequence[i] for i in s1])
    targets, names_all, generic = _pair_targets(rb, sequence, pairs, strand1)
    dst = np.vstack(targets)

    def assemble(tors):
        strand2 = build_strand(rb, [sequence[j] for j in s2], torsions=tors)
        src = np.vstack([np.array([strand2[L - 1 - k][n] for n in names_all[k]]) for k in range(L)])
        r, t = kabsch(src, dst)
        fit = src @ r.T + t
        return strand2, r, t, fit

    base_tors = [dict(aform) for _ in range(L)]
    strand2, r, t, fit = assemble(base_tors)

    def to_model(strand2, r, t):
        out = {}
        for idx, atoms in zip(s1, strand1):
            out[idx] = atoms
        for idx, atoms in zip(s2, strand2):
            out[idx] = {k: v @ r.T + t for k, v in atoms.items()}
        return out
    out = to_model(strand2, r, t)
    fit_rmsd = float(np.sqrt(np.mean(np.sum((fit - dst) ** 2, axis=1))))
    worst, dmin = _stem_quality(out, sequence, pairs)
    report = {'pair_fit_rmsd': fit_rmsd, 'noncanonical_pairs_approximated': generic, 'refined': False,
              'max_hbond_compression': round(worst, 3), 'min_contact': None if not np.isfinite(dmin) else round(dmin, 3)}
    if refine and L >= 2 and (worst > 0.5 or dmin < 2.5):
        from scipy.optimize import least_squares
        names = ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi')
        sd = np.array([max(rb.g.raw['aform_torsion_sd'][n], 4.0) for n in names])

        def tors_of(x):
            x = x.reshape(L, 6)
            return [{**aform, **{n: aform[n] + x[k, q] for q, n in enumerate(names)}} for k in range(L)]

        def resid(x):
            _, _, _, f = assemble(tors_of(x))
            return np.concatenate([(f - dst).ravel(), 0.5 * (x.reshape(L, 6) / sd).ravel()])
        sol = least_squares(resid, np.zeros(6 * L), method='trf', max_nfev=25, x_scale=5.0, diff_step=1e-3)
        strand2b, rb2, tb2, fitb = assemble(tors_of(sol.x))
        out_b = to_model(strand2b, rb2, tb2)
        worst_b, dmin_b = _stem_quality(out_b, sequence, pairs)
        if (worst_b <= 0.5 and dmin_b >= 2.5) or (min(dmin_b, 2.5) >= min(dmin, 2.5) and worst_b < worst):
            out = out_b
            report.update({'refined': True, 'pair_fit_rmsd': float(np.sqrt(np.mean(np.sum((fitb - dst) ** 2, axis=1)))),
                           'max_hbond_compression': round(worst_b, 3),
                           'min_contact': None if not np.isfinite(dmin_b) else round(dmin_b, 3),
                           'max_torsion_shift_deg': round(float(np.max(np.abs(sol.x))), 1)})
    base_i, base_j = pairs[0]
    loc = {}
    for r_, atoms in out.items():
        d = (1, r_ - base_i) if r_ in s1 else (2, base_j - r_)
        loc[d] = {k: v.copy() for k, v in atoms.items()}
    _STEM_CACHE[key] = (loc, dict(report), None, None)
    return out, report
