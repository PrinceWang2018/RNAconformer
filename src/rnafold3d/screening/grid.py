"""The receptor on a grid: the scoring function as a lookup.

Docking a compound means scoring tens of thousands of poses, and every pose costs a sum over receptor atoms.
Since every term depends only on the receptor and on what *kind* of atom the compound puts at a point, the
sum can be done once per point of a grid and looked up afterwards. That is the whole reason a screen of
100,000 compounds is possible on a laptop, and the reason the pocket has to be a pocket: the grid covers a
box around the site, not the molecule.

What is stored, per grid point:

* per radius of compound atom (a handful of values: 1.7 for oxygen, 1.9 for carbon …), the shape and
  contact terms — gauss1, gauss2, repulsion, hydrophobic, the two hydrogen-bond directions, metal;
* once for all atoms, the electrostatic potential (multiplied by the atom's charge at lookup) and the
  desolvation field (multiplied by |charge|);
* once for all atoms, the six independent entries of the symmetric matrix M(x) = Σ_b w_b n_b n_bᵀ over the
  nearby base planes. The stacking term of a compound atom is n·M(x)·n with n the normal of its own
  aromatic ring, which is *exactly* the sum over bases of w_b cos²(angle between the planes) — the
  direction survives the gridding instead of being averaged away.

Values between grid points come from trilinear interpolation; a pose that leaves the box is given a large
penalty, which is what keeps a search inside the pocket it was asked about.
"""
from __future__ import annotations

import numpy as np

from . import score as S

SPACING = 0.375                  # Å, the usual docking grid step
PAD = 4.0                        # Å, how far the box reaches beyond the pocket points
OUTSIDE = 10.0                   # score penalty per atom outside the box
SHAPE_CHANNELS = ('gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'hbond_donor', 'hbond_acceptor', 'metal')


class Grid:
    """Maps of the receptor terms over a box, and the lookup that scores a pose.

    Build it from a receptor and a box (`Grid.for_pocket` takes a pocket from `pockets.find`), then call
    `energy(xyz, lig)` with the coordinates of a pose. The result is the same number `score.score` gives,
    within the interpolation error — `tests` check that on real complexes.
    """

    def __init__(self, rec, lo, hi, spacing=SPACING, radii=(1.5, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2), cutoff=S.CUTOFF):
        self.spacing = float(spacing)
        self.lo = np.asarray(lo, float)
        self.radii = tuple(sorted(set(float(r) for r in radii)))
        self.cutoff = float(cutoff)
        shape = np.maximum(np.ceil((np.asarray(hi, float) - self.lo) / self.spacing).astype(int) + 1, 2)
        self.shape = tuple(int(x) for x in shape)
        self.hi = self.lo + (np.array(self.shape) - 1) * self.spacing
        self.maps = {}
        self._cache = {}
        self._build(rec)

    # ------------------------------------------------------------------------------------- building
    def _axes(self):
        return [self.lo[k] + np.arange(self.shape[k]) * self.spacing for k in range(3)]

    def _build(self, rec):
        n = self.shape
        chan = {(r, c): np.zeros(n, np.float32) for r in self.radii for c in SHAPE_CHANNELS}
        elec = np.zeros(n, np.float32)
        desolv = np.zeros(n, np.float32)
        stack = np.zeros(n + (6,), np.float32)
        ax = self._axes()
        near = rec.within((self.lo + self.hi) / 2, float(np.linalg.norm(self.hi - self.lo) / 2 + self.cutoff + 3))
        for a in near:
            p = rec.xyz[a]
            # every term is zero at and beyond the cutoff (the smooth cut below takes the last 2 Å down to
            # it, and the contact terms die out well before), so the stamp only has to reach that far
            sl, d = self._sub(ax, p, self.cutoff)
            if d is None:
                continue
            w = S._slope(d, self.cutoff - 2.0, self.cutoff)          # smooth cut at the box of influence
            if rec.charge[a]:
                elec[sl] += (S.COULOMB * rec.charge[a] / (4.0 * np.maximum(d, S.EPS_MIN) ** 2) * w).astype(np.float32)
            desolv[sl] += ((4 / 3 * np.pi * rec.radius[a] ** 3) * np.exp(-(d / 3.5) ** 2) / 100.0).astype(np.float32)
            if rec.metal[a]:
                met = np.exp(-((d - 2.1) / 0.5) ** 2) * w
            for r in self.radii:
                sd = d - r - rec.radius[a]                           # surface distance for this atom radius
                chan[(r, 'gauss1')][sl] += (np.exp(-(sd / 0.5) ** 2) * w).astype(np.float32)
                chan[(r, 'gauss2')][sl] += (np.exp(-((sd - 3.0) / 2.0) ** 2) * w).astype(np.float32)
                chan[(r, 'repulsion')][sl] += (np.where(sd < 0, sd * sd, 0.0)).astype(np.float32)
                if rec.hydrophobic[a]:
                    chan[(r, 'hydrophobic')][sl] += S._slope(sd, 0.5, 1.5).astype(np.float32)
                if rec.acceptor[a]:
                    chan[(r, 'hbond_donor')][sl] += S._slope(sd, -0.7, 0.0).astype(np.float32)
                if rec.donor[a]:
                    chan[(r, 'hbond_acceptor')][sl] += S._slope(sd, -0.7, 0.0).astype(np.float32)
                if rec.metal[a]:
                    chan[(r, 'metal')][sl] += met.astype(np.float32)
        for b in range(len(rec.ring_centre)):
            c, nb = rec.ring_centre[b], rec.ring_normal[b]
            sl, d = self._sub(ax, c, 9.0)
            if d is None:
                continue
            v = np.stack(np.meshgrid(ax[0][sl[0]] - c[0], ax[1][sl[1]] - c[1], ax[2][sl[2]] - c[2],
                                     indexing='ij'), axis=-1)
            perp = np.abs(v @ nb)
            lateral = np.sqrt(np.maximum(0.0, (v * v).sum(-1) - perp ** 2))
            wgt = np.exp(-((perp - 3.6) / 0.9) ** 2) * np.exp(-(lateral / 2.2) ** 2)
            outer = np.array([nb[0] * nb[0], nb[1] * nb[1], nb[2] * nb[2],
                              nb[0] * nb[1], nb[0] * nb[2], nb[1] * nb[2]])
            stack[sl] += (wgt[..., None] * outer).astype(np.float32)
        self.maps = {'shape': chan, 'electrostatic': elec, 'desolvation': desolv, 'stacking': stack}

    def _sub(self, ax, centre, reach):
        """The slice of the grid within `reach` of a point, and the distances of those points to it."""
        i0 = np.maximum(np.floor((centre - reach - self.lo) / self.spacing).astype(int), 0)
        i1 = np.minimum(np.ceil((centre + reach - self.lo) / self.spacing).astype(int) + 1, self.shape)
        if np.any(i1 <= i0):
            return None, None
        sl = tuple(slice(int(i0[k]), int(i1[k])) for k in range(3))
        d2 = ((ax[0][sl[0]] - centre[0]) ** 2)[:, None, None] + ((ax[1][sl[1]] - centre[1]) ** 2)[None, :, None] \
            + ((ax[2][sl[2]] - centre[2]) ** 2)[None, None, :]
        return sl, np.sqrt(d2)

    # ------------------------------------------------------------------------------------- lookup
    def interpolate(self, field, xyz):
        """Trilinear interpolation of a map (or of a stack of maps) at the given points."""
        p = (np.asarray(xyz, float) - self.lo) / self.spacing
        inside = np.all((p >= 0) & (p <= np.array(self.shape) - 1.001), axis=1)
        p = np.clip(p, 0, np.array(self.shape) - 1.001)
        i = np.floor(p).astype(int)
        f = p - i
        out = 0.0
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = ((1 - f[:, 0]) if dx == 0 else f[:, 0]) * ((1 - f[:, 1]) if dy == 0 else f[:, 1]) \
                        * ((1 - f[:, 2]) if dz == 0 else f[:, 2])
                    v = field[i[:, 0] + dx, i[:, 1] + dy, i[:, 2] + dz]
                    out = out + (w[:, None] * v if v.ndim > 1 else w * v)
        return out, inside

    def prepare_ligand(self, lig, weights=None):
        """What a compound needs from the grid, worked out once per compound.

        Every atom falls into a class — its radius and the four flags the maps distinguish (apolar, donor,
        acceptor, able to reach a metal) — and each class has one map with the weights already folded in,
        so a pose costs a single interpolation of one stacked array rather than one per term. The stack is
        cached on the grid and shared by every compound with the same classes, which in a library is nearly
        all of them. Also three atoms of its own ring per aromatic atom, from which the ring normal of a
        pose is a cross product rather than a fit.
        """
        w = dict(S.WEIGHTS, **(weights or {}))
        radius = np.array([min(self.radii, key=lambda r: abs(r - x)) for x in lig.radius])
        metal_ok = lig.acceptor | (lig.charge < -0.3)
        keys = [(float(radius[a]), bool(lig.hydrophobic[a]), bool(lig.donor[a]), bool(lig.acceptor[a]),
                 bool(metal_ok[a])) for a in range(len(lig.radius))]
        classes = sorted(set(keys))
        stack = self._stack(tuple(classes), w)
        tri = np.zeros((len(lig.radius), 3), int)
        for ring in lig.rings:
            ring = sorted(set(ring))
            if len(ring) >= 3 and lig.aromatic[ring].all():
                tri[ring] = [ring[0], ring[1], ring[2]]
        return {'radius': radius, 'tri': tri, 'lig': lig, 'weights': w, 'stack': stack,
                'class_of': np.array([classes.index(k) for k in keys]), 'n_classes': len(classes),
                'classes': {r: np.nonzero(radius == r)[0] for r in sorted(set(radius))},
                'aromatic': np.nonzero(lig.aromatic)[0]}

    def _stack(self, classes, w):
        """The stacked, weighted maps for a set of atom classes: [class maps …, electrostatic, desolvation,
        the six stacking channels]. Cached, because a library shares its atom classes."""
        key = (classes, tuple(sorted((k, round(float(v), 6)) for k, v in w.items())))
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        ch = self.maps['shape']
        out = np.zeros(self.shape + (len(classes) + 8,), np.float32)
        for c, (r, phob, don, acc, met) in enumerate(classes):
            m = (w['gauss1'] * ch[(r, 'gauss1')] + w['gauss2'] * ch[(r, 'gauss2')]
                 + w['repulsion'] * ch[(r, 'repulsion')])
            if phob:
                m = m + w['hydrophobic'] * ch[(r, 'hydrophobic')]
            if don:
                m = m + w['hbond'] * ch[(r, 'hbond_donor')]
            if acc:
                m = m + w['hbond'] * ch[(r, 'hbond_acceptor')]
            if met:
                m = m + w['metal'] * ch[(r, 'metal')]
            out[..., c] = m
        out[..., len(classes)] = w['electrostatic'] * self.maps['electrostatic']
        out[..., len(classes) + 1] = w['desolvation'] * self.maps['desolvation']
        out[..., len(classes) + 2:] = w['stacking'] * self.maps['stacking']
        if len(self._cache) > 8:
            self._cache.clear()
        self._cache[key] = out
        return out

    def energy(self, xyz, lig, prep=None, weights=None, breakdown=False):
        """The score of a pose from the maps — the same quantity as `score.score`, by lookup.

        One interpolation of the stacked maps gives every term at once. The ring normals of the compound
        are taken from the pose itself, so a turned branch stacks as it really would.
        """
        if breakdown:
            return self._energy_terms(xyz, lig, prep, weights)
        prep = prep if prep is not None and (weights is None or prep['weights'] == dict(S.WEIGHTS, **weights)) \
            else self.prepare_ligand(lig, weights)
        xyz = np.asarray(xyz, float)
        vals, inside = self.interpolate(prep['stack'], xyz)
        nc = prep['n_classes']
        total = float(vals[np.arange(len(xyz)), prep['class_of']].sum())
        total += float(vals[:, nc] @ lig.charge + vals[:, nc + 1] @ np.abs(lig.charge))
        arom = prep['aromatic']
        if len(arom):
            m = vals[arom, nc + 2:]
            n = ring_normals_fast(xyz, prep['tri'])[arom]
            total += float(np.sum(m[:, 0] * n[:, 0] ** 2 + m[:, 1] * n[:, 1] ** 2 + m[:, 2] * n[:, 2] ** 2
                                  + 2 * (m[:, 3] * n[:, 0] * n[:, 1] + m[:, 4] * n[:, 0] * n[:, 2]
                                         + m[:, 5] * n[:, 1] * n[:, 2])))
        total /= 1.0 + prep['weights']['rotatable'] * lig.n_rotatable
        return total + OUTSIDE * float(np.sum(~inside))

    def energies(self, poses, lig, prep):
        """The score of many poses of one compound in a single lookup.

        `poses` is (m, atoms, 3). The search asks for a dozen poses at a time — the pose itself and one per
        parameter for the gradient — and doing them together turns a dozen trips through numpy into one,
        which is most of the speed of the minimiser.
        """
        poses = np.asarray(poses, float)
        m, n = poses.shape[0], poses.shape[1]
        vals, inside = self.interpolate(prep['stack'], poses.reshape(-1, 3))
        nc = prep['n_classes']
        vals = vals.reshape(m, n, -1)
        total = vals[:, np.arange(n), prep['class_of']].sum(1)
        total += vals[:, :, nc] @ lig.charge + vals[:, :, nc + 1] @ np.abs(lig.charge)
        arom = prep['aromatic']
        if len(arom):
            mm = vals[:, arom, nc + 2:]
            nn = np.stack([ring_normals_fast(x, prep['tri'])[arom] for x in poses])
            total += (mm[:, :, 0] * nn[:, :, 0] ** 2 + mm[:, :, 1] * nn[:, :, 1] ** 2
                      + mm[:, :, 2] * nn[:, :, 2] ** 2
                      + 2 * (mm[:, :, 3] * nn[:, :, 0] * nn[:, :, 1] + mm[:, :, 4] * nn[:, :, 0] * nn[:, :, 2]
                             + mm[:, :, 5] * nn[:, :, 1] * nn[:, :, 2])).sum(1)
        total /= 1.0 + prep['weights']['rotatable'] * lig.n_rotatable
        return total + OUTSIDE * (~inside.reshape(m, n)).sum(1)

    def _energy_terms(self, xyz, lig, prep=None, weights=None):
        """The same score, term by term, from the unweighted channels: for a report, and for the tests that
        check the maps against `score.pair_terms`."""
        w = dict(S.WEIGHTS, **(weights or {}))
        prep = prep or self.prepare_ligand(lig, weights)
        xyz = np.asarray(xyz, float)
        parts = {t: 0.0 for t in S.TERMS}
        outside = 0.0
        for r, idx in prep['classes'].items():
            if not len(idx):
                continue
            p = xyz[idx]
            metal_ok = lig.acceptor[idx] | (lig.charge[idx] < -0.3)
            for ch, term, mask in (('gauss1', 'gauss1', None), ('gauss2', 'gauss2', None),
                                   ('repulsion', 'repulsion', None),
                                   ('hydrophobic', 'hydrophobic', lig.hydrophobic[idx]),
                                   ('hbond_donor', 'hbond', lig.donor[idx]),
                                   ('hbond_acceptor', 'hbond', lig.acceptor[idx]),
                                   ('metal', 'metal', metal_ok)):
                if mask is not None and not mask.any():
                    continue
                v, inside = self.interpolate(self.maps['shape'][(r, ch)], p)
                parts[term] += float(np.sum(v if mask is None else v[mask]))
            outside += float(np.sum(~self.interpolate(self.maps['shape'][(r, 'gauss1')], p)[1]))
        v, _in = self.interpolate(self.maps['electrostatic'], xyz)
        parts['electrostatic'] = float(np.sum(v * lig.charge))
        v, _in = self.interpolate(self.maps['desolvation'], xyz)
        parts['desolvation'] = float(np.sum(v * np.abs(lig.charge)))
        arom = prep['aromatic']
        if len(arom):
            m, _in = self.interpolate(self.maps['stacking'], xyz[arom])
            n = ring_normals_fast(xyz, prep['tri'])[arom]
            parts['stacking'] = float(np.sum(m[:, 0] * n[:, 0] ** 2 + m[:, 1] * n[:, 1] ** 2
                                             + m[:, 2] * n[:, 2] ** 2
                                             + 2 * (m[:, 3] * n[:, 0] * n[:, 1] + m[:, 4] * n[:, 0] * n[:, 2]
                                                    + m[:, 5] * n[:, 1] * n[:, 2])))
        total = sum(w[t] * parts[t] for t in S.TERMS) / (1.0 + w['rotatable'] * lig.n_rotatable)
        total += OUTSIDE * outside
        return total, {'raw': parts, 'weighted': {t: w[t] * parts[t] for t in S.TERMS}, 'outside': outside}

    @classmethod
    def for_pocket(cls, rec, pocket=None, centre=None, size=None, pad=PAD, **kw):
        """A grid over a pocket from `pockets.find` (its points), or over a box given by centre and size."""
        if pocket is not None:
            pts = np.asarray(pocket['points'] if isinstance(pocket, dict) else pocket, float)
            if not len(pts):
                raise ValueError('the pocket has no points: call pockets.find(..., keep_points=True)')
            lo, hi = pts.min(0) - pad, pts.max(0) + pad
        else:
            c = np.asarray(centre, float)
            s = np.asarray(size if size is not None else 20.0, float) * np.ones(3) / 2
            lo, hi = c - s, c + s
        return cls(rec, lo, hi, **kw)

    @property
    def centre(self):
        return (self.lo + self.hi) / 2

    @property
    def size(self):
        return self.hi - self.lo

    def __repr__(self):
        return (f'<Grid {self.shape[0]}x{self.shape[1]}x{self.shape[2]} at {self.spacing} Å, '
                f'{len(self.radii)} radii, {self.nbytes / 1e6:.0f} MB>')

    @property
    def nbytes(self):
        return sum(v.nbytes for v in self.maps['shape'].values()) + self.maps['electrostatic'].nbytes \
            + self.maps['desolvation'].nbytes + self.maps['stacking'].nbytes


def ring_normals_fast(xyz, tri):
    """Ring normals from three atoms of each ring: the pose's own normals, at the cost of a cross product.

    A ring is rigid, but the branch it sits on may have turned, so the normal has to come from the moved
    atoms. `tri` is what `Grid.prepare_ligand` worked out once (three ring atoms per aromatic atom)."""
    xyz = np.asarray(xyz, float)
    a = xyz[tri[:, 1]] - xyz[tri[:, 0]]
    b = xyz[tri[:, 2]] - xyz[tri[:, 0]]
    n = np.cross(a, b)
    norm = np.linalg.norm(n, axis=1)
    return n / np.where(norm > 1e-9, norm, 1.0)[:, None]
