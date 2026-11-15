"""Where could a small molecule bind? Pockets of an RNA structure, and how they differ between two conformations.

The method is the classic buriedness scan (LIGSITE / fpocket in spirit), applied to RNA:

1. The structure is put on a grid (0.8 Å). A grid point is *occupied* if it is within the van der Waals
   radius of an atom plus a probe radius (1.4 Å, a water): a compound atom cannot sit there.
2. Every free grid point is scanned along seven directions (the three axes and four cube diagonals). A
   direction counts as *enclosing* when the structure is met on both sides of the point within 12 Å. The
   count, 0 to 7, is the buriedness of the point: 0 in the bulk solvent, 7 inside a closed cavity.
3. Free points with a buriedness of at least 5 of 7 are kept and grouped into connected clusters. A
   cluster is a pocket; its volume is the number of points times the volume of a grid cell.

   The threshold decides what counts as a site, and it was chosen by measuring rather than by looking at
   pictures: `benchmark/screening/pocket_sites.py` takes thirty experimental RNA-ligand complexes, removes
   the compound and asks where the site it occupied appears in the ranking. At 5 of 7 the site is the
   top-ranked pocket in 21 of the 30 and is found at all in 24; at 6 of 7, in 9 and 12. Six gives tidier,
   better separated sites — and misses the shallow ones, which is what most RNA sites are. Hence 5, with
   `burial_min=6` for the deep pockets only.
4. Each pocket is described by the nucleotides lining it (any atom within 4.5 Å of one of its points),
   its mean buriedness, and how much of its lining is not in a Watson-Crick helix. RNA-binding small
   molecules bind loops, bulges and junctions, not the grooves of an A-form helix, so that last term is
   what separates a site worth screening from a stretch of major groove.

`compare` matches the pockets of two conformations by the nucleotides that line them, not by position in
space, so it does not depend on how the two structures were superposed. A pocket that is open in one
conformation and closed in the other is what a conformation-selective screen is looking for.

Nothing here is an affinity: it is geometry. A pocket in a predicted model is only as good as the model.

This module has no dependency beyond numpy and scipy, and takes plain coordinates, so it can be used on
any structure:

    from rnafold3d.screening import pockets
    found = pockets.find(xyz, elements=elements, residue_of=residue_of, flexible=flexible)
"""
from __future__ import annotations

import numpy as np

VDW = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'P': 1.80, 'S': 1.80, 'F': 1.47, 'CL': 1.75,
       'BR': 1.85, 'I': 1.98, 'MG': 1.73, 'NA': 2.27, 'K': 2.75, 'ZN': 1.39, 'MN': 1.61, 'H': 1.20}
DEFAULT_VDW = 1.70                 # an unknown element is treated as carbon
PROBE = 1.4                        # Å, water; the radius of the smallest useful ligand atom
SPACING = 0.8                      # Å, grid step
SCAN = 12.0                        # Å, how far a buriedness scan looks for the structure in one direction
BURIAL_MIN = 5                     # of the 7 direction pairs, for a grid point to belong to a pocket
CONTACT = 4.5                      # Å, atom to pocket point, for a nucleotide to count as lining the pocket
MIN_VOLUME = 100.0                 # Å³; smaller clusters hold no drug-like compound
WIDE_VOLUME = 1500.0               # Å³; above this a site is a groove or a cleft rather than a pocket
MAX_POINTS = 8_000_000             # grid points; above this the spacing is coarsened (a very long RNA)

# the three axes and four cube diagonals: 7 direction pairs, as in LIGSITE
DIRECTIONS = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1))


def radii(elements):
    """Van der Waals radii of a list of element symbols."""
    return np.array([VDW.get(str(e).strip().upper(), DEFAULT_VDW) for e in elements], float)


def _shift(a, d):
    """`a` moved by one lattice step `d`, the edge padded with False."""
    out = np.zeros_like(a)
    src = tuple(slice(max(0, -k), n - max(0, k)) for k, n in zip(d, a.shape))
    dst = tuple(slice(max(0, k), n - max(0, -k)) for k, n in zip(d, a.shape))
    out[dst] = a[src]
    return out


def _occupancy(xyz, rad, spacing, margin):
    """Boolean grid of the points a probe centre cannot reach, with the grid origin."""
    lo = xyz.min(0) - margin
    hi = xyz.max(0) + margin
    shape = np.maximum(np.ceil((hi - lo) / spacing).astype(int) + 1, 3)
    occ = np.zeros(tuple(shape), bool)
    cell = np.arange(3)
    for p, r in zip(xyz, rad):
        r = r + PROBE
        i0 = np.maximum(np.floor((p - r - lo) / spacing).astype(int), 0)
        i1 = np.minimum(np.ceil((p + r - lo) / spacing).astype(int) + 1, shape)
        if np.any(i1 <= i0):
            continue
        axes = [lo[k] + np.arange(i0[k], i1[k]) * spacing - p[k] for k in cell]
        d2 = (axes[0] ** 2)[:, None, None] + (axes[1] ** 2)[None, :, None] + (axes[2] ** 2)[None, None, :]
        occ[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] |= d2 <= r * r
    return occ, lo


def _buriedness(occ, spacing, scan=SCAN):
    """For every grid point, in how many of the 7 direction pairs the structure lies on both sides (0-7)."""
    out = np.zeros(occ.shape, np.uint8)
    for d in DIRECTIONS:
        steps = max(1, int(round(scan / (spacing * float(np.linalg.norm(d))))))
        both = None
        for sense in (1, -1):
            step = tuple(sense * k for k in d)
            seen = np.zeros_like(occ)
            cur = occ
            for _ in range(steps):
                cur = _shift(cur, step)
                seen |= cur
            both = seen if both is None else (both & seen)
        out += both
    return out


def lining_residues(points, xyz, residue_of, cutoff=CONTACT):
    """Nucleotides with an atom within `cutoff` of one of the pocket points, and how many atoms each has there."""
    from scipy.spatial import cKDTree
    if not len(points):
        return {}
    near = cKDTree(points).query_ball_point(xyz, cutoff, return_length=True)
    out = {}
    for a in np.nonzero(near)[0]:
        r = int(residue_of[a])
        out[r] = out.get(r, 0) + 1
    return out


def find(xyz, elements=None, radius=None, residue_of=None, flexible=None, spacing=SPACING,
         burial_min=BURIAL_MIN, min_volume=MIN_VOLUME, contact=CONTACT, keep_points=True, top=None):
    """Pockets of one structure, most promising first.

    xyz          (n_atoms, 3) coordinates of the heavy atoms
    elements     element symbol per atom (or `radius`, the van der Waals radius per atom)
    residue_of   residue index per atom; without it the lining is reported per atom index
    flexible     per residue, how much the nucleotide is outside a Watson-Crick helix: 1 for unpaired or
                 non-canonically paired, 0 inside a helix (used for the score, and nothing else)
    spacing      grid step in Å (coarsened automatically for a very long RNA)
    burial_min   how many of the 7 direction pairs must enclose a point (5 of 7; 6 for deep pockets only)
    min_volume   smallest pocket to report, in Å³
    keep_points  include the grid points of each pocket ('points'), which a docking grid is built on
    top          report at most this many pockets

    Each pocket is a dict with: volume (Å³), center, burial (0-1, mean over its points), enclosure
    (fraction of its points enclosed in all seven direction pairs), extent (diagonal of its bounding
    box, Å), residues (sorted indices of the nucleotides lining it), contacts (atoms near the
    pocket, per residue), loop_fraction (the contact-weighted share of the lining outside a helix, or
    None without `flexible`), score (see `_score`), wide (a groove rather than a pocket), and the grid
    spacing and points.
    """
    from scipy import ndimage
    xyz = np.asarray(xyz, float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or len(xyz) < 4:
        raise ValueError('pocket finding needs at least four atoms as (n, 3) coordinates')
    rad = np.asarray(radius, float) if radius is not None else (
        radii(elements) if elements is not None else np.full(len(xyz), DEFAULT_VDW))
    if residue_of is None:
        residue_of = np.arange(len(xyz))
    residue_of = np.asarray(residue_of, int)

    reach = float(rad.max() + PROBE)
    box = (xyz.max(0) - xyz.min(0) + 2 * (reach + 2 * spacing))
    spacing = float(max(spacing, (box.prod() / MAX_POINTS) ** (1 / 3)))     # a very long RNA gets a coarser grid
    margin = reach + 2 * spacing
    occ, origin = _occupancy(xyz, rad, spacing, margin)
    burial = _buriedness(occ, spacing)
    mask = (~occ) & (burial >= burial_min)
    labels, _n = ndimage.label(mask)                     # 6-connectivity: a diagonal leak does not merge pockets
    cell = spacing ** 3
    out = []
    for lab, window in enumerate(ndimage.find_objects(labels), 1):
        if window is None:
            continue
        sub = labels[window] == lab
        volume = int(sub.sum()) * cell
        if volume < min_volume:
            continue
        idx = np.argwhere(sub) + [w.start for w in window]
        pts = origin + idx * spacing
        b = burial[window][sub].astype(float)
        contacts = lining_residues(pts, xyz, residue_of, contact)
        res = sorted(contacts)
        loop = None
        if flexible is not None and contacts:
            f = np.asarray(flexible, float)
            w = np.array([contacts[r] for r in res], float)
            loop = float(np.clip(np.sum(w * f[res]) / w.sum(), 0.0, 1.0))
        p = {'volume': round(volume, 1), 'center': [round(float(v), 2) for v in pts.mean(0)],
             'burial': round(float(b.mean()) / len(DIRECTIONS), 3),
             'enclosure': round(float(np.mean(b >= len(DIRECTIONS))), 3),
             'extent': round(float(np.linalg.norm(pts.max(0) - pts.min(0))), 1),
             'min': [round(float(v), 2) for v in pts.min(0)],
             'max': [round(float(v), 2) for v in pts.max(0)],
             'residues': res, 'contacts': {int(r): int(c) for r, c in sorted(contacts.items())},
             'loop_fraction': None if loop is None else round(loop, 3),
             'wide': bool(volume > WIDE_VOLUME), 'spacing': round(spacing, 3), 'n_points': int(len(idx))}
        p['score'] = _score(p)
        if keep_points:
            p['points'] = pts
        out.append(p)
    out.sort(key=lambda p: (-p['score'], -p['volume']))
    if top:
        out = out[:int(top)]
    for i, p in enumerate(out):
        p['id'] = i + 1
    return out


def _score(p):
    """How promising a pocket is for a screen, 0 to 1: enclosed, large enough, and not a helix groove.

    score = burial x (0.4 + 0.6 x loop_fraction) x min(1, volume / 300 Å³)

    The three terms are the three ways a site fails: open to the solvent, too small for a compound, or a
    stretch of A-form groove, which is where docking programs put their false positives. Without pairing
    information the middle term is 1 (the pocket is judged on geometry alone). This is a ranking, not a
    measurement: two pockets within ~0.05 are not meaningfully different.
    """
    loop = p['loop_fraction']
    shape = 1.0 if loop is None else 0.4 + 0.6 * loop
    return round(float(p['burial'] * shape * min(1.0, p['volume'] / 300.0)), 3)


def volume_near(pockets, residues, xyz, residue_of, cutoff=CONTACT):
    """How much pocket volume of one structure lies against a given set of nucleotides.

    For every pocket, the share of its grid points within `cutoff` of an atom of those nucleotides, times
    its volume. This is what makes a site comparable between two conformations: the same nucleotides are
    asked how much open space they hold, so a pocket that splits in two, or merges with its neighbour, is
    still measured, which matching whole pockets to each other cannot do.

    Returns the total volume (Å³) and the pockets that contribute: [(index, volume, share)].
    """
    from scipy.spatial import cKDTree
    sel = np.isin(np.asarray(residue_of, int), np.asarray(sorted(residues), int))
    if not sel.any() or not pockets:
        return 0.0, []
    tree = cKDTree(np.asarray(xyz, float)[sel])
    total, parts = 0.0, []
    for i, p in enumerate(pockets):
        pts = p.get('points')
        if pts is None or not len(pts):
            raise ValueError('volume_near needs the pocket points: call find(..., keep_points=True)')
        share = float(np.mean(tree.query_ball_point(pts, cutoff, return_length=True) > 0))
        if share > 0.05:
            total += share * p['volume']
            parts.append((i, round(share * p['volume'], 1), round(share, 3)))
    parts.sort(key=lambda t: -t[1])
    return round(total, 1), parts


def compare(pockets_a, pockets_b, xyz_a, xyz_b, residue_of, cutoff=CONTACT, same_site=0.3):
    """The same sites in two conformations of one RNA: which pocket opens where.

    A site is the set of nucleotides lining a pocket. Each site of A is measured in both conformations with
    `volume_near` — the open volume held by those nucleotides — and the same is done for the sites of B that
    are not already covered by a site of A (Jaccard overlap below `same_site`). The sites are compared by
    their nucleotides, never by position in space, so the superposition of the two structures does not
    matter, and pockets that merge or split between the conformations are still compared.

    One entry per site, the largest change first:
      side            'A' or 'B': the conformation whose pocket defined the site
      a, b            index of the pocket that defined it, and the largest contributing pocket in the other
      residues        the nucleotides of the site
      volume_a/_b     open volume at those nucleotides in each conformation, Å³
      difference      (volume_a - volume_b) / max(volume_a, volume_b): +1 open in A only, -1 in B only
      overlap         Jaccard overlap of the two pockets' lining nucleotides (0 when the site is absent)
      open_in         'A', 'B' or 'both' (|difference| <= 0.25 counts as both)
    """
    sets_a = [set(p['residues']) for p in pockets_a]
    sets_b = [set(p['residues']) for p in pockets_b]
    sites = [('A', i, sets_a[i]) for i in range(len(pockets_a))]
    for j, sb in enumerate(sets_b):
        if all(not (sa & sb) or len(sa & sb) / len(sa | sb) < same_site for sa in sets_a):
            sites.append(('B', j, sb))
    rows = []
    for side, k, res in sites:
        va, pa = volume_near(pockets_a, res, xyz_a, residue_of, cutoff)
        vb, pb = volume_near(pockets_b, res, xyz_b, residue_of, cutoff)
        if max(va, vb) <= 0:
            continue
        diff = (va - vb) / max(va, vb)
        i = k if side == 'A' else (pa[0][0] if pa else None)
        j = k if side == 'B' else (pb[0][0] if pb else None)
        ov = 0.0
        if i is not None and j is not None:
            union = len(sets_a[i] | sets_b[j])
            ov = len(sets_a[i] & sets_b[j]) / union if union else 0.0
        rows.append({'side': side, 'a': i, 'b': j, 'residues': sorted(res),
                     'volume_a': va, 'volume_b': vb, 'overlap': round(float(ov), 3),
                     'difference': round(float(diff), 3),
                     'open_in': 'both' if abs(diff) <= 0.25 else ('A' if diff > 0 else 'B')})
    rows.sort(key=lambda r: (-abs(r['difference']), -max(r['volume_a'], r['volume_b'])))
    return rows
