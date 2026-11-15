"""More experimental data for an ensemble: SAXS profiles, FRET distances and cryo-EM maps.

Each data type follows the path of the chemical probing data (probing.py): read the file, compute what every
conformation predicts, score each conformation's agreement with the data, and reweight the ensemble by maximum
entropy so that its average fits the data (Hummer & Köfinger 2015; Bottaro et al. 2020). The reweighting is
generic (maxent()): every dataset contributes observables per conformation, data values, errors and, where the
data have an arbitrary scale, a fitted linear scale. Several datasets (SAXS + FRET + SHAPE, ...) can be fitted
together; theta sets how far the data may move the weights away from the prior.

SAXS   I(q) per conformation by the Debye formula over all heavy atoms, from a distance histogram per element
       pair (0.1 A bins) with Cromer-Mann atomic form factors minus a Gaussian excluded-volume term (Fraser et al.
       1978; solvent 0.334 e/A^3). No hydration shell: the fitted scale c and constant background b of
       I_exp ~ c I_calc + b absorb part of it, and the agreement is best judged at q < 0.25 A^-1. Agreement:
       reduced chi2 with the experimental errors (or 3 % of I where none are given).
FRET   distances between two nucleotides (residue numbers; the C4' atom unless another is named), given as a
       distance in A or as an efficiency with the Foerster radius R0. Without dye models the dye linkers are not
       represented: the default error is 5 A. Observables: the efficiency E = 1/(1 + (r/R0)^6) when R0 is given
       (the ensemble average of E is what is measured), the distance otherwise.
cryo-EM a density map (MRC/CCP4). Every conformation is placed in the map (principal axes of the density above
       the contour level against those of the model, 4 proper orientations, then a local search of rotation and
       shift), its density is simulated at the map resolution (Gaussians on the atoms) and compared by the
       correlation over the voxels above the contour; the ensemble average density is fitted for reweighting.
"""
from __future__ import annotations

import base64
import hashlib
import math
import re
import struct

import numpy as np

# Cromer-Mann coefficients (International Tables Vol. C): a1..a4, b1..b4, c
CROMER_MANN = {
    'C': (2.31000, 1.02000, 1.58860, 0.865000, 20.8439, 10.2075, 0.568700, 51.6512, 0.215600),
    'N': (12.2126, 3.13220, 2.01250, 1.16630, 0.005700, 9.89330, 28.9975, 0.582600, -11.529),
    'O': (3.04850, 2.28680, 1.54630, 0.867000, 13.2771, 5.70110, 0.323900, 32.9089, 0.250800),
    'P': (6.43450, 4.17910, 1.78000, 1.49080, 1.90670, 27.1570, 0.526000, 68.1645, 1.11490),
}
VOLUME = {'C': 16.44, 'N': 2.49, 'O': 9.13, 'P': 5.73}    # displaced solvent volumes (A^3), Svergun/CRYSOL
RHO_SOLVENT = 0.334
ELEMENTS = ('C', 'N', 'O', 'P')


# ------------------------------------------------------------------------------------------- reweighting
def maxent(datasets, K, theta=None, prior=None):
    """Maximum-entropy weights over K conformations for one or more datasets.

    datasets: [{'O': K x m predicted observables, 'y': m data, 'sigma': m errors,
                'scale': 'none' | 'linear' (y ~ a O + b, a >= 0) | 'multiplicative' (y ~ a O)}].
    Minimises sum chi2/2 - theta * S(w || w0). Returns weights, chi2/m per dataset before and after, the
    effective number of conformations and the fitted scales."""
    from scipy.optimize import minimize
    w0 = np.full(K, 1.0 / K) if prior is None else np.clip(np.asarray(prior, float), 1e-12, None)
    w0 = w0 / w0.sum()
    m_tot = sum(len(d['y']) for d in datasets)
    theta = float(theta) if theta else max(1.0, 0.1 * m_tot)
    D = []
    for d in datasets:
        O = np.asarray(d['O'], float)
        y = np.asarray(d['y'], float)
        s = np.asarray(d['sigma'], float)
        D.append((O, y, np.clip(s, 1e-9, None), d.get('scale', 'none')))

    def scale(mu, y, s, kind):
        if kind == 'none':
            return 1.0, 0.0
        wts = 1.0 / s ** 2
        if kind == 'multiplicative':
            a = float(np.sum(wts * mu * y) / max(1e-30, np.sum(wts * mu * mu)))
            return max(0.0, a), 0.0
        A = np.vstack([mu, np.ones_like(mu)]).T * np.sqrt(wts)[:, None]
        (a, b), *_ = np.linalg.lstsq(A, y * np.sqrt(wts), rcond=None)
        if a < 0:
            a, b = 0.0, float(np.sum(wts * y) / np.sum(wts))
        return float(a), float(b)

    def terms(w):
        out = []
        for O, y, s, kind in D:
            mu = w @ O
            a, b = scale(mu, y, s, kind)
            r = (y - a * mu - b) / s
            out.append((float(r @ r), a, b, r, s, O))
        return out

    def f(z):
        z = z - z.max()
        w = w0 * np.exp(z)
        w /= w.sum()
        tt = terms(w)
        S = -float(np.sum(w * np.log(w / w0)))
        F = 0.5 * sum(t[0] for t in tt) - theta * S
        g_w = theta * (np.log(w / w0) + 1)
        for c2, a, b, r, s, O in tt:
            g_w = g_w - (r / s * a) @ O.T
        return F, w * (g_w - float(w @ g_w))
    res = minimize(f, np.zeros(K), jac=True, method='L-BFGS-B', options={'maxiter': 800})
    z = res.x - res.x.max()
    w = w0 * np.exp(z)
    w /= w.sum()
    before, after = terms(w0), terms(w)
    S = -float(np.sum(w * np.log(w / w0)))
    return {'weights': [round(float(x), 6) for x in w], 'prior': [round(float(x), 6) for x in w0], 'theta': theta,
            'chi2_prior': [round(t[0] / len(d[1]), 4) for t, d in zip(before, D)],
            'chi2': [round(t[0] / len(d[1]), 4) for t, d in zip(after, D)],
            'scales': [[round(t[1], 6), round(t[2], 6)] for t in after], 'n_data': [len(d[1]) for d in D],
            'phi_eff': round(math.exp(S), 4), 'n_eff_conformations': round(math.exp(-float(np.sum(w * np.log(w)))), 2),
            'converged': bool(res.success)}


# ------------------------------------------------------------------------------------------- shared helpers
def _numbers(line):
    out = []
    for t in re.split(r'[\s,;]+', line.strip()):
        try:
            out.append(float(t))
        except ValueError:
            return None
    return out


def _heavy(doc, coords):
    el = np.array([str(e).upper()[:1] for e in doc['atoms']['elements']])
    keep = np.isin(el, ELEMENTS)
    return el[keep], np.asarray(coords, float)[:, keep, :]


# ------------------------------------------------------------------------------------------- SAXS
def parse_saxs(text):
    """q, I, sigma from a 2-4 column profile (.dat, ATSAS/BioSAXS); q in 1/nm is converted to 1/A."""
    rows = []
    for line in (text or '').splitlines():
        if not line.strip() or line.lstrip()[:1] in '#;!%':
            continue
        v = _numbers(line)
        if v and len(v) >= 2 and v[1] > 0 and v[0] > 0:
            rows.append(v[:3])
    if len(rows) < 10:
        raise ValueError('a SAXS profile needs at least 10 lines of q and I(q)')
    q = np.array([r[0] for r in rows])
    I = np.array([r[1] for r in rows])
    s = np.array([r[2] if len(r) > 2 and r[2] > 0 else 0.03 * r[1] for r in rows])
    notes = []
    unit = '1/A'
    if q.max() > 1.0:                       # 1/nm (q up to ~5 /nm)
        q = q / 10.0
        unit = '1/nm (converted)'
    keep = q <= 0.5
    if not keep.all():
        notes.append(f'{int((~keep).sum())} points above q = 0.5 1/A left out')
    q, I, s = q[keep], I[keep], s[keep]
    if len(rows[0]) < 3:
        notes.append('no errors in the file: 3 % of I(q) assumed')
    return {'q': q.round(6).tolist(), 'I': I.tolist(), 'sigma': s.tolist(), 'unit': unit, 'notes': notes,
            'q_range': [round(float(q.min()), 4), round(float(q.max()), 4)], 'n_points': int(len(q))}


def form_factor(e, q):
    a1, a2, a3, a4, b1, b2, b3, b4, c = CROMER_MANN[e]
    s2 = (q / (4 * math.pi)) ** 2
    f = a1 * np.exp(-b1 * s2) + a2 * np.exp(-b2 * s2) + a3 * np.exp(-b3 * s2) + a4 * np.exp(-b4 * s2) + c
    V = VOLUME[e]
    g = RHO_SOLVENT * V * np.exp(-q * q * V ** (2 / 3) / (4 * math.pi))
    return f - g


def saxs_profile(elements, X, q, bin_width=0.1):
    """Debye I(q) of one structure (heavy atoms) via element-pair distance histograms."""
    from scipy.spatial.distance import pdist
    q = np.asarray(q, float)
    F = {e: form_factor(e, q) for e in ELEMENTS}
    I = np.zeros(len(q))
    idx = {e: np.where(elements == e)[0] for e in ELEMENTS}
    for e in ELEMENTS:
        I += len(idx[e]) * F[e] ** 2
    dmax = float(np.max(pdist(X[::max(1, len(X) // 400)]))) + 20.0 if len(X) > 1 else 1.0
    edges = np.arange(0.0, dmax + bin_width, bin_width)
    r = 0.5 * (edges[1:] + edges[:-1])
    qr = np.outer(q, r)
    sinc = np.where(qr > 1e-9, np.sin(qr) / np.where(qr > 1e-9, qr, 1.0), 1.0)
    for a in range(len(ELEMENTS)):
        for b in range(a, len(ELEMENTS)):
            ea, eb = ELEMENTS[a], ELEMENTS[b]
            A, B = X[idx[ea]], X[idx[eb]]
            if not len(A) or not len(B):
                continue
            # every unordered atom pair once (pdist within an element, all cross pairs between two), counted twice
            d = pdist(A) if a == b else np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2).ravel()
            h, _ = np.histogram(d, bins=edges)
            I += 2.0 * F[ea] * F[eb] * (sinc @ h)
    return I


def saxs_observables(doc, coords, data):
    el, X = _heavy(doc, coords)
    q = np.array(data['q'])
    return np.array([saxs_profile(el, X[k], q) for k in range(len(X))])


# ------------------------------------------------------------------------------------------- FRET
def parse_fret(text, n):
    """FRET pairs, one per line: `i j distance [sd]` or `i j E=0.45 R0=54 [sd=...]`; options atom=C5' (both)
    or atom1=/atom2=. Residue numbers 1..n."""
    out, notes = [], []
    for ln, line in enumerate((text or '').splitlines(), 1):
        line = line.split('#')[0].strip()
        if not line:
            continue
        toks = re.split(r'[\s,;]+', line)
        kv = {t.split('=', 1)[0].lower(): t.split('=', 1)[1] for t in toks if '=' in t}
        plain = [t for t in toks if '=' not in t]
        try:
            i, j = int(plain[0]), int(plain[1])
        except (IndexError, ValueError):
            if ln == 1:
                continue                          # a header line
            raise ValueError(f'line {ln}: two residue numbers expected') from None
        if not (1 <= i <= n and 1 <= j <= n) or i == j:
            raise ValueError(f'line {ln}: residues must be two different numbers in 1..{n}')
        item = {'i': i, 'j': j, 'atom1': kv.get('atom1', kv.get('atom', "C4'")), 'atom2': kv.get('atom2', kv.get('atom', "C4'"))}
        if 'e' in kv:
            E = float(kv['e'])
            R0 = float(kv.get('r0', 54.0))
            if not 0 < E < 1:
                raise ValueError(f'line {ln}: an efficiency is between 0 and 1')
            item.update({'E': E, 'R0': R0, 'distance': round(R0 * (1 / E - 1) ** (1 / 6), 2),
                         'sd': float(kv.get('sd', 0.05))})
            if 'r0' not in kv:
                notes.append(f'line {ln}: R0 = 54 A assumed (Cy3/Cy5)')
        else:
            if len(plain) < 3:
                raise ValueError(f'line {ln}: a distance (A) or E=... is needed')
            item.update({'distance': float(plain[2]), 'sd': float(plain[3]) if len(plain) > 3 else float(kv.get('sd', 5.0))})
        out.append(item)
    if not out:
        raise ValueError('no FRET pairs found')
    return {'pairs': out, 'notes': sorted(set(notes)), 'n_points': len(out)}


def _atom_xyz(doc, coords, res, name):
    r = doc['residues'][res]
    names = doc['atoms']['names']
    for q in range(r['n_atoms']):
        if names[r['first_atom'] + q] == name:
            return np.asarray(coords, float)[:, r['first_atom'] + q, :]
    raise ValueError(f'residue {res + 1} has no atom {name}')


def fret_observables(doc, coords, data):
    """K x m: efficiency for pairs given as E (with R0), distance (A) otherwise."""
    cols = []
    for p in data['pairs']:
        d = np.linalg.norm(_atom_xyz(doc, coords, p['i'] - 1, p['atom1']) - _atom_xyz(doc, coords, p['j'] - 1, p['atom2']), axis=1)
        cols.append(1.0 / (1.0 + (d / p['R0']) ** 6) if 'E' in p else d)
    return np.array(cols).T


def fret_targets(data):
    y = np.array([p['E'] if 'E' in p else p['distance'] for p in data['pairs']], float)
    s = np.array([p['sd'] for p in data['pairs']], float)
    return y, s


# ------------------------------------------------------------------------------------------- cryo-EM
_MAPS = {}


def parse_map(raw, name='map', resolution=None, contour=None):
    """Read an MRC/CCP4 map (bytes or base64); keeps it in memory for fitting. Returns a summary with its id."""
    if isinstance(raw, str):
        raw = base64.b64decode(raw.split(',', 1)[-1])
    if len(raw) < 1024:
        raise ValueError('not an MRC/CCP4 map (file too short)')
    hdr = raw[:1024]
    endian = '<' if struct.unpack('<i', hdr[12:16])[0] in (0, 1, 2, 6, 12) else '>'
    nx, ny, nz, mode = struct.unpack(endian + '4i', hdr[:16])
    nxs, nys, nzs = struct.unpack(endian + '3i', hdr[16:28])
    mx, my, mz = struct.unpack(endian + '3i', hdr[28:40])
    cella = struct.unpack(endian + '3f', hdr[40:52])
    mapc, mapr, maps = struct.unpack(endian + '3i', hdr[64:76])
    nsym = struct.unpack(endian + 'i', hdr[92:96])[0]
    origin = struct.unpack(endian + '3f', hdr[196:208])
    dtype = {0: 'i1', 1: 'i2', 2: 'f4', 6: 'u2', 12: 'f2'}.get(mode)
    if dtype is None or min(nx, ny, nz) < 2:
        raise ValueError(f'unsupported MRC map (mode {mode})')
    off = 1024 + nsym
    count = nx * ny * nz
    data = np.frombuffer(raw[off:off + count * np.dtype(dtype).itemsize], dtype=endian + dtype).astype(np.float32)
    if len(data) != count:
        raise ValueError('the map is truncated')
    vol = data.reshape(nz, ny, nx)                     # sections, rows, columns
    order = [mapc - 1, mapr - 1, maps - 1]             # axis (x=0, y=1, z=2) of columns, rows, sections
    if sorted(order) != [0, 1, 2]:
        order = [0, 1, 2]
    vol = np.transpose(vol, (2, 1, 0))                 # columns, rows, sections
    inv = np.argsort(order)
    vol = np.transpose(vol, inv)                       # x, y, z
    starts = np.array([nxs, nys, nzs], float)[inv]
    grid = np.array([mx, my, mz], float)
    voxel = np.array(cella) / np.where(grid > 0, grid, 1)
    org = np.array(origin, float)
    if not np.any(org):
        org = starts * voxel
    if contour is None:
        contour = float(np.mean(vol) + 2.0 * np.std(vol))
    if resolution is None:
        resolution = round(float(max(4.0, 3.0 * voxel.mean())), 1)
    key = hashlib.sha1(raw[:4096] + str(len(raw)).encode()).hexdigest()[:16]
    _MAPS[key] = {'vol': vol, 'voxel': voxel, 'origin': org}
    while len(_MAPS) > 4:
        _MAPS.pop(next(iter(_MAPS)))
    return {'map': key, 'name': name, 'shape': list(vol.shape), 'voxel_A': [round(float(v), 3) for v in voxel],
            'resolution': float(resolution), 'contour': round(float(contour), 6),
            'n_points': int(np.sum(vol > contour)), 'notes': []}


def _map_points(m, contour):
    vol = m['vol']
    idx = np.argwhere(vol > contour)
    return idx * m['voxel'] + m['origin'], vol[vol > contour]


def _simulate(X, pts, resolution, tree=None):
    """Simulated density of atoms X at the points pts (Gaussian per atom, sigma = resolution / (pi sqrt 2))."""
    from scipy.spatial import cKDTree
    sig = resolution / (math.pi * math.sqrt(2.0))
    tp = tree if tree is not None else cKDTree(pts)
    sm = cKDTree(X).sparse_distance_matrix(tp, 3.0 * sig, output_type='coo_matrix')
    return np.bincount(sm.col, weights=np.exp(-sm.data ** 2 / (2 * sig * sig)), minlength=len(pts))


def _corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b) / den if den > 0 else 0.0


def _rot(v):
    from scipy.spatial.transform import Rotation
    return Rotation.from_rotvec(v).as_matrix()


def fit_in_map(X, pts, dens, resolution, sub=1500):
    """Rigid placement of model X in the map: principal axes, then a local search. Returns (X placed, cc)."""
    rng = np.random.default_rng(0)
    from scipy.spatial import cKDTree
    sel = rng.choice(len(pts), sub, replace=False) if len(pts) > sub else np.arange(len(pts))
    P, D = pts[sel], dens[sel]
    TP = cKDTree(P)
    cm, cx = pts.mean(0), X.mean(0)
    _, vm = np.linalg.eigh(np.cov((pts - cm).T))
    _, vx = np.linalg.eigh(np.cov((X - cx).T))
    best = (-2.0, None)
    for sx in (1, -1):
        for sy in (1, -1):
            R = vx @ np.diag([sx, sy, sx * sy]) @ vm.T          # proper rotation: model axes onto map axes
            if np.linalg.det(R) < 0:
                R = vx @ np.diag([sx, sy, -sx * sy]) @ vm.T
            Y = (X - cx) @ R + cm
            cc = _corr(_simulate(Y, P, resolution, TP), D)
            if cc > best[0]:
                best = (cc, (R, cm))
    cc, (R, t) = best
    step_a, step_t = 0.15, 2.0
    for _it in range(60):
        improved = False
        for k in range(6):
            for sgn in (1, -1):
                dv = np.zeros(6)
                dv[k] = sgn * (step_a if k < 3 else step_t)
                R2 = R @ _rot(dv[:3])
                t2 = t + dv[3:]
                c2 = _corr(_simulate((X - cx) @ R2 + t2, P, resolution, TP), D)
                if c2 > cc + 1e-4:
                    cc, R, t, improved = c2, R2, t2, True
        if not improved:
            step_a, step_t = step_a / 2, step_t / 2
            if step_t < 0.25:
                break
    return (X - cx) @ R + t, cc


def map_observables(doc, coords, data):
    """Per conformation: density simulated at the map points above the contour after rigid fitting, and the cc."""
    m = _MAPS.get(data['map'])
    if m is None:
        raise ValueError('the map is no longer loaded; load it again')
    pts, dens = _map_points(m, data['contour'])
    if len(pts) < 50:
        raise ValueError('fewer than 50 voxels above the contour level; lower it')
    rng = np.random.default_rng(1)
    if len(pts) > 6000:
        sel = rng.choice(len(pts), 6000, replace=False)
        pts, dens = pts[sel], dens[sel]
    _el, X = _heavy(doc, coords)
    O, cc = [], []
    for k in range(len(X)):
        Y, c = fit_in_map(X[k], pts, dens, data['resolution'])
        O.append(_simulate(Y, pts, data['resolution']))
        cc.append(round(_corr(O[-1], dens), 4))
    return np.array(O), dens, cc


# ------------------------------------------------------------------------------------------- one entry point
def fit(doc, coords, datasets, theta=None, prior=None):
    """Agreement of every conformation with every dataset, and the joint maximum-entropy reweighting.

    datasets: [{'kind': 'saxs'|'fret'|'map', ...parsed data}]. Returns per dataset the agreement per conformation
    (SAXS: reduced chi2 after fitting scale and background; FRET: reduced chi2; map: correlation) and, with more
    than one conformation, the reweighting."""
    K = len(doc['conformations'])
    out, md = [], []
    for d in datasets:
        kind = d.get('kind')
        if kind == 'saxs':
            O = saxs_observables(doc, coords, d)
            y, s = np.array(d['I']), np.array(d['sigma'])
            agr, curves = [], []
            for k in range(K):
                A = np.vstack([O[k], np.ones_like(O[k])]).T / s[:, None]
                (a, b), *_ = np.linalg.lstsq(A, y / s, rcond=None)
                r = (y - a * O[k] - b) / s
                agr.append(round(float(r @ r / len(y)), 4))
                curves.append([float(f'{v:.5g}') for v in a * O[k] + b])
            md.append({'O': O, 'y': y, 'sigma': s, 'scale': 'linear'})
            out.append({'kind': 'saxs', 'agreement': agr, 'measure': 'chi2', 'lower_is_better': True, 'curves': curves,
                        'rg': [round(float(_rg(doc, coords, k)), 2) for k in range(K)]})
        elif kind == 'fret':
            O = fret_observables(doc, coords, d)
            y, s = fret_targets(d)
            agr = [round(float(np.sum(((O[k] - y) / s) ** 2) / len(y)), 4) for k in range(K)]
            md.append({'O': O, 'y': y, 'sigma': s, 'scale': 'none'})
            out.append({'kind': 'fret', 'agreement': agr, 'measure': 'chi2', 'lower_is_better': True,
                        'predicted': O.round(4).tolist()})
        elif kind == 'map':
            O, dens, cc = map_observables(doc, coords, d)
            s = np.full(len(dens), max(1e-6, 0.1 * float(np.std(dens))))
            md.append({'O': O, 'y': dens, 'sigma': s, 'scale': 'linear'})
            out.append({'kind': 'map', 'agreement': cc, 'measure': 'cc', 'lower_is_better': False})
        else:
            raise ValueError(f'unknown data type {kind!r}')
    res = {'datasets': out}
    if K > 1 and md:
        res['reweighting'] = maxent(md, K, theta=theta, prior=prior)
    return res


def _rg(doc, coords, k):
    _el, X = _heavy(doc, coords)
    Y = X[k]
    return math.sqrt(float(np.mean(np.sum((Y - Y.mean(0)) ** 2, axis=1))))
