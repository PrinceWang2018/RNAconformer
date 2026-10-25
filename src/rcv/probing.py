"""Chemical probing data (SHAPE, DMS) for an ensemble.

* Reading reactivity profiles: two-column text (position, value; RNAstructure .shape, many lab formats), ShapeMapper
  .map / profile files (position, reactivity, stderr, ...), CSV with a header (a column named reactivity, value,
  shape or dms), or one value per line. Missing values: -999, NaN, NA or empty.
* Normalisation (optional): the 2-8 % rule of the Weeks lab - values above the boxplot outlier limit
  (Q3 + 1.5 IQR, at most the top 10 %) are set aside and every value is divided by the mean of the next 8 %.
  DMS: A and C only.
* Agreement of each conformation with the data: the area under the ROC curve of reactivity as a predictor of
  unpaired nucleotides (1 = perfect, 0.5 = no information), using the base pairs annotated in 3D.
* Ensemble reweighting (maximum entropy, Bayesian/maximum-entropy reweighting after Hummer & Köfinger 2015 and
  Bottaro et al. 2020): weights w_k of the conformations minimise chi2(w)/2 - theta * S(w || w0), where the
  predicted reactivity of nucleotide i is alpha * <u_i>_w + beta (u = 1 if unpaired in conformation k; alpha >= 0
  and beta fitted by least squares), S is the relative entropy to the prior w0 (uniform or given) and theta sets how
  far the data may move the weights. Reported: weights, chi2 before and after, the effective number of
  conformations (exp of the entropy) and the agreement of the reweighted ensemble.
"""
from __future__ import annotations

import math
import re

import numpy as np

MISSING = {'', 'nan', 'na', 'n/a', '-999', '-999.0', 'none', 'null'}


def _num(tok):
    t = tok.strip().lower()
    if t in MISSING:
        return None
    v = float(t)
    if not math.isfinite(v) or v <= -500:
        return None
    return v


def parse(text, n):
    """Reactivities for positions 1..n from text; returns (values with None for missing, format name, notes)."""
    lines = [ln for ln in re.split(r'\r?\n', text or '') if ln.strip() and not ln.lstrip().startswith(('#', ';', '>'))]
    if not lines:
        raise ValueError('the file has no data')
    split = lambda ln: [t for t in re.split(r'[\s,;]+', ln.strip()) if t != '']
    head = split(lines[0])
    notes = []
    vals = [None] * n

    def is_num(t):
        try:
            _num(t)
            return True
        except ValueError:
            return False
    col = None
    body = lines
    if not is_num(head[0]):                     # a header names the columns; data lines start with a position
        names = [t.lower() for t in head]
        for want in ('reactivity', 'norm_profile', 'shape', 'dms', 'value', 'profile', 'reactivity_profile'):
            col = next((k for k, nm in enumerate(names) if want in nm), None)
            if col is not None:
                break
        pos_col = next((k for k, nm in enumerate(names) if nm in ('pos', 'position', 'nucleotide', 'nt', 'seqnum', 'index')), 0)
        if col is None:
            raise ValueError(f'no reactivity column in the header ({", ".join(head)})')
        body = lines[1:]
        fmt = f'table with header (column "{head[col]}")'
        for ln in body:
            t = split(ln)
            if len(t) <= max(col, pos_col):
                continue
            try:
                p = int(float(t[pos_col]))
            except ValueError:
                continue
            if 1 <= p <= n:
                vals[p - 1] = _num(t[col])
        return vals, fmt, notes
    width = len(head)
    if width == 1:
        nums = [_num(split(ln)[0]) for ln in body]
        if len(nums) != n:
            notes.append(f'{len(nums)} values for {n} nucleotides; the first {min(n, len(nums))} are used')
        for i, v in enumerate(nums[:n]):
            vals[i] = v
        return vals, 'one value per line', notes
    # position + value (+ more columns: ShapeMapper .map has position, reactivity, stderr, nucleotide)
    fmt = 'position and value' if width == 2 else f'{width} columns (position, reactivity, ...)'
    outside = 0
    for ln in body:
        t = split(ln)
        if len(t) < 2:
            continue
        try:
            p = int(float(t[0]))
        except ValueError:
            continue
        if 1 <= p <= n:
            vals[p - 1] = _num(t[1])
        else:
            outside += 1
    if outside:
        notes.append(f'{outside} positions outside 1..{n} ignored')
    return vals, fmt, notes


def normalise(vals, seq, kind='SHAPE'):
    """2-8 % normalisation; returns (normalised values, factor)."""
    use = [i for i, v in enumerate(vals) if v is not None and (kind != 'DMS' or seq[i] in 'AC')]
    x = np.array([max(0.0, vals[i]) for i in use], float)
    if len(x) < 10:
        raise ValueError('too few values to normalise (need at least 10)')
    q1, q3 = np.percentile(x, [25, 75])
    limit = q3 + 1.5 * (q3 - q1)
    xs = np.sort(x)[::-1]
    n_out = min(int(np.sum(xs > limit)), int(0.1 * len(xs)))
    if n_out == 0:
        n_out = int(round(0.02 * len(xs)))
    nxt = xs[n_out:n_out + max(1, int(round(0.08 * len(xs))))]
    factor = float(nxt.mean()) if len(nxt) and nxt.mean() > 0 else 1.0
    out = [None if v is None or (kind == 'DMS' and seq[i] not in 'AC') else v / factor for i, v in enumerate(vals)]
    return out, factor


def unpaired_matrix(doc):
    """u[k][i] = 1 if nucleotide i is unpaired in conformation k (canonical pairs annotated in 3D, pseudoknots too)."""
    n = len(doc['sequence'])
    U = np.ones((len(doc['conformations']), n))
    for k, c in enumerate(doc['conformations']):
        for p in c.get('pairs', []):
            U[k, p[0]] = U[k, p[1]] = 0.0
    return U


def auc(scores, labels):
    """Area under the ROC curve (Mann-Whitney), labels 1 = positive."""
    s = np.asarray(scores, float)
    y = np.asarray(labels, int)
    pos, neg = s[y == 1], s[y == 0]
    if len(pos) == 0 or len(neg) == 0:
        return None
    order = np.argsort(s, kind='mergesort')
    ranks = np.empty(len(s))
    ranks[order] = np.arange(1, len(s) + 1)
    # ties: average ranks
    ss = s[order]
    i = 0
    while i < len(ss):
        j = i
        while j + 1 < len(ss) and ss[j + 1] == ss[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j + 2) / 2
        i = j + 1
    return float((ranks[y == 1].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def _mask(vals, seq, kind):
    return np.array([v is not None and (kind != 'DMS' or seq[i] in 'AC') for i, v in enumerate(vals)])


def agreement(doc, vals, kind='SHAPE'):
    seq = doc['sequence']
    U = unpaired_matrix(doc)
    m = _mask(vals, seq, kind)
    y = np.array([v if v is not None else 0.0 for v in vals])
    return [None if m.sum() < 5 else (lambda a: round(a, 4) if a is not None else None)(auc(y[m], U[k][m])) for k in range(len(U))]


def reweight(doc, vals, kind='SHAPE', theta=None, prior=None, sigma=0.3):
    """Maximum-entropy reweighting of the conformations to fit the reactivities."""
    from scipy.optimize import minimize
    seq = doc['sequence']
    U = unpaired_matrix(doc)
    m = _mask(vals, seq, kind)
    if m.sum() < 5:
        raise ValueError('too few nucleotides with data')
    y = np.array([vals[i] for i in range(len(seq)) if m[i]], float)
    y = np.clip(y, 0.0, 4.0)
    Um = U[:, m]
    K = len(U)
    w0 = np.full(K, 1.0 / K) if prior is None else np.asarray(prior, float) / np.sum(prior)
    w0 = np.clip(w0, 1e-12, None)
    w0 /= w0.sum()
    theta = float(theta) if theta else max(1.0, 0.1 * len(y))

    def fit_scale(mu):
        A = np.vstack([mu, np.ones_like(mu)]).T
        (alpha, beta), *_ = np.linalg.lstsq(A, y, rcond=None)
        if alpha < 0:                       # reactivity must not fall with unpairing
            alpha, beta = 0.0, float(y.mean())
        return float(alpha), float(beta)

    def chi2_of(w):
        mu = w @ Um
        a, b = fit_scale(mu)
        r = y - a * mu - b
        return float((r ** 2).sum() / sigma ** 2), a, b, mu

    def f(z):
        z = z - z.max()
        w = w0 * np.exp(z)
        w /= w.sum()
        c2, a, b, mu = chi2_of(w)
        S = -float(np.sum(w * np.log(w / w0)))
        F = 0.5 * c2 - theta * S
        r = y - a * mu - b
        g_w = -(r @ (a * Um.T)) / sigma ** 2 + theta * (np.log(w / w0) + 1)
        grad = w * (g_w - float(w @ g_w))
        return F, grad
    res = minimize(f, np.zeros(K), jac=True, method='L-BFGS-B', options={'maxiter': 500})
    z = res.x - res.x.max()
    w = w0 * np.exp(z)
    w /= w.sum()
    c0, a0, b0, mu0 = chi2_of(w0)
    c1, a1, b1, mu1 = chi2_of(w)
    S = -float(np.sum(w * np.log(w / w0)))
    # agreement (AUC) of the prior and the reweighted ensemble: expected unpairing vs data
    return {'weights': [round(float(x), 6) for x in w], 'prior': [round(float(x), 6) for x in w0],
            'chi2_prior': round(c0 / len(y), 4), 'chi2': round(c1 / len(y), 4), 'theta': theta, 'sigma': sigma,
            'alpha': round(a1, 4), 'beta': round(b1, 4), 'n_data': int(len(y)),
            'phi_eff': round(float(math.exp(S)), 4), 'n_eff_conformations': round(float(math.exp(-np.sum(w * np.log(w)))), 2),
            'auc_prior': auc(y, (mu0 > 0.5).astype(int)) if 0 < (mu0 > 0.5).sum() < len(mu0) else None,
            'auc': auc(y, (mu1 > 0.5).astype(int)) if 0 < (mu1 > 0.5).sum() < len(mu1) else None,
            'converged': bool(res.success)}
