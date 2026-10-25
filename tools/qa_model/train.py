"""Train the RCV confidence model (rnafold3d/qa.py) on the dataset of make_dataset.py.

Network: standardised features -> 48 tanh -> 24 tanh -> sigmoid (predicted lDDT per nucleotide); loss: weighted mean
squared error (every target weighs the same, natives half) + L2. Adam, numpy only. Grouped 5-fold cross-validation
by target chooses the number of epochs; the final model is trained on all targets and written as an .npz that
rnafold3d.qa loads.

    PYTHONPATH=src python tools/qa_model/train.py --data benchmark/runs/qa_train_v1/qa_dataset.npz \
        --out src/rnafold3d/data/qa_v1.npz
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import qa  # noqa: E402

H1, H2 = 48, 24


def load(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z['meta']))
    X, y, model = z['X'].astype(np.float64), z['y'].astype(np.float64), z['model']
    keep = np.isfinite(y)
    return X[keep], y[keep], model[keep], meta


def weights_for(model, meta):
    """Every target weighs the same; within a target every model the same; natives half."""
    target = np.array([meta[m]['target'] for m in model])
    w = np.ones(len(model))
    for t in np.unique(target):
        sel = target == t
        mods = np.unique(model[sel])
        for m in mods:
            s2 = model == m
            w[s2] = (0.5 if meta[m]['kind'] == 'native' else 1.0) / (len(mods) * s2.sum())
    return w / w.mean(), target


def init(F, rng):
    return {'W1': rng.normal(0, 1 / np.sqrt(F), (F, H1)), 'b1': np.zeros(H1),
            'W2': rng.normal(0, 1 / np.sqrt(H1), (H1, H2)), 'b2': np.zeros(H2),
            'W3': rng.normal(0, 1 / np.sqrt(H2), (H2, 1)), 'b3': np.zeros(1)}


def forward(p, Z):
    A1 = Z @ p['W1'] + p['b1']; T1 = np.tanh(A1)
    A2 = T1 @ p['W2'] + p['b2']; T2 = np.tanh(A2)
    out = 1 / (1 + np.exp(-(T2 @ p['W3'] + p['b3']).ravel()))
    return out, (T1, T2)


def fit(Z, y, w, epochs, rng, lr=3e-3, l2=1e-4, batch=512, Zv=None, yv=None, wv=None):
    p = init(Z.shape[1], rng)
    m = {k: np.zeros_like(v) for k, v in p.items()}
    v2 = {k: np.zeros_like(v) for k, v in p.items()}
    b1, b2, eps, t = 0.9, 0.999, 1e-8, 0
    curve = []
    n = len(y)
    for ep in range(epochs):
        order = rng.permutation(n)
        for s in range(0, n, batch):
            idx = order[s:s + batch]
            Zb, yb, wb = Z[idx], y[idx], w[idx]
            out, (T1, T2) = forward(p, Zb)
            sw = wb.sum()
            dA3 = (2 * wb * (out - yb) * out * (1 - out) / sw)[:, None]
            g = {'W3': T2.T @ dA3 + 2 * l2 * p['W3'], 'b3': dA3.sum(axis=0)}
            dA2 = (dA3 @ p['W3'].T) * (1 - T2 ** 2)
            g['W2'] = T1.T @ dA2 + 2 * l2 * p['W2']; g['b2'] = dA2.sum(axis=0)
            dA1 = (dA2 @ p['W2'].T) * (1 - T1 ** 2)
            g['W1'] = Zb.T @ dA1 + 2 * l2 * p['W1']; g['b1'] = dA1.sum(axis=0)
            t += 1
            for k in p:
                m[k] = b1 * m[k] + (1 - b1) * g[k]
                v2[k] = b2 * v2[k] + (1 - b2) * g[k] ** 2
                p[k] -= lr * (m[k] / (1 - b1 ** t)) / (np.sqrt(v2[k] / (1 - b2 ** t)) + eps)
        if Zv is not None:
            ov, _ = forward(p, Zv)
            curve.append(float(np.sum(wv * (ov - yv) ** 2) / wv.sum()))
    return p, curve


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--epochs', type=int, default=40)
    ap.add_argument('--folds', type=int, default=5)
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    X, y, model, meta = load(a.data)
    w, target = weights_for(model, meta)
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd < 1e-6] = 1.0
    Z = (X - mu) / sd
    rng = np.random.default_rng(a.seed)
    targets = np.unique(target)
    rng.shuffle(targets)
    folds = np.array_split(targets, a.folds)
    curves, cv_r, cv_rank = [], [], []
    for k, ft in enumerate(folds):
        va = np.isin(target, ft)
        p, curve = fit(Z[~va], y[~va], w[~va], a.epochs, np.random.default_rng(a.seed + k), Zv=Z[va], yv=y[va], wv=w[va])
        curves.append(curve)
        out, _ = forward(p, Z[va])
        cv_r.append(float(np.corrcoef(out, y[va])[0, 1]))
        # ranking within targets: global prediction vs global lDDT (decoys only)
        for t in ft:
            mods = [mm for mm in np.unique(model[target == t]) if meta[mm]['kind'] != 'native']
            if len(mods) < 3:
                continue
            pg = [out[model[va] == mm].mean() for mm in mods]
            tg = [meta[mm]['lddt'] for mm in mods]
            cv_rank.append(float(np.corrcoef(np.argsort(np.argsort(pg)), np.argsort(np.argsort(tg)))[0, 1]))
        print(f'fold {k + 1}: per-nucleotide r {cv_r[-1]:.3f}, val loss {curve[-1]:.4f} (best epoch {int(np.argmin(curve)) + 1})', flush=True)
    mean_curve = np.mean(curves, axis=0)
    best_ep = int(np.argmin(mean_curve)) + 1
    print(f'CV: per-nucleotide r {np.mean(cv_r):.3f} ± {np.std(cv_r):.3f}; within-target ranking Spearman {np.nanmean(cv_rank):.3f} '
          f'(n={len(cv_rank)}); best epochs {best_ep}', flush=True)
    p, _ = fit(Z, y, w, best_ep, np.random.default_rng(a.seed))
    info = {'created': datetime.date.today().isoformat(), 'feature_version': qa.FEATURE_VERSION, 'features': qa.feature_names(),
            'hidden': [H1, H2], 'epochs': best_ep, 'training': {'data': str(a.data), 'targets': int(len(targets)),
            'models': int(len(np.unique(model))), 'nucleotides': int(len(y))},
            'cv': {'per_nucleotide_r': round(float(np.mean(cv_r)), 4), 'per_nucleotide_r_sd': round(float(np.std(cv_r)), 4),
                   'ranking_spearman': round(float(np.nanmean(cv_rank)), 4), 'folds': a.folds}}
    np.savez(a.out, mean=mu, std=sd, feature_version=np.array(qa.FEATURE_VERSION), info=json.dumps(info), **p)
    print(json.dumps(info['cv']), '->', a.out)


if __name__ == '__main__':
    main()
