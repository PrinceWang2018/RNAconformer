"""Evaluate the RCV confidence model on the benchmark (models the network never saw, families excluded from training).

Reports, for the benchmark's dev and test splits separately:
  * per-nucleotide Pearson r between predicted and true lDDT (decoys);
  * within each target, Spearman correlation of the predicted global confidence with the true global lDDT;
  * top-1 selection: the conformer ranked first by the engine score, by the confidence model, and by a blend
    (z-scores, weight chosen on the dev split only), as median C4' RMSD and lDDT of the chosen model;
  * natives: predicted confidence of the experimental reference structures.

    PYTHONPATH=src python tools/qa_model/evaluate.py --eval benchmark/runs/v2_0_1_pilot/qa_eval.npz \
        --weights src/rnafold3d/data/qa_v1.npz --out docs/benchmark_v1/qa_v1_eval.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import qa  # noqa: E402


def spearman(a, b):
    if len(a) < 3:
        return None
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    if ra.std() == 0 or rb.std() == 0:
        return None
    return float(np.corrcoef(ra, rb)[0, 1])


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / (v.std() or 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--eval', required=True)
    ap.add_argument('--weights', required=True)
    ap.add_argument('--dataset', default=str(ROOT / 'benchmark' / 'data' / 'dataset_v1.json'))
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    zf = np.load(a.eval, allow_pickle=False)
    meta = json.loads(str(zf['meta']))
    X, y, model = zf['X'].astype(np.float64), zf['y'].astype(np.float64), zf['model']
    w = qa.load_weights(a.weights)
    pred = qa.forward(X, w)
    split = {s['sample_id']: s['split'] for s in json.loads(Path(a.dataset).read_text())['samples']}
    for k, m in enumerate(meta):
        sel = model == k
        m['pred'] = float(pred[sel].mean())
    out = {'weights': json.loads(str(w['info'])) if 'info' in w else None}
    # natives
    nat = [m['pred'] for m in meta if m['kind'] == 'native']
    dec = [m['pred'] for m in meta if m['kind'] != 'native']
    out['natives'] = {'mean_confidence': round(float(np.mean(nat)), 3), 'decoys_mean_confidence': round(float(np.mean(dec)), 3)}
    for sp in ('dev', 'test', 'all'):
        tg = sorted({m['target'] for m in meta if sp == 'all' or split.get(m['target']) == sp})
        dmask = np.array([meta[k]['kind'] != 'native' and (sp == 'all' or split.get(meta[k]['target']) == sp) for k in model])
        fin = dmask & np.isfinite(y)
        r = float(np.corrcoef(pred[fin], y[fin])[0, 1]) if fin.sum() > 2 else None
        rho, sel = [], {'engine': [], 'confidence': [], 'best': [], 'blend': {}}
        for t in tg:
            for mode in ('independent', 'oracle'):
                ms = [m for m in meta if m['target'] == t and m['mode'] == mode and m['kind'] == 'conformer' and m['rmsd'] is not None]
                if len(ms) < 2:
                    continue
                s = spearman([m['pred'] for m in ms], [m['lddt'] for m in ms])
                if s is not None:
                    rho.append(s)
                eng = min(ms, key=lambda m: m['rank'])
                con = max(ms, key=lambda m: m['pred'])
                best = min(ms, key=lambda m: m['rmsd'])
                sel['engine'].append(eng); sel['confidence'].append(con); sel['best'].append(best)
                for lam in (0.25, 0.5, 0.75):
                    zs = lam * z([m['pred'] for m in ms]) - (1 - lam) * z([m['score'] for m in ms])
                    sel['blend'].setdefault(lam, []).append(ms[int(np.argmax(zs))])
        summ = lambda L: {'median_rmsd': round(float(np.median([m['rmsd'] for m in L])), 2), 'mean_rmsd': round(float(np.mean([m['rmsd'] for m in L])), 2),
                          'median_lddt': round(float(np.median([m['lddt'] for m in L])), 3), 'n': len(L)} if L else None
        res = {'targets': len(tg), 'per_nucleotide_r': round(r, 3) if r is not None else None,
               'ranking_spearman_mean': round(float(np.mean(rho)), 3) if rho else None, 'ranking_cases': len(rho),
               'top1_engine': summ(sel['engine']), 'top1_confidence': summ(sel['confidence']), 'best_of_models': summ(sel['best']),
               'top1_blend': {str(k): summ(v) for k, v in sel['blend'].items()}}
        if sel['engine']:
            d = np.array([c['rmsd'] - e['rmsd'] for c, e in zip(sel['confidence'], sel['engine'])])
            res['confidence_vs_engine'] = {'better': int((d < -0.5).sum()), 'worse': int((d > 0.5).sum()), 'same': int((np.abs(d) <= 0.5).sum())}
        out[sp] = res
    txt = json.dumps(out, indent=1)
    print(txt)
    if a.out:
        Path(a.out).write_text(txt)


if __name__ == '__main__':
    main()
