"""Ranking/sampling diagnostics for a benchmark run: per target, the RMSD of every valid candidate
versus its model score (Spearman rho), top-1 vs best-candidate gap, and Rg statistics.

Usage: PYTHONPATH=src python benchmark/analyze_ranking.py --run-id dev_probe2 [--mode oracle]
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'src'))
from rnafold3d.metrics import compare, read_numbered_pdb   # noqa: E402


def spearman(x, y):
    if len(x) < 3:
        return None
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def rg(path):
    res = read_numbered_pdb(path)
    c = np.array([r.atoms["C4'"] for r in res.values() if "C4'" in r.atoms])
    return float(np.sqrt(np.mean(np.sum((c - c.mean(0)) ** 2, axis=1))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--mode', default='oracle')
    a = ap.parse_args()
    run = ROOT / 'runs' / a.run_id / a.mode
    out = []
    for tdir in sorted(p for p in run.iterdir() if p.is_dir()):
        sid = tdir.name
        ref = ROOT / 'data' / 'samples' / sid / 'reference.pdb'
        scores = {}
        if not (tdir / 'metadata' / 'scores.csv').exists() or not (tdir / 'candidates').exists():
            continue
        with (tdir / 'metadata' / 'scores.csv').open() as f:
            for row in csv.DictReader(f):
                if row['status'] == 'valid':
                    scores[int(row['trajectory'])] = float(row['total_score'])
        xs, ys, rgs = [], [], []
        for p in sorted((tdir / 'candidates').glob('traj_*.pdb')):
            k = int(p.stem.split('_')[1])
            m = compare(ref, p)
            if m['c4_rmsd'] is None or k not in scores:
                continue
            xs.append(scores[k])
            ys.append(m['c4_rmsd'])
            rgs.append(rg(p))
        if not xs:
            out.append({'sample_id': sid, 'n_valid': 0})
            continue
        i_top = int(np.argmin(xs))
        out.append({'sample_id': sid, 'n_valid': len(xs), 'top_score_rmsd': round(ys[i_top], 2),
                    'best_rmsd': round(min(ys), 2), 'median_rmsd': round(float(np.median(ys)), 2),
                    'spearman_score_rmsd': None if spearman(xs, ys) is None else round(spearman(xs, ys), 2),
                    'native_rg': round(rg(ref), 1), 'model_rg_median': round(float(np.median(rgs)), 1)})
    for r in out:
        print(json.dumps(r))


if __name__ == '__main__':
    main()
