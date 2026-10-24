"""Run the v2 engine on the frozen pilot dataset and evaluate against experimental references.

Modes
  independent : FASTA only. ViennaRNA predicts the 2D hypotheses (no information from the target).
  oracle      : the nested canonical 2D annotated from the reference 3D is given as input.
                ORACLE results measure 3D assembly given correct 2D and must never be reported
                as sequence-only predictions.

Predefined (before running) success criteria: best-of-K C4' RMSD <= 5 A and <= 10 A;
eRMSD <= 1.0 (native-like base arrangement); INF-WC reported as is. Failed targets count as
failures in every aggregate.

Usage:
  PYTHONPATH=src python benchmark/run_benchmark.py --run-id v2_dev --split dev --modes oracle,independent
"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')   # one BLAS thread per process; parallelism is across trajectories

import argparse
import csv
import datetime
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'src'))

from rnafold3d.predict import predict, PredictConfig, read_dbn, ENGINE, ENGINE_VERSION   # noqa: E402
from rnafold3d.metrics import compare                                                   # noqa: E402
from rnafold3d.ss import vienna_version                                                 # noqa: E402
from rnafold3d.search import ScoreWeights                                               # noqa: E402

THRESHOLDS = {'c4_rmsd_5A': 5.0, 'c4_rmsd_10A': 10.0, 'ermsd_1': 1.0}


def evaluate_run(sample_dir, run_dir):
    man = json.loads((run_dir / 'metadata' / 'manifest.json').read_text())
    rows = []
    for c in man['conformers']:
        m = compare(sample_dir / 'reference.pdb', run_dir / c['files']['pdb'])
        rows.append({'id': c['id'], 'rank': c['rank'], 'ss_id': c['secondary_structure_id'],
                     'total_score': c['total_score'], **m})
    return man, rows


def candidate_rmsds(sample_dir, run_dir):
    out = []
    cdir = run_dir / 'candidates'
    if cdir.exists():
        for p in sorted(cdir.glob('traj_*.pdb')):
            out.append(compare(sample_dir / 'reference.pdb', p)['c4_rmsd'])
    return [x for x in out if x is not None]


def summarize_target(meta, mode, man, rows, wall, cand=None):
    rm = [r['c4_rmsd'] for r in rows if r['c4_rmsd'] is not None]
    er = [r['ermsd'] for r in rows if r['ermsd'] is not None]
    inf = [r['inf_wc'] for r in rows]
    top = rows[0] if rows else None
    out = {'sample_id': meta['sample_id'], 'split': meta['split'], 'mode': mode, 'length': meta['length'],
           'stratum': meta['stratum'], 'family': meta['family'], 'n_models': len(rows),
           'trajectories': man['statistics']['trajectories'], 'valid_trajectories': man['statistics']['valid'],
           'top1_c4_rmsd': top['c4_rmsd'] if top else None, 'best_c4_rmsd': min(rm) if rm else None,
           'mean_c4_rmsd': round(float(np.mean(rm)), 3) if rm else None,
           'top1_ermsd': top['ermsd'] if top else None, 'best_ermsd': min(er) if er else None,
           'top1_inf_wc': top['inf_wc'] if top else None, 'best_inf_wc': max(inf) if inf else None,
           'top1_ss_id': top['ss_id'] if top else None, 'wall_time_s': round(wall, 1),
           'yield': bool(rows),
           'best_candidate_c4_rmsd': min(cand) if cand else None,
           'median_candidate_c4_rmsd': round(float(np.median(cand)), 3) if cand else None}
    return out


def aggregate(records):
    agg = {}
    for mode in sorted({r['mode'] for r in records}):
        rs = [r for r in records if r['mode'] == mode]
        n = len(rs)

        def med(key):
            v = [r[key] for r in rs if r[key] is not None]
            return round(float(np.median(v)), 2) if v else None

        def frac(key, thr, better='low'):
            ok = [r for r in rs if r[key] is not None and (r[key] <= thr if better == 'low' else r[key] >= thr)]
            return round(len(ok) / n, 3) if n else None
        agg[mode] = {
            'targets': n, 'yield': round(sum(r['yield'] for r in rs) / n, 3) if n else None,
            'median_top1_c4_rmsd': med('top1_c4_rmsd'), 'median_best_c4_rmsd': med('best_c4_rmsd'),
            'median_top1_ermsd': med('top1_ermsd'), 'median_best_ermsd': med('best_ermsd'),
            'median_top1_inf_wc': med('top1_inf_wc'), 'median_best_inf_wc': med('best_inf_wc'),
            'success_best_c4_le_5A': frac('best_c4_rmsd', 5.0), 'success_best_c4_le_10A': frac('best_c4_rmsd', 10.0),
            'success_top1_c4_le_10A': frac('top1_c4_rmsd', 10.0), 'success_best_ermsd_le_1': frac('best_ermsd', 1.0),
            'median_wall_time_s': med('wall_time_s'),
            'by_stratum': {},
        }
        for st in sorted({r['stratum'] for r in rs}):
            sub = [r for r in rs if r['stratum'] == st]
            v = [r['best_c4_rmsd'] for r in sub if r['best_c4_rmsd'] is not None]
            t = [r['top1_c4_rmsd'] for r in sub if r['top1_c4_rmsd'] is not None]
            agg[mode]['by_stratum'][st] = {'targets': len(sub), 'with_models': len(v),
                                           'median_best_c4_rmsd': round(float(np.median(v)), 2) if v else None,
                                           'median_top1_c4_rmsd': round(float(np.median(t)), 2) if t else None}
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    ap.add_argument('--split', default='all', choices=['dev', 'test', 'all'])
    ap.add_argument('--samples', default='')
    ap.add_argument('--modes', default='oracle,independent')
    ap.add_argument('--models', type=int, default=5)
    ap.add_argument('--trajectories', type=int, default=16)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--steps', type=int, default=0)
    ap.add_argument('--data', default=str(ROOT / 'data'))
    ap.add_argument('--set', default='', help='score weight overrides, e.g. stacking=0,compaction=0')
    a = ap.parse_args()
    data = Path(a.data)
    ds = json.loads((data / 'dataset_v1.json').read_text())
    samples = [s for s in ds['samples'] if a.split == 'all' or s['split'] == a.split]
    if a.samples:
        want = set(a.samples.split(','))
        samples = [s for s in samples if s['sample_id'] in want or s['pdb_id'] in want]
    samples.sort(key=lambda s: s['length'])
    run_root = ROOT / 'runs' / a.run_id
    run_root.mkdir(parents=True, exist_ok=True)
    overrides = {}
    for item in filter(None, a.set.split(',')):
        k, v = item.split('=')
        default = getattr(ScoreWeights(), k)
        overrides[k] = (v.lower() in ('1', 'true', 'yes')) if isinstance(default, bool) else float(v)
    weights = ScoreWeights(**overrides)
    cfg = {'engine': ENGINE, 'engine_version': ENGINE_VERSION, 'weight_overrides': overrides,
           'models': a.models, 'trajectories': a.trajectories,
           'seed': a.seed, 'steps': a.steps, 'workers': a.workers, 'modes': a.modes, 'split': a.split,
           'vienna': vienna_version(), 'python': platform.python_version(), 'thresholds': THRESHOLDS,
           'started': datetime.datetime.now().isoformat(timespec='seconds')}
    (run_root / 'config.json').write_text(json.dumps(cfg, indent=1))
    records, per_model = [], []
    for s in samples:
        sdir = data / 'samples' / s['sample_id']
        seq = (sdir / 'sequence.fasta').read_text().split('\n')[1].strip()
        for mode in a.modes.split(','):
            out = run_root / mode / s['sample_id']
            done = out / 'metadata' / 'manifest.json'
            t0 = time.time()
            if not done.exists():
                dbn = read_dbn(sdir / 'oracle_2d.dbn', len(seq)) if mode == 'oracle' else None
                pc = PredictConfig(n_models=a.models, trajectories=a.trajectories, seed=a.seed,
                                   workers=a.workers, anneal_steps=a.steps, save_candidates=True, weights=weights)
                try:
                    predict(seq, out, name=s['sample_id'], dbn=dbn, config=pc, log=None)
                except Exception as exc:
                    print('ERROR', s['sample_id'], mode, exc)
                    records.append({'sample_id': s['sample_id'], 'split': s['split'], 'mode': mode, 'length': s['length'],
                                    'stratum': s['stratum'], 'family': s['family'], 'n_models': 0, 'yield': False,
                                    'error': str(exc), **{k: None for k in ('top1_c4_rmsd', 'best_c4_rmsd', 'mean_c4_rmsd',
                                    'top1_ermsd', 'best_ermsd', 'top1_inf_wc', 'best_inf_wc', 'top1_ss_id')},
                                    'trajectories': a.trajectories, 'valid_trajectories': 0, 'wall_time_s': time.time() - t0})
                    continue
            wall = json.loads(done.read_text())['statistics']['wall_time_s']
            man, rows = evaluate_run(sdir, out)
            rec = summarize_target(s, mode, man, rows, wall, candidate_rmsds(sdir, out))
            records.append(rec)
            for r in rows:
                per_model.append({'sample_id': s['sample_id'], 'mode': mode, **{k: v for k, v in r.items() if k != 'pairs'}})
            print(f"{s['sample_id']:8s} {mode:11s} n={s['length']:3d} models={rec['n_models']} valid={rec['valid_trajectories']}/{rec['trajectories']} "
                  f"top1={rec['top1_c4_rmsd']} best={rec['best_c4_rmsd']} eRMSD(best)={rec['best_ermsd']} "
                  f"INF(top1)={rec['top1_inf_wc']} cand(best/median)={rec['best_candidate_c4_rmsd']}/{rec['median_candidate_c4_rmsd']} "
                  f"t={rec['wall_time_s']}s", flush=True)
    rep = run_root / 'reports'
    rep.mkdir(exist_ok=True)
    keys = list(records[0].keys()) if records else []
    keys = sorted({k for r in records for k in r}, key=lambda k: (k not in keys, keys.index(k) if k in keys else 0))
    with (rep / 'per_target.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(records)
    if per_model:
        with (rep / 'per_model.csv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(per_model[0].keys()))
            w.writeheader()
            w.writerows(per_model)
    agg = aggregate(records)
    (rep / 'aggregate.json').write_text(json.dumps({'config': cfg, 'aggregate': agg,
                                                    'finished': datetime.datetime.now().isoformat(timespec='seconds')}, indent=1))
    failures = [r for r in records if not r['yield']]
    with (rep / 'failures.csv').open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['sample_id', 'mode', 'length', 'valid_trajectories', 'trajectories', 'error'])
        for r in failures:
            w.writerow([r['sample_id'], r['mode'], r['length'], r.get('valid_trajectories'), r.get('trajectories'), r.get('error', '')])
    print(json.dumps(agg, indent=1))


if __name__ == '__main__':
    main()
