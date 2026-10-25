"""Features and labels for the RCV confidence model.

For every model of a run (engine trajectories in candidates/, exported conformers, and the reference itself):
per-nucleotide features (rnafold3d.qa.residue_features) and labels (all-atom lDDT per nucleotide against the
reference), plus per-model C4' RMSD, global lDDT, engine score and rank.

    # training set: decoys of tools/qa_model targets (benchmark/runs/qa_train_v1)
    PYTHONPATH=src python tools/qa_model/make_dataset.py --data benchmark/data_qa --run benchmark/runs/qa_train_v1 \
        --out benchmark/runs/qa_train_v1/qa_dataset.npz
    # evaluation set: the benchmark (never used for training or model selection)
    PYTHONPATH=src python tools/qa_model/make_dataset.py --data benchmark/data --run benchmark/runs/v2_0_1_pilot \
        --out benchmark/runs/v2_0_1_pilot/qa_eval.npz --conformers
"""
from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import qa                                   # noqa: E402
from rnafold3d.metrics import read_numbered_pdb, c4_rmsd   # noqa: E402


def one_target(args):
    sample_dir, run_dirs, conformers, natives = args
    sample_dir = Path(sample_dir)
    seq = (sample_dir / 'sequence.fasta').read_text().split('\n')[1].strip()
    ref = read_numbered_pdb(sample_dir / 'reference.pdb')
    rows = []
    jobs = []
    if natives:
        jobs.append(('native', 'reference', sample_dir / 'reference.pdb', None, None))
    for mode, rd in run_dirs:
        rd = Path(rd)
        man = rd / 'metadata' / 'manifest.json'
        if not man.exists():
            continue
        m = json.loads(man.read_text())
        if conformers:
            for c in m['conformers']:
                jobs.append(('conformer', mode, rd / c['files']['pdb'], c.get('total_score'), c.get('rank')))
        else:
            for p in sorted((rd / 'candidates').glob('traj_*.pdb')):
                jobs.append(('candidate', mode, p, None, None))
    for kind, mode, path, score, rank in jobs:
        try:
            mod = read_numbered_pdb(path)
            keys = sorted(mod)
            if len(keys) != len(seq):
                continue
            residues = [mod[k] for k in keys]
            X = qa.residue_features(seq, residues)
            lab = qa.lddt_per_residue(ref, mod) if kind != 'native' else {k: 1.0 for k in keys}
            y = np.array([lab.get(k, np.nan) for k in keys], np.float32)
            rmsd = c4_rmsd(ref, mod)[0] if kind != 'native' else 0.0
            rows.append({'target': sample_dir.name, 'kind': kind, 'mode': mode, 'path': str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
                         'score': score, 'rank': rank, 'rmsd': rmsd, 'lddt': float(np.nanmean(y)) if np.isfinite(y).any() else None,
                         'X': X.astype(np.float32), 'y': y})
        except Exception as exc:                       # a broken file must not stop the dataset
            print('skip', path, exc, flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--run', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--modes', default='independent,oracle')
    ap.add_argument('--conformers', action='store_true', help='use the exported conformers (with engine scores) instead of all trajectories')
    ap.add_argument('--no-natives', action='store_true')
    ap.add_argument('--workers', type=int, default=2)
    a = ap.parse_args()
    data, run = Path(a.data), Path(a.run)
    ds = json.loads((data / 'dataset_v1.json').read_text())
    tasks = []
    for s in ds['samples']:
        sd = data / 'samples' / s['sample_id']
        runs = [(m, run / m / s['sample_id']) for m in a.modes.split(',') if (run / m / s['sample_id'] / 'metadata' / 'manifest.json').exists()]
        if runs:
            tasks.append((str(sd), [(m, str(r)) for m, r in runs], a.conformers, not a.no_natives))
    print(f'{len(tasks)} targets with finished runs', flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, rr in enumerate(ex.map(one_target, tasks), 1):
            rows += rr
            if k % 10 == 0:
                print(f'  {k}/{len(tasks)} targets, {len(rows)} models', flush=True)
    X = np.concatenate([r['X'] for r in rows]) if rows else np.zeros((0, len(qa.feature_names())), np.float32)
    y = np.concatenate([r['y'] for r in rows]) if rows else np.zeros(0, np.float32)
    model_of = np.concatenate([np.full(len(r['y']), k) for k, r in enumerate(rows)]) if rows else np.zeros(0, int)
    meta = [{k: v for k, v in r.items() if k not in ('X', 'y')} for r in rows]
    np.savez_compressed(a.out, X=X, y=y, model=model_of, meta=json.dumps(meta), features=json.dumps(qa.feature_names()),
                        feature_version=qa.FEATURE_VERSION)
    print(f'{len(rows)} models, {len(y)} nucleotides ({np.isfinite(y).sum()} labelled) -> {a.out}')


if __name__ == '__main__':
    main()
