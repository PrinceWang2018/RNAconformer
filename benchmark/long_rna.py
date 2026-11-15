"""Long RNAs: domain assembly (rnafold3d/domains.py) against the whole-chain search, on library chains of 240-400 nt.

Targets: library chains of 240-400 nt (A/C/G/U only, >= 90 % of C4' atoms resolved), one per non-redundant
class, at most --n of them. Every target is predicted from its sequence alone with templates off (the point is
the search, not the template library), once with --domains auto and once with --domains off, with the same
seed and number of trajectories. Scored: TM-score of C4' atoms against the library coordinates (residue-index
correspondence, US-align's RNA d0, as template_selection.py), valid trajectories and wall time.

    python benchmark/long_rna.py run --n 8 --trajectories 16 --workers 16
    python benchmark/long_rna.py report
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'src'))
sys.path.insert(0, str(HERE))
RUNS = HERE / 'runs' / 'long_rna'


def targets(n):
    from rnafold3d.templates import load_library
    lib = load_library()
    out, seen = [], set()
    for k in range(lib.n_chains):
        seq = lib.chain_seq(k)
        if not 240 <= len(seq) <= 400 or set(seq) - set('ACGU') or lib.class_of[k] in seen:
            continue
        X = lib.chain_xyz(k)
        if np.isnan(X[:, 0]).mean() > 0.1:
            continue
        seen.add(lib.class_of[k])
        m = lib.meta[k]
        out.append({'id': f"{m['pdb']}_{m['chain']}", 'k': k, 'sequence': seq, 'title': m['title'][:60]})
    rng = np.random.default_rng(7)
    idx = sorted(rng.choice(len(out), min(n, len(out)), replace=False))
    return [out[i] for i in idx]


def c4_of(pdb):
    from rnafold3d.domains import read_model
    return np.array([r["C4'"] for r in read_model(pdb)])


def cmd_run(a):
    from rnafold3d.predict import predict, PredictConfig
    for t in targets(a.n):
        for mode in ('auto', 'off'):
            out = RUNS / mode / t['id']
            if (out / 'metadata' / 'manifest.json').exists():
                continue
            cfg = PredictConfig(n_models=5, trajectories=a.trajectories, workers=a.workers, seed=42, templates='off',
                                domains=mode)
            print(f"== {t['id']} {len(t['sequence'])} nt domains={mode}", flush=True)
            t0 = time.time()
            man = predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
            print('   ', man['statistics'], round(time.time() - t0, 1), 's', flush=True)


def cmd_report(a):
    from rnafold3d.templates import load_library
    from template_selection import tm_aligned
    lib = load_library()
    rows = []
    for t in targets(a.n):
        ref = lib.chain_xyz(t['k'])
        row = {'id': t['id'], 'n': len(t['sequence'])}
        for mode in ('auto', 'off'):
            run = RUNS / mode / t['id']
            mp = run / 'metadata' / 'manifest.json'
            if not mp.exists():
                continue
            man = json.loads(mp.read_text())
            tms = []
            for c in man['conformers']:
                X = c4_of(run / c['files']['pdb'])
                ok = ~np.isnan(ref[:, 0])
                tms.append(tm_aligned(X[ok], ref[ok], len(ref)))
            row[mode] = {'top1': tms[0] if tms else 0.0, 'best5': max(tms) if tms else 0.0,
                         'valid': man['statistics']['valid'], 'time': man['statistics']['wall_time_s'],
                         'domains': len((man['input'].get('domains') or {}).get('domains') or [])}
        rows.append(row)
        f = lambda r: f"{r['top1']:.3f}/{r['best5']:.3f} valid {r['valid']:2d} {r['time']:6.0f}s" if r else '-'
        print(f"{row['id']:8s} {row['n']:4d} nt  auto {f(row.get('auto'))} ({(row.get('auto') or {}).get('domains', 0)} domains)  off {f(row.get('off'))}")
    both = [r for r in rows if 'auto' in r and 'off' in r]
    if both:
        for key in ('top1', 'best5', 'valid', 'time'):
            print(f"mean {key}: domains {np.mean([r['auto'][key] for r in both]):.3f}  whole chain {np.mean([r['off'][key] for r in both]):.3f}")
    (RUNS / 'report.json').write_text(json.dumps(rows, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--n', type=int, default=8)
    r.add_argument('--trajectories', type=int, default=16)
    r.add_argument('--workers', type=int, default=16)
    p = sub.add_parser('report')
    p.add_argument('--n', type=int, default=8)
    a = ap.parse_args()
    {'run': cmd_run, 'report': cmd_report}[a.cmd](a)


if __name__ == '__main__':
    main()
