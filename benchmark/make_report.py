"""Write reports/benchmark_summary.md for a benchmark run (tables per split and mode).

Usage: python benchmark/make_report.py --run-id v2_0_pilot
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def fnum(x, nd=1):
    if x in (None, ''):
        return '—'
    return f'{float(x):.{nd}f}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    a = ap.parse_args()
    run = ROOT / 'runs' / a.run_id
    rows = list(csv.DictReader((run / 'reports' / 'per_target.csv').open()))
    null = json.loads((ROOT / 'runs' / 'null_baseline.json').read_text()) if (ROOT / 'runs' / 'null_baseline.json').exists() else {}
    agg = json.loads((run / 'reports' / 'aggregate.json').read_text())
    ds = {s['sample_id']: s for s in json.loads((ROOT / 'data' / 'dataset_v1.json').read_text())['samples']}
    lines = [f'# Benchmark summary — run `{a.run_id}`', '',
             f"Engine {agg['config']['engine']} {agg['config']['engine_version']}; {agg['config']['trajectories']} "
             f"trajectories and up to {agg['config']['models']} models per target and mode; seed {agg['config']['seed']}; "
             f"ViennaRNA {agg['config']['vienna']}.", '',
             'C4\' RMSD in Å over residues observed in the reference; **top-1** = rank-1 model, **best** = best of the '
             'exported models, **cand** = best of all valid trajectory models (sampling ceiling); **null** = median '
             'C4\' RMSD of 30 random arrangements with the correct (oracle) stems; INF-WC on canonical pairs; '
             'ViennaRNA F1 = MFE base-pair F1 vs. the reference nested pairs. A dash means no valid model.', '']
    for split in ('dev', 'test'):
        for mode in ('independent', 'oracle'):
            rs = [r for r in rows if r['split'] == split and r['mode'] == mode]
            if not rs:
                continue
            title = {'independent': 'sequence only (ViennaRNA 2D)', 'oracle': 'oracle 2D (from reference; upper bound)'}[mode]
            lines += [f'## {split} split — {title}', '',
                      '| target | nt | family | valid traj | top-1 | best | cand | null | eRMSD best | INF top-1 | Vienna F1 | time s |',
                      '|---|---|---|---|---|---|---|---|---|---|---|---|']
            for r in sorted(rs, key=lambda r: int(r['length'])):
                sid = r['sample_id']
                lines.append(f"| {sid} | {r['length']} | {r['family']} | {r['valid_trajectories']}/{r['trajectories']} | "
                             f"{fnum(r['top1_c4_rmsd'])} | {fnum(r['best_c4_rmsd'])} | {fnum(r.get('best_candidate_c4_rmsd'))} | "
                             f"{fnum(null.get(sid, {}).get('median'))} | {fnum(r['best_ermsd'], 2)} | {fnum(r['top1_inf_wc'], 2)} | "
                             f"{fnum(ds[sid]['independent_2d']['f1_vs_oracle_nested'], 2)} | {fnum(r['wall_time_s'], 0)} |")
            have = [r for r in rs if r['top1_c4_rmsd'] not in ('', None)]
            top = [float(r['top1_c4_rmsd']) for r in have]
            best = [float(r['best_c4_rmsd']) for r in have]
            nulls = [null[r['sample_id']]['median'] for r in have if r['sample_id'] in null]
            beat = sum(float(r['best_c4_rmsd']) < null[r['sample_id']]['median'] for r in have if r['sample_id'] in null)
            beat1 = sum(float(r['top1_c4_rmsd']) < null[r['sample_id']]['median'] for r in have if r['sample_id'] in null)
            lines += ['', f"Targets with ≥1 valid model: {len(have)}/{len(rs)}. "
                          f"Median top-1 {np.median(top):.1f} Å, median best-of-K {np.median(best):.1f} Å "
                          f"(null median {np.median(nulls):.1f} Å over the same targets). "
                          f"Top-1 better than null: {beat1}/{len(have)}; best-of-K better than null: {beat}/{len(have)}. "
                          f"Best-of-K ≤ 5 Å: {sum(b <= 5 for b in best)}/{len(rs)}; ≤ 10 Å: {sum(b <= 10 for b in best)}/{len(rs)}.", '']
    lines += ['## Aggregate (all targets; failures counted as failures)', '', '```json', json.dumps(agg['aggregate'], indent=1), '```', '']
    (run / 'reports' / 'benchmark_summary.md').write_text('\n'.join(lines))
    print('\n'.join(lines[:60]))


if __name__ == '__main__':
    main()
