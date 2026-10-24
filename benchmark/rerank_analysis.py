"""Compare top-1 selection rules offline on saved candidates (no re-sampling).

Rules: 'score' (engine total score), 'score_no_torsion' (total minus torsion-strain terms),
'medoid' (candidate with the lowest mean C4' RMSD to the other candidates of the same run),
'ss_medoid' (medoid among candidates of the lowest-energy 2D hypothesis), 'oracle_best' (ceiling)
and 'mean' (expected value of a random pick). Choose a rule on the dev split only.

Usage: PYTHONPATH=src python benchmark/rerank_analysis.py --run-id v2_0_pilot
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'src'))
from rnafold3d.metrics import read_numbered_pdb, c4_rmsd, kabsch_rmsd   # noqa: E402


def load(run_dir, sid):
    rows = {int(r['trajectory']): r for r in csv.DictReader((run_dir / 'metadata' / 'scores.csv').open())}
    val = json.loads((run_dir / 'metadata' / 'validation.json').read_text())
    ref = read_numbered_pdb(ROOT / 'data' / 'samples' / sid / 'reference.pdb')
    cands = []
    for p in sorted((run_dir / 'candidates').glob('traj_*.pdb')):
        k = int(p.stem.split('_')[1])
        r = rows[k]
        pred = read_numbered_pdb(p)
        loop_tor = sum(c.get('torsion_energy', 0.0) or 0.0 for c in (val[str(k)]['closures'] or []))
        cands.append({'k': k, 'ss': r['ss_id'], 'score': float(r['total_score']),
                      'torsion': float(r['torsion'] or 0) + loop_tor,
                      'rmsd': c4_rmsd(ref, pred)[0],
                      'c4': np.array([pred[i].atoms["C4'"] for i in sorted(pred)])})
    return cands


def pick(cands, rule):
    if rule == 'score':
        return min(cands, key=lambda c: c['score'])
    if rule == 'score_no_torsion':
        return min(cands, key=lambda c: c['score'] - c['torsion'])
    if rule in ('medoid', 'ss_medoid'):
        pool = cands
        if rule == 'ss_medoid':
            ss = sorted({c['ss'] for c in cands})[0]
            pool = [c for c in cands if c['ss'] == ss] or cands
        if len(pool) == 1:
            return pool[0]
        m = np.array([[kabsch_rmsd(a['c4'], b['c4']) for b in pool] for a in pool])
        return pool[int(np.argmin(m.sum(axis=1)))]
    raise ValueError(rule)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run-id', required=True)
    a = ap.parse_args()
    run = ROOT / 'runs' / a.run_id
    ds = {s['sample_id']: s for s in json.loads((ROOT / 'data' / 'dataset_v1.json').read_text())['samples']}
    rules = ('score', 'score_no_torsion', 'medoid', 'ss_medoid')
    table = []
    for mode in ('independent', 'oracle'):
        for tdir in sorted((run / mode).iterdir()):
            if not (tdir / 'candidates').exists():
                continue
            cands = [c for c in load(tdir, tdir.name) if c['rmsd'] is not None]
            if not cands:
                continue
            row = {'sample_id': tdir.name, 'split': ds[tdir.name]['split'], 'mode': mode, 'n': len(cands),
                   'oracle_best': min(c['rmsd'] for c in cands), 'mean': float(np.mean([c['rmsd'] for c in cands]))}
            for rule in rules:
                row[rule] = pick(cands, rule)['rmsd']
            table.append(row)
    out = {}
    for split in ('dev', 'test'):
        for mode in ('independent', 'oracle'):
            rs = [r for r in table if r['split'] == split and r['mode'] == mode]
            if rs:
                out[f'{split}/{mode}'] = {k: round(float(np.median([r[k] for r in rs])), 2)
                                          for k in ('oracle_best', 'mean') + rules}
                out[f'{split}/{mode}']['targets'] = len(rs)
    for r in table:
        print(json.dumps({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()}))
    print(json.dumps(out, indent=1))
    (run / 'reports').mkdir(exist_ok=True)
    (run / 'reports' / 'rerank_analysis.json').write_text(json.dumps({'per_target': table, 'medians': out}, indent=1))


if __name__ == '__main__':
    main()
