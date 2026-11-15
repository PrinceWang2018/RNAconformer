"""Fit the weight of the base-base potential (rnafold3d/basepot.py) in the ranking score, offline.

Uses a development-benchmark run made with --candidates (every valid trajectory model saved, junction
geometry on as by default, so the weight is fitted with those restraints present). For every target and weight w,
model 1 is the candidate with the lowest total_score + w * basepot; the reported number is the mean TM-score of
model 1 over the targets (and of the candidate ranked first among those of template hypotheses, which the engine
always ranks first). Template targets keep their template ranking.

    python benchmark/rerank_basepot.py --label base075
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'src'))
sys.path.insert(0, str(HERE))


def read_pdb(p):
    res = {}
    for line in Path(p).read_text().splitlines():
        if line.startswith('ATOM'):
            k = int(line[22:26])
            res.setdefault(k, {})[line[12:16].strip()] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return [res[k] for k in sorted(res)]


def main():
    import dev_templates as B
    from rnafold3d import basepot as BP
    from rnafold3d.annotate import dotbracket_to_pairs
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', default='base075')
    ap.add_argument('--weights', default='0,0.05,0.1,0.2,0.3,0.5,0.75,1,1.5,2,3')
    a = ap.parse_args()
    cache = B.RUNS / a.label / 'rerank_basepot.json'
    data = json.loads(cache.read_text()) if cache.exists() else {}
    for t in B.targets():
        run = B.RUNS / a.label / t['id']
        if t['id'] in data or not (run / 'candidates').is_dir():
            continue
        rows = {int(r['trajectory']): r for r in csv.DictReader((run / 'metadata' / 'scores.csv').open())}
        hyps = {}
        lines = (run / 'input' / 'secondary_structures.dbn').read_text().splitlines()
        for k in range(0, len(lines), 3):
            hid = lines[k].split('|')[1]
            hyps[hid] = (lines[k], lines[k + 2])
        cands = []
        for p in sorted((run / 'candidates').glob('traj_*.pdb')):
            k = int(p.stem.split('_')[1])
            r = rows[k]
            head, dbn = hyps[r['ss_id']]
            nested = ''.join(c if c in '().' else '.' for c in dbn)
            pairs = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(nested)]
            e = BP.score(t['sequence'], read_pdb(p), pairs)
            tm = B.tm(p, t['reference'])[0]
            cands.append({'k': k, 'ss': r['ss_id'], 'template': 'template' in head, 'total': float(r['total_score']),
                          'bp': e, 'tm': tm})
        data[t['id']] = cands
        print(t['id'], len(cands), flush=True)
    cache.write_text(json.dumps(data))
    ws = [float(x) for x in a.weights.split(',')]
    print(f"{'w':>5s} {'top1 all':>8s} {'no-tmpl':>8s} {'tmpl':>6s}  better/worse vs w=0 (no-template targets)")
    base = {}
    for w in ws:
        top, top_nt, top_t, diff = [], [], [], []
        for tid, cs in data.items():
            if not cs:
                continue
            tm_c = [c for c in cs if c['template']]
            pool = tm_c if tm_c else cs
            best = min(pool, key=lambda c: c['total'] + w * c['bp'])
            top.append(best['tm'])
            (top_t if tm_c else top_nt).append(best['tm'])
            if w == 0:
                base[tid] = best['tm']
            elif not tm_c:
                diff.append(best['tm'] - base[tid])
        print(f'{w:5.2f} {np.mean(top):8.4f} {np.mean(top_nt):8.4f} {np.mean(top_t):6.4f}  '
              f'{sum(d > 0.005 for d in diff)}/{sum(d < -0.005 for d in diff)}')
    orc = np.mean([max(c['tm'] for c in cs) for cs in data.values() if cs])
    rnd = np.mean([np.mean([c['tm'] for c in cs]) for cs in data.values() if cs])
    print(f'oracle (best candidate) {orc:.4f}; random candidate {rnd:.4f}')
    # rank correlation of the potential with TM-score within targets
    rs = []
    from scipy.stats import spearmanr
    for cs in data.values():
        if len(cs) >= 5:
            rs.append(spearmanr([c['bp'] for c in cs], [-c['tm'] for c in cs])[0])
    print('mean within-target Spearman(basepot, -TM):', round(float(np.nanmean(rs)), 3), 'n =', len(rs))
    rs = []
    for cs in data.values():
        if len(cs) >= 5:
            rs.append(spearmanr([c['total'] for c in cs], [-c['tm'] for c in cs])[0])
    print('mean within-target Spearman(total score, -TM):', round(float(np.nanmean(rs)), 3))


if __name__ == '__main__':
    main()
