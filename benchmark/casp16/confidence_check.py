"""Blind check of the confidence model (rnafold3d.qa) on the CASP16 models.

The confidence model predicts the all-atom lDDT of every nucleotide from the coordinates of a model alone. Here
its prediction is compared with the lDDT the models really have against the experimental structures. No CASP16
target was used to train or to choose the model (training: tools/qa_model, 132 BGSU chains), and nothing in the
engine is tuned on this table.

    python benchmark/casp16/confidence_check.py --runs benchmark/casp16/runs_v076

Output: a table on stdout and benchmark/casp16/results/confidence_check.json. One core, a few minutes.
"""
from __future__ import annotations

import argparse
import difflib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import qa                                   # noqa: E402
from rnafold3d.structure_io import read_pdb                # noqa: E402

BANDS = ((0.0, 0.5, 'very low'), (0.5, 0.65, 'low'), (0.65, 0.8, 'medium'), (0.8, 1.01, 'high'))


class _Res:
    def __init__(self, atoms):
        self.atoms = atoms


def _chain(path):
    return max(read_pdb(path)[0].values(), key=lambda ch: len(ch.residues))


def reference(tid, seq):
    """{index in the target sequence: residue of the experimental structure} (resolved residues only)."""
    rr = list(_chain(HERE / 'cache' / 'references' / f'{tid}.pdb').residues)
    rseq = ''.join(r.base[-1] for r in rr)
    out = {}
    for a, b, n in difflib.SequenceMatcher(None, rseq, seq, autojunk=False).get_matching_blocks():
        for k in range(n):
            out[b + k] = _Res({nm: np.asarray(x, float) for nm, x in rr[a + k].atoms.items()})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', default=str(HERE / 'runs_v076'))
    ap.add_argument('--out', default=str(HERE / 'results' / 'confidence_check.json'))
    a = ap.parse_args()
    runs = Path(a.runs)
    rows, P, T, TPL = [], [], [], []
    for t in json.loads((HERE / 'targets.json').read_text()):
        mp = runs / t['id'] / 'metadata' / 'manifest.json'
        if not mp.exists():
            continue
        man, seq = json.loads(mp.read_text()), t['sequence']
        ref = reference(t['id'], seq)
        per = []
        for c in man['conformers']:
            mr = list(_chain(runs / t['id'] / c['files']['pdb']).residues)
            p = qa.predict(seq, mr)
            if p is None:
                continue
            true = qa.lddt_per_residue(ref, {i: _Res({nm: np.asarray(x, float) for nm, x in r.atoms.items()})
                                             for i, r in enumerate(mr)})
            idx = sorted(true)
            pv, tv = np.array([p['per_residue'][i] for i in idx]), np.array([true[i] for i in idx])
            P.append(pv), T.append(tv), TPL.append(np.full(len(idx), bool(c.get('template'))))
            per.append({'rank': c['rank'], 'confidence': round(float(pv.mean()), 4), 'lddt': round(float(tv.mean()), 4),
                        'r_per_nucleotide': round(float(np.corrcoef(pv, tv)[0, 1]), 3) if pv.std() > 0 and tv.std() > 0 else None,
                        'template': c.get('template') or None})
        if not per:
            continue
        rows.append({'target': t['id'], 'length': len(seq), 'resolved': len(ref), 'models': per})
        print(f"{t['id']:7s} {len(seq):4d} nt  model 1: confidence {per[0]['confidence']:.2f}, lDDT {per[0]['lddt']:.2f} "
              f"({'template' if per[0]['template'] else 'no template'}) | confidence " + ' '.join(f"{m['confidence']:.2f}" for m in per)
              + ' | lDDT ' + ' '.join(f"{m['lddt']:.2f}" for m in per), flush=True)
    P, T, TPL = np.concatenate(P), np.concatenate(T), np.concatenate(TPL)
    g = np.array([[m['confidence'], m['lddt']] for r in rows for m in r['models']])
    g1 = np.array([[r['models'][0]['confidence'], r['models'][0]['lddt']] for r in rows])
    summary = {
        'runs': runs.name, 'targets': len(rows), 'models': len(g), 'nucleotides': int(len(P)),
        'model_1_r_over_targets': round(float(np.corrcoef(g1[:, 0], g1[:, 1])[0, 1]), 3),
        'per_model_r': round(float(np.corrcoef(g[:, 0], g[:, 1])[0, 1]), 3),
        'per_model_bias': round(float(np.mean(g[:, 0] - g[:, 1])), 3),
        'per_nucleotide_r': round(float(np.corrcoef(P, T)[0, 1]), 3),
        'per_nucleotide_mean_abs_error': round(float(np.abs(P - T).mean()), 3),
        'lddt_model_1': round(float(g1[:, 1].mean()), 3),
        'lddt_highest_confidence': round(float(np.mean([max(r['models'], key=lambda m: m['confidence'])['lddt'] for r in rows])), 3),
        'lddt_best_of_5': round(float(np.mean([max(m['lddt'] for m in r['models']) for r in rows])), 3),
        'bands': []}
    print(f"\n{summary['targets']} targets, {summary['models']} models, {summary['nucleotides']} nucleotides ({runs.name})")
    print(f"model 1, over the targets: r = {summary['model_1_r_over_targets']:.3f}; all models: r = {summary['per_model_r']:.3f}, "
          f"bias (confidence - lDDT) {summary['per_model_bias']:+.3f}")
    print(f"per nucleotide: r = {summary['per_nucleotide_r']:.3f}, mean absolute error {summary['per_nucleotide_mean_abs_error']:.3f}")
    print(f"mean lDDT of model 1 {summary['lddt_model_1']:.3f}; of the model with the highest confidence "
          f"{summary['lddt_highest_confidence']:.3f}; of the best of 5 {summary['lddt_best_of_5']:.3f}")
    print(f"{'band':22s} {'nucleotides':>11s} {'true lDDT':>10s} {'from a template':>26s} {'without template':>26s}")
    for lo, hi, name in BANDS:
        s = (P >= lo) & (P < hi)
        row = {'band': name, 'from': lo, 'to': min(hi, 1.0), 'nucleotides': int(s.sum()),
               'lddt': round(float(T[s].mean()), 3) if s.any() else None}
        cells = []
        for key, sel in (('template', TPL), ('no_template', ~TPL)):
            k = s & sel
            row[key] = {'nucleotides': int(k.sum()), 'lddt': round(float(T[k].mean()), 3) if k.any() else None}
            cells.append(f"{int(k.sum()):6d} nt, lDDT {T[k].mean():.3f}" if k.any() else f"{0:6d} nt")
        summary['bands'].append(row)
        print(f"{name + f' ({lo:.2f}-{min(hi, 1.0):.2f})':22s} {row['nucleotides']:11d} {row['lddt'] if row['lddt'] is not None else float('nan'):10.3f} "
              f"{cells[0]:>26s} {cells[1]:>26s}")
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({'summary': summary, 'targets': rows}, indent=1))
    print(f'written {a.out}')


if __name__ == '__main__':
    main()
