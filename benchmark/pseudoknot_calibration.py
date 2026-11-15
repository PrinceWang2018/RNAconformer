"""Calibration of the pseudoknot penalty (rnafold3d/pseudoknots.py) on the library structures.

For every library chain of 40-200 nt (A/C/G/U only, >= 5 annotated pairs; one chain per non-redundant class),
the nested ViennaRNA MFE and the structure with one crossing helix (fold_with_pseudoknot) are compared with the
canonical base pairs annotated in the chain's 3D structure, for several penalties: base-pair F1, PPV and
sensitivity summed over all chains, chains whose prediction changes, gains or loses correct pairs.

    python benchmark/pseudoknot_calibration.py --penalties 2,3,4,5,6,8
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))


def main():
    import RNA
    from rnafold3d.templates import load_library
    from rnafold3d.annotate import dotbracket_to_pairs
    from rnafold3d import pseudoknots as PK
    ap = argparse.ArgumentParser()
    ap.add_argument('--penalties', default='2,3,4,5,6,8')
    ap.add_argument('--out', default=str(ROOT / 'benchmark' / 'runs' / 'pseudoknot_calibration.json'))
    a = ap.parse_args()
    pens = [float(x) for x in a.penalties.split(',')]
    lib = load_library()
    tot = {p: [0, 0, 0, 0, 0, 0] for p in pens}          # tp, predicted, true, changed, better, worse
    base = [0, 0, 0]
    n = 0
    for k in range(lib.n_chains):
        seq = lib.chain_seq(k)
        if not 40 <= len(seq) <= 200 or set(seq) - set('ACGU'):
            continue
        true = set(lib.chain_pairs(k))
        if len(true) < 5:
            continue
        n += 1
        md = RNA.md()
        md.noLP = 1
        s, _e = RNA.fold_compound(seq, md).mfe()
        mfe = {(i - 1, j - 1) for i, j in dotbracket_to_pairs(s)}
        tp0 = len(mfe & true)
        base = [base[0] + tp0, base[1] + len(mfe), base[2] + len(true)]
        for p in pens:
            r = PK.fold_with_pseudoknot(seq, penalty=p)
            pred = mfe if r is None else set(r['nested']) | set(r['pk'])
            tp = len(pred & true)
            z = tot[p]
            z[0] += tp
            z[1] += len(pred)
            z[2] += len(true)
            if r is not None:
                z[3] += 1
                z[4] += tp > tp0
                z[5] += tp < tp0
    f1 = lambda tp, pr, tr: 2 * tp / (pr + tr)
    rows = {'chains': n, 'mfe': {'f1': f1(*base), 'ppv': base[0] / base[1], 'sensitivity': base[0] / base[2]}}
    print(f'{n} chains; nested MFE: F1 {f1(*base):.4f} PPV {base[0] / base[1]:.3f} sensitivity {base[0] / base[2]:.3f}')
    for p, z in tot.items():
        rows[str(p)] = {'f1': f1(*z[:3]), 'ppv': z[0] / z[1], 'sensitivity': z[0] / z[2], 'changed': z[3],
                        'better': z[4], 'worse': z[5]}
        print(f'penalty {p:4.1f}: F1 {f1(*z[:3]):.4f} PPV {z[0] / z[1]:.3f} sensitivity {z[0] / z[2]:.3f} '
              f'changed {z[3]} (more correct pairs {z[4]}, fewer {z[5]})')
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    main()
