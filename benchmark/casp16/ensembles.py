"""Do the five models of a method describe an RNA that moves? CASP16 targets with several experimental states.

Three CASP16 RNAs were solved in more than one conformation:
  R1203  HIV-1 RRE stem-loop II, two crystal structures (TM-score 0.56 between them),
  R1255  SARS-CoV-2 SL5 (124 nt), 10 cryo-EM conformations of a rotating four-way junction,
  R1256  bat coronavirus HKU5 SL5 (127 nt), 30 cryo-EM conformations.
For every method (RCV and the CASP16 groups of compare.py) the five submitted models are scored against every
experimental conformation (US-align, as in compare.py). Reported per target:
  coverage   mean over the experimental conformations of the best TM-score among the five models,
  states     how many experimental conformations are matched best by different models (1 = one model fits all),
  spread     1 - mean pairwise TM-score of the five models, next to the spread of the experimental ensemble.
And for all 15 single-state targets of compare.py: the spread of the five models of each method.

    python benchmark/casp16/ensembles.py run --workers 16      # RCV on R1255 and R1256 (R1203 comes from compare.py run)
    python benchmark/casp16/ensembles.py evaluate              # -> results/ensembles.json, results/ENSEMBLES.md
"""
from __future__ import annotations

import argparse
import itertools
import json
import shutil
import statistics as st
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / 'src'))
import compare  # noqa: E402

MULTI = [
    {'id': 'R1203', 'family': 'HIV-1 RRE stem-loop II', 'refs': ['R1203v1', 'R1203v2'],
     'sequence': next(t['sequence'] for t in compare.TARGETS if t['id'] == 'R1203')},
    {'id': 'R1255', 'family': 'SARS-CoV-2 SL5', 'refs': [f'R1255_{i}' for i in range(1, 11)],
     'sequence': 'GGACACGAGUAACUCGUCUAUCUGCUGCAGGCUGCUUACGGUUUCGUCCGUGUUGCAGCCGAUCAUCAGAACAUCUAGGUUUCGUCCGGGUGU'
                 'UACCGAAAGGUCAGAUGGAGAGCCUUGUCCC'},
    {'id': 'R1256', 'family': 'Bat coronavirus HKU5 SL5', 'refs': [f'R1256_{i}' for i in range(1, 31)],
     'sequence': 'GGAGCAUCGUGUCUCAAGUGCUUCACGGUCACAAUAUACCGUUUCGUCGGGUGCGUGGCAAUUCGGUGCACAUCAUGUCUUUCGUGGCUGGUG'
                 'UGGUCAAGGUGCGACAAGUAUAGAGCAGAGCUCC'},
]
METHODS = {'rcv': 'RCV 0.7.3', '304': 'AlphaFold 3', '052': 'trRosettaRNA', '481': 'Vfold', '183': 'GuangzhouRNA'}


def references(c, t):
    """Experimental conformations of a target as cleaned PDB files (RNA atoms, chain A)."""
    out = c / 'references_multi' / t['id']
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(c / 'casp16.targets_RNA.tgz') as tf:
        for r in t['refs']:
            dest = out / f'{r}.pdb'
            if not dest.exists():
                raw = tf.extractfile(f'./Targets_RNA/{r}.pdb').read().decode()
                chain = next(l[21] for l in raw.splitlines() if l.startswith('ATOM'))
                dest.write_text(compare.clean_reference(raw, chain))
    return [out / f'{r}.pdb' for r in t['refs']]


def group_models(c, t, g, tmp):
    """The five models a CASP16 group submitted (paths of temporary files), [] if missing."""
    tgz = c / 'predictions' / f"{t['id']}.tar.gz"
    out = []
    with tarfile.open(tgz) as tf:
        for k in range(1, 6):
            try:
                data = tf.extractfile(f"{t['id']}/{t['id']}TS{g}_{k}").read()
            except KeyError:
                continue
            p = tmp / f"{t['id']}_{g}_{k}.pdb"
            p.write_bytes(data)
            out.append(p)
    return out


def rcv_models(run, tid):
    man_p = run / tid / 'metadata' / 'manifest.json'
    if not man_p.exists():
        return [], None
    man = json.loads(man_p.read_text())
    return [run / tid / c['files']['pdb'] for c in man['conformers']], man


def spread(models):
    """1 - mean pairwise TM-score (residue-index superposition) of a set of models of one sequence."""
    tms = [compare.score(a, b)[0] for a, b in itertools.combinations(models, 2)]
    return round(1 - st.mean(tms), 3) if tms else None


def cmd_run(a):
    from rnafold3d.predict import predict, PredictConfig
    for t in MULTI:
        if t['id'] == 'R1203':
            continue
        out = Path(a.runs) / t['id']
        if (out / 'metadata' / 'manifest.json').exists() and not a.force:
            print('skip (done)', t['id'])
            continue
        if out.exists():
            shutil.rmtree(out)
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=a.seed, workers=a.workers,
                            templates=a.templates, template_cutoff=a.template_cutoff)
        print(f"== {t['id']} ({len(t['sequence'])} nt)", flush=True)
        man = predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
        print('   ', man['statistics'], flush=True)


def cmd_evaluate(a):
    c = Path(a.cache)
    tmp = c / 'ensemble_tmp'
    tmp.mkdir(exist_ok=True)
    res = {'multi': [], 'single_state_spread': {}}
    for t in MULTI:
        refs = references(c, t)
        ref_spread = spread(refs)
        row = {'target': t['id'], 'family': t['family'], 'length': len(t['sequence']), 'conformations': len(refs),
               'experimental_spread': ref_spread, 'methods': {}}
        for m in METHODS:
            if m == 'rcv':
                models, man = rcv_models(Path(a.runs), t['id'])
            else:
                models = group_models(c, t, m, tmp)
            if not models:
                continue
            mat = [[compare.score(p, r)[0] for r in refs] for p in models]       # models x conformations
            best = [max(mat[k][j] for k in range(len(models))) for j in range(len(refs))]
            arg = {max(range(len(models)), key=lambda k: mat[k][j]) for j in range(len(refs))}
            if m == 'rcv' and man:
                released = compare._release_dates()
                used = sorted({c.get('template') for c in man['conformers'] if c.get('template')})
                row['rcv_valid'] = f"{man['statistics']['valid']} of {man['statistics']['trajectories']}"
                row['rcv_templates'] = [f"{u} ({' '.join(released.get(p.split('_')[0].upper(), '?') for p in u.split('+'))})"
                                        for u in used]
            row['methods'][m] = {'coverage': round(st.mean(best), 3), 'worst_state': round(min(best), 3),
                                 'model1_best_state': round(max(mat[0]), 3), 'states_by_different_models': len(arg),
                                 'spread': spread(models), 'n_models': len(models),
                                 'tm_matrix': [[round(x, 3) for x in r] for r in mat]}
            print(f"{t['id']} {METHODS[m]:14s} coverage {row['methods'][m]['coverage']:.3f}  worst "
                  f"{row['methods'][m]['worst_state']:.3f}  spread {row['methods'][m]['spread']}  "
                  f"(experimental {ref_spread})  models used {len(arg)}", flush=True)
        res['multi'].append(row)
    # spread of the five models on the 15 single-state targets
    for m in METHODS:
        vals = {}
        for t in compare.TARGETS:
            models = rcv_models(Path(a.runs), t['id'])[0] if m == 'rcv' else group_models(c, t, m, tmp)
            if len(models) >= 2:
                vals[t['id']] = spread(models)
        res['single_state_spread'][m] = {'per_target': vals, 'median': round(st.median(vals.values()), 3) if vals else None}
        print(f"spread of 5 models, 15 targets, {METHODS[m]}: median {res['single_state_spread'][m]['median']}", flush=True)
    shutil.rmtree(tmp, ignore_errors=True)
    out = HERE / 'results'
    (out / 'ensembles.json').write_text(json.dumps(res, indent=1))
    (out / 'ENSEMBLES.md').write_text(report(res))


def report(res):
    L = ['# RNAs that move: CASP16 targets with several experimental conformations', '',
         'Generated by `benchmark/casp16/ensembles.py evaluate`. Each method\'s five models are scored against every '
         'experimental conformation (US-align TM-score, as in `compare.py`). *Coverage*: mean over the conformations of '
         'the best TM-score among the five models. *Worst state*: the conformation fitted least well. *Spread*: '
         '1 − mean pairwise TM-score of the five models (0 = five copies of one structure); the experimental '
         'ensemble\'s own spread is given for comparison. *Models used*: how many different models are the best fit '
         'for at least one conformation.', '']
    for row in res['multi']:
        tp = row.get('rcv_templates')
        L += [f"## {row['target']} — {row['family']}, {row['length']} nt, {row['conformations']} experimental "
              f"conformations (spread {row['experimental_spread']:.2f})", '',
              (f"RCV's templates: {', '.join(tp)}." if tp else "RCV found no template.")
              + (f" Valid RCV trajectories: {row['rcv_valid']}." if row.get('rcv_valid') else ''), '',
              '| Method | Coverage | Worst state | Model 1, best state | Spread of the 5 models | Models used |',
              '|---|---|---|---|---|---|']
        for m, v in row['methods'].items():
            L.append(f"| {METHODS[m]} | {v['coverage']:.3f} | {v['worst_state']:.3f} | {v['model1_best_state']:.3f} | "
                     f"{v['spread']:.2f} | {v['states_by_different_models']} of {v['n_models']} |")
        L.append('')
    L += ['## Spread of the five models on the 15 single-state targets (median)', '',
          '| Method | Median spread |', '|---|---|']
    for m, v in res['single_state_spread'].items():
        L.append(f"| {METHODS[m]} | {v['median']:.2f} |")
    L += ['', '## Notes', '',
          '* RCV ran with the CASP16 settings of `compare.py` (5 models, 40 trajectories, templates released before '
          '2024-05-01). For R1255 and R1256 the library holds earlier cryo-EM structures of the same RNAs (PDB 8UYS, '
          '8UYJ, 8UYM and 8QO5, released December 2023 to March 2024, 4.7–7.3 Å), which RCV used as templates (listed '
          'per target above); the CASP16 groups could use them too. Few valid trajectories (R1256) mostly mean failed loop closures on these coarse templates.',
          '* Spread is only useful when the alternatives are real: on a rigid RNA, five different models mean four '
          'wrong ones. RCV builds its five models from different secondary-structure and template hypotheses on '
          'purpose; the CASP16 groups were asked for five ranked guesses of one structure.']
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['run', 'evaluate'])
    ap.add_argument('--cache', default=str(HERE / 'cache'))
    ap.add_argument('--runs', default=str(HERE / 'runs_v073_tmpl'))
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--trajectories', type=int, default=40)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--templates', choices=['auto', 'off'], default='auto')
    ap.add_argument('--template-cutoff', default='2024-05-01')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    {'run': cmd_run, 'evaluate': cmd_evaluate}[a.command](a)


if __name__ == '__main__':
    main()
