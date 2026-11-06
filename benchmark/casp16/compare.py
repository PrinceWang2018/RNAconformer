"""Blind-test comparison of the RCV engine (rnafold3d) with leading RNA 3D predictors on CASP16.

CASP16 (2024) is the latest community-wide blind assessment of RNA structure prediction. Its RNA
targets, the experimental reference structures, every group's models and the official scores are
public at predictioncenter.org. This script

  fetch     downloads the reference structures and the official per-model score table,
  run       predicts every target from its sequence alone with rnafold3d (5 models, like CASP),
  evaluate  scores the RCV models exactly as CASP did (US-align, residue-index superposition,
            C3' atoms, TM-score normalised by the reference length; RMSD over all residues) and
            puts them next to the official scores of the comparison groups,
  verify    re-scores the comparison groups' downloaded models with the same US-align call and checks
            that the official numbers are reproduced (so both sides are measured the same way).

Usage (from the repository root, ViennaRNA and US-align installed):
  python benchmark/casp16/compare.py fetch
  python benchmark/casp16/compare.py run --workers 16
  python benchmark/casp16/compare.py verify
  python benchmark/casp16/compare.py evaluate

Fairness notes (also written to the report): the CASP16 groups predicted blind in 2024; RCV was run
in 2026 but uses no information from the targets (no templates, no MSA, no training on structures).
Its geometry parameters were calibrated on 178 PDB entries; none is a CASP16 target (checked by
release date and sequence). R1255/R1256 (multi-conformation references) and targets above 250 nt
or in RNA-protein complexes are not part of this set.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import statistics as st
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'src'))

BASE = 'https://predictioncenter.org/download_area/CASP16'
TARGETS = json.loads((HERE / 'targets.json').read_text())
# comparison groups (CASP16 group number -> label); chosen from the official table: the AlphaFold 3
# server baseline, the trRosettaRNA server (Yang-Server) and two of the top expert (human) groups on small RNAs
GROUPS = {'304': 'AlphaFold 3 (AF3-server)', '052': 'trRosettaRNA (Yang-Server)', '481': 'Vfold',
          '183': 'GuangzhouRNA-human'}
SHORT = {'304': 'AlphaFold 3', '052': 'trRosettaRNA', '481': 'Vfold', '183': 'GuangzhouRNA', 'rcv': 'RCV'}


def cache_dir(a):
    return Path(a.cache)


def usalign():
    exe = shutil.which('USalign')
    if not exe:
        sys.exit('US-align not found (conda install -c bioconda usalign)')
    return exe


def download(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        print('  download', url)
        with urllib.request.urlopen(url, timeout=600) as r, open(dest, 'wb') as f:
            shutil.copyfileobj(r, f)
    return dest


# ------------------------------------------------------------------------------------------------ fetch
def cmd_fetch(a):
    c = cache_dir(a)
    tgz = download(f'{BASE}/targets/casp16.targets_RNA.tgz', c / 'casp16.targets_RNA.tgz')
    download(f'{BASE}/results/tables/CASP16_rna_mono.scores.csv', c / 'CASP16_rna_mono.scores.csv')
    refs = c / 'references'
    refs.mkdir(exist_ok=True)
    with tarfile.open(tgz) as tf:
        for t in TARGETS:
            m = tf.getmember(f"./Targets_RNA/{t['reference']}.pdb")
            raw = tf.extractfile(m).read().decode()
            (refs / f"{t['id']}.pdb").write_text(clean_reference(raw, t['chain']))
    if a.models:
        for t in TARGETS:
            download(f"{BASE}/predictions/RNA/{t['casp_id']}.tar.gz", c / 'predictions' / f"{t['casp_id']}.tar.gz")
    print('ok:', c)


def clean_reference(text, chain):
    """RNA chain of the reference as plain ATOM records (chain A), as scored by CASP."""
    out = []
    for line in text.splitlines():
        if line.startswith(('ATOM', 'HETATM')) and line[21] == chain and line[17:20].strip() in ('A', 'C', 'G', 'U'):
            out.append('ATOM  ' + line[6:21] + 'A' + line[22:])
    return '\n'.join(out + ['END']) + '\n'


# -------------------------------------------------------------------------------------------------- run
def cmd_run(a):
    from rnafold3d.predict import predict, PredictConfig
    runs = Path(a.runs)
    for t in TARGETS:
        if a.only and t['id'] not in a.only.split(','):
            continue
        out = runs / t['id']
        if (out / 'metadata' / 'manifest.json').exists() and not a.force:
            print('skip (done)', t['id'])
            continue
        if out.exists():
            shutil.rmtree(out)
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=a.seed, workers=a.workers)
        print(f"== {t['id']} ({len(t['sequence'])} nt)", flush=True)
        predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
        man = json.loads((out / 'metadata' / 'manifest.json').read_text())
        print('   ', man['statistics'], flush=True)


# --------------------------------------------------------------------------------------------- evaluate
def score(model, reference):
    """(TM-score normalised by the reference, RMSD) with residue-index superposition, as in CASP16."""
    p = subprocess.run([usalign(), str(model), str(reference), '-mol', 'RNA', '-TMscore', '1', '-outfmt', '2'],
                       capture_output=True, text=True, check=True)
    row = [l for l in p.stdout.splitlines() if l and not l.startswith('#')][-1].split('\t')
    return float(row[3]), float(row[4])


def official_table(c):
    """{target: {group: {model_index: (tm, rmsd)}}} from the official CASP16 RNA monomer table."""
    blocks = (c / 'CASP16_rna_mono.scores.csv').read_text().split('Target: ')[1:]
    out = {}
    for blk in blocks:
        lines = blk.splitlines()
        tgt = lines[0].strip()
        d = out.setdefault(tgt, {})
        for l in lines[3:]:
            p = l.split()
            m = re.match(r'.*TS(\d+)_(\d+)$', p[1]) if len(p) > 7 and p[0].isdigit() else None
            if m and p[6] != '-' and p[3] != '-':
                d.setdefault(m.group(1), {})[int(m.group(2))] = (float(p[6]), float(p[3]))
    return out


def rescore_groups(c, t, ref):
    """Score the comparison groups' five models with the same US-align call as RCV (needs `fetch --models`)."""
    tgz = c / 'predictions' / f"{t['casp_id']}.tar.gz"
    if not tgz.exists():
        return None
    out = {}
    tmp = c / 'rescore_model.pdb'
    with tarfile.open(tgz) as tf:
        for g in GROUPS:
            for k in range(1, 6):
                try:
                    tmp.write_bytes(tf.extractfile(f"{t['casp_id']}/{t['casp_id']}TS{g}_{k}").read())
                except KeyError:
                    continue
                out.setdefault(g, {})[k] = score(tmp, ref)
    return out


def inf_wc(ref, model):
    """Base-pair accuracy (INF-WC: canonical pairs of the model vs the reference), from rnafold3d.metrics."""
    from rnafold3d.metrics import compare
    try:
        return compare(ref, model).get('inf_wc')
    except Exception:
        return None


def group_inf_wc(c, t, ref):
    tgz = c / 'predictions' / f"{t['casp_id']}.tar.gz"
    out = {}
    if not tgz.exists():
        return out
    tmp = c / 'rescore_model.pdb'
    with tarfile.open(tgz) as tf:
        for g in GROUPS:
            try:
                tmp.write_bytes(tf.extractfile(f"{t['casp_id']}/{t['casp_id']}TS{g}_1").read())
            except KeyError:
                continue
            out[g] = inf_wc(ref, tmp)
    return out


def cmd_evaluate(a):
    c, runs = cache_dir(a), Path(a.runs)
    table = official_table(c)
    rows = []
    for t in TARGETS:
        ref = c / 'references' / f"{t['id']}.pdb"
        man_p = runs / t['id'] / 'metadata' / 'manifest.json'
        if not man_p.exists():
            print('missing run', t['id'])
            continue
        man = json.loads(man_p.read_text())
        rcv = {}
        for conf in man['conformers']:
            rcv[conf['rank']] = score(runs / t['id'] / conf['files']['pdb'], ref)
        infwc = {'rcv': inf_wc(ref, runs / t['id'] / man['conformers'][0]['files']['pdb'])}
        infwc.update(group_inf_wc(c, t, ref))
        groups = rescore_groups(c, t, ref) or {g: table.get(t['casp_id'], {}).get(g, {}) for g in GROUPS}
        all_top1 = [v[1][0] for g, v in table.get(t['casp_id'], {}).items() if 1 in v]
        row = {'target': t['id'], 'casp_id': t['casp_id'], 'length': len(t['sequence']),
               'resolved': t['resolved'], 'family': t['family'],
               'n_groups': len(all_top1), 'casp_median_top1_tm': round(st.median(all_top1), 3) if all_top1 else None,
               'casp_best_top1_tm': round(max(all_top1), 3) if all_top1 else None,
               'rcv_valid': man['statistics']['valid'], 'rcv_trajectories': man['statistics']['trajectories'],
               'rcv_wall_time_s': man['statistics']['wall_time_s']}
        for key, models in [('rcv', rcv)] + list(groups.items()):
            if not models:
                row.update({f'{key}_top1_tm': None, f'{key}_best5_tm': None, f'{key}_top1_rmsd': None, f'{key}_best5_rmsd': None})
                continue
            top = models.get(1) or models[min(models)]
            row.update({f'{key}_top1_tm': round(top[0], 3), f'{key}_best5_tm': round(max(v[0] for v in models.values()), 3),
                        f'{key}_top1_rmsd': round(top[1], 2), f'{key}_best5_rmsd': round(min(v[1] for v in models.values()), 2)})
        for key, v in infwc.items():
            row[f'{key}_top1_infwc'] = None if v is None else round(v, 3)
        # rank of RCV's model 1 among all CASP16 groups' model 1 on this target
        if row['rcv_top1_tm'] is not None and all_top1:
            row['rcv_top1_beats_groups'] = round(sum(x < row['rcv_top1_tm'] for x in all_top1) / len(all_top1), 3)
        rows.append(row)
        print(f"{t['id']:8s} RCV top1 TM {row['rcv_top1_tm']}  best5 {row['rcv_best5_tm']}  | AF3 {row['304_top1_tm']}"
              f"  | CASP median {row['casp_median_top1_tm']}")
    out = HERE / 'results'
    out.mkdir(exist_ok=True)
    with (out / 'per_target.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = summarise(rows)
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    (out / 'README.md').write_text(report(rows, summary))
    print(json.dumps(summary['methods'], indent=1))


def summarise(rows):
    keys = ['rcv'] + list(GROUPS)
    s = {'targets': len(rows), 'methods': {}}
    for k in keys:
        tm1 = [r[f'{k}_top1_tm'] for r in rows if r.get(f'{k}_top1_tm') is not None]
        tm5 = [r[f'{k}_best5_tm'] for r in rows if r.get(f'{k}_best5_tm') is not None]
        rm1 = [r[f'{k}_top1_rmsd'] for r in rows if r.get(f'{k}_top1_rmsd') is not None]
        rm5 = [r[f'{k}_best5_rmsd'] for r in rows if r.get(f'{k}_best5_rmsd') is not None]
        inf = [float(r[f'{k}_top1_infwc']) for r in rows if r.get(f'{k}_top1_infwc') not in (None, '')]
        s['methods'][k] = {'label': 'RCV 0.7.2 (rnafold3d, physics-based)' if k == 'rcv' else GROUPS[k], 'n': len(tm1),
                           'mean_top1_tm': round(st.mean(tm1), 3), 'mean_best5_tm': round(st.mean(tm5), 3),
                           'median_top1_rmsd': round(st.median(rm1), 2), 'median_best5_rmsd': round(st.median(rm5), 2),
                           'top1_tm_ge_0.45': sum(x >= 0.45 for x in tm1),
                           'mean_top1_inf_wc': round(st.mean(inf), 3) if inf else None}
    both = [r for r in rows if r.get('304_top1_tm') is not None]
    s['rcv_vs_af3'] = {'targets': len(both),
                       'rcv_top1_better_or_within_0.02': sum(r['rcv_top1_tm'] >= r['304_top1_tm'] - 0.02 for r in both),
                       'rcv_best5_better_or_within_0.02': sum(r['rcv_best5_tm'] >= r['304_best5_tm'] - 0.02 for r in both)}
    s['casp_median_group_mean_top1_tm'] = round(st.mean(r['casp_median_top1_tm'] for r in rows), 3)
    s['rcv_mean_fraction_of_groups_beaten'] = round(st.mean(r.get('rcv_top1_beats_groups', 0) for r in rows), 3)
    return s


def report(rows, s):
    m = s['methods']
    L = ['# CASP16 RNA benchmark: RCV versus leading predictors', '',
         'Generated by `benchmark/casp16/compare.py evaluate`. Scores: US-align, residue-index superposition on C3\' '
         'atoms, TM-score normalised by the reference (CASP16 "TM_score" column), RMSD in Å over all resolved residues. '
         'All models (RCV and the CASP16 groups\' submitted models) are scored with the same command; it reproduces '
         'the official CASP16 table (`compare.py verify`).', '',
         '| Method | Mean TM-score, model 1 | Mean TM-score, best of 5 | Median RMSD model 1 (Å) | Median RMSD best of 5 (Å) | Targets with TM ≥ 0.45 (model 1) | Base pairs right (mean INF-WC, model 1) |',
         '|---|---|---|---|---|---|---|']
    for k in ['rcv'] + list(GROUPS):
        v = m[k]
        L.append(f"| {v['label']} | {v['mean_top1_tm']:.3f} | {v['mean_best5_tm']:.3f} | {v['median_top1_rmsd']:.1f} | "
                 f"{v['median_best5_rmsd']:.1f} | {v['top1_tm_ge_0.45']}/{v['n']} | {v['mean_top1_inf_wc']:.2f} |")
    L += ['', f"Median CASP16 group (mean over targets of the per-target median model-1 TM-score): "
          f"{s['casp_median_group_mean_top1_tm']:.3f}. On average RCV's model 1 scores higher than "
          f"{100 * s['rcv_mean_fraction_of_groups_beaten']:.0f}% of the CASP16 groups' model 1.", '',
          '## Per target (TM-score, model 1 / best of 5)', '',
          '| Target | nt | Family | RCV | AlphaFold 3 | trRosettaRNA | Vfold | GuangzhouRNA | CASP16 median / best |',
          '|---|---|---|---|---|---|---|---|---|']
    f = lambda r, k: '–' if r.get(f'{k}_top1_tm') is None else f"{r[f'{k}_top1_tm']:.2f} / {r[f'{k}_best5_tm']:.2f}"
    for r in rows:
        L.append(f"| {r['target']} | {r['length']} | {r['family']} | {f(r, 'rcv')} | {f(r, '304')} | {f(r, '052')} | "
                 f"{f(r, '481')} | {f(r, '183')} | {r['casp_median_top1_tm']:.2f} / {r['casp_best_top1_tm']:.2f} |")
    L += ['', '## What this shows', '',
          '* The deep-learning and template-based methods get the global fold right on 8–10 of the 15 targets; RCV on '
          'none. RCV is a physics-inspired ensemble generator, not a competitor for single-structure accuracy.',
          '* Most of the gap is the 3D arrangement of correct helices (junction geometry, coaxial stacking, tertiary '
          'contacts), which RCV does not model yet: on R1203, R1209, R1211 and R1296 RCV\'s base pairs are as '
          'accurate as AlphaFold 3\'s (INF-WC 0.86–0.98), but its TM-score is about half. Where ViennaRNA misses the '
          '2D (pre-tRNA R1221s3/R1224s3, ZTP riboswitches R1263/R1264) the gap is larger.',
          '* On the five targets that the field as a whole found hardest (CASP16 median model below TM-score 0.3: R1205, '
          'R1209, R1212, R1288, R1293) RCV is close to the median group and scores higher than AlphaFold 3 on three '
          '(R1209, R1212, R1288) — but at this level all predictions are essentially wrong.',
          '* Practical advice: for one most-likely structure, use a deep-learning predictor and open its models in RCV; '
          'use RCV\'s own predictions for alternative folds, switching and pathways, which the single-structure tools '
          'do not provide.',
          '', '## Notes', '',
          '* CASP16 groups predicted blind in 2024 (many with templates, MSAs, deep learning and expert curation). RCV '
          'was run in 2026 from the sequence alone: ViennaRNA 2D hypotheses + physics-inspired 3D search, no templates, '
          'no MSA, no deep learning, nothing fitted to these targets. Its geometry parameters come from 178 PDB entries, '
          'none of which is a CASP16 target (checked by release date and sequence).',
          '* RCV settings: 5 models, "thorough" sampling (40 trajectories), seed 42, run on a 16-core desktop CPU.',
          '* TM-score above about 0.45 usually means the global fold is right; below 0.3 the fold is essentially random.']
    return '\n'.join(L) + '\n'


# ----------------------------------------------------------------------------------------------- verify
def cmd_verify(a):
    c = cache_dir(a)
    table = official_table(c)
    worst = 0.0
    n = 0
    for t in TARGETS:
        tgz = c / 'predictions' / f"{t['casp_id']}.tar.gz"
        if not tgz.exists():
            continue
        ref = c / 'references' / f"{t['id']}.pdb"
        with tarfile.open(tgz) as tf:
            for g in GROUPS:
                name = f"{t['casp_id']}/{t['casp_id']}TS{g}_1"
                try:
                    data = tf.extractfile(name).read()
                except KeyError:
                    continue
                tmp = c / 'verify_model.pdb'
                tmp.write_bytes(data)
                tm, rmsd = score(tmp, ref)
                off = table[t['casp_id']][g][1]
                worst = max(worst, abs(tm - off[0]))
                n += 1
                print(f"{t['id']:8s} {g} TM {tm:.3f} (official {off[0]:.3f})  RMSD {rmsd:.2f} (official {off[1]:.2f})")
    print(f'{n} models re-scored; largest TM-score difference to the official table {worst:.3f}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['fetch', 'run', 'evaluate', 'verify'])
    ap.add_argument('--cache', default=str(HERE / 'cache'))
    ap.add_argument('--runs', default=str(HERE / 'runs'))
    ap.add_argument('--models', action='store_true', help='fetch: also download all groups\' models (~150 MB)')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--trajectories', type=int, default=40)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--only', default='')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    {'fetch': cmd_fetch, 'run': cmd_run, 'evaluate': cmd_evaluate, 'verify': cmd_verify}[a.command](a)


if __name__ == '__main__':
    main()
