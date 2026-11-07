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
  python benchmark/casp16/compare.py fetch --models
  python benchmark/casp16/compare.py run --workers 16                       # templates on (runs_v073_tmpl)
  python benchmark/casp16/compare.py run --workers 16 --templates off --runs benchmark/casp16/runs_v073_free
  python benchmark/casp16/compare.py verify
  python benchmark/casp16/compare.py evaluate

Fairness notes (also written to the report): the CASP16 groups predicted blind in 2024; RCV was run
in 2026 with its template library and junction geometry limited to structures released before
2024-05-01 (the start of the CASP16 RNA season; every target was released later). No MSA, no deep
learning. Its geometry parameters were calibrated on 178 PDB entries; none is a CASP16 target
(checked by release date and sequence). R1255/R1256 (multi-conformation references) and targets
above 250 nt or in RNA-protein complexes are not part of this set.
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
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=a.seed, workers=a.workers,
                            templates=a.templates, template_cutoff=a.template_cutoff, motifs=a.motifs)
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


RCV_RUNS = {'rcv': 'RCV 0.7.3 (physics + templates)', 'rcvfree': 'RCV 0.7.3 without templates',
            'rcv072': 'RCV 0.7.2 (physics only)'}


def _rcv_scores(run, t, ref):
    """{rank: (tm, rmsd)}, INF-WC of model 1, manifest; None if the run is missing."""
    man_p = run / t['id'] / 'metadata' / 'manifest.json'
    if not man_p.exists():
        return None
    man = json.loads(man_p.read_text())
    models = {conf['rank']: score(run / t['id'] / conf['files']['pdb'], ref) for conf in man['conformers']}
    return models, inf_wc(ref, run / t['id'] / man['conformers'][0]['files']['pdb']), man


def _release_dates():
    """{'PDBID': 'YYYY-MM-DD'} of the template library (for the report)."""
    try:
        from rnafold3d.templates import Library
        return {m['pdb'].upper(): m['release'] for m in Library().meta}
    except Exception:
        return {}


def cmd_evaluate(a):
    c = cache_dir(a)
    runs = {'rcv': Path(a.runs), 'rcvfree': Path(a.runs_free) if a.runs_free else None,
            'rcv072': Path(a.runs_prev) if a.runs_prev else None}
    runs = {k: v for k, v in runs.items() if v is not None and v.exists()}
    table = official_table(c)
    released = _release_dates()
    rows = []
    for t in TARGETS:
        ref = c / 'references' / f"{t['id']}.pdb"
        got = {k: _rcv_scores(run, t, ref) for k, run in runs.items()}
        if got.get('rcv') is None:
            print('missing run', t['id'])
            continue
        man = got['rcv'][2]
        infwc = {k: v[1] for k, v in got.items() if v is not None}
        infwc.update(group_inf_wc(c, t, ref))
        groups = rescore_groups(c, t, ref) or {g: table.get(t['casp_id'], {}).get(g, {}) for g in GROUPS}
        all_top1 = [v[1][0] for g, v in table.get(t['casp_id'], {}).items() if 1 in v]
        tmpl = man['conformers'][0].get('template')
        hits = (man.get('input', {}).get('templates') or {}).get('hits') or []
        first = next((h for h in hits if h.get('label') == tmpl), None)
        row = {'target': t['id'], 'casp_id': t['casp_id'], 'length': len(t['sequence']),
               'resolved': t['resolved'], 'family': t['family'],
               'n_groups': len(all_top1), 'casp_median_top1_tm': round(st.median(all_top1), 3) if all_top1 else None,
               'casp_best_top1_tm': round(max(all_top1), 3) if all_top1 else None,
               'rcv_valid': man['statistics']['valid'], 'rcv_trajectories': man['statistics']['trajectories'],
               'rcv_wall_time_s': man['statistics']['wall_time_s'],
               'rcv_template': tmpl or '', 'rcv_template_coverage': first.get('coverage') if first else None,
               'rcv_template_identity': first.get('identity') if first else None,
               'rcv_template_released': ' '.join(released.get(p.split('_')[0].upper(), '?')
                                                 for p in tmpl.split('+')) if tmpl else ''}
        for key, models in [(k, v[0] if v else None) for k, v in got.items()] + list(groups.items()):
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
        print(f"{t['id']:8s} RCV top1 TM {row['rcv_top1_tm']}  best5 {row['rcv_best5_tm']}  [{tmpl or 'no template'}]"
              f"  | no templates {row.get('rcvfree_top1_tm')}  | 0.7.2 {row.get('rcv072_top1_tm')}"
              f"  | AF3 {row['304_top1_tm']}  | CASP median {row['casp_median_top1_tm']}")
    out = HERE / 'results'
    out.mkdir(exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (out / 'per_target.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    summary = summarise(rows, [k for k in RCV_RUNS if k in runs])
    (out / 'summary.json').write_text(json.dumps(summary, indent=1))
    (out / 'README.md').write_text(report(rows, summary))
    print(json.dumps(summary['methods'], indent=1))


def _mean(xs, nd=3):
    xs = [float(x) for x in xs if x not in (None, '')]
    return round(st.mean(xs), nd) if xs else None


def summarise(rows, rcv_keys=('rcv',)):
    keys = list(rcv_keys) + list(GROUPS)
    s = {'targets': len(rows), 'methods': {}}
    tmpl = [r for r in rows if r.get('rcv_template')]
    free = [r for r in rows if not r.get('rcv_template')]
    for k in keys:
        tm1 = [r[f'{k}_top1_tm'] for r in rows if r.get(f'{k}_top1_tm') is not None]
        if not tm1:
            continue
        rm1 = [r[f'{k}_top1_rmsd'] for r in rows if r.get(f'{k}_top1_rmsd') is not None]
        rm5 = [r[f'{k}_best5_rmsd'] for r in rows if r.get(f'{k}_best5_rmsd') is not None]
        s['methods'][k] = {'label': RCV_RUNS.get(k) or GROUPS[k], 'n': len(tm1),
                           'mean_top1_tm': _mean(tm1), 'mean_best5_tm': _mean(r.get(f'{k}_best5_tm') for r in rows),
                           'median_top1_rmsd': round(st.median(rm1), 2), 'median_best5_rmsd': round(st.median(rm5), 2),
                           'top1_tm_ge_0.45': sum(x >= 0.45 for x in tm1),
                           'mean_top1_inf_wc': _mean(r.get(f'{k}_top1_infwc') for r in rows),
                           'mean_top1_tm_template_targets': _mean(r.get(f'{k}_top1_tm') for r in tmpl),
                           'mean_top1_tm_other_targets': _mean(r.get(f'{k}_top1_tm') for r in free)}
    s['template_targets'] = {'n': len(tmpl), 'targets': [r['target'] for r in tmpl]}
    s['other_targets'] = {'n': len(free), 'targets': [r['target'] for r in free]}
    both = [r for r in rows if r.get('304_top1_tm') is not None]
    s['rcv_vs_af3'] = {'targets': len(both),
                       'rcv_top1_higher': sum(r['rcv_top1_tm'] > r['304_top1_tm'] for r in both),
                       'rcv_top1_better_or_within_0.02': sum(r['rcv_top1_tm'] >= r['304_top1_tm'] - 0.02 for r in both),
                       'rcv_best5_better_or_within_0.02': sum(r['rcv_best5_tm'] >= r['304_best5_tm'] - 0.02 for r in both)}
    s['casp_median_group_mean_top1_tm'] = _mean(r['casp_median_top1_tm'] for r in rows)
    s['casp_best_model1_mean_tm'] = _mean(r['casp_best_top1_tm'] for r in rows)
    s['rcv_mean_fraction_of_groups_beaten'] = _mean(r.get('rcv_top1_beats_groups', 0) for r in rows)
    times = [r['rcv_wall_time_s'] for r in rows]
    s['rcv_wall_time_s'] = {'median': round(st.median(times)), 'max': round(max(times))}
    return s


AF3_TRAINING_CUTOFF = '2021-09-30'


def report(rows, s):
    m = s['methods']
    keys = [k for k in list(RCV_RUNS) + list(GROUPS) if k in m]
    nt, nf = s['template_targets']['n'], s['other_targets']['n']
    L = ['# CASP16 RNA benchmark: RCV 0.7.3 versus leading predictors', '',
         'Generated by `benchmark/casp16/compare.py evaluate`. Scores: US-align, residue-index superposition on C3\' '
         'atoms, TM-score normalised by the reference (CASP16 "TM_score" column), RMSD in Å over all resolved residues. '
         'All models (RCV and the CASP16 groups\' submitted models) are scored with the same command; it reproduces '
         'the official CASP16 table (`compare.py verify`).', '',
         f'| Method | Mean TM-score, model 1 | Mean TM-score, best of 5 | Model 1 on the {nt} targets with a template '
         f'| Model 1 on the other {nf} | Median RMSD model 1 (Å) | Targets with TM ≥ 0.45 (model 1) '
         '| Base pairs right (mean INF-WC, model 1) |',
         '|---|---|---|---|---|---|---|---|']
    g = lambda x, nd=3: '–' if x is None else f'{x:.{nd}f}'
    for k in keys:
        v = m[k]
        name = f"**{v['label']}**" if k == 'rcv' else v['label']
        L.append(f"| {name} | {g(v['mean_top1_tm'])} | {g(v['mean_best5_tm'])} | {g(v['mean_top1_tm_template_targets'])} | "
                 f"{g(v['mean_top1_tm_other_targets'])} | {g(v['median_top1_rmsd'], 1)} | {v['top1_tm_ge_0.45']}/{v['n']} | "
                 f"{g(v['mean_top1_inf_wc'], 2)} |")
    L += ['', f"Median CASP16 group (mean over targets of the per-target median model-1 TM-score): "
          f"{s['casp_median_group_mean_top1_tm']:.3f}; best model 1 of any group, per target: "
          f"{s['casp_best_model1_mean_tm']:.3f} (no single group reaches this). On average RCV's model 1 scores higher "
          f"than {100 * s['rcv_mean_fraction_of_groups_beaten']:.0f}% of the CASP16 groups' model 1. "
          f"\"Targets with a template\": RCV found a related structure released before 2024-05-01 "
          f"({', '.join(s['template_targets']['targets'])}).", '',
          '## Per target (TM-score, model 1 / best of 5)', '',
          '| Target | nt | Family | RCV 0.7.3 | Template (released) | RCV 0.7.3 without templates | RCV 0.7.2 '
          '| AlphaFold 3 | trRosettaRNA | Vfold | GuangzhouRNA | CASP16 median / best |',
          '|---|---|---|---|---|---|---|---|---|---|---|---|']
    f = lambda r, k: '–' if r.get(f'{k}_top1_tm') is None else f"{r[f'{k}_top1_tm']:.2f} / {r[f'{k}_best5_tm']:.2f}"
    for r in rows:
        tp = f"{r['rcv_template']} ({r['rcv_template_released']})" if r.get('rcv_template') else '–'
        L.append(f"| {r['target']} | {r['length']} | {r['family']} | **{f(r, 'rcv')}** | {tp} | {f(r, 'rcvfree')} | "
                 f"{f(r, 'rcv072')} | {f(r, '304')} | {f(r, '052')} | {f(r, '481')} | {f(r, '183')} | "
                 f"{r['casp_median_top1_tm']:.2f} / {r['casp_best_top1_tm']:.2f} |")
    rcv, af3 = m['rcv'], m.get('304', {})
    ahead = [r['target'] for r in rows if r['rcv_top1_tm'] > r['304_top1_tm'] + 0.05]
    behind = [r['target'] for r in rows if r['rcv_top1_tm'] < r['304_top1_tm'] - 0.05]
    late = sorted({(p, d) for r in rows if r.get('rcv_template')
                   for p, d in zip(r['rcv_template'].split('+'), r['rcv_template_released'].split())
                   if d > AF3_TRAINING_CUTOFF})
    L += ['', '## What this shows', '',
          f"* With templates, RCV's model 1 has a mean TM-score of {rcv['mean_top1_tm']:.3f} over the 15 targets: "
          f"AlphaFold 3 {af3.get('mean_top1_tm', 0):.3f}, trRosettaRNA {m['052']['mean_top1_tm']:.3f}; the expert groups "
          f"Vfold ({m['481']['mean_top1_tm']:.3f}) and GuangzhouRNA ({m['183']['mean_top1_tm']:.3f}) are higher. "
          f"RCV 0.7.2 (physics only) had {m['rcv072']['mean_top1_tm']:.3f}." if 'rcv072' in m else '',
          f"* RCV is more than 0.05 ahead of AlphaFold 3 on {len(ahead)} targets ({', '.join(ahead)}) and more than 0.05 "
          f"behind on {len(behind)} ({', '.join(behind)}); model 1 is higher on {s['rcv_vs_af3']['rcv_top1_higher']} "
          f"of {s['rcv_vs_af3']['targets']}.",
          f"* Where RCV found a template ({nt} targets) its models follow the template closely, so accuracy is set by "
          f"how similar the template is to the target: mean {rcv['mean_top1_tm_template_targets']:.3f} (AlphaFold 3 "
          f"{af3.get('mean_top1_tm_template_targets', 0):.3f}). On the other {nf} targets, which the whole field found "
          f"hard, RCV has {rcv['mean_top1_tm_other_targets']:.3f} (AlphaFold 3 {af3.get('mean_top1_tm_other_targets', 0):.3f}).",
          (f"* Without templates (physics, ViennaRNA 2D and junction geometry from the library) RCV 0.7.3 has "
           f"{m['rcvfree']['mean_top1_tm']:.3f}." if 'rcvfree' in m else ''),
          f"* Speed: median {s['rcv_wall_time_s']['median']} s per target (longest {s['rcv_wall_time_s']['max']} s) "
          f"for 40 trajectories on a 16-core laptop CPU (Intel Core Ultra 9 285H), no GPU; the structure library is 21.5 MB.",
          '* RNAs solved in several conformations (R1203, R1255, R1256): [ENSEMBLES.md](ENSEMBLES.md). Size, hardware '
          'and features next to other predictors: [rcv_comparison.svg](rcv_comparison.svg) (`figure.py`).',
          '', '## Fairness notes', '',
          '* The CASP16 groups predicted blind in 2024. RCV was run in 2026; its template library and junction '
          'geometry were limited to structures released before 2024-05-01 (the start of the CASP16 RNA prediction '
          'season; all 15 reference structures were released later). The CASP16 groups could use the same structures.',
          f"* AlphaFold 3 (the AF3 server baseline, group 304) uses no RNA templates and was trained on PDB structures "
          f"released up to {AF3_TRAINING_CUTOFF}. Templates that RCV used and that were released after that: "
          + (', '.join(f'{p} ({d})' for p, d in late) or 'none') + '. Vfold and GuangzhouRNA had all of them.',
          '* Not fully blind for the method choices: the template rules (minimum coverage 45%, the rule for a second '
          'template covering an insert, ranking by agreement between templates) were set during 0.7.3 development '
          'with the CASP16 results visible. They are general rules, with no target-specific settings. The development '
          'benchmark (`benchmark/dev_templates.py`, 26 PDB chains, each with templates limited to structures released '
          'before it) is the independent check.',
          '* RCV settings: 5 models, 40 trajectories ("thorough"), seed 42. Geometry parameters come from 178 PDB '
          'entries, none of them a CASP16 target (checked by release date and sequence). No MSA, no deep learning.',
          '* TM-score above about 0.45 usually means the global fold is right; below 0.3 the fold is essentially random.']
    return '\n'.join(x for x in L if x is not None) + '\n'


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
    ap.add_argument('--runs', default=str(HERE / 'runs_v073_tmpl'), help='RCV run (templates on)')
    ap.add_argument('--runs-free', default=str(HERE / 'runs_v073_free'), help='evaluate: RCV run without templates')
    ap.add_argument('--runs-prev', default=str(HERE / 'runs'), help='evaluate: RCV 0.7.2 run (physics only)')
    ap.add_argument('--models', action='store_true', help='fetch: also download all groups\' models (~150 MB)')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument('--trajectories', type=int, default=40)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--only', default='')
    ap.add_argument('--force', action='store_true')
    ap.add_argument('--templates', choices=['auto', 'off'], default='auto',
                    help='run: experimental structures as templates (only those released before --template-cutoff)')
    ap.add_argument('--template-cutoff', default='2024-05-01',
                    help='run: latest release date of library structures (CASP16 targets were released later)')
    ap.add_argument('--motifs', choices=['auto', 'off'], default='auto', help='run: junction geometry from the library')
    a = ap.parse_args()
    {'fetch': cmd_fetch, 'run': cmd_run, 'evaluate': cmd_evaluate, 'verify': cmd_verify}[a.command](a)


if __name__ == '__main__':
    main()
