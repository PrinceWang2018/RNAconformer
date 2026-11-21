"""Template-based modelling benchmark: 3D models for library chains that have a template (larger than the
26-chain pilot of dev_templates.py, where only six chains have one).

Queries are taken from a run of benchmark/template_selection.py: chains released on or after its --since date
for which the search finds a template among the structures released before them (their own PDB entry excluded).
Each is predicted blind in the same way (templates and junction geometry limited to earlier structures) and the
models are compared with the chain's own experimental C4' atoms: US-align, residue-index superposition, TM-score
normalised by the experimental chain, as in dev_templates.py (one atom per residue is given to US-align under
the name C3', for the model and the reference alike). PDB entries of the CASP16 targets are left out, so that
CASP16 stays a blind check for anything tuned here.

    python benchmark/template_models.py run --label fit4 --workers 16
    python benchmark/template_models.py run --label fit1 --workers 16 --set template_fit_attempts=1
    python benchmark/template_models.py report --labels fit1,fit4
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'src'))
RUNS = HERE / 'runs' / 'template_models'
SELECTION = HERE / 'runs' / 'template_selection' / 'v074.json'


def targets(a):
    from rnafold3d.templates import load_library
    lib = load_library()
    idx = {f"{m['pdb']}_{m['chain']}": k for k, m in enumerate(lib.meta)}
    casp = {t['pdb'].upper() for t in json.loads((HERE / 'casp16' / 'targets.json').read_text()) if t.get('pdb')}
    rows = [r for r in json.loads(Path(a.selection).read_text())['queries']
            if r['structure'] and a.min_len <= r['length'] <= a.max_len and r['query'].split('_')[0].upper() not in casp]
    rows = rows[::max(1, len(rows) // a.n)][:a.n]
    out = []
    for r in rows:
        k = idx[r['query']]
        m = lib.meta[k]
        out.append({'id': r['query'], 'index': k, 'pdb': m['pdb'], 'release': m['release'], 'sequence': lib.chain_seq(k),
                    'title': m['title']})
    return lib, out


def write_trace(path, seq, xyz):
    lines = [f"ATOM  {i:5d}  C3'   {b if b in 'ACGU' else 'N'} A{i:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00           C"
             for i, (b, x) in enumerate(zip(seq, xyz), 1) if not np.isnan(x[0])]
    Path(path).write_text('\n'.join(lines) + '\nTER\nEND\n')


def model_c4(path, n):
    xyz = np.full((n, 3), np.nan)
    for line in Path(path).read_text().splitlines():
        if line.startswith('ATOM') and line[12:16].strip() == "C4'":
            i = int(line[22:26]) - 1
            if 0 <= i < n:
                xyz[i] = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
    return xyz


def tm(model, ref):
    exe = shutil.which('USalign') or str(Path(sys.executable).parent / 'USalign')
    p = subprocess.run([exe, str(model), str(ref), '-mol', 'RNA', '-TMscore', '1', '-outfmt', '2'],
                       capture_output=True, text=True, check=True)
    return float([l for l in p.stdout.splitlines() if l and not l.startswith('#')][-1].split('\t')[3])


def confidence(path, seq):
    """Predicted lDDT of a model (the confidence model of rnafold3d.qa, from the coordinates alone)."""
    from rnafold3d import qa
    from rnafold3d.structure_io import read_pdb
    chain = max(read_pdb(path)[0].values(), key=lambda ch: len(ch.residues))
    p = qa.predict(seq, chain.residues)
    return None if p is None else round(float(p['global']), 4)


def _extra(a):
    out = {}
    for kv in filter(None, (a.set or '').split(',')):
        k, v = kv.split('=', 1)
        out[k] = {'true': True, 'false': False}.get(v.lower(), int(v) if v.lstrip('-').isdigit() else v)
    return out


def cmd_run(a):
    from rnafold3d.predict import predict, PredictConfig
    _lib, ts = targets(a)
    print(f'{len(ts)} targets', flush=True)
    for t in ts:
        out = RUNS / a.label / t['id']
        if (out / 'metadata' / 'manifest.json').exists():
            continue
        if out.exists():
            shutil.rmtree(out)
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=a.seed, workers=a.workers,
                            template_cutoff=t['release'], template_exclude=(t['pdb'],), **_extra(a))
        print(f"== {t['id']} ({len(t['sequence'])} nt, before {t['release']})", flush=True)
        man = predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
        print('   ', man['statistics'], [h['label'] for h in man['input']['templates'].get('hits', [])], flush=True)


def evaluate(a, label, lib, ts):
    cache = RUNS / label / 'scores.json'
    rows = json.loads(cache.read_text()) if cache.exists() else {}
    tmp = RUNS / label / '_traces'
    tmp.mkdir(parents=True, exist_ok=True)
    for t in ts:
        run = RUNS / label / t['id']
        mp = run / 'metadata' / 'manifest.json'
        if t['id'] in rows or not mp.exists():
            continue
        man = json.loads(mp.read_text())
        n = len(t['sequence'])
        write_trace(tmp / 'ref.pdb', t['sequence'], lib.chain_xyz(t['index']))
        models = []
        for c in man['conformers']:
            write_trace(tmp / 'model.pdb', t['sequence'], model_c4(run / c['files']['pdb'], n))
            fit = c.get('template_fit') or {}
            models.append({'tm': tm(tmp / 'model.pdb', tmp / 'ref.pdb'), 'template': c.get('template'),
                           'fit_accepted': fit.get('accepted'), 'released': fit.get('released_residues', 0),
                           'confidence': confidence(run / c['files']['pdb'], t['sequence'])})
        rows[t['id']] = {'length': n, 'models': models, 'valid': man['statistics']['valid'],
                         'templates': [h['label'] for h in man['input'].get('templates', {}).get('hits', [])],
                         'time': man['statistics']['wall_time_s']}
    shutil.rmtree(tmp, ignore_errors=True)
    cache.write_text(json.dumps(rows, indent=1))
    return rows


def cmd_report(a):
    lib, ts = targets(a)
    labels = a.labels.split(',')
    res = {lab: evaluate(a, lab, lib, ts) for lab in labels}
    common = [t['id'] for t in ts if all(t['id'] in res[lab] and res[lab][t['id']]['models'] for lab in labels)]
    print(f'{len(common)} targets with models in every run ({", ".join(labels)})')
    top = {lab: np.array([res[lab][q]['models'][0]['tm'] for q in common]) for lab in labels}
    best = {lab: np.array([max(m['tm'] for m in res[lab][q]['models']) for q in common]) for lab in labels}
    print(f"{'run':12s} {'model 1':>8s} {'best of 5':>10s} {'model 1 from a template':>24s} {'rejected fits':>14s} {'rescued fits':>13s} {'time, s':>8s}")
    for lab in labels:
        ms = [m for q in common for m in res[lab][q]['models'] if m['template']]
        print(f"{lab:12s} {top[lab].mean():8.3f} {best[lab].mean():10.3f} "
              f"{sum(bool(res[lab][q]['models'][0]['template']) for q in common):24d} "
              f"{sum(m['fit_accepted'] is False for m in ms):9d}/{len(ms):<4d} {sum(bool(m['released']) for m in ms):13d} "
              f"{sum(res[lab][q]['time'] for q in common):8.0f}")
    for lab in labels:                    # would the confidence model choose a better first model?
        rows = [res[lab][q]['models'] for q in common]
        if not all(m.get('confidence') is not None for ms in rows for m in ms):
            continue
        first = np.array([ms[0]['tm'] for ms in rows])
        line = f"{lab}: model 1 {first.mean():.3f}; model with the highest confidence {np.mean([max(ms, key=lambda m: m['confidence'])['tm'] for ms in rows]):.3f}"
        for gap in (0.05, 0.1, 0.15, 0.2):
            pick = []
            for ms in rows:
                b = max(ms, key=lambda m: m['confidence'])
                pick.append(b['tm'] if b['confidence'] - ms[0]['confidence'] >= gap else ms[0]['tm'])
            d = np.array(pick) - first
            line += f"; ahead by >= {gap}: {np.mean(pick):.3f} ({int((d > 0.02).sum())} better, {int((d < -0.02).sum())} worse)"
        print(line)
    if len(labels) == 2:
        rng = np.random.default_rng(1)
        for name, v in (('model 1', top), ('best of 5', best)):
            d = v[labels[1]] - v[labels[0]]
            boot = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(20000)])
            print(f'{name}: {labels[1]} - {labels[0]} = {d.mean():+.4f} (95% interval {np.percentile(boot, 2.5):+.4f} to '
                  f'{np.percentile(boot, 97.5):+.4f}); better by > 0.02 in {int((d > 0.02).sum())}, worse in {int((d < -0.02).sum())}')
        d = top[labels[1]] - top[labels[0]]
        for i in np.argsort(-np.abs(d))[:12]:
            if abs(d[i]) > 0.02:
                q = common[i]
                print(f"   {q:10s} {res[labels[0]][q]['length']:4d} nt  {top[labels[0]][i]:.3f} -> {top[labels[1]][i]:.3f}  "
                      f"{res[labels[0]][q]['models'][0]['template']} -> {res[labels[1]][q]['models'][0]['template']}")
    (RUNS / 'report.json').write_text(json.dumps({lab: res[lab] for lab in labels}, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['run', 'report'])
    ap.add_argument('--label', default='run')
    ap.add_argument('--labels', default='')
    ap.add_argument('--selection', default=str(SELECTION))
    ap.add_argument('-n', type=int, default=48, help='number of targets (evenly spaced over the eligible queries)')
    ap.add_argument('--min-len', type=int, default=40)
    ap.add_argument('--max-len', type=int, default=140)
    ap.add_argument('--trajectories', type=int, default=16)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--set', default='', help='PredictConfig fields, key=value,...')
    a = ap.parse_args()
    {'run': cmd_run, 'report': cmd_report}[a.command](a)


if __name__ == '__main__':
    main()
