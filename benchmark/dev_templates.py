"""Development benchmark for template-based and template-free modelling (pilot dataset, 26 chains).

Every target is predicted blind with respect to its own structure: templates and junction motifs are
limited to library structures released before the target's PDB release date. Scores: US-align TM-score
(residue-index superposition, normalised by the reference, as in compare.py for CASP16).

    python benchmark/dev_templates.py run  --label tmpl --templates auto --trajectories 16 --workers 16
    python benchmark/dev_templates.py report --labels notmpl,tmpl
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'src'))
SAMPLES = HERE / 'data' / 'samples'
RUNS = HERE / 'runs' / 'dev_templates'


def release_dates():
    from rnafold3d.templates import load_library
    lib = load_library()
    dates = {}
    for m in lib.meta:
        dates.setdefault(m['pdb'].upper(), m['release'])
    return dates


def targets(only=None):
    out = []
    for d in sorted(SAMPLES.iterdir()):
        meta = json.loads((d / 'metadata.json').read_text())
        if only and d.name not in only:
            continue
        seq = ''.join(l.strip() for l in (d / 'sequence.fasta').read_text().splitlines() if not l.startswith('>'))
        out.append({'id': d.name, 'pdb': meta['pdb_id'], 'family': meta.get('family', ''), 'sequence': seq,
                    'reference': d / 'reference.pdb'})
    return out


def tm(model, ref):
    exe = shutil.which('USalign') or str(Path(sys.executable).parent / 'USalign')
    p = subprocess.run([exe, str(model), str(ref), '-mol', 'RNA', '-TMscore', '1', '-outfmt', '2'],
                       capture_output=True, text=True, check=True)
    row = [l for l in p.stdout.splitlines() if l and not l.startswith('#')][-1].split('\t')
    return float(row[3]), float(row[4])


def cmd_run(a):
    from rnafold3d.predict import predict, PredictConfig
    dates = release_dates()
    import urllib.request
    for t in targets(a.only.split(',') if a.only else None):
        out = RUNS / a.label / t['id']
        if (out / 'metadata' / 'manifest.json').exists() and not a.force:
            continue
        if out.exists():
            shutil.rmtree(out)
        cutoff = dates.get(t['pdb'].upper())
        if not cutoff:
            q = '{entry(entry_id:"%s"){rcsb_accession_info{initial_release_date}}}' % t['pdb']
            req = urllib.request.Request('https://data.rcsb.org/graphql', data=json.dumps({'query': q}).encode(),
                                         headers={'Content-Type': 'application/json'})
            cutoff = json.loads(urllib.request.urlopen(req, timeout=60).read())['data']['entry'][
                'rcsb_accession_info']['initial_release_date'][:10]
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=a.seed, workers=a.workers,
                            templates=a.templates, template_cutoff=cutoff, motifs=a.motifs)
        print(f"== {t['id']} ({len(t['sequence'])} nt, cutoff {cutoff})", flush=True)
        man = predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
        print('   ', man['statistics'], [h['label'] for h in man['input']['templates'].get('hits', [])], flush=True)


def evaluate(label):
    rows = {}
    for t in targets():
        run = RUNS / label / t['id']
        mp = run / 'metadata' / 'manifest.json'
        if not mp.exists():
            continue
        man = json.loads(mp.read_text())
        tms = [tm(run / c['files']['pdb'], t['reference'])[0] for c in man['conformers']]
        rows[t['id']] = {'top1': tms[0] if tms else 0.0, 'best5': max(tms) if tms else 0.0,
                         'templates': [h['label'] for h in man['input'].get('templates', {}).get('hits', [])],
                         'source1': man['conformers'][0].get('template', 'physics') if man['conformers'] else '',
                         'time': man['statistics']['wall_time_s']}
    return rows


def cmd_report(a):
    labels = a.labels.split(',')
    res = {l: evaluate(l) for l in labels}
    ids = sorted(set().union(*[set(r) for r in res.values()]))
    print(f"{'target':10s} " + ' '.join(f'{l:>22s}' for l in labels))
    for i in ids:
        cells = []
        for l in labels:
            r = res[l].get(i)
            cells.append(f"{r['top1']:.2f}/{r['best5']:.2f} {r['source1'][:8]:>8s} {r['time']:5.0f}s" if r else ' ' * 22)
        print(f'{i:10s} ' + ' '.join(cells))
    for l in labels:
        v = list(res[l].values())
        if v:
            print(f"{l}: n={len(v)} mean top1 {sum(x['top1'] for x in v) / len(v):.3f} best5 "
                  f"{sum(x['best5'] for x in v) / len(v):.3f}  with template: {sum(bool(x['templates']) for x in v)}")
    (RUNS / 'report.json').write_text(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--label', required=True)
    r.add_argument('--templates', default='auto')
    r.add_argument('--motifs', default='auto')
    r.add_argument('--trajectories', type=int, default=16)
    r.add_argument('--workers', type=int, default=16)
    r.add_argument('--seed', type=int, default=42)
    r.add_argument('--only', default='')
    r.add_argument('--force', action='store_true')
    p = sub.add_parser('report')
    p.add_argument('--labels', required=True)
    a = ap.parse_args()
    {'run': cmd_run, 'report': cmd_report}[a.cmd](a)


if __name__ == '__main__':
    main()
