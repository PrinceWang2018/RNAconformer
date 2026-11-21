"""Upper bound for the template search: for every query of benchmark/template_selection.py, the structural
similarity (US-align, sequence-independent, C4' atoms, TM-score normalised by the query) of every library
chain to the query's own experimental structure. One core, resumable; the matrix is the ground truth that
search and ranking rules are judged against (which chain was the best template that existed at the time).

    python benchmark/template_oracle.py scan  --selection benchmark/runs/template_selection/v074.json
    python benchmark/template_oracle.py report --selection benchmark/runs/template_selection/v074.json
"""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'benchmark' / 'casp16'))

from oracle_templates import write_trace, usalign_many          # noqa: E402
from rnafold3d.templates import load_library                    # noqa: E402

OUT = ROOT / 'benchmark' / 'runs' / 'template_selection' / 'oracle'


def chain_index(lib):
    return {f"{m['pdb']}_{m['chain']}": k for k, m in enumerate(lib.meta)}


def cmd_scan(a):
    lib = load_library()
    idx = chain_index(lib)
    lens = np.diff(lib.off)
    queries = json.loads(Path(a.selection).read_text())['queries']
    queries.sort(key=lambda r: bool(r['structure']))             # queries without a hit first
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        names = {}
        for k in range(lib.n_chains):
            if lens[k] <= a.max_length and write_trace(tmp / f'c{k}.pdb', lib.chain_seq(k), lib.chain_xyz(k)) >= 15:
                names[f'c{k}.pdb'] = k
        for n, r in enumerate(queries):
            dest = OUT / f"{r['query']}.npy"
            if dest.exists():
                continue
            k = idx[r['query']]
            L = int(lens[k])
            use = [nm for nm, j in names.items() if lens[j] >= 0.4 * L and j != k]
            res = usalign_many(tmp, use, tmp / f'c{k}.pdb')
            tm = np.full(lib.n_chains, np.nan, np.float32)
            for nm, (t, _rmsd, _n) in res.items():
                tm[names[nm]] = t
            np.save(dest, tm)
            if n % 20 == 0:
                print(f'{n}/{len(queries)} {r["query"]} {L} nt best {np.nanmax(tm):.3f}', flush=True)


def allowed_tm(lib, r, tm):
    m = lib.meta[chain_index(lib)[r['query']]]
    ok = lib.allowed(m['release'], (m['pdb'],))
    return np.where(ok & ~np.isnan(tm), tm, -1.0)


def cmd_report(a):
    lib = load_library()
    queries = json.loads(Path(a.selection).read_text())['queries']
    idx = chain_index(lib)
    rows = []
    for r in queries:
        f = OUT / f"{r['query']}.npy"
        if not f.exists():
            continue
        tm = allowed_tm(lib, r, np.load(f))
        best = int(np.argmax(tm))
        hits = r['structure']
        first = float(tm[idx[f"{hits[0]['pdb']}_{hits[0]['chain']}"]]) if hits else None
        found = max((float(tm[idx[f"{h['pdb']}_{h['chain']}"]]) for h in hits), default=None)
        rows.append({'query': r['query'], 'length': r['length'], 'oracle': float(tm[best]),
                     'oracle_chain': f"{lib.meta[best]['pdb']}_{lib.meta[best]['chain']}", 'first': first, 'found': found})
    hit = [x for x in rows if x['first'] is not None]
    no = [x for x in rows if x['first'] is None]
    print(f'{len(rows)} queries scanned: {len(hit)} with a hit, {len(no)} without')
    if hit:
        print(f"with a hit   : best chain in the library {np.mean([x['oracle'] for x in hit]):.3f}, first-ranked template "
              f"{np.mean([x['first'] for x in hit]):.3f}, best template found {np.mean([x['found'] for x in hit]):.3f} (structural TM)")
        print(f"               first-ranked template more than 0.10 below the best in the library: "
              f"{sum(x['oracle'] - x['first'] > 0.10 for x in hit)} of {len(hit)}")
    if no:
        o = np.array([x['oracle'] for x in no])
        print(f'without a hit: best chain in the library {o.mean():.3f}; >= 0.45: {(o >= 0.45).sum()} of {len(o)}; '
              f'>= 0.55: {(o >= 0.55).sum()}; >= 0.70: {(o >= 0.70).sum()}')
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('command', choices=['scan', 'report'])
    ap.add_argument('--selection', default=str(ROOT / 'benchmark' / 'runs' / 'template_selection' / 'v074.json'))
    ap.add_argument('--max-length', type=int, default=1200)
    ap.add_argument('--out', default='')
    a = ap.parse_args()
    {'scan': cmd_scan, 'report': cmd_report}[a.command](a)
