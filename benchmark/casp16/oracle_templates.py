"""How good a template does the library hold for every CASP16 target, and does the search find it?

For each target the experimental structure is superposed (US-align, sequence-independent, C4' atoms, TM-score
normalised by the target) on every library chain released before the cut-off. That gives the best template
the library could offer (an upper bound no search can exceed without a larger library), next to the templates
the engine's search returns from the sequence alone and the TM-score those have against the target.

    python benchmark/casp16/oracle_templates.py scan  [--only R1271,R1203]   # slow part, one core, cached
    python benchmark/casp16/oracle_templates.py report

Output: benchmark/casp16/results/oracle_templates.json and a table on stdout. Reading it: "best in library"
well above "first found" means the search or the ranking loses accuracy; both low means no template exists.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d.templates import load_library, search, rank_templates, sw_scan      # noqa: E402

TARGETS = json.loads((HERE / 'targets.json').read_text())
CUTOFF = '2024-05-01'
OUT = HERE / 'results' / 'oracle_templates.json'
CACHE = HERE / 'cache' / 'oracle'


def write_trace(path, seq, xyz):
    """One C3'-named atom per residue (US-align's RNA atom), residues without coordinates left out."""
    lines = []
    for i, (b, x) in enumerate(zip(seq, xyz), 1):
        if np.isnan(x[0]):
            continue
        lines.append(f"ATOM  {i % 100000:5d}  C3'   {b if b in 'ACGU' else 'N'} A{i % 10000:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00           C")
    Path(path).write_text('\n'.join(lines) + '\nTER\nEND\n')
    return len(lines)


def reference_trace(tid):
    """Sequence and C4' coordinates of the reference structure (resolved residues only)."""
    seq, xyz = [], []
    for line in (HERE / 'cache' / 'references' / f'{tid}.pdb').read_text().splitlines():
        if line.startswith(('ATOM', 'HETATM')) and line[12:16].strip() == "C4'":
            seq.append(line[17:20].strip()[-1])
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return ''.join(seq), np.array(xyz)


def usalign_many(folder, names, ref):
    """TM-score normalised by `ref` of every structure in `names` (files in folder) against ref."""
    lst = Path(folder) / 'list.txt'
    lst.write_text('\n'.join(names) + '\n')
    r = subprocess.run(['USalign', '-dir1', str(folder) + '/', str(lst), str(ref), '-outfmt', '2', '-mol', 'RNA'],
                       capture_output=True, text=True)
    out = {}
    for line in r.stdout.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        out[Path(f[0].split(':')[0]).name] = (float(f[3]), float(f[4]), int(f[10]))      # TM by ref, RMSD, aligned
    return out


def cmd_scan(a):
    lib = load_library()
    lens = np.diff(lib.off)
    allowed = lib.allowed(CUTOFF)
    CACHE.mkdir(parents=True, exist_ok=True)
    only = {x for x in a.only.split(',') if x}
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        names = {}
        for k in np.flatnonzero(allowed):
            k = int(k)
            if lens[k] > a.max_length:
                continue
            nm = f'c{k}.pdb'
            if write_trace(tmp / nm, lib.chain_seq(k), lib.chain_xyz(k)) >= 15:
                names[nm] = k
        for t in TARGETS:
            tid = t['id']
            dest = CACHE / f'{tid}.json'
            if (only and tid not in only) or (dest.exists() and not a.force):
                continue
            seq, xyz = reference_trace(tid)
            ref = tmp / f'ref_{tid}.pdb'
            write_trace(ref, seq, xyz)
            L = len(seq)
            use = [nm for nm, k in names.items() if 0.4 * L <= lens[k]]
            res = usalign_many(tmp, use, ref)
            rows = sorted(({'chain_index': names[nm], 'pdb': lib.meta[names[nm]]['pdb'], 'chain': lib.meta[names[nm]]['chain'],
                            'release': lib.meta[names[nm]]['release'], 'length': int(lens[names[nm]]),
                            'title': lib.meta[names[nm]]['title'][:70], 'tm': tm, 'rmsd': rmsd, 'aligned': n}
                           for nm, (tm, rmsd, n) in res.items()), key=lambda r: -r['tm'])
            dest.write_text(json.dumps({'target': tid, 'length': L, 'chains_compared': len(res), 'top': rows[:25],
                                        'all': {str(r['chain_index']): r['tm'] for r in rows}}, indent=1))
            print(f'{tid}: {len(res)} chains; best {rows[0]["pdb"]}_{rows[0]["chain"]} TM {rows[0]["tm"]:.3f}', flush=True)


def cmd_report(a):
    lib = load_library()
    out = []
    print(f'{"target":8s} {"nt":>4s}  {"best in library":28s} {"first found by the search":34s} {"best found":18s}')
    for t in TARGETS:
        tid = t['id']
        f = CACHE / f'{tid}.json'
        if not f.exists():
            continue
        o = json.loads(f.read_text())
        seq = t['sequence']
        hits = search(seq, lib, cutoff=CUTOFF, max_hits=8, scan=sw_scan(seq, lib), structure=True)
        hits = rank_templates(seq, hits, lib, 3)
        tm_of = lambda h: o['all'].get(str(h['chain_index']))
        best = o['top'][0]
        first = hits[0] if hits else None
        found_best = max(hits, key=lambda h: tm_of(h) or 0) if hits else None
        row = {'target': tid, 'length': o['length'], 'best_in_library': {k: best[k] for k in ('pdb', 'chain', 'tm', 'release', 'title')},
               'top5': [{k: r[k] for k in ('pdb', 'chain', 'tm', 'length', 'title')} for r in o['top'][:5]],
               'found': [{'pdb': h['pdb'], 'chain': h['chain'], 'tm': tm_of(h), 'identity': h['identity'], 'coverage': h['coverage'],
                          'rank_score': h['rank_score'], 'found_by': h['found_by']} for h in hits]}
        out.append(row)
        fs = f"{first['pdb']}_{first['chain']} {tm_of(first) or float('nan'):.3f} (id {first['identity']:.2f} cov {first['coverage']:.2f})" if first else '–'
        bs = f"{found_best['pdb']}_{found_best['chain']} {tm_of(found_best) or float('nan'):.3f}" if found_best else '–'
        print(f"{tid:8s} {o['length']:4d}  {best['pdb'] + '_' + best['chain']:10s} {best['tm']:.3f} {best['title'][:10]:10s}  {fs:34s} {bs:18s}")
    OUT.write_text(json.dumps(out, indent=1))
    b = [r['best_in_library']['tm'] for r in out]
    f1 = [(r['found'][0]['tm'] or 0.0) if r['found'] else 0.0 for r in out]
    print(f'\nmean TM of the best library chain {np.mean(b):.3f}; of the first template found {np.mean(f1):.3f} '
          f'(0 where none); targets with a library chain >= 0.45: {sum(x >= 0.45 for x in b)} of {len(b)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('command', choices=['scan', 'report'])
    ap.add_argument('--only', default='')
    ap.add_argument('--max-length', type=int, default=1200, help='scan: skip library chains longer than this')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    {'scan': cmd_scan, 'report': cmd_report}[a.command](a)
