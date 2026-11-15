"""Template-selection benchmark: leave-one-out over the built-in structure library (no 3D modelling).

Every library chain released on or after --since (40-250 nt) is a query; templates are limited to chains
released before it (its own entry excluded). For every hit the benchmark computes what the template is
worth: the TM-score-like agreement between the template's C4' atoms and the query's own experimental C4'
atoms over the aligned residues, normalised by the query length (unaligned residues count as zero). This
measures the search and ranking rules on hundreds of RNAs in minutes, independent of CASP16.

    python benchmark/template_selection.py run --out benchmark/runs/template_selection/v074.json
    python benchmark/template_selection.py report benchmark/runs/template_selection/v074.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

FEATURES = ('pdb', 'chain', 'evalue', 'identity', 'coverage', 'pairs_mapped', 'pairs_compatible', 'template_pairs',
            'found_by', 'resolution', 'method', 'template_length', 'consensus', 'rank_score')


def tm_aligned(A, B, L):
    """TM-score of two matched C4' sets after weighted superposition, normalised by L."""
    if len(A) < 5:
        return 0.0
    d0 = max(0.5, 0.6 * (max(L, 13) - 0.5) ** 0.5 - 2.5)          # RNA d0 of US-align
    w = np.ones(len(A))
    best = 0.0
    for _ in range(12):
        a0 = (A * w[:, None]).sum(0) / w.sum()
        b0 = (B * w[:, None]).sum(0) / w.sum()
        u, _s, vt = np.linalg.svd(((A - a0) * w[:, None]).T @ (B - b0))
        R = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
        d = np.linalg.norm((A - a0) @ R - (B - b0), axis=1)
        w = 1.0 / (1.0 + (d / d0) ** 2)
        best = max(best, float(w.sum() / L))
    return best


def one(k):
    from rnafold3d import templates as T
    lib = T.load_library()
    m = lib.meta[k]
    q = lib.chain_seq(k)
    if set(q) - set('ACGU'):
        return None
    X = lib.chain_xyz(k)
    allowed = lib.allowed(m['release'], (m['pdb'],))
    scan = np.where(allowed, T.sw_scan(q, lib), -1)
    profile = T.pair_profile(q)
    out = {'query': f"{m['pdb']}_{m['chain']}", 'length': len(q), 'release': m['release'], 'title': m['title'][:60]}
    for label, structure in (('sequence', False), ('structure', True)):
        hits = T.search(q, lib, cutoff=m['release'], exclude=(m['pdb'],), max_hits=8, scan=scan,
                        structure=structure, profile=profile)
        hits = T.rank_templates(q, hits, lib, n=len(hits)) if hits else []
        rows = []
        for h in hits:
            Y = lib.chain_xyz(h['chain_index'])
            aln = [(a, b) for a, b in h['aligned'] if not np.isnan(X[a, 0])]
            tm = tm_aligned(np.array([X[a] for a, _ in aln]), np.array([Y[b] for _, b in aln]), len(q)) if aln else 0.0
            rows.append({**{f: h.get(f) for f in FEATURES}, 'tm': round(tm, 4)})
        out[label] = rows
    return out


def cmd_run(a):
    import multiprocessing as mp
    from rnafold3d.templates import load_library
    lib = load_library()
    ks = [k for k, m in enumerate(lib.meta) if (m['release'] or '') >= a.since and a.min_len <= m['length'] <= a.max_len
          and np.isnan(lib.chain_xyz(k)[:, 0]).mean() < 0.2]
    if a.limit:
        ks = ks[::max(1, len(ks) // a.limit)]
    print(f'{len(ks)} queries', flush=True)
    with mp.get_context('fork').Pool(a.workers) as pool:
        rows = [r for r in pool.imap(one, ks, chunksize=4) if r]
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps({'since': a.since, 'library': lib.version, 'queries': rows}))
    print('wrote', a.out)
    summary(rows)


def summary(rows, key=None):
    """Mean value of the first-ranked template and of the best hit; queries without a hit count as 0."""
    n = len(rows)
    print(f"{'search':10s} {'with hit':>9s} {'top-1 TM':>9s} {'best TM':>8s} {'top-1 TM (hit)':>15s} {'top-1 >= 0.45':>14s}")
    for label in ('sequence', 'structure'):
        top, best, hit = [], [], 0
        for r in rows:
            h = r[label]
            if key and h:
                h = sorted(h, key=key)
            hit += bool(h)
            top.append(h[0]['tm'] if h else 0.0)
            best.append(max(x['tm'] for x in h) if h else 0.0)
        wh = [t for t, r in zip(top, rows) if r[label]]
        print(f"{label:10s} {hit:5d}/{n:<4d} {np.mean(top):9.3f} {np.mean(best):8.3f} {np.mean(wh) if wh else 0:15.3f} "
              f"{sum(t >= 0.45 for t in top):9d}")


def cmd_report(a):
    for f in a.files:
        d = json.loads(Path(f).read_text())
        print(f, d['library'], 'since', d['since'])
        summary(d['queries'])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run')
    r.add_argument('--since', default='2021-01-01')
    r.add_argument('--min-len', type=int, default=40)
    r.add_argument('--max-len', type=int, default=250)
    r.add_argument('--limit', type=int, default=0, help='use about this many queries (evenly spaced)')
    r.add_argument('--workers', type=int, default=14)
    r.add_argument('--out', default=str(ROOT / 'benchmark' / 'runs' / 'template_selection' / 'run.json'))
    r.set_defaults(fn=cmd_run)
    p = sub.add_parser('report')
    p.add_argument('files', nargs='+')
    p.set_defaults(fn=cmd_report)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
