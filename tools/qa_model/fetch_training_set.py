"""Training structures for the RCV confidence model (tools/qa_model), disjoint from the benchmark.

Source: the BGSU RNA 3D Hub non-redundant list (Leontis & Zirbel; resolution cut-off 4.0 A), one representative
chain per equivalence class. Kept: single-chain representatives, polyribonucleotide, 15-120 nt, canonical
sequence A/C/G/U only, X-ray/EM <= 3.5 A or NMR.

Excluded, so that the benchmark (benchmark/manifest_v1.json) stays a fair test:
  1. every equivalence class that contains a benchmark entry;
  2. every chain with >= 80 % sequence identity to a benchmark sequence (semi-global alignment, identities over the
     shorter sequence; short RNAs reach ~70 % against unrelated long sequences by chance, so related families are
     caught by rule 3);
  3. every chain whose entry title names a benchmark family (keyword list below).
Each chain is downloaded alone from the RCSB ModelServer (first model only).

Usage:
    python tools/qa_model/fetch_training_set.py --out benchmark/data_qa
"""
from __future__ import annotations

import argparse
import csv
import datetime
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NR_URL = 'http://rna.bgsu.edu/rna3dhub/nrlist/download/current/4.0A/csv'
GRAPHQL = 'https://data.rcsb.org/graphql'
MODEL_URL = 'https://models.rcsb.org/v1/{pdb}/atoms?auth_asym_id={chain}&model_nums=1&encoding=cif'

# benchmark families (manifest_v1.json), in the words entry titles use
FAMILY_WORDS = (
    'uucg', 'sarcin', 'ricin', 'tar ', 'tar-', 'tar rna', 'trans-activation', 'mmtv', 'mouse mammary', 'malachite',
    'guanidine', 'fluoride', 'sam-', 'sam ', 's-adenosyl', 'adomet', 'zmp', 'ztp', 'group ii', 'group-ii',
    'adenine', 'purine', 'guanine', 'xpt', 'pbux', '5-hydroxytryptophan', '5ht', 'xrrna', 'zika', 'flavivir',
    'exoribonuclease', 'l1 ligase', 'ligase ribozyme', 'trna', 'transfer rna', 'tpp', 'thiamine', 'thi-box', 'thim',
    'c-di-gmp', 'cyclic di-gmp', 'cyclic-di-gmp', 'c-di-amp', 'cyclic di-amp', 'ydao', 'group i', 'group-i', 'p4-p6',
    'p4p6', 'tetrahymena', 'rnase p', 'ribonuclease p', 'lysine', 'm-box', 'mgta', 'ykok', 'cobalamin', 'adocbl',
    'b12', 'lariat', 'capping ribozyme')
IDENTITY_MAX = 0.8


def fetch(url, data=None, timeout=60, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'} if data else {})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 + 3 * k)


def identity(a, b):
    """Identities / length of the shorter sequence in a semi-global alignment (free end gaps)."""
    if len(a) > len(b):
        a, b = b, a
    n, m = len(a), len(b)
    NEG = -10 ** 9
    # score, identities; match +1, mismatch -1, gap -2; gaps at the ends of b are free
    prev = [(0, 0)] * (m + 1)
    for i in range(1, n + 1):
        cur = [(-2 * i, 0)] + [(NEG, 0)] * m
        for j in range(1, m + 1):
            same = a[i - 1] == b[j - 1]
            d = (prev[j - 1][0] + (1 if same else -1), prev[j - 1][1] + same)
            u = (prev[j][0] - 2, prev[j][1])
            left = (cur[j - 1][0] - 2, cur[j - 1][1])
            cur[j] = max(d, u, left)
        prev = cur
    return max(prev)[1] / max(1, n)


def entry_info(ids):
    """RCSB data API: per entry the title, method, resolution and per chain (auth id) type, length, sequence."""
    q = '''query($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id struct { title }
      rcsb_entry_info { resolution_combined experimental_method }
      polymer_entities { entity_poly { type rcsb_sample_sequence_length pdbx_seq_one_letter_code_can }
        rcsb_polymer_entity_container_identifiers { auth_asym_ids } } } }'''
    out = {}
    for k in range(0, len(ids), 80):
        body = json.dumps({'query': q, 'variables': {'ids': ids[k:k + 80]}}).encode()
        d = json.loads(fetch(GRAPHQL, body))
        for e in (d.get('data') or {}).get('entries') or []:
            info = e.get('rcsb_entry_info') or {}
            res = info.get('resolution_combined') or [None]
            chains = {}
            for p in e.get('polymer_entities') or []:
                ep = p.get('entity_poly') or {}
                for ch in (p.get('rcsb_polymer_entity_container_identifiers') or {}).get('auth_asym_ids') or []:
                    chains[ch] = {'type': ep.get('type'), 'length': ep.get('rcsb_sample_sequence_length'),
                                  'sequence': (ep.get('pdbx_seq_one_letter_code_can') or '').replace('\n', '')}
            out[e['rcsb_id'].upper()] = {'title': ((e.get('struct') or {}).get('title') or ''), 'method': info.get('experimental_method'),
                                         'resolution': res[0], 'chains': chains}
        print(f'  entry info {min(k + 80, len(ids))}/{len(ids)}', flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=str(ROOT / 'benchmark' / 'data_qa'))
    ap.add_argument('--min-length', type=int, default=15)
    ap.add_argument('--max-length', type=int, default=120)
    ap.add_argument('--max-resolution', type=float, default=3.5)
    a = ap.parse_args()
    out = Path(a.out)
    (out / 'raw').mkdir(parents=True, exist_ok=True)
    nr_path = out / 'nrlist_4.0A.csv'
    if not nr_path.exists():
        nr_path.write_bytes(fetch(NR_URL))
    classes = list(csv.reader(io.StringIO(nr_path.read_text())))
    bench = json.loads((ROOT / 'benchmark' / 'manifest_v1.json').read_text())
    bench_ids = {t['pdb_id'].upper() for t in bench['targets']}
    ds = json.loads((ROOT / 'benchmark' / 'data' / 'dataset_v1.json').read_text())
    bench_seqs = []
    for s in ds['samples']:
        p = ROOT / 'benchmark' / 'data' / 'samples' / s['sample_id'] / 'sequence.fasta'
        if p.exists():
            bench_seqs.append((s['sample_id'], p.read_text().split('\n')[1].strip()))
    reps, skipped = [], {}
    for row in classes:
        if len(row) < 3:
            continue
        cls, rep, members = row[0], row[1], row[2]
        if '+' in rep:
            skipped[cls] = 'multi-chain representative'
            continue
        member_ids = {m.split('|')[0].upper() for m in members.replace('+', ',').split(',') if m}
        if member_ids & bench_ids:
            skipped[cls] = 'class contains a benchmark entry'
            continue
        pdb, _model, chain = rep.split('|')
        reps.append((cls, pdb.upper(), chain))
    cache = out / 'entry_info.json'
    info = json.loads(cache.read_text()) if cache.exists() else {}
    missing = sorted({p for _, p, _ in reps} - set(info))
    if missing:
        info.update(entry_info(missing))
        cache.write_text(json.dumps(info))
    kept = []
    for cls, pdb, chain in reps:
        e = info.get(pdb)
        c = e and e['chains'].get(chain)
        why = None
        if not c:
            why = 'no entry information'
        elif c['type'] != 'polyribonucleotide':
            why = f"type {c['type']}"
        elif not (a.min_length <= (c['length'] or 0) <= a.max_length):
            why = f"length {c['length']}"
        elif set(c['sequence']) - set('ACGU'):
            why = 'non-ACGU canonical sequence'
        elif e['method'] != 'NMR' and (e['resolution'] is None or e['resolution'] > a.max_resolution):
            why = f"resolution {e['resolution']}"
        else:
            title = ' ' + e['title'].lower() + ' '
            word = next((w for w in FAMILY_WORDS if w in title), None)
            if word:
                why = f'benchmark family keyword "{word.strip()}"'
            else:
                best = max(((identity(c['sequence'], s), sid) for sid, s in bench_seqs), default=(0, ''))
                if best[0] >= IDENTITY_MAX:
                    why = f'{best[0]:.0%} identity to benchmark {best[1]}'
        if why:
            skipped[cls] = why
            continue
        kept.append({'class': cls, 'pdb_id': pdb, 'chain': chain, 'length': c['length'], 'sequence': c['sequence'],
                     'title': e['title'], 'method': e['method'], 'resolution': e['resolution']})
    # one chain per distinct sequence
    seen, unique = set(), []
    for k in sorted(kept, key=lambda r: (r['method'] == 'NMR', r['resolution'] or 9.9)):
        if k['sequence'] in seen:
            continue
        seen.add(k['sequence'])
        unique.append(k)
    print(f'{len(classes)} classes, {len(unique)} chains kept', flush=True)
    for i, k in enumerate(unique):
        f = out / 'raw' / f"{k['pdb_id']}_{k['chain']}.cif"
        if not f.exists():
            try:
                f.write_bytes(fetch(MODEL_URL.format(pdb=k['pdb_id'], chain=k['chain'])))
            except Exception as exc:
                k['download_error'] = str(exc)
        if (i + 1) % 25 == 0:
            print(f'  downloaded {i + 1}/{len(unique)}', flush=True)
    reasons = {}
    for why in skipped.values():
        key = why.split(' ')[0] if why.startswith(('length', 'type', 'resolution')) else why.split(' to benchmark')[0] if 'identity' in why else why
        reasons[key] = reasons.get(key, 0) + 1
    (out / 'training_chains.json').write_text(json.dumps({
        'created': datetime.date.today().isoformat(), 'source': NR_URL, 'filters': {
            'length': [a.min_length, a.max_length], 'max_resolution': a.max_resolution, 'identity_max': IDENTITY_MAX,
            'family_words': FAMILY_WORDS, 'benchmark_ids': sorted(bench_ids)},
        'skipped_summary': dict(sorted(reasons.items(), key=lambda x: -x[1])), 'chains': unique}, indent=1))
    print(json.dumps(dict(sorted(reasons.items(), key=lambda x: -x[1])[:15]), indent=1))


if __name__ == '__main__':
    sys.exit(main())
