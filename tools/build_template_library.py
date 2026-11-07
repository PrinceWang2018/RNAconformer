"""Build the RNA structure library used for template-based modelling and junction geometry.

Source: the BGSU RNA 3D Hub non-redundant list ("all" resolutions, one equivalence class per distinct
RNA chain; Leontis & Zirbel 2012) and coordinates from the RCSB PDB ModelServer (only the listed chain).

Per class the library keeps the representative chain and, when the representative was released after
`--split-date` (default 2024-05-01, before the CASP16 prediction season), also the best-ranked member
released before that date, so that benchmarks can restrict templates to structures that existed at the
time (template search filters by release date).

Stored per residue (compact, ~32 bytes before compression): sequence letter, C4' coordinates (int16,
0.05 A units relative to the chain centroid), P, C1', the glycosidic N and C2/C4/C6 as offsets from the
residue's C4' (int8, 0.1 A), backbone/glycosidic torsions (int8, 1.4 degree) and the canonical base pairs
(cis Watson-Crick and G-U wobble, rnafold3d.annotate).

    python tools/build_template_library.py --cache ~/.cache/rcv_templates --out src/rnafold3d/data/pdb_rna_v2.npz
    python tools/build_template_library.py --repack old_schema1.npz --out src/rnafold3d/data/pdb_rna_v2.npz
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime
import gzip
import io
import json
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d.structure_io import parse_cif, _mmcif_model          # noqa: E402
from rnafold3d.annotate import _pair_score                          # noqa: E402
from rnafold3d.geomutil import dihedral                              # noqa: E402

NR_URL = 'http://rna.bgsu.edu/rna3dhub/nrlist/download/current/all/csv'
NR_RELEASE_URL = 'http://rna.bgsu.edu/rna3dhub/nrlist/release/current'
MODEL_SERVER = 'https://models.rcsb.org/v1/{pdb}/atoms?auth_asym_id={chain}&model_nums={model}&encoding=cif'
GRAPHQL = 'https://data.rcsb.org/graphql'
ATOMS = ('P', "C4'", "C1'", 'NG', 'C2', 'C4', 'C6')       # NG = N9 (purines) / N1 (pyrimidines)
MISSING = -32768
SCALE = 20.0                                             # int16 units per angstrom (0.05 A), C4' atoms
REL_SCALE = 10.0                                         # int8 units per angstrom (0.1 A), atoms relative to C4'
TOR_SCALE = 127.0 / 180.0                                # int8 units per degree (1.4 deg)
MIN_LEN = 15


def fetch(url, data=None, tries=4, timeout=120):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'Accept-Encoding': 'gzip', 'User-Agent': 'rcv-template-builder',
                                                                  **({'Content-Type': 'application/json'} if data else {})})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                raw = r.read()
                if r.headers.get('Content-Encoding') == 'gzip':
                    raw = gzip.decompress(raw)
                return raw
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(2 + 3 * k)


def parse_nr(text):
    classes = []
    for row in csv.reader(io.StringIO(text)):
        if len(row) < 3:
            continue
        members = [m for m in row[2].split(',') if m]
        classes.append({'id': row[0], 'rep': row[1], 'members': members})
    return classes


def release_dates(pdb_ids, cache):
    path = cache / 'release_dates.json'
    known = json.loads(path.read_text()) if path.exists() else {}
    todo = sorted(set(pdb_ids) - set(known))
    for k in range(0, len(todo), 300):
        chunk = todo[k:k + 300]
        q = '{entries(entry_ids:%s){rcsb_id rcsb_accession_info{initial_release_date} exptl{method} ' \
            'rcsb_entry_info{resolution_combined} struct{title}}}' % json.dumps(chunk)
        d = json.loads(fetch(GRAPHQL, json.dumps({'query': q}).encode()))
        for e in d['data']['entries'] or []:
            res = (e.get('rcsb_entry_info') or {}).get('resolution_combined') or [None]
            known[e['rcsb_id']] = {'release': (e['rcsb_accession_info']['initial_release_date'] or '')[:10],
                                   'method': ','.join(x['method'] for x in e.get('exptl') or []),
                                   'resolution': res[0], 'title': (e.get('struct') or {}).get('title') or ''}
        for pid in chunk:
            known.setdefault(pid, {'release': '', 'method': '', 'resolution': None, 'title': ''})
        print(f'  release dates {min(k + 300, len(todo))}/{len(todo)}', flush=True)
        path.write_text(json.dumps(known))
    return known


def chain_text(pdb, model, chain, cache):
    p = cache / 'chains' / f'{pdb}_{model}_{chain}.cif.gz'
    if p.exists():
        return gzip.decompress(p.read_bytes()).decode(errors='replace')
    raw = fetch(MODEL_SERVER.format(pdb=pdb, chain=chain, model=model))
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(gzip.compress(raw))
    return raw.decode(errors='replace')


class _R:
    __slots__ = ('index', 'base', 'atoms')

    def __init__(self, index, base, atoms):
        self.index, self.base, self.atoms = index, base, atoms


def canonical_pairs_fast(res):
    """Canonical pairs (cWW A-U/G-C, G-U wobble) with a C1' neighbour search; residues: list of _R."""
    idx = [k for k, r in enumerate(res) if "C1'" in r.atoms and r.base in 'ACGU']
    if len(idx) < 2:
        return []
    pts = np.array([res[k].atoms["C1'"] for k in idx])
    cands = []
    for a, b in cKDTree(pts).query_pairs(12.5):
        ri, rj = res[idx[a]], res[idx[b]]
        if abs(ri.index - rj.index) < 4:
            continue
        ok, typ, score = _pair_score(ri, rj)
        if ok:
            cands.append((score, ri.index, rj.index))
    cands.sort()
    used, pairs = set(), []
    for _, i, j in cands:
        if i in used or j in used:
            continue
        used.update((i, j))
        pairs.append((min(i, j), max(i, j)))
    return sorted(pairs)


def torsions(prev, r, nxt):
    """alpha, beta, gamma, delta, epsilon, zeta, chi of residue r (degrees or None)."""
    a = r.atoms
    g = lambda d, k: d.atoms.get(k) if d is not None else None

    def dh(*pts):
        return None if any(p is None for p in pts) else float(dihedral(*pts))
    linked_prev = prev is not None and g(prev, "O3'") is not None and a.get('P') is not None and \
        np.linalg.norm(g(prev, "O3'") - a['P']) < 2.0
    linked_next = nxt is not None and a.get("O3'") is not None and g(nxt, 'P') is not None and \
        np.linalg.norm(a["O3'"] - g(nxt, 'P')) < 2.0
    alpha = dh(g(prev, "O3'"), a.get('P'), a.get("O5'"), a.get("C5'")) if linked_prev else None
    beta = dh(a.get('P'), a.get("O5'"), a.get("C5'"), a.get("C4'"))
    gamma = dh(a.get("O5'"), a.get("C5'"), a.get("C4'"), a.get("C3'"))
    delta = dh(a.get("C5'"), a.get("C4'"), a.get("C3'"), a.get("O3'"))
    eps = dh(a.get("C4'"), a.get("C3'"), a.get("O3'"), g(nxt, 'P')) if linked_next else None
    zeta = dh(a.get("C3'"), a.get("O3'"), g(nxt, 'P'), g(nxt, "O5'")) if linked_next else None
    if r.base in 'AG':
        chi = dh(a.get("O4'"), a.get("C1'"), a.get('N9'), a.get('C4'))
    else:
        chi = dh(a.get("O4'"), a.get("C1'"), a.get('N1'), a.get('C2'))
    return [alpha, beta, gamma, delta, eps, zeta, chi]


def process_chain(pdb, model, chain, cache):
    text = chain_text(pdb, model, chain, cache)
    cif = parse_cif(text)
    chains, _ = _mmcif_model(cif, rna_only=True, model_id=None, label=f'{pdb}.{chain}')
    ch = next((c for c in chains.values() if c.auth_asym_id == chain), None)
    if ch is None or not ch.residues:
        return None
    seq = ''.join(b if b in 'ACGU' else 'N' for b in ch.sequence)
    # residue index = label_seq_id (1-based position in the full polymer sequence)
    if not seq or max(r.index for r in ch.residues) > len(seq):
        seq = ''.join(r.base for r in ch.residues)
        for k, r in enumerate(ch.residues, 1):
            r.index = k
    n = len(seq)
    rmap = {r.index: r for r in ch.residues}
    xyz = np.full((n, len(ATOMS), 3), np.nan, np.float32)
    tor = np.full((n, 7), np.nan, np.float32)
    for k in range(1, n + 1):
        r = rmap.get(k)
        if r is None:
            continue
        for a, nm in enumerate(ATOMS):
            if nm == 'NG':
                nm = 'N9' if r.base in 'AG' else 'N1'
            v = r.atoms.get(nm)
            if v is not None:
                xyz[k - 1, a] = v
        t = torsions(rmap.get(k - 1), r, rmap.get(k + 1))
        tor[k - 1] = [np.nan if v is None else v for v in t]
    pairs = canonical_pairs_fast([_R(r.index, r.base, r.atoms) for r in ch.residues])
    pairs = [(i - 1, j - 1) for i, j in pairs]
    observed = int(np.sum(~np.isnan(xyz[:, 1, 0])))
    return {'seq': seq, 'xyz': xyz, 'tor': tor, 'pairs': pairs, 'observed': observed}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cache')
    ap.add_argument('--out', default=str(ROOT / 'src' / 'rnafold3d' / 'data' / 'pdb_rna_v2.npz'))
    ap.add_argument('--repack', help='convert an existing schema-1 library file instead of building')
    ap.add_argument('--split-date', default='2024-05-01')
    ap.add_argument('--max-length', type=int, default=5000)
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--limit', type=int, default=0, help='first N classes only (testing)')
    a = ap.parse_args()
    if a.repack:
        return repack(a.repack, a.out)
    if not a.cache:
        ap.error('--cache is required to build')
    cache = Path(a.cache).expanduser()
    cache.mkdir(parents=True, exist_ok=True)
    nr_path = cache / 'nr_all.csv'
    if not nr_path.exists():
        nr_path.write_bytes(fetch(NR_URL))
    nr_release = ''
    try:
        page = fetch(NR_RELEASE_URL).decode(errors='replace')
        import re
        m = re.search(r'Release\s+([0-9.]+)', page)
        nr_release = m.group(1) if m else ''
    except Exception:
        pass
    classes = parse_nr(nr_path.read_text())
    if a.limit:
        classes = classes[:a.limit]
    print(f'{len(classes)} equivalence classes (BGSU NR all, release {nr_release or "current"})')
    ids = {m.split('|')[0] for c in classes for ife in c['members'] for m in ife.split('+')}
    dates = release_dates(ids, cache)
    jobs = []
    for c in classes:
        rank = {}
        chosen = [(0, c['members'][0])]
        rep_date = dates.get(c['members'][0].split('|')[0], {}).get('release', '')
        if not rep_date or rep_date >= a.split_date:
            for k, ife in enumerate(c['members'][1:], 1):
                d = dates.get(ife.split('|')[0], {}).get('release', '')
                if d and d < a.split_date:
                    chosen.append((k, ife))
                    break
        for k, ife in chosen:
            for part in ife.split('+'):
                pdb, model, chain = part.split('|')
                key = (pdb, model, chain)
                if key in rank:
                    continue
                rank[key] = k
                jobs.append({'class': c['id'], 'rank': k, 'pdb': pdb, 'model': model, 'chain': chain,
                             **dates.get(pdb, {'release': '', 'method': '', 'resolution': None, 'title': ''})})
    print(f'{len(jobs)} chains to process')
    results = [None] * len(jobs)

    def work(k):
        j = jobs[k]
        try:
            return k, process_chain(j['pdb'], j['model'], j['chain'], cache), None
        except Exception as exc:
            return k, None, f'{type(exc).__name__}: {exc}'
    errors = 0
    with cf.ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(work, k) for k in range(len(jobs))]
        for n_done, f in enumerate(cf.as_completed(futs), 1):
            k, res, err = f.result()
            if err:
                errors += 1
                print('  error', jobs[k]['pdb'], jobs[k]['chain'], err, flush=True)
            results[k] = res
            if n_done % 200 == 0:
                print(f'  {n_done}/{len(jobs)} chains', flush=True)
    # assemble
    seqs, xyzs, tors, pairs, meta, offsets, origins = [], [], [], [], [], [0], []
    for j, r in zip(jobs, results):
        if r is None or r['observed'] < MIN_LEN or len(r['seq']) > a.max_length:
            continue
        xyz = r['xyz']
        ok = ~np.isnan(xyz[:, :, 0])
        origin = np.nanmean(xyz[:, 1], axis=0)
        q = np.where(ok[:, :, None], np.round((xyz - origin) * SCALE), MISSING)
        if np.nanmax(np.abs(np.where(ok[:, :, None], q, 0))) >= 32767:
            continue                                     # chain too large for the int16 grid
        t = np.where(np.isnan(r['tor']), MISSING, np.round(r['tor'] * 100))
        base = offsets[-1]
        seqs.append(r['seq'])
        xyzs.append(q.astype(np.int16))
        tors.append(t.astype(np.int16))
        pairs += [(base + i, base + jj) for i, jj in r['pairs']]
        origins.append(origin.astype(np.float32))
        offsets.append(base + len(r['seq']))
        meta.append({'pdb': j['pdb'], 'model': j['model'], 'chain': j['chain'], 'class': j['class'], 'rank': j['rank'],
                     'release': j['release'], 'method': j['method'], 'resolution': j['resolution'],
                     'title': j['title'][:120], 'length': len(r['seq']), 'observed': r['observed'],
                     'n_pairs': len(r['pairs'])})
    info = {'schema': 'rnafold3d.templates/1', 'created': datetime.date.today().isoformat(),
            'source': f'BGSU RNA 3D Hub non-redundant list (all resolutions) {nr_release}; coordinates RCSB PDB ModelServer',
            'split_date': a.split_date, 'atoms': list(ATOMS), 'scale_per_angstrom': SCALE, 'missing': MISSING,
            'torsions': ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi'], 'min_observed': MIN_LEN,
            'chains': meta}
    out = Path(a.out)
    save(out, ''.join(seqs), offsets, np.concatenate(xyzs), np.concatenate(tors), pairs, origins, info)
    print(f'{len(meta)} chains, {offsets[-1]} residues, {len(pairs)} pairs, {errors} download/parse errors -> '
          f'{out} ({out.stat().st_size / 1e6:.1f} MB)')


def save(out, seq, offsets, xyz_q, tor_q, pairs, origins, info):
    """Compact storage (schema 2): C4' on the int16 grid of the chain; the other atoms as int8 offsets from
    their residue's C4' (0.1 A); torsions as int8 (1.4 deg). About 24 bytes per residue before compression."""
    c4 = xyz_q[:, 1].astype(np.int16)
    ok = xyz_q[:, :, 0] != MISSING
    rel = np.round((xyz_q.astype(np.float64) - c4[:, None, :]) / SCALE * REL_SCALE)
    rel_ok = ok & ok[:, 1:2] & np.all(np.abs(rel) <= 127, axis=2)
    rel = np.where(rel_ok[:, :, None], rel, -128).astype(np.int8)
    rel = np.delete(rel, 1, axis=1)                      # C4' itself is the reference
    t_ok = tor_q != MISSING
    tor8 = np.where(t_ok, np.clip(np.round(tor_q / 100.0 * TOR_SCALE), -127, 127), -128).astype(np.int8)
    info = {**info, 'schema': 'rnafold3d.templates/2', 'rel_scale_per_angstrom': REL_SCALE,
            'tor_scale_per_degree': TOR_SCALE}
    np.savez_compressed(out, seq=np.frombuffer(seq.encode(), np.uint8), offsets=np.array(offsets, np.int64),
                        c4=c4, rel=rel, tor8=tor8, pairs=np.array(pairs, np.int32).reshape(-1, 2),
                        origin=np.array(origins, np.float32), info=np.frombuffer(json.dumps(info).encode(), np.uint8))


def repack(src, out):
    """Convert a schema-1 library (int16 coordinates of all atoms) to the compact schema 2."""
    z = np.load(src)
    info = json.loads(bytes(z['info']).decode())
    p = z['pairs']
    save(Path(out), bytes(z['seq']).decode(), z['offsets'], z['xyz'], z['tor'], [tuple(x) for x in p], z['origin'], info)
    print(f'{src} -> {out} ({Path(out).stat().st_size / 1e6:.1f} MB)')


if __name__ == '__main__':
    main()
