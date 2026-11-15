"""Benchmark 2 of the compound screen: does the score put known binders above look-alikes?

A redocking benchmark (`redock.py`) asks whether the search can find the right pose. It cannot tell you
whether the *ranking* is any good, because every compound in it is a binder. The question a screen actually
has to answer is the other one: given one target and a pile of compounds, do the ones that bind come out on
top? That is this benchmark, and under the rule of the plan (`docs/NEXT_STEPS.md` §1.4) it is **the only
benchmark the weights of the scoring function may ever be tuned against** — never the selectivity cases of
benchmark 3, which are too few to fit anything and are the thing being claimed.

    python benchmark/screening/enrichment.py fetch                     # actives and decoys from PubChem
    python benchmark/screening/enrichment.py run --jobs 12 --out runs/enrichment.json
    python benchmark/screening/enrichment.py report runs/enrichment.json
    python benchmark/screening/enrichment.py tune runs/enrichment.json  # weights, cross-validated

How it is built, and what each choice costs:

* **Six targets**, each an experimental RNA structure with its own pocket: the FMN, TPP, adenine, guanine
  and preQ1 riboswitches and the bacterial ribosomal A site. The pocket is the box around where the
  experimental ligand sat, with the ligand removed — the same easy case as benchmark 1, and for the same
  reason (a predicted pocket adds a second, unmeasured error to every number here).
* **Actives are hand-collected from the literature**, listed below with what they are. They are compounds
  reported to bind that RNA, not a curated public set: there is no RNA equivalent of DUD-E, and building
  one is a project of its own. Between four and twelve per target is few enough that one wrong assignment
  moves the number, so each one is named rather than hidden behind a file.
* **Decoys are property-matched from PubChem**: a random sample of the database, matched to each active on
  molecular weight, lipophilicity, polar surface area, rotatable bonds and hydrogen-bond counts, with
  everything 2D-similar to any active of that target removed. Property matching is what stops the score
  from winning on molecular weight alone. Twenty-five decoys per active, not the 50 to 100 the plan asks
  for: the search costs seconds per compound, and the ratio is what the run time is spent on. What
  property matching does not give is *experimentally inactive* compounds —
  a matched random compound is assumed inactive, which is the standard assumption of this kind of
  benchmark and is wrong for a few percent of them.
* Every compound is docked with the same settings, and the raw value of every term is stored for every
  kept pose. That is what makes tuning possible without redocking: a change of weights is a dot product
  over the stored terms, so a weight search costs seconds instead of days. Its limit is that it can only
  explore weights whose poses the original run already found — which is why `tune` reports how far it
  moved and the final weights have to be confirmed by a fresh run.

A random ranking scores AUC 0.5. The plan's line is that below about 0.65 the screen should not be offered
at all.

Compounds are downloaded from PubChem (3D conformers, public domain) and structures from RCSB, into
~/.rcv/benchmarks/enrichment; nothing is added to the repository.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from redock import CACHE as PDB_CACHE, URL as PDB_URL                  # noqa: E402

CACHE = Path.home() / '.rcv' / 'benchmarks' / 'enrichment'
PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
PROPS = ('MolecularWeight', 'XLogP', 'TPSA', 'RotatableBondCount', 'HeavyAtomCount',
         'HBondDonorCount', 'HBondAcceptorCount')
MATCH = ('MolecularWeight', 'XLogP', 'TPSA', 'RotatableBondCount', 'HBondDonorCount', 'HBondAcceptorCount')
CID_MAX = 160_000_000          # PubChem's compound ids run to about this; the sample is uniform over it
POOL = 24_000                  # background rows to sample; an aminoglycoside needs a big pool to match
SIMILARITY = 85                # a decoy 2D-similar to any active of the target at this % is thrown out

# One entry per target: the structure, the residue name of the experimental ligand that defines the pocket,
# and the actives with what each one is. Everything here is the literature, not a computation.
TARGETS = [
    dict(key='fmn', name='FMN riboswitch', pdb='3F2Q', ligand='FMN',
         what='F. nucleatum FMN riboswitch, 2.95 Å; the target of ribocil',
         actives=[('FMN', 'the natural ligand'),
                  ('riboflavin', 'the unphosphorylated vitamin, binds more weakly'),
                  ('roseoflavin', 'a natural flavin analogue antimicrobial'),
                  ('ribocil', 'a synthetic non-flavin binder found by a cell-based screen')]),
    dict(key='tpp', name='TPP riboswitch', pdb='2GDI', ligand='TPP',
         what='E. coli thiM TPP riboswitch, 2.05 Å',
         actives=[('thiamine pyrophosphate', 'the natural ligand'),
                  ('thiamine monophosphate', 'one phosphate short, binds weakly'),
                  ('thiamine', 'the vitamin itself, weaker again'),
                  ('pyrithiamine', 'an antimetabolite; its pyrophosphate binds the riboswitch'),
                  ('oxythiamine', 'an antimetabolite analogue')]),
    dict(key='adenine', name='adenine riboswitch', pdb='1Y26', ligand='ADE',
         what='V. vulnificus add adenine riboswitch, 2.1 Å',
         actives=[('adenine', 'the natural ligand'),
                  ('2-aminopurine', 'binds the purine riboswitch family'),
                  ('2,6-diaminopurine', 'binds both the adenine and the guanine switch'),
                  ('2-fluoroadenine', 'an adenine analogue used in the specificity studies')]),
    dict(key='guanine', name='guanine riboswitch', pdb='1U8D', ligand='HPA',
         what='B. subtilis xpt-pbuX guanine riboswitch, 1.95 Å',
         actives=[('guanine', 'the natural ligand'),
                  ('hypoxanthine', 'the ligand of this entry'),
                  ('xanthine', 'binds weakly'),
                  ('2,6-diaminopurine', 'binds both purine switches'),
                  ('6-chloroguanine', 'a guanine analogue from the specificity studies'),
                  ('7-deazaguanine', 'guanine without N7')]),
    dict(key='preq1', name='preQ1 riboswitch', pdb='3Q50', ligand='PRF',
         what='class II preQ1 riboswitch, 2.3 Å',
         actives=[('preQ1', 'the natural ligand, 7-aminomethyl-7-deazaguanine'),
                  ('preQ0', 'the nitrile precursor, binds well'),
                  ('7-deazaguanine', 'the shared core'),
                  ('guanine', 'binds weakly')]),
    dict(key='asite', name='ribosomal A site', pdb='1J7T', ligand='PAR',
         what='16S decoding-site model duplex with paromomycin, 2.5 Å; the aminoglycoside target',
         actives=[('paromomycin', 'the ligand of this entry'),
                  ('neomycin', 'neomycin B'),
                  ('ribostamycin', 'neomycin without ring D'),
                  ('neamine', 'the two-ring core'),
                  ('kanamycin A', '4,6-disubstituted'),
                  ('kanamycin B', '4,6-disubstituted'),
                  ('tobramycin', '3-deoxykanamycin B'),
                  ('apramycin', 'a bicyclic aminoglycoside'),
                  ('geneticin', 'G418'),
                  ('sisomicin', '4,6-disubstituted, unsaturated'),
                  ('netilmicin', 'N-ethylsisomicin'),
                  ('hygromycin B', 'binds the A site at an overlapping position')]),
]


# --------------------------------------------------------------------------------------- fetching

def _get(url, tries=4, pause=0.3):
    """One PubChem request, with the pauses its usage policy asks for (5 requests a second)."""
    for k in range(tries):
        try:
            time.sleep(pause)
            with urllib.request.urlopen(url, timeout=90) as r:
                return r.read()
        except Exception as err:                                      # noqa: BLE001 - retry, then give up
            if k == tries - 1:
                raise
            time.sleep(1.5 * (k + 1))
            del err
    return b''


def _properties(cids):
    """The property table for a list of CIDs, in batches."""
    out = {}
    cids = [int(c) for c in cids]
    for k in range(0, len(cids), 180):
        batch = cids[k:k + 180]
        url = f'{PUG}/compound/cid/{",".join(map(str, batch))}/property/{",".join(PROPS)}/JSON'
        try:
            rows = json.loads(_get(url))['PropertyTable']['Properties']
        except Exception:                                             # noqa: BLE001 - a batch may 404
            continue
        for row in rows:
            if 'MolecularWeight' in row and 'XLogP' in row:
                out[int(row['CID'])] = {p: float(row[p]) for p in PROPS if p in row}
    return out


def _resolve(name):
    """A compound name to its CID and properties, as PubChem resolves it."""
    url = f'{PUG}/compound/name/{urllib.parse.quote(name)}/property/{",".join(PROPS)}/JSON'
    rows = json.loads(_get(url))['PropertyTable']['Properties']
    row = rows[0]
    return int(row['CID']), {p: float(row[p]) for p in PROPS if p in row}


def _similar(cid, threshold=SIMILARITY, limit=4000):
    """The CIDs PubChem calls 2D-similar to this one: the compounds a decoy must not be."""
    url = (f'{PUG}/compound/fastsimilarity_2d/cid/{cid}/cids/JSON'
           f'?Threshold={threshold}&MaxRecords={limit}')
    try:
        return set(json.loads(_get(url, pause=0.6))['IdentifierList']['CID'])
    except Exception:                                                 # noqa: BLE001 - no neighbours is fine
        return set()


def _sdf(cids, path):
    """Download 3D conformers for these CIDs into one SDF. Compounds without a 3D record are skipped."""
    text = []
    for k in range(0, len(cids), 40):
        batch = ','.join(str(c) for c in cids[k:k + 40])
        try:
            data = _get(f'{PUG}/compound/cid/{batch}/SDF?record_type=3d').decode('utf-8', 'replace')
        except Exception:                                             # noqa: BLE001
            continue
        if '$$$$' in data:
            text.append(data)
    path.write_text(''.join(text), encoding='utf-8')
    return sum(t.count('$$$$') for t in text)


def _pdb(pdb):
    """The structure, from the redocking benchmark's cache or from RCSB."""
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    out = PDB_CACHE / f'{pdb}.pdb'
    if not out.exists():
        with urllib.request.urlopen(PDB_URL.format(pdb), timeout=90) as r:
            out.write_bytes(r.read())
    return out


def _pool(seed=0, size=POOL):
    """A random sample of PubChem with its properties, cached: the background the decoys come from."""
    path = CACHE / f'pool_{seed}_{size}.json'
    if path.exists():
        rows = json.loads(path.read_text(encoding='utf-8'))
        return {int(k): v for k, v in rows.items()}
    rng = np.random.default_rng(seed)
    out, tried = {}, 0
    while len(out) < size and tried < 6 * size:
        want = [int(c) for c in rng.integers(1, CID_MAX, 360)]
        tried += len(want)
        out.update(_properties(want))
        print(f'\r  background pool: {len(out)} compounds of {size} ({tried} ids tried)', end='', flush=True)
    print()
    CACHE.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out), encoding='utf-8')
    return out


def _match(active, pool, exclude, taken, ratio, scale):
    """The `ratio` compounds of the pool closest to this active in properties, and not similar to it."""
    cids, dist = [], []
    for cid, row in pool.items():
        if cid in exclude or cid in taken or any(p not in row for p in MATCH):
            continue
        d = sum(((row[p] - active[p]) / scale[p]) ** 2 for p in MATCH)
        cids.append(cid)
        dist.append(d)
    if not cids:
        return [], 0.0
    idx = np.argsort(dist)[:ratio]
    return [cids[i] for i in idx], float(np.sqrt(np.mean([dist[i] for i in idx])))


def fetch(args):
    CACHE.mkdir(parents=True, exist_ok=True)
    pool = _pool(args.pool_seed, args.pool)
    scale = {p: (np.std([r[p] for r in pool.values() if p in r]) or 1.0) for p in MATCH}
    manifest = {'similarity': SIMILARITY, 'ratio': args.ratio, 'pool': len(pool),
                'pool_seed': args.pool_seed, 'targets': []}
    for t in TARGETS:
        if args.only and t['key'] not in args.only.split(','):
            continue
        print(f'{t["name"]} ({t["pdb"]})')
        _pdb(t['pdb'])
        actives, exclude = [], set()
        for name, what in t['actives']:
            try:
                cid, prop = _resolve(name)
            except Exception as err:                                  # noqa: BLE001
                print(f'  {name}: NOT RESOLVED ({err})')
                continue
            actives.append({'name': name, 'cid': cid, 'what': what, 'properties': prop})
            print(f'  active  CID {cid:>12d}  MW {prop.get("MolecularWeight", 0):6.1f}  '
                  f'logP {prop.get("XLogP", float("nan")):5.1f}  {name}')
        for a in actives:                                             # no decoy may look like any active
            exclude |= _similar(a['cid'])
        exclude |= {a['cid'] for a in actives}
        decoys, taken, gaps = [], set(), []
        for a in actives:
            if any(p not in a['properties'] for p in MATCH):
                print(f'  {a["name"]}: no property row, no decoys matched to it')
                continue
            got, gap = _match(a['properties'], pool, exclude, taken, args.ratio, scale)
            taken |= set(got)
            gaps.append(gap)
            decoys += [{'cid': c, 'matched_to': a['name']} for c in got]
        print(f'  {len(actives)} actives, {len(decoys)} decoys, {len(exclude)} similar compounds excluded, '
              f'property gap {np.mean(gaps) if gaps else 0:.2f} standard deviations')
        n_a = _sdf([a['cid'] for a in actives], CACHE / f'{t["key"]}_actives.sdf')
        n_d = _sdf([d['cid'] for d in decoys], CACHE / f'{t["key"]}_decoys.sdf')
        print(f'  3D conformers: {n_a} of {len(actives)} actives, {n_d} of {len(decoys)} decoys')
        manifest['targets'].append(dict(t, actives=actives, decoys=decoys,
                                        property_gap=round(float(np.mean(gaps)) if gaps else 0.0, 3),
                                        n_actives_3d=n_a, n_decoys_3d=n_d))
    (CACHE / 'manifest.json').write_text(json.dumps(manifest, indent=1), encoding='utf-8')
    print(f'manifest written to {CACHE / "manifest.json"}')


# --------------------------------------------------------------------------------------- running

_CTX = {}


def _setup(target, exhaustiveness, seed, n_poses, box):
    """Receptor and grid for one target, kept in a module global so forked workers inherit them."""
    from rnafold3d.screening import grid
    from rnafold3d.screening.receptor import read_complex
    rec, ligands = read_complex(_pdb(target['pdb']))
    want = [m for m in ligands if target['ligand'] in m.properties.get('residues', [])]
    if not want:
        raise SystemExit(f'{target["pdb"]}: the ligand {target["ligand"]} is not in the entry')
    native = want[0]
    span = float((native.xyz.max(0) - native.xyz.min(0)).max())
    _CTX.clear()
    _CTX.update(rec=rec, grid=grid.Grid.for_pocket(rec, centre=native.xyz.mean(0), size=span + 2 * box),
                exhaustiveness=exhaustiveness, seed=seed, n_poses=n_poses, native=native)
    return _CTX


def _dock_one(job):
    """Dock one compound and keep the raw terms of every pose, which is what the tuning reads."""
    which, index, mol = job
    from rnafold3d.screening import dock, score
    rec = _CTX['rec']
    row = {'set': which, 'index': index, 'name': mol.name, 'atoms': len(mol),
           'rotatable': mol.n_rotatable, 'charge': round(float(mol.charge.sum()), 2),
           'properties': {k: v for k, v in mol.properties.items() if k in ('cid', 'matched_to', 'label')}}
    t0 = time.time()
    try:
        poses = dock.dock(mol, _CTX['grid'], receptor=rec, exhaustiveness=_CTX['exhaustiveness'],
                          seed=_CTX['seed'], n_poses=_CTX['n_poses'])
    except Exception as err:                                          # noqa: BLE001 - one bad compound
        row['error'] = f'{type(err).__name__}: {err}'
        return row
    if not poses:
        row['error'] = 'no pose'
        return row
    row['poses'] = []
    for p in poses:
        raw = score.pair_terms(p.xyz, mol, rec)
        exact = score.score(p.xyz, mol, rec)
        row['poses'].append({'score': round(p.score, 3),
                             'internal': round(p.score - exact, 3),
                             'terms': {k: round(v, 5) for k, v in raw.items()}})
    row['score'] = row['poses'][0]['score']
    row['seconds'] = round(time.time() - t0, 1)
    return row


def _compounds(target, which, limit=None):
    """The prepared compounds of one set, with their manifest labels attached."""
    from rnafold3d.screening import chem
    path = CACHE / f'{target["key"]}_{which}.sdf'
    if not path.exists():
        raise SystemExit(f'{path} is missing: run "fetch" first')
    by_cid = {}
    for row in target.get('actives' if which == 'actives' else 'decoys', []):
        by_cid.setdefault(int(row['cid']), row)
    out = []
    for mol in chem.read(path):
        if mol.properties.get('error') or not len(mol):
            continue
        try:
            cid = int(mol.name)
        except (TypeError, ValueError):
            cid = None
        row = by_cid.get(cid, {})
        mol.properties.update(cid=cid, label=row.get('name', str(cid)),
                              matched_to=row.get('matched_to', ''))
        mol.name = row.get('name', str(cid))
        out.append(mol)
        if limit and len(out) >= limit:
            break
    return out


def run(args):
    import multiprocessing as mp
    manifest = json.loads((CACHE / 'manifest.json').read_text(encoding='utf-8'))
    results, t_start = [], time.time()
    for target in manifest['targets']:
        if args.only and target['key'] not in args.only.split(','):
            continue
        ctx = _setup(target, args.exhaustiveness, args.seed, args.n_poses, args.box)
        jobs = [('actives', k, m) for k, m in enumerate(_compounds(target, 'actives', args.limit))]
        jobs += [('decoys', k, m) for k, m in enumerate(_compounds(target, 'decoys', args.limit))]
        n_act = sum(1 for j in jobs if j[0] == 'actives')
        print(f'{target["name"]}: {n_act} actives and {len(jobs) - n_act} decoys, '
              f'grid {ctx["grid"].shape} ({ctx["grid"].nbytes / 1e6:.0f} MB), '
              f'{args.exhaustiveness}, {args.jobs} workers', flush=True)
        t0 = time.time()
        if args.jobs > 1:
            with mp.get_context('fork').Pool(args.jobs) as pool:
                rows = []
                for r in pool.imap_unordered(_dock_one, jobs, chunksize=1):
                    rows.append(r)
                    print(f'\r  {len(rows)}/{len(jobs)} docked, {time.time() - t0:.0f} s',
                          end='', flush=True)
        else:
            rows = [_dock_one(j) for j in jobs]
        print()
        bad = [r for r in rows if 'error' in r]
        if bad:
            print(f'  {len(bad)} compounds could not be docked ' + ', '.join(
                sorted({r['error'].split(':')[0] for r in bad})))
        out = dict(target, results=rows, seconds=round(time.time() - t0, 1))
        results.append(out)
        _report_target(out)
    payload = {'benchmark': 'enrichment', 'exhaustiveness': args.exhaustiveness, 'seed': args.seed,
               'n_poses': args.n_poses, 'box': args.box, 'similarity': manifest['similarity'],
               'ratio': manifest['ratio'], 'pool': manifest['pool'],
               'weights': dict(_weights_used()), 'minutes': round((time.time() - t_start) / 60, 1),
               'targets': results}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=1), encoding='utf-8')
        print(f'written to {args.out}')
    _report(payload)
    return payload


def _weights_used():
    from rnafold3d.screening import score
    return dict(score.WEIGHTS)


# --------------------------------------------------------------------------------------- scoring

def _matrix(target, terms=None):
    """The stored terms of one target as arrays: one row per compound, pose-minimised later.

    Returns `(labels, is_active, raw, internal, rotatable, failed)` where `raw` is a list of (poses, terms)
    arrays — one per compound — so a new set of weights is a dot product away.
    """
    from rnafold3d.screening import score
    terms = tuple(terms or score.TERMS)
    labels, active, raw, internal, rot, failed = [], [], [], [], [], []
    for r in target['results']:
        labels.append(r.get('properties', {}).get('label') or r.get('name'))
        active.append(r['set'] == 'actives')
        rot.append(r.get('rotatable', 0))
        if 'poses' not in r:
            raw.append(np.zeros((0, len(terms))))
            internal.append(np.zeros(0))
            failed.append(True)
            continue
        raw.append(np.array([[p['terms'][t] for t in terms] for p in r['poses']], float))
        internal.append(np.array([p['internal'] for p in r['poses']], float))
        failed.append(False)
    return labels, np.array(active), raw, internal, np.array(rot, float), np.array(failed)


def _scores(raw, internal, rot, weights, w_rot, terms=None):
    """The best score of every compound under these weights (a failed compound scores +inf)."""
    from rnafold3d.screening import score
    terms = tuple(terms or score.TERMS)
    w = np.array([weights[t] for t in terms], float)
    out = np.full(len(raw), np.inf)
    for k, (m, intern) in enumerate(zip(raw, internal)):
        if not len(m):
            continue
        out[k] = float(np.min(m @ w / (1.0 + w_rot * rot[k]) + intern))
    return out


def auc(scores, active):
    """The probability that an active scores better than a decoy (Mann-Whitney U, ties count a half)."""
    a, d = scores[active], scores[~active]
    if not len(a) or not len(d):
        return float('nan')
    finite = np.isfinite(a)
    worst = np.max(scores[np.isfinite(scores)], initial=0.0) + 1.0
    a = np.where(finite, a, worst)
    d = np.where(np.isfinite(d), d, worst)
    better = (a[:, None] < d[None, :]).sum() + 0.5 * (a[:, None] == d[None, :]).sum()
    return float(better / (len(a) * len(d)))


def enrichment(scores, active, fraction):
    """The enrichment factor at the top `fraction` of the list: how many times better than chance."""
    n = len(scores)
    k = max(1, int(round(fraction * n)))
    order = np.argsort(np.where(np.isfinite(scores), scores, np.inf), kind='stable')
    hits = int(active[order[:k]].sum())
    expected = fraction * active.sum()
    return float(hits / expected) if expected else float('nan'), hits, k


def _report_target(target, weights=None, quiet=False):
    from rnafold3d.screening import score
    labels, active, raw, internal, rot, failed = _matrix(target)
    w = dict(score.WEIGHTS, **(weights or {}))
    s = _scores(raw, internal, rot, w, w['rotatable'])
    value = auc(s, active)
    ef1, h1, k1 = enrichment(s, active, 0.01)
    ef5, h5, k5 = enrichment(s, active, 0.05)
    if not quiet:
        order = np.argsort(np.where(np.isfinite(s), s, np.inf), kind='stable')
        top = ', '.join(f'{labels[i]}{"*" if active[i] else ""} {s[i]:.1f}' for i in order[:6])
        print(f'  AUC {value:.3f}   EF1% {ef1:5.1f} ({h1}/{k1})   EF5% {ef5:4.1f} ({h5}/{k5})   '
              f'{int(active.sum())} actives in {len(s)} compounds, {int(failed.sum())} undockable')
        print(f'  best six: {top}   (* = active)')
    return {'key': target['key'], 'name': target['name'], 'auc': round(value, 4),
            'ef1': round(ef1, 2), 'ef5': round(ef5, 2), 'n': len(s), 'n_actives': int(active.sum()),
            'undockable': int(failed.sum())}


def _report(payload, weights=None):
    rows = []
    for target in payload['targets']:
        print(f'{target["name"]} ({target["pdb"]}, {target["what"]})')
        rows.append(_report_target(target, weights))
    good = [r for r in rows if np.isfinite(r['auc'])]
    if good:
        print(f'\n{len(good)} targets, {payload["exhaustiveness"]}, '
              f'{payload["ratio"]} property-matched decoys per active')
        print(f'  mean AUC  {np.mean([r["auc"] for r in good]):.3f}  '
              f'(worst {min(r["auc"] for r in good):.3f} {min(good, key=lambda r: r["auc"])["key"]}, '
              f'best {max(r["auc"] for r in good):.3f} {max(good, key=lambda r: r["auc"])["key"]})')
        print(f'  mean EF1% {np.mean([r["ef1"] for r in good]):.1f}   '
              f'mean EF5% {np.mean([r["ef5"] for r in good]):.1f}')
        print('  a random ranking is AUC 0.500 and EF 1.0; the plan asks for at least 0.65')
    return rows


def report(args):
    payload = json.loads(Path(args.file).read_text(encoding='utf-8'))
    weights = json.loads(Path(args.weights).read_text(encoding='utf-8'))['weights'] if args.weights else None
    _report(payload, weights)


# --------------------------------------------------------------------------------------- tuning

TUNABLE = ('stacking', 'electrostatic', 'desolvation', 'metal')
FACTOR = 4.0                   # a tuned weight may move by this factor, and may not change sign


def _soft_auc(scores, active, temperature=0.7):
    """A differentiable stand-in for the AUC: the mean over active-decoy pairs of how clearly the active
    wins. The AUC itself is a step function of the weights and cannot be optimised directly."""
    a, d = scores[active], scores[~active]
    if not len(a) or not len(d):
        return 0.0
    worst = np.max(scores[np.isfinite(scores)], initial=0.0) + 1.0
    a = np.where(np.isfinite(a), a, worst)
    d = np.where(np.isfinite(d), d, worst)
    return float(np.mean(1.0 / (1.0 + np.exp((a[:, None] - d[None, :]) / temperature))))


def _prepared(payload, keys=None):
    out = []
    for target in payload['targets']:
        if keys and target['key'] not in keys:
            continue
        labels, active, raw, internal, rot, _failed = _matrix(target)
        out.append((target['key'], active, raw, internal, rot, labels))
    return out


def _fit(data, tunable, base, temperature):
    """Maximise the mean soft AUC over the targets in `data` by moving the tunable weights only."""
    from scipy.optimize import minimize
    x0 = np.array([base[t] for t in tunable], float)
    bounds = [tuple(sorted((v / FACTOR, v * FACTOR))) for v in x0]

    def cost(x):
        w = dict(base, **{t: float(v) for t, v in zip(tunable, x)})
        vals = [_soft_auc(_scores(raw, intern, rot, w, w['rotatable']), active, temperature)
                for _k, active, raw, intern, rot, _l in data]
        return -float(np.mean(vals))

    res = minimize(cost, x0, method='L-BFGS-B', bounds=bounds,
                   options={'maxiter': 200, 'eps': 1e-3})
    return dict(base, **{t: float(v) for t, v in zip(tunable, res.x)}), -float(res.fun)


def tune(args):
    from rnafold3d.screening import score
    payload = json.loads(Path(args.file).read_text(encoding='utf-8'))
    base = dict(score.WEIGHTS)
    tunable = tuple(args.terms.split(',')) if args.terms else TUNABLE
    data = _prepared(payload)
    keys = [k for k, *_ in data]
    print(f'tuning {", ".join(tunable)} on {len(data)} targets '
          f'(each weight within a factor of {FACTOR:.0f} of its default, sign fixed)')
    start = np.mean([auc(_scores(r, i, t, base, base['rotatable']), a) for _k, a, r, i, t, _l in data])
    fitted, soft = _fit(data, tunable, base, args.temperature)
    inside = np.mean([auc(_scores(r, i, t, fitted, fitted['rotatable']), a) for _k, a, r, i, t, _l in data])
    print(f'  default weights: mean AUC {start:.3f}')
    print(f'  fitted on all six: mean AUC {inside:.3f} (soft {soft:.3f})')
    for t in tunable:
        print(f'    {t:14s} {base[t]:+8.4f} → {fitted[t]:+8.4f}  '
              f'({fitted[t] / base[t]:.2f}× its default)')
    # Leave one target out: fit on the others, measure on the one held back. This is the number that says
    # whether the weights transfer, and with six targets it is the only honest one.
    held = {}
    for key in keys:
        rest = [d for d in data if d[0] != key]
        this = [d for d in data if d[0] == key]
        w, _ = _fit(rest, tunable, base, args.temperature)
        before = auc(_scores(this[0][2], this[0][3], this[0][4], base, base['rotatable']), this[0][1])
        after = auc(_scores(this[0][2], this[0][3], this[0][4], w, w['rotatable']), this[0][1])
        held[key] = (before, after, w)
        print(f'  held out {key:8s} AUC {before:.3f} → {after:.3f}')
    loo = np.mean([v[1] for v in held.values()])
    print(f'  cross-validated mean AUC {np.mean([v[0] for v in held.values()]):.3f} → {loo:.3f}')
    if loo <= np.mean([v[0] for v in held.values()]) + 0.005:
        print('  the fit does not transfer: keep the default weights')
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(
            {'weights': fitted, 'tunable': list(tunable), 'factor': FACTOR,
             'mean_auc_default': round(float(start), 4), 'mean_auc_fitted': round(float(inside), 4),
             'mean_auc_cross_validated': round(float(loo), 4),
             'held_out': {k: [round(v[0], 4), round(v[1], 4)] for k, v in held.items()},
             'source': str(args.file)}, indent=1), encoding='utf-8')
        print(f'written to {args.out}')


# --------------------------------------------------------------------------------------- cli

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fetch', help='resolve the actives, pick property-matched decoys, download 3D')
    f.add_argument('--ratio', type=int, default=25, help='decoys per active (default 25)')
    f.add_argument('--pool', type=int, default=POOL, help='background compounds to sample from PubChem')
    f.add_argument('--pool-seed', type=int, default=0)
    f.add_argument('--only', help='comma-separated target keys')
    r = sub.add_parser('run', help='dock every compound of every target')
    r.add_argument('--out')
    r.add_argument('--only')
    r.add_argument('--exhaustiveness', default='quick', choices=('quick', 'standard', 'thorough'))
    r.add_argument('--seed', type=int, default=1)
    r.add_argument('--n-poses', type=int, default=5)
    r.add_argument('--box', type=float, default=6.0, help='Å around the experimental ligand')
    r.add_argument('--jobs', type=int, default=max(1, 4), help='worker processes')
    r.add_argument('--limit', type=int, help='only this many compounds per set (for a smoke run)')
    p = sub.add_parser('report', help='AUC and enrichment of a finished run')
    p.add_argument('file')
    p.add_argument('--weights', help='a JSON file from "tune", to report that ranking instead')
    t = sub.add_parser('tune', help='fit the RNA weights to the run, with leave-one-target-out')
    t.add_argument('file')
    t.add_argument('--terms', help=f'comma-separated, default {",".join(TUNABLE)}')
    t.add_argument('--temperature', type=float, default=0.7)
    t.add_argument('--out')
    a = ap.parse_args(argv)
    {'fetch': fetch, 'run': run, 'report': report, 'tune': tune}[a.cmd](a)


if __name__ == '__main__':
    main()
