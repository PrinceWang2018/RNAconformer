"""Benchmark 1 of the compound screen: can the scoring function put a ligand back where it belongs?

A redocking benchmark is the first thing a docking score has to pass, and the only one that does not depend
on any claim about affinity: take an experimental RNA-ligand complex, remove the ligand, dock it back into
its own structure, and measure how far the best-scored pose is from where the crystal or the NMR model put
it. Below 2.5 Å is the usual line for "the right pose".

    python benchmark/screening/redock.py fetch                       # download the complexes (once)
    python benchmark/screening/redock.py run --out runs/redock.json
    python benchmark/screening/redock.py report runs/redock.json

What this is and is not:

* The set here is **hand-picked**, not the curated HARIBOSS set the plan asks for (`docs/NEXT_STEPS.md`
  §1.4): about thirty well-known RNA-ligand complexes — riboswitches, aptamers, a viral element, a
  ribosomal fragment — chosen to span sizes, charges and flexibility. It is enough to see whether the
  scoring function is in the right place and to catch a regression; it is not enough to publish a number.
* Redocking into the ligand's own structure is the easy case: the pocket is already the right shape. It
  measures the scoring function and the search, not what a screen will do against a predicted model, and
  it says nothing about enrichment (benchmark 2) or selectivity (benchmark 3).
* The RMSD allows atoms of the same element to swap (`dock.rmsd(..., elements)`), so a compound flipped by
  its own symmetry counts as the same pose. That is a lower bound of the proper graph-symmetry value.

The structures are downloaded from RCSB into ~/.rcv/benchmarks/rna_ligand (not into the repository) and
are the property of their depositors.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

CACHE = Path.home() / '.rcv' / 'benchmarks' / 'rna_ligand'
URL = 'https://files.rcsb.org/download/{}.pdb'

# pdb id, ligand residue name, what it is. Hand-picked; see the module text.
SET = [
    ('1Y26', 'ADE', 'adenine riboswitch (add), adenine'),
    ('1U8D', 'HPA', 'guanine riboswitch (xpt), hypoxanthine'),
    ('4FE5', 'HPA', 'guanine riboswitch, hypoxanthine (second entry)'),
    ('3GX5', 'SAM', 'SAM-II riboswitch, S-adenosylmethionine'),
    ('2QWY', 'SAM', 'SAM-III riboswitch, S-adenosylmethionine'),
    ('4NYB', '2QC', 'SAM-V riboswitch, SAM analogue'),
    ('2GDI', 'TPP', 'TPP riboswitch (thiM), thiamine pyrophosphate'),
    ('3D2G', 'TPP', 'TPP riboswitch, thiamine pyrophosphate (second entry)'),
    ('3IRW', 'C2E', 'c-di-GMP riboswitch, cyclic di-GMP (46 atoms, the largest here)'),
    ('3SKI', 'GNG', 'riboswitch with a guanosine analogue'),
    ('3Q50', 'PRF', 'preQ1 riboswitch class II, preQ1'),
    ('5BTP', 'AMZ', 'riboswitch with a purine analogue'),
    ('4TS2', '38E', 'riboswitch with a synthetic ligand'),
    ('2B57', '6AP', 'RNA with 2,6-diaminopurine'),
    ('1F1T', 'ROS', 'malachite green aptamer, malachite green (rigid, cationic)'),
    ('1FMN', 'FMN', 'FMN aptamer, flavin mononucleotide (NMR)'),
    ('1EHT', 'TEP', 'theophylline aptamer, theophylline (NMR)'),
    ('1O15', 'TEP', 'theophylline aptamer, theophylline (second entry)'),
    ('1AJU', 'ARG', 'HIV-1 TAR, argininamide (NMR)'),
    ('1ARJ', 'ARG', 'HIV-1 TAR, argininamide (NMR, second entry)'),
    ('1LVJ', 'PMZ', 'A-site RNA with an aminoglycoside-like ligand (NMR)'),
    ('1J7T', 'PAR', 'A-site RNA, paromomycin (+5, 9 rotatable bonds)'),
    ('2ET4', 'NMY', 'A-site RNA, neomycin B'),
    ('1EI2', 'NMY', 'neomycin aptamer, neomycin B (NMR)'),
    ('2KGP', 'MIX', 'neomycin riboswitch, ribostamycin (NMR)'),
    ('1NEM', 'NEB', 'neomycin aptamer (NMR); deposited as linked sugars'),
    ('1TOB', 'TOA', 'tobramycin aptamer (NMR); deposited as linked sugars'),
    ('2TOB', 'TOA', 'tobramycin aptamer (NMR, second entry); linked sugars'),
    ('1BYJ', 'GE1', 'A-site RNA, paromomycin (NMR); deposited as linked sugars'),
    ('1QD3', 'BDG', 'RNA aptamer with an aminoglycoside; linked sugars'),
]


def fetch(ids=None, force=False):
    """Download the complexes into the cache. Needs network access, once."""
    import urllib.request
    CACHE.mkdir(parents=True, exist_ok=True)
    got, missing = [], []
    for pdb, lig, what in SET:
        if ids and pdb not in ids:
            continue
        out = CACHE / f'{pdb}.pdb'
        if out.exists() and not force:
            got.append(pdb)
            continue
        try:
            with urllib.request.urlopen(URL.format(pdb), timeout=60) as r:
                data = r.read()
            if len(data) < 1000:
                raise ValueError('short file')
            out.write_bytes(data)
            got.append(pdb)
            print(f'  {pdb}  {len(data) / 1000:6.0f} kB  {what}')
        except Exception as err:                                  # noqa: BLE001 - report and carry on
            missing.append((pdb, str(err)))
            print(f'  {pdb}  FAILED: {err}')
    print(f'{len(got)} of {len(SET)} complexes in {CACHE}')
    return got, missing


def one(pdb, ligname, exhaustiveness='standard', seed=1, box=6.0, blind=False, n_poses=9):
    """Redock one complex. Returns a record, or None when the entry cannot be used."""
    from rnafold3d.screening import dock, grid, score
    from rnafold3d.screening.receptor import read_complex
    path = CACHE / f'{pdb}.pdb'
    if not path.exists():
        return None
    rec, ligands = read_complex(path)
    want = [m for m in ligands if ligname in m.properties.get('residues', [m.properties.get('residue')])]
    if not want:
        return {'pdb': pdb, 'ligand': ligname, 'error': 'ligand not found in the entry',
                'found': sorted('+'.join(m.properties.get('residues', [])) for m in ligands)}
    lig = want[0]
    t0 = time.time()
    if blind:
        from rnafold3d.screening import pockets as P
        found = P.find(rec.xyz, elements=rec.element, residue_of=rec.residue)
        if not found:
            return {'pdb': pdb, 'ligand': ligname, 'error': 'no pocket found'}
        g = grid.Grid.for_pocket(rec, pocket=found[0])
        site = float(np.linalg.norm(np.asarray(found[0]['center']) - lig.xyz.mean(0)))
    else:
        span = float((lig.xyz.max(0) - lig.xyz.min(0)).max())
        g = grid.Grid.for_pocket(rec, centre=lig.xyz.mean(0), size=span + 2 * box)
        site = 0.0
    build = time.time() - t0
    t0 = time.time()
    poses = dock.dock(lig, g, receptor=rec, exhaustiveness=exhaustiveness, seed=seed, n_poses=n_poses)
    rmsds = [dock.rmsd(p.xyz, lig.xyz, lig.elements) for p in poses]
    return {'pdb': pdb, 'ligand': '+'.join(lig.properties.get('residues', [ligname])),
            'atoms': len(lig), 'rotatable': lig.n_rotatable,
            'charge': round(float(lig.charge.sum()), 2), 'receptor_atoms': len(rec),
            'crystal_score': round(score.score(lig.xyz, lig, rec), 3),
            'scores': [round(p.score, 3) for p in poses], 'rmsd': [round(r, 2) for r in rmsds],
            'best_rmsd': round(min(rmsds), 2) if rmsds else None,
            'top_rmsd': round(rmsds[0], 2) if rmsds else None,
            'blind': bool(blind), 'site_distance': round(site, 1),
            'seconds': round(time.time() - t0, 1), 'grid_seconds': round(build, 1),
            'grid': {'shape': list(g.shape), 'spacing': g.spacing, 'mb': round(g.nbytes / 1e6, 1)}}


def run(args):
    out = []
    todo = [(p, l, w) for p, l, w in SET if not args.only or p in args.only.split(',')]
    for pdb, lig, what in todo:
        r = one(pdb, lig, args.exhaustiveness, args.seed, args.box, args.blind)
        if r is None:
            print(f'{pdb}: not downloaded (run "fetch" first)')
            continue
        if 'error' in r:
            print(f'{pdb} {lig}: {r["error"]} {r.get("found", "")}')
        else:
            print(f'{pdb} {r["ligand"]:14s} {r["atoms"]:3d} atoms {r["rotatable"]:2d}rot q{r["charge"]:+.0f}  '
                  f'top {r["top_rmsd"]:5.2f} Å  best {r["best_rmsd"]:5.2f} Å  score {r["scores"][0]:7.2f} '
                  f'(crystal {r["crystal_score"]:7.2f})  {r["seconds"]:5.1f} s   {what}')
        r['description'] = what
        out.append(r)
    payload = {'set': 'hand-picked RNA-ligand complexes', 'n': len(out), 'blind': bool(args.blind),
               'exhaustiveness': args.exhaustiveness, 'seed': args.seed, 'results': out}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=1), encoding='utf-8')
        print(f'written to {args.out}')
    summarise(payload)


def summarise(payload):
    ok = [r for r in payload['results'] if r.get('top_rmsd') is not None]
    if not ok:
        print('no results')
        return
    top = np.array([r['top_rmsd'] for r in ok])
    best = np.array([r['best_rmsd'] for r in ok])
    print(f'\n{len(ok)} complexes, {payload["exhaustiveness"]}'
          f'{", blind (pocket finder)" if payload.get("blind") else ", into the ligand\'s own site"}')
    print(f'  best-scored pose:  median {np.median(top):.2f} Å, '
          f'under 2.5 Å in {int((top < 2.5).sum())} of {len(ok)} ({100 * (top < 2.5).mean():.0f}%)')
    print(f'  any kept pose:     median {np.median(best):.2f} Å, '
          f'under 2.5 Å in {int((best < 2.5).sum())} of {len(ok)} ({100 * (best < 2.5).mean():.0f}%)')
    found = np.array([r['scores'][0] <= r['crystal_score'] for r in ok])
    print(f'  the search found a pose scoring at least as well as the experimental one in '
          f'{int(found.sum())} of {len(ok)}')
    print(f'  total {sum(r["seconds"] for r in ok) / 60:.1f} min of docking, '
          f'{sum(r["grid_seconds"] for r in ok) / 60:.1f} min of grid building')
    hard = sorted(ok, key=lambda r: -r['top_rmsd'])[:5]
    print('  worst: ' + ', '.join(f'{r["pdb"]} {r["top_rmsd"]:.1f} Å' for r in hard))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    f = sub.add_parser('fetch', help='download the complexes from RCSB')
    f.add_argument('--force', action='store_true')
    r = sub.add_parser('run', help='redock them')
    r.add_argument('--out')
    r.add_argument('--only', help='comma-separated PDB ids')
    r.add_argument('--exhaustiveness', default='standard', choices=('quick', 'standard', 'thorough'))
    r.add_argument('--seed', type=int, default=1)
    r.add_argument('--box', type=float, default=6.0, help='Å added around the ligand on each side')
    r.add_argument('--blind', action='store_true', help='dock into the best pocket instead of the known site')
    p = sub.add_parser('report', help='summarise a run')
    p.add_argument('file')
    a = ap.parse_args(argv)
    if a.cmd == 'fetch':
        fetch(force=a.force)
    elif a.cmd == 'run':
        run(a)
    else:
        summarise(json.loads(Path(a.file).read_text(encoding='utf-8')))


if __name__ == '__main__':
    main()
