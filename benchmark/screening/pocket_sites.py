"""Does the pocket finder find the sites where compounds actually bind?

Pocket finding (`rnafold3d/screening/pockets.py`) is a geometric analysis, and the only honest way to check
it is against sites that are known: take the experimental RNA-ligand complexes of the redocking benchmark,
**remove the compound**, run the pocket finder on the RNA alone, and ask where in its ranking the site the
compound occupied appears.

    python benchmark/screening/pocket_sites.py                  # uses the complexes fetched by redock.py
    python benchmark/screening/pocket_sites.py --burial-min 5 --out runs/pocket_sites.json

A site counts as found when a pocket's lining covers at least half of the nucleotides the compound touches
(an atom within 4.5 Å). The ranking is the pocket score, which without base-pairing information is geometry
alone — the same number the Details tab shows, except that the tab also knows which nucleotides are in a
helix and uses it.

The limits are the same as the redocking benchmark's: these are the compound's *own* structures, so the
pocket is already open and shaped. On a predicted model, or on a conformation that does not bind, nothing
here applies.
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from redock import CACHE, SET                                        # noqa: E402

CONTACT = 4.5


def one(pdb, ligname, burial_min=6, top=10):
    from rnafold3d.screening import pockets as P
    from rnafold3d.screening.receptor import read_complex
    from scipy.spatial import cKDTree
    path = CACHE / f'{pdb}.pdb'
    if not path.exists():
        return None
    rec, ligands = read_complex(path)
    want = [m for m in ligands if ligname in m.properties.get('residues', [m.properties.get('residue')])]
    if not want:
        return {'pdb': pdb, 'error': 'ligand not found'}
    lig = want[0]
    near = cKDTree(lig.xyz).query_ball_point(rec.xyz, CONTACT, return_length=True)
    site = sorted({int(rec.residue[a]) for a in np.nonzero(near)[0]})
    if len(site) < 3:
        return {'pdb': pdb, 'error': 'the compound touches almost nothing'}
    t0 = time.time()
    found = P.find(rec.xyz, elements=rec.element, residue_of=rec.residue, burial_min=burial_min,
                   keep_points=False)
    cover = [len(set(p['residues']) & set(site)) / len(site) for p in found]
    rank = next((k + 1 for k, c in enumerate(cover) if c >= 0.5), None)
    centre = [float(np.linalg.norm(np.asarray(p['center']) - lig.xyz.mean(0))) for p in found]
    return {'pdb': pdb, 'ligand': '+'.join(lig.properties.get('residues', [ligname])), 'atoms': len(lig),
            'site': site, 'site_size': len(site), 'pockets': len(found),
            'rank': rank, 'best_cover': round(max(cover), 2) if cover else 0.0,
            'cover_of_top': round(cover[0], 2) if cover else 0.0,
            'nearest_centre': round(min(centre), 1) if centre else None,
            'volumes': [p['volume'] for p in found[:top]], 'seconds': round(time.time() - t0, 2)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--burial-min', type=int, default=6, choices=(5, 6, 7))
    ap.add_argument('--only')
    ap.add_argument('--out')
    a = ap.parse_args(argv)
    rows = []
    for pdb, lig, what in SET:
        if a.only and pdb not in a.only.split(','):
            continue
        r = one(pdb, lig, a.burial_min)
        if r is None:
            print(f'{pdb}: not downloaded (run "redock.py fetch" first)')
            continue
        if 'error' in r:
            print(f'{pdb}: {r["error"]}')
        else:
            print(f'{pdb} {r["ligand"]:14s} site {r["site_size"]:2d} nt  {r["pockets"]:2d} pockets  '
                  f'rank {str(r["rank"] or "-"):>3s}  best cover {r["best_cover"]:.2f}  '
                  f'top-pocket cover {r["cover_of_top"]:.2f}  {r["seconds"]:4.1f} s   {what}')
        rows.append(r)
    ok = [r for r in rows if 'error' not in r and r.get('pockets') is not None]
    if ok:
        rank = [r['rank'] for r in ok]
        print(f'\n{len(ok)} complexes, sites of at least half the contacting nucleotides, '
              f'burial {a.burial_min} of 7')
        for k in (1, 3, 5):
            print(f'  site among the top {k} pockets: {sum(1 for x in rank if x and x <= k)} of {len(ok)}')
        print(f'  site found at all:            {sum(1 for x in rank if x)} of {len(ok)}')
        print(f'  median pockets per structure: {np.median([r["pockets"] for r in ok]):.0f}')
        print(f'  median coverage of the best pocket: {np.median([r["best_cover"] for r in ok]):.2f}')
        missed = [r['pdb'] for r in ok if not r['rank']]
        if missed:
            print('  not found: ' + ', '.join(missed))
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(json.dumps({'burial_min': a.burial_min, 'results': rows}, indent=1),
                               encoding='utf-8')
        print(f'written to {a.out}')


if __name__ == '__main__':
    main()
