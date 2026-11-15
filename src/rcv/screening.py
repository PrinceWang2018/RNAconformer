"""Compound screening in the viewer: a library against two conformations of one RNA.

The engine side is `rnafold3d.screening` (pockets, the scoring function, its grid, the search and the
two-state screen); this file is what the window talks to. It holds the parsed libraries of this session,
turns what the user picked in the window into a pair of boxes on the two conformations, and runs the screen
as an ordinary RCV job so that progress, cancelling and the job list work as they do for a prediction.

What the window may say, and what it may not: the numbers here are a **ranking**, with the measured limits
of `docs/SCREENING_STATUS.md` shown next to them. Nothing is an affinity, nothing is in kcal/mol, and no
compound is said to block anything. The screen is behind *experimental* in the interface until the three
benchmarks of the plan (`docs/NEXT_STEPS.md` §1.4) are all run and reported.
"""
from __future__ import annotations

import secrets
import threading
import time
from pathlib import Path

import numpy as np

from . import rcvfile
from .jobs import Job, plural

MAX_LIBRARY = 200_000          # compounds parsed from one file
PREVIEW_ROWS = 300             # rows the job's status carries; the full table is a separate request
LIBRARY_LIMIT = 8              # parsed libraries kept per session

_LIBRARIES = {}
_LOCK = threading.Lock()


# ------------------------------------------------------------------------------------ libraries

def parse_library(text, name='', ph=7.0, limit=MAX_LIBRARY, fmt=None):
    """Parse a library file (SDF, MOL, MOL2) into prepared compounds, kept under a handle.

    Returns what the window shows before a run: how many compounds were read, how many were rejected and
    why, the spread of their size and charge, and a handle to start a screen with. Rejected records are
    counted by reason rather than dropped silently — a library that is half unreadable is a mistake the
    user has to see before waiting for a run.
    """
    from rnafold3d.screening import chem
    text = text.decode('utf-8', 'replace') if isinstance(text, (bytes, bytearray)) else str(text)
    kind = (fmt or ('mol2' if '@<TRIPOS>' in text[:4000] else 'sdf')).lower()
    reader = chem.read_mol2 if kind == 'mol2' else chem.read_sdf
    mols, rejected = [], {}
    for mol in reader(text, limit=limit):
        why = mol.properties.get('error') if mol.properties else None
        if why or not len(mol):
            why = why or 'no atoms'
            rejected[why] = rejected.get(why, 0) + 1
            continue
        if ph is not None and ph != 7.0:
            mol = chem.prepare(mol.elements, mol.xyz, mol.bonds, hydrogens=mol.hydrogens,
                               name=mol.name, properties=mol.properties, ph=ph)
        mols.append(mol)
    handle = secrets.token_hex(5)
    rec = {'handle': handle, 'name': name or 'library', 'molecules': mols, 'format': kind,
           'rejected': rejected, 'created': time.time(), 'ph': ph}
    with _LOCK:
        _LIBRARIES[handle] = rec
        for old in sorted(_LIBRARIES.values(), key=lambda r: r['created'])[:-LIBRARY_LIMIT]:
            _LIBRARIES.pop(old['handle'], None)
    return library_summary(handle)


def library_summary(handle):
    from rnafold3d.screening.screen import _weight
    rec = _LIBRARIES.get(handle)
    if rec is None:
        raise KeyError('this library is no longer loaded: choose the file again')
    mols = rec['molecules']
    atoms = np.array([len(m) for m in mols]) if mols else np.zeros(0)
    charge = np.array([float(m.charge.sum()) for m in mols]) if mols else np.zeros(0)
    rot = np.array([m.n_rotatable for m in mols]) if mols else np.zeros(0)
    return {'handle': handle, 'name': rec['name'], 'format': rec['format'], 'ph': rec['ph'],
            'count': len(mols), 'rejected': rec['rejected'],
            'rejected_total': int(sum(rec['rejected'].values())),
            'atoms': [int(atoms.min()), int(np.median(atoms)), int(atoms.max())] if len(atoms) else [],
            'weight': [round(float(np.min([_weight(m) for m in mols])), 1),
                       round(float(np.max([_weight(m) for m in mols])), 1)] if mols else [],
            'charge': [round(float(charge.min()), 1), round(float(charge.max()), 1)] if len(charge) else [],
            'rotatable': int(rot.max()) if len(rot) else 0,
            'examples': [{'name': m.name or f'compound {i + 1}', 'formula': m.formula, 'atoms': len(m),
                          'charge': round(float(m.charge.sum()), 1), 'rotatable': m.n_rotatable}
                         for i, m in enumerate(mols[:6])]}


def library_molecules(handle):
    rec = _LIBRARIES.get(handle)
    if rec is None:
        raise KeyError('this library is no longer loaded: choose the file again')
    return rec['molecules']


# ------------------------------------------------------------------------------------ the site

def site_of(doc, coords, k, spec, burial_min=None):
    """Turn what the user picked into a box on conformation k, and say what it came from.

    `spec` is one of
        {'mode': 'pocket', 'pocket': 1}            a pocket of `pockets.find`, by its rank (1-based)
        {'mode': 'residues', 'residues': [...]}    the current selection in the main window
        {'mode': 'box', 'centre': [x,y,z], 'size': 24}
    """
    from rnafold3d.screening import pockets as P
    from rnafold3d.screening.screen import Site
    from .insight import _flexible
    xyz = np.asarray(coords, float)[k]
    mode = (spec or {}).get('mode', 'pocket')
    if mode == 'box':
        site = Site(centre=tuple(float(c) for c in spec['centre']), size=float(spec.get('size', 24.0)),
                    label='a box the user gave')
        return site, {'mode': mode}
    if mode == 'residues':
        res = [int(r) for r in spec.get('residues', [])]
        if not res:
            raise ValueError('no nucleotides are selected: select the site in the main window, or pick a pocket')
        atoms = np.isin(np.asarray(doc['atoms']['residue'], int), res)
        site = Site.from_atoms(xyz[atoms], label='the selected nucleotides', residues=res)
        return site, {'mode': mode, 'residues': res, 'n_atoms': int(atoms.sum())}
    found = P.find(xyz, elements=doc['atoms']['elements'], residue_of=doc['atoms']['residue'],
                   flexible=_flexible(doc, k), **({} if burial_min is None else {'burial_min': int(burial_min)}))
    if not found:
        raise ValueError('no pocket was found in this conformation: select the site in the main window instead')
    want = int((spec or {}).get('pocket', 1))
    pocket = next((p for p in found if p['id'] == want), found[0])
    site = Site.from_pocket(pocket)
    return site, {'mode': 'pocket', 'pocket': pocket['id'], 'volume': pocket['volume'],
                  'score': pocket['score'], 'residues': pocket['residues'], 'wide': pocket['wide'],
                  'n_pockets': len(found)}


def matching_site(doc, coords, k, site, residues):
    """The same site on the other conformation: the box around the nucleotides that lined the first one.

    Not the pockets of the second conformation — the question a two-state screen asks is what happens to
    *this* site when the RNA changes shape, and the pocket may have closed there entirely (which is the
    interesting answer, and it has no pocket to match).
    """
    from rnafold3d.screening.screen import Site
    xyz = np.asarray(coords, float)[k]
    res = [int(r) for r in (residues or site.residues)]
    if not res:
        return Site(centre=site.centre, size=site.size, label='the same box'), {'mode': 'box'}
    atoms = np.isin(np.asarray(doc['atoms']['residue'], int), res)
    if not atoms.any():
        return Site(centre=site.centre, size=site.size, label='the same box'), {'mode': 'box'}
    other = Site.from_atoms(xyz[atoms], label='the same nucleotides', residues=res)
    return other, {'mode': 'residues', 'residues': res}


# ------------------------------------------------------------------------------------ the job

class ScreenJob(Job):
    """Dock every compound of a library into one or two conformations and rank them by selectivity."""

    kind = 'screen'
    STAGES = (('library', 'Library'), ('site', 'Pocket'), ('grid', 'Grids'),
              ('dock', 'Docking'), ('rank', 'Ranking'))

    def __init__(self, state, params):
        entry = state.entries.get(params.get('document') or '')
        if entry is None:
            raise ValueError('unknown ensemble')
        self.entry = entry
        self.params = dict(params or {})
        library_summary(self.params.get('library') or '')     # fail now, not after the grids are built
        doc = entry.doc
        name = doc.get('name') or 'ensemble'
        self.a = int(params.get('a', 0))
        self.b = None if params.get('b') is None else int(params['b'])
        n = len(doc['conformations'])
        if not 0 <= self.a < n or (self.b is not None and not 0 <= self.b < n):
            raise ValueError('conformation index out of range')
        if self.b == self.a:
            self.b = None                      # A against itself has no selectivity to report
        titles = [c.get('label') or f'conformation {i + 1}' for i, c in enumerate(doc['conformations'])]
        both = f'{titles[self.a]} vs {titles[self.b]}' if self.b is not None else titles[self.a]
        super().__init__(state, f'Compound screen: {name}, {both}')
        self.rows = []
        self.summary = {}
        self.sites = {}
        self.library = None

    # ---------------------------------------------------------------- reporting
    def info(self):
        out = super().info()
        out['result'] = dict(self.result or {}, rows=self.rows[:PREVIEW_ROWS], n_rows=len(self.rows))
        return out

    def table(self, limit=None, offset=0, poses=False):
        """The full table, for the window's results view and the exports."""
        from rnafold3d.screening.screen import sort_key
        rows = sorted(self.rows, key=sort_key)
        rows = rows[int(offset):int(offset) + int(limit)] if limit else rows[int(offset):]
        if not poses:
            rows = [{k: v for k, v in r.items() if k != 'poses'} for r in rows]
        return {'rows': rows, 'total': len(self.rows), 'summary': self.result, 'sites': self.sites,
                'parameters': self.parameters()}

    def parameters(self):
        p = self.params
        return {'a': self.a, 'b': self.b, 'exhaustiveness': p.get('exhaustiveness', 'quick'),
                'n_poses': int(p.get('n_poses', 3)), 'keep_poses': int(p.get('keep_poses', 200)),
                'ph': p.get('ph', 7.0), 'jobs': self.n_jobs, 'library': (self.library or {}).get('name'),
                'burial_min': p.get('burial_min'), 'weights': 'measured defaults'}

    @property
    def n_jobs(self):
        from rnafold3d.screening.screen import workers
        return workers(self.params.get('jobs'), per_worker_mb=120)

    # ---------------------------------------------------------------- the run
    def run(self):
        from rnafold3d.screening.receptor import Receptor
        from rnafold3d.screening.screen import Screen, run as run_screen, summarise
        doc = self.entry.doc
        coords = rcvfile.doc_coords(doc)
        p = self.params

        self.stage('library', 'active')
        self.library = library_summary(p['library'])
        mols = library_molecules(p['library'])
        if int(p.get('limit') or 0):
            mols = mols[:int(p['limit'])]
        if not mols:
            raise ValueError('the library has no readable compounds')
        self.stage('library', 'done', f'{plural(len(mols), "compound")}')
        self.note(f'{plural(len(mols), "compound")} from {self.library["name"]}'
                  + (f'; {plural(self.library["rejected_total"], "record")} could not be read'
                     if self.library['rejected_total'] else '.'))

        self.stage('site', 'active')
        site_a, info_a = site_of(doc, coords, self.a, p.get('site'), p.get('burial_min'))
        rec_a = Receptor.from_document(doc, coords, self.a)
        self.sites = {'a': dict(site_a.as_dict(), **info_a)}
        rec_b = site_b = None
        if self.b is not None:
            site_b, info_b = matching_site(doc, coords, self.b, site_a, info_a.get('residues'))
            rec_b = Receptor.from_document(doc, coords, self.b)
            self.sites['b'] = dict(site_b.as_dict(), **info_b)
        self.stage('site', 'done', site_a.label or 'a box')
        if site_a.cropped:
            self.note('The pocket is wider than the largest box a screen builds a grid for, so only its '
                      'centre is screened. A site this wide is usually several clefts merged into one: '
                      'select the nucleotides of the one you mean in the main window instead.', 'warning')

        self.stage('grid', 'active')
        screen = Screen(rec_a, site_a, rec_b, site_b, exhaustiveness=p.get('exhaustiveness', 'quick'),
                        seed=int(p.get('seed', 1)), n_poses=int(p.get('n_poses', 3)))
        screen.prepare()
        self.stage('grid', 'done', f'{screen.nbytes / 1e6:.0f} MB')

        self.stage('dock', 'active')
        jobs, total = self.n_jobs, len(mols)
        self.note(f'{plural(total, "compound")} on {plural(jobs, "core")}, '
                  f'{p.get("exhaustiveness", "quick")} search.')
        t0 = time.time()

        def progress(done, _total, row):
            self.rows.append(row)
            self.progress = 0.1 + 0.85 * done / max(1, total)
            rate = done / max(1e-6, time.time() - t0)
            self.eta_s = (total - done) / rate if rate > 0 else None
            self.stage('dock', 'active', f'{done} of {total}, {rate:.1f} per second')
            if done % 50 == 0 or done == total:
                self.result = dict(summarise(self.rows, top=0), done=done, total=total,
                                   per_second=round(rate, 2))

        self.rows = []
        run_screen(screen, mols, jobs=jobs, keep_poses=int(p.get('keep_poses', 200)),
                   progress=progress, cancelled=lambda: self.cancelled)
        self.stage('dock', 'done', f'{plural(len(self.rows), "compound")} in '
                                   f'{(time.time() - t0) / 60:.1f} min')

        self.stage('rank', 'active')
        self.result = dict(summarise(self.rows, top=20), done=len(self.rows), total=total,
                           per_second=round(len(self.rows) / max(1e-6, time.time() - t0), 2),
                           two_state=self.b is not None, sites=self.sites,
                           parameters=self.parameters())
        failed = self.result.get('failed', 0)
        if failed:
            self.note(f'{plural(failed, "compound")} could not be placed in either state.', 'warning')
        self.note('These are rankings, not affinities: see the measured limits of the scoring function '
                  'before you order anything.', 'warning')
        self.stage('rank', 'done', f'{plural(len(self.rows), "compound")} ranked')


# ------------------------------------------------------------------------------------ exports

def to_csv(job):
    """The results as a CSV table, in the ranked order the window shows."""
    rows = job.table()['rows']
    two = job.b is not None
    head = ['rank', 'name', 'formula', 'weight', 'atoms', 'rotatable', 'charge', 'score_a']
    head += ['score_b', 'selectivity', 'worse'] if two else []
    head += ['strain_a'] + (['strain_b'] if two else []) + ['seconds', 'note']
    out = [','.join(head)]
    for i, r in enumerate(rows):
        cells = [i + 1, r.get('name', ''), r.get('formula', ''), r.get('weight', ''), r.get('atoms', ''),
                 r.get('rotatable', ''), r.get('charge', ''), r.get('score_a', '')]
        cells += [r.get('score_b', ''), r.get('selectivity', ''), r.get('worse', '')] if two else []
        cells += [r.get('strain_a', '')] + ([r.get('strain_b', '')] if two else [])
        cells += [r.get('seconds', ''), r.get('error', '')]
        out.append(','.join('' if c is None else f'"{c}"' if isinstance(c, str) and (',' in c or '"' in c)
                            else str(c) for c in cells))
    return '\n'.join(out) + '\n'


def to_sdf(job, state='a', limit=200):
    """The kept poses as an SDF, best first: what goes into a modelling program or a figure."""
    from rnafold3d.screening import chem
    mols = library_molecules(job.params['library'])
    rows = [r for r in job.table(poses=True)['rows'] if (r.get('poses') or {}).get(state)][:limit]
    out = []
    for rank, r in enumerate(rows):
        mol = mols[int(r['index'])]
        props = {'RCV_rank': rank + 1, 'RCV_score_a': r.get('score_a'), 'RCV_score_b': r.get('score_b'),
                 'RCV_selectivity': r.get('selectivity'), 'RCV_state': state.upper(),
                 'RCV_strain': r.get(f'strain_{state}'), 'RCV_note': 'a ranking score, not an affinity'}
        out.append((mol, np.asarray(r['poses'][state], float),
                    {k: v for k, v in props.items() if v is not None}))
    if not out:
        return ''
    return chem.write_sdf([m for m, _x, _p in out], xyz=[x for _m, x, _p in out],
                          properties=[p for _m, _x, p in out])


def pose_for(job, index, state='a'):
    """One pose for the main window to draw: its atoms, bonds and the nucleotides it touches."""
    from rnafold3d.screening.receptor import Receptor
    from scipy.spatial import cKDTree
    row = next((r for r in job.rows if int(r['index']) == int(index)), None)
    if row is None or not (row.get('poses') or {}).get(state):
        raise KeyError('no pose was kept for this compound: it is outside the number the run kept')
    mol = library_molecules(job.params['library'])[int(index)]
    xyz = np.asarray(row['poses'][state], float)
    k = job.a if state == 'a' else job.b
    doc = job.entry.doc
    coords = rcvfile.doc_coords(doc)
    rec = Receptor.from_document(doc, coords, k)
    near = cKDTree(xyz).query_ball_point(rec.xyz, 4.5, return_length=True)
    residues = sorted({int(rec.residue[i]) for i in np.nonzero(near)[0]})
    return {'name': row.get('name'), 'index': int(index), 'state': state, 'conformation': k,
            'elements': list(mol.elements), 'xyz': [[round(float(c), 3) for c in p] for p in xyz],
            'bonds': [[int(i), int(j)] for i, j, _o in mol.bonds], 'residues': residues,
            'score': row.get(f'score_{state}'), 'selectivity': row.get('selectivity'),
            'strain': row.get(f'strain_{state}')}


def library_from_path(path, ph=7.0, limit=MAX_LIBRARY):
    p = Path(path)
    if not p.is_file():
        raise ValueError(f'not a file: {path}')
    if p.suffix.lower() not in ('.sdf', '.mol', '.mol2', '.sd'):
        raise ValueError('a compound library is an SDF, MOL or MOL2 file')
    return parse_library(p.read_text(encoding='utf-8', errors='replace'), name=p.name, ph=ph, limit=limit)
