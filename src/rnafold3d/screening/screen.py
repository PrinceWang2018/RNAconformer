"""The two-state screen: which compounds prefer conformation A over conformation B.

This is the part of the screen that is not docking. Docking one compound into one pocket is `dock.py`;
what a conformational switch asks for is the *difference* between two pockets of the same RNA:

    selectivity = score(B) − score(A)

positive when the compound fits A better, which is what stabilises A and blocks the transition to B. The
number is a **ranking**, not an affinity, and not a probability that anything will happen at the bench: it
picks a few hundred compounds out of a hundred thousand to try. Two compounds with the same selectivity are
ranked by the worse of their two scores, because a compound that binds neither state selectively is of no
use — that is `sort_key` below, and it is the order every list in RCV shows.

What the score cannot see, and what the caller has to keep in view:

* The RNA is rigid in both states. A compound that would open a pocket that neither conformation has is
  invisible here, and one that fits a pocket only because the model has it slightly too wide scores too
  well. The two-state difference cancels some of this, since both states come from the same engine.
* `strain` is the compound's own internal cost in the kept pose. A good score with a high strain is a pose
  holding the compound in a shape it would not take by itself; the table shows it for that reason.
* A compound that cannot be placed at all (too big for the box, no pose found) gets `None` for that state,
  and is reported rather than dropped: a compound that fits A and *cannot* be placed in B at all is the
  strongest selectivity there is, and is also the most likely to be an artefact of a too-small box.

The weights of the scoring function are the ones measured in `docs/SCREENING_STATUS.md`, and under the
plan's rule they may only ever be tuned against the enrichment benchmark, never against the selectivity
cases this file produces.
"""
from __future__ import annotations

import math
import multiprocessing as mp
import os
import time
from dataclasses import dataclass, field

import numpy as np

from . import dock as D
from . import grid as G
from . import score as S

PAD = 6.0                       # Å of box around a pocket: room for a compound to lie across its mouth
MAX_ATOMS = 120                 # a compound bigger than this is not a small molecule; it is skipped
MAX_SIZE = 30.0                 # Å, the widest box a screen will build a grid for: a 30 Å box is already
                                # 80³ points per channel (about 65 MB of grid), and the search cost
                                # grows with the volume. The
                                # open clefts of an RNA can form one pocket spanning the whole molecule
                                # (`pockets.find` marks those `wide`); screening it as one box is both slow
                                # and pointless, so the box is cropped to the centre of the pocket and says
                                # so in `cropped`, for the interface to warn about.


@dataclass
class Site:
    """Where to dock: the centre and the edge of the box, with what it was derived from."""
    centre: tuple
    size: float
    label: str = ''
    residues: tuple = ()
    cropped: bool = False       # the pocket was wider than MAX_SIZE and only its centre is screened

    def __post_init__(self):
        if self.size > MAX_SIZE:
            self.size, self.cropped = MAX_SIZE, True

    @classmethod
    def from_pocket(cls, pocket, pad=PAD):
        """The box around a pocket of `pockets.find`: its own bounding box, or its points if it kept them."""
        if 'points' in pocket and len(pocket['points']):
            return cls.from_atoms(pocket['points'], pad, label=str(pocket.get('id', '')),
                                  residues=pocket.get('residues', ()))
        lo, hi = np.asarray(pocket['min'], float), np.asarray(pocket['max'], float)
        return cls(centre=tuple((lo + hi) / 2), size=float((hi - lo).max()) + 2 * pad,
                   label=str(pocket.get('id', '')), residues=tuple(pocket.get('residues', ())))

    @classmethod
    def from_atoms(cls, xyz, pad=PAD, label='', residues=()):
        xyz = np.asarray(xyz, float).reshape(-1, 3)
        lo, hi = xyz.min(0), xyz.max(0)
        return cls(centre=tuple((lo + hi) / 2), size=float((hi - lo).max()) + 2 * pad, label=label,
                   residues=tuple(residues))

    def as_dict(self):
        return {'centre': [round(float(c), 2) for c in self.centre], 'size': round(float(self.size), 1),
                'label': self.label, 'residues': list(self.residues), 'cropped': bool(self.cropped)}


@dataclass
class Result:
    """One compound against both states."""
    name: str = ''
    index: int = 0
    atoms: int = 0
    rotatable: int = 0
    charge: float = 0.0
    formula: str = ''
    weight: float = 0.0
    score_a: float = None
    score_b: float = None
    strain_a: float = 0.0
    strain_b: float = 0.0
    seconds: float = 0.0
    error: str = ''
    poses: dict = field(default_factory=dict)
    properties: dict = field(default_factory=dict)

    @property
    def selectivity(self):
        """score(B) − score(A): positive when the compound prefers state A. None if a state has no pose."""
        if self.score_a is None or self.score_b is None:
            return None
        return self.score_b - self.score_a

    @property
    def worse(self):
        """The weaker of the two scores (the one that decides whether it binds anything at all)."""
        have = [s for s in (self.score_a, self.score_b) if s is not None]
        return max(have) if have else None

    def as_dict(self, poses=False):
        out = {'name': self.name, 'index': self.index, 'atoms': self.atoms, 'rotatable': self.rotatable,
               'charge': round(self.charge, 2), 'formula': self.formula, 'weight': round(self.weight, 1),
               'score_a': _round(self.score_a), 'score_b': _round(self.score_b),
               'selectivity': _round(self.selectivity), 'worse': _round(self.worse),
               'strain_a': round(self.strain_a, 2), 'strain_b': round(self.strain_b, 2),
               'seconds': round(self.seconds, 2), 'error': self.error}
        if self.properties:
            out['properties'] = self.properties
        if poses and self.poses:
            out['poses'] = {k: [[round(float(c), 3) for c in p] for p in v]
                            for k, v in self.poses.items()}
        return out


def _round(x, places=2):
    return None if x is None else round(float(x), places)


def sort_key(row):
    """The order the plan asks for: most selective first, ties broken by the worse of the two scores.

    Works on a `Result` or on the dict it produces. A compound that could not be placed in one of the two
    states has no selectivity and sorts after every compound that has one, however well it scored.
    """
    sel = row.selectivity if isinstance(row, Result) else row.get('selectivity')
    worse = row.worse if isinstance(row, Result) else row.get('worse')
    return (0 if sel is not None else 1, -(sel or 0.0), worse if worse is not None else math.inf)


class Screen:
    """Two receptors, two boxes, one scoring function. Build it, `prepare()` it, then `one()` per compound.

    Both states are usually two conformations of the same molecule, so they have the same atoms in
    different places; nothing here requires that, and a single state (``receptor_b=None``) gives a plain
    ranking by score with no selectivity, which is what the enrichment benchmark and a single-pocket screen
    need.
    """

    def __init__(self, receptor_a, site_a, receptor_b=None, site_b=None, exhaustiveness='quick',
                 seed=1, n_poses=3, weights=None, spacing=0.375):
        self.rec = {'a': receptor_a, 'b': receptor_b}
        self.site = {'a': site_a, 'b': site_b}
        self.exhaustiveness = exhaustiveness
        self.seed = int(seed)
        self.n_poses = int(n_poses)
        self.weights = dict(S.WEIGHTS, **(weights or {}))
        self.spacing = float(spacing)
        self.grid = {}

    @property
    def states(self):
        return [k for k in ('a', 'b') if self.rec[k] is not None and self.site[k] is not None]

    def prepare(self):
        """Build the grid of every state (seconds, and the memory the screen then holds)."""
        for k in self.states:
            if k not in self.grid:
                self.grid[k] = G.Grid.for_pocket(self.rec[k], centre=self.site[k].centre,
                                                 size=self.site[k].size, spacing=self.spacing)
        return self

    @property
    def nbytes(self):
        return sum(g.nbytes for g in self.grid.values())

    def one(self, mol, index=0, keep_poses=False):
        """Dock one compound into every state. Never raises: a failure lands in `Result.error`."""
        t0 = time.time()
        row = Result(name=mol.name or f'compound {index + 1}', index=index, atoms=len(mol),
                     rotatable=mol.n_rotatable, charge=float(mol.charge.sum()), formula=mol.formula,
                     weight=_weight(mol), properties=dict(mol.properties or {}))
        row.properties.pop('error', None)
        if not len(mol):
            row.error = 'no atoms'
            return row
        if len(mol) > MAX_ATOMS:
            row.error = f'{len(mol)} atoms: too big for a small-molecule screen'
            return row
        self.prepare()
        for k in self.states:
            try:
                poses = D.dock(mol, self.grid[k], receptor=self.rec[k], weights=self.weights,
                               exhaustiveness=self.exhaustiveness, seed=self.seed, n_poses=self.n_poses)
            except Exception as err:                       # noqa: BLE001 - one compound must not stop a run
                row.error = f'{type(err).__name__}: {err}'
                continue
            if not poses:
                continue
            setattr(row, f'score_{k}', float(poses[0].score))
            setattr(row, f'strain_{k}', float(S.strain(poses[0].xyz, mol)))
            if keep_poses:
                row.poses[k] = poses[0].xyz.tolist()
        if row.score_a is None and row.score_b is None and not row.error:
            row.error = 'no pose fits the box'
        row.seconds = time.time() - t0
        return row


# --------------------------------------------------------------------------------- running a library

_WORK = {}


def _init(args):
    """Worker start-up: build the screen here (the spawn path, Windows) or inherit it (fork)."""
    if args is not None:
        _WORK['screen'] = Screen(*args[0], **args[1]).prepare()


def _job(task):
    index, mol, keep = task
    return _WORK['screen'].one(mol, index, keep).as_dict(poses=keep)


def _fork_ok():
    return 'fork' in mp.get_all_start_methods() and os.name != 'nt'


def workers(n=None, per_worker_mb=None):
    """How many processes to use: what was asked for, what the machine has, and what memory allows."""
    have = os.cpu_count() or 1
    n = have if n in (None, 0) else int(n)
    n = max(1, min(n, have))
    if per_worker_mb and not _fork_ok():
        # Without fork (Windows) every worker builds and holds its own copy of the grids, so the limit is
        # memory rather than cores. Where the free memory cannot be read, assume a modest 2 GB is free.
        free_mb = 2048.0
        try:
            import psutil                                   # optional dependency
            free_mb = psutil.virtual_memory().available / 1e6
        except Exception:                                   # noqa: BLE001 - the default above stands
            pass
        n = max(1, min(n, int(free_mb * 0.6 / per_worker_mb)))
    return n


def run(screen, molecules, jobs=1, keep_poses=0, progress=None, cancelled=None):
    """Screen a library. Yields nothing; returns the rows, in the order the compounds came in.

    `keep_poses` keeps the coordinates of the best pose of each state for the first *n* compounds **in the
    order they arrive**, which is what a run that cannot hold a hundred thousand poses in memory does; the
    caller re-docks the short list afterwards if it wants poses for the compounds that came out on top.
    `progress(done, total, row)` is called after every compound, and `cancelled()` stops the run.
    """
    molecules = list(molecules)
    total = len(molecules)
    tasks = [(k, m, k < keep_poses) for k, m in enumerate(molecules)]
    rows = []
    if jobs <= 1:
        screen.prepare()
        for task in tasks:
            if cancelled is not None and cancelled():
                break
            rows.append(screen.one(task[1], task[0], task[2]).as_dict(poses=task[2]))
            if progress:
                progress(len(rows), total, rows[-1])
        return rows
    if _fork_ok():
        # Fork, so the grids are built once and the workers read the parent's copy — the difference between
        # 60 MB and 60 MB per core on a laptop. Python warns that forking a process with threads in it (the
        # viewer runs an HTTP server) can deadlock; the workers here allocate and compute and touch nothing
        # the parent's threads hold, which is the case the warning does not apply to. Windows has no fork
        # and takes the branch below, where each worker builds its own grids.
        screen.prepare()
        _WORK['screen'] = screen
        ctx, args = mp.get_context('fork'), None
    else:
        ctx = mp.get_context('spawn')
        args = ((screen.rec['a'], screen.site['a'], screen.rec['b'], screen.site['b']),
                {'exhaustiveness': screen.exhaustiveness, 'seed': screen.seed, 'n_poses': screen.n_poses,
                 'weights': screen.weights, 'spacing': screen.spacing})
    with ctx.Pool(jobs, initializer=_init, initargs=(args,)) as pool:
        try:
            for row in pool.imap(_job, tasks, chunksize=1):
                rows.append(row)
                if progress:
                    progress(len(rows), total, row)
                if cancelled is not None and cancelled():
                    pool.terminate()
                    break
        finally:
            _WORK.pop('screen', None)
    return rows


def _weight(mol):
    """Molecular weight of the heavy atoms plus the hydrogens counted onto them."""
    mass = {'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'P': 30.974, 'S': 32.06,
            'CL': 35.45, 'BR': 79.904, 'I': 126.904, 'B': 10.81, 'SI': 28.085, 'SE': 78.97}
    total = float(np.sum(mol.hydrogens) * mass['H'])
    for e in mol.elements:
        total += mass.get(str(e).upper(), 0.0)
    return total


def summarise(rows, top=10):
    """What a report says about a finished run: counts, the spread of the scores, and the short list."""
    ok = [r for r in rows if r.get('selectivity') is not None]
    one = [r for r in rows if r.get('selectivity') is None and r.get('worse') is not None]
    sel = np.array([r['selectivity'] for r in ok], float) if ok else np.zeros(0)
    out = {'compounds': len(rows), 'scored': len(ok) + len(one), 'one_state_only': len(one),
           'failed': sum(1 for r in rows if r.get('error')),
           'seconds': round(sum(r.get('seconds', 0.0) for r in rows), 1),
           'top': sorted([r for r in rows if r.get('worse') is not None], key=sort_key)[:top]}
    if len(sel):
        out.update(selectivity_median=round(float(np.median(sel)), 2),
                   selectivity_max=round(float(sel.max()), 2),
                   selectivity_min=round(float(sel.min()), 2))
    return out
