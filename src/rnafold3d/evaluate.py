"""Conformation evaluator: comparable energy terms for ANY all-atom RNA model of one sequence.

The evaluator does not depend on how a model was produced (rnafold3d engine, NMR ensemble,
another predictor), so a viewer can show how the energy changes between conformations.
Every term is documented here and in docs/RCV.md; none of them is a free energy of the 3D state.

Terms (kcal/mol-scaled):
    dg2d     Turner-2004 nearest-neighbour free energy of the nested canonical base pairs annotated in
             the 3D coordinates (ViennaRNA eval_structure, 37 C). Thermodynamic, covers 2D only.
    elec     Debye-Hueckel phosphate-phosphate repulsion: q_eff = -0.5 e per phosphate (counterion
             condensation), relative permittivity 78.5, Debye length 7.9 A (~150 mM monovalent salt).
    stack3d  -1.0 per base-base stack that is not a canonical helix step (loop and tertiary stacking,
             which dg2d does not count). Geometric detection, heuristic weight.
    steric   +1.0 per severe heavy-atom clash (< 2.2 A, more than 3 bonds apart; rnafold3d.checks).
    total    dg2d + elec + stack3d + steric  ("evaluator energy"; heuristic composite).
Descriptors: n_pairs, n_pk_pairs, n_stacks, rg (heavy atoms), overlaps_per_1000, valid.
"""
from __future__ import annotations

import numpy as np

from .annotate import base_plane, canonical_pairs, nested_subset, pairs_to_dotbracket
from .build import atom_order
from .checks import validate_model

COULOMB = 332.0637            # kcal mol^-1 A e^-2
Q_EFF = -0.5
EPSILON = 78.5
DEBYE_LENGTH = 7.9
STACK_WEIGHT = -1.0
CLASH_WEIGHT = 1.0

TERMS = {
    'dg2d': {'label': 'ΔG 2D (Turner)', 'unit': 'kcal/mol',
             'description': 'ViennaRNA nearest-neighbour free energy of the nested canonical pairs found in '
                            'this model (37 °C). Thermodynamic estimate of the secondary structure only.'},
    'elec': {'label': 'Electrostatics', 'unit': 'kcal/mol',
             'description': 'Debye–Hückel phosphate repulsion; q = −0.5 e, ε = 78.5, Debye length 7.9 Å '
                            '(~150 mM monovalent salt, no Mg²⁺).'},
    'stack3d': {'label': 'Loop/tertiary stacking', 'unit': 'kcal/mol',
                'description': '−1.0 per base–base stack outside canonical helix steps (heuristic weight).'},
    'steric': {'label': 'Steric clashes', 'unit': 'kcal/mol',
               'description': '+1.0 per severe heavy-atom clash (< 2.2 Å, > 3 bonds apart).'},
}
TOTAL = {'key': 'eval_total', 'label': 'Evaluator energy', 'unit': 'kcal/mol', 'lower_is_better': True,
         'terms': ['dg2d', 'elec', 'stack3d', 'steric'],
         'description': 'Heuristic composite: ΔG 2D + electrostatics + loop/tertiary stacking + clash penalty. '
                        'Comparable between conformations of one sequence; not a free energy of the 3D state.'}
DESCRIPTORS = {
    'n_pairs': {'label': 'Canonical base pairs', 'unit': '', 'lower_is_better': False},
    'n_pk_pairs': {'label': 'Pseudoknotted pairs', 'unit': '', 'lower_is_better': None},
    'n_stacks': {'label': 'Base stacks', 'unit': '', 'lower_is_better': False},
    'rg': {'label': 'Radius of gyration', 'unit': 'Å', 'lower_is_better': None},
    'overlaps_per_1000': {'label': 'Atom overlaps / 1000 atoms', 'unit': '', 'lower_is_better': True},
}


class _Res:
    __slots__ = ('index', 'base', 'atoms')

    def __init__(self, index, base, atoms):
        self.index, self.base, self.atoms = index, base, atoms


def annotate_pairs(sequence, residues):
    """Canonical pairs of a model. residues: list (per position) of {atom: xyz}.

    Returns (pairs, nested, pseudoknotted, dot_bracket); pairs are 1-based (i, j, type).
    """
    objs = [_Res(k + 1, sequence[k], residues[k]) for k in range(len(sequence)) if residues[k]]
    pairs = canonical_pairs(objs)
    nested, removed = nested_subset([(i, j) for i, j, _ in pairs])
    dbn = pairs_to_dotbracket(len(sequence), nested + removed) if pairs else '.' * len(sequence)
    return pairs, nested, removed, dbn


def find_stacks(sequence, residues):
    """Base-base stacks from ring geometry: list of 1-based (i, j)."""
    planes = []
    for k, b in enumerate(sequence):
        p = base_plane(residues[k], b) if residues[k] and b in 'ACGU' else None
        planes.append(p)
    idx = [k for k, p in enumerate(planes) if p is not None]
    if len(idx) < 2:
        return []
    cen = np.array([planes[k][0] for k in idx])
    nor = np.array([planes[k][1] for k in idx])
    out = []
    for a in range(len(idx)):
        d = cen[a + 1:] - cen[a]
        dist = np.linalg.norm(d, axis=1)
        for off in np.nonzero(dist < 5.5)[0]:
            b = a + 1 + off
            cosang = abs(float(nor[a] @ nor[b]))
            if cosang < np.cos(np.radians(35.0)):
                continue
            v1 = abs(float(d[off] @ nor[a]))
            v2 = abs(float(d[off] @ nor[b]))
            vert = 0.5 * (v1 + v2)
            lat = np.sqrt(max(0.0, dist[off] ** 2 - min(v1, v2) ** 2))
            if 2.9 <= vert <= 4.3 and lat <= 2.8:
                out.append((idx[a] + 1, idx[b] + 1))
    return out


def _vienna_eval(sequence, dbn):
    try:
        import RNA
    except ImportError:
        return None
    try:
        return float(RNA.fold_compound(sequence).eval_structure(dbn))
    except Exception:          # malformed / unsupported structure: report as unavailable
        return None


def electrostatics(residues):
    p = np.array([r['P'] for r in residues if r and 'P' in r])
    if len(p) < 2:
        return 0.0
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    iu = np.triu_indices(len(p), 1)
    r = np.maximum(d[iu], 2.0)
    return float(COULOMB * Q_EFF * Q_EFF / EPSILON * np.sum(np.exp(-r / DEBYE_LENGTH) / r))


def radius_of_gyration(residues):
    x = np.array([v for r in residues if r for v in r.values()])
    if len(x) == 0:
        return 0.0
    return float(np.sqrt(np.mean(np.sum((x - x.mean(axis=0)) ** 2, axis=1))))


def evaluate_many(sequence, residue_sets, validations=None, workers=1):
    """Evaluate several conformations, in parallel worker processes when `workers` > 1."""
    validations = validations or [None] * len(residue_sets)
    if workers <= 1 or len(residue_sets) < 2:
        return [evaluate(sequence, r, v) for r, v in zip(residue_sets, validations)]
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    ctx = mp.get_context('spawn')            # safe from threaded servers; workers re-import the package
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
        return list(pool.map(_evaluate_args, [(sequence, r, v) for r, v in zip(residue_sets, validations)]))


def _evaluate_args(args):
    return evaluate(*args)


def evaluate(sequence, residues, validation=None):
    """Evaluate one conformation. residues: list (per sequence position, may be empty dicts) of {atom: xyz}.

    `validation` may pass a precomputed rnafold3d.checks report to avoid recomputing it.
    Returns {'values': {term: value}, 'pairs': [...], 'pk_pairs': [...], 'stacks': [...], 'dbn': str,
             'validation': summary}.
    """
    seq = sequence.upper().replace('T', 'U')
    std = [{k: v for k, v in (r or {}).items() if b in 'ACGU' and k in atom_order(b)} for b, r in zip(seq, residues)]
    pairs, nested, removed, dbn = annotate_pairs(seq, std)
    nested_dbn = pairs_to_dotbracket(len(seq), nested) if nested else '.' * len(seq)
    dg2d = _vienna_eval(seq, nested_dbn) if set(seq) <= set('ACGU') else None
    partner = {}
    for i, j in nested:
        partner[i], partner[j] = j, i
    stacks = find_stacks(seq, std)
    helix_steps = 0
    loop_stacks = 0
    for i, j in stacks:
        if j == i + 1 and i in partner and j in partner and partner[j] == partner[i] - 1:
            helix_steps += 1
        elif not (i in partner and partner[i] == j):
            loop_stacks += 1
    if validation is None:
        validation = validate_model(seq, std, max_report=20) if all(std) else None
    severe = (validation.get('severe_clashes') or 0) if validation else 0
    overlaps = None
    if validation:
        overlaps = validation.get('heavy_atom_overlaps_per_1000_atoms', validation.get('overlaps_per_1000'))
    values = {
        'dg2d': None if dg2d is None else round(dg2d, 2),
        'elec': round(electrostatics(std), 2),
        'stack3d': round(STACK_WEIGHT * loop_stacks, 2),
        'steric': round(CLASH_WEIGHT * severe, 2),
        'n_pairs': len(pairs),
        'n_pk_pairs': len(removed),
        'n_stacks': len(stacks),
        'rg': round(radius_of_gyration(std), 2),
        'overlaps_per_1000': overlaps,
    }
    values['eval_total'] = round(sum(values[k] or 0.0 for k in TOTAL['terms']), 2)
    summary = None
    if validation is not None:
        summary = {'passed': validation.get('passed'), 'issue_counts': validation.get('issue_counts', {}),
                   'severe_clashes': severe, 'overlaps_per_1000': overlaps,
                   'sugar_puckers': validation.get('sugar_puckers')}
    return {'values': values, 'pairs': [[i, j, t] for i, j, t in pairs],
            'pk_pairs': [[i, j] for i, j in removed], 'stacks': [[i, j] for i, j in stacks],
            'dbn': dbn, 'validation': summary, 'helix_steps': helix_steps}
