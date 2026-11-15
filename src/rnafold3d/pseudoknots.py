"""Pseudoknots in the search: a 2D hypothesis with one crossing helix, folded jointly with the nested part.

ViennaRNA folds nested structures only. Here every helix of >= MIN_BP stacked canonical pairs that is stable on
its own (Turner 2004 duplex energy, initiation included) is tried as a pseudoknot helix: its residues are made
unavailable, the rest is folded by ViennaRNA (MFE with those residues unpaired), and the helix is kept if it
crosses the resulting structure and

    E(nested part with the helix residues unpaired) + E(helix) + PENALTY  <  E(nested MFE).

This is the shape of the simple pseudoknot models (one crossing helix, a constant initiation penalty, the
loops it spans not further scored), enough for H-type pseudoknots, kissing loops and loop-to-strand pairs. The
penalty was calibrated on the library structures (benchmark/pseudoknot_calibration.py): over 1,509 chains of
40-200 nt it maximises the base-pair F1 against the pairs annotated in 3D (MFE 0.624, with pseudoknots 0.631;
139 chains change, 75 gain correct pairs, 23 lose some; penalty 3: 0.628, 5: 0.629, 6: 0.629).

In 3D the gain is small: the crossing helix is often predicted correctly but rarely forms in the models, and on
benchmark/dev_templates.py the hypothesis adds nothing on top of the motif fragments while slowing the runs that
have one. It is therefore off by default (PredictConfig.pseudoknots = 'auto' or --pseudoknots auto to use it).

The crossing pairs cannot be part of the fold tree: they become base-pair restraints in the 3D search (C1'-C1'
distance and the central hydrogen bond), as crossing pairs from templates already do. The nested hypothesis
stays as well, so a wrong pseudoknot costs trajectories, not the prediction.
"""
from __future__ import annotations

MIN_BP = 4
PENALTY = 4.0               # kcal/mol: pseudoknot initiation (calibrated, see module doc)
MIN_SPAN = 12               # a crossing helix spans at least this many residues
CANON = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}


def helix_candidates(seq, min_bp=MIN_BP, min_span=MIN_SPAN, max_energy=-PENALTY):
    """Maximal stacked runs of canonical pairs (>= min_bp), or their longest 3'-trimmed version, whose duplex
    energy is at most max_energy: [(energy, [(i, j)...])], most stable first."""
    import RNA
    n = len(seq)
    can = lambda i, j: 0 <= i < j < n and seq[i] + seq[j] in CANON
    out = []
    for i in range(n):
        for j in range(n - 1, i + min_span, -1):
            if not can(i, j) or can(i - 1, j + 1):
                continue
            run = []
            a, b = i, j
            while can(a, b) and b - a > 3:
                run.append((a, b))
                a, b = a + 1, b - 1
            for end in range(len(run), min_bp - 1, -1):
                sub = run[:end]
                (a0, b0), (a1, b1) = sub[0], sub[-1]
                e = RNA.fold_compound(seq[a0:a1 + 1] + '&' + seq[b1:b0 + 1]).eval_structure('(' * end + ')' * end)
                if e <= max_energy:
                    out.append((round(float(e), 2), sub))
                    break
    out.sort(key=lambda x: (x[0], x[1][0]))
    return out


def _crosses(helix, nested):
    a, c = helix[0]
    return any(a < i < c < j or i < a < j < c for i, j in nested)


def fold_with_pseudoknot(seq, penalty=PENALTY, temperature=37.0, no_lonely_pairs=True, max_candidates=150,
                         reactivities=None):
    """The best structure with one crossing helix if it beats the nested MFE once `penalty` is paid:
    {'nested': [(i, j)], 'pk': [(i, j)], 'energy', 'mfe', 'helix_energy', 'nested_energy'} (kcal/mol), or None.

    reactivities: SHAPE pseudo-energies as in ss.predict_structures (applied to the nested part)."""
    import RNA
    from .annotate import dotbracket_to_pairs
    md = RNA.md()
    md.temperature = float(temperature)
    md.noLP = 1 if no_lonely_pairs else 0

    def compound():
        fc = RNA.fold_compound(seq, md)
        if reactivities is not None:
            from .ss import SHAPE_SLOPE, SHAPE_INTERCEPT
            fc.sc_add_SHAPE_deigan([-999.0] + [(-999.0 if v is None else max(0.0, float(v))) for v in reactivities],
                                   SHAPE_SLOPE, SHAPE_INTERCEPT)
        return fc
    _s, mfe_e = compound().mfe()
    best = None
    for e_h, helix in helix_candidates(seq, max_energy=-penalty)[:max_candidates]:
        if best is not None and mfe_e + e_h + penalty >= best['energy']:
            break                                     # sorted by helix energy: no later helix can win
        fc = compound()
        for i, j in helix:
            fc.hc_add_up(i + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
            fc.hc_add_up(j + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        s, e = fc.mfe()
        nested = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(s)]
        if not _crosses(helix, nested):
            continue                                  # a nested alternative, not a pseudoknot
        total = e + e_h + penalty
        if total < mfe_e and (best is None or total < best['energy']):
            best = {'nested': nested, 'pk': list(helix), 'energy': round(total, 2), 'mfe': round(mfe_e, 2),
                    'helix_energy': e_h, 'nested_energy': round(e, 2)}
    return best


def hypothesis(seq, temperature=37.0, reactivities=None):
    """A 2D hypothesis (as in ss.predict_structures) with a pseudoknot, or None."""
    from .annotate import pairs_to_dotbracket
    r = fold_with_pseudoknot(seq, temperature=temperature, reactivities=reactivities)
    if r is None:
        return None
    dbn = list(pairs_to_dotbracket(len(seq), r['nested'], offset=0))
    for i, j in r['pk']:
        dbn[i], dbn[j] = '[', ']'
    return {'id': '', 'dbn': ''.join(dbn), 'pairs': r['nested'], 'energy_kcal_mol': r['energy'],
            'source': 'vienna_pseudoknot', 'restraints': r['pk'],
            'pseudoknot_pairs': [(i + 1, j + 1) for i, j in r['pk']], 'pseudoknot_treatment': 'base-pair restraints',
            'pseudoknot': {'helix_energy_kcal_mol': r['helix_energy'], 'nested_energy_kcal_mol': r['nested_energy'],
                           'penalty_kcal_mol': PENALTY, 'mfe_kcal_mol': r['mfe']}}


def helix_restraints(seq, pk_pairs, weight=1.0):
    """C4'-C4' distances of an ideal A-form stem for the crossing helix, as distance restraints
    ({'i','j','d','tol','w'}, the format of template restraints): they pull the two strands into a helix in the
    least-squares fit and in the restrained loop closure, before the base-pair restraints set the pairs."""
    from .build import ResidueBuilder, load_geometry, build_stem
    import numpy as np
    runs = []
    for i, j in sorted(pk_pairs):
        if runs and runs[-1][-1] == (i - 1, j + 1):
            runs[-1].append((i, j))
        else:
            runs.append([(i, j)])
    out = {'i': [], 'j': [], 'd': [], 'tol': [], 'w': []}
    rb = ResidueBuilder(load_geometry())
    for run in runs:
        if len(run) < 2:
            continue
        atoms, _rep = build_stem(rb, seq, run, refine=False)
        res = sorted(atoms)
        X = np.array([atoms[r]["C4'"] for r in res])
        for a in range(len(res)):
            for b in range(a + 1, len(res)):
                d = float(np.linalg.norm(X[a] - X[b]))
                out['i'].append(res[a])
                out['j'].append(res[b])
                out['d'].append(round(d, 2))
                out['tol'].append(round(0.5 + 0.03 * d, 2))
                out['w'].append(weight)
    return out
