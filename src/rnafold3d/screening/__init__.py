"""Compound screening against a conformational ensemble.

    pockets.py    where a small molecule could bind one conformation, and how that site differs in another
    chem.py       a compound as the scoring function needs it: atoms, rings, rotatable bonds, charges
    receptor.py   the RNA as the scoring function needs it: typed atoms, phosphate charges, base planes
    score.py      the scoring function — contact and shape terms, stacking, electrostatics, desolvation
    grid.py       the same function precomputed on a grid, so a pose costs one lookup
    dock.py       the search: position, orientation and torsions, by Monte Carlo with local minimisation

A screen is what these add up to and it does not exist yet: no libraries, no two-state selectivity score,
no window in the interface (`docs/NEXT_STEPS.md` §1, steps 3 onwards). What is here is usable from the
command line — `rnafold3d pockets`, `rnafold3d dock` — and measured by `benchmark/screening/redock.py`,
with the numbers in `docs/SCREENING_STATUS.md`.

**Nothing here predicts an affinity.** The score ranks poses and compounds; it is not in kcal/mol, and the
weights of its RNA-specific terms are set rather than fitted.
"""
from .pockets import compare, find, lining_residues, volume_near

__all__ = ['find', 'compare', 'lining_residues', 'volume_near', 'chem', 'receptor', 'score', 'grid', 'dock']
