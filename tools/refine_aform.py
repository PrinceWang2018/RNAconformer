"""Refine the calibrated A-form torsion set so that long ideal duplexes are self-consistent.

Minimizes the Watson-Crick pair-fit RMSD of a 14-bp mixed-sequence duplex built with
uniform torsions, with a penalty (in units of the calibrated SD) for departing from the
calibrated helix means. Writes 'aform_torsions_refined' into the geometry file.
Uses only calibration-derived quantities (no benchmark structures).
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from rnafold3d.build import ResidueBuilder, Geometry, build_stem

path = Path(sys.argv[1] if len(sys.argv) > 1 else 'src/rnafold3d/data/geometry_v1.json')
raw = json.loads(path.read_text())
names = ('alpha', 'beta', 'gamma', 'epsilon', 'zeta', 'chi')
mean = np.array([raw['aform_torsions'][k] for k in names])
sd = np.array([max(raw['aform_torsion_sd'][k], 3.0) for k in names])
seq = 'GCAUGGCUACGUCA' + 'UGACGUAGCCAUGC'
pairs = [(k, 27 - k) for k in range(14)]


def objective(x):
    raw2 = dict(raw)
    raw2['aform_torsions'] = {**raw['aform_torsions'], **{k: float(v) for k, v in zip(names, x)}}
    rb = ResidueBuilder(Geometry(raw2))
    _, rep = build_stem(rb, seq, pairs)
    return rep['pair_fit_rmsd'] ** 2 + 0.05 * np.sum(((x - mean) / sd) ** 2), rep['pair_fit_rmsd']


f0, r0 = objective(mean)
res = minimize(lambda x: objective(x)[0], mean, method='Nelder-Mead',
               options={'maxiter': 600, 'xatol': 0.05, 'fatol': 1e-5})
f1, r1 = objective(res.x)
print('pair-fit RMSD 14 bp: before %.3f  after %.3f' % (r0, r1))
print('shift (deg):', dict(zip(names, np.round(res.x - mean, 2))))
raw['aform_torsions_refined'] = {**raw['aform_torsions'], **{k: round(float(v), 2) for k, v in zip(names, res.x)}}
raw['aform_refinement'] = {'method': 'Nelder-Mead on 14-bp duplex pair-fit RMSD + 0.05*sum(((x-mean)/sd)^2)',
                           'pair_fit_rmsd_before': round(r0, 3), 'pair_fit_rmsd_after': round(r1, 3)}
path.write_text(json.dumps(raw, indent=1))
