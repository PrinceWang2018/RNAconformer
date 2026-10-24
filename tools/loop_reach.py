"""Empirical reach of an n-nucleotide loop: distance O3'(a) -> C5'(b) when building
forward from residue a through n residues with torsions drawn from the calibrated
loop-conformer library. Stores percentiles in the geometry file (key 'loop_reach')."""
import json
import sys
from pathlib import Path

import numpy as np

from rnafold3d.build import Geometry, ResidueBuilder

path = Path(sys.argv[1] if len(sys.argv) > 1 else 'src/rnafold3d/data/geometry_v1.json')
raw = json.loads(path.read_text())
g = Geometry(raw)
rb = ResidueBuilder(g)
lib = g.loop_library
rng = np.random.default_rng(7)
names = ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi')


def sample():
    row = lib[rng.integers(len(lib))]
    t = dict(zip(names, row[:7]))
    return t, ('C2endo' if row[3] > 115 else 'C3endo')


out = {'n': [], 'p50': [], 'p95': [], 'p995': [], 'max': []}
for n in range(0, 41):
    d = []
    for _ in range(400 if n <= 12 else 150):
        prev = rb.unit('A', 'C3endo', -160.0)
        t_prev, _ = sample()
        for k in range(n + 1):
            t, puck = sample()
            P, O5, C5, C4, C3 = rb.forward_link(prev, t_prev['epsilon'], t_prev['zeta'], t['alpha'], t['beta'], t['gamma'], puck)
            if k == n:
                d.append(np.linalg.norm(C5 - rb.unit('A', 'C3endo', -160.0)["O3'"]))
                break
            prev = rb.place_unit_forward('A', puck, t['chi'], C5, C4, C3)
            t_prev = t
    d = np.array(d)
    out['n'].append(n)
    for k, q in (('p50', 50), ('p95', 95), ('p995', 99.5)):
        out[k].append(round(float(np.percentile(d, q)), 2))
    out['max'].append(round(float(d.max()), 2))
    print(n, out['p50'][-1], out['p95'][-1], out['p995'][-1], out['max'][-1])
raw['loop_reach'] = {'definition': "distance O3'(a) to C5'(b) across n unpaired residues, library torsions",
                     **out}
path.write_text(json.dumps(raw, indent=1))
