"""Null reference: C4' RMSD of random fold-tree arrangements with the correct (oracle) stems.

For each target, 30 random conformations are built from the oracle 2D with library torsions and
all loops built explicitly (no search, no compaction, chain-break/clash ignored). This gives the RMSD
expected from correct helices in a random arrangement; predictions should be judged against it.
Usage: PYTHONPATH=src python benchmark/null_baseline.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'src'))
from rnafold3d.fold import FoldModel                       # noqa: E402
from rnafold3d.predict import read_dbn                     # noqa: E402
from rnafold3d.annotate import dotbracket_to_pairs         # noqa: E402
from rnafold3d.metrics import read_numbered_pdb, kabsch_rmsd   # noqa: E402

out = {}
ds = json.loads((ROOT / 'data' / 'dataset_v1.json').read_text())
for s in sorted(ds['samples'], key=lambda s: s['length']):
    sdir = ROOT / 'data' / 'samples' / s['sample_id']
    seq = (sdir / 'sequence.fasta').read_text().split('\n')[1].strip()
    pairs = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(read_dbn(sdir / 'oracle_2d.dbn', len(seq)))]
    ref = read_numbered_pdb(sdir / 'reference.pdb')
    m = FoldModel(seq, pairs)
    m.set_loops_active(True)
    rng = np.random.default_rng(0)
    vals = []
    for _ in range(30):
        m.sample_torsions(rng)
        m.build()
        keys = [k for k in sorted(ref) if "C4'" in ref[k].atoms and m.present[m.idx[(k - 1, "C4'")]]]
        vals.append(kabsch_rmsd([ref[k].atoms["C4'"] for k in keys], [m.xyz[m.idx[(k - 1, "C4'")]] for k in keys]))
    out[s['sample_id']] = {'length': s['length'], 'median': round(float(np.median(vals)), 2),
                           'p10': round(float(np.percentile(vals, 10)), 2), 'min': round(float(min(vals)), 2)}
    print(s['sample_id'], out[s['sample_id']], flush=True)
(ROOT / 'runs' / 'null_baseline.json').write_text(json.dumps(out, indent=1))
