"""Build the knowledge-based base-base potential (src/rnafold3d/data/basepot_v1.npz) from the structure library.

Counted: every library chain (all released before 2024-05-01, the library's split date, so no CASP16 target is
in it) that is not homologous to a development-benchmark target (benchmark/data/samples; homology as for the
junction geometry: a significant Smith-Waterman hit, E <= 0.01). See rnafold3d/basepot.py for the features.

    python tools/build_base_potential.py
"""
from __future__ import annotations

import datetime
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))


def main():
    from rnafold3d import basepot as BP
    from rnafold3d.templates import load_library, homologs, sw_scan
    from rnafold3d.motifs import nested_greedy
    lib = load_library()
    excluded = np.zeros(lib.n_chains, bool)
    dev = []
    for d in sorted((ROOT / 'benchmark' / 'data' / 'samples').iterdir()):
        seq = ''.join(l.strip() for l in (d / 'sequence.fasta').read_text().splitlines() if not l.startswith('>'))
        dev.append(d.name)
        excluded |= homologs(seq, lib, sw_scan(seq, lib))
    counts = np.zeros((BP.N_TYPES, 2, BP.N_D, BP.H_BINS), np.int64)
    used = 0
    for k in range(lib.n_chains):
        if excluded[k]:
            continue
        seq = lib.chain_seq(k)
        if len(seq) < 20:
            continue
        atoms = {nm: lib.chain_xyz(k, nm) for nm in ('C2', 'C4', 'C6')}
        cen, nrm = BP.base_frames(seq, atoms)
        hid = BP.helix_id(len(seq), nested_greedy(lib.chain_pairs(k)))
        t, s, d, h = BP.pair_features(seq, cen, nrm, hid)
        np.add.at(counts, (t, s, d, h), 1)
        used += 1
    energy = BP.potential_from_counts(counts)
    info = {'created': datetime.date.today().isoformat(), 'library': lib.version, 'chains_counted': used,
            'chains_excluded_homologous_to_dev': int(excluded.sum()), 'dev_targets': dev,
            'pairs_counted': int(counts.sum()), 'd_min': BP.D_MIN, 'd_max': BP.D_MAX, 'd_step': BP.D_STEP,
            'h_bins': BP.H_BINS, 'h_cut': BP.H_CUT, 'min_sep': BP.MIN_SEP}
    np.savez_compressed(BP.DATA, energy=energy.astype(np.float32), counts=counts.astype(np.int32),
                        info=np.frombuffer(json.dumps(info).encode(), np.uint8))
    print(json.dumps(info, indent=1))


if __name__ == '__main__':
    main()
