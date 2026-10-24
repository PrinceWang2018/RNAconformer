"""Ensemble packages: the stable data contract between prediction and the future viewer (WP5).

A package is a directory with
    metadata/manifest.json      (schema 'rnafold3d.ensemble/1', see schemas/ensemble_manifest.schema.json)
    conformations/conf_###.pdb  (+ .cif)  chain A, residues numbered 1..N along the sequence
Packages are produced by `rnafold3d predict` or by `import_ensemble` (external PDB/mmCIF files),
so a viewer never needs to know which engine (if any) produced the conformers.
"""
from __future__ import annotations

import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

from .build import atom_order
from .checks import validate_model
from .metrics import kabsch_rmsd
from .structure_io import read_structure, write_pdb, write_mmcif

SCHEMA = 'rnafold3d.ensemble/1'


def _chain_to_residues(chain, sequence):
    """Map an RNA chain onto 1..N of `sequence`; returns residues list and mapping report."""
    rm = chain.residue_map()
    seq_chain = chain.sequence
    if seq_chain != sequence:
        raise ValueError('chain sequence differs from the reference sequence; provide matching structures')
    residues = []
    for i, b in enumerate(sequence, 1):
        r = rm.get(i)
        residues.append({} if r is None else {k: v for k, v in r.atoms.items() if k in atom_order(b)})
    missing = [i for i, r in enumerate(residues, 1) if not r]
    return residues, missing


def import_ensemble(paths, out_dir, name='imported', chain_id=None, sequence=None, validate=True):
    """Package external structures (one conformer per file, first model) into an ensemble directory."""
    out = Path(out_dir)
    (out / 'conformations').mkdir(parents=True, exist_ok=True)
    (out / 'metadata').mkdir(parents=True, exist_ok=True)
    conformers = []
    c4 = []
    for k, p in enumerate(paths, 1):
        chains, meta = read_structure(p)
        if chain_id:
            ch = chains[chain_id]
        else:
            ch = max(chains.values(), key=lambda c: len(c.residues))
        seq = sequence or ch.sequence
        sequence = seq
        residues, missing = _chain_to_residues(ch, seq)
        cid = f'conf_{k:03d}'
        write_pdb(out / 'conformations' / f'{cid}.pdb', residues, seq,
                  remarks=[f'IMPORTED FROM {Path(p).name} CHAIN {ch.label_asym_id}'])
        write_mmcif(out / 'conformations' / f'{cid}.cif', residues, seq, entry_id=cid.upper(), title=f'{name} {cid}')
        val = validate_model(seq, residues) if validate and not missing else None
        conformers.append({
            'id': cid, 'rank': k, 'files': {'pdb': f'conformations/{cid}.pdb', 'mmcif': f'conformations/{cid}.cif'},
            'source': {'file': str(Path(p).name), 'sha256': hashlib.sha256(Path(p).read_bytes()).hexdigest(),
                       'chain': ch.label_asym_id, 'title': meta.get('title'), 'method': meta.get('method')},
            'missing_residues': missing,
            'validation': None if val is None else {k2: val[k2] for k2 in ('passed', 'severe_clashes', 'issue_counts',
                                                                          'sugar_puckers')},
        })
        c4.append({i: r["C4'"] for i, r in enumerate(residues) if "C4'" in r})
    # pairwise C4' RMSD matrix over residues common to all conformers
    common = sorted(set.intersection(*[set(d) for d in c4])) if c4 else []
    n = len(c4)
    mat = [[0.0] * n for _ in range(n)]
    if len(common) >= 3:
        arrs = [np.array([d[i] for i in common]) for d in c4]
        for a in range(n):
            for b in range(a + 1, n):
                mat[a][b] = mat[b][a] = round(kabsch_rmsd(arrs[a], arrs[b]), 3)
    manifest = {
        'schema': SCHEMA,
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'software': {'name': 'rnafold3d', 'engine': 'external-import'},
        'input': {'name': name, 'sequence': sequence, 'length': len(sequence or ''),
                  'sha256': hashlib.sha256((sequence or '').encode()).hexdigest(),
                  'secondary_structure_mode': 'not_applicable'},
        'chain': {'id': 'A', 'numbering': '1..N matches input sequence positions', 'atoms': 'heavy atoms only'},
        'conformers': conformers,
        'pairwise_c4_rmsd': {'ids': [c['id'] for c in conformers], 'matrix': mat, 'n_common_residues': len(common)},
        'notes': ['Imported ensemble: no prediction scores; validation is the rnafold3d independent checker.'],
    }
    (out / 'metadata' / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    return manifest


def load_package(path):
    """Read and minimally check an ensemble package (for viewers and tests)."""
    p = Path(path)
    man = json.loads((p / 'metadata' / 'manifest.json').read_text())
    if man.get('schema') != SCHEMA:
        raise ValueError(f'unsupported manifest schema {man.get("schema")!r}')
    for c in man['conformers']:
        for f in c['files'].values():
            if not (p / f).exists():
                raise FileNotFoundError(p / f)
    return man
