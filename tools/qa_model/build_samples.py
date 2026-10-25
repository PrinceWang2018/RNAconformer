"""Training samples for the confidence model from the chains of fetch_training_set.py.

For every downloaded chain that folds on itself (>= 4 nested canonical pairs) and is >= 90 % observed this writes
<out>/samples/<PDB>_<chain>/ with sequence.fasta, reference.pdb (normalised like benchmark/build_dataset.py),
oracle_2d.dbn and metadata.json, and <out>/dataset_v1.json in the format benchmark/run_benchmark.py reads, so the
decoys are made with:

    PYTHONPATH=src python benchmark/run_benchmark.py --data benchmark/data_qa --run-id qa_train_v1 \
        --modes independent,oracle --trajectories 8 --models 8 --workers 6
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d.structure_io import read_mmcif                                         # noqa: E402
from rnafold3d.annotate import canonical_pairs, nested_subset, pairs_to_dotbracket    # noqa: E402
from rnafold3d.build import atom_order                                                # noqa: E402


def reference_pdb(pid, chain, residues):
    lines = [f'REMARK 250 NORMALIZED REFERENCE {pid} CHAIN {chain}; RESIDUE NUMBER = label_seq_id',
             'REMARK 250 MODIFIED NUCLEOTIDES RENAMED TO THEIR CANONICAL PARENTS; NON-STANDARD ATOMS DROPPED']
    serial = 1
    for r in residues:
        for name in atom_order(r.base):
            if name not in r.atoms:
                continue
            x, y, z = r.atoms[name]
            el = 'O' if name.startswith('OP') else name[0]
            nm = name.ljust(3) if len(name) < 4 else name
            lines.append(f'ATOM  {serial:5d}  {nm:<3s} {r.base:>3s} A{r.index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}')
            serial += 1
    return '\n'.join(lines + ['TER', 'END']) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=str(ROOT / 'benchmark' / 'data_qa'))
    ap.add_argument('--min-pairs', type=int, default=4)
    a = ap.parse_args()
    data = Path(a.data)
    chains = json.loads((data / 'training_chains.json').read_text())['chains']
    index, excluded = [], {}
    for c in chains:
        sid = f"{c['pdb_id']}_{c['chain']}"
        f = data / 'raw' / f'{sid}.cif'
        try:
            if not f.exists():
                raise FileNotFoundError('not downloaded')
            parsed, _ = read_mmcif(f)
            if len(parsed) != 1:
                raise ValueError(f'{len(parsed)} chains in the file')
            ch = next(iter(parsed.values()))
            seq = ch.sequence
            if set(seq) - set('ACGU'):
                raise ValueError('non-ACGU sequence')
            rmap = ch.residue_map()
            if len(rmap) < 0.9 * len(seq):
                raise ValueError(f'{len(rmap)} of {len(seq)} residues observed')
            if sum("C4'" in r.atoms for r in ch.residues) < 0.9 * len(seq):
                raise ValueError('C4\' coverage below 90 %')
            for r in ch.residues:
                if r.base != seq[r.index - 1]:
                    raise ValueError(f'residue {r.index} differs from the sequence')
            pairs = [(i, j) for i, j, _ in canonical_pairs(ch.residues)]
            nested, removed = nested_subset(pairs)
            if len(nested) < a.min_pairs:
                raise ValueError(f'{len(nested)} nested canonical pairs (does not fold on itself)')
            sdir = data / 'samples' / sid
            sdir.mkdir(parents=True, exist_ok=True)
            (sdir / 'sequence.fasta').write_text(f'>{sid}\n{seq}\n')
            (sdir / 'reference.pdb').write_text(reference_pdb(c['pdb_id'], c['chain'], ch.residues))
            (sdir / 'oracle_2d.dbn').write_text(f'>{sid}|oracle_2d|nested canonical pairs annotated from reference 3D\n{seq}\n{pairs_to_dotbracket(len(seq), nested)}\n')
            meta = {'sample_id': sid, 'pdb_id': c['pdb_id'], 'chain': c['chain'], 'split': 'train', 'family': c['class'],
                    'stratum': 'train', 'length': len(seq), 'title': c['title'], 'method': c['method'],
                    'resolution': c['resolution'], 'nested_pairs': len(nested), 'pseudoknot_pairs': len(removed)}
            (sdir / 'metadata.json').write_text(json.dumps(meta, indent=1))
            index.append(meta)
        except Exception as exc:
            excluded[sid] = str(exc)[:160]
    (data / 'dataset_v1.json').write_text(json.dumps({'dataset_version': 'qa_train_v1', 'created': datetime.date.today().isoformat(),
                                                      'samples': index, 'excluded': excluded}, indent=1))
    lengths = sorted(s['length'] for s in index)
    print(f'{len(index)} samples ready, {len(excluded)} excluded; lengths {lengths[:1]}..{lengths[-1:]} median {lengths[len(lengths) // 2] if lengths else None}')
    reasons = {}
    for v in excluded.values():
        k = v.split(' (')[0] if 'nested' in v else v.split(':')[0][:40]
        reasons[k] = reasons.get(k, 0) + 1
    print(json.dumps(dict(sorted(reasons.items(), key=lambda x: -x[1])[:10]), indent=1))


if __name__ == '__main__':
    main()
