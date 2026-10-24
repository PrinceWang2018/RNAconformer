"""Build the frozen pilot benchmark dataset (WP1) from benchmark/manifest_v1.json.

For every target this writes benchmark/data/samples/<PDB>_<chain>/:
    sequence.fasta            full deposited polymer sequence (canonical parents of modified residues)
    reference.pdb             normalized coordinates: chain A, residue number = label_seq_id (1..N)
    residue_map.json          label_seq_id -> author numbering, residue names, modification flags
    oracle_2d.dbn             nested canonical pairs annotated from the reference (ORACLE: uses 3D)
    oracle_2d_full.dbn        all annotated canonical pairs incl. pseudoknots (for reference only)
    independent_2d.dbn        ViennaRNA MFE from sequence alone (independent of the 3D structure)
    metadata.json             provenance, QC, checksums, strata and 2D comparison
and benchmark/data/dataset_v1.json (index) and dataset_v1.csv (inventory).

Usage: PYTHONPATH=src python benchmark/build_dataset.py [--download]
"""
import argparse
import csv
import datetime
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'src'))

from rnafold3d.structure_io import read_mmcif, write_pdb          # noqa: E402
from rnafold3d.annotate import canonical_pairs, nested_subset, pairs_to_dotbracket  # noqa: E402
from rnafold3d.build import atom_order                             # noqa: E402

DEV_SET = {'2KOC', '480D', '1EHZ', '1Y26', '1KXK', '1GID'}


def stratum(n, strata):
    for name, (lo, hi) in strata.items():
        if lo <= n <= hi:
            return name
    return 'out_of_range'


def f1(pred, ref):
    pred, ref = set(pred), set(ref)
    if not pred and not ref:
        return 1.0
    tp = len(pred & ref)
    p = tp / len(pred) if pred else 0.0
    r = tp / len(ref) if ref else 0.0
    return 0.0 if tp == 0 else 2 * p * r / (p + r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default=str(ROOT / 'manifest_v1.json'))
    ap.add_argument('--raw', default=str(ROOT / 'data' / 'raw'))
    ap.add_argument('--out', default=str(ROOT / 'data'))
    ap.add_argument('--download', action='store_true')
    a = ap.parse_args()
    manifest = json.loads(Path(a.manifest).read_text())
    raw = Path(a.raw)
    out = Path(a.out)
    import RNA
    index = []
    rows = []
    for t in manifest['targets']:
        pid, chain = t['pdb_id'].upper(), t['chain']
        sid = f'{pid}_{chain}'
        path = raw / f'{pid}.cif'
        rec = {'sample_id': sid, 'pdb_id': pid, 'chain': chain, 'family': t['family'], 'source': t['source'],
               'split': 'dev' if pid in DEV_SET else 'test'}
        try:
            if not path.exists():
                if not a.download:
                    raise FileNotFoundError(f'{path} missing (run with --download)')
                url = f'https://files.rcsb.org/download/{pid}.cif'
                data = urllib.request.urlopen(url, timeout=60).read()
                path.write_bytes(data)
            chains, meta = read_mmcif(path)
            if chain not in chains:
                raise ValueError(f'chain {chain} not found (available: {sorted(chains)})')
            ch = chains[chain]
            seq = ch.sequence
            n = len(seq)
            bad = sorted(set(seq) - set('ACGU'))
            if bad:
                raise ValueError(f'non-ACGU symbols in canonical sequence: {bad}')
            rmap = ch.residue_map()
            observed = sorted(rmap)
            if observed and (observed[0] < 1 or observed[-1] > n):
                raise ValueError('label_seq_id outside 1..N')
            for r in ch.residues:
                if r.base != seq[r.index - 1]:
                    raise ValueError(f'residue {r.index} base {r.base} differs from sequence {seq[r.index - 1]}')
            c4_cov = sum("C4'" in rmap[i].atoms for i in observed) / n
            sdir = out / 'samples' / sid
            sdir.mkdir(parents=True, exist_ok=True)
            (sdir / 'sequence.fasta').write_text(f'>{sid}\n{seq}\n')
            # normalized reference: canonical atom names only, numbering = label_seq_id
            residues, seqs, nums = [], [], []
            lines = [f'REMARK 250 NORMALIZED REFERENCE {pid} CHAIN {chain} (label_asym_id); RESIDUE NUMBER = label_seq_id',
                     'REMARK 250 MODIFIED NUCLEOTIDES RENAMED TO THEIR CANONICAL PARENTS; NON-STANDARD ATOMS DROPPED']
            serial = 1
            for r in ch.residues:
                allowed = set(atom_order(r.base))
                for name in atom_order(r.base):
                    if name not in r.atoms or name not in allowed:
                        continue
                    x, y, z = r.atoms[name]
                    el = 'O' if name.startswith('OP') else name[0]
                    nm = name.ljust(3) if len(name) < 4 else name
                    lines.append(f'ATOM  {serial:5d}  {nm:<3s} {r.base:>3s} A{r.index:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}')
                    serial += 1
            lines += ['TER', 'END']
            (sdir / 'reference.pdb').write_text('\n'.join(lines) + '\n')
            rmap_json = {str(r.index): {'auth_asym_id': r.auth_asym_id, 'auth_seq_id': r.auth_seq_id,
                                        'ins_code': r.ins_code, 'comp_id': r.comp_id, 'parent': r.base,
                                        'modified': r.modified, 'n_atoms': len(r.atoms)} for r in ch.residues}
            (sdir / 'residue_map.json').write_text(json.dumps(rmap_json, indent=1))
            # oracle 2D from the reference (canonical pairs only)
            pairs = canonical_pairs(ch.residues)
            pp = [(i, j) for i, j, _ in pairs]
            nested, removed = nested_subset(pp)
            oracle = pairs_to_dotbracket(n, nested)
            oracle_full = pairs_to_dotbracket(n, pp)
            (sdir / 'oracle_2d.dbn').write_text(f'>{sid}|oracle_2d|nested canonical pairs annotated from reference 3D\n{seq}\n{oracle}\n')
            (sdir / 'oracle_2d_full.dbn').write_text(f'>{sid}|oracle_2d_full|all canonical pairs incl. pseudoknots\n{seq}\n{oracle_full}\n')
            md = RNA.md()
            md.noLP = 1
            mfe, e = RNA.fold_compound(seq, md).mfe()
            (sdir / 'independent_2d.dbn').write_text(f'>{sid}|independent_2d|ViennaRNA {RNA.__version__} MFE noLP 37C {e:.2f} kcal/mol\n{seq}\n{mfe}\n')
            from rnafold3d.annotate import dotbracket_to_pairs
            ind_pairs = dotbracket_to_pairs(mfe)
            meta_out = {
                **rec, 'title': meta['title'], 'method': meta['method'], 'resolution': meta['resolution'],
                'length': n, 'observed_residues': len(observed), 'missing_residues': [i for i in range(1, n + 1) if i not in rmap],
                'c4_coverage': round(c4_cov, 3), 'modified_residues': [(r.index, r.comp_id) for r in ch.residues if r.modified],
                'stratum': stratum(n, manifest['length_strata']),
                'source_file': f'raw/{pid}.cif', 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                'retrieved_from': f'https://files.rcsb.org/download/{pid}.cif', 'prepared': datetime.date.today().isoformat(),
                'oracle_2d': {'pairs': len(pp), 'nested_pairs': len(nested), 'pseudoknot_pairs_removed': len(removed),
                              'method': 'rnafold3d.annotate.canonical_pairs (cWW + GU wobble, H-bond geometry)'},
                'independent_2d': {'method': f'ViennaRNA {RNA.__version__} MFE, noLP, 37 C, Turner 2004',
                                   'energy_kcal_mol': round(e, 2), 'pairs': len(ind_pairs),
                                   'f1_vs_oracle_nested': round(f1(ind_pairs, nested), 3),
                                   'f1_vs_oracle_full': round(f1(ind_pairs, pp), 3)},
                'review_status': 'automatic QC passed; manual review pending',
            }
            (sdir / 'metadata.json').write_text(json.dumps(meta_out, indent=1))
            index.append(meta_out)
            rows.append({'sample_id': sid, 'split': rec['split'], 'length': n, 'stratum': meta_out['stratum'],
                         'observed': len(observed), 'c4_coverage': meta_out['c4_coverage'],
                         'modified': len(meta_out['modified_residues']), 'oracle_pairs': len(nested),
                         'pk_pairs': len(removed), 'vienna_f1': meta_out['independent_2d']['f1_vs_oracle_nested'],
                         'family': t['family'], 'status': 'ready', 'error': ''})
        except Exception as exc:
            rows.append({'sample_id': sid, 'split': rec['split'], 'status': 'excluded', 'error': str(exc)[:200]})
            print('EXCLUDED', sid, exc)
    (out / 'dataset_v1.json').write_text(json.dumps({'dataset_version': manifest['dataset_version'],
                                                     'created': datetime.date.today().isoformat(),
                                                     'dev_set': sorted(DEV_SET), 'samples': index}, indent=1))
    keys = ['sample_id', 'split', 'length', 'stratum', 'observed', 'c4_coverage', 'modified', 'oracle_pairs',
            'pk_pairs', 'vienna_f1', 'family', 'status', 'error']
    with (out / 'dataset_v1.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})
    print(f'{sum(r["status"] == "ready" for r in rows)} ready, {sum(r["status"] != "ready" for r in rows)} excluded')


if __name__ == '__main__':
    main()
