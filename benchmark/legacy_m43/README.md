# RNAFold3D M4.3 — Experimental Dataset Population and Benchmark Execution

## Scope
This release combines the M4.2 dataset collection/annotation scripts with a **reference-preparation and benchmark-execution driver**. It does not claim that any reference has been downloaded, annotated, manually reviewed, or benchmarked unless the execution report records it.

## Install and execute

```bash
pip install numpy biopython
PYTHONPATH=src python run_m43.py --dataset dataset --manifest candidate_manifest.json --online
```

Offline (uses only `dataset/raw/<PDB>.cif` already downloaded):

```bash
PYTHONPATH=src python run_m43.py --dataset dataset --manifest candidate_manifest.json
```

Predictions may be evaluated by adding `--predictions predictions`, where files reside under `predictions/<SAMPLE_ID>/*.pdb`. The reference is standardized to chain `A` and `label_seq_id`; prediction PDB files **must** use this same mapping. The program does not claim to run the RNAFold3D prediction engine or invent secondary structures. Produce input `.dbn` annotations and prediction PDBs separately; for independent-2D, provenance must be independent of the target 3D structure.

## Generated files

- `dataset/inventory.csv`: candidates and failures
- `dataset/samples/<SAMPLE_ID>/reference/reference.pdb`: normalized experimental coordinates (if available)
- `dataset/samples/<SAMPLE_ID>/reference/residue_mapping.json`: experimental author-to-label mapping
- `dataset/samples/<SAMPLE_ID>/reference/reference_qc.json`: basic reference QC and source checksum
- `dataset/samples/<SAMPLE_ID>/benchmark.json`: per-conformer C4' RMSD (only if predictions exist)
- `dataset/m43_execution_report.json`: overall execution and missing data status

## Scientific limitations

- The M4.1 extractor excludes modified nucleotides and does not validate all polymer sequence gaps. A sample is not accepted into a publication benchmark without manual curation.
- The standardization writes a reference PDB for scoring; original mmCIF must be retained. Legacy PDB residue numbering/atom field limitations apply.
- Current metrics cover aligned C4' RMSD, not eRMSD or INF.
- Running the scoring stage with an experimental reference is not the same as having run a complete sequence-to-PDB prediction.
- The provided PDB candidates are **unverified** until downloaded and inspected.
