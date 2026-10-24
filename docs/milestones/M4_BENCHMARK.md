# M4 — Benchmark and Scientific Validation Protocol (v0.1)

## Status
This deliverable is a **benchmark harness**, not an experimentally validated predictor. The new executable computes residue-matched, rigidly aligned C4′ RMSD and Best-of-K over supplied PDB files. It does not fabricate real RNA results.

## Reproducible protocol
1. Assemble a nonredundant reference set of experimentally determined single-chain RNAs, targeting 30–50, 51–100, 101–150 and 151–200 nt. Record PDB ID, chain, experimental method, resolution when available, ligands, metal ions, modifications, and secondary structure provenance.
2. Exclude RNA structures or homologous motifs used to tune any reconstruction templates or parameters from held-out testing. Never derive the predictor's input secondary structure from a hidden test reference without labeling it an oracle-secondary-structure experiment.
3. Freeze the predictor commit, configuration, ion conditions, random seeds, CPU/GPU model, and wall-clock budget. Run multiple independent seeds per RNA.
4. Export only validated PDB conformers as primary predictions; separately report invalid/rejected outputs and diagnostic debug structures.
5. Run `python -m rnafold3d.benchmark --reference experimental.pdb --predictions prediction_*.pdb --output benchmark.json` with matching chain IDs and residue numbering. For mismatched numbering, preprocess explicitly and document the mapping.
6. Report generation rate, fraction passing independent atomistic validation, C4′ RMSD Best-of-K, top-ranked RMSD, pairwise ensemble diversity, runtime and memory. Stratify by length, junction count, pseudoknots, and modified bases.
7. Implement eRMSD and interaction-network fidelity (INF) with established definitions/tools before making comparative claims. C4′ RMSD alone is insufficient for RNA structural accuracy.
8. Compare against available baseline tools under matched secondary-structure input, compute budget and ensemble size. Report failures and missing outputs, not only successes.

## Output contract
The JSON report includes a per-conformer C4′ RMSD or explicit failure reason, coverage, best-of-K, median RMSD and success-at-K threshold. A successful RMSD comparison is **not** evidence of chemically valid geometry.

## Scientific caveats
- PDB parser uses first model and first acceptable alternate location, standard RNA A/C/G/U, and exact (chain, residue number, insertion code) correspondence.
- Chain/residue mapping must be standardized for fair comparison. Residues without C4′ cannot be evaluated.
- Structural ensembles are not Boltzmann populations; uncalibrated energy rankings cannot be interpreted as probabilities.
- The M3 reconstruction pipeline may reject most arbitrary coarse-grained candidates. Zero exported PDBs is a legitimate and informative outcome.
- No genuine experimental RNA benchmark has been executed in this package; only synthetic regression fixtures are tested.

## M4 acceptance gates
- At least 30 curated experimental RNA targets spanning complexity and length, with provenance and frozen hold-out policy.
- Independent validation of bond/angle/dihedral, chirality, sugar pucker, clashes and backbone breaks; record all failure categories.
- Complete end-to-end per-target log from input sequence/dot-bracket to final conformer files, including zero-output cases.
- Baseline comparison and confidence intervals for performance statistics.
- Machine-readable benchmark summary, versioned dataset manifest and reproducible command logs.
