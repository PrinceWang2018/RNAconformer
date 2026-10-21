# RNAFold3D-Physics: Project Progress and Next Steps

**Status date:** 21 October 2026  
**Scope:** Physics-inspired RNA 3D conformation generation for RNAs up to 200 nucleotides, followed by a future multi-conformation visualization platform.  
**Document type:** Engineering status and validation roadmap

## 1. Project Objective

The long-term goal is a user-friendly RNA 3D visualization and analysis application that can load both internally generated ensembles and user-supplied PDB/mmCIF ensembles. The upstream RNAFold3D-Physics prediction engine is intended to produce several structurally distinct, geometrically valid RNA conformations from an RNA sequence and secondary structure. The engine avoids deep learning and prioritizes reliability, reproducibility, and interoperability over algorithmic novelty.

The visualization application has **not yet been implemented**. Current development concerns the upstream prediction and evaluation pipeline.

## 2. Architecture at a Glance

```text
RNA sequence + secondary structure (dot-bracket)
    -> secondary-structure graph (stems, loops, junctions)
    -> idealized A-form stem initialization
    -> rigid-body stem assembly and loop closure
    -> coarse-grained P/C4'/N guide coordinates
    -> geometry checks and physics-inspired scoring
    -> Monte Carlo / annealing / replica-exchange sampling
    -> clustering and ensemble selection
    -> experimental all-atom reconstruction and refinement
    -> independent stereochemical checks
    -> PDB candidates + metadata
    -> benchmark against experimental RNA structures
    -> [future] multi-conformation visualization software
```

## 3. Milestone Summary

| Milestone | Implemented prototype capabilities | Important remaining limitations |
|---|---|---|
| M0 — Prior-method audit | Reviewed published/official descriptions of Ernwin, RNAJP, IsRNA/IsRNA2, and TiRNA. | Full source-level, function-by-function audit and licensing review remain incomplete. |
| M1 — Geometry and topology | Sequence/dot-bracket parsing; stem/loop/junction identification; A-form-inspired three-site guides; stem rigid-body placement; virtual C4' loop closure; geometry filtering. | Pseudoknot assembly incomplete; virtual-chain closure is not chemical backbone closure. |
| M2.1 — Energy and geometry | Decomposed bonded, pairing, orientation, stacking, steric and screened-electrostatic scores; loop three-site guides; local geometry calibration. | Parameters are heuristic and uncalibrated; energy units are not physical free energies. |
| M2.2–M2.4 — Sampling and ensembles | Monte Carlo annealing, replica-exchange prototype, loop moves, trajectory snapshots, RMSD clustering and representative selection. | Limited loop degrees of freedom; no demonstrated equilibrium sampling or validated structural accuracy. |
| M3.1–M3.3 — All-atom reconstruction | Template fitting to P/C4'/N anchors; rigid-chain and selected torsion refinements; experimental sugar-pucker adjustment; PDB writing; independent basic stereochemical checks. | Complete RNA backbone chemistry, coupled sugar/base geometry, robust 200-nt refinement, and mmCIF export are not established. |
| M4–M4.3 — Benchmark framework | Reference PDB parsing, C4' RMSD/Best-of-K calculations, dataset manifest and annotation QC, residue mapping, reference normalization and execution scripts. | No confirmed real experimental benchmark runs have yet been completed in the reported environment. |

**Reported testing:** Earlier iterations reported increasing unit-test counts, reaching 87 tests for M3.3 and 91 for the M4 benchmark framework; M4.2 reported 7 annotation-tool tests. These are *reported prototype test results*, not independently reproduced in this status document, and do not establish predictive accuracy.

## 4. What Has Actually Been Demonstrated

- Software components and CLI workflows were developed for the stages above, with synthetic examples such as a 9-nt hairpin and a 24-nt multi-stem construct.
- The pipeline can represent and score coarse-grained RNA conformations and collect multiple candidate states.
- Prototype reconstruction can attempt to place full-atom nucleotide templates and reject candidates that fail basic checks.
- Benchmark code can compute aligned C4' RMSD and ensemble metrics when appropriate reference and predicted structures are available.

**Not yet demonstrated:** end-to-end generation of validated, chemically correct all-atom conformations for representative 50–200-nt experimental RNAs; comparative accuracy versus established RNA predictors; calibrated thermodynamic populations; robust performance across complex junctions and pseudoknots.

## 5. Current Scientific and Engineering Risks

1. **Reference data gap:** The M4.3 run reported 11 candidate PDB entries but zero successfully retrieved experimental coordinate files and zero real benchmark evaluations. Dataset population must be completed before making accuracy claims.
2. **Input leakage:** Keep *Oracle-2D* (secondary structure derived from target experimental coordinates) separate from *Independent-2D* (independently obtained secondary structure). Never mix their results.
3. **All-atom geometry:** Fitting whole nucleotide templates to three anchors does not guarantee valid O3'–P links, torsions, sugar puckers, stereochemistry, or absence of clashes.
4. **Sampling coverage:** Existing stem rigid-body and restricted loop moves may miss important conformational basins.
5. **Scoring validity:** Heuristic coarse-grained energies should not be reported in kcal/mol or interpreted as equilibrium populations.
6. **Sequence/structure correspondence:** Experimental author and label residue numbering, chain breaks, modified bases, and missing atoms require explicit mapping and auditing.
7. **Scale:** Existing atomistic refiners have had prototype size limits far below the 200-nt project target.

## 6. Immediate Priority: Populate a Real Benchmark Dataset

**Recommended first deliverable:** a manually reviewed pilot dataset of approximately 10–20 experimental RNA chains, followed by a broader 30–200-nt benchmark stratified by length and topology.

For every candidate:

1. Retrieve and preserve the original experimental mmCIF/PDB and provenance (accession, chain, experimental method, resolution where applicable, citation and retrieval date).
2. Select the biologically relevant RNA chain; identify ligands, ions, modified nucleotides, missing residues and alternate conformations.
3. Create a stable mapping from experimental `auth_*`/`label_*` identifiers to prediction indices; do not silently renumber away missing residues.
4. Require sequence consistency and sufficient C4' coordinate coverage; log exclusions and reasons.
5. Curate secondary structure separately for Oracle-2D and Independent-2D. Record method, source, reviewer and checksums.
6. Partition the benchmark by size (30–49, 50–99, 100–149, 150–200 nt), number of stems, junction complexity and pseudoknot status.
7. Freeze a versioned dataset manifest before prediction runs.

Potential sources: RNA-Puzzles standardized datasets, RCSB PDB, and RNA 3D Hub. Respect source licenses and provenance requirements.

## 7. End-to-End Benchmark Protocol

**Inputs:** frozen sequence and 2D annotation; no reference 3D coordinates supplied to the predictor.  
**Predictions:** multiple independent random seeds; fixed compute budget; save all valid candidates, not only the best.  
**Primary measures:**

- **Yield:** fraction of targets producing at least one all-atom PDB passing independent geometry checks.
- **Best-of-K C4' RMSD:** minimum aligned RMSD among the first K candidates, with coverage reported.
- **Top-1 RMSD:** RMSD of the engine's own top-ranked candidate.
- **Success@K:** proportion of targets with a candidate below a prespecified RMSD threshold.
- **Ensemble diversity:** pairwise RMSD and cluster counts at declared thresholds.
- **Runtime:** CPU/GPU, wall-clock time, memory, seed and hardware configuration.
- **Failure modes:** missing atoms, broken O3'–P bonds, stereochemical violations, clashes, incomplete loops, and unsuccessful reconstruction.

**Secondary measures to add:** RNA-aware eRMSD, interaction-network fidelity (INF), and external stereochemical validation. Compare against suitable versions of TiRNA, IsRNA, SimRNA, and/or Ernwin under as comparable inputs and compute budgets as possible.

Report Oracle-2D and Independent-2D separately. Predefine RMSD thresholds and treatment of incomplete predictions. Do not omit failed targets from aggregate success rates.

## 8. Engineering Work Packages

### WP1 — Experimental data ingestion (highest priority)
- Run M4.3 acquisition in an environment with reliable RCSB access.
- Review 10–20 real RNA chains and freeze `benchmark_manifest_v1.json`.
- Add checksums, provenance, chain mappings, exclusion reasons and independent secondary-structure annotations.

**Exit criterion:** at least 10 fully reviewed targets with valid reference coordinates, residue maps and documented 2D inputs.

### WP2 — End-to-end smoke tests
- Run the full prediction pipeline on a small hairpin and a multi-stem experimental RNA.
- Record whether any candidate survives all-atom reconstruction and independent validation.
- Preserve failure artifacts and diagnostic logs.

**Exit criterion:** reproducible, end-to-end runs with honest pass/fail reports; no silent fallback to invalid PDBs.

### WP3 — Repair the most consequential failures
- If atomistic yield is low, prioritize chain-aware covalent geometry and coupled sugar/backbone refinement before expanding sampling.
- If structures are valid but inaccurate, prioritize topology coverage, junction moves and energy calibration.
- If candidate generation is good but Top-1 is poor, improve ranking separately from search.

**Exit criterion:** measurable improvement on the frozen pilot benchmark without using reference structures to tune individual test cases.

### WP4 — Independent benchmarking and reproducibility
- Implement eRMSD and INF or integrate validated implementations.
- Add external RNA stereochemical checking.
- Compare against established tools and publish scripts, exact settings, failure counts and confidence intervals.

**Exit criterion:** a versioned benchmark report with real measurements and transparent limitations.

### WP5 — Output contract for future visualization
- Standardize PDB/mmCIF numbering, atom naming, chain IDs and residue correspondence across all conformers.
- Write a machine-readable `manifest.json` containing scores, validation flags, provenance, cluster IDs and seeds.
- Permit external PDB ensembles to use the same metadata schema, without requiring RNAFold3D predictions.

**Exit criterion:** an ensemble package specification that a future viewer can consume without knowledge of the prediction engine internals.

## 9. Proposed Deliverable Layout

```text
benchmark_v1/
  manifest.json
  samples/<sample_id>/
    reference.cif
    reference_normalized.pdb
    residue_map.json
    sequence.fasta
    oracle_2d.dbn
    independent_2d.dbn
    annotation_provenance.json
  runs/<run_id>/
    config.json
    conformations/conf_001.pdb
    conformations/conf_002.pdb
    scores.csv
    validation.json
    manifest.json
  reports/
    per_target.csv
    aggregate.json
    failures.csv
    benchmark_summary.md
```

## 10. Recommended Next Decision

**Do not begin the visualization application yet.** First, populate and review the real experimental dataset, execute a small end-to-end benchmark, and establish whether the current all-atom reconstruction can produce even a modest number of chemically acceptable conformers. This evidence should determine whether the next investment goes into atomistic refinement, sampling, or ranking.

The final product vision remains unchanged: a general-purpose RNA multi-conformation visualization and analysis platform, able to handle both predicted and user-imported structure ensembles. A stable, validated structural-data interface is the bridge between the two projects.

---

*Status interpretation:* “Implemented” in this report denotes functionality described in previous development iterations, not a claim of independent source-code audit or experimental validation. Reported tests and prototype outputs should be reproduced from the archived project packages before public release.
