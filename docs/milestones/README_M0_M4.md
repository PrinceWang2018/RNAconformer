# RNAFold3D-Physics — M0/M1 Starter

A minimal, tested, independent RNA secondary-structure parsing and geometry kernel. Intended for RNA up to 200 nt; no machine learning, structure fragments or upstream source code.

## Install and run

```bash
pip install -e .
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' --output example.json
python -m unittest discover -s tests -v
```

## Implemented
- Sequence and multi-bracket dot-bracket validation (including crossing base pairs)
- Pair map and consecutive base-pair stem detection
- Unpaired runs and preliminary structural-element adjacency
- Simplified A-form-like helical guide coordinates
- Proper rigid transformations and a conservative endpoint distance feasibility bound
- JSON export and six passing unit tests

## Not implemented
- Complete junction topology graph, precise loop classification, actual atomic A-form RNA geometry, loop closure, physical energies, Monte Carlo, full-atom reconstruction, or validated PDB outputs.

**Important:** Helix points are abstract guide coordinates, not physical atom positions. The endpoint reachability check is only a necessary upper-bound filter, not a sufficient chain-closure test. `T` is normalized to `U`; modified nucleotides are not supported. Bracket types support pseudoknot pair maps but not full pseudoknot topology modeling.

See `M0_AUDIT.md` for source-audit status and follow-up work.

## M1 extension (2026-10-21)

Added `topology.py` and `geometry.py`:
- Planar secondary-structure tree with external, hairpin, bulge, internal-loop, and multibranch-junction nodes, per-arm unpaired segments, and parent-child stem edges.
- Crossing base pairs are **detected but not decomposed**: `supported_planar=false`; `--strict` rejects them. This is intentional, not pseudoknot support.
- Analytic A-form-inspired P/C4'/N three-site **guide coordinates** and stem endpoint frames; proper rigid transformations.
- These sites are synthetic scaffolds, **not experimentally calibrated nucleotide geometries**, and cannot be labeled as validated atom coordinates or exported as all-atom PDBs.
- No loop closure, force field, energy minimization, or physical conformational ensemble yet.

Run: `PYTHONPATH=src python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' --output demo.json`

Test: `PYTHONPATH=src python -m unittest discover -s tests -v`

The earlier `core.py` `elements` field is retained for backward compatibility; use `topology.junctions` for the authoritative planar topology classification.

## M1.3 + M1.4 extension

`assembly.py` generates multiple seeded, topology-aware rigid-body stem arrangements. It initializes child stems near the parent stem's interior end and applies stochastic rigid translations/rotations, with a simple loop-reachability and inter-stem clash penalty. Scores are **heuristic geometry penalties, not physical free energies**.

`closure.py` implements FABRIK inverse kinematics for virtual C4' guide chains with fixed endpoints and equal virtual link lengths. Every attempted loop returns an explicit `converged` flag, endpoint residual and link-length residual. Unreachable endpoints are reported rather than silently connected. Unpaired terminal tails are not yet modeled; zero-unpaired inter-stem joins are checked. Loop orientation, sugar pucker, atomic torsions, steric clashes, base stacking, and chirality are **not** enforced.

Run a multi-candidate geometry demonstration:

```bash
PYTHONPATH=src python -m rnafold3d --sequence AAAAAAAAAAAAAAAAAAAAAAAA --structure '((..((...))..((...))..))' --assemble 8 --seed 42 --proposals 150 --output m1_assembly.json
PYTHONPATH=src python -m unittest discover -s tests -v
```

These JSON structures are **not PDB predictions**. A reported closure is a mathematical virtual-link closure only, not a chemically valid RNA backbone. Candidates with failed closure must be rejected or repaired in a later constrained optimizer. Pseudoknots are explicitly unsupported by this planar assembly stage.

## M1.5 — Geometry Validation & Candidate Filtering

The new `rnafold3d.validation` module performs **guide-coordinate** checks:
- Attempts C4' virtual-loop closure and records failure reasons and residuals.
- Reports unmodeled residues (including terminal tails); strict filtering rejects incomplete chains.
- Checks nonlocal guide-site clashes between different stems / loop residues, excluding close sequence neighbors and uncalibrated within-stem geometry.
- Checks virtual C4' link continuity outside individual stems.
- Ranks candidates by existing heuristic placement score, then greedily removes near-duplicates using Kabsch-aligned C4' RMSD.
- Returns machine-readable reports, including cluster membership and rejection reasons.

Example:

```bash
python -m rnafold3d --sequence AAAAAAAAAAAAAAAAAAAAAAAA \
  --structure '((..((...))..((...))..))' --assemble 12 --proposals 100 \
  --filter --max-selected 8 --rmsd-threshold 3.0 --seed 42 \
  --output m1_5_results.json
```

`--allow-incomplete` is an **exploratory option only**: terminal tails are still missing from coordinate output. Strict filtering can legitimately return zero candidates; no synthetic success is reported. Thresholds are provisional guide-site diagnostics, not experimentally validated RNA atomic geometry. M1.5 does not output PDB or certify molecular validity. M2 should add physically calibrated potentials and proper backbone geometry.


## M2.1 — Experimental Physics Energy Engine

`rnafold3d.energy` implements an independently testable **guide-coordinate surrogate score** with six components: adjacent C4 guide bonds, secondary-structure base-anchor pairing restraints, pair-facing orientation, direction-sensitive intra-stem stacking, nonlocal C4 steric repulsion, and screened phosphate-guide repulsion. See `EnergyConfig` for tunable parameters and `EnergyBreakdown` for decomposed scores.

**Scientific limitation:** The existing M1 P/C4/N sites are synthetic geometric guides; base-anchor directions are not physical base-plane normals. Energy values are **arbitrary model units, not kcal/mol or free energies**. The pairing term is a secondary-structure restraint, and the stacking term is a simple orientation proxy. Electrostatics omit explicit ions and solvent. This is a tested software foundation, not a validated physical force field. In particular, score comparisons across RNA sequences, lengths, or conformations with different modeled residue counts are not scientifically justified.

Example (stem-only diagnostic scores, not complete RNA conformer energies):

```bash
python -m rnafold3d --sequence AAAAAAAAAAAAAAAAAAAAAAAA \
  --structure '((..((...))..((...))..))' --assemble 4 \
  --score-energy --output m2_1_example.json
```

The optional `m2_1_energy` field is emitted for each assembly and scores **stem sites only**, because the M1 loop solver currently emits C4-only coordinates. The next step is to reconstruct complete P/C4/N loop sites and calibrate/test the geometric representation before physically meaningful conformational optimization.

## M2.1 bridge — complete internal-loop three-site guides

`rnafold3d.loop_sites.build_full_chain_guides(model, assembly)` adds P/C4/N **synthetic guide sites** for internal unpaired runs with successful FABRIK C4-chain closure. It preserves existing stem coordinates exactly, tracks missing residues and failed closures, and reports maximum cross-boundary C4 link-length error. Terminal tails are intentionally not extrapolated. Failed/incomplete candidates receive `m2_1_energy: null` rather than a misleading partial-chain ranking.

CLI assembly JSON now includes `full_chain_guides` and, with `--score-energy`, a complete-chain experimental score only when `full_chain_guides.complete` is true. All energy values remain arbitrary surrogate units. No atomistic backbone continuity, sugar pucker, base plane orientation, chirality, or steric validity is guaranteed.

Example:

```bash
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' \
  --assemble 2 --seed 42 --score-energy --output results.json
```

## M2.1 Geometry Calibration (guide-coordinate consistency)

Run with `--calibrate-guides` to regularize loop-only P/C4'/N local frames and
include `m2_1_geometry_calibration` in each assembly JSON record. The method
preserves every C4' guide coordinate and all Stem sites; it adjusts only the
P and N guide offsets of interior Loop residues with both adjacent residues
present. It checks local frame degeneracy, C4' neighbor distances and
Stem–Loop interface distances. See `src/rnafold3d/calibration.py`.

**Scientific limitation:** This is NOT experimental calibration and does not
implement RNA backbone torsion-angle distributions, ribose puckering,
phosphodiester bond stereochemistry, or real all-atom structure validation.
The constants are configurable heuristic guide targets, not force-field values.
This phase can flag bad geometry but cannot repair bad C4' paths or guarantee
that a physically valid atomistic RNA chain can be reconstructed.

Example:

```bash
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' \
  --assemble 2 --seed 42 --calibrate-guides --score-energy \
  --output calibrated.json
```

## M2.2 — Experimental conformational sampling

`src/rnafold3d/sampling.py` adds seeded Metropolis annealing and parallel-tempering-style replica exchange. The move is a whole-stem rigid transformation; each proposal regenerates and closes loops and is rejected if the complete three-site chain cannot be reconstructed. Energies remain **arbitrary guide-model units**, and temperature is a **dimensionless search control**, not Kelvin. The output is **not** an equilibrium thermodynamic ensemble or an atomistic PDB.

```bash
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' --assemble 2 --sample anneal --sample-steps 100 --output anneal.json
python -m rnafold3d --sequence AAAAAAAAAAAAAAAAAAAAAAAA --structure '((..((...))..((...))..))' --assemble 2 --sample remc --sample-steps 100 --output remc.json
```

Output `m2_2_sampling` includes best guide coordinates, component energies, acceptance and invalid-proposal counts, and traces. Candidate errors are reported instead of silently accepted. The original M1/M2.1 JSON output remains available.

### Important current limitations

* A single isolated stem has only global rigid-body moves and therefore cannot meaningfully change its internal energy; multi-stem RNAs are the intended sampling test.
* Loop coordinates are rebuilt by the existing FABRIK solver, not sampled with chemically accurate RNA backbone torsions.
* The sampler currently rejects incomplete or failed closures but does not impose full atomistic stereochemistry or an experimentally calibrated RNA free-energy potential.
* Pseudoknot topology and user-specified tertiary contacts remain unsupported by this search pipeline.
* Replica exchange swaps states between search temperatures, but because the approximate score and closure generator are not a calibrated statistical-mechanical model, it must not be used to infer physical populations.

## M2.3: energy ranking and conformer ensemble selection

`ensemble.py` aggregates M2.2 **best and final** sampled states across initial assemblies, rejects incomplete/nonfinite or geometrically invalid C4' guide chains, ranks by experimental model energy, and greedily selects structurally distinct representatives by Kabsch-aligned RMSD. The output JSON key `m2_3_ensemble` includes all inspected candidates, rejection reasons, scores, cluster identifiers, and full P/C4'/N guide coordinates for selected representatives.

Example:

```bash
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' \
  --assemble 8 --sample anneal --sample-steps 100 \
  --select-ensemble --max-selected 5 --rmsd-threshold 3.0 \
  --seed 42 --output ensemble.json
```

**Scientific limitations:** The energy values are arbitrary model scores. This is *not* a calibrated free-energy landscape, equilibrium ensemble, or full-atom RNA PDB predictor. It can return zero selected conformers under strict geometry gates. M2.2 currently retains only best and final snapshots per run; a later milestone should implement trajectory snapshot retention to increase conformational coverage. Selected coordinates are synthetic three-site guides and must not be presented as experimentally valid PDB atom positions.


## M2.4 — Ensemble Sampling Improvement

`ensemble_sampling.py` introduces reversible axial rotation of internal-loop three-site guides, preserving C4 virtual link lengths and fixed stem anchors; captures periodic snapshots from annealing or replica-exchange replicas; aggregates multiple initial assemblies and selects RMSD-diverse candidates. The loop move is **not** a chemical torsion move. No thermodynamic ensemble or all-atom RNA geometry is claimed.

```bash
python -m rnafold3d --sequence GGGAAACCC --structure '(((...)))' --assemble 3 --pool-sampling anneal --sample-steps 50 --snapshot-interval 5 --loop-move-probability 0.7 --seed 42 --output m24.json
```

For replica exchange, set `--pool-sampling remc`. Snapshots and selection reports appear under `m2_4_pool`. The output includes coordinates for selected representatives only, while the audit lists all candidates.

## M3 experimental atomistic export

See [M3_README.md](M3_README.md). This adapter requires external complete nucleotide templates and fails closed when basic linkage or clash validation fails. It does not yet provide stereochemical refinement or mmCIF export.

## M3.3: Independent stereochemical validation

The M3 atomistic exporter now records `stereochemistry` for each reconstructed conformer in `manifest.json` and rejects export by default if either basic or new stereochemical checks fail. See `M3_3_STATUS.md` for coverage and limitations. These are provisional diagnostic checks, not calibrated experimental validation.


## M4 benchmark harness
See `M4_BENCHMARK.md`. Run `python -m rnafold3d.benchmark --reference ref.pdb --predictions conf_*.pdb --output report.json`. This does not imply experimental validation.
