# RNAFold3D-Physics: Algorithm Design Specification

**Document ID:** RFP-ADS-001  
**Version:** 0.1 (implementation baseline)  
**Date:** 2026-10-21  
**Status:** Proposed; not a validated predictor  
**Target:** RNA molecules of 20–200 nucleotides; priority evaluation range 50–150 nt

## 1. Purpose and scope

Develop a reproducible, non-deep-learning RNA tertiary-conformation generator that accepts a sequence and secondary structure and exports a *diverse ensemble* of geometrically validated, all-atom RNA PDB/mmCIF files. The prediction engine is an upstream component of a future RNA ensemble visualization application. It must remain usable independently of that application.

**Primary requirements:** correct chemical connectivity; preservation of specified secondary-structure pairs; multiple distinct candidate conformations; deterministic reruns given identical seeds and versions; traceable scoring and validation. **Not guaranteed:** native-fold identification, calibrated equilibrium populations, general performance for every 200-nt RNA, or correct ligand/Mg²⁺-dependent folds without corresponding conditions.

**Non-goals for v1:** deep learning, fragment-template lookup from known RNA tertiary structures, molecular dynamics at explicit-solvent/all-atom resolution, interactive visualization, prediction of unknown secondary structure, and automatic inference of noncanonical tertiary contacts.

## 2. Input/output contract

### 2.1 Inputs

- RNA sequence in FASTA or plain text; uppercase A/C/G/U; T is rejected or explicitly normalized under an opt-in flag.
- Dot-bracket structure of exactly the same nucleotide length. Standard parentheses supported in v1; extended bracket alphabets and crossing pairs are a planned extension.
- Optional: temperature (K), monovalent salt concentration (M), Mg²⁺ concentration (M), desired number of models K (default 10), random seed, sampling budget, topology-search settings, and explicit experimentally known tertiary restraints.
- Secondary structure is an *input constraint*, not a predicted result. Contradictory or impossible inputs must yield actionable errors.

### 2.2 Outputs

```text
run_<id>/
  input/sequence.fasta
  input/secondary_structure.dbn
  conformations/conf_001.pdb
  conformations/conf_001.cif
  conformations/conf_002.pdb
  ...
  metadata/manifest.json
  metadata/scores.csv
  metadata/clusters.json
  metadata/validation.json
  config.json
  logs/run.log
```

All exported conformers must have consistent sequence, chain IDs, residue numbering, atom names, and model ordering. Export PDB only if format limits permit; mmCIF is the authoritative lossless format. A rejected conformer must never be silently exported as a valid prediction.

`manifest.json` must record input hashes, software/parameter versions, seed, sampling budget, ionic assumptions, conformer IDs, relative file paths, clustering method, ranking method, and validation status. Energy values are model-dependent scores, **not** probabilities or comparable absolute free energies across different parameter sets.

## 3. Prior-art algorithm analysis and design decisions

The following table separates established algorithmic ideas from our engineering choices. Details that vary by implementation/version require confirmation against the original articles and source releases before reproduction.

| Method | Structural representation | Junction handling | Loop closure/continuity | Stem orientation search | Escape from local minima | Design takeaway |
|---|---|---|---|---|---|---|
| **Ernwin** | Coarse structural elements; helices and loops assembled in a graph | Samples relative arrangements of coarse elements using structural priors | Element connectivity and closure are enforced/assessed during assembly; not an atomistic loop-closure solver | MCMC moves on element geometry | Resampling and Monte Carlo acceptance | Reuse graph abstraction, not structure-derived motif priors in strict template-free mode |
| **RNAJP** | Coarse-grained RNA with junction-sensitive interactions | Specialized helix/junction and loop interaction terms | Backbone connectivity and junction geometry are modeled during refinement; exact closure implementation must be checked in source | Conformational search guided by junction interactions | Junction rearrangement/relaxation and repeated sampling | Add explicit junction topology hypotheses and rearrangement moves |
| **IsRNA / IsRNA2** | Multi-site nucleotide coarse graining (IsRNA2 described with five sites per nucleotide) | Emergent from force field and secondary constraints | Continuous backbone geometry during simulation; does not by itself guarantee perfect all-atom closure | Thermal dynamics across replica temperatures | Replica-exchange molecular dynamics | Benchmark refined coarse-grained scoring and replica exchange; do not duplicate complex force field in v1 |
| **TiRNA** | Three-site nucleotide coarse graining (P, C4′, base site N1/N9) | Enhanced stem-orientation exploration, effective ionic effects | Connected backbone during sampling; final atomic geometry requires reconstruction/validation | Stem-angle sampling and conformational Monte Carlo | Simulated annealing / replica-exchange Monte Carlo | Use three-site representation and temperature-exchange sampling as practical baseline |

**Prior-art caveat:** this is a design-oriented synthesis, not a verified source-code audit. In particular, precise move sets, loop closure mechanisms, parameter provenance, and ion treatments must be checked against each software version. Avoid claiming originality for rigid-stem assembly, multi-resolution modeling, or ensemble generation alone.

## 4. High-level pipeline

```text
Sequence + dot-bracket
    -> validate and parse
    -> secondary-structure graph (stems, loops, bulges, junctions)
    -> analytical A-form stem construction
    -> discrete topology hypotheses (coaxial stacks, relative stem order)
    -> fast geometric feasibility pruning
    -> loop/junction backbone construction and closure
    -> three-site coarse-grained Monte Carlo / replica exchange
    -> validity filtering + ensemble clustering
    -> all-atom reconstruction + restrained refinement
    -> final validation + PDB/mmCIF + metadata
```

A strict gate is required between topology generation and physical refinement. Geometrically impossible topologies should be rejected early, while soft scores must not replace chain-continuity validation.

## 5. Internal data structures

```python
@dataclass(frozen=True)
class BasePair:
    i: int  # zero-based
    j: int
    pair_type: str  # AU, GC, GU, other

@dataclass
class Stem:
    id: str
    pairs: list[BasePair]
    pose: np.ndarray  # 4x4 rigid transform, local -> global
    local_sites: np.ndarray  # (n_sites, 3)

@dataclass
class Loop:
    id: str
    residue_indices: list[int]
    left_anchor: int | None
    right_anchor: int | None
    kind: str  # hairpin/internal/bulge/junction-linker

@dataclass
class Conformer:
    sites: np.ndarray  # (N, 3, 3): P, C4prime, base site
    energy_terms: dict[str, float]
    topology_id: str
    seed: int
    valid: bool
    validation: dict
```

Represent the secondary structure as a graph of stems (nodes) and loop/junction connections (edges or hyperedges). Maintain a residue-level adjacency map for every covalent phosphodiester link. A pseudoknot-aware graph is required for later extended dot-bracket support; do not silently reinterpret crossing pairs as nested pairs.

## 6. Geometric construction

### 6.1 A-form stems

Generate idealized A-form RNA duplexes from standard chemical geometry, including helical rise, twist, radius, strand polarity, and canonical Watson–Crick base-pair orientation. Keep these as *versioned physical/geometric parameters* with documented units and sources, rather than fitting to a test RNA. GU wobble pairs require distinct geometry. All atoms/sites should be expressed in a local stem frame, then transformed using SE(3) rigid poses.

### 6.2 Junction hypotheses

For each junction, enumerate a bounded set of coaxial-stacking pairings and noncoaxial arrangements. Sample relative stem rotations and translations while preserving anchor directions. Use a beam of candidate partial assemblies rather than exhaustive enumeration, and deduplicate near-identical rigid arrangements.

### 6.3 Early feasibility pruning

Reject a proposed arrangement when (a) covalent linker end-to-end separation exceeds a conservative contour-length upper bound; (b) stem volumes clash severely; (c) anchor orientations imply impossible backbone geometry under a conservative reachable-set approximation; or (d) required secondary pairs are geometrically irreconcilable. A distance-only upper bound is necessary but **not sufficient** for closure. Any uncertain test must be permissive to avoid discarding potentially feasible folds.

### 6.4 Loop and junction closure

Construct loop backbones in internal coordinates with constrained bond lengths/angles and flexible torsions. Optimize to satisfy both endpoint positions **and endpoint frames**. Start with cyclic coordinate descent (CCD) and a damped least-squares closure objective; consider kinematic closure (KIC) as a later alternative. Multiple random starts are mandatory. Validate every consecutive phosphate/sugar linkage and reject closures that exceed specified tolerances. Do not use a distance penalty alone as proof of successful closure.

## 7. Three-site coarse-grained model

Each nucleotide has three sites: phosphate P, sugar C4′, and base representative N9 (purines) or N1 (pyrimidines). Because one base site alone does not uniquely define the base plane, derive a local base frame from neighboring backbone geometry and/or add auxiliary virtual orientation vectors for pairing and stacking calculations. Do not pretend that three positions by themselves fully resolve arbitrary base orientations.

### 7.1 Energy decomposition

\[
E = E_{\mathrm{bonded}} + E_{\mathrm{excluded}} + E_{\mathrm{pair}} + E_{\mathrm{stack}} + E_{\mathrm{electrostatic}} + E_{\mathrm{junction}} + E_{\mathrm{restraints}}.
\]

- **Bonded:** harmonic or tabulated bonds/angles, periodic torsional terms, chirality constraints.
- **Excluded volume:** short-range repulsion; omit or scale bonded near-neighbor exclusions consistently.
- **Pairing:** secondary-structure-specific distance and orientation restraints; canonical AU/GC and GU wobble geometry differentiated.
- **Stacking:** anisotropic base stacking using separation, normal alignment, and lateral offset; include intrahelical and potential tertiary stacks without indiscriminate attraction.
- **Electrostatics:** screened phosphate repulsion (Debye–Hückel baseline), with explicit dielectric, temperature, and ionic-strength assumptions. Mg²⁺-specific binding is out of scope for v1.
- **Junction:** soft coaxial stacking and orientation preferences; avoid imposing a single preferred fold.
- **Restraints:** optional user-provided contact/distance restraints, labeled separately from physics-only terms.

All energies must use documented units (e.g., kcal/mol and Å) and a single temperature convention. Parameter values should be tracked in a machine-readable, versioned file. Avoid treating arbitrarily weighted pseudo-energies as thermodynamically exact.

## 8. Conformation sampling

### 8.1 Move set

1. Rigid stem rotation/translation (with linker-aware proposals).
2. Junction angle and coaxial-stack reassignment.
3. Loop torsion crankshaft/pivot moves followed by closure.
4. Local nucleotide/base-orientation perturbations.
5. Occasional large topology reassembly moves to cross barriers.

### 8.2 Monte Carlo acceptance

For symmetric proposals at fixed temperature, accept a proposed move with probability `min(1, exp(-(E_new-E_old)/(k_B*T)))`. For asymmetric proposal distributions, include the Hastings proposal-ratio correction. Reject geometry-invalid candidates before evaluating acceptance. Use annealing for initialization and fixed-temperature chains for comparable sampling.

### 8.3 Replica exchange

Run multiple replicas at increasing temperatures. Periodically attempt swaps of configurations between neighboring temperatures using the correct exchange acceptance criterion based on inverse temperatures and energies. Track exchange rates, acceptance rates, energy traces, and independent-seed reproducibility. Adaptive ladder tuning is optional, but temperature ladders must be frozen for production comparisons.

```python
for topology in feasible_topologies:
    x = build_and_close(topology)
    if not geometry_valid(x):
        continue
    replicas = initialize_replicas(x, temperature_ladder)
    for sweep in range(num_sweeps):
        for r in replicas:
            proposal = propose_move(r.structure)
            proposal = repair_loop_closure(proposal)
            if geometry_valid(proposal):
                metropolis_hastings_update(r, proposal)
        if sweep % exchange_interval == 0:
            attempt_neighbor_replica_swaps(replicas)
        if sweep % snapshot_interval == 0:
            save_valid_snapshots(replicas)
```

**Important:** if the search includes explicit diversity bonuses, those alter the sampled distribution. Keep novelty-driven exploration separate from unbiased temperature-based sampling and final physical ranking.

## 9. Ensemble selection and ranking

Apply a strict validity gate first. Cluster valid candidates using aligned backbone RMSD, RNA-specific eRMSD, or contact-map similarity; report which distance was used. Select medoids or best-scoring members from distinct clusters rather than simply the K lowest energies. Report cluster occupancy as a **sampling statistic**, not a physical state probability. Provide separate `physics_score`, `geometry_metrics`, `cluster_id`, and `rank` fields. Optional secondary statistical scorers must be labeled and disabled in strict no-structure-data mode.

## 10. All-atom reconstruction and refinement

Reconstruct sugar rings, bases, phosphate oxygens, and full RNA backbone from coarse-grained sites, preserving residue identity and chain order. Prefer an independently validated reconstruction library/tool if its license and redistribution terms permit; otherwise implement a constrained internal-coordinate builder with nucleotide stereochemistry checks. Perform restrained all-atom minimization with a documented force field and solvent approximation. Recheck the originally supplied secondary pairs and all covalent geometry after minimization. If reconstruction fails, retain the coarse-grained model internally and mark the all-atom export as failed; **never label a pseudoatom PDB as an all-atom model**.

## 11. Validation and acceptance criteria

**Per-conformer hard checks:** residue count and sequence identity; finite coordinates; no chain breaks; stereochemistry and bond geometry within specified tolerances; severe steric clash count below a predefined threshold; secondary-pair satisfaction above a predefined threshold; successful atom reconstruction; consistent residue and atom names.

**Per-run checks:** K distinct valid models when feasible; documented shortfall otherwise; all files readable by at least two independent structure parsers; stable results for a fixed seed/environment; no missing metadata; successful round-trip from mmCIF to parsed coordinates.

**Scientific benchmarks:** stratify by 20–50, 50–100, 100–150, and 150–200 nt; hairpins, 3-way and 4-way junctions, and pseudoknots (only once supported). Report Best-of-K RMSD, Top-1 RMSD, eRMSD, INF, valid-conformer rate, ensemble diversity, wall time, and success at fixed compute budgets. Evaluate ranking separately from candidate generation. Use held-out experimental structures **only for evaluation**, not as templates or parameter-fitting targets in strict mode.

**Release gates (initial engineering targets, to calibrate):** 100% correct parser behavior on fixture suite; 0 silently broken chains in exported models; 100% valid metadata schemas; deterministic fixed-seed smoke tests; at least three nonredundant valid conformers on designated small benchmark cases. Do not promise native-like RMSD before empirical benchmarking.

## 12. Software modules and public interfaces

```text
rnafold3d/
  io/             # FASTA, dot-bracket, config, PDB/mmCIF
  graph/          # stems, loops, junctions, pseudoknots
  geometry/       # A-form builder, SE(3), anchors, closure
  topology/       # enumeration, reachability pruning
  forcefield/     # bonded, pairing, stacking, electrostatics
  sampling/       # MC, annealing, replica exchange
  reconstruction/ # all-atom generation and refinement
  ensemble/       # validation, clustering, ranking
  schemas/        # JSON Schema, parameter versions
  tests/
```

CLI contract:

```bash
rnafold3d predict \
  --fasta input.fa \
  --dbn input.dbn \
  --models 10 \
  --seed 42 \
  --config physics-v1.yaml \
  --out run_001
rnafold3d validate --input run_001
rnafold3d cluster --input run_001 --metric ermsd
```

Suggested implementation: Python 3.11+, NumPy, SciPy, Numba for hot loops, optional C++ acceleration after profiling, BioPython/Gemmi for file interoperability, pytest and JSON Schema for testing. Keep the computational core independent of the future visualization UI.

## 13. Milestones

| Milestone | Deliverables | Exit criterion |
|---|---|---|
| M0: reference audit | Versioned paper/source notes for Ernwin, RNAJP, IsRNA2, TiRNA; licensing inventory | Every borrowed mechanism traceable to a source |
| M1: geometry | Parser, structure graph, A-form stems, topology generator, closure solver | Small test cases produce continuous coarse-grained chains |
| M2: sampling | Three-site force field, MC/annealing, optional replica exchange | Multiple valid distinct candidates with deterministic logs |
| M3: atomistic output | Reconstruction, restrained minimization, PDB/mmCIF, metadata | All-atom files pass independent validators |
| M4: benchmark | Held-out 20–200 nt suite, baselines, ablations | Accuracy/diversity/runtime results reproducible at equal budget |

Implement and validate M1 on short hairpins and 3-way junctions before extending to 200 nt. Pseudoknots and strong Mg²⁺-dependent systems require explicit feature flags until independently tested.

## 14. Risk register

| Risk | Mitigation |
|---|---|
| Wrong global topology despite low score | Diverse topology hypotheses; cluster-level reporting; Best-of-K evaluation |
| Loop closure failure | Multiple starts, frame-aware closure, early feasibility pruning, hard validation |
| Oversimplified ionic physics | Document salt assumptions; avoid claiming accurate Mg²⁺ binding |
| Poor base-plane orientation in three-site model | Auxiliary local frames or orientation sites; validate stacking geometry |
| Expensive search near 200 nt | Hierarchical search, adaptive local refinement, fixed compute budgets |
| Reconstruction produces bad stereochemistry | Validated atom builder, restrained minimization, independent validation |
| Data leakage in benchmark | Freeze parameters and split targets before final evaluation |
| License incompatibility | Audit external code and parameter licenses before integration |

## 15. Originality and implementation policy

The prediction engine prioritizes **reliability > ensemble diversity > computational efficiency > novelty**. Modest possible contributions include conservative loop-reachability pruning, reproducible topology-first sampling, explicit geometric acceptance gates, and interoperable ensemble metadata. None should be called novel until compared against prior work. The future visualization/analysis application is a separate research and product deliverable and is deliberately outside this specification.

## 16. Reference starting points (verify versions before coding)

- Ernwin source: https://github.com/ViennaRNA/ernwin
- IsRNA server: https://rna.physics.missouri.edu/IsRNA/index.html
- TiRNA source: https://github.com/Tan-group/TiRNA
- SimRNAweb v2.0: https://doi.org/10.1093/nar/gkae356
- RCSB PDB (benchmark experimental structures): https://www.rcsb.org/

**Reference audit task:** locate primary Ernwin, RNAJP, IsRNA2, and TiRNA methods publications; record exact algorithmic mechanisms, equations, parameters, versions, and software licenses in `docs/prior_art_audit.md` before copying implementation details. This specification intentionally does not present unverified method-specific details as exact reproductions.
