# RNAfold3D v2 engine — validation and modification report

**Date:** 2026-10-22 · **Engine:** `rnafold3d-v2-foldtree` 2.0.1 (current) · **Geometry parameters:** `geometry_v1.json` 1.0.0

**Summary.** The archived pipeline could not turn a sequence into a structure (0/26 benchmark
targets). The reworked engine produces validated all-atom ensembles from sequence alone for
26/26 targets (14–192 nt). Small RNAs are accurate (hairpins and TAR 1.8–4 Å); larger multi-branched
RNAs reach random-arrangement accuracy in the top-ranked model and modestly better in best-of-5.
Tertiary physics, pseudoknots and ranking — not chemistry — now limit accuracy.

## 1. Why the prediction engine was reworked

Reproducing the archived M0–M4.3 packages showed:

| Check | Finding |
|---|---|
| Test suite | 91/91 legacy tests pass (reproduced). |
| Sequence-only input | Not supported: the legacy CLI requires a dot-bracket; no 2D prediction. |
| All-atom output | Never produced from a sequence: the M3 exporter needs nucleotide template files that were never supplied. |
| Guide geometry | The analytic P/C4'/N guide sites are not chemically consistent with real nucleotides (e.g. P–C4' 3.0 Å vs ~3.9 Å; base-site pairing target 6 Å vs ~9 Å), so rigid template fitting cannot give bonded chains. |
| Benchmark data | M4.3 retrieved 0 references; no real evaluation had been run. |
| Legacy baseline (this work) | Given ideal templates and oracle 2D, the legacy M1–M3 path yields **0 PDBs on 26/26 targets** (§5.1). |

Per WP3 of the roadmap ("if atomistic yield is low, prioritise chain-aware covalent geometry"),
the engine was rebuilt so that covalent chemistry is correct by construction, and sequence-only
input was added. Legacy modules are kept unchanged and still pass their tests.

## 2. v2 method

1. **2D hypotheses (sequence only).** ViennaRNA 2.7 (Turner 2004, 37 °C, no lonely pairs):
   MFE plus up to two structurally distinct suboptimals from deterministic Wuchty enumeration
   (base-pair distance ≥ 15 % of MFE pairs). Trajectories: 60 % MFE, rest split. A user
   dot-bracket may replace this step (pseudoknotted and non-canonical pairs are dropped with a warning).
2. **Calibrated local geometry** (`tools/calibrate_geometry.py`): 178 high-resolution
   (≤ 1.6 Å) RNA-only crystal structures, 4108 complete nucleotides, **disjoint from all benchmark
   targets and families** (title-keyword and ID exclusion, chains ≤ 50 nt). Derived: sugar units
   for C3'-endo/C2'-endo, planar base templates, phosphate group, bond/angle statistics,
   Watson–Crick/wobble pair frames, A-form torsions (refined for duplex self-consistency:
   14-bp pair-fit 0.69 → 0.20 Å), a 1158-residue single-nucleotide loop conformer library and an
   empirical loop-reach table. Derived values agree with published standards (O3'–P 1.602 Å,
   C3'–O3'–P 120.2°, O3'–P–O5' 104.0°; A-form α/β/γ/δ/ε/ζ/χ = −67/173/56/79/−151/−72/−163°).
3. **All-atom builder** (`build.py`): NeRF placement along P–O5'–C5'–C4'–C3'–O3' with rigid
   sugar/base/phosphate units. Built stems reproduce tRNA-Phe helices to 0.5–1.0 Å all-atom RMSD.
   Stems with a compressed pair (any H-bond < 2.4 Å, e.g. next to G·U wobbles) or a contact
   < 2.5 Å get a regularised least-squares refinement of strand-2 torsions (v2.0.1; ~12 % of stems).
4. **Fold tree** (`fold.py`): stems are rigid bodies; a minimum spanning tree of short linkers
   carries the global arrangement; remaining linkers and hairpins are closure loops.
   Torsion rotations keep every bond and angle exact (verified to 1e-13°).
5. **Search** (`search.py`): staged simulated annealing — (i) stems + tree linkers, short
   (≤ 4 nt) closure loops explicit; (ii) remaining loops grown from the conformer library with
   clash/reach checks and closed by CCD; (iii) joint annealing with every loop explicit.
   Every move is followed by CCD re-closure of the loops it perturbs (closure-preserving moves),
   sterics ramp soft → hard, a chain-break restraint holds loops closed. (iv) exact CCD closure,
   (v) clash relief of loops/tails/χ.
6. **Score** (dimensionless, kT-like; *not* a free energy): coarse-site and all-atom
   (bond-graph-aware, Verlet list) sterics, chain-break, compaction prior
   Rg ≤ 5.5·N^⅓ Å for branched topologies only (Hyeon et al. 2006), coaxial helix-end
   stacking, screened phosphate repulsion (Debye length 10 Å), torsion prior from the conformer
   library. Ranking adds the ViennaRNA ΔΔG of the model's 2D hypothesis (in kT).
7. **Independent validation** (`checks.py`, fail-closed): complete heavy atoms, bonds within
   max(6σ, 0.08 Å), angles within max(6σ, 8°), O3'–P 1.50–1.72 Å for every link, D-ribose
   chirality at C1'–C4', no non-bonded heavy-atom pair < 2.2 Å. Real structures pass
   (2KOC, 480D, 1GID); failures on 1EHZ/1Y26 are genuine (modified bases; a 2.08 Å P–O5' bond).
8. **Ensemble selection:** greedy score-ordered C4'-RMSD diversity (threshold 0.15·Rg₀, 2–6 Å).
9. **Output contract (WP5):** `metadata/manifest.json` (`schemas/ensemble_manifest.schema.json`),
   PDB + mmCIF per conformer, scores/clusters/validation files; `rnafold3d import` packages
   external PDB/mmCIF ensembles in the same format (tested on the 2KOC NMR ensemble).

## 3. Benchmark design (WP1)

`benchmark/manifest_v1.json`: 26 experimental single chains, 14–192 nt (8 RNA-Puzzles targets),
prepared by `benchmark/build_dataset.py` with checksums, residue maps (label/auth), modified-residue
mapping, oracle 2D (canonical pairs annotated from the reference, nested subset) and independent 2D
(ViennaRNA MFE from sequence only). **Dev split** (used for all development decisions):
2KOC, 480D, 1EHZ, 1Y26, 1KXK, 1GID. **Test split:** the other 20, evaluated once with frozen settings.
Success thresholds fixed in advance: best-of-K C4' RMSD ≤ 5 Å / ≤ 10 Å, eRMSD ≤ 1.0.

## 4. Results

Two complete runs were made with identical settings: **v2.0.0** (the frozen engine, run once on all
targets) and **v2.0.1** (after failure analysis: stem refinement for compressed pairs, BLAS
thread fix; see `CHANGELOG.md`). v2.0.1 is the current engine; both are reported. Reports are copied
to `docs/benchmark_v1/`.

### 4.0 v2.0.1 (current engine, `benchmark/runs/v2_0_1_pilot/`)

| | Legacy M1–M3 | v2.0.1 sequence only | v2.0.1 oracle 2D |
|---|---|---|---|
| Targets with ≥ 1 valid all-atom model | 0 / 26 | **26 / 26** | 25 / 26 |
| Valid trajectories (all targets) | — | 159 / 208 | 149 / 208 |
| Median top-1 / best-of-5 C4' RMSD | — | 21.2 / 17.5 Å | 18.0 / 14.3 Å |
| Best-of-5 ≤ 5 Å / ≤ 10 Å (of 26) | — | 3 / 4 | 3 / 5 |
| Median top-1 INF-WC | — | 0.77 | 0.97 |
| Median wall time per target (8 traj., 8 logical CPUs) | — | 125 s | 133 s |

Held-out **test split** (20 targets): sequence only — models for 20/20, median top-1 21.7 Å, best-of-5
19.3 Å, null 19.4 Å, top-1 better than null 8/20, best-of-5 better than null 14/20; oracle 2D — 19/20,
top-1 19.4 Å, best-of-5 14.4 Å, null 19.1 Å, top-1 better than null 11/19, best-of-5 16/19.

**Run-to-run variability.** Results are deterministic for a given engine and seed (unchanged
targets reproduced to the last digit). Between v2.0.0 and v2.0.1, 25/52 target–mode results changed
because some stems changed; per-target best-of-5 values moved by up to ±6 Å in both directions.
With 8 trajectories, per-target numbers are therefore dominated by sampling noise; only aggregate
trends and large differences are meaningful.

**Ranking (v2.0.1).** Dev split median top-1 under alternative rules: engine score 11.9 Å vs 12.7–13.2 Å
(sequence only) and 12.5 vs 12.4–13.1 Å (oracle) — no rule is clearly better, ranking unchanged.
Sampling ceiling on the test split (best valid candidate): 18.0 / 14.4 Å vs top-1 21.7 / 19.4 Å.

### 4.1 v2.0.0 (frozen engine, `benchmark/runs/v2_0_pilot/`)

Budget: 8 trajectories, ≤ 5 exported models per target and mode, seed 42, 8 workers. Full tables:
`benchmark/runs/v2_0_pilot/reports/benchmark_summary.md`; per model: `per_model.csv`.

| | Legacy M1–M3 | v2 sequence only | v2 oracle 2D |
|---|---|---|---|
| Targets with ≥ 1 valid all-atom model | 0 / 26 | **25 / 26** | 23 / 26 |
| Valid trajectories (all targets) | — | 152 / 208 | 140 / 208 |
| Median top-1 C4' RMSD (26 targets*) | — | 19.3 Å | 16.1 Å |
| Median best-of-5 C4' RMSD (26 targets*) | — | 16.6 Å | 14.1 Å |
| Best-of-5 ≤ 5 Å / ≤ 10 Å | — | 3 / 5 of 26 | 3 / 5 of 26 |
| Median top-1 INF-WC (base pairs) | — | 0.77 | 0.97 |
| Median wall time per target (8 traj., 8 CPUs) | — | 102 s | 107 s |

\*medians over targets with a valid model; failures are listed separately and count as failures in the
success rates.

**Held-out test split (20 targets, evaluated once):**

| | sequence only | oracle 2D |
|---|---|---|
| Targets with ≥ 1 valid model | 19 / 20 | 17 / 20 |
| Median top-1 / best-of-5 C4' RMSD | 21.0 / 17.6 Å | 18.6 / 16.1 Å |
| Random-arrangement null† (same targets) | 19.1 Å | 18.4 Å |
| Top-1 better than null | 5 / 19 | 9 / 17 |
| Best-of-5 better than null | 13 / 19 | 14 / 17 |
| Best-of-5 ≤ 10 Å | 3 / 20 | 3 / 20 |

†median C4' RMSD of 30 random fold-tree arrangements with the correct (oracle) stems and no search
(`benchmark/null_baseline.py`), median over the same targets.

**Interpretation (both runs).**
* The rework fixes the engineering failure: the engine now turns a sequence into valid, closed,
  stereochemically checked all-atom ensembles for 14–192-nt RNAs (25/26 targets in v2.0.0, 26/26 in
  v2.0.1), where the legacy pipeline produced nothing.
* Hairpins and small motifs are predicted well from sequence alone (UUCG 1.8 Å, HIV-1 TAR 3.4 Å,
  sarcin–ricin 4.4 Å best-of-5; 3.5 Å with oracle 2D).
* For multi-branched RNAs the search improves on random arrangements in best-of-5 (e.g. tRNA 10.8 vs
  20.6 Å, P4-P6 23.1 vs 36.8 Å, M-box 17.6 vs 25.7 Å), but the **top-ranked model is often no better
  than random** (5/19 better than null on the test split). Global folds of riboswitches/ribozymes,
  which depend on tertiary contacts, pseudoknots and ligands, are not recovered.
* Ranking: on the dev split no alternative rule (score without torsion strain, medoid/consensus)
  beat the engine score (median top-1 12.2 vs 12.4–12.8 Å sequence-only; 11.6 vs 12.2–13.5 Å oracle),
  so ranking was left unchanged (`reports/rerank_analysis.json`). On the test split the engine score
  was also best except for a medoid restricted to the MFE hypothesis in sequence-only mode
  (19.9 vs 21.0 Å) — reported, not adopted. The sampling ceiling (best valid candidate) is
  17.4 / 15.4 Å vs 21.0 / 18.6 Å for top-1: better ranking alone is worth ~3.5 Å in median.
* Sequence-only vs oracle 2D: INF-WC 0.77 vs 0.97; median best-of-5 RMSD 16.6 vs 14.1 Å.

## 5. Evidence tables

### 5.1 Legacy pipeline baseline (WP2 "before")

`benchmark/legacy_baseline.py` → `benchmark/runs/legacy_baseline/legacy_baseline.json`. The legacy
M1–M2.4 guide pipeline (8 assemblies, pooled annealing) followed by M3 template export, given oracle 2D
and ideal nucleotide templates generated by the v2 builder (its best case):

| Outcome | Targets |
|---|---|
| No complete guide chain (terminal tails/loops cannot be modelled) | 24 / 26 |
| Complete guide chains, atomistic candidates built | 2 / 26 (2KOC: 3, 1ANR: 7 candidates) |
| Candidates passing the legacy checks (exported PDBs) | 0 (max O3'–P "bond" 10.4–10.8 Å) |
| Candidates passing the v2 validator | 0 |

### 5.2 / 5.3 Per-target tables (dev and test splits, both modes)

`docs/benchmark_v1/v2_0_1_pilot/benchmark_summary.md` (current engine) and
`docs/benchmark_v1/v2_0_pilot/benchmark_summary.md` (frozen v2.0.0): per target valid trajectories,
top-1 / best-of-5 / best-candidate C4' RMSD, null RMSD, eRMSD, INF-WC, ViennaRNA F1 and wall time.
Raw data: `per_target.csv`, `per_model.csv`, `failures.csv`, `rerank_analysis.json`;
dataset index `docs/benchmark_v1/dataset_v1.csv`; baselines `legacy_baseline.json`, `null_baseline.json`.

## 6. Limitations

* **No tertiary-contact physics.** Scores contain sterics, chain closure, helix-end stacking, a
  compaction prior and screened electrostatics; there are no base–base tertiary interactions
  (A-minor motifs, tetraloop–receptor, kissing loops, ribose zippers) and no Mg²⁺. Global folds
  of riboswitches and large ribozymes are therefore usually not recovered (§4).
* **Pseudoknots.** The fold tree is nested-only. Crossing pairs of a user 2D can be applied as
  base-pair distance restraints (`PredictConfig(pseudoknot_restraints=True)`), but this is
  EXPERIMENTAL: on the dev target 1Y26 restrained runs did not recover the kissing pairs and
  produced invalid models; it is off by default. Sequence-only mode never predicts pseudoknots.
* **Non-canonical pairs** in a user 2D are dropped with a warning (only A-U, G-C, G-U helices are built).
* **Ranking** is weakly informative: top-1 is usually not the best of the sampled candidates.
  Scores are dimensionless model units, not free energies; cluster sizes are not populations.
* **Sampling budget and variance.** Benchmarks used 8 trajectories per target and mode (pilot
  budget); wall time ≈ 15 s (14 nt) to ≈ 6 min (192 nt) per 8 trajectories on 8 logical CPUs.
  Per-target results vary by several Å between runs; large RNAs still yield only 0–5 valid models
  out of 8 trajectories (e.g. 4P9R oracle 0/8 in v2.0.1).
* **Validation** is heavy-atom only (no hydrogens, no MolProbity); chirality, bonds, angles,
  linkage and severe clashes are hard gates, rotamer quality is not assessed.
* **2D input dominates.** Sequence-only accuracy is bounded by ViennaRNA MFE accuracy
  (median base-pair F1 ≈ 0.8 on this set; 0 for the H-type pseudoknot 1RNK).
* The pilot set (26 chains) is small; numbers have wide uncertainty and should not be compared
  to published methods without matched inputs and budgets.

## 7. Recommended next steps

1. **Tertiary physics (largest expected gain):** base–base stacking and hydrogen-bond terms for
   unpaired nucleotides (coarse base frames already available), junction-family priors for
   3-/4-way junctions (Lescoute & Westhof), and a pseudoknot-capable fold tree (crossing stems as
   jumps with two closure loops) instead of distance restraints.
2. **Ranking:** a statistical-potential re-ranker trained only on non-benchmark structures, or an
   ensemble-consensus (centroid) selection; report Spearman(score, RMSD) on the dev split.
3. **Throughput:** numba/C++ for CCD and neighbour lists (spec §12) to afford 50–100 trajectories
   for 150–200 nt; enable restrained all-atom minimisation (e.g. OpenMM, AMBER OL3) as an optional
   polish step behind the same validator.
4. **Benchmark:** extend to ≥ 30 reviewed targets with RNA-Puzzles-style all-atom RMSD and
   external validators (MolProbity), and compare against an established tool at matched 2D input.
5. **Viewer (now unblocked):** the ensemble package (`schemas/ensemble_manifest.schema.json`,
   `rnafold3d import` for external PDB/mmCIF ensembles) is stable enough to start the
   visualization application against it.
