# RNAfold3D v2 engine — validation and modification report

**Date:** 2026-11-20 · **Engine:** `rnafold3d-v2-foldtree` 2.5.0 (current; repeated fit onto the template, §4.00000; modified nucleotides and Z-RNA since 2.4.0; motif fragments, §4.0000; structure-aware template search, §4.000; templates and junction geometry, §4.00) ·
**Geometry parameters:** `geometry_v1.json` 1.0.0 · **Structure library:** `pdb_rna_v2.npz` (3,964 PDB RNA chains)

**Summary.** The archived pipeline could not turn a sequence into a structure (0/26 benchmark
targets). The reworked engine produces validated all-atom ensembles from sequence alone for
26/26 targets (14–192 nt). Small RNAs are accurate (hairpins and TAR 1.8–4 Å); larger multi-branched
RNAs reach random-arrangement accuracy in the top-ranked model and modestly better in best-of-5.
Tertiary physics, pseudoknots and ranking — not chemistry — now limit accuracy. v2.1.0 adds related
experimental structures (templates) and junction geometry from a built-in PDB library: mean TM-score of model 1
on the 26 chains 0.214 → 0.271, and on the blind CASP16 targets 0.251 → 0.490 (AlphaFold 3: 0.480). v2.2.0 adds a
template search that also aligns by secondary structure and keeps models from different templates: 0.271 → 0.284
on the 26 chains; CASP16 model 1 unchanged (0.489), best of 5 0.511 → 0.524 (within the noise of 15 targets). v2.3.0 grafts
hairpin-loop fragments from the library: 0.284 → 0.291 on the 26 chains (0.216 → 0.227 without a template), CASP16
model 1 0.489 → 0.495 (within the noise); pseudoknot hypotheses and domain assembly are options, off (§4.0000). v2.5.0
repeats the fit onto the template with the residues at fault released instead of discarding it: on 48 chains
with a template model 1 0.582 → 0.605, CASP16 model 1 0.495 → 0.501 (§4.00000; v2.4.0 added modified nucleotides
and left unmodified runs unchanged).

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
thread fix; see `CHANGELOG.md`); both are reported in §4.0–4.1. Reports are copied to `docs/benchmark_v1/`.
The engines v2.1.0 (templates and junction geometry), v2.2.0 and v2.3.0 were measured on the same 26 chains with a
development benchmark that limits every target to structures released before it (§4.00, §4.000).

### 4.00000 v2.5.0 (current engine): the fit onto the template is repeated; where the remaining accuracy is

**Diagnosis** (development data; CASP16 kept as the blind check).

* *A template model is as good as its template and its fit.* The TM-score of a model equals the TM-score of
  its template under the alignment, as long as the final Cartesian fit onto the template atoms is accepted.
  In 0.7.5 the fit was discarded whenever the fitted model failed the geometry checks — 73 of 150 template
  models on the 48-chain set below — and those models kept the looser shape of the torsion-space search (CASP16
  R1271: restraints satisfied 0.39, 8.6–9.0 Å from the template atoms, TM 0.39–0.41, while the one model with an
  accepted fit had 0.65).
* *Upper bound of the library* (`benchmark/template_oracle.py`, 916 queries of `template_selection.py`; US-align,
  sequence-independent, each query against every chain released before it). 655 queries with a hit: the best
  chain in the library 0.752, the first-ranked template 0.639, the best template found 0.694; the first-ranked
  one is more than 0.10 below the best chain for 225 queries. 261 queries without a hit: best chain 0.471;
  ≥ 0.45 for 124, ≥ 0.70 for 21. CASP16 (`casp16/oracle_templates.py`): best chain 0.573 on average against 0.392
  for the first template found (0 where none); a chain ≥ 0.45 exists for 10 of 15 targets.
  Of the difference between the structural and the alignment-based TM-score of the first template (0.664
  against 0.594, 500 queries) about 0.047 is not real (pairs with identical sequence and full coverage lose as
  much); the real alignment headroom is about 0.025, in hits covering less than 90%.

**Adopted: `template_fit_attempts` = 4.** After a fit that fails validation, the residues with the problems and
their neighbours (one more on each side per attempt) are released from the template and the fit is repeated,
until the model passes or more than 35% of the chain would be released; models whose fit stays rejected rank
after the models with an accepted fit. Same seeds, the step on and off:

| | one fit (0.7.5) | 2.5.0 | difference (95% interval, bootstrap over targets) |
|---|---|---|---|
| 48 chains with a template (`template_models.py`, `fit1`, `fit4`): rejected fits | 73 of 150 | 8 of 138 | |
| – mean TM-score, model 1 | 0.582 | **0.605** | +0.023 (+0.007 to +0.043); 8 better by > 0.02, 2 worse |
| – mean TM-score, best of 5 | 0.613 | 0.622 | +0.009 (−0.001 to +0.023) |
| CASP16, 15 targets (`runs_v076_fit1`, `runs_v076`): model 1 | 0.495 | **0.501** | +0.006 (−0.000 to +0.017) |
| – best of 5 | 0.527 | 0.527 | +0.000 |
| – right global fold (TM ≥ 0.45), model 1 | 8 | 9 | |
| – median RMSD of model 1 | 7.8 Å | 6.0 Å | |
| – run time, 15 targets, 16 cores | 68 min | 68 min | |
| 6 chains with a template of the 26-chain set (`dev_templates.py`, `fit1_076`, `fit4_076`): model 1 | 0.508 | 0.497 | −0.011 (1ANR −0.03, 5KPY −0.02) |

The 48 chains are library chains of 40–140 nt released after their templates, CASP16 entries left out; each is
predicted with templates and junction geometry limited to earlier structures and its own entry excluded.
`runs_v076_fit1` reproduces `runs_v075` to the last digit. On CASP16 the whole gain is R1271 (0.391 → 0.474);
with 15 targets it is not established. Runs without a template are unchanged. A side effect on the confidence
model is described in `docs/V07_VALIDATION.md` §F.

**Measured and not adopted** (each on development data):

| Idea | Result |
|---|---|
| Rank templates by the share of the 2D energy they explain | first template 0.586 → 0.55–0.57 |
| Free partnerless or mismatched residues, or extend helices, in the template's 2D hypothesis | base-pair F1 0.7515 → 0.754 |
| Structure-only search for queries without a hit (half of 70 such queries have a library chain ≥ 0.45) | accepted candidates average 0.2: no specificity |
| Helix-aware gap penalties in the alignment | better for tRNAs (one template 0.556 → 0.649), neutral overall (0.5964 → 0.5958; first hits 0.6106 → 0.6076–0.6120) |
| Template ranking with a gap term, resolution, length; cross-validated linear ranker on 13 features | 0.6105 → +0.002 (best single rule), +0.0046 (ranker); the oracle is 0.653 |
| Structure alignment replaces the sequence alignment on a tie as well | 228 queries: first template unchanged (0.406 / 0.442); 44 hits better, 5 worse, first templates −0.035 in sum |
| Choose model 1 by the confidence model when another model leads by ≥ 0.1 | 48 chains, on top of 2.5.0: 0.605 → 0.605 |
| Five templates instead of three (`max_templates=5`, run `fit4_t5`) | 48 chains, on top of 2.5.0: model 1 +0.0001 (−0.0002 to +0.0005), best of 5 +0.0003 (−0.003 to +0.003) |

Run times of the three 48-chain runs (3,516 s, 5,055 s and 2,960 s for `fit1`, `fit4`, `fit4_t5`) were taken
while the computer was in other use and do not compare; the two CASP16 runs, taken one after the other, do.

**What limits accuracy now.** With a template: its choice (0.639 against 0.694 for the best one found, 0.752 for
the best in the library) — the rules tried above do not capture it. Without a hit: 124 of 261 queries have a
structural relative in the library that sequence and 2D alignment do not find. Without any relative: sampling
of the tertiary arrangement (§4.0000).

### 4.0000 v2.3.0: motif fragments; pseudoknots, domains and a base–base potential measured

What changed (details in `CHANGELOG.md`, package 0.7.5): hairpin loops of 3–12 nt take torsions and C4′
distances from library loops with a matching sequence, and grafted junction strands keep their torsions (on).
Three more items were built and measured with the feature on and off on the same 26 chains (16 trajectories,
seed 42; runs `ctl075`, `frag075`, `pk075`, `all075` of `benchmark/dev_templates.py`) and are off by default:

| Dev benchmark, mean TM-score | model 1, 26 chains | model 1, 20 without template | best of 5 | run time, 26 chains |
|---|---|---|---|---|
| 2.2.0 code, features off (`ctl075`) | 0.284 | 0.216 | 0.303 | 2,108 s |
| + fragments (`frag075`, the 2.3.0 default) | **0.291** | **0.227** | **0.307** | 2,092 s |
| + pseudoknot hypothesis (`pk075`) | 0.287 | 0.219 | 0.303 | 2,273 s |
| + both (`all075`) | 0.291 | 0.225 | 0.306 | 2,277 s |

* **Pseudoknot hypothesis** (`--pseudoknots auto`): 2D base-pair F1 on 1,509 library chains 0.624 → 0.631
  (`benchmark/pseudoknot_calibration.py`, penalty 4 kcal/mol best of 3–6); correct crossing pairs on the dev chains
  where checkable (1RNK 6/7, 4ENC 6/6, 3D0U 6/9), but the crossing helix rarely forms in the 3D models — the
  stated limit "pseudoknots" stands, now at the 3D-building step rather than the 2D step.
* **Domain assembly** (`--domains auto`, `benchmark/long_rna.py`, six library chains of 255–393 nt, templates off):
  model 1 0.091 against 0.104 for the whole-chain search, 2.5 against 3.3 valid trajectories of 16, 1.8× the time.
* **Base–base potential** (`basepot.py`): offline re-ranking of the saved candidates (`benchmark/rerank_basepot.py`)
  gains nothing at a usable weight. On the chains without a template, random choice among the candidates 0.220,
  engine ranking 0.221, the best candidate 0.254: sampling, not ranking, is the limit.
* **Speed**: numpy-level changes only; a 2-trajectory profile of 1Y26 82 → 76 s (the library motif index), the
  trajectory inner loop unchanged.
* **Template search** (`benchmark/template_selection.py`): unchanged code, unchanged numbers (0.428 / 0.461).

CASP16, 15 targets (`runs_v075`, results in `benchmark/casp16/results/`), paired against 2.2.0 (`compare.py
paired`, bootstrap): model 1 0.489 → **0.495**, difference +0.006 (95% interval −0.016 to +0.026); best of 5 0.524 →
0.527, +0.003 (−0.007 to +0.016). The six targets without a template: model 1 0.263 → 0.281 (R1209 +0.09, R1296
+0.04, R1205 −0.02, R1293 −0.02). Template targets move because the fragments change the trajectories of the
uncovered parts (R1203 +0.05, R1271 −0.11 in model 1, best of 5 unchanged). Not established with 15 targets.

### 4.000 v2.2.0: structure-aware template search

What changed (details in `CHANGELOG.md`, package 0.7.4): library chains with a weak sequence alignment (E-value
up to 10) are aligned by sequence and secondary structure together (pairing probabilities of the query against
the base pairs of the template) and kept when ≥ 60% of the query is covered and ≥ 90% of the covered template
pairs can form; the same alignment extends sequence hits, except across an uncovered internal part of ≥ 20 nt;
the kept models take the best model of each template in turn; SHAPE reactivities and a user alignment
(RNAalifold consensus) can enter the 2D hypotheses. The template ranking rule is unchanged: no alternative was
better on the template benchmark.

| | 2.1.0 | 2.2.0 |
|---|---|---|
| Template search, 916 library chains, own class excluded (`benchmark/template_selection.py`): mean TM-score of the first-ranked template | 0.416 | **0.428** |
| … of the best template found | 0.444 | **0.461** |
| … chains with a template | 642 | 655 |
| Development benchmark, 26 chains: mean TM-score, model 1 | 0.271 | **0.284** |
| … best of 5 | 0.292 | **0.302** |
| … targets with a template | 5 | 7 |
| CASP16, 15 targets: mean TM-score, model 1 | 0.490 | 0.489 |
| … best of 5 | 0.511 | 0.524 |

CASP16, paired over targets (`compare.py paired`, bootstrap): model 1 difference −0.000 (95% interval −0.001 to
+0.000), the first model comes from the same template on every target; best of 5 +0.012 (−0.001 to +0.036),
from R1271 (0.50 → 0.66) and R1203 (0.49 → 0.52), not established with 15 targets. The search parameters were
set on the template benchmark; an intermediate version that mis-aligned a non-homologous insert in a shared
scaffold (R1203) led to the uncovered-part rule, so CASP16 is not blind for that rule. Still open: without a
template the 3D arrangement of helices is as in 2.1.0.

### 4.00 v2.1.0: templates and junction geometry

What changed (details in `CHANGELOG.md`, package 0.7.3): a search of the built-in PDB library (Smith–Waterman,
E-values, ranking by 3D agreement between hits); template 2D, C4′ distance restraints, torsions and atom
positions for template hypotheses; loop closure restrained to the template loop shape; Cartesian refinement
onto template atoms; junction geometry from non-homologous library structures as soft restraints for the other
hypotheses (left out of the ranking score); Cartesian regularisation of models that fail validation only
through local strain.

**Development benchmark** (`benchmark/dev_templates.py`; the 26 pilot chains, 16 trajectories, seed 42;
templates and junction geometry limited to structures released before each target, so no target is its own
template). Mean US-align TM-score (normalised by the reference, as in CASP16):

| | physics only (as 2.0) | + junction geometry | + templates (default) |
|---|---|---|---|
| Mean TM-score, model 1 | 0.214 | 0.230 | **0.271** |
| Mean TM-score, best of 5 | 0.244 | 0.254 | **0.292** |
| Targets with a template | – | – | 5 of 26 |

Junction geometry: better on 12 targets, worse on 5 (by ≤ 0.05). Templates: 1EHZ tRNA-Phe 0.42 → 0.89 (486D),
3MXH 0.25 → 0.75 (3IWN), 5KPY 0.20 → 0.39 (1Y27); 1ANR (29-nt NMR structure) 0.21 → 0.16 and 4QLM 0.23 → 0.20.

**CASP16** (`benchmark/casp16/`, 15 blind targets, templates released before 2024-05-01): mean TM-score of
model 1 0.490 with templates (AlphaFold 3 0.480, trRosettaRNA 0.493, Vfold 0.534, GuangzhouRNA 0.526;
2.0.2: 0.251). The template rules were set with the CASP16 results visible; the development benchmark
above is the independent check. Limitations: a template is only as good as its similarity to the target (tRNA
R1271: closest earlier tRNA conformation 0.50); coarse cryo-EM templates can make loop closure fail (R1256: 6
of 40 trajectories valid); without a template the 3D arrangement of helices is still the weak point.

### 4.0 v2.0.1 (`benchmark/runs/v2_0_1_pilot/`)

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

`docs/benchmark_v1/v2_0_1_pilot/benchmark_summary.md` (engine v2.0.1) and
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
