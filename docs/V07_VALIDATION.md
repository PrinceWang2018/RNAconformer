# RCV 0.7 — validation of the new analyses

This report checks the analyses added in RCV 0.7 on the RCV benchmark (`benchmark/manifest_v1.json`): the
topology-aware ensemble map, the transition analysis (basins, routes, kinetics, bottleneck), the φ-value scan, the
mechanism morph, probing-data reweighting and the optional confidence model. Every number below is produced by
two scripts and stored as JSON next to this file:

    PYTHONPATH=src python tools/validate07.py --out docs/benchmark_v1/v07_validation.json        # A-E, ~1 min
    PYTHONPATH=src python tools/qa_model/evaluate.py --eval benchmark/runs/v2_0_1_pilot/qa_eval.npz \
        --weights src/rnafold3d/data/qa_v1.npz --out docs/benchmark_v1/qa_v1_eval.json          # F

**Ensembles.** For each of the 26 benchmark targets the engine's exported models of both modes of the pilot run
(`benchmark/runs/v2_0_1_pilot`, *independent* and *oracle*, 5 + 5 models) form one ensemble of up to 10
conformations of one sequence (14–192 nt). These are the ensembles a user would open in RCV. Energies are ViennaRNA
Turner 2004 energies at 37 °C and standard salt. Run times were measured on a laptop-class CPU (WSL2).

## A. Topology-aware ensemble map vs classical MDS

Map distances should follow the C4′ RMSD (stress, Spearman correlation, trustworthiness and continuity of
k-neighbourhoods, Venna & Kaski 2001). Map neighbours should also share a topology (abstract shape) and a 2D state
(purity: the share of each point's k nearest map neighbours with the same label). Medians over the targets:

| | map | classical MDS |
|---|---|---|
| stress (lower is better) | **0.197** | 0.233 |
| Spearman correlation of distances | 0.870 | **0.933** |
| trustworthiness | 0.889 | **0.931** |
| continuity | 0.909 | **0.946** |
| 2D-state purity of neighbours (25 ensembles with ≥ 2 states) | **0.70** | 0.45 |
| *9 ensembles with several topologies:* topology purity | **0.90** | 0.72 |
| *9 ensembles with several topologies:* trustworthiness | 0.84 | **0.90** |

The map has lower stress than MDS in 22 of 25 ensembles. It keeps topologies and 2D states together much better.
The price is a small loss of neighbourhood preservation, because a conformation can no longer sit next to a similar
one of another topology. That trade is the design goal: secondary structure is a constraint, not a feature that
competes with RMSD. The separation between 2D states inside a disk (0.1 × the median RMSD) was chosen on these
data: 0 lowers the 2D-state purity to 0.63, and 0.25 lowers trustworthiness and continuity (0.877 / 0.877) for the
same purity.

## B. Transitions between the two most different 2D states

For each ensemble with at least two 2D states, A and B are the two states with the largest base-pair distance
(25 pairs).

**Basins.** In 28 of the 50 end points, the model's 2D structure is not a local minimum: it lies a median 1.15
kcal/mol above the minimum it relaxes into. Measuring barriers from such a point overstates or understates them,
so RCV 0.7 measures kinetics between basins. In 1 pair (480D) A and B relax into the same basin, and the analysis
reports that no barrier separates them. The most extreme case is 4GXY: the direct path between the model structures
climbs 74.3 kcal/mol, because model B is 72.5 kcal/mol above its basin (a strained model, not a real barrier).
Between the basins the direct path climbs 16.4 kcal/mol, the best route 12.0, and the network gives 11.8.

**Routes and kinetics** (24 pairs in different basins; medians):

| | barrier (kcal/mol) |
|---|---|
| direct path between the basins (findpath) | 6.1 |
| best single route through intermediates (minimax) | 5.9 |
| effective barrier of the network, RT·ln(k₀·MFPT) | 5.5 |

The network is faster than the direct path between the basins in 19 of 24 pairs: parallel routes and
intermediates speed up the switch. In 3 pairs (5TPY, 3MXH, 1U9S) it is slower, because the RNA spends 93–100 % of
the switching time in wells deeper than the basin of A. These wells are an A-side minimum (3MXH), or a state about
as stable as B on the way to it (5TPY, 1U9S). The analysis reports them as traps. The occupancy behind this rule
(the expected time in each node before B is reached) was checked against a stochastic simulation of the same
network (0.9986 vs 0.9986, `tests/test_rcv07.py`).

All 25 analyses gave identical results when repeated (fixed sampling seed). They take a median of 0.05 s, at most
0.5 s for the largest network (3MXH, 30 nodes), including the 3D module analysis.

## C. φ-value scan: are the kinetic hot spots at the cooperative unit?

The cooperative unit (base pairs that open together on the climb to the highest saddle) is predicted from the
energy profile alone. The φ-value scan independently evaluates all 3 × length point mutations. Over the 21 pairs
with a unit:

* the unit nucleotides are more sensitive than the rest in **21 of 21** pairs: median of the per-pair mean |ΔΔG‡|
  is 3.3 vs 0.33 kcal/mol;
* the unit's median sensitivity percentile is **0.885** (0.5 would be no relation);
* a median of 50 % of the most sensitive positions fall into the unit (as many positions as the unit has), against
  7.9 % expected by chance.

The scan takes a median of 0.12 s and at most 1.1 s (4GXY, 172 nt, 516 mutants).

## D. Mechanism morph: does the chain stay intact?

Chain distortion is the largest relative change of a consecutive C4′–C4′ distance along a 25-frame morph from A
to B (0 = every virtual bond keeps its length). Medians over 25 pairs:

| morph | chain distortion |
|---|---|
| linear interpolation of coordinates | 0.81 |
| rigid modules (all turning together) | 0.21 — lower than linear in 25 of 25 |
| staged mechanism morph (modules in the order of their base-pair events) | 0.20 — lower than linear in 24 of 25 |

Linear interpolation shrinks virtual bonds by up to 81 % in the typical case, because helices pass through each
other. The module-based morphs keep the chain within about 20 %.

## E. Probing data: does maximum-entropy reweighting recover a known mixture?

For each ensemble, synthetic SHAPE-like profiles were generated from a known mixture: 70 % of the weight on the
conformations of one 2D state and 30 % spread over the others. Noise was 0.1, 0.2 and 0.4 (reactivity scale ~1).
Each profile was fitted from a uniform prior, where that state holds a median 30 % of the weight. Medians over 25
ensembles:

| noise | share of the state after fitting (true 0.70, prior 0.30) | correlation of fitted and true weights |
|---|---|---|
| 0.1 | 0.60 | 0.999 |
| 0.2 | 0.62 | 0.995 |
| 0.4 | 0.58 | 0.991 |

Reweighting moves about three quarters of the way from the prior to the truth. It stops short by design: the
maximum-entropy penalty (θ) keeps the weights close to the prior unless the data demand more. The result is stable
up to a noise of 0.4. These profiles follow the same linear reactivity model the fit assumes, so this test checks
the method and is optimistic for real data.

## F. Confidence model (predicted per-nucleotide lDDT)

**Training data** (`tools/qa_model`). 132 RNA chains of the BGSU non-redundant set (release 4.0 Å; 15–120 nt), with
the RCV engine's models of both modes (1,121 models) and the native structures. Labels are the per-nucleotide
all-atom lDDT against the native structure. The chains are kept apart from the benchmark: no shared equivalence
class, below 80 % sequence identity, and no benchmark family named in the entry title.

**Network.** 41 hand-made features per nucleotide → 48 → 24 → 1 (tanh, sigmoid), 3,217 parameters, numpy, a
fraction of a second per ensemble. Five-fold cross-validation grouped by chain gives a per-nucleotide correlation
r = 0.817 ± 0.018.

**Benchmark** (26 targets, never used for training or model selection):

| | dev (6 targets) | test (20 targets) |
|---|---|---|
| per-nucleotide r (predicted vs true lDDT) | 0.835 | 0.754 |
| within-target ranking, mean Spearman | 0.47 | 0.10 |
| top-1 model by engine score: median RMSD (Å) | **12.8** | 21.5 |
| top-1 model by confidence: median RMSD (Å) | 14.5 | **20.3** |
| confidence better / worse / same than the engine's pick | 2 / 4 / 4 | 15 / 9 / 13 |

Native structures score 0.86 on average and engine models 0.49. The model reliably tells which **parts** of a
model are trustworthy (per-nucleotide r 0.75–0.84 on unseen families). It does not reliably pick the best **model**:
better on the test split, worse on the dev split. RCV 0.7 therefore keeps the engine's ranking and shows confidence
as a per-nucleotide colour scheme and a per-model metric. Re-ranking will be reconsidered when a model wins on
both splits.

## Limitations

* All kinetics are at the secondary-structure level (Turner energies, findpath saddles, Arrhenius rates with an
  assumed prefactor k₀ = 10⁶ s⁻¹). They have not yet been compared with measured folding or switching rates, so
  barriers are upper bounds and times are order-of-magnitude estimates. A comparison with experimental riboswitch
  kinetics is the natural next validation.
* The ensembles come from one prediction engine. Their 2D diversity (1–3 topologies, 1–5 states) is that of the
  engine, not of the RNAs in solution.
* Probing recovery was tested on synthetic data only (section E).
* The cooperative unit and the φ-value scan (section C) are both derived from the same energy model, so their
  agreement is a consistency check, not experimental evidence. In-silico φ-values are a prediction to be tested by
  mutagenesis.
