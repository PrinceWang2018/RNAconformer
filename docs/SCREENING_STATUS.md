# Compound screening: what works, what is measured, what is not

This file says plainly how far the screening part of RCV has got and what the numbers behind it are. It is
the file the plan (`docs/NEXT_STEPS.md` §1.4) asks for, and the rule it states is kept here: **a screening
function that cannot be shown to work is worse than none, because the user cannot tell.**

Status as of 0.7.5 (2026-11-14): steps 1 to 5 of the plan exist — pocket finding, a scoring function, its
grid, a search, the two-state screen, and the screening window — with two of the three benchmarks run. In
one line each: redocking puts the right pose first in 10 of 30 experimental complexes (median 4.24 Å,
below what the plan asks); ranking known binders against property-matched decoys gives a mean AUC of
**0.864** over six targets, clearing the plan's 0.65 on five of them and scoring exactly chance on the
sixth. What is still missing: benchmark 3 (the published conformation-selective compounds, the thing the screen is
actually for), prepared libraries, and the speed the hundred-thousand-compound target assumes. The screen is
marked **experimental** in the interface and shows the numbers below next to its results. Nothing in RCV
reports an affinity, nothing is in kcal/mol, and nothing claims a compound will bind or block anything.

---

## 1. What exists

| Piece | Where | What it does |
|---|---|---|
| Pocket finding | `rnafold3d/screening/pockets.py` | buriedness scan, pockets with volume, lining and a score; the difference between two conformations |
| Compounds | `screening/chem.py` | SDF/MOL/MOL2/PDB readers, rings, aromaticity, rotatable bonds, torsion tree, charges at pH 7, SDF writer |
| Receptor | `screening/receptor.py` | the RNA typed from its atom names, one −1 per phosphate, a plane per base; splitting a PDB entry into RNA and compounds |
| Scoring function | `screening/score.py` | gauss1, gauss2, repulsion, hydrophobic, hbond (Vina's form, Apache-2.0, as a reference), plus stacking against the base planes, screened electrostatics, desolvation, metal contacts, and the torsion-entropy division |
| Grid | `screening/grid.py` | every term precomputed on a 0.375 Å grid, one interpolation per pose; stacking kept directional by storing Σ w·nnᵀ per point |
| Search | `screening/dock.py` | Monte Carlo with local minimisation over position, orientation and torsions; deterministic for a seed; poses rescored exactly |
| Command line | `rnafold3d pockets`, `dock`, `screen` | pockets of a structure; docking into a pocket; a library against two conformations |
| Two-state screen | `screening/screen.py` | a library against two conformations: selectivity = score(B) − score(A) at the site of A, measured in both; one worker per core |
| Window | `rcv/screening.py`, `static/js/screening.js` | *Tools ▸ Compound screening…*: target, library, parameters, run, ranked table, pose back in the 3D view |
| Benchmark 1 | `benchmark/screening/redock.py` | redocking of experimental RNA–ligand complexes |
| Benchmark 2 | `benchmark/screening/enrichment.py` | known binders against property-matched decoys, on six targets |

The scoring function is not a potential of mean force and the number it produces is not an energy. It ranks.

## 2. Benchmark 1 — redocking (the one that is run)

Take an experimental RNA–ligand complex, remove the compound, dock it back into its own structure, and
measure how far the best-scored pose is from the experimental one. Under 2.5 Å is the usual line for "the
right pose". The RMSD lets atoms of the same element swap, so a compound flipped by its own symmetry counts
as the same pose (a lower bound of the proper graph-symmetry RMSD).

The set is **hand-picked**, not the curated HARIBOSS set the plan asks for: about thirty well-known
complexes — riboswitches (adenine, guanine, SAM I/II/III/V, TPP, preQ1, c-di-GMP, glycine), aptamers
(malachite green, theophylline, FMN, tobramycin, neomycin), the HIV-1 TAR element with argininamide, and
A-site RNA with aminoglycosides. It spans 10 to 50 heavy atoms, 0 to 10 rotatable bonds and charges from
−2 to +5. It is enough to see whether the function is in the right place and to catch a regression. It is
not enough to publish a number, and it is the easy case: the pocket already has the right shape, because
the compound shaped it.

Run it with:

    python benchmark/screening/redock.py fetch                 # downloads from RCSB into ~/.rcv
    python benchmark/screening/redock.py run --out benchmark/runs/screening/redock_standard.json

### Results (0.7.5, `--exhaustiveness standard --seed 1`, 91 minutes on one core)

`benchmark/runs/screening/redock_standard.json`, re-run after the charge-perception fix below.

| | best-scored pose | any of the kept poses |
|---|---|---|
| median RMSD | **4.24 Å** | 3.50 Å |
| under 2.5 Å | **10 of 30 (33%)** | 11 of 30 (37%) |

**The charge-perception fix moved individual cases a long way and the median not at all.** While building
benchmark 2 it turned out that any aromatic nitrogen with three single bonds was being read as a
pyridinium, so theophylline came out +2 and riboflavin +1 — a whole unit of charge wrong in front of a
polyanion. A nitrogen is now charged only when it also holds a double bond inside its ring. Redocked with
the corrected charges, against the same run with the old ones:

| | fixed | broke |
|---|---|---|
| | 1F1T malachite green 9.39 → **0.32 Å**, 2GDI TPP 10.69 → **2.80 Å** | 4TS2 0.66 → 6.84 Å, 3D2G 2.25 → 8.84 Å, 1J7T 2.14 → 8.51 Å, 2TOB 2.50 → 5.87 Å |
| median | | 3.97 → 4.24 Å, under 2.5 Å 11 → 10 of 30 |

The fix is right on its own terms — those two compounds were being docked as the wrong molecule — and it
changed the summary by less than the spread between neighbouring cases. That is the honest reading of a
30-case set: it has no power to resolve a change of this size, and a difference of one case in thirty is
not evidence either way.

**This is below what the plan asks for** (a median under 2.5 Å and the right pose in the top five for half
the cases) and below what published RNA docking programs report. The function is therefore offered as
**experimental** only: the screening window and the command line both print these numbers next to their
results, and nothing in RCV calls the score an affinity.

What the numbers say about *why*:

* In **24 of 30** the search found a pose scoring at least as well as the experimental one. So in most of
  the failures the search is doing its job and the **scoring function prefers the wrong pose** — the ranking
  is what needs work, not the sampling.
* The five worst cases, redocked again at `thorough` (`redock_thorough_worst.json`, measured with the
  earlier charge perception — two of the five, 2GDI and 1F1T, are now solved at `standard`), split exactly
  along that line:

  | case | standard | thorough | what it was |
  |---|---|---|---|
  | 2GDI thiamine pyrophosphate | 10.69 Å | **1.20 Å** | sampling: the experimental pose scored better all along |
  | 1F1T malachite green | 9.39 Å | **0.32 Å** | sampling: likewise (score −11.0 against −9.6) |
  | 5BTP purine analogue | 8.04 Å | 9.11 Å | ranking: our pose scores −8.5, the experimental one −5.6 |
  | 1ARJ argininamide (NMR) | 6.64 Å | 6.86 Å | ranking: −5.8 against −2.1 |
  | 1BYJ paromomycin (NMR) | 8.53 Å | 8.14 Å | ranking: −8.1 against −6.2 |

  More search fixes the cases where the experimental pose already scored better, and does nothing for the
  cases where the function prefers a wrong pose — which is the whole argument for benchmark 2 being the
  next piece of work rather than a faster search. Cost of `thorough`: 2 to 21 minutes per compound.
* It works where the site is a tight pocket around a small rigid compound: 2B57 0.27 Å, 1F1T 0.32 Å,
  3SKI 0.61 Å, 1Y26 0.87 Å, 4NYB 0.89 Å, 1FMN 1.56 Å, 1U8D 1.71 Å, 3Q50 1.71 Å, 1O15 1.86 Å, 4FE5 1.91 Å.
  Every one of them is under 20 heavy atoms with at most 6 rotatable bonds, except FMN.
* It fails where the compound is big and flexible in a shallow site — every aminoglycoside and every
  SAM entry is between 3.9 and 9.0 Å — and on the two NMR entries whose experimental pose scores
  *positive* under these radii (1LVJ +19.7, 2KGP +8.6): close contacts in those models that the repulsion
  term cannot accept. Four cases never found a pose scoring as well as the experimental one, so there the
  search itself ran out: 3IRW c-di-GMP (46 atoms, 12 rotatable bonds, −4.0 against −11.0), 2QWY and 3GX5
  SAM, 1NEM neomycin.
* Cost: 10 s for a rigid 10-atom compound, 3–10 minutes for 40 atoms with 9–12 rotatable bonds. That is
  100 to 1000 times slower than the screen's budget needs (§6).

The weights are not tuned against this — the plan's rule is that only the enrichment benchmark of §3 may
tune them, and when it was allowed to (§3), it did not earn a change. What this run is for is to know
where the function stands and to catch a regression in it.

## 2b. Do the pockets cover the sites compounds use?

`benchmark/screening/pocket_sites.py` takes the same thirty complexes, **removes the compound**, runs the
pocket finder on the RNA alone and asks where the site the compound occupied appears in the ranking. A site
counts as found when a pocket's lining covers at least half of the nucleotides the compound touches
(within 4.5 Å).

| enclosed on | site = top pocket | site in the top 3 | site found at all | pockets per structure (median) |
|---|---|---|---|---|
| **5 of 7** (default) | 21 / 30 | 23 / 30 | 24 / 30 | 2 |
| 6 of 7 (*deep pockets*) | 9 / 30 | 12 / 30 | 12 / 30 | 1 |

This is what set the default threshold. Six of seven gives tidier, better separated pockets and misses the
shallow sites, which is what most RNA sites are — the aptamers and the NMR structures in the set have
almost none at six. The six missed at five of seven (1Y26, 1FMN, 1EHT, 1AJU, 1LVJ, 2KGP) are surface sites
with almost no enclosure at all.

The caveat is the same as everywhere here: these are the compounds' own structures, so the site is already
open and shaped. On a predicted model, or on the conformation of a switch that does not bind, the numbers
say nothing — which is exactly why the two-conformation comparison exists.

    python benchmark/screening/pocket_sites.py [--burial-min 5|6|7] [--out runs/pocket_sites.json]

## 3. Benchmark 2 — enrichment (the one that says whether the ranking is worth anything)

Redocking says whether the search can find the right pose. It cannot say whether the *ranking* is any good,
because every compound in it is a binder. This benchmark asks the question a screen actually has to answer:
given one target and a pile of compounds, do the ones that bind come out on top?

`benchmark/screening/enrichment.py`. Six targets, each an experimental structure with the pocket taken from
where its own ligand sat, the ligand removed:

* **Actives**: compounds reported in the literature to bind that RNA, collected by hand and named in the
  file — there is no RNA equivalent of DUD-E. Between four and twelve per target, which is few enough that
  one wrong assignment moves the number.
* **Decoys**: 25 per active, drawn from a random sample of 30,000 PubChem compounds and matched to each
  active on molecular weight, lipophilicity, polar surface area, rotatable bonds and hydrogen-bond counts,
  with everything PubChem calls 2D-similar (≥ 85%) to any active of that target thrown out. Property
  matching is what stops a score from winning on molecular weight alone. What it cannot give is compounds
  *measured* to be inactive: a matched random compound is assumed inactive, which is the standard
  assumption of this kind of benchmark and is wrong for a few percent of them.

    python benchmark/screening/enrichment.py fetch            # actives and decoys, from PubChem and RCSB
    python benchmark/screening/enrichment.py run --jobs 12 --out benchmark/runs/screening/enrichment_quick.json
    python benchmark/screening/enrichment.py report benchmark/runs/screening/enrichment_quick.json
    python benchmark/screening/enrichment.py tune  benchmark/runs/screening/enrichment_quick.json

### Results (0.7.5, `--exhaustiveness quick --seed 1`)

| target | structure | actives | decoys | AUC | EF 1% | EF 5% | property gap (σ) |
|---|---|---|---|---|---|---|---|
| FMN riboswitch | 3F2Q | 4 | 95 | **0.855** | 25 | 10.0 | 0.74 |
| TPP riboswitch | 2GDI | 4 | 100 | **0.495** | 0 | 0.0 | 0.46 |
| adenine riboswitch | 1Y26 | 4 | 98 | **0.997** | 25 | 20.0 | 0.58 |
| guanine riboswitch | 1U8D | 6 | 148 | **0.871** | 0 | 6.7 | 0.64 |
| preQ1 riboswitch | 3Q50 | 4 | 99 | **0.992** | 25 | 15.0 | 0.79 |
| ribosomal A site | 1J7T | 12 | 160 | **0.972** | 17 | 10.0 | 3.38 |
| **mean** | | 34 | 700 | **0.864** | 15.3 | 10.3 | |

49 minutes on 12 cores, nothing undockable. A random ranking scores AUC 0.500 and EF 1.0; the plan asked
for at least 0.65, and five of the six targets clear it.

Read the numbers with their limits in view:

* **TPP is a failure, not a near miss.** AUC 0.495 is exactly chance, and not one active lands in the top
  5%. TPP is the most charged ligand in the set (a diphosphate, formal −2 at pH 7) binding a site that
  holds its phosphates through a pair of Mg²⁺; the score there is dominated by terms this function treats
  crudely — a 6-point metal term and a screened-Coulomb approximation with no explicit ions. Compounds with
  loose phosphates outscore the real ligand. Anything with a metal-bridged phosphate site should be
  assumed to rank no better than random.
* **The A site's 0.972 is the least trustworthy good number.** Its decoys sit 3.38σ from its actives in
  property space (every other target is under 0.8σ), because PubChem holds few 3D records for compounds as
  large and polar as the aminoglycosides — the decoy pool simply does not contain their match. Part of
  that AUC is the score preferring big polycations, which on this target happens to be the right answer.
  It also has only 160 decoys rather than the 300 the ratio asks for, for the same reason.
* **EF 1% rests on one or two compounds.** With four actives in ~100 compounds the top 1% is a single
  slot, so EF 1% is 25 or 0 and nothing in between. The AUC is the number to read; the EFs are reported
  because a screen is used at its top end, not averaged.
* Five of six targets are nucleobase or flavin pockets, where stacking against a base plane does most of
  the work. The set is not evidence about pockets of another kind.

The plan's rule is kept: **the weights of the scoring function may only ever be tuned against this
benchmark**, never against the selectivity cases of §4. `tune` searches them by reweighting the stored
terms of every kept pose (a dot product, not a re-run) and reports leave-one-target-out: fit on five
targets, measure on the sixth.

**It was run, and the weights were left alone.** Allowing the four RNA-specific terms (stacking,
electrostatic, desolvation, metal) to move within a factor of 4 of their defaults, with signs fixed:

| weights | mean AUC | cross-validated |
|---|---|---|
| defaults, as shipped | 0.864 | 0.864 |
| fitted on all six targets | 0.888 | 0.874 |

A cross-validated gain of +0.010 is not a reason to change a scoring function. Three signs point the same
way: three of the four fitted weights ran to the edge of the allowed box and the fourth close to it
(stacking 4.00×, desolvation 0.25×, metal 0.25×, electrostatic 3.14×), which is a fit leaning on a
boundary rather than finding an optimum; the
gain is concentrated in two targets (TPP +0.090, A site +0.025) while two others get worse (guanine
−0.043, preQ1 −0.010); and with 34 actives in total that spread is within what one reassigned active would
move. On top of that, reweighting *stored* poses only explores weights near the ones the poses were docked
with — a stacking term four times heavier would have found different poses, and would almost certainly
cost redocking accuracy (§2), which the surrogate cannot see. The shipped weights stay as they are, and
the tuning run is recorded here rather than applied.

## 4. Benchmark 3 — selectivity (not run)

Compounds known to stabilise one state of a switch (ribocil on the FMN riboswitch, the published TPP
riboswitch ligands, risdiplam-type SMN2 modifiers) should rank in the top percent by selectivity against
the alternative state. There are fewer than twenty such cases: this will be an anecdote, not a benchmark,
and the manual has to say so.

## 5. What the weights are, and that they are not fitted

    gauss1 −0.0356   gauss2 −0.00516   repulsion 0.840   hydrophobic −0.0351   hbond −0.587
    stacking −0.250  electrostatic 0.050   desolvation 0.110   metal −0.700   rotatable 0.0585

The first five and the last are AutoDock Vina's, used as a reference (Apache-2.0). The four RNA terms are
**set, not fitted**: they put a face-to-face ring stack at about −1.4, a salt bridge at about −1 and a
buried charge at about +1, which is the scale of a hydrogen bond in the same units. Benchmark 2 is what
would fit them; it now exists, it was run against them, and it did not justify moving any of the four
(§3). They remain set by hand.

## 6. Known limits

* The RNA is rigid. A screen against two conformations is the answer this tool gives to that, not induced fit.
* No water, no explicit Mg²⁺ geometry beyond a distance term, no entropy beyond the torsion count.
* Tautomers, protonation beyond the groups listed in `chem.ionised`, and ring conformers are not enumerated:
  a compound is docked as the file draws it. RDKit is the intended answer and is not a dependency.
* The search is Python; a compound takes seconds, not the 0.3 s the 100,000-compound target needs. The
  analytic gradient and a compiled inner loop are the next performance step.
* Pocket finding is validated only on the compounds' own structures (§2b), which is the easy case.
