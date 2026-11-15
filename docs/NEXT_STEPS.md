# Planned functions, not yet implemented

Written 2026-11-12, after package 0.7.4; step 1 of §1 (pocket finding) was implemented in 0.7.5 and is
marked as done in §1.7. Everything else here is a plan: no code for it exists in the tree.
Each item says what it should do, where it would go, what it needs, and how it would be judged. Order
inside each section is roughly the order I would implement in.

The one function in this file that was explicitly asked for next, and is specified in most detail, is
**conformation-selective compound screening** (§1).

---

## 1. Conformation-selective compound screening

**What it is.** The user has an ensemble with two conformations, A and B, that interconvert (a riboswitch,
a thermometer, a pulled state). The question is which compounds of a library bind A much better than B, so
that A is stabilised and the transition to B is blocked (or the reverse). The output is a ranked list of
compounds with a selectivity score, not an affinity prediction: this is a ranking tool to pick a few hundred
compounds out of 100,000 for an experiment, and the documentation has to say so plainly.

Target scale: ~100,000 compounds on a laptop CPU, with a GPU option left for later. This is the scale that
decides the whole design: everything below follows from "100,000 compounds must finish overnight on 16 cores".

**Why it fits RCV.** RCV already has the two things a conformation-selective screen needs and that docking
programs do not have: an ensemble of 3D conformations of one RNA with populations, and the free-energy
difference between the states (`switch.py`, `energies.py`). A screen that scores a compound against A and
against B and reports the difference is a few hundred lines on top of a docking core, and nothing else in
the field packages it for a bench scientist.

### 1.1 The sub-window

A screen is a long job with many parameters, so it gets its own window (like the start page, not like a
dialog): **Tools ▸ Compound screening…**, and a button in the Switch tab next to the chosen A and B.
`src/rcv/static/js/screening.js`, served at `/_screen`, opened by `windows.py` as a second webview;
falls back to a full-page view in the single-window build. The window is not part of the main layout,
so the main interface does not grow (the user's constraint: keep the interface as it is).

Four cards, top to bottom:

1. **Target.** The ensemble, the two conformations (A and B, prefilled from the Switch tab), and the
   pocket. Pocket choice: *whole RNA* (blind), *the pockets found automatically* (list, each with volume
   and the nucleotides lining it, clicking one shows it in the main window's 3D view), or *the current
   selection in the main window*. A screen of the whole RNA at 100,000 compounds is too slow, so if the
   user leaves it blind the window says how long it will take and offers the pocket list instead.
2. **Library.** A file picker (SDF, SMILES, MOL2, or a prepared `.rcvlib` pack), a recent-libraries list,
   and *Download a standard library…* (§1.5). Shows the compound count, how many were parsed, how many
   were rejected and why. No library ships in the installer.
3. **Parameters.** Exhaustiveness (*quick* / *standard* / *thorough*), the number of poses kept per
   compound and state, the selectivity threshold for the short list, the ligand conditions
   (pH for protonation, Mg²⁺), and whether to keep poses for all compounds or only the short list
   (100,000 compounds × 2 states × poses is gigabytes; default is to keep poses for the top 2,000 only).
4. **Run.** Start, pause, resume, stop. A progress row with compounds done, compounds per second,
   time left, and a live histogram of the selectivity score. Partial results are usable while it runs:
   the table fills from the top.

Results, in the same window: a sortable table (compound id, score for A, score for B, selectivity,
pocket, strain, a 2D depiction of the compound), filters (selectivity, molecular weight, logP, rings,
PAINS-like substructures), *Show in 3D* which loads the pose into the main window as a ligand with the
contacting nucleotides selected, and export (CSV, SDF with poses, a `.rcv` with the best poses as scenes).
A report page for the lab: the ten best compounds, their contacts, the pocket, and the parameters used.

### 1.2 The scoring core

Written in `src/rnafold3d/screening/`, importable and usable from the command line without the GUI
(`rnafold3d screen --ensemble x.rcv --a conf_1 --b conf_2 --library lib.sdf --out run/`), because a
100,000-compound screen belongs on a workstation or a cluster as often as on the laptop.

Pipeline per compound:

1. **Prepare** (once per compound): parse, add hydrogens at the chosen pH, enumerate 1–8 ring-conformer
   and tautomer variants, find the rotatable bonds, compute the atom types and partial charges.
2. **Pocket grids** (once per state, not per compound): the pocket of A and the pocket of B are each put
   on a 0.375 Å grid holding, per atom type, the sum of the pairwise terms with the RNA — this is what
   makes 100,000 compounds feasible, since the expensive part is then a grid lookup.
3. **Search**: rigid-body placement plus torsion sampling of the compound in the grid, by Monte Carlo with
   local gradient minimisation; a few hundred evaluations for *quick*, a few thousand for *thorough*.
   Deterministic for a given seed, like the rest of the engine.
4. **Score** (the part that has to be RNA-specific, not borrowed from protein docking): van der Waals,
   directional hydrogen bonds, base stacking as an explicit term with a face-to-face geometry (the
   dominant interaction of RNA-binding small molecules, and the term protein scoring functions get worst),
   electrostatics with a distance-dependent dielectric and explicit backbone phosphate charges, the
   desolvation of the phosphate backbone, Mg²⁺-mediated contacts, and a torsional-entropy penalty.
5. **Selectivity**: score(B) − score(A) for the best pose in each, plus the strain (how much the pocket
   of A has to be deformed, which flags scores that come from a pocket that would not stay open).
   Two compounds with the same selectivity are ranked by the worse of their two scores, since a compound
   that binds neither state selectively is useless.

Parallelism: one worker per core over compounds, shared read-only grids, a results file appended in
compound order so a stopped run resumes. Memory stays flat (grids plus one compound per worker).

Budget: 100,000 compounds × 2 states must fit in a night on 16 cores. That is ~0.3 s per compound per
state at *standard*, which a grid-based search reaches for a 25-atom compound with 5 rotatable bonds,
but only if the pocket is a pocket and not the whole RNA. If the prototype misses this, the fallback is a
two-stage screen: a cheap filter (pharmacophore and shape match against the pocket, ~5 ms per compound)
over all 100,000, then the full search on the best 10,000.

### 1.3 Pocket finding — implemented in 0.7.5

What was built follows this section, with two decisions made during the work: sites are detected at 6 of 7
enclosing directions (at 5 the clefts of an RNA merge into one cavity spanning the molecule), and two
conformations are compared by measuring the open volume held by the *same nucleotides* in both, rather than
by matching whole pockets — a pocket that merges with its neighbour in the other conformation has no match,
but its nucleotides still hold a measurable volume. The ranking by the difference between A and B is in the
Details tab; the screen will take the pockets from the same module. What is described below and was *not*
built: the pocket list is not yet validated against known RNA-ligand complexes.

Original plan — `screening/pockets.py`: a grid flood-fill for buried volume (the classic LIGSITE/fpocket approach applied
to RNA), then keep pockets of 100–1,500 Å³, score them by burial, by how many nucleotides are
non-canonically paired or unpaired (RNA ligands bind loops, bulges and junctions, not helices), and by how
much they differ between A and B. A pocket that exists in A and is closed in B is exactly what a
conformation-selective screen wants, so the list is sorted by that difference by default. This is also
useful on its own in the Details tab, so `pockets.py` is written as a standalone analysis first and the
screen uses it.

### 1.4 Validation — the part that decides whether this ships

A screening function that cannot be shown to work is worse than no screening function, because the user
cannot tell. Three benchmarks, in `benchmark/screening/`:

1. **Pose reproduction.** The RNA–ligand complexes in the PDB with a resolution better than 3.0 Å and a
   drug-like ligand (roughly 150–250 complexes; the HARIBOSS set is the curated starting point). Redock
   the ligand into its own structure and measure the RMSD of the best-scored pose. Target: a median below
   2.5 Å and the right pose in the top 5 for half the cases. Published RNA docking programs reach about
   this; being clearly worse means the scoring function is not ready.
2. **Enrichment.** For the RNA targets with many known actives (the FMN, TPP and ZTP riboswitches,
   the HIV-1 TAR element, the pre-miR-21 and MALAT1 constructs from the published screens), mix the actives
   into 50–100 property-matched decoys per active and report the enrichment at 1% and the AUC. A random
   ranking is AUC 0.5; below about 0.65 the screen should not be offered.
3. **Selectivity, the one that matters and the one with almost no data.** Compounds known to stabilise one
   state of a switch (ribocil on the FMN riboswitch, the published TPP riboswitch ligands, risdiplam-type
   SMN2 splice modifiers as a hard case) should rank in the top percent by selectivity against the
   alternative state. There will be fewer than twenty such cases, so this is an anecdote rather than a
   benchmark; the manual must say that selectivity is the least validated of the three.

Until all three exist and are reported in `docs/SCREENING_STATUS.md`, the function stays behind
*Advanced settings ▸ experimental* with the measured numbers shown in the window itself.

### 1.5 Libraries

Nothing ships in the installer. `tools/fetch_library.py` downloads and prepares a library into
`~/.rcv/libraries/NAME.rcvlib`, with the compound records, the precomputed conformers and an index:

| Library | Compounds | Prepared size | Why |
|---|---|---|---|
| Approved and clinical-stage drugs | ~4,000 | ~60 MB | repurposing, the first screen anyone should run |
| An RNA-focused subset of a commercial catalogue | ~20,000 | ~300 MB | the compounds most likely to bind RNA |
| A general diversity set | ~100,000 | ~1.5 GB | the stated target scale |
| User's own file | any | — | SDF, SMILES, MOL2 |

The download link is the user's to provide (his own hosting or a Zenodo deposit); the tool takes a URL and
a checksum from a small manifest that ships with RCV, so a library can be added without a new release.
A `.rcvlib` is a plain zip with a documented layout, so a lab can prepare one and pass it around.

### 1.6 Dependencies, and what to avoid

The honest options for compound handling are RDKit (a 100+ MB dependency, permissive licence, and the only
real choice for tautomers, protonation and 2D depiction) or writing a minimal SDF/SMILES reader plus a
conformer generator (small, no dependency, and a long tail of chemistry bugs). Recommendation: RDKit, as an
optional dependency that the screening window asks to install, so the base installer stays its current
size and a user who never screens never downloads it. The scoring core itself must not depend on RDKit —
it takes prepared atoms, types and charges — so the prepared-library path works without it.

Not to be done: no reuse of a GPL docking program's code or parameters (AutoDock Vina's scoring function
is Apache-2.0 and may be used as a reference, rDock and most others are GPL and would change RCV's
licence); no claim of affinity prediction in kcal/mol in the interface; no claim that a compound "will"
block a transition.

### 1.7 Steps

1. **Done (0.7.5).** `screening/pockets.py` + the Details-tab pocket list, with the pocket-difference score
   (`compare` / `volume_near`: the same nucleotides measured in both conformations, rather than matching whole
   pockets, which fails when a pocket merges with its neighbour). Sites are detected at 6 of 7 enclosing
   directions, not the usual 5: at 5 the open clefts of an RNA merge into one cavity spanning the molecule.
   Validated in `benchmark/screening/pocket_sites.py` on thirty experimental complexes with the compound
   removed: the site it occupied is the top-ranked pocket in 21 of 30 and found in 24 (at 6 of 7: 9 and 12,
   which is why the default threshold is 5). Still the easy case — the compound shaped the site it is
   measured in.
2. **Done (0.7.5), with benchmark 1 run.** The scoring function (`screening/score.py`), its grid
   (`grid.py`), the compound and receptor preparation (`chem.py`, `receptor.py`) and the search
   (`dock.py`), with `rnafold3d dock` on top and `benchmark/screening/redock.py` measuring it. Numbers and
   limits: `docs/SCREENING_STATUS.md`. What the work changed about the plan: the stacking term is stored on
   the grid as the matrix Σ w·nnᵀ so that the angle between the compound's ring and the base survives the
   gridding; the search runs on the grid but every kept pose is rescored from the atoms; the finite-difference
   step of the minimiser has to be about a tenth of a grid cell, not a tiny number, or it follows the kinks of
   the interpolation. Still missing here: analytic gradients and a compiled inner loop, without which a
   compound takes seconds rather than the 0.3 s the 100,000-compound budget assumes.
3. **Done (0.7.5).** Benchmark 2 (enrichment): `benchmark/screening/enrichment.py`, six RNA targets with
   actives collected from the literature and 25 property-matched decoys per active drawn from PubChem,
   with everything 2D-similar to an active removed. **Mean AUC 0.864** (34 actives, 700 decoys; the plan
   asked for 0.65), five targets between 0.855 and 0.997 and the TPP riboswitch at 0.495 — exactly chance,
   and the one site in the set that holds its ligand through bridging Mg²⁺. Full table, and why the A
   site's 0.972 is the least trustworthy of the good numbers: `docs/SCREENING_STATUS.md` §3. The weights
   were tuned against this benchmark, as the rule allows, and **left unchanged**: +0.010 cross-validated
   with three of the four fitted weights pinned to the edge of their allowed range is not a result. Two things the
   work changed about the plan: the raw
   value of every term is stored for every kept pose, so a change of weights is a dot product rather than
   a re-run, which is what makes `tune` affordable; and tuning is judged leave-one-target-out, because
   with six targets an in-sample gain means nothing. A charge bug found by it (an aromatic nitrogen with
   three single bonds was taken for a pyridinium, so theophylline came out +2 and riboflavin +1) was
   fixed, and benchmark 1 re-run.
4. **Done (0.7.5).** The screening window (`src/rcv/static/js/screening.js`, `src/rcv/screening.py`,
   `/api/screen/*`): target, library, parameters, run, and the ranked table with *Show in 3D*. It is a
   page of its own rather than a second OS window, which is the same thing for the user and does not need
   a second webview; the main layout is untouched. Manual: `docs/RCV.md` §8.7.
5. **Done (0.7.5)** except benchmark 3. The two-state screen is `screening/screen.py` and
   `rnafold3d screen`: selectivity = score(B) − score(A) at the site of A measured in both conformations,
   ties broken by the worse of the two scores. Benchmark 3 (the published conformation-selective
   compounds) is still not run, and until it is, the window says so.
6. `tools/fetch_library.py`, the `.rcvlib` format and the manifest.
7. The experimental flag, removed only if benchmark 3 and a faster search support it.

A rough size: 2,500–4,000 lines of Python, 800–1,200 lines of JavaScript, and the benchmarks are the
larger half of the work.

---

## 2. Engine accuracy (rnafold3d)

Each one is measured on the same two benchmarks as 0.7.4 (`benchmark/template_selection.py`,
`benchmark/dev_templates.py`) with CASP16 kept as the blind check.

**Status (0.7.5, engine 2.3.0; numbers in `CHANGELOG.md` and `docs/V2_STATUS.md` §4.0000).** Fragments: done and
on — the only item that moved the dev benchmark (no template 0.216 → 0.227). Base–base potential: built and
measured, not in the ranking (no gain: the limit is sampling). Pseudoknots: in the search as an option, off — the
2D is right more often but the helix rarely forms in 3D. Speed: numpy-level only; a compiled inner loop is still
open. Long RNAs: domain assembly as an experimental option, off — worse than the whole-chain search on six chains
of 255–393 nt. More data types: done (SAXS, FRET, cryo-EM maps in the viewer). Complexes: waiting, as planned.

* **Motif-level 3D fragment templates.** The library currently gives whole-chain templates and junction
  geometries. Between them sits the fragment: a hairpin loop, a bulge, a k-turn, a tetraloop–receptor,
  taken from any library chain whose sequence matches the motif, and grafted onto the fold tree. This is
  what FARFAR2 and Vfold do, and it is the single most likely gain for targets without a template — which
  is where RCV is weakest (CASP16 without a template: 0.26).
* **A knowledge-based base–base potential.** Statistics over the library of the relative position and
  orientation of base pairs by type and separation, as a smooth potential in the ranking score. Replaces
  hand-set tertiary terms with measured ones; the risk is double-counting with the junction restraints,
  so it must be fitted with those restraints present.
* **Pseudoknots in the search, not only in the analysis.** 0.7.4 finds pseudoknot helices in a structure
  but does not fold them: the 2D hypotheses come from ViennaRNA, which cannot. A pseudoknot-capable 2D
  (an exact algorithm for the simple classes, or a heuristic plus the template pairs) would fix the
  targets where the global fold is wrong because a long-range pair is missing.
* **Speed.** The fold-tree inner loop in numba or as a small C extension. Not accuracy, but it buys
  trajectories, and trajectories are accuracy at fixed quality: the number of valid trajectories is
  already the limit on the coarse-template targets (R1256: 6 of 40).
* **Complexes.** Two chains (a duplex, an RNA–RNA kissing loop) — the engine is single-chain throughout,
  so this is a larger change than it sounds and should wait.
* **Long RNAs.** A domain finder that cuts a 500-nt RNA into independently folded domains at the junctions
  the 2D makes obvious, folds each, and reassembles. The current limit is practical, not principled.
* **More data types for reweighting.** Beyond SHAPE/DMS: SAXS profiles, FRET distances, cryo-EM maps as
  restraints, each with the same shape as the probing path that 0.7.4 added.

## 3. Interface and retention

* ~~**Project library.**~~ *Done in 0.7.5 (File ▸ Project library, `rcv/library.py`).* The one thing that makes a tool sticky: a folder RCV watches, with every run,
  structure and `.rcv` the user has opened, searchable, with a thumbnail and the parameters, so coming
  back after a month is one click rather than hunting for a file.
* **Template alerts.** The library in the installer is a snapshot. A check against the PDB for new
  structures related to an RNA the user has worked on, and an offer to re-predict, turns RCV into
  something that gives the user news.
* ~~**Recipes.**~~ *Done in 0.7.5 (Tools ▸ Recipes, `rcv/recipes.py`, four built-in).* Named parameter sets for a whole analysis ("riboswitch report", "pull and compare"), saved
  and shared as a small file, so a lab standardises how it uses RCV.
* ~~**A plain-language console.**~~ *Done in 0.7.5 (`console/plain.js`, `help plain`).* The console already does everything; a natural-language front end over it
  (local, pattern-based, not a model) lowers the entry cost for the scientists who will never type
  `conditions 25, 150 mM`.
* 
* 

## 4. Layout

0.7.4 tightened the chrome (menu bar, toolbar, panel heads, timeline, tab padding) and widened the side
panels slightly above 2,200 px. **Done in 0.7.5** (`static/js/windows.js`, user guide §4.3):

* ~~Draggable splitters on both side panels, with the widths saved per screen size.~~ Drag the inner edge
  (left 190–520 px, right 312–760 px, the 3D view keeps ≥ 360 px); double-click = default; saved per
  `screen.width × screen.height`; *Window ▸ Default panel widths*.
* ~~A collapsed-panel rail.~~ A hidden side panel folds to a 36-px rail; its icons open the panel (on the right,
  on that tab, plugin tabs included).
* ~~A *presentation* mode.~~ `F5` / Window menu: 3D only, full screen, no menu bar, large labels, title and
  legend; `Esc`/`F5` or the button that appears on mouse movement leaves it.

