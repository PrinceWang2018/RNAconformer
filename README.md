# RCV — RNA Conformation Viewer

![RCV](docs/images/rcv-wordmark.svg)

**Predict, explore and present RNA conformations on your own computer.** RCV combines the rnafold3d
prediction engine (physics-inspired, **no deep learning in the prediction**) with an interactive desktop viewer
and analyses of how an RNA switches between its conformations. Free and open source (GPL-3.0).

**Author:** Zixiang Wang · princewang2018@gmail.com

* **Start page: the two ways to use RCV** *(0.7.1)* — **1 · predict RNA conformations from sequence(s)**: paste a
  sequence (empty = demo sequence) or a **multi-sequence FASTA for batch prediction** (live stages, time left,
  notes; choose which result to open), settings under *Advanced settings*; **2 · view experimental or predicted
  structures of one RNA**: open NMR ensembles, MD snapshots, models of any predictor, earlier runs or `.rcv` files
  (several files open together as one ensemble), or **fetch a PDB entry**.
* **A desktop program, not a web page** — mouse tools (select, pull, paint, label, measure), **undo / redo /
  clear**, RCV's own right-click menus, no accidental text selection; the command-line tools are in the
  **Tools** menu.
* **One file per RNA: `.rcv`** — every conformation with energies, 2D structures, validation, provenance and
  your view, 10–14× smaller than the same models as PDB; open format (ZIP + JSON), `rcv convert/info/extract`.
* **Pull to switch** — drag a residue: the structure stretches and switches to the nearest stable computed
  conformation; *Explore* runs the engine with the pull as a temporary bias to find new stable ones.
* **2D window linked to 3D** — the secondary structure of the shown conformation, coloured like the 3D view,
  with linked hover/selection, morphing in step with the 3D morph.
* **Switch analysis** — Boltzmann populations of the ensemble's 2D states (any temperature) and a mutation scan
  that ranks the point mutations shifting the switch towards A or B, with one-click *Predict in 3D*.
* **Pathways: how A turns into B** *(0.7)* — refolding paths and barriers between the basins of A and B (ViennaRNA
  findpath), routes through intermediates ranked by **transition-path-theory flux**, the mean first-passage time
  to switch, traps and a feasibility verdict; the **cooperative unit** of base pairs that must open together before
  the saddle; **rigid modules** and hinges between the two 3D models with the **rate-limiting coordinated motion**;
  an **in-silico φ-value scan** of every point mutation (which ones speed the switch up or slow it down); and a
  **mechanism morph** that plays the motion in order. Validated on the benchmark ensembles:
  [docs/V07_VALIDATION.md](docs/V07_VALIDATION.md).
* **Topology-aware ensemble map** *(0.7)* — every secondary-structure topology in its own disk, conformations placed
  by RMSD inside, routes drawn on top; quality reported against classical MDS.
* **Conditions, ligands and probing data** *(0.7)* — temperature, Na⁺ and ligand binding to an aptamer (switch
  dose–response); **SHAPE/DMS** profiles: agreement of every conformation and maximum-entropy **reweighting** of
  the ensemble.
* **Confidence model** *(0.7, optional)* — a small neural network (numpy, runs on any laptop) predicts the local
  accuracy of every nucleotide (lDDT) and colours the structure by it.
* **Command console** — an RNA-ensemble selection language (`helix H2`, `switching`, `rmsf > 3`, `within 5 of …`)
  and familiar molecular-viewer commands (select, show/hide, color, label, zoom/orient, distance/angle/dihedral,
  align/fit, scenes, png, fetch, scripts), plus a sequence bar, measurements and saved scenes.
* **Publication output** — figures (PNG/TIFF/PDF, vector energy panel), **multi-conformation figure panels**
  with 2D drawings, a **plot studio** (clustered RMSD matrix, RMSD to a reference, RMSF, funnel, energy profile;
  journal sizes in mm, vector PDF/SVG), **annotations** (labels, region brackets, notes), **colour palettes**
  (colour-blind-safe sets, your own, with a colour-vision check), **molecular surfaces** and per-object
  **transparency**, videos including **transition A → B** (directly or through intermediates), an
  **interactive single-file HTML page** for supplements and reviewers, and an auto-generated **methods text**.
  Figure text uses the bundled open-licence font Liberation Sans (SIL OFL), embedded in PDFs.
* **In your language, your layout** — interface in English, 简体中文, 繁體中文, 日本語, 한국어, Español, Français,
  Deutsch, Italiano, Português (Brasil) and Русский (figures stay in English); every panel and window can be
  shown, hidden or maximised, `H` keeps only the 3D view, `F11` full screen.
* **Local and extensible** — engine on your CPU cores, rendering on your GPU; JavaScript and Python plugins
  with a template generator and a **[developer guide](docs/PLUGINS.md)**.

**Download:** installers for Windows, macOS and Linux (.deb, .rpm) — Python, engine, ViennaRNA, demo data and help included — at
**[www.rnaconfviewer.com](https://www.rnaconfviewer.com/#download)**. From the source code:

```bash
pip install -e .            # or: export PYTHONPATH=src
rcv                         # opens the app (python -m rcv works too); rcv --demo for example data
```

User guide: **[docs/RCV.md](docs/RCV.md)** · plugins: **[docs/PLUGINS.md](docs/PLUGINS.md)** · installers: **[packaging/README.md](packaging/README.md)** ·
engine status and benchmark: **[docs/V2_STATUS.md](docs/V2_STATUS.md)**.

![RCV viewer](docs/images/rcv_viewer.png)

## Prediction engine (rnafold3d)

```
FASTA ──► ViennaRNA 2D hypotheses (MFE + distinct suboptimals; thermodynamic, Turner 2004)
      ──► rigid all-atom A-form stems (calibrated geometry, exact covalent chemistry)
      ──► fold-tree search in torsion space (annealing; closure-preserving CCD moves)
      ──► CCD loop closure + clash relief
      ──► independent all-atom validation (fail-closed)
      ──► ranked, structurally distinct ensemble: PDB + mmCIF + manifest.json
```

**Status (v2.0.1, 2026-10-22):** valid all-atom ensembles from sequence alone for 26/26 pilot
benchmark targets (14–192 nt); accurate for hairpins and small motifs (≈2–4 Å), near random-arrangement
accuracy for large multi-branched RNAs (no tertiary-contact physics yet).
Details, benchmark results and limitations: **[docs/V2_STATUS.md](docs/V2_STATUS.md)**, tables in
`docs/benchmark_v1/`, and
`RNAFold3D_Project_Progress_and_Next_Steps.md`. Earlier milestone notes (M0-M4, guide-site
engine) are in `docs/milestones/`; the original archives (`*.zip`) are kept unchanged.

## Install

Python >= 3.10 with `numpy`, `scipy`; `ViennaRNA` (python module `RNA`) for sequence-only
prediction (`conda install -c bioconda viennarna` or `pip install ViennaRNA`; an external dependency under its own licence: this source
repository does not contain it, the desktop installers bundle it with credit as its licence allows, see
`packaging/THIRD-PARTY-NOTICES.txt`). Optional: `threadpoolctl`.

```bash
pip install -e .            # or: export PYTHONPATH=src
```

## Use

```bash
# sequence only (2D predicted by ViennaRNA)
python -m rnafold3d predict --fasta my_rna.fa --models 5 --seed 42 --out runs/my_rna --workers 8

# with a known secondary structure (dot-bracket; pseudoknotted pairs are dropped with a warning)
python -m rnafold3d predict --fasta my_rna.fa --dbn my_rna.dbn --models 5 --out runs/my_rna_2d

# independent validation of any PDB/mmCIF
python -m rnafold3d validate runs/my_rna/conformations/conf_001.pdb

# package an external ensemble (e.g. NMR models) in the same format the viewer will read
python -m rnafold3d import model_*.pdb --out ensembles/my_nmr --name my_nmr

# 2D hypotheses only
python -m rnafold3d ss --fasta my_rna.fa

# RCV app with a run folder, an .rcv file, an NMR/multi-model file, or several models of one RNA
rcv runs/my_rna
rcv models/*.pdb --chain A                     # several files of one RNA = one ensemble (wildcards work on Windows too)
rcv convert runs/my_rna -o my_rna.rcv          # one compact file with every conformation
rcv info my_rna.rcv --methods                   # what it contains (+ methods paragraph)
rcv extract my_rna.rcv -o pdb/                  # back to PDB files
rcv batch my_rnas.fasta                         # predict every sequence (one .rcv each)
rcv html my_rna.rcv -o my_rna.html              # interactive single-file web page
rcv export-document runs/my_rna --out my_rna.json   # ensemble document without the GUI
```

Output package (`runs/my_rna/`):

```
input/sequence.fasta, input/secondary_structures.dbn
conformations/conf_001.pdb, conf_001.cif, ...      chain A, residues 1..N, heavy atoms
metadata/manifest.json   schema rnafold3d.ensemble/1 (schemas/ensemble_manifest.schema.json)
metadata/scores.csv      every trajectory: status, score terms, cluster
metadata/clusters.json   greedy C4'-RMSD clusters (sampling statistics, not populations)
metadata/validation.json per-trajectory validation and loop-closure reports
config.json, logs/run.log
```

Runs are deterministic for a given seed and package version (independent of worker count).

## Benchmark

```bash
python benchmark/build_dataset.py              # 26 experimental targets, 14-192 nt (benchmark/manifest_v1.json)
python benchmark/run_benchmark.py --run-id v2 --split test --modes independent,oracle
python benchmark/legacy_baseline.py            # legacy M1-M3 pipeline on the same targets
```

`independent` = FASTA only; `oracle` = 2D annotated from the reference structure (upper bound for
3D assembly; never to be reported as sequence-only accuracy).

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests     # engine + app (+ JS units if Node.js is installed)
python tools/viewer_smoketest.py                          # 42 GUI scenarios in headless Edge/Chrome, with screenshots
```

## Layout

```
src/rnafold3d/
  v2 engine: structure_io, geomutil, build, fold, loops, search, checks, ss, predict, metrics,
             annotate, ensemble_io, cli; data/geometry_v1.json (calibrated parameters)
  evaluate.py  model-independent conformation energy terms (used by the app)
src/rcv/
  app.py (entry point, CLI tools, `engine` passthrough), server.py (local API), jobs.py (predictions, batches,
  pull explorations), data.py (ensemble documents, 2D layouts), rcvfile.py (.rcv format), switch.py (2D switch
  analysis), methods.py, htmlexport.py, fonts.py (PDF font subsets), desktop.py (window shells, native
  dialogs), paths.py, templates/ (plugin templates),
  static/ (GUI: HTML/CSS/ES modules incl. console/, panels/, export/, viewer/; three.js r160 vendored under
  its MIT license; Liberation Sans under the SIL OFL 1.1; built-in plugins)
packaging/  PyInstaller spec, Windows (Inno Setup), macOS (.app/.dmg) and Linux (.deb/.rpm) build scripts, icons
  legacy M1-M4 guide-site modules (core, topology, assembly, closure, energy, sampling, atomistic, ...)
tools/      calibrate_geometry.py, refine_aform.py, loop_reach.py (parameter derivation);
            viewer_smoketest.py + viewer_harness/ (GUI test scenarios)
benchmark/  manifest, dataset builder, runner, legacy baseline, data/ (raw mmCIF git-ignored)
schemas/    ensemble manifest JSON schema (viewer contract)
tests/      legacy suite + test_v2_engine.py + test_app.py + test_pull_engine.py + test_rcv05.py + test_rcv06.py
            (+ js/app_unit.mjs, js/rcv05_unit.mjs, js/rcv06_unit.mjs)
docs/       RCV.md (user guide), PLUGINS.md (plugin developer guide), V2_STATUS.md, images/
examples/   viewer_demo/ (TAR NMR ensemble, TAR and tRNA predictions), viewer_plugins/ (Python plugin example),
            scripts/ (console script example)
```
