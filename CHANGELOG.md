# Changelog

## Package 0.7.0 — 2026-10-26 — RCV · RNA Conformation Viewer

From pictures of conformations to the mechanism between them: how A turns into B, how fast, through which
intermediates, which base pairs and which coordinated 3D motions limit the switch, and which mutations, salt or
ligand change it. The new analyses run on the user's computer in under a second and are validated on the
benchmark ensembles (`docs/V07_VALIDATION.md`).

### Added
- **Pathway tab: transition from conformation A to B** (`rcv/transition.py`, `js/panels/pathway.js`; user guide
  8.1).
  - *Basins:* A and B first relax into their local minima of the 2D landscape (steepest descent), and kinetics are
    measured between the basins. A and B in one basin are reported as fluctuations of one fold.
  - *Routes:* the direct path between the basins (ViennaRNA findpath), plus a network of intermediates: the basins
    of the ensemble's other 2D states and local minima from Boltzmann samples (fixed seed, reproducible).
  - *Kinetics:* transition path theory gives committors and the routes ranked by reactive flux. The **mean
    first-passage time** gives the switching time, the half-life of A and the effective barrier; k₀ is adjustable.
  - *Traps:* wells deeper than A that hold the RNA for most of the switching time, found from the expected time
    spent in each node.
  - *Display:* a plain-language verdict, the routes on the ensemble map, and the energy profile (relaxation, route,
    direct path, saddle) with hover showing the base pairs that open or form.
- **Bottleneck and coordinated motions:** the **cooperative unit** (base pairs that open together on the climb to
  the highest saddle) and the **nucleus** of B. In 3D, **rigid modules** (helix-seeded, merged within 2 Å), with
  their turn, shift, screw axis and time window, plus **hinges**. The module that carries the unit is named the
  **rate-limiting coordinated motion**. Colour by mechanism, event order or module.
- **Mechanism morph** (`js/mechmorph.js`): A turns into B with the rigid modules moving in the order of their
  base-pair events. Each module keeps its shape, and a distance-geometry refinement keeps the chain connected.
  Chain distortion is 0.20 vs 0.81 for linear interpolation. *Mechanism video…* records it.
- **φ-value scan:** all 3 × length point mutations along the route, with ΔΔG‡, ΔΔG and φ. Mutations that speed up
  or slow down the switch, colour by kinetic sensitivity, copy as a table.
- **Topology-aware ensemble map** (`rcv/landscape.py`, `js/landscape.js`; Energy tab ▸ *Map*, also in the Pathway
  tab, the plot studio and the HTML page):
  - One disk per secondary-structure topology (abstract shape), labelled with a small 2D drawing; inside each disk,
    distances follow C4′ RMSD and 2D states stay together.
  - The caption gives trustworthiness, continuity and topology purity, compared with MDS.
- **Conditions** (`rcv/conditions.py`, `js/conditions.js`; Tools menu): temperature, Na⁺ (ViennaRNA salt
  correction) and a ligand that binds an aptamer motif (theophylline preset, or your own motif and Kd).
  - Kept per ensemble and saved in `.rcv` files.
  - Used by the Switch and Pathway tabs; the Switch tab adds a *Ligand response* section (2D-state populations
    over ligand concentration).
- **Probing data** (`rcv/probing.py`, `js/probing.js`; Tools ▸ *Load probing data*):
  - Reads SHAPE and DMS reactivities (`.shape`, ShapeMapper `.map`, CSV, plain lists) with 2–8 % normalisation.
  - Adds the colour scheme *Reactivity* and the metric *SHAPE agreement (AUC)*.
  - **Maximum-entropy reweighting** of the ensemble (metric *SHAPE data weight*, χ², effective number of
    conformations); optional pseudo-energies for the 2D-state populations.
  - Saved in `.rcv` files.
- **Optional confidence model** (`rnafold3d/qa.py`, weights `data/qa_v1.npz`, training in `tools/qa_model/`):
  - Predicts the per-nucleotide lDDT of every conformation in the background: a small neural network with 3,217
    parameters, numpy only, no GPU, never part of structure prediction.
  - Adds the metric *Confidence (predicted lDDT)* and the colour scheme *Confidence (per nucleotide)*; switch it
    off in Preferences.
  - Per-nucleotide r = 0.75–0.84 on benchmark families it never saw. The engine's ranking is kept, because the
    model does not reliably pick a better top model.
- **Publication output:**
  - The plot studio adds *Ensemble map*, *Refolding path* and *φ-value scan*.
  - The interactive HTML page shows the map, the pathway and the mechanism morph.
  - Methods text covers the new analyses with exact settings and a list of **references** to cite.
- Console: `pathway A, B`, `phiscan`, `mechanism [staged|rigid]`, `conditions …`, `probing …`. Command line:
  `rcv pathway`, `rcv probing`, `rcv confidence`; conditions options `--temperature --sodium --ligand --kd --conc`.
- Local API: `/api/transition`, `/api/transition/scan`, `/api/landscape`, `/api/probing/parse`,
  `/api/probing/fit`, `/api/confidence`, `/api/switch/dose`; `/api/switch` takes conditions. Plugin events
  `pathway`, `conditions`, `probing` (`docs/PLUGINS.md`).
- All new interface text in the 10 interface languages.
- Validation: `docs/V07_VALIDATION.md`, `tools/validate07.py` (map vs MDS, basins and kinetics, φ hot spots vs the
  cooperative unit, morph strain, probing recovery, confidence model on the benchmark).
- Tests: `tests/test_rcv07.py` (MFPT and time occupancy against hand-solved and simulated networks, basins, traps,
  φ-values, rigid modules, conditions, probing, confidence), `tests/js/rcv07_unit.mjs`, and new GUI scenarios
  (`pathway1`, `pathwaybasin`, `probing1`, `conditions1`, `confidence1`, `htmlpage07`, `i18n07`, `plots07`, …).

### Changed
- The Energy tab's map views are *Funnel*, *Map* (topology-aware) and *MDS*.
- `js/ui.js` can be imported outside a browser (the console commands' unit tests load it in Node).

## Package 0.6.1 — 2026-10-25 — RCV · RNA Conformation Viewer

### Fixed
- **Distance tracker: tracked distances are saved.** They are stored in `.rcv` files (and in workspaces, for the
  shown ensemble) and come back when the file is opened, also the distance shown in the energy panel; before, they
  were lost on closing the ensemble or reloading the window. They are also kept when a pull exploration extends the
  ensemble.
- **Distance tracker: removing a distance (✕) removes its metric** from the energy panel, the Details table, the
  plot studio and the CSV export; if the energy panel showed it, it returns to the default metric.
- Distance tracker: adding and removing distances can be undone and redone (`Ctrl+Z` / `Ctrl+Y`), like measurements.
- A tracked distance now belongs only to the ensemble it was made on; it no longer adds an empty metric to the
  other open ensembles.

### Added
- **Plugin API version 3** (`js/plugins.js`, `docs/PLUGINS.md` sections 2.7 and 2.9; versions 1 and 2 unchanged):
  `api.metrics.register({ ..., ensemble })` for metrics of one ensemble, `api.metrics.unregister(key[, ensemble])`,
  and `api.session.register({ save, load })` / `api.session.changed()` for plugin state that is saved in `.rcv`
  files and workspaces (like measurements) and is part of undo/redo.
- GUI smoke test `trackersave` (remove, undo/redo, `.rcv` round trip, extended ensemble).

## Package 0.6.0 — 2026-10-24 — RCV · RNA Conformation Viewer

### Fixed
- **"Try the demo ensembles" failed with `501 Unsupported method ('{}GET')`**: the local server now reads every
  request body (also of requests that do not use it), so the next request on a kept-alive connection is not
  corrupted. Regression test in `tests/test_rcv05.py`.
- The old logo in the window's title bar and the taskbar: the app window now has its own browser profile
  (`--user-data-dir`, `--class=RCV`, no extensions), a multi-size `favicon.ico` (16–256 px), PNG icons up to
  512 px and a web-app manifest, so the RCV logo is shown sharply on high-DPI screens.
- Dialog widths: a style bug kept every dialog at 560 px; wide dialogs (plot studio, palettes, figure panels)
  now get their width.
- Figure panels: the size of the 2D-structure dots in the exported file now equals the preview (drawing unit
  scaled with the export resolution); tall 2D layouts are turned to fit their box.
- **Plot studio: the pairwise RMSD matrix did not draw** ("L is not a function"): a local variable hid the
  line-drawing helper used by the colour-bar ticks.
- **Installed app could not load ViennaRNA on some PCs** ("DLL initialization routine failed"): PyInstaller had
  bundled an old `MSVCP140.dll` found on the build machine's PATH (from a Java JDK). The build now takes each
  Visual C++ runtime DLL from Windows when that copy is newer.
- **Installed app on Windows in Chinese and other languages**: frozen builds ignore `PYTHONIOENCODING`, so the
  engine wrote its progress in the system code page (e.g. GBK) and failed on characters such as `Å`. The app now
  runs in Python UTF-8 mode with unbuffered output: progress arrives line by line and non-ASCII names work.

### Changed
- **Start page simplified**: one box with *Predict conformations*, *Open files…*, *Fetch PDB…* and *Demo
  ensembles*; name, conformations, sampling, cores, dot-bracket, seed and links are under *Advanced settings*.
  The words "(no deep learning)" were removed from the start page.
- **Figure font**: all figure text (images, panels, plots, videos, PDF) uses **Liberation Sans**, bundled under the
  SIL Open Font License 1.1 (metric-compatible with the Helvetica/Arial style journals ask for; no proprietary
  font is used or required). PDFs embed subsets of the font (`rcv/fonts.py`, fontTools; standard PDF fonts
  without fontTools); PDF text supports rotation (vertical axis titles).
- **Right-panel charts and tables at publication level**: axis lines with outward ticks, rotated axis titles with
  units, figure font, an export button per chart (PDF/SVG/PNG/TIFF/copy) and three-rule (booktabs) tables.
- Representation is a drop-down in the toolbar; labels hide on narrower windows so the toolbar stays on one row.
- About dialog, citation, installers and documents name the author: **Zixiang Wang (princewang2018@gmail.com)**.

### Added
- **Desktop behaviour** (`js/nativefeel.js`): RCV's own right-click menus everywhere (residue menu, chart export,
  table copy, cut/copy/paste in fields) instead of the browser menu; no text selection outside text fields; no
  image dragging, middle-click scrolling, page zoom, reload, print, find, view-source or back/forward — these keys
  do RCV things (`Ctrl+F` select, `Ctrl+P` export image, `Ctrl` + wheel / `Ctrl+±` zoom the 3D view). *Developer
  mode* (Preferences) restores the browser tools for plugin authors.
- **Mouse tools** (`js/tools3d.js`): Select `V`, Pull `G`, **Paint** `B` (residues or atoms, palette or any
  colour), **Label** `T`, Measure `M` — in the toolbar and the View menu.
- **Annotations** (`js/annotations.js`): residue labels, region brackets over a selection and free notes/titles;
  drag to move, double-click to edit text, colour, size, arrow and box; follow morphs; drawn into images, panels
  and videos; saved in `.rcv` files.
- **Undo / redo / clear** (`js/history.js`): toolbar buttons, Edit menu, `Ctrl+Z` / `Ctrl+Y`; per-ensemble
  history of annotations, paint, colours, display, labels, measurements, selection, scheme, representation, style,
  surface and transparency (console commands included); *Clear* by kind or everything (undoable).
- **Tools menu** (`js/tools.js`) with the command-line tools: convert to `.rcv`, inspect `.rcv`, extract PDB (ZIP or
  multi-model), batch from FASTA, **Fetch from the PDB** dialog (`Ctrl+Shift+F`, also on the start page), HTML
  page, figure panels, plot studio. Server endpoints `/api/rcv/info`, `/api/rcv/extract`.
- **Plot studio** (`js/plots.js`, `js/cluster.js`, `js/export/chartexport.js`): pairwise RMSD matrix (UPGMA
  clustering with dendrogram, colour bar, values), RMSD to a reference (first, best, shown, previous, mean),
  per-residue RMSF with pairing track and displacement, energy-vs-RMSD funnel coloured by 2D state, energy
  profile, vector 2D structure; journal sizes in mm, font size in pt, PDF/SVG/PNG/TIFF at 300–1200 dpi, copy.
- **Colour palettes** (`js/palettes.js`, `js/palettecheck.js`): RCV, Okabe–Ito, Tol bright, Tol muted, Pastel,
  Grayscale and user palettes (edit, duplicate, JSON copy/paste); the palette drives all colour schemes, paint,
  annotations and plots; colour-vision check (deuteranopia, protanopia, tritanopia; closest pair in OKLab ΔE).
- **Molecular surface** (`js/viewer/surface.js`, `js/surfacectl.js`): Gaussian surface of all residues or the
  selection, coloured like the residues, rebuilt during morphs and video frames (`U`, toolbar, `show surface`).
- **Per-object transparency**: cartoon, bases, atoms, surface and overlay (View ▸ Surface and transparency…,
  `set transparency, 0.4, surface`).
- **Plugin development**: guide `docs/PLUGINS.md` (also shown inside RCV), *New plugin from template…* (JavaScript
  or Python, `/api/plugins/new`), *Open plugin folder*; plugin API version 2 adds `commands`, `annotations`,
  `palette`, `display.surface/transparency/color` and `export.chart`.
- **Interface languages** (`js/i18n.js`, `static/i18n/*.json`): English, 简体中文, 繁體中文, 日本語, 한국어,
  Español, Français, Deutsch, Italiano, Português (Brasil) and Русский — about 1,100 interface strings each.
  Globe button in the menu bar and on the start page, and in Preferences; switches at once and is remembered;
  follows the system language until chosen. Figures, charts, exports, the console and the methods text stay in
  English on purpose.
- **Window layout** (`js/windows.js`, **Window** menu): show or hide every part (side panels, 2D window, console,
  sequence bar, toolbar, playback bar, status bar, 3D title and legend); maximise the side panels, the 2D window
  and the console (title-bar button or double-click; `Esc` restores); `H` only the 3D view; `F11` full screen;
  reset layout. The layout is remembered.
- **Help inside RCV**: Help ▸ *User guide* (`F1`, with contents and pictures) and *What is new in this version*
  show the documentation in the app, also in installed copies.
- **Installers**: Windows setup in nine languages with Start-menu entries for RCV, *RCV with demo ensembles*, the
  user guide, *What is new* and the plugin guide (as web pages, `packaging/make_docs.py`), licence and
  `THIRD-PARTY-NOTICES.txt` (ViennaRNA credit and licence, and the other bundled software); macOS disk image
  with an *Applications* link and the user guide (macOS 14 or later); **Linux packages** `.deb` (Debian, Ubuntu,
  Mint) and `.rpm` (Fedora, RHEL, Rocky, AlmaLinux, openSUSE) with the `rcv` command, a menu entry (*RCV with demo
  ensembles*, *RCV User Guide*), icons, the `.rcv` file type and the user guide (`packaging/linux/`). The GitHub
  workflow `.github/workflows/installers.yml` builds the Windows x64, macOS Apple Silicon, macOS Intel and Linux
  x86_64 installers, test-installs the Linux packages on Ubuntu, Debian, Fedora and AlmaLinux (and, for a version
  tag, attaches everything to the release behind the download links of www.rnaconfviewer.com).
- Tests: `tests/test_rcv06.py`, `tests/test_rcv06_lang.py`, `tests/js/rcv06_unit.mjs`, `tests/js/i18n_unit.mjs`
  (translation files and matcher, help documents); GUI smoke test extended to 42 scenarios (window layout, all
  languages, user guide).

## Package 0.5.0 — 2026-10-23 — RCV · RNA Conformation Viewer

### Changed
- Renamed to **RCV (RNA Conformation Viewer)** with a new logo (hairpin in a viewfinder with a second,
  translucent pose). Python package `rcv`, command `rcv`, distribution `rna-conformation-viewer`, installers
  `RCV.exe` / `RCV.app`, user folders `RCV`. Settings, plugins, projects, sessions and documents of
  RNAconfSwitch 0.4 are still read.
- Document schema `rcv.document/1` (older schemas accepted); documents carry 2D layouts (ViennaRNA NAView) and
  source-file checksums.
- `Ctrl+S` saves an `.rcv` file; the JSON session became the *workspace* (`Ctrl+Shift+S`).
- Toolbar labels shortened so the toolbar fits on one row.

### Added
- **`.rcv` file format** (`rcv/1`, ZIP + JSON, `rcvfile.py`): all conformations of one RNA with energies, terms,
  2D structures and layouts, validation, ensemble statistics, provenance, engine manifest, methods text,
  thumbnail and the viewer state; coordinates as delta-coded, byte-shuffled 0.001 Å integers with SHA-256
  check — 10–14× smaller than the same models as PDB. Predictions write `<name>.rcv` next to their PDB files
  and reopen from it instantly. CLI: `rcv convert`, `rcv info`, `rcv extract`. Installers register `.rcv`.
- **2D structure window** (`panels/ss2d.js`, `ss2dlayout.js`): corner window with the 2D drawing of the shown
  conformation, coloured like 3D, linked hover/selection, morphing in step with 3D (pairs fade out/in), new
  pairs highlighted, pseudoknot arcs, three sizes, four corners; fallback radial layout without ViennaRNA.
- **Transition video A → B**: direct or through intermediates (lowest-RMSD route), holds, boomerang, slow
  turn, overlay with A/B energies and the metric along the route, morphing 2D inset; right-click entry in the
  conformation list.
- **Start page**: *Predict* with an empty box runs the demo sequence (HIV-1 TAR); **FASTA batch prediction**
  (`BatchJob`, `rcv batch`): per-sequence validation, one stage per sequence, overall time estimate, results
  table with a required choice of the prediction to open, `batch_summary.csv`; *Load FASTA…* and drop.
- **Command console** (`console/`): RCV's own command language with ~45 commands (select, show/hide per
  residue, color, bg, label, zoom/center/orient/turn/spin, conf/morph/play, distance/angle/dihedral,
  align/rms/fit, pairs, info, view save/recall, ls/delete, set, png/save/video, load, fetch, run scripts),
  an RNA-ensemble selection language (resi/pos/resn/chain, paired/unpaired/pk, helix, switching, changed,
  rmsf/moved thresholds, within … of, named selections), history, tab completion, multi-line paste; the
  command log is stored with each ensemble.
- **Sequence bar**, **measurements** (measure mode `M`, live values during morphs, plotted as metrics),
  **named selections** and **scenes** (camera fly-to), per-residue representations, custom colours and labels.
- **Switch tab** (`switch.py`, `panels/switch.js`): 2D states of the ensemble with ViennaRNA energies,
  full-ensemble probabilities and Boltzmann shares at any temperature; single-mutation scan between two states
  (ΔΔG heat map, two-state shares, strongest shifts with full-ensemble probabilities), *Predict in 3D* for a
  mutant, *Colour 3D by mutation sensitivity*, copy as table.
- **Figure panels** (`export/panels.js`): several conformations in one figure, lettered, with captions, 2D
  drawings beside each model and a shared legend; journal widths, PNG/TIFF/PDF.
- **Interactive HTML export** (`htmlexport.py`, `js/lite.js`, `rcv html`): single self-contained page (viewer,
  morphs, energies, 2D, methods) that opens from a file in any current browser.
- **Methods text & citation** (`methods.py`): methods paragraph from the recorded provenance, citation, command
  log; stored in `.rcv` files.
- **Right-click residue menu** in the 3D view (select residue/pair/helix/neighbours, sticks, spheres, hide, label,
  colour, zoom, measure from here), every action logged as a console command.
- Colour schemes **Structure elements**, **Pairing frequency**, **Switching residues**, **Mutation sensitivity**.
- Tests: `tests/test_rcv05.py`, `tests/js/rcv05_unit.mjs`; GUI smoke test extended to 26 scenarios; example console script `examples/scripts/tar_switch_figure.rcvs`.

## Package 0.4.0 — 2026-10-23 — RNAconfSwitch (engine 2.0.1 + pull restraints)

### Added
- **Renamed to RNAconfSwitch**; new application package `rnaconfswitch` (command `rnaconfswitch`,
  `python -m rnaconfswitch`); the engine keeps its module name `rnafold3d` and its `rnafold3d` command.
- **Local desktop application**: native window through pywebview (WebView2/WKWebView) or an Edge/Chrome app
  window; closing the window quits; native Open/Save dialogs; `rnaconfswitch engine ...` runs the bundled engine
  (used by frozen builds); user data, logs, plugins and a projects folder per operating system.
- **Start page**: predict from sequence (validation, sampling depth, CPU cores, time estimate) or open
  structures; recent predictions; demo button; system line (CPU cores, GPU, ViennaRNA).
- **Prediction jobs with progress**: stages, progress bar, time left, plain-language notes, cancel; the
  viewer opens automatically when done.
- **Pull to switch**: elastic drag of a residue, live preview of the nearest stable conformation (view-aware
  matching on the unmoved residues, stability preference among near-equal matches), anchored snap on release,
  result card with Undo and **Explore new conformations**.
- **Engine pull restraints**: flat-bottom distance restraints as a removable bias (`ScoreWeights.pull`,
  `Scorer(pull_restraints=...)`), unbiased relaxation stage, seed secondary structures, `--job-json`;
  explorations are merged into the open ensemble ("pull" badge) with unbiased scores.
- **Packaging**: PyInstaller spec, Windows build script + Inno Setup installer, macOS `.app`/`.dmg` script,
  icons; verified on Linux including spawn-started engine workers in the frozen app.
- Jobs dock, Preferences ▸ projects folder, `docs/RNAconfSwitch.md`, `packaging/README.md`,
  `tests/test_app.py`, `tests/test_pull_engine.py`, pull tests in `tests/js/app_unit.mjs`, smoke-test scenarios
  welcome / pull / predict / explore.

### Changed
- Ensemble documents use schema `rnaconfswitch.document/1` (older `rnafold3d.viewer/1` documents and
  `.rnaview.json` sessions still open); sessions are saved as `.rnaconf.json`.
- Residues are matched by position when all conformations share the same sequence (lets NMR models numbered
  17–45 merge with engine models numbered 1–29).
- Engine worker processes use `spawn` outside Linux (macOS app bundles); `RNAFOLD3D_MP_START` overrides.
- The project library scans the projects folder, demo data and `./runs` only (fast start page).

## Package 0.3.0 — 2026-10-22 — RNAconformer Viewer (engine unchanged at 2.0.1)

### Added
- **RNAconformer Viewer** (`python -m rnafold3d view`): local desktop app (stdlib server + Edge/Chrome app
  window, works on WSL) for several conformations of one RNA. Trackball 3D view (cartoon tube + base slabs,
  sticks, ball & stick, spheres), styles with screen-space ambient occlusion and outlines (Studio,
  Illustrative, Glossy, Depth cue, Plain), validated colour schemes (nucleotide, helix families, pairing,
  RMSF flexibility, displacement, position), ensemble overlay, flexibility putty, residue selection and labels.
- Conformation switching with eased morphs, timeline scrubbing and tours; orders by rank, energy,
  similarity path or input; hover previews; superposition on selected residues.
- Energy panel that follows the conformation: value and change vs the previous conformation, energy
  profile with moving marker, funnel plot and MDS similarity map, additive term breakdown, two-sided
  base-pair arc diagram.
- Publication output: PNG (DPI), TIFF (Deflate, DPI, alpha), PDF (vector energy panel and text, Symbol
  font for Greek), journal size presets, transparent backgrounds; MP4/WebM/PNG-frame videos with energy
  overlay; PDB (current/all), CSV energy table, sessions.
- Plugins: JavaScript interface plugins (panels, menu items, metrics, colour schemes, 3D overlays) and
  Python analysis plugins; built-in Distance tracker and Similarity matrix; example Python plugin.
- In-GUI prediction from sequence with progress, run library and server-side folder browser.
- `rnafold3d.evaluate`: model-independent energy terms (ΔG 2D via ViennaRNA, Debye–Hückel electrostatics,
  loop/tertiary stacking, clash penalty) and descriptors for any RNA model.
- `structure_io.read_models`: every model of PDB/mmCIF files (path or text).
- `rnafold3d export-view`, demo data (`examples/viewer_demo`), `tests/test_viewer.py`,
  `tests/js/viewer_unit.mjs`, `tools/viewer_smoketest.py`, `docs/VIEWER.md`.

### Changed
- Engine logs per-trajectory progress (`progress k/N`) as trajectories finish; results are re-sorted by
  trajectory index, so outputs are unchanged.
- Package version 0.3.0; GUI files are shipped as package data.

## 2.0.1 (2026-10-22)
- Stems: if an ideal A-form stem contains a compressed pair (any Watson-Crick/wobble H-bond < 2.4 Å) or a
  non-bonded contact < 2.5 Å (seen next to G·U wobbles and at stem ends), strand-2 torsions are refined
  by regularised least squares toward the calibrated pair geometry (~12 % of stems in the pilot set).
  Found by failure analysis of the v2.0 pilot run (3MXH: 0/8 valid models in both modes).
- Parallel runs: all BLAS/LAPACK libraries (including SciPy's bundled OpenBLAS, loaded lazily by
  `scipy.optimize`) are imported before being limited to one thread per worker, and
  `OPENBLAS/OMP/MKL_NUM_THREADS=1` are set for CLI/benchmark processes. Without this, the stem
  refinement oversubscribed the CPU (load ~30 on 8 cores; one target took > 5 h instead of ~2 min).
- User 2D input: non-canonical pairs are dropped with a warning (they produced invalid models).
- CLI: input errors are reported without tracebacks (exit code 2).
- Removed dead code; pseudoknot restraints documented as experimental and off by default.

## 2.0.0 (2026-10-21)
- New v2 engine: sequence-only input via ViennaRNA, calibrated all-atom builder, fold-tree torsion
  search with closure-preserving CCD moves, staged loop closure, independent validator, ensemble
  package (manifest schema `rnafold3d.ensemble/1`), `predict/validate/ss/import` CLI.
- Benchmark dataset v1 (26 targets), runner, legacy/null baselines.
- Legacy M0–M4 modules unchanged.
