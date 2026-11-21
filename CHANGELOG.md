# Changelog

## Package 0.7.6 — 2026-11-21 — RCV · RNA Conformation Viewer (engine rnafold3d 2.5.0)

Modified nucleotides in the prediction from sequence: 2D and 3D structures of RNA with m6A, pseudouridine,
A-to-I editing, oxidative damage, 2′-O-methyl, LNA, phosphorothioate and more (`docs/RCV.md` §3.3). Z-RNA:
left-handed helices marked in every displayed frame, and optionally predicted (§3.4). An audit of the basic
functions (2026-11-20): cancelling stops the whole computation, more PDB entries open, ligands, ions and
partner chains are shown next to the RNA, the 2D structure can be exported as dot-bracket, CT or BPSEQ, and
errors are reported in words.

### Added
- **Catalogue of 97 modified nucleotides** (`src/rnafold3d/modifications.py`): methylations (m6A, m1A, m5C,
  m7G, m2,2G …), pseudouridine and m1Ψ, dihydrouridine, A-to-I (inosine) and C-to-U editing, damage (8-oxo-G,
  8-oxo-A, O6-methyl-G), thio and halogen bases, 27 hypermodified tRNA nucleosides (t6A, i6A, queuosine,
  wybutosine, mnm5s2U …), sugar and backbone modifications (Nm, 2′-F, LNA, deoxy, phosphorothioate) and
  synthetic analogues (2-aminopurine, 2,6-diaminopurine, 5-bromo-U …). Each has its MODOMICS short name,
  aliases, PDB chemical component, parent base, effect on pairing, the atoms to build and a note.
- **Sequence notation**: `GGC[m6A]CU[Y]A*G` (brackets for a modified nucleotide, `*` for a phosphorothioate),
  in the start page, FASTA files, batches and on the command line; `rnafold3d predict --modifications
  POS:CODE,...`; `rnafold3d modifications [--json]` lists the catalogue.
- **2D with modifications**: ViennaRNA's measured parameters (m6A, Ψ, inosine, dihydrouridine, purine,
  7-deaza-A); generated soft constraints with estimated corrections for the others, including new pairs that
  only the modification allows (8-oxo-G·A, m2,2G·A, k2C·A); hard constraints keeping Watson-Crick-blocking
  modifications unpaired. They work together with probing data and alignments, and in `evaluate`.
- **3D with modifications, inside the search** (`fold.py`, `loops.py`, `search.py`):
  - sugar pucker preferences in the loop conformers drawn for non-helical residues (D and dN → C2′-endo; Nm,
    2′-F, s2U, Ψ → C3′-endo; LNA locked C3′-endo); helical residues keep A-form;
  - syn/anti: a chi prior with a syn well for 8-oxo-G, 8-oxo-A and m8A (1.5–3 kT) and an anti bias for LNA and
    2′-F, used by the annealing, the loop-closure torsion energy and the chi refinement;
  - modification atoms as *riders* of the model, placed from their base, sugar or phosphate frame before every
    contact evaluation and scored with the other atoms; removed atoms (O2′ of dN, N2 of I, O2 of s2U) left out;
  - interactions of the modification itself (`modinteract.py`, score term `modification`): hydrogen bonds of
    the polar groups it adds and of the donors it creates (Ψ N1-H, 8-oxo N7-H), with parent-atom roles updated
    for methylation; the Ψ N1-H water bridge to its own or the preceding phosphate (not for m1Ψ); a stacking
    change for residues outside helices (Ψ, m1Ψ, m5C, m6A, i6A/t6A, wybutosine, 5-halo-U stronger; D weaker);
  - a local relaxation of the modified residues and their neighbourhood at full steric weight;
  - the finished models get the final atoms (rotamers against the whole model); the manifest reports per
    modified residue its pucker and syn/anti against the preference, and the contacts before and after relaxation.
  Pairs a modification allows become helix restraints. PDB output writes HETATM and MODRES records, mmCIF the
  one-letter sequence with (CCD) and a `_chem_comp` loop with the parent. Unmodified sequences sample exactly as
  before.
- **Start page ▸ Ψ Modified nucleotides…**: an editor window — clickable sequence map, searchable list per base
  with pairing badges and notes, phosphorothioate linkages, presets (m6A at DRACH motifs, every U → m1Ψ / Ψ)
  and a summary table. The status line counts the modifications.
- **Viewer**: *Colour by ▸ Modified nucleotides* with a legend by kind, the modification in the hover tip,
  `select modified` and `resn m6A` in the console, a paragraph in the methods text; modified residues of opened
  PDB/mmCIF files are recognised by their PDB component.
- **AlphaFold 3 (local)** passes modifications with a PDB component in its input JSON.
- `GET /api/modifications`, `POST /api/modifications/parse`; `tests/test_v076.py`.

- **Z-RNA identification** (`src/rnafold3d/zrna.py`, `src/rcv/static/js/zrna.js`), on by default: regions of
  at least three stacked Watson-Crick pairs with left-handed two-step twist and syn purines / anti pyrimidines,
  found on every displayed frame (conformations, morphs, scrubbing, video, pulling) and drawn with their own
  theme — halo, dashed axis and badge in 3D, halo and tag in 2D, also in image, panel and video exports. *View ▸
  Z-RNA regions* switches it; `select zrna` in the console; `zrna_regions` per model in the manifest. Finds the
  Z duplexes of PDB 2GXB and 1T4X (12/12 models) and nothing in the 3964 chains of the template library.
- **Z-RNA prediction** (*Advanced settings ▸ Predict Z-RNA (left-handed helices) when likely*, off by default;
  `rnafold3d predict --zrna auto --zrna-salt M --zrna-binder`): an estimated A→Z free energy per helix segment
  (dinucleotide steps, salt, junctions, 8-substituted purines, Z-binding protein); below +3 kcal/mol a hypothesis
  with that segment built as a Z-form stem (`build_stem(..., form='Z')`, torsions fitted to 1T4X and 2GXB, 1.1–1.4 Å
  from them) and the junction pairs opened, with a quarter more trajectories. The estimate and the candidates
  are in the manifest (`input.zrna`) and the methods text.
- **m8G and 8-Br-G** in the modification catalogue (syn guanosines; Z-RNA promoters).
- `tests/test_zrna.py`, `tests/js/zrna_unit.mjs` (Python and JavaScript detectors give the same regions),
  `tests/data/zrna_2gxb_duplex.pdb`; smoke scenarios `zrna`, `zrnashot`.

- **Microsoft Store package (MSIX)** for Windows: `packaging/windows/build_msix.ps1` packs the PyInstaller
  build with a manifest template (`packaging/windows/msix/AppxManifest.xml.in`: full-trust desktop app, `.rcv`
  file association, identity from Partner Center), tile and Store icons at every scale (`msix/Assets/`, drawn by
  `make_assets.py`) and a `resources.pri`; the *Build installers* workflow builds it as the artifact
  `store-msix-x64` (identity from the repository variables `MSIX_IDENTITY_NAME`, `MSIX_PUBLISHER`,
  `MSIX_PUBLISHER_DISPLAY_NAME`). Submission steps: `packaging/README.md`.
- `RCV.exe` carries a Windows version resource (publisher, product, description, version).

- **What is bound to the RNA is shown** (audit C5; `rcv/data.py` `context_of`, `static/js/viewer/context.js`,
  `contextctl.js`). Opening a PDB entry kept the RNA chain and dropped everything else, so a riboswitch appeared
  without its ligand and an RNA–protein complex as bare RNA. Now ligands are drawn as sticks, ions as spheres
  and protein, DNA and other RNA chains as thin traces (one point per residue), named on hover. They move with
  their conformation (superposition, morphs, scrubbing, mechanism playback, figure panels, videos), are stored
  in `.rcv` files and shown on the interactive HTML page. One switch: *View ▸ Bound molecules*. They are context
  only: no analysis uses them. In an ensemble only what every conformation has is shown; hydrogens and water are
  left out; above 6000 atoms or 30 000 residues the parts nearest to the RNA are kept. AlphaFold 3 results with
  partner chains show them the same way. Smoke scenarios `bound`, `boundshot`.
- **Export of the secondary structure** (requested): *Export ▸ Secondary structure (dot-bracket, CT, BPSEQ)…*,
  the download button of the 2D window, and `rcv ss INPUT [-o FILE] [--format dbn|ct|bpseq] [--which
  distinct|all] [--conformation N]`. Dot-bracket is the ViennaRNA form (`>name`, sequence, structure) read by
  RNAfold, RNAplot, forna, VARNA and RNAstructure; CT is the connect table of mfold and RNAstructure (the last
  column keeps the residue numbers of the file); BPSEQ. The shown conformation, every different 2D structure once
  (with the conformations that share it) or all; preview, copy and save. Window and command line write the
  same text (`rcv/ss2dexport.py`, `static/js/export/secondary.js`; compared in the tests).
- `tests/test_v076_audit.py` (51 tests), `tests/js/audit_unit.mjs`, `tests/js/notes_unit.mjs`; smoke scenarios
  `audit76`, `langstart` and `langnotes`.
### Changed
- **Sequences copied from a database record can be pasted as they are** (audit C1): position numbers, spaces
  and a closing `//` are ignored (`1 ggcagatctg agcctgggag …`), in the start page, FASTA files, batches and
  `rnafold3d predict --sequence`. Digits inside the brackets of a modified nucleotide (`[m2,2G]`) are kept.
- **A known 2D structure is checked before the run** (audit C2): length, unmatched and unclosed brackets with
  their position, characters that are not dot-bracket, and hairpin loops shorter than 3 nt (which made every
  trajectory fail after minutes). The same check in the start page, the batch list, the server and the engine.
- **Cancel** (audit C3): a run older than a minute asks first, because what was computed is discarded and a
  cancelled run cannot be continued; the button shows *Cancelling…* until the engine has stopped.
- **Notes of a running prediction can be read** (audit C4): the list follows the newest note only while it is
  scrolled to the bottom; it is refilled only when a note arrives, so a selection in it survives.
- **Start page wording** (audit C6): the empty box no longer says twice that the demo sequence will be used.
- **The window says when something is missing or breaks** (audit B1–B3, B5): a lost connection to the engine is
  stated in the progress card, which keeps the last known state and reconnects (before: an empty card after one
  failed request); without WebGL 2 a page explains the usual reasons instead of a blank window; a lost graphics
  context is announced and the view redrawn when it is back; errors that no code handles are shown in a note and
  written to the log. In browser-window mode RCV quits ten minutes after its window is gone, unless a
  computation is running; a computer that slept does not count.
- **Local server** (audit B6): the session token is checked before the body of a request is read; a body that
  is not a JSON object, and numbers given as words, are answered with status 400 and a sentence naming the value.
- **Messages for files that cannot be read** (audit B4): a truncated ATOM line, a coordinate that is not a
  number and an incomplete mmCIF atom table are reported with the file, the line or atom and what is wrong; a
  file without RNA says what it holds (“protein (129 residues) and 2 ligands or ions (MG, SAM)”). The engine's
  command line reports impossible values (`--models 0`, `--seed -1`, `--modifications abc:m6A`) in one line.
- Mirror-image RNA (L-nucleotides) opens with a note that the analyses are made for natural D-RNA, and is left
  out when the template library is rebuilt.
- CASP16 result files, figure and tables regenerated from the 0.7.6 runs (`benchmark/casp16/results/`).
- **Engine 2.5.0: the fit onto the template is repeated instead of given up** (`predict.py`,
  `PredictConfig.template_fit_attempts`, default 4). The last stage of a template model moves every atom onto
  the template. When the result failed the geometry checks, the fit was discarded and the model kept the looser
  shape of the torsion-space search — in about half of the template models. Now the residues with the problems
  and their neighbours (one more on each side per attempt) are released from the template and the fit is
  repeated, up to 35% of the chain; models whose fit is still rejected rank after those with an accepted fit.
  Measured with the step on and off, same seeds:
  - 48 library chains with a template, released after their templates (`benchmark/template_models.py`):
    rejected fits 73 of 150 → 8 of 138; mean TM-score of model 1 0.582 → 0.605 (+0.023, 95% interval +0.007 to
    +0.043; 8 chains better by more than 0.02, 2 worse), best of 5 0.613 → 0.622 (+0.009, −0.001 to +0.023).
  - CASP16, 15 targets (`runs_v076` against `runs_v076_fit1`, which reproduces 0.7.5 exactly): model 1
    0.495 → 0.501 (+0.006, −0.000 to +0.017; R1271 0.39 → 0.47, the others unchanged), best of 5 0.527 → 0.527;
    right global fold on 9 of 15 targets (8), median RMSD of model 1 7.8 → 6.0 Å; 68 minutes for the 15
    targets with and without the step.
  - The six chains with a template of the 26-chain development set: model 1 0.508 → 0.497 (two chains −0.02 to
    −0.03; within the noise of six).
  - Five templates instead of three (`max_templates=5`, run `fit4_t5`): model 1 +0.0001, best of 5 +0.0003 on
    the 48 chains; three stay the default.
  Runs without a template are unchanged; template runs differ from 0.7.5 for the same seed.
  `template_fit_attempts=1` gives the former behaviour.
- **Confidence model checked blind on CASP16** (`benchmark/casp16/confidence_check.py`, `docs/V07_VALIDATION.md`
  §F): colour bands in the right order, the *high* band somewhat optimistic; models fitted onto a distant
  template are over-rated (the model reads coordinates only). Stated in the manual (§8.5). A correction with
  the template's rank score and choosing model 1 by confidence were tested and not adopted. The model itself is
  unchanged.
- **Variants (Switch tab)**: one entry that cannot be read no longer stops the whole list; the others are
  analysed and the unusable ones are listed with the reason. Positions are the residue numbers RCV shows (a
  structure whose file starts at residue 17 was addressed by counting from 1 while the table showed the file's
  numbers). `34G>U` and an `n.`/`r.` prefix are read; a variant entered twice counts once.
- **Target sites (Details tab) about 100 times faster, and every window on long RNAs.** The opening energy of
  all windows now comes from one pass over the ensemble (the unpaired probabilities of RNAplfold with window and
  pair span set to the full length) instead of one constrained partition function per window. The numbers are
  the same (difference below 0.0001 kcal/mol on RNAs of 76–3000 nt, at 37 °C/1 M and 25 °C/0.15 M Na⁺): 300 nt
  take 0.2 s instead of 20 s, 1500 nt 15 s instead of about 40 min. RNAs longer than 400 nt are no longer
  sampled at every n-th window. The former method remains as a fallback for a ViennaRNA without that interface.
- The Pull tool stands out in the toolbar: accent-coloured, its label shown at every window width, and a small
  dot marks it until it is first used.
- Base-pair annotation accepts a pair whose non-central hydrogen bond is missing only for three-bond pairs, so
  inosine–C (two hydrogen bonds) is annotated.
- Pairs across a chain break (O3′–P not bonded) are annotated even when the residues are fewer than four apart,
  and the 2D energy of such a model is the co-folding energy of its strands. The RNA chain option accepts
  several chains joined, e.g. `A+B`, to show both strands of a duplex.

### Fixed
- **`rcv batch` went silent after 60 notes** (audit F1). The command printed the notes it had not shown yet by
  counting them, but a job reports only its newest 60, so a batch of more than a few sequences stopped printing
  although it went on and finished correctly. Every note is printed now. Ctrl+C stops the running prediction
  and its worker processes before the command ends (exit code 130; measured: all processes gone after 0.2 s);
  before, an interrupt that reached only the command left the engine running.
- **Lines of the start page that stayed English in the other ten languages** (audit F2). The status line under
  the sequence, the time estimate, the progress card of a running prediction (per cent, time left, stages and
  their details, AlphaFold 3 steps, “Failed:”), the list and the results table of a batch, the running jobs in
  the viewer and the variants table of the Switch tab (headings, effect, explanation, “Not analysed” with its
  reasons) had no translation, or only for one wording of a line that is put together from parts: “3D sampling
  · 12 of 32 trajectories” came out half translated. Such lines are now built from separately translated
  parts, and 66 texts were added to each language file. The smoke scenario `langstart` walks these states in
  German and Japanese and fails on English text left in both.
- **The notes of a running prediction were English in every language** (audit F5), although the guide listed
  them as translated. The lines under the progress bar (“Sampling 32 independent trajectories on 4 CPU
  cores.”, “30 of 32 trajectories passed all-atom validation; …”, the template and 2D-structure lines, the
  notes of a pull exploration, a batch and an AlphaFold 3 run, and the engine's five warnings) now follow the
  interface language: 63 texts were added to each language file (with F8). A note tells the page which parts of it are
  data (the results folder, the title of a template structure, the reason a trajectory was rejected, a command
  line, a system error) and which name stands in front of it (the sequence of a batch); those parts are shown
  as they are. Lines printed by AlphaFold 3 itself and errors of other programs stay English. `rcv batch` on
  the command line prints the same English text as before. A test reads every note text from the source and
  fails when one has no entry in the language files.
- **The note naming a template never appeared** (audit F7). The engine's line about a template says in brackets
  how it was found (`… pairs 46/46 [sequence] 'title'`); RCV's reader of that line did not expect the brackets,
  so “Template 8UYJ_A (title): sequence identity 100%, covers 100% of the sequence, E-value …” was shown only
  for the additional template of a part the first one left uncovered. Checked on the 15 CASP16 logs: 28
  template lines, 28 notes (before: none of them).
- **Internal labels in the notes about 2D structures** (audit F8). With probing data, an alignment, a
  pseudoknot hypothesis or Z-RNA the note read “ss1 (vienna_mfe_probing, …)”, “alignment_consensus”,
  “vienna_pseudoknot”, “Z-RNA form of ss1 (vienna_mfe)”. They now read “minimum free energy with probing
  data”, “suboptimal with probing data”, “consensus of your alignment”, “with a pseudoknot” and “Z-RNA form of
  ss1 (minimum free energy)”, in every interface language. A test reads the labels from the engine's source and
  fails when one has no words.
- **Traditional Chinese: “passed” was written 透過** (“by means of”) in two entries (audit F9); now 通過 and
  通過全原子驗證.
- **The GUI smoke test lost a scenario now and then to a browser that never loaded the test page** (audit
  F10; six times on 2026-11-21, each ending after the full time limit without a line from the page). Under WSL
  the Windows browser kept its temporary profile on the Linux disk, which Windows reaches through
  `\\wsl.localhost`. Measured with 1200 bare browser starts: 7 of 800 with the profile there never loaded a
  page (five of them within the same 30 seconds, while starts with the profile on the Windows disk went on
  working), 0 of 400 with the profile on the Windows disk; a start took 5.6 s against 1.1 s. Turning the
  browser's crash reporter off did not help (2 of 400). `tools/viewer_smoketest.py` now keeps the profiles in
  the Windows temporary folder when it drives a Windows browser from WSL and removes them afterwards. As a second
  safeguard the test page reports that it has loaded; a scenario whose browser has not reported within 30 s
  (`--start-timeout`) is started once more, and the report names it. The move brought one new problem, found
  in the first full run and removed: with its profile on the fast disk the browser stood still for up to 13 s
  about 45 s after its start, which pushed the longest scenario (`langstart`, 49 s) to 59–85 s in six runs
  and once past its 90 s limit. The test browser is now started with its background work switched off
  (component updates, background networking, sync, timer throttling); `langstart` then took 45 s in three
  runs. Which of these switches is the one that matters was not separated. Full smoke test: 540 s, before
  about 730 s.
- **“2 min” and “40 s” kept the English unit letters inside translated lines** (audit F6), next to “1 小時 5
  分鐘” for longer times. Short durations are now written like long ones in Japanese, Korean, Russian and both
  Chinese variants (German, Spanish, French, Italian and Portuguese keep the international symbols). The entry
  for a bare unit takes a number only, so a sentence that ends in “ s” is not touched.
- **Cancel left the computation running** (audit A1). Cancelling a prediction stopped the engine process but
  not its worker processes, which went on using every chosen core until their trajectories ended (minutes),
  also after RCV was closed. The engine now ends its workers when it is told to stop, every worker ends by
  itself when the engine is gone, and RCV stops the whole process tree (`jobs.stop_process_tree`; a launcher
  script of a local AlphaFold 3 with its children at once). Measured with four workers: all gone 0.2 s after
  Cancel.
- **Pull exploration of an RNA with modified nucleotides lost the modifications** (audit A2): the options of
  the exploration replaced those of the sequence notation instead of being added to them.
- **Structures that did not open** (audit A3, A4; `structure_io.py`): entries whose models are told apart by
  the alternate-conformer label (1ZEV), mirror-image RNA (0A/0C/0G/0U; 2G32, 2GPM, 2GQ4–2GQ7), residues whose
  component declares a DNA parent while they are 2′-substituted RNA (3P4A), components without a declared
  parent (the base is read from the atoms), and PDB files that give the parent in MODRES or the residue in
  SEQRES only. Of two alternate positions of an atom the first in the file is used (before: only the one
  labelled A, so residues with conformers B and C lost atoms). A free base or another hetero group without a
  sugar is no longer taken for a residue of the chain. All 327 cached PDB entries of the benchmarks open
  (320 before).
- **“No valid model” was reported as a fragment of JSON** (audit A5): the message now says that no trajectory
  passed the all-atom validation and what to try; other engine failures show the engine's own sentence.
- **Batch predictions ignored the Z-RNA option** of the start page (audit A6); probing data and alignments,
  which belong to one sequence, are said not to be used in a batch instead of being dropped silently.
- **Seed 0 became seed 42**, and a word in a number field stopped the request with a Python error (audit A7):
  numbers are read in one place, with their range and a sentence for the user.
- **`rcv extract` could write outside the chosen folder** when an `.rcv` file carried conformation names with
  path separators (audit A8): names are reduced to letters, digits, `_ . -` and made unique.
- **Modified nucleotides checked against the wwPDB chemical component dictionary** (`modifications.COMPONENTS`,
  `tests/test_components.py`, reference data in `tests/data/pdb_components.json`). Every modified residue RCV
  writes under a PDB component name now is that component: same atoms, same bonds, same configuration at every
  stereocentre, and the dictionary's atom names (before, 14 of the 25 components were written with RCV's own
  names, e.g. `CM6` for the methyl that 6MZ calls `C9`, so other programs could not match the atoms).
  61 entries now have a component (25 before), among them LNA (LCA/LCC/LCG/TLN), 2′-fluoro
  (AF2/CFZ/GF2/UFT), m1Ψ (B8H), ac4C (4AC), s2U (SUR), lysidine (M3X), mnm5s2U (U8U), mcm5s2U (70U), 8-bromo-G,
  2,6-diaminopurine; 67 components are recognised when a structure is opened. Entries without a nucleotide
  component of exactly that chemistry keep the parent name. Consequences: more modifications are recognised in
  opened PDB entries, and more are passed to AlphaFold 3 (which takes them as PDB components).
  - queuosine is written as 56B (QUO, used before, is an analogue with a saturated ring; still read) and m7G as
    G7M (7MG is the 8-hydro form; still read);
  - **queuosine** was built as the wrong diastereomer (trans-diol, either hand at the ring carbon): now the
    (1S,4S,5R)-cis-diol with the ring geometry of the component;
  - **t6A, ms2t6A, m6t6A**: the side-chain rotamer search could turn L-threonine into D- or allo-threonine;
    a stereocentre now takes its second substituent relative to the first, so no rotamer can invert it;
  - **LNA**: the methylene bridge was 2.1 Å from C4′ (a bond is 1.56 Å). O2′ and the methylene are now fitted
    together to the bicyclic sugar (bonds within 0.05 Å on a C3′-endo sugar);
  - **lysidine**: the carboxyl carbon was named `CL` and written as chlorine;
  - both hydroxyls of queuosine can donate hydrogen bonds in the search (one before).
- **Z-RNA prediction: a left-handed model was always ranked last**, also when the engine's own estimate said
  Z is the major form (r(CG)₆ with a Z-binding protein: ΔG(A→Z) −1.5 kcal/mol, 92 % Z, yet models 4 and 5).
  The 3D score is calibrated on A-form helices and scores a Z stem about 50 units worse. The best Z model is
  now set level with the best model of its parent structure, so the estimated ΔG(A→Z), already part of the 2D
  energy, decides: model 1 when Z is favoured, next to the first models when it is a minor form
  (`helix_form_score_excluded` and `z_fraction_estimate` in the manifest; `docs/RCV.md` §3.4).

### Removed
- **Check for updates** (0.7.5): *Help ▸ Check for updates…*, the daily check at start-up, its Preferences
  option, `/api/update`, `src/rcv/update.py`, `static/js/update.js`, `tests/test_update.py` and its texts. RCV no
  longer contacts any server on its own and never points to or opens an installer; Store copies are updated by
  the Store, the other installers from the website. (Store policy does not allow apps to update themselves, and
  programs that fetch installers are what antivirus heuristics look for.) The Inno Setup script no longer sets
  `AppUpdatesURL`.

## Package 0.7.5 — 2026-11-15 — RCV · RNA Conformation Viewer (engine rnafold3d 2.3.0)

First public release since 0.7.3: it also contains everything listed under 0.7.4, which was not released on its own.

Conformation-selective compound screening (`docs/NEXT_STEPS.md` §1), steps 1 to 5: where a compound could
bind, how those sites differ between two conformations, a scoring function and a docking search for RNA, and
the screen itself — a compound library ranked by how much better it fits one conformation than the other, in
its own window and on the command line. It is marked *experimental* and shows what it was measured at,
because one of its three benchmarks is still missing (`docs/SCREENING_STATUS.md`).

### Added
- **Pocket finding** (`src/rnafold3d/screening/pockets.py`, a new `rnafold3d.screening` package). A buriedness
  scan of the free space around the atoms on a 0.8 Å grid (the LIGSITE approach: a point counts as enclosed in
  a direction when the RNA closes the space on both sides within 12 Å, along the three axes and the four cube
  diagonals), clustered into sites of at least 100 Å³ enclosed on at least five of the seven directions. Each site has its volume, how enclosed it is, the
  nucleotides lining it, how much of that lining lies outside a Watson-Crick helix, and a score that ranks the
  sites (enclosed × outside-a-helix × size). `find()` takes plain coordinates and depends only on numpy and
  scipy; a 76-nt tRNA takes 0.2 s, a 158-nt structure 1 s.
- **`benchmark/screening/pocket_sites.py`**: thirty experimental RNA-compound complexes with the compound
  removed; the site it occupied is the top-ranked pocket in 21 of 30 and among those found in 24. This is
  what set the threshold: at 6 of 7 enclosing directions it is 9 and 12, because most RNA sites are shallow.
- **`volume_near()` and `compare()`**: the same site measured in two conformations — the open volume held by
  the nucleotides that line a pocket — so pockets that merge or split between the conformations are still
  compared, and the superposition of the two structures does not matter. A site open in one conformation and
  closed in the other is what a conformation-selective compound would use.
- **Details tab ▸ *Pockets: where a compound could bind*** (closed *advanced* section, as in 0.7.4): the pocket
  table of the shown conformation, *Compare with* a second conformation and the sites that differ between the
  two, *Colour 3D by pocket*, copy as a table. A row selects the nucleotides lining the pocket in 3D and 2D.
  *only deep pockets (6 of 7)* keeps the closed sites only. Manual: `docs/RCV.md` §8.6.
- Server: `POST /api/insight/pockets` (`conformation`, `other`, `burial_min`); it does not need ViennaRNA.
- Tests: `tests/test_v075.py` (15 tests, among them the cavity of a hollow sphere against its analytic volume);
  GUI smoke test: the `insightpockets` scenario (62 in total).
- Interface texts of the new section in all 10 languages.
- **Window layout** (`static/js/windows.js`; `docs/NEXT_STEPS.md` §4). **Panel widths**: drag the inner edge of
  either side panel (left 190–520 px, right 312–760 px; the 3D view keeps at least 360 px), double-click for the
  default; widths are remembered per screen size; *Window ▸ Default panel widths*. **Rails**: a hidden side
  panel folds to a 36-px strip of icons instead of disappearing; one click opens it, on the right on that tab.
  **Presentation mode** (`F5`, Window menu): only the 3D view, full screen, no menu bar, large labels, title and
  legend; `Esc`/`F5` or an *Exit presentation* button (shown when the mouse moves) leaves it. Texts in all 10
  languages; GUI smoke test: `layout1`, `railshot`, `presentshot`.
- **Project library** (*File ▸ Project library…*, `Ctrl+L`, or the start page; `src/rcv/library.py`,
  `static/js/library.js`). Every prediction run, `.rcv` file and structure that has been opened, plus what is
  found in the watched folders, with a thumbnail, the parameters and a search over names and sequences; nothing
  is copied or moved, the index points to the files where they are (`<user data>/library/index.json`). Manual:
  `docs/RCV.md` §11.
- **Plain language in the console** (`static/js/console/plain.js`). A line that is not a command is read as
  plain English and shown as the commands it stands for ("colour by flexibility" → `color flexibility`), which
  are then recorded in the command log, so the log stays a reproducible script. Local patterns only: no language
  model, nothing is sent anywhere; English, like the commands. `help plain`. Manual: `docs/RCV.md` §9.1.
- **Recipes** (*Tools ▸ Recipes…*; `src/rcv/recipes.py`, `static/js/recipes.js`). A named analysis — prediction
  settings plus the steps that follow — in a small `.rcvrecipe` file (JSON) that a lab can share, so everyone
  looks at an RNA the same way. Four built-in recipes; user recipes in `<user data>/recipes/`. Server:
  `/api/recipes`. Manual: `docs/RCV.md` §9.2.
- **Check for updates** (*Help ▸ Check for updates…*; `src/rcv/update.py`, `static/js/update.js`). Finds the newest
  release and the installer for this computer (Windows, macOS Apple Silicon / Intel — also an Intel build under
  Rosetta —, Linux `.deb` / `.rpm` by how RCV was installed or the distribution) and asks: *Download* opens it in the
  web browser, *Skip this version*, *Later*, *Release notes*. Nothing is installed automatically. Installed copies
  check once a day at start-up (off in the update window or Preferences, or with `RCV_NO_UPDATE_CHECK=1`); source
  runs only when asked. Source: the GitHub release of the website repository, or the website's `stats.json` when
  GitHub is unreachable. Texts in all 10 languages; tests: `tests/test_update.py` (10), GUI smoke test `update1` (the
  smoke test turns the start-up check off; 66 scenarios in total).
- **Zenodo DOI** ([10.5281/zenodo.24163503](https://doi.org/10.5281/zenodo.24163503)): in the software citation
  of *Export › Methods text*, a DOI badge and a *How to cite* section in the README, and `CITATION.cff` (GitHub's
  *Cite this repository*).

- **A scoring function for RNA and a docking search** (step 2 of the plan), as a library and two commands:
  - `screening/chem.py`: compounds from SDF, MOL, MOL2 and PDB entries — rings, aromaticity, rotatable bonds,
    the torsion tree, and the charge of the groups that are ionised at pH 7 (a carboxylate is −1 whether the
    file drew it protonated or not, a thiazolium is +1, the N9 of a guanine is not); an SDF writer for poses.
  - `screening/receptor.py`: the RNA typed from its atom names (which nitrogen donates, which carbon is
    apolar), one −1 per phosphate, and a plane per base; `read_complex` splits a PDB entry into the RNA and
    the compounds, keeping modified nucleotides in the chain and joining an aminoglycoside's sugars into one
    compound.
  - `screening/score.py`: AutoDock Vina's shape and contact terms (Apache-2.0, used as a reference) plus the
    three an RNA score needs — stacking against the base planes with the angle between them, screened
    electrostatics with explicit phosphate charges, and the cost of burying a charge.
  - `screening/grid.py`: the whole function precomputed on a 0.375 Å grid, so a pose costs one interpolation.
    Stacking stays directional on the grid by storing Σ w·nnᵀ per point, which gives back cos² of the angle
    between the compound's ring and each base exactly.
  - `screening/dock.py`: Monte Carlo with local minimisation over position, orientation and torsions,
    deterministic for a seed, poses rescored from the atoms rather than from the grid.
  - `rnafold3d pockets` and `rnafold3d dock` (poses out as SDF); manual: `docs/RCV.md` §17.1.
  - `benchmark/screening/redock.py`: thirty experimental RNA-ligand complexes downloaded from RCSB, each
    ligand docked back into its own structure, reported as RMSD to the experimental pose. **Result: median
    4.24 Å, under 2.5 Å in 10 of 30** — below what the plan asks of a scoring function, which is why the
    screen is marked experimental and shows this number in its own window. In 24 of the 30 the
    search found a pose scoring at least as well as the experimental one, so what needs work is the ranking,
    not the sampling. It redocks small rigid compounds well (2B57 0.27 Å, 1F1T 0.32 Å, 3SKI 0.61 Å,
    1Y26 0.87 Å, 4NYB 0.89 Å) and large flexible ones in shallow sites badly (every aminoglycoside and
    every SAM entry is 3.9 to 9.0 Å).
  - `docs/SCREENING_STATUS.md`: what is measured, what is not, and the weights that are set rather than fitted.
  - Tests: `tests/test_screening.py` (no network).

- **The screen itself** (steps 3 to 5 of the plan): a compound library against two conformations.
  - **Tools ▸ Compound screening…** (also a button in the Switch tab): its own page, so the main layout is
    unchanged — target (the two conformations and a pocket, or the nucleotides selected in the main window),
    library, parameters, run with live progress, and a sortable, filterable table that exports as CSV or as an
    SDF of poses. *Show in 3D* draws the pose in the main window with the nucleotides it touches selected.
    Four example compounds ship with RCV so a first run takes a minute; no library does. Manual: `docs/RCV.md`
    §8.7. Texts in all 10 languages; GUI smoke test: `screenwindow` (a real run of four compounds against two
    conformations) and `screenpage`.
  - `screening/screen.py` and `rnafold3d screen`: the two-state screen. **Selectivity = score(B) − score(A)**
    for the best pose in each state, with state B measured at the nucleotides that line the site in state A —
    not at a pocket of its own, because a site that has closed in B has no pocket to match. Compounds with the
    same selectivity are ranked by the worse of their two scores. One worker per core over compounds, the grids
    built once and shared; a site wider than 30 Å is cropped to its centre and says so.
  - `src/rcv/screening.py` and `POST /api/screen` (`/library`, `/site`, `/<id>/results`, `/csv`, `/sdf`,
    `/pose`): the libraries of a session, the site the user picked, and the run as an ordinary RCV job, so
    progress and cancelling work as they do for a prediction.
  - **Benchmark 2 — enrichment** (`benchmark/screening/enrichment.py`): six RNA targets (the FMN, TPP, adenine,
    guanine and preQ1 riboswitches and the ribosomal A site) with actives collected from the literature and 25
    property-matched decoys per active drawn from PubChem, everything 2D-similar to an active removed. The raw
    value of every term is kept for every pose, so `tune` searches the weights by reweighting stored poses
    (a dot product) and reports leave-one-target-out, never in-sample only.
    **Result: mean AUC 0.864 over the six targets** (34 actives, 700 decoys, 49 minutes on 12 cores; a random
    ranking is 0.500 and the plan asks for 0.65). Per target: adenine 0.997, preQ1 0.992, A site 0.972,
    guanine 0.871, FMN 0.855, and **TPP 0.495, which is exactly chance** — the one ligand in the set whose
    site holds a diphosphate through two Mg²⁺, i.e. the terms this function treats most crudely. The A site's
    0.972 is flattered by its decoys sitting 3.38σ from its actives in property space, because PubChem has
    few 3D records for compounds as large and polar as the aminoglycosides. The screening window shows the
    mean AUC next to its results.
    **The weights were tuned against this benchmark and then left alone**: fitting the four RNA terms within a
    factor of 4 of their defaults gained +0.024 AUC in sample but only +0.010 cross-validated, with three of
    the four weights pinned to the edge of the allowed range and two of the six targets made worse. That is not a
    change a 34-active set can justify, so the shipped weights are still the ones set by hand.
  - Tests: `tests/test_v075_screening.py` (7 tests: library parsing, the site on both conformations, a whole
    run through the routes, the exports and the pose).

- **Engine accuracy, rnafold3d 2.3.0** (`docs/NEXT_STEPS.md` §2). Every item was measured on the 26 dev chains
  (`benchmark/dev_templates.py`, 16 trajectories, seed 42) with each feature on and off. What is on by default
  is what improved the numbers; the rest is there as an option, with the measurement that says why it is off.
  - **Motif fragments** (`fragments.py`, on by default; `--fragments off`). Hairpin loops of 3–12 nt take their
    backbone torsions and C4′ distances from library loops with a matching sequence (loop bases weighted double,
    the closing pair included, homologues of the query excluded), sampled per trajectory; the strands of a grafted
    junction now carry their torsions as well. Grafted loops are closed in the fragment's shape. **Dev chains
    without a template: top model TM-score 0.216 → 0.227; all 26: 0.284 → 0.291** (best of five 0.303 → 0.307);
    no change in run time. Better on 1Y26, 480D, 4L81, 4ENC, 1F1T, 4QLM (+0.03 to +0.05), worse on 2OIU (−0.04),
    1U9S, 3V7E (−0.02).
  - **Pseudoknots in the search** (`pseudoknots.py`, off by default; `--pseudoknots auto`). A 2D with one
    crossing helix, folded jointly with the nested part (each stable helix of ≥ 4 pairs is blocked, the rest
    refolded by ViennaRNA, and the helix kept when it crosses the result and wins by a 4 kcal/mol initiation
    penalty). The penalty was calibrated on 1,509 library chains (`benchmark/pseudoknot_calibration.py`): base-pair
    F1 0.624 → 0.631. In 3D the crossing helix becomes A-form distance and base-pair restraints, its loops are
    built explicitly from the start. On the dev chains it is found in 8 of 26 and is correct where it can be
    checked (1RNK 6 of 7 pairs, 4ENC 6/6, 3D0U 6/9), but the helix rarely forms in the models: top model 0.284 →
    0.287 alone, nothing on top of the fragments (0.291 both ways), and runs that have one take up to twice as
    long — so it stays off.
  - **Long RNAs** (`domains.py`, experimental, off; `--domains auto`). Chains of ≥ 240 nt are cut at the exterior
    loop and at the first multiloop of long helices into domains of ≤ 160 nt, each predicted as its own chain,
    and the models become a multi-part template for the whole-chain search. On six library chains of 255–393 nt
    without templates (`benchmark/long_rna.py`) it was **worse** than the whole-chain search: top model 0.091
    against 0.104, valid trajectories 2.5 against 3.3 of 16, 1.8× the time. Both are near the level of
    unrelated structures at that length.
  - **Knowledge-based base–base potential** (`basepot.py`, `data/basepot_v1.npz`, 9 kB; built by
    `tools/build_base_potential.py` from 3,315 library chains with homologues of the dev chains removed): base
    centroid distance × relative orientation by nucleotide pair and sequence separation. **Not in the ranking**:
    refitting the ranking on the saved candidates (`benchmark/rerank_basepot.py`) changes nothing up to a weight of
    3 and gains 0.007 only at a weight of 40. Random choice among the candidates gives 0.220, the engine's ranking
    0.221 and the best candidate 0.254 — what limits the targets without a template is sampling, not ranking.
  - **Speed.** No compiler or numba in the build environment, so numpy only: the motif index of the library is
    built 10× faster (6.2 → 0.6 s per run, `motifs.nested_greedy`); the clash term avoids square roots it does not
    need. A 2-trajectory profile of 1Y26 goes from 82 to 76 s; the trajectory inner loop itself is not faster.
    A compiled inner loop is still to do.
  - `rnafold3d predict --fragments / --pseudoknots / --domains auto|off`; `benchmark/dev_templates.py --set
    key=value` (any engine option) and `--candidates`.
  - Tests: `tests/test_v075_engine.py` (14 tests).
- **SAXS, FRET and cryo-EM maps** (`src/rcv/expdata.py`, `static/js/expdata.js`; *Tools ▸ Load SAXS, FRET or
  cryo-EM data…*). The same path as the probing data: what each conformation predicts, its agreement, and a
  maximum-entropy reweighting of the ensemble, several datasets jointly. SAXS: Debye profile from all heavy atoms
  with form factors minus displaced solvent, fitted scale and background, χ² and Rg; FRET: distances or
  efficiencies with R₀; maps: MRC/CCP4, each conformation placed in the map and its simulated density correlated
  with the map above the contour. Results in a closed section of the Details tab (curve, metrics, data weights);
  saved in `.rcv` files; Methods text paragraph. Server: `POST /api/expdata/parse`, `/api/expdata/fit`. Texts in
  all 10 languages; GUI smoke test: `expdata1`, `expdatashot`. Manual: `docs/RCV.md` §8.4.

### Fixed
- **Charge perception**: an aromatic nitrogen with three bonds was treated as a cation whenever its ring had no
  pyrrole position, which made every N-methylated ring amide a cation — theophylline came out +2, riboflavin
  and FMN +1 too high. It is now charged only when it also holds a double bond inside the ring (a pyridinium, a
  thiazolium), with the old rule kept for files that bring no bond orders at all (a compound out of a PDB
  entry). An amine on an aromatic ring is left neutral (an aniline is not protonated at pH 7). Benchmark 1 was
  re-run with the corrected charges: two cases that were being docked as the wrong molecule came right (1F1T
  malachite green 9.39 → 0.32 Å, 2GDI thiamine pyrophosphate 10.69 → 2.80 Å), four others moved the other way
  (4TS2, 3D2G, 1J7T, 2TOB), and the median went 3.97 → 4.24 Å. The fix is right on its own terms; a 30-case
  set cannot resolve a change this small.
- `rnafold3d pockets --burial-min` defaulted to 6 while everything else used the measured 5.

### Not done yet
- Benchmark 3 (selectivity): the published conformation-selective compounds, the thing the screen is for and
  the one with almost no data. Until it is run the screen stays marked *experimental*, and every score in it
  is a ranking, not an affinity, with nothing in kcal/mol anywhere.
- Compound libraries: `tools/fetch_library.py`, the `.rcvlib` format and its manifest. Nothing ships but the
  four example compounds.
- Tautomers, ring conformers and protonation beyond the groups in `chem.ionised` — a compound is docked as its
  file draws it.
- Metals and phosphates: with a 6-point metal term and no explicit ions, a site that holds its ligand through
  bridging Mg²⁺ ranks no better than random (the TPP riboswitch, AUC 0.495). An explicit-ion model is what
  benchmark 2 says is missing.
- Engine: complexes (two chains) as the plan says; a compiled fold-tree inner loop; making pseudoknot helices
  actually form in 3D; a domain assembly that beats the whole-chain search.
- Speed: a compound takes seconds, not the 0.3 s the 100,000-compound target needs. The analytic gradient and
  a compiled inner loop are the next performance step.

## Package 0.7.4 — 2026-11-12 — RCV · RNA Conformation Viewer (engine rnafold3d 2.2.0; not released on its own, part of 0.7.5)

Six advanced analyses inside the existing tabs (the interface is otherwise unchanged), and an engine that finds
more distant templates, keeps models from different templates, and can use your own probing data or alignment.

### Added
- **Advanced analyses, as closed sections at the bottom of the Switch, Pathway and Details tabs** (marked
  *advanced*; nothing is computed until a section is opened and its button pressed; rows select the nucleotides
  in 3D and 2D; tables copy as tab-separated text). Manual: `docs/RCV.md` §8.6.
  - *Switch ▸ Experiments that tell A from B*: disrupt-and-rescue double mutants (shift ≥ 1 kcal/mol, rescue
    ≥ 50%), the segments where SHAPE/DMS reactivity should differ, and the nucleotide pairs whose distance
    differs most between the 3D models of the two states (FRET, PELDOR/DEER).
  - *Switch ▸ Variants*: paste substitutions (`G34U`, `A12C+G30C`); change of the pairing profile (the measure
    of riboSNitch studies), pairs gained and lost, ΔΔG and the A share, with an effect class.
  - *Pathway ▸ Folding during transcription*: the most stable structure of every nascent transcript, pairs kept
    or temporary, how much of A and B has formed at each length, and the lengths where the structure rearranges.
  - *Details ▸ Structural motifs*: hairpin loop types, bulges and internal loops with the bend between helices,
    junctions with coaxial stacks, pseudoknot helices, non-canonical pairs, base triples and A-minor contacts
    from the atom positions.
  - *Details ▸ Target sites*: the free energy to open every 14–25-nt window, the share of conformations in which
    it is unpaired, GC content and the antisense sequence, for antisense oligonucleotides, siRNAs and probes.
  - *Details ▸ Compare with another ensemble* of the same RNA open in RCV: superposition weighted towards the
    part that stays the same, displacement per nucleotide, regions that moved, RMSD, TM-score, local agreement
    and base-pairing differences.
  - Server: `POST /api/insight/<experiments|variants|cotranscription|motifs|targets|compare>` (`rcv/insight.py`).
- **Structure-aware template search** (`rnafold3d/templates.py`, on by default; `PredictConfig.template_structure`).
  Library chains whose sequence alignment is weak (E-value up to 10) are aligned by sequence and secondary
  structure together (the query's pairing probabilities from ViennaRNA against the template's base pairs) and
  kept when the alignment covers ≥ 60% of the query and ≥ 90% of the template pairs it covers can form in the
  query. The same alignment extends sequence hits that stopped where the sequences diverge; it is skipped when
  the hit leaves an internal part of ≥ 20 nt uncovered (an insert in a scaffold).
- **Models from different templates.** The kept models now take the best model of each template in turn before
  a second model of the same template, so the five models cover the templates found.
- **Your own data in the prediction** (start page ▸ *Advanced settings*; `rnafold3d predict --shape FILE
  --alignment FILE`): SHAPE reactivities as pseudo-energies for the 2D hypotheses (Deigan, 1.8 and −0.6
  kcal/mol), and an alignment of homologs whose consensus structure (RNAalifold) becomes one of the hypotheses.
- `benchmark/template_selection.py`: leave-one-out benchmark of the template search on 916 library chains;
  `benchmark/casp16/compare.py paired`: per-target differences between two runs with bootstrap intervals.
- Tests: `tests/test_v074.py` (19 tests); GUI smoke test: three new scenarios (61 in total).
- Interface texts of the new sections in all 10 languages.

### Fixed
- An advanced section closed itself again when the tab was drawn a moment after the click (the Pathway tab
  finishing its analysis, for instance).

### Changed
- **Benchmarks, 0.7.4 against 0.7.3.**
  - Template search (916 chains, each searched against the library without its own class): mean TM-score of
    the first-ranked template 0.416 → 0.428, of the best template found 0.444 → 0.461; chains with a template
    642 → 655.
  - Development benchmark (26 PDB chains, templates released before each one): mean TM-score of model 1
    0.271 → 0.284, best of 5 0.292 → 0.302; chains with a template 5 → 7.
  - CASP16 (15 targets): model 1 unchanged (0.489; 0.7.3: 0.490 — the first model comes from the same template
    on every target); best of 5 0.511 → 0.524, from the tRNA R1271 (0.50 → 0.66) and R1203 (0.49 → 0.52). The
    95% bootstrap interval of the mean best-of-5 difference is −0.001 to +0.036: with 15 targets this is not
    established. 1–9 minutes per target on the 16-core laptop CPU.
- CASP16 report, figures and ensemble page (`benchmark/casp16/results/`) regenerated with 0.7.4.
- `docs/NEXT_STEPS.md`: the planned functions that are not implemented, in detail — conformation-selective
  compound screening (its own sub-window, scoring core, pocket finding, validation benchmarks and libraries),
  engine accuracy, retention features and the remaining layout work.

- **A more compact window.** Menu bar 38 → 32 px, toolbar 46 → 36 px with smaller controls, panel headers
  36 → 30 px, playback bar 68 → 54 px, narrower tab and list padding, side panels 268/360 → 244/352 px (the
  panels hold the same content width as before, the padding is what shrank). Above 2,200 px — 27-inch monitors
  and up — the side panels are 252/372 px, where there is room for them.

### Not changed
- The rule that ranks templates (identity × coverage × 3D agreement): no alternative tried on the template
  benchmark was better.
- Without a template, the 3D arrangement of helices is modelled as in 0.7.3.

## Package 0.7.3 — 2026-11-07 — RCV · RNA Conformation Viewer (engine rnafold3d 2.1.0)

3D interactions in the engine: related experimental structures (templates) and junction geometry from known
structures, without deep learning, in a 21 MB library. AlphaFold 3 can be used as a second engine (your own
installation) and AlphaFold Server results open directly.

### Added
- **Template-based modelling with a built-in structure library.** `rnafold3d/data/pdb_rna_v2.npz` holds 3,964 RNA
  chains from the PDB (one per non-redundant class of the BGSU RNA 3D Hub, release 4.64): sequence, coarse
  coordinates (C4′ plus P, C1′ and four base atoms), torsions and base pairs, 21.5 MB, no network needed. The
  query is aligned to all chains in 1–3 s (exact Smith–Waterman with BLASTN scoring and E-values, collinear
  multi-segment alignment). Hits must have E ≤ 10⁻³ (≤ 10⁻² when they cover ≥ 60% of the query) and cover ≥ 45%;
  they are ranked by identity × coverage × 3D agreement with the other hits, which moves templates in unusual
  states down. When the best template leaves a part of ≥ 20 nt uncovered (an insert in a scaffold), a template
  for that part is added. A template gives the 2D (the rest folded by ViennaRNA with the template pairs
  enforced), C4′ distance restraints, starting torsions and atom positions. Development benchmark (26 PDB
  chains, templates released before each one; 5 had a template): mean TM-score of model 1 0.271, e.g. tRNA-Phe 1EHZ
  0.89.
- **Cartesian refinement** (`rnafold3d/cartesian.py`). After the fold-tree search, every atom is pulled onto
  the superposed template atoms (each template part on its own) while bond lengths, angles, ring geometry and
  chirality stay those of the ideal model and atom overlaps are pushed apart (distance restraints, L-BFGS,
  about 1 s); a polish restores ideal covalent geometry. A template model reproduces its template to 0.2–0.4 Å
  on the restrained atoms (R1263: TM-score 0.93, the template itself 0.93).
- **Restrained loop closure.** Loops of template trajectories are closed by a Levenberg–Marquardt fit of the
  loop torsions to the closure, the template restraints and the template torsions together (CCD alone changed a
  few torsions by up to 180° and lost the loop shape).
- **Junction geometry from known structures (template-free).** Every internal loop, bulge and multi-way junction
  of the 2D gets, per trajectory, the geometry of a library junction with the same helix count and strand
  lengths (homologues of the query excluded), as soft C4′ distance restraints; the ensemble samples junction
  arrangements seen in experimental structures. The restraints only guide the sampling: their energy is left out
  of the score that ranks the models (`guide_energy_excluded` in the manifest). Development benchmark (26 PDB
  chains, structures released before each one): mean TM-score of model 1 0.230 with junction geometry,
  0.214 without (better on 12 chains, worse on 5).
- **Regularisation of strained models.** Trajectories that fail validation only through local strain
  (closure-junction angles, a stretched O3′–P bond, a few severe overlaps) are repaired in Cartesian space
  (atoms move 0.05–0.1 Å) instead of being discarded.
- **AlphaFold 3 as a second engine** (start page, card 1, *Engine: AlphaFold 3 (local)*): RCV writes the
  AlphaFold 3 input, runs your installation (Python + `run_alphafold.py` + model parameters, optional
  databases, or a custom command for Docker/WSL/clusters), shows its log and opens all sampled models as one
  ensemble coloured by pLDDT. AlphaFold 3 is not included in RCV.
- **AlphaFold Server results** (card 2): open the downloaded `.zip` (or a local AlphaFold 3 output folder); all
  models open as one ensemble, best ranking score first, with pLDDT, ranking score, pTM and fraction disordered.
- Start page: *Use related experimental structures (templates)* in the advanced settings (on by default);
  `rnafold3d predict --templates off`, `--template-cutoff YYYY-MM-DD`, `--template-exclude`.
- Details tab: the template of every conformation, with a link to the PDB entry; progress notes name the
  templates found.
- `tools/build_template_library.py` builds the library from the BGSU list and the RCSB ModelServer
  (`--split-date` keeps a pre-CASP16 member of every class; `--repack` converts older files).

### Changed
- **CASP16 benchmark** (`benchmark/casp16/`): RCV is run with templates (only structures released before
  2024-05-01, the start of the CASP16 prediction season) and reported next to 0.7.2.
  Mean TM-score of model 1 over the 15 targets: 0.490 (AlphaFold 3 0.480, trRosettaRNA 0.493,
  Vfold 0.534, GuangzhouRNA 0.526; 0.7.2: 0.251); right global fold (TM ≥ 0.45) on 8 of 15
  targets (AlphaFold 3: 8; 0.7.2: none). 1–7 minutes per target on a 16-core laptop CPU without GPU, at most
  1.3 GB of memory; the same models with 4 cores (deterministic), in about twice the time.
- New: `ensembles.py` (CASP16 RNAs solved in several conformations: how well the five models of each method
  cover them), `timing.py` (time and memory per RNA), `figure.py` (four-panel comparison of size, hardware,
  speed and RNAs that move), `benchmark/dev_templates.py` (development benchmark with per-target date cutoffs).
- Engine 2.1.0: runs differ from 2.0.2 for the same seed (templates and junction geometry are on by default;
  `--templates off` plus `PredictConfig(motifs='off')` give the 2.0 behaviour).

### Fixed
- CASP16 target R1296 is HIV-1 RRE stem-loop II with G34U (R1209's RNA), not a tRNA: label corrected in
  `benchmark/casp16/targets.json`.

## Package 0.7.2 — 2026-11-06 — RCV · RNA Conformation Viewer (engine rnafold3d 2.0.2)

Bug fixes and a blind-test comparison with the leading RNA 3D structure predictors (CASP16).

### Fixed
- **Stems next to G-U wobble pairs were sometimes built with clashing atoms, so whole folds were lost.** The
  rigid A-form stem builder could squeeze the pair next to a wobble (G26-C41 O6-N4 at 2.05 Å in CASP16 R1205),
  overlap the two bases of a terminal wobble (1.5 Å) or push stacked neighbours together. Stems stay rigid during
  the search, so every trajectory on such a secondary structure failed validation and the hypothesis silently
  dropped out of the ensemble — even when it was the most stable (MFE) fold. Affected: about 5% of stems that
  contain a G-U pair (81 of 1600 random test stems). The refinement now adjusts both strands, pushes compressed
  pair hydrogen bonds and close contacts apart explicitly and checks stacked neighbours as the validator does:
  0 of 1600 stems clash, and the pair geometry fits better (0.15–0.2 Å instead of 1.4 Å RMSD in the failing
  cases). R1205: 18 of 20 trajectories valid (was 4), 5 models (was 4), MFE fold included. Test:
  `tests/test_v072.py`.
- **Command-line tools no longer print Python tracebacks for ordinary mistakes.** `rcv info`, `extract`,
  `convert`, `html`, `batch`, `pathway`, `probing`, `confidence` and `export-document` report a missing file, a
  file that is not an .rcv file, or a folder given to `rcv info`/`extract` (with a hint to `rcv convert`) in one
  line and exit with code 2. `RCV_DEBUG=1` shows the full traceback.

### Added
- **CASP16 benchmark** (`benchmark/casp16/`): RCV predicts the 15 single-chain CASP16 RNA targets up to 250 nt
  from sequence alone, and every model — RCV's and the submitted models of AlphaFold 3, trRosettaRNA, Vfold and
  GuangzhouRNA — is scored with the same US-align call as the official CASP16 assessment (`compare.py verify`
  reproduces the official table). Results: `benchmark/casp16/results/`, summary in the README.

## Package 0.7.1 — 2026-10-27 — RCV · RNA Conformation Viewer

### Changed
- **Start page: the two ways to use RCV, side by side.** Card 1, *Predict RNA conformations from sequence(s)*: the
  sequence box (one sequence, or several in FASTA format), *Predict conformations*, *Load a FASTA file…* and the
  advanced settings. Card 2, *View experimental or predicted structures of one RNA*: a large *Open files…* target
  (files can be dropped on it), *Fetch PDB…*, *Demo ensembles*, *Browse folders…* and *All prediction runs…*. The
  new text is in all 10 interface languages.
- **Command line: structure files given together open as one ensemble.** `rcv models/*.pdb` (or
  `rcv a.pdb b.cif …`) opens one ensemble with one conformation per file (one per model for multi-model files);
  before, every argument opened as an ensemble of its own. `.rcv` files, prediction folders and folders of
  structures still open on their own. Files of different RNAs given together open as one ensemble per RNA, with a
  note; `--separate` opens every file on its own.
- **Wildcards work on Windows.** The Command Prompt and PowerShell pass `models\*.pdb` on unexpanded, so RCV
  expands it itself, for `rcv`, `rcv convert`, `rcv html`, `rcv pathway`, `rcv probing`, `rcv confidence`,
  `rcv export-document`, `rnafold3d import` and `rnafold3d validate`. A pattern that matches nothing gives a short
  message.
- Structure files are read in file-name order everywhere (`model_2` before `model_10`), whatever order the shell
  passes them in, as already for folders and the Open dialog; the first is the reference for superposition.
  Duplicate paths are read once.

### Added
- `DifferentRNAError` (`rcv.data`) for structures that are not all the same RNA; `rcv.app.expand_paths` and
  `rcv.app.ensemble_groups`. Test: `tests/test_app.py::TestCli::test_structure_files_open_as_one_ensemble`.

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
