// Start page: the two ways to use RCV (1 predict from sequence(s) with the RCV engine or a local AlphaFold 3,
// 2 view structures of one RNA, including AlphaFold Server results); progress of local predictions.
import { el, slug, loadPref, savePref, readFileText, fmt } from './util.js';
import { icon } from './icons.js';
import { toast } from './ui.js';
import { watchJob, cancelJob, fmtDuration } from './jobs.js';
import { af3SettingsDialog, af3Configured, af3Summary, loadAf3Settings, AF3_SERVER } from './af3.js';
import * as modseq from './modseq.js';

const SAMPLING = { fast: [2, 6, 'Fast'], standard: [4, 8, 'Standard'], thorough: [8, 16, 'Thorough'] };
const trajSeconds = n => 45 + n + 0.004 * n * n;          // same rough model as the job runner

export const DEMO = { name: 'HIV1_TAR_demo', sequence: 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC', label: 'HIV-1 TAR, 29 nt' };

// FASTA (several records) or one bare sequence -> [{ name, sequence, structure }] (mirrors rcv.jobs.parse_fasta).
// A sequence may hold modified nucleotides in brackets (GGC[m6A]CU, modseq.js): `sequence` keeps them as written,
// `length` and the checks use the parent bases, `mods` lists them.
export function parseFastaRecords(text) {
  const recs = [];
  let cur = null;
  for (const raw of String(text || '').split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith(';') || line.startsWith('#') || line === '//') continue;
    if (line.startsWith('>')) { cur = { name: line.slice(1).trim().split(/\s+/)[0] || '', sequence: '', structure: '' }; recs.push(cur); continue; }
    if (!cur) { cur = { name: '', sequence: '', structure: '' }; recs.push(cur); }
    const tok = line.split(/\s+/)[0];
    if (/^[().\[\]{}<>]+$/.test(tok) && cur.sequence) cur.structure += tok;
    else cur.sequence += modseq.cleanLine(line);
  }
  const seen = new Map();
  return recs.filter(r => r.sequence).map((r, k) => {
    const base = slug(r.name || `seq${k + 1}`);
    const n = seen.get(base) || 0;
    seen.set(base, n + 1);
    const p = modseq.parse(r.sequence);
    const len = p.seq.length;
    const bad = p.bad;
    const problem = p.problems.length ? p.problems[0] : bad.length ? `letters ${bad.join('')}` : len < 8 ? 'shorter than 8 nt' : len > 400 ? 'longer than 400 nt'
      : r.structure && r.structure.length !== len ? 'dot-bracket length differs' : r.structure ? dotBracketProblem(r.structure, len) : null;
    return { ...r, name: n ? `${base}_${n + 1}` : base, problem, length: len, mods: p.mods };
  });
}

// Why a dot-bracket cannot be used, in words, or null when it is fine (mirrors rnafold3d.annotate.check_dotbracket).
// A line made of parts, joined by " · ". Each part is its own text node, so it is translated on its own whatever stands
// next to it (the interface languages translate whole text nodes; see i18n.js).
function setParts(node, parts) {
  const shown = parts.filter(Boolean);
  node.replaceChildren(...shown.flatMap((p, i) => (i ? [' · ', el('span', { text: p })] : [el('span', { text: p })])));
}

export function dotBracketProblem(ss, n) {
  if (n && ss.length !== n) return `dot-bracket has ${ss.length} characters, the sequence ${n} nt`;
  const close = { ')': '(', ']': '[', '}': '{', '>': '<' };
  const stacks = { '(': [], '[': [], '{': [], '<': [] };
  for (let k = 0; k < ss.length; k++) {
    const ch = ss[k];
    if (stacks[ch]) stacks[ch].push(k + 1);
    else if (close[ch]) {
      const i = stacks[close[ch]].pop();
      if (i === undefined) return `dot-bracket: "${ch}" at position ${k + 1} has no opening bracket`;
      if (k + 1 - i < 4) return `dot-bracket: the pair ${i}–${k + 1} closes a loop of fewer than 3 nt`;
    } else if (!'.-_,:'.includes(ch)) return `dot-bracket: "${ch}" at position ${k + 1} is not allowed (use . ( ) [ ] { } < >)`;
  }
  const open = Object.entries(stacks).filter(([, st]) => st.length).map(([o, st]) => [st[st.length - 1], o]).sort((a, b) => a[0] - b[0])[0];
  return open ? `dot-bracket: "${open[1]}" at position ${open[0]} is never closed` : null;
}

export function parseSequenceInput(text) {
  const lines = String(text || '').split(/\r?\n/);
  let name = '';
  const parts = [];
  for (const l of lines) {
    if (l.startsWith('>')) { if (!name) name = l.slice(1).trim().split(/\s+/)[0] || ''; continue; }
    if (/^\s*\/\/\s*$/.test(l)) continue;             // the end mark of a GenBank sequence block
    parts.push(l);
  }
  const raw = modseq.cleanLine(parts.join(''));
  const p = modseq.parse(raw);
  const bad = p.bad;
  // seq: parent bases (what is folded); raw: as written, with [modifications] (what is sent to the engine)
  return { name, seq: p.seq, raw, bad, mods: p.mods, problems: p.problems, unchecked: p.unknown };
}

export class Welcome {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    this.job = null;
    this.stopWatch = null;
    this.prefs = loadPref('predictForm', { models: 8, sampling: 'standard' });
    this._build();
  }

  // ---------------------------------------------------------------- layout
  // A compact start page: two cards side by side, prediction settings folded away.
  _build() {
    const app = this.app;
    this.viewerBtn = el('button', { class: 'btn', onclick: () => app.showPage('viewer') }, 'Open viewer', icon('next'));
    this.themeBtn = el('button', { class: 'btn ghost icon-only', title: 'Light / dark interface', 'aria-label': 'Toggle light or dark interface',
      onclick: () => app.setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark') });
    const top = el('header', { class: 'wl-top' },
      el('div', { class: 'brand' }, el('img', { class: 'logo', src: 'img/icon.svg?v=0.6', alt: '', draggable: 'false' }), el('span', { class: 'brand-name', text: 'RCV' }), el('span', { class: 'brand-sub', text: 'RNA Conformation Viewer' })),
      el('span', { class: 'grow' }), this.viewerBtn,
      el('button', { class: 'btn ghost', title: 'User guide (F1)', onclick: () => import('./dialogs.js').then(m => m.userGuideDialog(app)) }, icon('book'), el('span', { text: 'User guide' })),
      el('button', { class: 'btn ghost small lang-btn', title: 'Interface language', 'aria-label': 'Interface language', onclick: e => app.languageMenu(e.currentTarget) }),
      this.themeBtn);
    this.cards = el('div', { class: 'wl-cards' }, this._mainCard(), this._viewCard());
    this.progress = el('section', { class: 'wl-progress', hidden: true });
    this.batchBox = el('section', { class: 'wl-batch', hidden: true });
    this.recent = el('section', { class: 'wl-recent' });
    this.sys = el('div', { class: 'wl-sys' });
    const main = el('main', { class: 'wl-main' },
      el('div', { class: 'wl-hero' }, el('h1', { text: 'Predict, explore and present RNA conformations' }),
        el('p', { text: 'Paste a sequence to predict its conformations, or open structures you already have.' })),
      this.batchBox, this.cards, this.progress, this.recent);
    this.root.replaceChildren(top, main, el('footer', { class: 'wl-foot' }, this.sys));
    app._syncLanguage?.();
  }

  // The two ways to use RCV, side by side: 1 predict from sequence(s), 2 view structures of one RNA.
  _cardHead(n, title, sub, alt = false) {
    return el('div', { class: 'wl-card-head' }, el('span', { class: `wl-badge${alt ? ' alt' : ''}`, translate: 'no', text: n }),
      el('div', {}, el('h2', { text: title }), el('p', { text: sub })));
  }

  _mainCard() {
    const app = this.app, p = this.prefs;
    this.nameIn = el('input', { class: 'input', placeholder: 'from the FASTA header', 'aria-label': 'Project name', maxlength: 60 });
    this.seqIn = el('textarea', { class: 'input wl-seq', rows: 4, spellcheck: 'false', 'aria-label': 'RNA sequence',
      placeholder: 'Paste an RNA sequence or FASTA (several sequences = batch)…' });
    this.seqStatus = el('div', { class: 'wl-seqstatus' });
    this.modelsIn = el('input', { class: 'input', type: 'number', min: 1, max: 30, value: p.models ?? 8, 'aria-label': 'Conformations to keep' });
    this.samplingSeg = el('div', { class: 'seg', role: 'radiogroup', 'aria-label': 'Sampling depth' },
      ...Object.entries(SAMPLING).map(([k, v]) => el('button', { type: 'button', 'data-s': k, text: v[2] })));
    this.sampling = p.sampling || 'standard';
    this.samplingSeg.addEventListener('click', e => { const s = e.target.closest('button')?.dataset.s; if (s) { this.sampling = s; this._estimate(); } });
    this.coresIn = el('input', { type: 'range', min: 1, max: 8, value: 7, 'aria-label': 'CPU cores' });
    this.coresLabel = el('span', { class: 'muted' });
    this.ssIn = el('input', { class: 'input mono', placeholder: 'e.g. ((((((...((((......))))))))))', 'aria-label': 'Known secondary structure' });
    this.seedIn = el('input', { class: 'input', type: 'number', min: 0, step: 1, value: 42, 'aria-label': 'Random seed' });
    this.templatesIn = el('input', { type: 'checkbox', checked: p.templates !== false, 'aria-label': 'Use related experimental structures' });
    this.templatesIn.addEventListener('change', () => { this.prefs.templates = this.templatesIn.checked; savePref('predictForm', this.prefs); });
    // Z-RNA prediction (off by default: it adds trajectories); salt and Z-binding protein shown only when on
    this.zrnaIn = el('input', { type: 'checkbox', checked: !!p.zrna, 'aria-label': 'Predict left-handed Z-RNA helices' });
    this.zSaltIn = el('input', { class: 'input wl-zsalt', type: 'number', min: 5, max: 6000, step: 5, value: p.zrnaSaltMM || 150, 'aria-label': 'Na+ in mM for Z-RNA' });
    this.zBinderIn = el('input', { type: 'checkbox', checked: !!p.zrnaBinder, 'aria-label': 'Z-binding protein present' });
    this.zOpts = el('span', { class: 'wl-zopts' }, el('span', { class: 'muted', text: 'Na⁺ (mM)' }), this.zSaltIn,
      el('label', { class: 'check', title: 'ADAR1 Zα, ZBP1 or PKZ bound: they stabilise left-handed RNA strongly.' }, this.zBinderIn, el('span', { text: 'Z-binding protein (Zα)' })));
    const zSave = () => { this.prefs.zrna = this.zrnaIn.checked; this.prefs.zrnaSaltMM = +this.zSaltIn.value || 150; this.prefs.zrnaBinder = this.zBinderIn.checked; this.zOpts.hidden = !this.zrnaIn.checked; savePref('predictForm', this.prefs); this._estimate(); };
    for (const c of [this.zrnaIn, this.zSaltIn, this.zBinderIn]) c.addEventListener('change', zSave);
    this.zOpts.hidden = !this.zrnaIn.checked;
    // a recipe (0.7.5): its prediction settings fill the form, its steps run on the result
    this.recipe = null;
    this.recipeSel = el('select', { class: 'select', 'aria-label': 'Recipe' }, el('option', { value: '', text: 'None' }));
    this.recipeSel.addEventListener('change', async () => {
      const r = this.recipeSel.value ? await app.recipes.find(this.recipeSel.value) : null;
      this.setRecipe(r, { fromSelect: true });
    });
    // optional experimental / evolutionary information for the 2D structures (0.7.4)
    const filePick = (key, title) => {
      const input = el('input', { type: 'file', hidden: true, accept: '.txt,.shape,.dat,.csv,.tsv,.map,.react,.fa,.fasta,.afa,.aln,.sto,.stk,text/plain' });
      const btn = el('button', { class: 'btn small ghost wl-file', type: 'button', title });
      const show = () => { btn.replaceChildren(this[key] ? `${this[key].name}  ✕` : 'Choose file…'); btn.classList.toggle('on', !!this[key]); };
      input.addEventListener('change', async () => { const f = input.files[0]; input.value = ''; if (f) { this[key] = { name: f.name, text: (await f.text()).slice(0, 4e6) }; show(); } });
      btn.addEventListener('click', () => { if (this[key]) { this[key] = null; show(); } else input.click(); });
      show();
      return el('span', { class: 'wl-filepick' }, btn, input);
    };
    this.probingPick = filePick('probingFile', 'Reactivities (SHAPE, DMS): one value per nucleotide, or "position value" lines. They guide the 2D structures of the prediction as pseudo-energies. Click again to remove the file.');
    this.alignPick = filePick('alignmentFile', 'Aligned homologous sequences (FASTA with gaps, Stockholm or Clustal) that include this sequence. Base pairs supported by covariation become the first 2D structure. Click again to remove the file.');
    this.estimate = el('div', { class: 'wl-estimate' });
    this.predictBtn = el('button', { class: 'btn primary wl-go', onclick: () => this.startPrediction() }, icon('flask'), el('span', { text: 'Predict conformations' }));
    // engine: the built-in RCV engine (physics + templates) or a local AlphaFold 3 installation
    this.engine = p.engine === 'alphafold3' ? 'alphafold3' : 'rcv';
    this.engineSeg = el('div', { class: 'seg', role: 'radiogroup', 'aria-label': 'Prediction engine' },
      el('button', { type: 'button', 'data-e': 'rcv', title: 'Built-in engine: physics-based 3D sampling guided by known structures; runs on the CPU' }, 'RCV engine'),
      el('button', { type: 'button', 'data-e': 'alphafold3', title: 'Your own AlphaFold 3 installation (Linux, NVIDIA GPU)' }, 'AlphaFold 3 (local)'));
    this.engineSeg.addEventListener('click', e => { const v = e.target.closest('button')?.dataset.e; if (v) { this.engine = v; this.prefs.engine = v; savePref('predictForm', this.prefs); this._engineUi(); this._estimate(); } });
    this.af3Info = el('span', { class: 'muted af3-summary' });
    this.af3Row = el('div', { class: 'wl-af3', hidden: true }, icon('info'), this.af3Info, el('span', { class: 'grow' }),
      el('button', { class: 'btn small', type: 'button', onclick: () => af3SettingsDialog(app, sv => { this.af3 = sv; this._engineUi(); this._estimate(); }) }, 'AlphaFold 3 settings…'));
    this.batchList = el('div', { class: 'wl-batchlist', hidden: true });
    this.fastaIn = el('input', { type: 'file', accept: '.fa,.fasta,.fna,.txt,.seq', hidden: true });
    this.fastaIn.addEventListener('change', () => this.fastaIn.files[0] && this.loadFastaFile(this.fastaIn.files[0]));
    this.seqIn.addEventListener('dragover', e => { e.preventDefault(); this.seqIn.classList.add('over'); });
    this.seqIn.addEventListener('dragleave', () => this.seqIn.classList.remove('over'));
    this.seqIn.addEventListener('drop', async e => {
      e.preventDefault(); e.stopPropagation(); this.seqIn.classList.remove('over');
      const f = e.dataTransfer.files?.[0];
      if (f) await this.loadFastaFile(f);
    });
    for (const c of [this.seqIn, this.modelsIn, this.coresIn, this.ssIn]) c.addEventListener('input', () => this._estimate());
    this.seqIn.addEventListener('keydown', e => { if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) this.startPrediction(); });
    this.adv = el('details', { class: 'wl-adv' }, el('summary', { text: 'Advanced settings' }),
      el('div', { class: 'grid2' },
        el('div', { class: 'field' }, el('label', { text: 'Name' }), this.nameIn),
        el('div', { class: 'field' }, el('label', { text: 'Conformations to keep' }), this.modelsIn)),
      el('div', { class: 'grid2' },
        el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Sampling' }), this.samplingSeg),
        el('div', { class: 'field' }, el('span', { class: 'flabel' }, 'CPU cores ', this.coresLabel), this.coresIn)),
      el('div', { class: 'grid2' },
        el('div', { class: 'field' }, el('label', { text: 'Known secondary structure (dot-bracket)' }), this.ssIn),
        el('div', { class: 'field' }, el('label', { text: 'Random seed' }), this.seedIn)),
      el('div', { class: 'grid2' },
        el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Probing data (SHAPE, DMS), optional' }), this.probingPick),
        el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Alignment of homologs, optional' }), this.alignPick)),
      el('div', { class: 'grid2' },
        el('div', { class: 'field' }, el('label', { text: 'Recipe: its steps run on the result' }),
          el('div', { class: 'row', style: { gap: '6px' } }, this.recipeSel,
            el('button', { class: 'btn ghost icon-only', type: 'button', title: 'Recipes…', 'aria-label': 'Recipes…', onclick: () => import('./recipes.js').then(m => m.recipesDialog(app, this.recipe?.id)) }, icon('book')))),
        el('div')),
      el('label', { class: 'check wl-templates', title: 'The built-in library holds about 4000 RNA chains from the PDB (21 MB). Related structures (sequence identity, covering at least about half of the sequence) guide the 3D model; switch off for a physics-only prediction.' },
        this.templatesIn, el('span', { text: 'Use related experimental structures (templates) when the library has one' })),
      el('div', { class: 'wl-zrna' },
        el('label', { class: 'check', title: 'Adds a hypothesis with one helix built left-handed (Z-RNA: syn purines, zig-zag backbone) when alternating purine-pyrimidine pairs (CG, CA/UG repeats) could turn under these conditions: high salt, 8-substituted purines (m8G, 8-Br-G) or a Z-binding protein. About a quarter more run time.' },
          this.zrnaIn, el('span', { text: 'Predict Z-RNA (left-handed helices) when likely' })), this.zOpts),
      el('div', { class: 'wl-adv-links' },
        el('button', { class: 'linkbtn', type: 'button', onclick: () => this.useDemo() }, `Fill in the demo sequence (${DEMO.label})`)));
    return el('section', { class: 'wl-card wl-main-card' },
      this._cardHead('1', 'Predict RNA conformations from sequence(s)', 'Paste one sequence, or several in FASTA format. Runs on this computer.'),
      this.seqIn, el('div', { class: 'wl-statusrow' }, this.seqStatus, el('span', { class: 'grow' }), this.estimate), this.batchList, this.fastaIn,
      el('div', { class: 'wl-engine' }, el('span', { class: 'flabel', text: 'Engine' }), this.engineSeg), this.af3Row,
      el('div', { class: 'wl-actions' },
        this.predictBtn,
        el('button', { class: 'btn', onclick: () => this.chooseFasta() }, icon('open'), 'Load a FASTA file…'),
        this.modBtn = el('button', { class: 'btn ghost', type: 'button', onclick: () => this.editModifications(),
          title: 'Modified nucleotides (m6A, Ψ, inosine from A-to-I editing, 8-oxo-G, 2\'-O-methyl, LNA …): choose them on the sequence. They change base pairing in the 2D structures and are built into the 3D models. In the sequence they are written in brackets, e.g. GGC[m6A]CU.' },
          el('span', { class: 'wl-modmark', translate: 'no', text: 'Ψ' }), 'Modified nucleotides…')),
      this.adv);
  }

  _viewCard() {
    const app = this.app;
    this.fileIn = el('input', { type: 'file', multiple: true, accept: '.rcv,.pdb,.ent,.cif,.mmcif,.json,.zip', hidden: true });
    this.fileIn.addEventListener('change', () => this.fileIn.files.length && app.openLocalFiles([...this.fileIn.files]));
    const chooseFiles = async () => {
      if (app.native?.available()) { const paths = await app.native.chooseFiles(); if (paths.length) await app.openNativePaths(paths); }
      else this.fileIn.click();
    };
    // a big target for files; a drop on it reaches the window's drop handler, which opens the files
    const drop = el('button', { class: 'wl-drop', type: 'button', onclick: chooseFiles, title: 'Open an .rcv file, PDB/mmCIF files or a prediction folder' },
      icon('open'), el('b', { text: 'Open files…' }),
      el('span', { text: 'PDB/mmCIF, multi-model NMR or .rcv files, a prediction folder, or an AlphaFold Server result (.zip). You can also drop them here.' }));
    drop.addEventListener('dragover', () => drop.classList.add('over'));
    drop.addEventListener('dragleave', e => { if (!drop.contains(e.relatedTarget)) drop.classList.remove('over'); });
    drop.addEventListener('drop', () => drop.classList.remove('over'));
    this.demoBtn = el('button', { class: 'btn ghost', onclick: () => app.openDemo(), title: 'HIV-1 TAR NMR ensemble and two predictions' }, icon('cube'), 'Demo ensembles');
    return el('section', { class: 'wl-card wl-view-card' },
      this._cardHead('2', 'View experimental or predicted structures of one RNA', 'Several files or a multi-model file open together as one ensemble.', true),
      drop, this.fileIn,
      el('div', { class: 'wl-actions' },
        el('button', { class: 'btn', onclick: () => app.fetchDialog(), title: 'Download a structure from the Protein Data Bank' }, icon('download'), 'Fetch PDB…'),
        this.demoBtn),
      el('div', { class: 'wl-adv-links' },
        el('button', { class: 'linkbtn', type: 'button', onclick: () => app.openDialog('server') }, 'Browse folders…'),
        el('button', { class: 'linkbtn', type: 'button', onclick: () => app.openDialog('library') }, 'All prediction runs…'),
        el('button', { class: 'linkbtn', type: 'button', title: 'Predict on the AlphaFold Server, download the result .zip and open it here',
          onclick: () => import('./dialogs.js').then(m => m.openExternal(app, AF3_SERVER)) }, 'AlphaFold Server…')));
  }

  // Engine choice: AlphaFold 3 needs a local installation; sampling options apply to the RCV engine only.
  _engineUi() {
    const af3 = this.engine === 'alphafold3';
    for (const b of this.engineSeg.children) b.classList.toggle('on', b.dataset.e === this.engine);
    this.af3Row.hidden = !af3;
    this.af3Info.textContent = af3Summary(this.af3);
    this.af3Row.classList.toggle('warn', af3 && !af3Configured(this.af3));
    for (const c of [this.samplingSeg, this.modelsIn, this.ssIn]) c.closest('.field')?.classList.toggle('disabled', af3);
    this.templatesIn.closest('label')?.classList.toggle('disabled', af3);
    this.zrnaIn.closest('.wl-zrna')?.classList.toggle('disabled', af3);
  }

  // ---------------------------------------------------------------- state
  async refresh() {
    const app = this.app;
    const sys = app.info?.system;
    if (sys) {
      this.coresIn.max = sys.cpus;
      if (!this._coresSet) { this.coresIn.value = this.prefs.cores ? Math.min(this.prefs.cores, sys.cpus) : Math.max(1, sys.cpus - 1); this._coresSet = true; }
    }
    this.viewerBtn.hidden = !app.ensembles.size;
    this.themeBtn.replaceChildren(icon(document.documentElement.dataset.theme === 'dark' ? 'sun' : 'moon'));
    this.demoBtn.hidden = !app.info?.demo;
    if (this.af3 === undefined) { this.af3 = null; loadAf3Settings(app).then(s => { this.af3 = s; this._engineUi(); this._estimate(); }); }
    this._engineUi();
    if (!this._advInit) { this._advInit = true; if (this.prefs.advancedOpen) this.adv.open = true; this.adv.addEventListener('toggle', () => { this.prefs.advancedOpen = this.adv.open; savePref('predictForm', this.prefs); }); }
    const gpu = app.viewer?.gpuName();
    this.sys.replaceChildren(
      el('span', {}, el('span', { class: `dot ${app.api.online ? '' : 'off'}` }), app.api.online ? 'Runs locally on this computer' : 'Engine not reachable'),
      sys ? el('span', { text: `${sys.cpus} CPU cores` }) : null,
      gpu ? el('span', { text: `GPU: ${gpu}` }) : null,
      sys ? el('span', { text: sys.vienna ? `ViennaRNA ${sys.vienna}` : 'ViennaRNA missing — sequence prediction needs it' }) : null,
      sys ? el('span', { class: 'grow' }) : null,
      sys ? el('span', { class: 'muted', title: sys.projects_dir, text: `Projects: ${sys.projects_dir}` }) : null,
      el('span', { class: 'muted', text: `RCV ${app.info?.version || ''}` }));
    if (!modseq.getCatalogue() && app.api.online) modseq.loadCatalogue(app.api).then(c => { if (c) this._estimate(); });
    this._estimate();
    if (this.job) this._showProgress(true);
    await this._fillRecipes();
    await this._recent();
  }

  async _fillRecipes() {
    let list = [];
    try { list = await this.app.recipes.list(true); } catch { /* offline */ }
    const cur = this.recipe?.id || '';
    this.recipeSel.replaceChildren(el('option', { value: '', text: 'None' }),
      ...list.filter(r => !r.error).map(r => el('option', { value: r.id, text: r.name })));
    if (this.recipe && !list.some(r => r.id === cur)) this.recipeSel.append(el('option', { value: cur || '_unsaved', text: this.recipe.name }));
    this.recipeSel.value = cur || (this.recipe ? '_unsaved' : '');
  }

  // Chosen here or in Tools > Recipes ("Use for the next prediction"); null clears it.
  setRecipe(r, { fromSelect = false } = {}) {
    this.recipe = r || null;
    const p = r?.prediction || {};
    if (p.models) this.modelsIn.value = p.models;
    if (p.sampling) this.sampling = p.sampling;
    if (p.seed !== undefined && p.seed !== null) this.seedIn.value = p.seed;
    if (p.templates) this.templatesIn.checked = p.templates !== 'off';
    if (!fromSelect) { this.adv.open = true; this._fillRecipes(); }
    this._estimate();
  }

  _estimate() {
    const records = parseFastaRecords(this.seqIn.value);
    this.batchList.hidden = records.length < 2;
    if (records.length >= 2) return this._estimateBatch(records);
    const af3 = this.engine === 'alphafold3';
    this.predictBtn.lastChild.textContent = af3 ? 'Predict with AlphaFold 3' : this.seqIn.value.trim() ? 'Predict conformations' : 'Predict the demo sequence';
    const { seq, bad, name, mods, problems } = parseSequenceInput(this.seqIn.value);
    if (name && !this.nameIn.value) this.nameIn.placeholder = name;
    for (const b of this.samplingSeg.children) b.classList.toggle('on', b.dataset.s === this.sampling);
    const cores = +this.coresIn.value;
    this.coresLabel.textContent = `· ${cores} of ${this.coresIn.max}`;
    const models = Math.max(1, +this.modelsIn.value || 8);
    const [per, low] = SAMPLING[this.sampling];
    const traj = Math.max(low, per * models);
    const n = seq.length;
    const ss = this.ssIn.value.trim();
    let msg = '', cls = '';
    const more = [];                                   // further parts of the status line
    if (!n) msg = `No sequence yet: Predict runs the demo sequence (${DEMO.label}). T is read as U; position numbers and spaces are ignored.`;
    else if (problems.length) { msg = `${problems[0][0].toUpperCase()}${problems[0].slice(1)}. Modified nucleotides are written in brackets, e.g. GGC[m6A]CU.`; cls = 'bad'; }
    else if (bad.length) { msg = `Unsupported letters: ${bad.join(' ')}. Only A, C, G, U (T) are modelled.`; cls = 'bad'; }
    else if (n < 8) { msg = `${n} nt — at least 8 nt are needed.`; cls = 'bad'; }
    else if (n > 400) { msg = `${n} nt — above the 400-nt limit of the engine.`; cls = 'bad'; }
    else {
      const comp = ['A', 'C', 'G', 'U'].map(b => `${b} ${(seq.match(new RegExp(b, 'g')) || []).length}`).join(' · ');
      msg = `${n} nt · ${comp}`;
      if (n > 200) more.push('above the tested 200-nt range');
      if (mods.length) more.push(`${mods.length} modified: ${modseq.summary(mods)}`);
      cls = n > 200 ? 'warn' : 'ok';
    }
    const ssBad = ss && n && !af3 ? dotBracketProblem(ss, n) : null;
    if (ssBad) { more.push(ssBad); cls = 'bad'; }
    setParts(this.seqStatus, [msg, ...more]);
    this.seqStatus.className = `wl-seqstatus ${cls}`;
    const empty = !this.seqIn.value.trim();
    const ok = (empty || (n >= 8 && n <= 400 && !bad.length && !problems.length && !ssBad)) && this.app.api.online && !this.job;
    this.modBtn?.classList.toggle('on', mods.length > 0);
    this.predictBtn.disabled = !ok;
    const nn = empty ? DEMO.sequence.length : n;
    if (af3) {
      const seeds = this.af3?.seeds || 1;
      this.estimate.textContent = `AlphaFold 3 · ${seeds * 5} models (${seeds} seed${seeds > 1 ? 's' : ''} × 5) · on the GPU`;
      return;
    }
    const zx = this.zrnaIn?.checked ? Math.max(2, Math.round(0.25 * traj)) : 0;      // Z-RNA: up to this many extra
    this.estimate.textContent = nn >= 8 && !bad.length ? `${traj}${zx ? ` (+${zx} Z-RNA)` : ''} trajectories · about ${fmtDuration(trajSeconds(nn) * Math.ceil((traj + zx) / cores) + 15 + (zx ? 10 : 0))} on ${cores} core${cores > 1 ? 's' : ''}` : '';
  }

  _estimateBatch(records) {
    for (const b of this.samplingSeg.children) b.classList.toggle('on', b.dataset.s === this.sampling);
    const cores = +this.coresIn.value;
    this.coresLabel.textContent = `· ${cores} of ${this.coresIn.max}`;
    const models = Math.max(1, +this.modelsIn.value || 8);
    const [per, low] = SAMPLING[this.sampling];
    const traj = Math.max(low, per * models);
    const bad = records.filter(r => r.problem);
    const lens = records.map(r => r.length);
    const zx = this.zrnaIn?.checked ? Math.max(2, Math.round(0.25 * traj)) : 0;      // Z-RNA: up to this many extra
    const total = records.reduce((s, r) => s + trajSeconds(r.length) * Math.ceil((traj + zx) / cores) + 25, 0);
    const unused = this.probingFile || this.alignmentFile ? ' · probing data and alignments apply to one sequence: not used in a batch' : '';
    this.seqStatus.textContent = bad.length ? `${records.length} sequences · ${bad.length} cannot be predicted (see the list)` : `${records.length} sequences (batch) · ${Math.min(...lens)}–${Math.max(...lens)} nt · one .rcv file per sequence${unused}`;
    this.seqStatus.className = `wl-seqstatus ${bad.length ? 'bad' : 'ok'}`;
    this.batchList.replaceChildren(...records.slice(0, 200).map((r, i) => el('div', { class: `wl-bitem ${r.problem ? 'bad' : ''}` },
      el('span', { class: 'n', text: String(i + 1) }), el('b', { text: r.name }), el('span', { class: 'muted', text: `${r.length} nt${r.mods?.length ? ` · ${r.mods.length} modified` : ''}${r.structure ? ' · with 2D' : ''}` }),
      el('span', { class: 'grow' }), r.problem ? el('span', { class: 'state-bad', text: r.problem }) : el('span', { class: 'state-ok' }, icon('check')))));
    this.predictBtn.lastChild.textContent = `Predict ${records.length} sequences`;
    this.predictBtn.disabled = !!bad.length || records.length > 100 || !this.app.api.online || !!this.job;
    this.estimate.textContent = `${records.length} × ${traj}${zx ? ` (+${zx} Z-RNA)` : ''} trajectories · about ${fmtDuration(total)} on ${cores} core${cores > 1 ? 's' : ''}` + (records.length > 100 ? ' · at most 100 sequences per batch' : '');
  }

  // ---------------------------------------------------------------- input helpers
  // Modified nucleotides: a sub-window to choose them on the sequence; the result is written back as [code] tokens.
  async editModifications() {
    const recs = parseFastaRecords(this.seqIn.value);
    if (recs.length > 1) { toast('Choose modifications for one sequence at a time. In a batch, write them in brackets in each sequence, e.g. GGC[m6A]CU.', 'info', 'Modified nucleotides', 7000); return; }
    const m = await import('./modifications.js');
    await m.modificationsDialog(this.app, {
      text: this.seqIn.value.trim() ? this.seqIn.value : `>${DEMO.name}\n${DEMO.sequence}`,
      onApply: text => { this.seqIn.value = text; this._estimate(); },
    });
  }

  useDemo() {
    this.seqIn.value = `>${DEMO.name}\n${DEMO.sequence}`;
    this.nameIn.value = '';
    this._estimate();
  }

  async chooseFasta() {
    const app = this.app;
    if (app.native?.available()) {
      const [path] = await app.native.chooseFiles('fasta');
      if (!path) return;
      try { this._setFasta(await app.native.readText(path), path.split(/[\\/]/).pop()); } catch (err) { toast(err.message, 'error', 'FASTA'); }
    } else this.fastaIn.click();
  }

  async loadFastaFile(file) {
    if (file.size > 20 * 1024 * 1024) { toast('This file is larger than 20 MB.', 'error', 'FASTA'); return; }
    this._setFasta(await readFileText(file), file.name);
  }

  _setFasta(text, name) {
    this.seqIn.value = text.trim();
    const recs = parseFastaRecords(text);
    this._estimate();
    toast(recs.length > 1 ? `${recs.length} sequences read from ${name}: they will be predicted one after another.` : `Sequence read from ${name}.`, 'info', 'FASTA');
    this.seqIn.focus();
  }

  async _recent() {
    const app = this.app;
    this.recent.replaceChildren();
    if (!app.api.online) return;
    let lib;
    try { lib = await app.api.get('/api/library'); } catch { return; }
    const items = lib.packages.slice(0, 8);
    if (!items.length) return;
    const list = el('div', { class: 'wl-recent-list' });
    for (const p of items) {
      const date = p.created ? p.created.replace('T', ' ').slice(0, 16) : '';
      list.append(el('button', { class: 'wl-recent-item', title: p.path, onclick: () => app.openServerPaths([p.path]) },
        icon('cube'), el('div', { class: 'meta' }, el('b', { text: p.name }), el('span', {}, `${p.length ?? '?'} nt · ${p.n_conformations} conformation${p.n_conformations === 1 ? '' : 's'}`, date ? ' · ' : '', date ? el('span', { class: 'nowrap', text: date }) : ''))));
    }
    this.recent.append(el('div', { class: 'section-title' }, 'Recent predictions', el('span', { class: 'grow' }),
      el('button', { class: 'linkbtn', onclick: () => app.openLibrary() }, 'Project library…')), list);
  }

  // ---------------------------------------------------------------- prediction
  // The random seed as entered: 0 is a seed like any other; an empty or unusable field gives the default 42.
  _seed() {
    const v = String(this.seedIn.value).trim(), n = Math.trunc(Number(v));
    return v !== '' && Number.isFinite(n) && n >= 0 && n <= 2147483647 ? n : 42;
  }

  _zrnaParams() {
    return { zrna: 'auto', zrna_salt: (+this.zSaltIn.value || 150) / 1000, zrna_binder: this.zBinderIn.checked };
  }

  async startPrediction() {
    const app = this.app;
    if (!this.seqIn.value.trim()) {
      this.useDemo();
      toast(`No sequence was entered, so the demo sequence (${DEMO.label}) is predicted. Paste your own sequence or a FASTA file to predict it instead.`, 'info', 'Demo sequence', 7000);
    }
    const records = parseFastaRecords(this.seqIn.value);
    const af3 = this.engine === 'alphafold3';
    if (af3 && !af3Configured(this.af3)) { af3SettingsDialog(app, sv => { this.af3 = sv; this._engineUi(); this._estimate(); }); return; }
    if (af3 && records.length > 1) { toast('AlphaFold 3 runs one sequence at a time here. Predict the sequences one by one, or use the RCV engine for a batch.', 'error', 'AlphaFold 3'); return; }
    if (records.length > 1) return this.startBatch(records);
    const { raw, name: fastaName } = parseSequenceInput(this.seqIn.value);
    const name = this.nameIn.value.trim() || fastaName || 'rna';
    const params = { name: slug(name), sequence: raw, structure: this.ssIn.value.trim(), models: +this.modelsIn.value || 8,
      sampling: this.sampling, cores: +this.coresIn.value, seed: this._seed(), templates: this.templatesIn.checked ? 'auto' : 'off',
      ...(this.zrnaIn.checked && !af3 ? this._zrnaParams() : {}),
      ...(af3 ? { engine: 'alphafold3' } : { probing_text: this.probingFile?.text || '', alignment_text: this.alignmentFile?.text || '' }) };
    this.prefs = { ...this.prefs, models: params.models, sampling: params.sampling, cores: params.cores };
    savePref('predictForm', this.prefs);
    try {
      this.job = await app.api.post('/api/predict', params);
    } catch (err) { toast(err.message, 'error', 'Cannot start the prediction'); return; }
    this._showProgress(true);
    this._render(this.job);
    this.stopWatch = watchJob(app.api, this.job.id, {
      onUpdate: j => { this.job = j; this._render(j); },
      onDone: async j => {
        this.job = j; this._render(j);
        try { await app.loadServerDocument(j.document.id); } catch (err) { toast(err.message, 'error'); }
        this.job = null;
        this._showProgress(false);
        app.showPage('viewer');
        if (this.recipe?.steps?.length) { const r = this.recipe; setTimeout(() => app.recipes.run(r), 300); }   // the chosen recipe's steps
        toast(j.result?.engine === 'alphafold3' ? `${j.document.n_conformations} AlphaFold 3 models of ${j.document.name}, best ranking score first. Colour by Confidence to see the pLDDT.`
          : `${j.document.n_conformations} conformation${j.document.n_conformations === 1 ? '' : 's'} of ${j.document.name}. Drag a residue with the Pull tool to switch between them.`, 'ok', 'Prediction finished', 7000);
      },
      onFail: j => { this.job = j; this._render(j); },
    });
  }

  async startBatch(records) {
    const app = this.app;
    const params = { records: records.map(r => ({ name: r.name, sequence: r.sequence, structure: r.structure })), models: +this.modelsIn.value || 8,
      sampling: this.sampling, cores: +this.coresIn.value, seed: this._seed(), templates: this.templatesIn.checked ? 'auto' : 'off',
      ...(this.zrnaIn.checked ? this._zrnaParams() : {}) };
    this.prefs = { ...this.prefs, models: params.models, sampling: params.sampling, cores: params.cores };
    savePref('predictForm', this.prefs);
    try {
      this.job = await app.api.post('/api/batch', params);
    } catch (err) { toast(err.message, 'error', 'Cannot start the batch'); return; }
    this.batch = null;
    this._showProgress(true);
    this._render(this.job);
    this.stopWatch = watchJob(app.api, this.job.id, {
      onUpdate: j => { this.job = j; this._render(j); },
      onDone: j => {
        this.job = null; this.batch = j; this._showProgress(false); this._estimate(); this._batchResults();
        const items = j.result?.items || [];
        toast(`${items.filter(it => it.status === 'done').length} of ${items.length} sequences predicted. Choose the one to open${app.page === 'welcome' ? ' below' : ' on the start page (Ctrl+H)'}.`, 'ok', 'Batch finished', 8000);
      },
      onFail: j => {
        this.job = j; this._render(j);
        if ((j.result?.items || []).some(it => it.status === 'done')) { this.batch = j; this._batchResults(); }
      },
    });
  }

  // Batch finished: the user chooses which prediction to open (the others stay in the batch folder).
  _batchResults() {
    const app = this.app, j = this.batch;
    if (!j) { this.batchBox.hidden = true; return; }
    const items = j.result?.items || [];
    let chosen = null;
    const openBtn = el('button', { class: 'btn primary', disabled: true, onclick: async () => {
      const it = items[chosen];
      if (!it?.rcv) return;
      await app.openServerPaths([it.rcv]);
    } }, icon('cube'), 'Open selected in the viewer');
    const rows = items.map((it, i) => {
      const ok = it.status === 'done' && it.rcv;
      const radio = el('input', { type: 'radio', name: 'batch-pick', disabled: ok ? null : true, 'aria-label': `Choose ${it.name}` });
      radio.addEventListener('change', () => { chosen = i; openBtn.disabled = false; for (const r of tb.children) r.classList.toggle('on', r === tr); });
      const tr = el('tr', { class: ok ? '' : 'muted', ondblclick: () => { if (ok) { radio.checked = true; chosen = i; openBtn.click(); } } },
        el('td', {}, radio), el('td', {}, el('b', { text: it.name })), el('td', { class: 'num', text: String(it.length) }),
        el('td', { class: 'num', text: ok ? String(it.n_conformations ?? '–') : '–' }), el('td', { class: 'num', text: it.best_score !== null && it.best_score !== undefined ? fmt(it.best_score, 1) : '–' }),
        el('td', { class: 'num', text: it.dg_mfe !== null && it.dg_mfe !== undefined ? fmt(it.dg_mfe, 1) : '–' }), el('td', { class: 'num', text: it.distinct_2d ?? '–' }),
        el('td', { class: 'num', text: it.elapsed_s ? `${Math.round(it.elapsed_s)} s` : '–' }),
        el('td', {}, ok ? el('span', { class: 'state-ok', text: 'done' })
          : it.error ? el('span', { class: 'state-bad' }, el('span', { text: 'failed' }), ': ', el('span', { text: it.error.slice(0, 60), translate: 'no' }))
            : el('span', { class: 'state-bad', text: it.status })));
      tr.addEventListener('click', e => { if (ok && e.target !== radio) { radio.checked = true; radio.dispatchEvent(new Event('change')); } });
      return tr;
    });
    const tb = el('tbody', {}, ...rows);
    const done = items.filter(it => it.status === 'done').length;
    this.batchBox.hidden = false;
    this.batchBox.replaceChildren(el('div', { class: 'wl-pg-card' },
      el('div', { class: 'wl-pg-head' }, el('h2', { text: `Batch finished · ${done} of ${items.length} sequences predicted` }), el('span', { class: 'grow' }),
        el('button', { class: 'btn ghost', onclick: () => { this.batch = null; this.batchBox.hidden = true; } }, 'Dismiss')),
      el('p', { class: 'desc', text: 'Choose the prediction to show in the viewer. Every sequence was saved as an .rcv file in the batch folder, so the others can be opened later (File ▸ Open, or Recent predictions).' }),
      el('div', { class: 'wl-table' }, el('table', { class: 'table' },
        el('thead', {}, el('tr', {}, el('th', { text: '' }), el('th', { text: 'Sequence' }), el('th', { class: 'num', text: 'nt' }), el('th', { class: 'num', text: 'Conformations' }),
          el('th', { class: 'num', text: 'Best score' }), el('th', { class: 'num', text: 'ΔG MFE' }), el('th', { class: 'num', text: '2D states' }), el('th', { class: 'num', text: 'Time' }), el('th', { text: 'Status' }))), tb)),
      el('div', { class: 'wl-row end' },
        el('span', { class: 'muted grow', text: j.result?.out ? `Folder: ${j.result.out}` : '' }),
        app.native?.available() && j.result?.out ? el('button', { class: 'btn', onclick: () => app.native.reveal(j.result.out) }, icon('folder'), 'Show folder') : null,
        openBtn)));
    this.batchBox.scrollIntoView({ block: 'nearest' });
  }

  _showProgress(on) {
    this.progress.hidden = !on;
    this.cards.hidden = on;
    this.recent.hidden = on;
  }

  _render(j) {
    const steps = el('ol', { class: 'wl-steps' }, ...(j.stages || []).map(s => el('li', { class: `st-${s.state}` },
      el('span', { class: 'st-icon' }, s.state === 'done' ? icon('check') : s.state === 'active' ? el('span', { class: 'spinner small' }) : s.state === 'error' ? icon('warn') : el('span', { class: 'st-dot' })),
      el('div', {}, el('b', { text: s.label }), s.detail ? el('span', { text: s.detail }) : null))));
    const pct = Math.round((j.progress || 0) * 100);
    const finished = ['done', 'failed', 'cancelled'].includes(j.status);
    const cancelling = !finished && (j.cancelling || this._cancelling === j.id);
    const meta = j.lost && !finished ? 'The connection to the RCV engine was lost; the last known state is shown. RCV keeps trying to reconnect. If RCV itself was closed, the run has stopped: start RCV again.'
      : cancelling ? 'Cancelling… stopping the computation'
      : j.status === 'running' || j.status === 'queued'
        ? `${pct}% · ${j.eta_s ? 'about ' + fmtDuration(j.eta_s) + ' left' : 'estimating time'} · ${fmtDuration(j.elapsed_s)} elapsed`
        : j.status === 'done' ? `Finished in ${fmtDuration(j.elapsed_s)} — opening the viewer…` : j.status === 'cancelled' ? 'Cancelled' : null;
    const actions = el('div', { class: 'wl-row end' });
    const back = () => { this.stopWatch?.(); this.job = null; this._cancelling = null; this._showProgress(false); this._estimate(); };
    if (!finished) {
      if (j.lost) actions.append(el('button', { class: 'btn', onclick: back, title: 'Leave this page. A run that is still going on continues in the background.' }, 'Back to the start page'));
      actions.append(el('button', { class: 'btn', disabled: cancelling || j.lost ? true : null, onclick: async () => {
        if (await cancelJob(this.app, this.job || j)) { this._cancelling = j.id; if (this.job) this._render(this.job); }
      } }, cancelling ? 'Cancelling…' : 'Cancel'));
    } else if (j.status !== 'done') {
      actions.append(el('button', { class: 'btn primary', onclick: back }, 'Back'));
    }
    // The card is built once per run and its parts are refreshed in place. The notes list stays in the page: it is
    // refilled only when a note was added (a selection in it survives the refresh every second), and it follows the
    // newest note only while the user has not scrolled up to read.
    let c = this._card;
    if (!c || c.job !== j.id || !c.root.isConnected) {
      c = this._card = { job: j.id, sig: '', head: el('div', { class: 'wl-pg-head' }), steps: el('div'), bar: el('div'), meta: el('div'),
        notes: el('ul', { class: 'wl-notes', 'aria-live': 'polite' }), actions: el('div') };
      c.root = el('div', { class: 'wl-pg-card' }, c.head, c.steps, c.bar, c.meta,
        el('div', { class: 'section-title', style: { marginTop: '14px' } }, 'Notes'), c.notes, c.actions);
      this.progress.replaceChildren(c.root);
    }
    c.head.replaceChildren(el('h2', { text: j.kind === 'batch' ? `Batch prediction · ${j.title || ''}` : `Predicting conformations · ${j.title || ''}` }),
      el('span', { class: 'muted', text: j.kind === 'batch' ? `${(j.result?.items || []).length} sequences, one after another, on ${j.result?.cores} CPU cores` : j.result?.engine === 'alphafold3' ? `AlphaFold 3 · ${j.result.seeds} seed(s) × 5 samples` : j.result?.cores ? `${j.result.trajectories} trajectories on ${j.result.cores} CPU cores` : '' }));
    c.steps.replaceChildren(steps);
    c.bar.replaceChildren(el('div', { class: `progress big ${j.status === 'failed' || (j.lost && !finished) ? 'bad' : ''}` }, el('i', { style: { width: `${pct}%` } })));
    // a failure: "Failed:" and the reason are two texts, so the first is translated also when the reason is not
    c.meta.replaceChildren(meta !== null ? el('div', { class: `wl-pg-meta${j.lost && !finished ? ' lost' : ''}`, text: meta })
      : el('div', { class: 'wl-pg-meta' }, el('span', { text: 'Failed:' }), ' ', el('span', { text: j.error || '' })));
    c.actions.replaceChildren(actions);
    const sig = `${(j.notes || []).length}|${(j.notes || []).at(-1)?.t ?? ''}|${(j.notes || []).at(-1)?.text ?? ''}`;
    if (sig !== c.sig && j.notes) {
      const follow = c.notes.scrollHeight - c.notes.scrollTop - c.notes.clientHeight < 28;
      const top = c.notes.scrollTop;
      // lead: a name in front of the text (the sequence of a batch); keep: data inside the text (a path, a title).
      // Neither is translated (i18n.js, data-i18n-keep).
      c.notes.replaceChildren(...j.notes.map(n => {
        const lead = n.lead && n.text.startsWith(n.lead) ? n.lead : '';
        return el('li', { class: `nt-${n.level}` },
          el('span', { class: 'nt-t', text: `${Math.floor(n.t / 60)}:${String(Math.floor(n.t % 60)).padStart(2, '0')}` }), icon(n.level === 'error' ? 'warn' : n.level === 'warning' ? 'warn' : 'info'),
          el('span', {}, lead ? el('span', { translate: 'no', text: lead }) : null,
            el('span', { text: n.text.slice(lead.length), 'data-i18n-keep': n.keep?.length ? JSON.stringify(n.keep) : null })));
      }));
      c.notes.scrollTop = follow ? c.notes.scrollHeight : top;
      c.sig = sig;
    }
    if (finished) this._cancelling = null;
    this._estimate();
  }
}
