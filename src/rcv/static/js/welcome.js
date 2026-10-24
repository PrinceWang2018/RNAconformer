// Start page: predict an ensemble from sequence, or open structures; progress of local predictions.
import { el, slug, loadPref, savePref, readFileText, fmt } from './util.js';
import { icon } from './icons.js';
import { toast } from './ui.js';
import { watchJob, fmtDuration } from './jobs.js';

const SAMPLING = { fast: [2, 6, 'Fast'], standard: [4, 8, 'Standard'], thorough: [8, 16, 'Thorough'] };
const trajSeconds = n => 45 + n + 0.004 * n * n;          // same rough model as the job runner

export const DEMO = { name: 'HIV1_TAR_demo', sequence: 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC', label: 'HIV-1 TAR, 29 nt' };

// FASTA (several records) or one bare sequence -> [{ name, sequence, structure }] (mirrors rcv.jobs.parse_fasta).
export function parseFastaRecords(text) {
  const recs = [];
  let cur = null;
  for (const raw of String(text || '').split(/\r?\n/)) {
    const line = raw.trim();
    if (!line || line.startsWith(';') || line.startsWith('#')) continue;
    if (line.startsWith('>')) { cur = { name: line.slice(1).trim().split(/\s+/)[0] || '', sequence: '', structure: '' }; recs.push(cur); continue; }
    if (!cur) { cur = { name: '', sequence: '', structure: '' }; recs.push(cur); }
    const tok = line.split(/\s+/)[0];
    if (/^[().\[\]{}<>]+$/.test(tok) && cur.sequence) cur.structure += tok;
    else cur.sequence += line.replace(/\s+/g, '').toUpperCase().replace(/T/g, 'U');
  }
  const seen = new Map();
  return recs.filter(r => r.sequence).map((r, k) => {
    const base = slug(r.name || `seq${k + 1}`);
    const n = seen.get(base) || 0;
    seen.set(base, n + 1);
    const bad = [...new Set(r.sequence.replace(/[ACGU]/g, ''))];
    const problem = bad.length ? `letters ${bad.join('')}` : r.sequence.length < 8 ? 'shorter than 8 nt' : r.sequence.length > 400 ? 'longer than 400 nt'
      : r.structure && r.structure.length !== r.sequence.length ? 'dot-bracket length differs' : null;
    return { ...r, name: n ? `${base}_${n + 1}` : base, problem };
  });
}

export function parseSequenceInput(text) {
  const lines = String(text || '').split(/\r?\n/);
  let name = '';
  const parts = [];
  for (const l of lines) {
    if (l.startsWith('>')) { if (!name) name = l.slice(1).trim().split(/\s+/)[0] || ''; continue; }
    parts.push(l);
  }
  const raw = parts.join('').replace(/\s+/g, '').toUpperCase().replace(/T/g, 'U');
  const bad = [...new Set(raw.replace(/[ACGU]/g, ''))];
  return { name, seq: raw, bad };
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
  // A compact start page: one box for a sequence (or FASTA), four buttons, settings folded away.
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
    this.cards = el('div', { class: 'wl-cards' }, this._mainCard());
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

  _mainCard() {
    const app = this.app, p = this.prefs;
    this.nameIn = el('input', { class: 'input', placeholder: 'from the FASTA header', 'aria-label': 'Project name', maxlength: 60 });
    this.seqIn = el('textarea', { class: 'input wl-seq', rows: 4, spellcheck: 'false', 'aria-label': 'RNA sequence',
      placeholder: 'Paste an RNA sequence or FASTA (several sequences = batch)…  Empty: the demo sequence is used.' });
    this.seqStatus = el('div', { class: 'wl-seqstatus' });
    this.modelsIn = el('input', { class: 'input', type: 'number', min: 1, max: 30, value: p.models ?? 8, 'aria-label': 'Conformations to keep' });
    this.samplingSeg = el('div', { class: 'seg', role: 'radiogroup', 'aria-label': 'Sampling depth' },
      ...Object.entries(SAMPLING).map(([k, v]) => el('button', { type: 'button', 'data-s': k, text: v[2] })));
    this.sampling = p.sampling || 'standard';
    this.samplingSeg.addEventListener('click', e => { const s = e.target.closest('button')?.dataset.s; if (s) { this.sampling = s; this._estimate(); } });
    this.coresIn = el('input', { type: 'range', min: 1, max: 8, value: 7, 'aria-label': 'CPU cores' });
    this.coresLabel = el('span', { class: 'muted' });
    this.ssIn = el('input', { class: 'input mono', placeholder: 'e.g. ((((((...((((......))))))))))', 'aria-label': 'Known secondary structure' });
    this.seedIn = el('input', { class: 'input', type: 'number', value: 42, 'aria-label': 'Random seed' });
    this.estimate = el('div', { class: 'wl-estimate' });
    this.predictBtn = el('button', { class: 'btn primary wl-go', onclick: () => this.startPrediction() }, icon('flask'), el('span', { text: 'Predict conformations' }));
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
    // opening structures
    this.fileIn = el('input', { type: 'file', multiple: true, accept: '.rcv,.pdb,.ent,.cif,.mmcif,.json', hidden: true });
    this.fileIn.addEventListener('change', () => this.fileIn.files.length && app.openLocalFiles([...this.fileIn.files]));
    const chooseFiles = async () => {
      if (app.native?.available()) { const paths = await app.native.chooseFiles(); if (paths.length) await app.openNativePaths(paths); }
      else this.fileIn.click();
    };
    this.demoBtn = el('button', { class: 'btn ghost', onclick: () => app.openDemo(), title: 'HIV-1 TAR NMR ensemble and two predictions' }, icon('cube'), 'Demo ensembles');
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
      el('div', { class: 'wl-adv-links' },
        el('button', { class: 'linkbtn', type: 'button', onclick: () => this.chooseFasta() }, 'Load a FASTA file…'),
        el('button', { class: 'linkbtn', type: 'button', onclick: () => this.useDemo() }, `Fill in the demo sequence (${DEMO.label})`),
        el('button', { class: 'linkbtn', type: 'button', onclick: () => app.openDialog('server') }, 'Browse folders…'),
        el('button', { class: 'linkbtn', type: 'button', onclick: () => app.openDialog('library') }, 'All prediction runs…')));
    return el('section', { class: 'wl-card wl-main-card' },
      this.seqIn, el('div', { class: 'wl-statusrow' }, this.seqStatus, el('span', { class: 'grow' }), this.estimate), this.batchList, this.fastaIn, this.fileIn,
      el('div', { class: 'wl-actions' },
        this.predictBtn,
        el('button', { class: 'btn', onclick: chooseFiles, title: 'Open an .rcv file, PDB/mmCIF files or a prediction folder' }, icon('open'), 'Open files…'),
        el('button', { class: 'btn', onclick: () => app.fetchDialog(), title: 'Download a structure from the Protein Data Bank' }, icon('download'), 'Fetch PDB…'),
        this.demoBtn),
      this.adv,
      el('div', { class: 'wl-hint', text: 'Tip: drop .rcv, PDB/mmCIF or FASTA files anywhere in this window.' }));
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
    this._estimate();
    if (this.job) this._showProgress(true);
    await this._recent();
  }

  _estimate() {
    const records = parseFastaRecords(this.seqIn.value);
    this.batchList.hidden = records.length < 2;
    if (records.length >= 2) return this._estimateBatch(records);
    this.predictBtn.lastChild.textContent = this.seqIn.value.trim() ? 'Predict conformations' : 'Predict the demo sequence';
    const { seq, bad, name } = parseSequenceInput(this.seqIn.value);
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
    if (!n) msg = `Paste a sequence or load a FASTA file (several sequences run as a batch); T is read as U. Empty: the demo sequence (${DEMO.label}) is used.`;
    else if (bad.length) { msg = `Unsupported letters: ${bad.join(' ')}. Only A, C, G, U (T) are modelled.`; cls = 'bad'; }
    else if (n < 8) { msg = `${n} nt — at least 8 nt are needed.`; cls = 'bad'; }
    else if (n > 400) { msg = `${n} nt — above the 400-nt limit of the engine.`; cls = 'bad'; }
    else {
      const comp = ['A', 'C', 'G', 'U'].map(b => `${b} ${(seq.match(new RegExp(b, 'g')) || []).length}`).join(' · ');
      msg = `${n} nt · ${comp}${n > 200 ? ' · above the tested 200-nt range' : ''}`;
      cls = n > 200 ? 'warn' : 'ok';
    }
    if (ss && n && ss.length !== n) { msg += ` · dot-bracket has ${ss.length} characters`; cls = 'bad'; }
    this.seqStatus.textContent = msg;
    this.seqStatus.className = `wl-seqstatus ${cls}`;
    const empty = !this.seqIn.value.trim();
    const ok = (empty || (n >= 8 && n <= 400 && !bad.length && !(ss && ss.length !== n))) && this.app.api.online && !this.job;
    this.predictBtn.disabled = !ok;
    const nn = empty ? DEMO.sequence.length : n;
    this.estimate.textContent = nn >= 8 && !bad.length ? `${traj} trajectories · about ${fmtDuration(trajSeconds(nn) * Math.ceil(traj / cores) + 15)} on ${cores} core${cores > 1 ? 's' : ''}` : '';
  }

  _estimateBatch(records) {
    for (const b of this.samplingSeg.children) b.classList.toggle('on', b.dataset.s === this.sampling);
    const cores = +this.coresIn.value;
    this.coresLabel.textContent = `· ${cores} of ${this.coresIn.max}`;
    const models = Math.max(1, +this.modelsIn.value || 8);
    const [per, low] = SAMPLING[this.sampling];
    const traj = Math.max(low, per * models);
    const bad = records.filter(r => r.problem);
    const lens = records.map(r => r.sequence.length);
    const total = records.reduce((s, r) => s + trajSeconds(r.sequence.length) * Math.ceil(traj / cores) + 25, 0);
    this.seqStatus.textContent = bad.length ? `${records.length} sequences · ${bad.length} cannot be predicted (see the list)` : `${records.length} sequences (batch) · ${Math.min(...lens)}–${Math.max(...lens)} nt · one .rcv file per sequence`;
    this.seqStatus.className = `wl-seqstatus ${bad.length ? 'bad' : 'ok'}`;
    this.batchList.replaceChildren(...records.slice(0, 200).map((r, i) => el('div', { class: `wl-bitem ${r.problem ? 'bad' : ''}` },
      el('span', { class: 'n', text: String(i + 1) }), el('b', { text: r.name }), el('span', { class: 'muted', text: `${r.sequence.length} nt${r.structure ? ' · with 2D' : ''}` }),
      el('span', { class: 'grow' }), r.problem ? el('span', { class: 'state-bad', text: r.problem }) : el('span', { class: 'state-ok' }, icon('check')))));
    this.predictBtn.lastChild.textContent = `Predict ${records.length} sequences`;
    this.predictBtn.disabled = !!bad.length || records.length > 100 || !this.app.api.online || !!this.job;
    this.estimate.textContent = `${records.length} × ${traj} trajectories · about ${fmtDuration(total)} on ${cores} core${cores > 1 ? 's' : ''}` + (records.length > 100 ? ' · at most 100 sequences per batch' : '');
  }

  // ---------------------------------------------------------------- input helpers
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
      el('button', { class: 'linkbtn', onclick: () => app.openDialog('library') }, 'All projects…')), list);
  }

  // ---------------------------------------------------------------- prediction
  async startPrediction() {
    const app = this.app;
    if (!this.seqIn.value.trim()) {
      this.useDemo();
      toast(`No sequence was entered, so the demo sequence (${DEMO.label}) is predicted. Paste your own sequence or a FASTA file to predict it instead.`, 'info', 'Demo sequence', 7000);
    }
    const records = parseFastaRecords(this.seqIn.value);
    if (records.length > 1) return this.startBatch(records);
    const { seq, name: fastaName } = parseSequenceInput(this.seqIn.value);
    const name = this.nameIn.value.trim() || fastaName || 'rna';
    const params = { name: slug(name), sequence: seq, structure: this.ssIn.value.trim(), models: +this.modelsIn.value || 8,
      sampling: this.sampling, cores: +this.coresIn.value, seed: +this.seedIn.value || 42 };
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
        toast(`${j.document.n_conformations} conformation${j.document.n_conformations === 1 ? '' : 's'} of ${j.document.name}. Drag a residue with the Pull tool to switch between them.`, 'ok', 'Prediction finished', 7000);
      },
      onFail: j => { this.job = j; this._render(j); },
    });
  }

  async startBatch(records) {
    const app = this.app;
    const params = { records: records.map(r => ({ name: r.name, sequence: r.sequence, structure: r.structure })), models: +this.modelsIn.value || 8,
      sampling: this.sampling, cores: +this.coresIn.value, seed: +this.seedIn.value || 42 };
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
        el('td', {}, ok ? el('span', { class: 'state-ok', text: 'done' }) : el('span', { class: 'state-bad', text: it.error ? `failed: ${it.error.slice(0, 60)}` : it.status })));
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
    const meta = j.status === 'running' || j.status === 'queued'
      ? `${pct}% · ${j.eta_s ? 'about ' + fmtDuration(j.eta_s) + ' left' : 'estimating time'} · ${fmtDuration(j.elapsed_s)} elapsed`
      : j.status === 'done' ? `Finished in ${fmtDuration(j.elapsed_s)} — opening the viewer…` : j.status === 'cancelled' ? 'Cancelled' : `Failed: ${j.error || ''}`;
    const notes = el('ul', { class: 'wl-notes', 'aria-live': 'polite' }, ...(j.notes || []).map(n => el('li', { class: `nt-${n.level}` },
      el('span', { class: 'nt-t', text: `${Math.floor(n.t / 60)}:${String(Math.floor(n.t % 60)).padStart(2, '0')}` }), icon(n.level === 'error' ? 'warn' : n.level === 'warning' ? 'warn' : 'info'), el('span', { text: n.text }))));
    const actions = el('div', { class: 'wl-row end' });
    if (!finished) {
      actions.append(el('button', { class: 'btn', onclick: () => this.app.api.post(`/api/jobs/${j.id}/cancel`, {}) }, 'Cancel'));
    } else if (j.status !== 'done') {
      actions.append(el('button', { class: 'btn primary', onclick: () => { this.stopWatch?.(); this.job = null; this._showProgress(false); this._estimate(); } }, 'Back'));
    }
    this.progress.replaceChildren(el('div', { class: 'wl-pg-card' },
      el('div', { class: 'wl-pg-head' }, el('h2', { text: j.kind === 'batch' ? `Batch prediction · ${j.title || ''}` : `Predicting conformations · ${j.title || ''}` }),
        el('span', { class: 'muted', text: j.kind === 'batch' ? `${(j.result?.items || []).length} sequences, one after another, on ${j.result?.cores} CPU cores` : j.result?.cores ? `${j.result.trajectories} trajectories on ${j.result.cores} CPU cores` : '' })),
      steps,
      el('div', { class: `progress big ${j.status === 'failed' ? 'bad' : ''}` }, el('i', { style: { width: `${pct}%` } })),
      el('div', { class: 'wl-pg-meta', text: meta }),
      el('div', { class: 'section-title', style: { marginTop: '14px' } }, 'Notes'),
      notes, actions));
    notes.scrollTop = notes.scrollHeight;
    this._estimate();
  }
}
