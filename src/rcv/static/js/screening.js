// Compound screening (experimental): a compound library against two conformations of one RNA, ranked by
// how much better each compound fits conformation A than conformation B.
//
// This is a page of its own (Tools ▸ Compound screening…, and the button in the Switch tab), not a part of
// the main layout: a screen is a long job with many parameters, and the main window stays as it is. Four
// cards — target, library, parameters, run — then the ranked table, with "Show in 3D" putting a pose back
// into the main window as an overlay with the nucleotides it touches selected.
//
// Every number here is a ranking. The header says what the scoring function was measured at and links to
// the status document; nothing in this file calls anything an affinity, and no compound is said to block a
// transition. Server side: rcv/screening.py; engine: rnafold3d.screening.
import { el, fmt, fmtSigned, saveBlob, loadPref, savePref } from './util.js';
import { icon } from './icons.js';
import { toast, field, select, askConfirm } from './ui.js';

// What the scoring function is measured at (docs/SCREENING_STATUS.md). Shown in the window, as the plan
// (docs/NEXT_STEPS.md §1.4) requires while the screen is experimental: one place, so it cannot go stale
// in three.
export const MEASURED = {
  redock: 'redocking 30 experimental complexes: median 4.24 Å, the right pose first in 10 of 30',
  pockets: 'pocket finding: the site of a known compound is the top pocket in 21 of 30 structures',
  enrichment: 'known binders against property-matched decoys, 6 targets: mean AUC 0.864 '
    + '(0.50 is random), worst 0.50, best 1.00',
  warning: 'A ranking to choose a few compounds to test, not an affinity and not a prediction that a '
    + 'compound will bind or block anything.',
};

const EXH = { quick: 'Quick', standard: 'Standard', thorough: 'Thorough' };
const POLL = 1200;

export function openScreening(app) {
  if (!app.ens) { toast('Open or predict an ensemble first.', 'error', 'Compound screening'); return null; }
  if (!app.api.online) { toast('Screening needs the local RCV engine.', 'error', 'Compound screening'); return null; }
  app.screening = app.screening || new ScreenWindow(app);
  app.screening.show();
  return app.screening;
}

class ScreenWindow {
  constructor(app) {
    this.app = app;
    this.prefs = loadPref('screen', {});
    this.library = null;          // the parsed library summary from the server
    this.job = null;              // the running or finished job
    this.rows = [];
    this.sort = { key: 'selectivity', dir: -1 };
    this.filters = { selectivity: null, weight: null, text: '' };
    this.root = document.getElementById('screenpage');
    if (!this.root) {
      this.root = el('div', { class: 'welcome scr-page', id: 'screenpage', hidden: true });
      document.body.append(this.root);
    }
    this._build();
    app.on('ensemble', () => { if (this.app.page === 'screen') this._fillConformations(); });
  }

  show() {
    this._fillConformations();
    this.app.showPage('screen');
    this._pockets();
  }

  close() { this.app.showPage('viewer'); }

  // ---------------------------------------------------------------- layout
  _build() {
    const app = this.app;
    const top = el('header', { class: 'wl-top' },
      el('div', { class: 'brand' }, icon('flask'), el('span', { class: 'brand-name', text: 'Compound screening' }),
        el('span', { class: 'badge warn', text: 'experimental' })),
      el('span', { class: 'grow' }),
      el('button', { class: 'btn ghost', onclick: () => import('./dialogs.js').then(m => m.userGuideDialog(app)) }, icon('book'), 'User guide'),
      el('button', { class: 'btn', onclick: () => this.close() }, icon('prev'), 'Back to the viewer'));
    this.note = el('div', { class: 'scr-note' },
      el('p', { text: MEASURED.warning }),
      el('ul', {}, ...[MEASURED.redock, MEASURED.enrichment, MEASURED.pockets]
        .map(t => el('li', { text: t }))));
    this.cards = el('div', { class: 'scr-cards' },
      this._targetCard(), this._libraryCard(), this._paramCard(), this._runCard());
    this.results = el('section', { class: 'wl-card scr-results', hidden: true });
    const main = el('main', { class: 'wl-main scr-main' },
      el('div', { class: 'wl-hero' }, el('h1', { text: 'Which compounds prefer one conformation over the other' }),
        el('p', { text: 'Pick the two conformations and the site, choose a compound file, and the screen docks every compound into both and ranks them by the difference.' })),
      this.note, this.cards, this.results);
    this.root.replaceChildren(top, main);
  }

  _head(n, title, sub) {
    return el('div', { class: 'wl-card-head' }, el('span', { class: 'wl-badge', translate: 'no', text: n }),
      el('div', {}, el('h2', { text: title }), el('p', { text: sub })));
  }

  // ---------------------------------------------------------------- 1 target
  _targetCard() {
    this.selA = select([], '0', { 'aria-label': 'Conformation A' });
    this.selB = select([], '', { 'aria-label': 'Conformation B' });
    this.selA.addEventListener('change', () => { this._pockets(); this._estimate(); });
    this.selB.addEventListener('change', () => this._estimate());
    this.siteMode = select([
      { value: 'pocket', label: 'a pocket found automatically' },
      { value: 'residues', label: 'the nucleotides selected in the main window' },
    ], this.prefs.siteMode || 'pocket', { 'aria-label': 'Where to dock' });
    this.siteMode.addEventListener('change', () => { this._siteUI(); this._estimate(); });
    this.selPocket = select([], '1', { 'aria-label': 'Pocket' });
    this.selPocket.addEventListener('change', () => this._siteInfo());
    this.burial = select([{ value: '5', label: 'sites (enclosed on 5 of 7 sides)' },
      { value: '6', label: 'only deep pockets (6 of 7)' }], String(this.prefs.burial || 5), { 'aria-label': 'How enclosed' });
    this.burial.addEventListener('change', () => this._pockets());
    this.pocketRow = el('div', { class: 'grid2' }, field('Pocket', this.selPocket), field('Detect', this.burial));
    this.siteInfo = el('div', { class: 'desc' });
    return el('section', { class: 'wl-card' },
      this._head(1, 'Target', 'two conformations of the open ensemble, and the site to dock into'),
      el('div', { class: 'grid2' },
        field('Conformation A (the one to stabilise)', this.selA),
        field('Conformation B (the one to block)', this.selB, 'without a second conformation the list is ranked by score alone')),
      field('Site', this.siteMode), this.pocketRow, this.siteInfo);
  }

  _fillConformations() {
    const ens = this.app.ens;
    if (!ens) return;
    const names = ens.conformations.map((c, i) => ({ value: String(i), label: c.label || `Conformation ${i + 1}` }));
    const keepA = this.selA.value, keepB = this.selB.value;
    this.selA.replaceChildren(...names.map(o => el('option', { value: o.value, text: o.label })));
    this.selB.replaceChildren(el('option', { value: '', text: '— none: rank by score only —' }),
      ...names.map(o => el('option', { value: o.value, text: o.label })));
    // prefilled from what the main window shows: A is the conformation on screen, B the next one
    const shown = this.app.state?.conf ?? 0;
    this.selA.value = keepA && Number(keepA) < names.length ? keepA : String(shown);
    this.selB.value = keepB && keepB !== '' && Number(keepB) < names.length ? keepB
      : (names.length > 1 ? String((shown + 1) % names.length) : '');
    this._siteUI();
  }

  _siteUI() {
    this.pocketRow.hidden = this.siteMode.value !== 'pocket';
    this._siteInfo();
  }

  async _pockets() {
    if (this.siteMode.value !== 'pocket' || !this.app.ens) return;
    this.selPocket.replaceChildren(el('option', { value: '1', text: 'looking for pockets…' }));
    try {
      const res = await this.app.api.post('/api/insight/pockets', {
        document: this.app.ens.id, conformation: Number(this.selA.value), burial_min: Number(this.burial.value),
      });
      this.selPocket.replaceChildren(...res.pockets.map(p => el('option', {
        value: String(p.pocket),
        text: `${p.pocket}. ${fmt(p.volume, 0)} Å³, ${p.residues.length} nt${p.wide ? ', wide' : ''}`,
      })));
      if (!res.pockets.length) this.selPocket.replaceChildren(el('option', { value: '1', text: 'no pocket found' }));
      this._siteInfo();
    } catch (err) { this.siteInfo.replaceChildren(el('span', { class: 'bad', text: err.message })); }
  }

  siteSpec() {
    if (this.siteMode.value === 'residues') {
      return { mode: 'residues', residues: [...(this.app.state?.selection || [])] };
    }
    return { mode: 'pocket', pocket: Number(this.selPocket.value || 1) };
  }

  async _siteInfo() {
    const app = this.app;
    if (!app.ens) return;
    this.siteInfo.replaceChildren(el('span', { class: 'muted', text: 'measuring the site…' }));
    try {
      const s = await app.api.post('/api/screen/site', {
        document: app.ens.id, conformation: Number(this.selA.value),
        site: this.siteSpec(), burial_min: Number(this.burial.value),
      });
      this.site = s;
      const parts = [`a box of ${fmt(s.size, 0)} Å`, s.residues?.length ? `${s.residues.length} nucleotides lining it` : null,
        s.volume ? `${fmt(s.volume, 0)} Å³ of open space` : null].filter(Boolean);
      this.siteInfo.replaceChildren(el('span', { text: parts.join(' · ') }),
        ...(s.cropped ? [el('div', { class: 'bad', text: 'This site is wider than the largest box a screen grids. Only its centre is screened — select the nucleotides you mean in the main window instead.' })] : []));
      this._estimate();
    } catch (err) {
      this.site = null;
      this.siteInfo.replaceChildren(el('span', { class: 'bad', text: err.message }));
    }
  }

  // ---------------------------------------------------------------- 2 library
  _libraryCard() {
    this.libBtn = el('button', { class: 'btn', onclick: () => this._chooseLibrary() }, icon('open'), 'Choose a compound file…');
    this.demoBtn = el('button', { class: 'btn ghost', title: 'Four small compounds that come with RCV, to see what a run does', onclick: () => this.loadExample() }, icon('flask'), 'Use the example file');
    this.libInfo = el('div', { class: 'desc', text: 'An SDF, MOL or MOL2 file with 3D coordinates. Nothing ships with RCV: a compound is docked as the file draws it, so the file has to be prepared (3D, hydrogens, one conformer per record).' });
    this.libFile = el('input', { type: 'file', accept: '.sdf,.sd,.mol,.mol2', hidden: true });
    this.libFile.addEventListener('change', () => this._uploadLibrary());
    return el('section', { class: 'wl-card' },
      this._head(2, 'Library', 'the compounds to try'),
      el('div', { class: 'row' }, this.libBtn, this.demoBtn), this.libInfo, this.libFile);
  }

  async _chooseLibrary() {
    const app = this.app;
    if (app.native?.available()) {
      const paths = await app.native.chooseFiles('compounds').catch(() => null);
      if (!paths || !paths.length) return;
      await this._load(() => app.api.raw('POST', '/api/screen/library', null, { 'X-Path': encodeURIComponent(paths[0]) }));
      return;
    }
    this.libFile.click();
  }

  // The four compounds that ship with RCV (static/examples/demo-compounds.sdf): purines and xanthines,
  // small enough that a run finishes while you watch it. Not a library — a way to see what the window does.
  async loadExample() {
    await this._load(async () => {
      const text = await (await fetch('examples/demo-compounds.sdf')).arrayBuffer();
      return this.app.api.raw('POST', '/api/screen/library', text,
        { 'Content-Type': 'application/octet-stream', 'X-Name': encodeURIComponent('demo-compounds.sdf') });
    });
  }

  _uploadLibrary() {
    const f = this.libFile.files[0];
    if (!f) return;
    this._load(async () => this.app.api.raw('POST', '/api/screen/library', await f.arrayBuffer(),
      { 'Content-Type': 'application/octet-stream', 'X-Name': encodeURIComponent(f.name) }));
  }

  async _load(request) {
    this.libInfo.replaceChildren(el('span', { class: 'muted', text: 'reading the compounds…' }));
    try {
      const res = await (await request()).json();
      this.library = res;
      const bad = res.rejected_total
        ? el('div', { class: 'bad', text: `${res.rejected_total} records could not be read: ` + Object.entries(res.rejected).map(([w, n]) => `${n} × ${w}`).join(', ') })
        : '';
      this.libInfo.replaceChildren(
        el('div', {}, el('b', { text: res.name }), ` · ${res.count} compounds`),
        el('div', { class: 'muted', text: res.count ? `${res.atoms[0]}–${res.atoms[2]} heavy atoms (median ${res.atoms[1]}), ${res.weight[0]}–${res.weight[1]} Da, charge ${fmtSigned(res.charge[0], 0)} to ${fmtSigned(res.charge[1], 0)} at pH 7, up to ${res.rotatable} rotatable bonds` : '' }),
        el('div', { class: 'muted small', text: res.examples.map(e => e.name).join(', ') }), bad);
      this._estimate();
    } catch (err) {
      this.library = null;
      this.libInfo.replaceChildren(el('span', { class: 'bad', text: err.message }));
    }
  }

  // ---------------------------------------------------------------- 3 parameters
  _paramCard() {
    this.exh = select(Object.entries(EXH).map(([v, l]) => ({ value: v, label: l })), this.prefs.exh || 'quick', { 'aria-label': 'Search depth' });
    this.poses = el('input', { class: 'input', type: 'number', min: 1, max: 9, value: this.prefs.poses || 3, 'aria-label': 'Poses per compound' });
    this.keep = el('input', { class: 'input', type: 'number', min: 0, max: 5000, value: this.prefs.keep ?? 200, 'aria-label': 'Compounds whose pose is kept' });
    this.cores = el('input', { class: 'input', type: 'number', min: 1, max: 64, value: this.prefs.cores || Math.max(1, (navigator.hardwareConcurrency || 4) - 1), 'aria-label': 'Cores' });
    this.limit = el('input', { class: 'input', type: 'number', min: 0, max: 200000, value: this.prefs.limit || 0, 'aria-label': 'Compounds at most' });
    for (const c of [this.exh, this.poses, this.cores, this.limit]) c.addEventListener('change', () => this._estimate());
    this.estimate = el('div', { class: 'desc' });
    return el('section', { class: 'wl-card' },
      this._head(3, 'Parameters', 'how hard to search, and how much to keep'),
      el('div', { class: 'grid2' },
        field('Search', this.exh, 'quick is for a first look; thorough is minutes per compound'),
        field('Poses kept per compound and state', this.poses),
        field('Compounds whose pose is kept', this.keep, 'the poses of the whole library would be gigabytes'),
        field('Cores', this.cores),
        field('At most this many compounds (0 = all)', this.limit)),
      el('div', { class: 'desc', text: 'Compounds are ionised at pH 7 (carboxylate, phosphate, amine, amidine; see the status document for what is not modelled). The RNA is rigid in both conformations.' }));
  }

  _estimate() {
    const n = Math.min(this.library?.count || 0, Number(this.limit.value) || Infinity);
    const states = this.selB.value === '' ? 1 : 2;
    const per = { quick: 4, standard: 18, thorough: 90 }[this.exh.value] * (this.library?.atoms?.[1] ? Math.max(0.4, this.library.atoms[1] / 20) : 1);
    const cores = Math.max(1, Number(this.cores.value) || 1);
    const secs = n * states * per / cores;
    this.estimate.replaceChildren(el('span', {
      text: n ? `about ${human(secs)} for ${n} compounds × ${states} state${states > 1 ? 's' : ''} on ${cores} core${cores > 1 ? 's' : ''} — a rough guess from the compound size, not a measurement.`
        : 'Choose a compound file to see how long a run would take.',
    }));
    if (this.startBtn) this.startBtn.disabled = !n || !this.app.ens;
    return secs;
  }

  // ---------------------------------------------------------------- 4 run
  _runCard() {
    this.startBtn = el('button', { class: 'btn primary', disabled: true, onclick: () => this._start() }, icon('play'), 'Start the screen');
    this.stopBtn = el('button', { class: 'btn', hidden: true, onclick: () => this._stop() }, icon('x'), 'Stop');
    this.bar = el('div', { class: 'scr-bar' }, el('span', { style: { width: '0%' } }));
    this.status = el('div', { class: 'desc' });
    this.stages = el('div', { class: 'scr-stages' });
    this.hist = el('div', { class: 'scr-hist' });
    return el('section', { class: 'wl-card' },
      this._head(4, 'Run', 'progress, and the table fills as it goes'),
      this.estimate, el('div', { class: 'row' }, this.startBtn, this.stopBtn), this.bar, this.status,
      this.stages, this.hist);
  }

  async _start() {
    const app = this.app;
    if (!this.library) { toast('Choose a compound file first.', 'error', 'Compound screening'); return; }
    if (this.siteMode.value === 'residues' && !(this.app.state?.selection?.size)) {
      toast('Select the nucleotides of the site in the main window, or choose a pocket.', 'error', 'Compound screening');
      return;
    }
    const body = {
      document: app.ens.id, a: Number(this.selA.value), b: this.selB.value === '' ? null : Number(this.selB.value),
      site: this.siteSpec(), burial_min: Number(this.burial.value), exhaustiveness: this.exh.value,
      n_poses: Number(this.poses.value), keep_poses: Number(this.keep.value), jobs: Number(this.cores.value),
      limit: Number(this.limit.value) || null, library: this.library.handle,
    };
    this.prefs = { siteMode: this.siteMode.value, burial: Number(this.burial.value), exh: this.exh.value,
      poses: Number(this.poses.value), keep: Number(this.keep.value), cores: Number(this.cores.value),
      limit: Number(this.limit.value) };
    savePref('screen', this.prefs);
    try {
      this.job = await app.api.post('/api/screen', body);
      this.startBtn.hidden = true;
      this.stopBtn.hidden = false;
      this.rows = [];
      this._poll();
    } catch (err) { toast(err.message, 'error', 'Compound screening'); }
  }

  async _stop() {
    if (!this.job) return;
    if (!(await askConfirm('Stop the screen?', 'The compounds already done stay in the table.', 'Stop'))) return;
    try { await this.app.api.post(`/api/jobs/${this.job.id}/cancel`); } catch { /* it may have just finished */ }
  }

  async _poll() {
    if (!this.job) return;
    let info = null;
    try { info = await this.app.api.get(`/api/screen/${this.job.id}`); } catch (err) { this.status.textContent = err.message; return; }
    this.job = info;
    const r = info.result || {};
    const pct = Math.round(100 * (info.progress || 0));
    this.bar.firstChild.style.width = `${pct}%`;
    const bits = [`${r.done || 0} of ${r.total || '?'} compounds`];
    if (r.per_second) bits.push(`${fmt(r.per_second, 2)} per second`);
    if (info.eta_s) bits.push(`about ${human(info.eta_s)} left`);
    if (r.failed) bits.push(`${r.failed} could not be placed`);
    this.status.replaceChildren(el('span', { text: bits.join(' · ') }),
      ...(info.error ? [el('div', { class: 'bad', text: info.error })] : []));
    this.stages.replaceChildren(...(info.stages || []).map(s => el('span', { class: `scr-stage ${s.state}` },
      el('b', { text: s.label }), s.detail ? ` ${s.detail}` : '')));
    if ((r.rows || []).length) { this.rows = r.rows; this._renderResults(); this._histogram(); }
    if (['running', 'queued'].includes(info.status)) { setTimeout(() => this._poll(), POLL); return; }
    this.startBtn.hidden = false;
    this.stopBtn.hidden = true;
    if (info.status === 'done') {
      await this._loadAll();
      toast(`${this.rows.length} compounds ranked.`, 'ok', 'Compound screening');
    } else if (info.status === 'cancelled') {
      await this._loadAll();
      toast('Stopped. The compounds already done are in the table.', 'warning', 'Compound screening');
    }
  }

  // The spread of the selectivity scores so far, as twenty bars: a run where every compound scores the same
  // is a run whose site does not differ between the two conformations, and that is worth seeing early.
  _histogram() {
    const vals = this.rows.map(r => r.selectivity).filter(v => v != null);
    if (vals.length < 8) { this.hist.replaceChildren(); return; }
    const lo = Math.min(...vals), hi = Math.max(...vals), span = hi - lo || 1;
    const bins = new Array(20).fill(0);
    for (const v of vals) bins[Math.min(19, Math.floor(20 * (v - lo) / span))]++;
    const top = Math.max(...bins);
    this.hist.replaceChildren(
      el('div', { class: 'scr-bars' }, ...bins.map((n, i) => el('span', {
        style: { height: `${Math.round(100 * n / top)}%` },
        title: `${n} compounds, selectivity ${fmt(lo + span * i / 20, 1)} to ${fmt(lo + span * (i + 1) / 20, 1)}`,
      }))),
      el('div', { class: 'muted small', text: `selectivity ${fmtSigned(lo, 1)} to ${fmtSigned(hi, 1)} over ${vals.length} compounds` }));
  }

  async _loadAll() {
    try {
      const res = await this.app.api.get(`/api/screen/${this.job.id}/results`);
      this.rows = res.rows || [];
      this.summary = res.summary || {};
      this._renderResults();
      this._histogram();
    } catch (err) { toast(err.message, 'error', 'Compound screening'); }
  }

  // ---------------------------------------------------------------- results
  _renderResults() {
    const two = this.job && this.job.result && this.job.result.two_state !== false && this.selB.value !== '';
    const rows = this._filtered();
    this.results.hidden = false;
    const cols = [['name', 'Compound'], ['weight', 'Da'], ['atoms', 'Atoms'], ['rotatable', 'Rot.'],
      ['charge', 'Charge'], ['score_a', 'Score A'], ...(two ? [['score_b', 'Score B'], ['selectivity', 'Selectivity']] : []),
      ['strain_a', 'Strain']];
    const head = el('tr', {}, el('th', { text: '#' }), ...cols.map(([k, label]) => el('th', {
      class: `sortable${this.sort.key === k ? ' sorted' : ''}`, onclick: () => this._sortBy(k),
    }, label, this.sort.key === k ? (this.sort.dir < 0 ? ' ▾' : ' ▴') : '')), el('th', { text: '' }));
    const body = rows.slice(0, 400).map((r, i) => el('tr', { class: r.error ? 'muted' : '' },
      el('td', { class: 'num', text: String(i + 1) }),
      el('td', { title: r.formula || '', text: r.name }),
      el('td', { class: 'num', text: fmt(r.weight, 0) }),
      el('td', { class: 'num', text: String(r.atoms) }),
      el('td', { class: 'num', text: String(r.rotatable) }),
      el('td', { class: 'num', text: fmtSigned(r.charge, 0) }),
      el('td', { class: 'num', text: r.score_a == null ? '—' : fmt(r.score_a, 2) }),
      ...(two ? [el('td', { class: 'num', text: r.score_b == null ? '—' : fmt(r.score_b, 2) }),
        el('td', { class: `num ${r.selectivity > 0 ? 'good' : ''}`, text: r.selectivity == null ? '—' : fmtSigned(r.selectivity, 2) })] : []),
      el('td', { class: 'num', title: 'how much the pose distorts the compound itself', text: fmt(r.strain_a, 1) }),
      el('td', {}, r.error ? el('span', { class: 'muted small', text: r.error })
        : el('button', { class: 'btn small ghost', onclick: () => this._show3D(r) }, icon('cube'), '3D'))));
    const table = el('table', { class: 'table scr-tbl' }, el('thead', {}, head), el('tbody', {}, ...body));
    this.results.replaceChildren(
      el('div', { class: 'wl-card-head' }, el('div', {}, el('h2', { text: 'Ranked compounds' }),
        el('p', { text: two ? 'Selectivity is score(B) − score(A): positive means the compound prefers conformation A. Ties go to the compound with the better of its two worse scores.' : 'Ranked by score in the one conformation.' })),
        el('span', { class: 'grow' }),
        el('button', { class: 'btn small ghost', onclick: () => this._copy() }, icon('copy'), 'Copy as table'),
        el('button', { class: 'btn small ghost', onclick: () => this._download('csv') }, icon('download'), 'CSV'),
        el('button', { class: 'btn small ghost', onclick: () => this._download('sdf') }, icon('download'), 'Poses (SDF)')),
      this._filterRow(rows.length),
      el('div', { class: 'scr-tblwrap' }, table),
      rows.length > 400 ? el('p', { class: 'muted small', text: `showing the first 400 of ${rows.length}; the CSV has all of them` }) : '',
      el('p', { class: 'muted small', text: MEASURED.warning }));
  }

  _filterRow(n) {
    const selIn = el('input', { class: 'input small', type: 'number', step: 0.5, placeholder: 'any', value: this.filters.selectivity ?? '', 'aria-label': 'Minimum selectivity' });
    const wIn = el('input', { class: 'input small', type: 'number', step: 10, placeholder: 'any', value: this.filters.weight ?? '', 'aria-label': 'Largest molecular weight' });
    const txt = el('input', { class: 'input small', placeholder: 'name contains…', value: this.filters.text, 'aria-label': 'Name filter' });
    const apply = () => {
      this.filters = { selectivity: selIn.value === '' ? null : Number(selIn.value),
        weight: wIn.value === '' ? null : Number(wIn.value), text: txt.value.trim() };
      this._renderResults();
    };
    for (const c of [selIn, wIn, txt]) c.addEventListener('change', apply);
    return el('div', { class: 'scr-filters' }, el('span', { class: 'muted', text: `${n} shown` }),
      field('Selectivity at least', selIn), field('Weight at most (Da)', wIn), field('Name', txt));
  }

  _filtered() {
    const f = this.filters;
    let rows = this.rows.filter(r => (f.selectivity == null || (r.selectivity ?? -Infinity) >= f.selectivity)
      && (f.weight == null || (r.weight ?? 0) <= f.weight)
      && (!f.text || String(r.name || '').toLowerCase().includes(f.text.toLowerCase())));
    const k = this.sort.key, d = this.sort.dir;
    rows = rows.slice().sort((a, b) => {
      const x = a[k], y = b[k];
      if (x == null && y == null) return 0;
      if (x == null) return 1;
      if (y == null) return -1;
      return typeof x === 'string' ? d * x.localeCompare(y) : d * (x - y);
    });
    return rows;
  }

  _sortBy(k) {
    this.sort = { key: k, dir: this.sort.key === k ? -this.sort.dir : (k === 'name' ? 1 : -1) };
    this._renderResults();
  }

  _copy() {
    const rows = this._filtered();
    const head = ['rank', 'compound', 'Da', 'atoms', 'rotatable', 'charge', 'score_A', 'score_B', 'selectivity', 'strain'];
    const text = [head.join('\t'), ...rows.map((r, i) => [i + 1, r.name, fmt(r.weight, 0), r.atoms, r.rotatable,
      r.charge, r.score_a ?? '', r.score_b ?? '', r.selectivity ?? '', r.strain_a].join('\t'))].join('\n');
    navigator.clipboard?.writeText(text).then(() => toast(`${rows.length} compounds copied as a tab-separated table.`, 'ok'));
  }

  async _download(kind) {
    try {
      const r = await this.app.api.raw('GET', `/api/screen/${this.job.id}/${kind}`);
      const blob = await r.blob();
      await saveBlob(blob, `screen-${this.job.id}.${kind}`);
    } catch (err) { toast(err.message, 'error', kind.toUpperCase()); }
  }

  // Put one pose back into the main window: the compound as an overlay, the nucleotides it touches selected.
  async _show3D(row) {
    const app = this.app;
    try {
      const pose = await app.api.post(`/api/screen/${this.job.id}/pose`, { index: row.index, state: 'a' });
      showPose(app, pose);
      this.close();
      toast(`${pose.name}: the pose is drawn in conformation ${pose.conformation + 1} and the ${pose.residues.length} nucleotides it touches are selected.`, 'ok', 'Compound screening');
    } catch (err) { toast(err.message, 'error', 'Show in 3D'); }
  }
}

// The pose as a viewer overlay (points and bonds, projected over the 3D view) plus the contacting
// nucleotides selected. Clear it with View ▸ Clear, like the other overlays.
export function showPose(app, pose) {
  const COLOR = { C: '#444c56', N: '#2a78d6', O: '#d2402a', P: '#e08b1f', S: '#c2a52b', F: '#3aa655', CL: '#3aa655', BR: '#9a5b2e' };
  app.goTo(pose.conformation, { animate: false, force: true });
  if (pose.residues?.length) app.select(pose.residues, 'set');
  app.viewer.addOverlay('screen-pose', () => {
    const out = [];
    for (const [i, j] of pose.bonds) out.push({ type: 'line', a: pose.xyz[i], b: pose.xyz[j], color: '#555', width: 3 });
    for (let i = 0; i < pose.xyz.length; i++) {
      out.push({ type: 'point', p: pose.xyz[i], color: COLOR[String(pose.elements[i]).toUpperCase()] || '#777', r: 4 });
    }
    return out;
  });
  app.viewer.needsRender = true;
}

export function clearPose(app) { app.viewer?.removeOverlay('screen-pose'); }

function human(s) {
  if (!isFinite(s)) return 'a long time';
  if (s < 90) return `${Math.round(s)} s`;
  if (s < 5400) return `${Math.round(s / 60)} min`;
  if (s < 86400) return `${(s / 3600).toFixed(1)} h`;
  return `${(s / 86400).toFixed(1)} days`;
}
