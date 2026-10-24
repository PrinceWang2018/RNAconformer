// RCV (RNA Conformation Viewer) — application controller.
import { el, $, $$, Emitter, clamp, easeInOut, fmt, fmtVal, loadPref, savePref, saveBlob, readFileText, slug, luminance, debounce, ensureFigureFont } from './util.js';
import { icon } from './icons.js';
import { Ensemble, toPDB, DOC_SCHEMAS } from './ensemble.js';
import { Viewer, STYLES } from './viewer/viewer.js';
import { SCHEMES, BACKGROUNDS, residueColors, legend as schemeLegend } from './colors.js';
import { tooltip, tipRow, toast, menuBar, popupMenu, busy } from './ui.js';
import { EnergyPanel } from './panels/energy.js';
import { EnsembleList, ConformationList, Timeline } from './panels/conformations.js';
import { DetailsPanel } from './panels/details.js';
import { PluginHost } from './plugins.js';
import { openDialog, exportImageDialog, exportVideoDialog, preferencesDialog, pluginsDialog, shortcutsDialog, aboutDialog, methodsDialog, panelsDialog, pluginGuideDialog, newPluginDialog, userGuideDialog, openWebsite } from './dialogs.js';
import { exportFigure, composeFigure } from './export/figure.js';
import { Welcome } from './welcome.js';
import { PullController } from './pullui.js';
import { JobsDock } from './jobs.js';
import { Native } from './native.js';
import { SecondaryInset } from './panels/ss2d.js';
import { SequenceBar } from './panels/seqbar.js';
import { ConsolePanel } from './console/panel.js';
import { Measurements } from './measure.js';
import { Scenes } from './scenes.js';
import { SwitchPanel } from './panels/switch.js';
import { fetchDialog, toolsMenuItems } from './tools.js';
import { setPdfFontProvider } from './export/encoders.js';
import { installNativeFeel } from './nativefeel.js';
import { Palettes, paletteDialog } from './palettes.js';
import { Annotations } from './annotations.js';
import { MouseTools } from './tools3d.js';
import { History, clearItems } from './history.js';
import { plotStudio } from './plots.js';
import { SurfaceController, transparencyDialog } from './surfacectl.js';
import { composePanels, exportPanels } from './export/panels.js';
import { recordVideo, framePlan } from './export/video.js';
import { initLanguage, setLanguage, currentLanguage, LANGUAGES, languageName, untranslated, t as tr } from './i18n.js';
import { WindowManager } from './windows.js';

// ------------------------------------------------------------------ server API
class Api {
  constructor() {
    const m = /token=([A-Za-z0-9_-]+)/.exec(location.hash);
    let tok = m ? m[1] : null;
    try {
      if (tok) sessionStorage.setItem('rcv.token', tok); else tok = sessionStorage.getItem('rcv.token');
    } catch { /* storage unavailable */ }
    this.token = tok;
    this.online = false;
    if (m) history.replaceState(null, '', location.pathname + location.search);
  }
  async req(method, path, body) {
    const opt = { method, headers: { 'X-Viewer-Token': this.token || '' } };
    if (body !== undefined) { opt.headers['Content-Type'] = 'application/json'; opt.body = JSON.stringify(body); }
    const r = await fetch(path, opt);
    let data = null;
    try { data = await r.json(); } catch { /* non-JSON */ }
    if (!r.ok) throw new Error((data && data.error) || `${r.status} ${r.statusText}`);
    return data;
  }
  // raw request (binary bodies / responses)
  async raw(method, path, body, headers = {}) {
    const r = await fetch(path, { method, body, headers: { 'X-Viewer-Token': this.token || '', ...headers } });
    if (!r.ok) {
      let msg = `${r.status} ${r.statusText}`;
      try { msg = (await r.json()).error || msg; } catch { /* non-JSON */ }
      throw new Error(msg);
    }
    return r;
  }
  get(p) { return this.req('GET', p); }
  post(p, b) { return this.req('POST', p, b ?? {}); }
  del(p) { return this.req('DELETE', p); }
}

const hexToLin = hex => {
  const n = parseInt(hex.slice(1), 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255].map(v => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); });
};

// ------------------------------------------------------------------ app
class App extends Emitter {
  constructor() {
    super();
    this.api = new Api();
    this.prefs = loadPref('prefs', {});
    this.ensembles = new Map();
    this.ens = null;
    this.order = [];
    this.display = null;
    this.morph = null;
    this.playing = false;
    this.morphDuration = this.prefs.morphDuration ?? 1.2;
    this.holdTime = this.prefs.holdTime ?? 0.8;
    this.loop = this.prefs.loop ?? true;
    this.state = {
      conf: 0, prev: null, metric: null,
      orderMode: this.prefs.orderMode || 'rank',
      colorScheme: this.prefs.colorScheme || 'nucleotide',
      uniformColor: this.prefs.uniformColor || '#5b7a99',
      selection: new Set(),
      ghosts: false, putty: false,
    };
    this._colorCache = new Map();
    this.pluginTabs = [];
    // export functions of this window's modules (for plugins and automated tests)
    this.lib = { exportFigure, composeFigure, composePanels, exportPanels, recordVideo, framePlan };
    this.i18n = { untranslated, t: tr, current: currentLanguage, languages: LANGUAGES };
  }

  savePrefs() {
    Object.assign(this.prefs, { orderMode: this.state.orderMode, colorScheme: this.state.colorScheme, morphDuration: this.morphDuration, holdTime: this.holdTime, loop: this.loop,
      style: this.viewer.settings.style, representation: this.viewer.settings.representation, background: this.viewer.settings.background, labelMode: this.viewer.labelMode });
    savePref('prefs', this.prefs);
  }

  // ---------------------------------------------------------------- startup
  async init() {
    await initLanguage();             // interface language (translates the page as it is built)
    ensureFigureFont();
    installNativeFeel(this);
    this.palettes = new Palettes(this);
    this.native = new Native(this.api);
    window.__native = this.native;
    this.setTheme(this.prefs.theme || 'system', false);
    this.viewer = new Viewer({
      canvas: $('#gl'), container: $('#viewport'), overlay: $('#overlay'), labels: $('#labels'),
      onPick: (p, e) => this._onPick(p, e), onHover: (p, e) => this._onHover(p, e), onDblPick: p => { if (p) this.viewer.centerOn(this.viewer.labelAnchor(p.residue)); },
    });
    const vs = this.viewer.settings;
    vs.style = this.prefs.style && STYLES[this.prefs.style] ? this.prefs.style : 'studio';
    vs.representation = this.prefs.representation || 'cartoon';
    vs.background = this.prefs.background || '#ffffff';
    vs.quality = this.prefs.quality || 'high';
    vs.outlineColor = this.prefs.outlineColor || '#141414';
    vs.elementColors = this.prefs.elementColors ?? true;
    this.viewer.labelMode = this.prefs.labelMode || 'selection';
    if (this.prefs.projection) this.viewer.setProjection(this.prefs.projection);
    this._buildChrome();
    this.wm = new WindowManager(this);
    this.energyPanel = new EnergyPanel(this, this.tabBody('energy'));
    this.details = new DetailsPanel(this, this.tabBody('details'));
    this.switchPanel = new SwitchPanel(this, this.tabBody('switch'));
    new EnsembleList(this, $('#ensemble-list'));
    new ConformationList(this, $('#conf-list'));
    this.timeline = new Timeline(this, $('#tl-track'), $('#tl-pos'));
    this._bindKeys();
    this._bindDrop();
    this._applyBackground();
    this.registerSessionPart('display', {
      save: ens => { const d = ens && ens.display; if (!d) return null; return { hidden: [...d.hidden], sticks: [...d.sticks], spheres: [...d.spheres], labels: [...d.labels], colors: [...d.colors], atomColors: [...(d.atomColors || [])], named: [...d.named].map(([k, v]) => [k, [...v]]) }; },
      load: (ens, v) => {
        const d = this.displayOf(ens);
        d.hidden = new Set(v?.hidden || []); d.sticks = new Set(v?.sticks || []); d.spheres = new Set(v?.spheres || []);
        d.labels = new Map(v?.labels || []); d.colors = new Map(v?.colors || []); d.atomColors = new Map(v?.atomColors || []); d.named = new Map((v?.named || []).map(([k, s]) => [k, new Set(s)]));
        d.version++;
        this.applyDisplay();
        this.emit('named');
      },
    });
    this.measure = new Measurements(this);
    this.annotations = new Annotations(this);
    this.surface = new SurfaceController(this);
    this.scenes = new Scenes(this);
    this.ss2d = new SecondaryInset(this, $('#viewport'));
    this.seqbar = new SequenceBar(this, $('#viewport'));
    this.console = new ConsolePanel(this, $('#console-host'));
    this.wm.decorate();
    this._bindContextMenu();
    this.pullCtl = new PullController(this);
    this.tools = new MouseTools(this);
    this.history = new History(this);
    $('#tool-host').replaceWith(this.tools.root);
    this.viewer.controls.interceptor = e => this.pullCtl.onPointerDown(e) || this.tools.onPointerDown(e);
    this.jobsDock = new JobsDock(this, $('#jobs-dock'));
    this.welcome = new Welcome(this, $('#welcome'));
    this._bindLifecycle();
    this.plugins = new PluginHost(this);
    let info = null;
    try {
      info = await this.api.get('/api/info');
      this.api.online = true;
      this.info = info;
    } catch (err) {
      this.api.online = false;
      toast('The viewer server is not reachable. Sessions and dropped files still open; predictions and folder browsing need the server (python -m rnafold3d view).', 'error', 'Offline');
    }
    this._status();
    // vector PDF figures embed the bundled figure font (Liberation Sans, subset by the local server)
    if (info) setPdfFontProvider(async (regular, bold) => { const r = await this.api.post('/api/pdf-font', { regular, bold }); return r.available ? r : null; });
    if (info) {
      await this.plugins.loadAll(info.plugins || []);
      const docs = info.documents || [];
      for (const d of docs) { try { await this.loadServerDocument(d.id, { select: false }); } catch (err) { toast(err.message, 'error', d.name); } }
      if (docs.length) {
        this.selectEnsemble([...this.ensembles.keys()][0]);
        if (docs[0].session?.view) this.applyViewState(docs[0].session.view);
      }
    }
    this.rebuildMenus();
    this.emit('documents');
    this.emit('ensemble');
    this.showPage(this.ensembles.size ? 'viewer' : 'welcome');
    document.body.dataset.ready = '1';
  }

  // ---------------------------------------------------------------- pages & lifecycle
  showPage(page) {
    this.page = page;
    document.body.dataset.page = page;
    $('#welcome').hidden = page !== 'welcome';
    $('#app').hidden = page !== 'viewer';
    if (page === 'welcome') this.welcome.refresh();
    else { this.viewer.resize(); this.viewer.needsRender = true; }
    this.emit('page', page);
  }

  _bindLifecycle() {
    // closing the window quits the app (the local engine process stops shortly after)
    window.addEventListener('pagehide', () => {
      try { navigator.sendBeacon(`/api/bye?token=${encodeURIComponent(this.api.token || '')}`, ''); } catch { /* ignore */ }
    });
    window.addEventListener('beforeunload', e => {
      const running = this.jobsDock?.running || (this.welcome?.job && !['done', 'failed', 'cancelled'].includes(this.welcome.job.status));
      if (running) { e.preventDefault(); e.returnValue = 'A computation is still running.'; }
    });
  }

  openDialog(tab) { return openDialog(this, tab); }

  async openDemo() {
    const done = busy('Opening the demo ensembles…');
    try {
      const list = await this.api.post('/api/open-demo', {});
      let first = null;
      for (const s of list) { const e = await this.loadServerDocument(s.id, { select: false }); first = first || e; }
      if (first) this.selectEnsemble(first.id);
    } catch (err) { toast(err.message, 'error', 'Demo'); } finally { done(); }
  }

  async openNativePaths(paths) {
    const sessions = paths.filter(p => /\.json$/i.test(p) && !/manifest\.json$/i.test(p));
    for (const p of sessions) {
      try { const r = await this.api.post('/api/read-json', { path: p }); await this.loadSessionText(r.text, p.split(/[\\/]/).pop()); }
      catch (err) { toast(err.message, 'error', 'Workspace'); }
    }
    const rcvs = paths.filter(p => /\.rcv$/i.test(p));
    for (const p of rcvs) await this.openServerPaths([p]);
    const rest = paths.filter(p => !sessions.includes(p) && !rcvs.includes(p));
    if (rest.length) await this.openServerPaths(rest);
  }

  // Swap an ensemble for its extended version (new conformations from a pull exploration), keeping the view.
  async replaceDocument(oldEnsId, newDocId) {
    const doc = await this.api.get(`/api/documents/${newDocId}`);
    const old = this.ensembles.get(oldEnsId);
    const ens = new Ensemble(doc, { serverId: newDocId });
    for (const m of this._clientMetrics || []) this._computeClientMetric(ens, m);
    ens.lastMetric = old?.lastMetric;
    // same RNA, more conformations: keep selections, colours, labels, measurements, scenes and the command log
    for (const key of ['display', 'measurements', 'scenes', 'commandLog']) if (old && old[key] !== undefined) ens[key] = old[key];
    const entries = [...this.ensembles.entries()].map(([id, e]) => (id === oldEnsId ? [ens.id, ens] : [id, e]));
    this.ensembles = new Map(entries);
    if (!this.ensembles.has(ens.id)) this.ensembles.set(ens.id, ens);
    if (!old || this.ens === old) {
      const curId = old?.conformations[this.state.conf]?.id;
      const k = Math.max(0, ens.conformations.findIndex(c => c.id === curId));
      this.selectEnsemble(ens.id, { keepView: true, conf: k, selection: [...this.state.selection] });
    }
    this.emit('documents');
    return ens;
  }

  // ---------------------------------------------------------------- chrome
  tabBody(id) { return this._tabBodies[id]; }

  _buildChrome() {
    // right tabs
    this._tabBodies = {};
    this._tabs = [];
    // tabs that do not fit fade out at the right edge; the mouse wheel scrolls them
    const tabs = $('#tabs');
    const fitTabs = () => tabs.classList.toggle('overflow', tabs.scrollWidth > tabs.clientWidth + 2 && tabs.scrollLeft + tabs.clientWidth < tabs.scrollWidth - 2);
    new ResizeObserver(fitTabs).observe(tabs);
    new MutationObserver(fitTabs).observe(tabs, { childList: true, subtree: true, characterData: true });
    tabs.addEventListener('scroll', fitTabs);
    tabs.addEventListener('wheel', e => { if (Math.abs(e.deltaY) > Math.abs(e.deltaX)) { tabs.scrollLeft += e.deltaY; e.preventDefault(); } }, { passive: false });
    this.addTab('energy', 'Energy');
    this.addTab('switch', 'Switch');
    this.addTab('details', 'Details');
    this.showTab('energy');
    // toolbar icons
    const iconFor = { fullscreen: 'fullscreen', focus: 'maximize', home: 'home', open: 'open', predict: 'flask', 'reset-view': 'reset', spin: 'spin', 'export-image': 'image', 'export-video': 'film', prev: 'prev', next: 'next', play: 'play', measure: 'ruler', ss2d: 'ss2d', console: 'terminal', undo: 'undo', redo: 'redo', clear: 'eraser' };
    for (const b of $$('[data-act]')) {
      if (!b.firstChild || b.classList.contains('icon-only')) b.replaceChildren(icon(iconFor[b.dataset.act] || 'info'));
      b.addEventListener('click', () => this.action(b.dataset.act));
    }
    // representation
    $('#sel-rep').addEventListener('change', e => this.setRepresentation(e.target.value));
    // colour scheme & style
    this.refreshSchemes();
    $('#sel-color').addEventListener('change', e => this.setColorScheme(e.target.value));
    const st = $('#sel-style');
    for (const [k, s] of Object.entries(STYLES)) st.append(el('option', { value: k, text: s.label }));
    st.value = this.viewer.settings.style;
    st.addEventListener('change', e => this.setStyle(e.target.value));
    $('#btn-bg').addEventListener('click', e => this._backgroundMenu(e.currentTarget));
    $('#tg-ghosts').addEventListener('click', () => this.toggleGhosts());
    $('#tg-putty').addEventListener('click', () => this.togglePutty());
    $('#tg-surface').addEventListener('click', () => this.surface.set(!this.surface.enabled));
    $('#btn-lang').addEventListener('click', e => this.languageMenu(e.currentTarget));
    this._syncLanguage();
    $('#btn-theme').addEventListener('click', () => this.setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'));
    $('#sel-order').value = this.state.orderMode;
    $('#sel-order').addEventListener('change', e => this.setOrder(e.target.value));
    const dur = $('#tl-dur'), hold = $('#tl-hold'), loop = $('#tl-loop');
    dur.value = this.morphDuration; hold.value = this.holdTime; loop.checked = this.loop;
    dur.addEventListener('input', () => { this.morphDuration = +dur.value; this.savePrefs(); });
    hold.addEventListener('input', () => { this.holdTime = +hold.value; this.savePrefs(); });
    loop.addEventListener('change', () => { this.loop = loop.checked; this.savePrefs(); });
    this._syncToolbar();
  }

  addTab(id, title) {
    const btn = el('button', { class: 'tab', role: 'tab', text: title, onclick: () => this.showTab(id) });
    const body = el('div', { class: 'tab-body', role: 'tabpanel', hidden: true });
    $('#tabs').append(btn);
    $('#tab-bodies').append(body);
    this._tabBodies[id] = body;
    this._tabs.push({ id, btn, body });
    return body;
  }
  showTab(id) {
    for (const t of this._tabs) { t.btn.classList.toggle('on', t.id === id); t.body.hidden = t.id !== id; t.btn.setAttribute('aria-selected', t.id === id ? 'true' : 'false'); }
    this.emit('tab', id);
  }
  addPluginTab({ id, title, render }) {
    const body = this.addTab(id, title);
    try { render(body); } catch (err) { body.append(el('div', { class: 'warnbox', text: err.message })); }
    return body;
  }

  rebuildMenus() {
    const vs = this.viewer.settings;
    const has = !!this.ens;
    const menus = [
      { label: 'File', items: () => [
        { label: 'Start page', kbd: 'Ctrl+H', action: () => this.showPage('welcome') },
        { label: 'New prediction from sequence…', kbd: 'Ctrl+N', action: () => this.action('predict') },
        { sep: true },
        { label: 'Open .rcv file or structures…', kbd: 'Ctrl+O', action: () => openDialog(this) },
        { label: 'Open prediction run…', action: () => openDialog(this, 'library'), disabled: !this.api.online },
        { label: 'Browse folders…', action: () => openDialog(this, 'server'), disabled: !this.api.online },
        { sep: true },
        { label: 'Save as .rcv file…', kbd: 'Ctrl+S', action: () => this.saveRcv(), disabled: !has || !this.api.online },
        { label: 'Save workspace (all open ensembles)…', kbd: 'Ctrl+Shift+S', action: () => this.saveWorkspace(), disabled: !has },
        { label: 'Open .rcv file or workspace…', action: () => this.openSessionFile() },
        { label: 'Fetch from the PDB…', kbd: 'Ctrl+Shift+F', action: () => this.fetchDialog(), disabled: !this.api.online },
        { sep: true },
        { label: 'Save current conformation as PDB…', action: () => this.exportPDB(false), disabled: !has },
        { label: 'Save all conformations (multi-model PDB)…', action: () => this.exportPDB(true), disabled: !has },
        { label: 'Save energy table (CSV)…', action: () => this.exportCSV(), disabled: !has },
        { sep: true },
        { label: 'Close ensemble', action: () => this.ens && this.closeEnsemble(this.ens.id), disabled: !has },
      ] },
      { label: 'Edit', items: () => [
        { label: 'Undo', kbd: 'Ctrl+Z', action: () => this.history.undo(), disabled: !this.history?.canUndo },
        { label: 'Redo', kbd: 'Ctrl+Y', action: () => this.history.redo(), disabled: !this.history?.canRedo },
        { sep: true }, { head: 'Clear' },
        ...clearItems(this),
        { sep: true }, { head: 'Selection' },
        { label: 'Select all', kbd: 'Ctrl+A', action: () => this.select([...Array(this.ens.nRes).keys()], 'set'), disabled: !has },
        { label: 'Invert selection', action: () => this.select([...Array(this.ens.nRes).keys()].filter(r => !this.state.selection.has(r)), 'set'), disabled: !has },
        { label: 'Select switching residues', action: () => this.console.run('select switching', { echo: this.console.visible }), disabled: !has },
        { label: 'Select paired residues', action: () => this.console.run('select paired', { echo: this.console.visible }), disabled: !has },
        { sep: true }, { head: 'Annotate' },
        { label: 'Label tool (click the structure)', kbd: 'T', action: () => this.tools.set('label'), disabled: !has },
        { label: 'Add a title / note on the view', action: () => { this.tools.set('label'); this.annotations.addNote('Title'); }, disabled: !has },
        { label: 'Paint tool', kbd: 'B', action: () => this.tools.set('paint'), disabled: !has },
      ] },
      { label: 'View', items: () => [
        { head: 'Representation' },
        ...[['cartoon', 'Cartoon', '1'], ['cartoon+sticks', 'Cartoon + atoms', '2'], ['sticks', 'Sticks', '3'], ['ballstick', 'Ball & stick', '4'], ['spacefill', 'Spheres', '5']]
          .map(([id, l, k]) => ({ label: l, kbd: k, checked: vs.representation === id, action: () => this.setRepresentation(id) })),
        { sep: true }, { head: 'Style' },
        ...Object.entries(STYLES).map(([id, s]) => ({ label: s.label, checked: vs.style === id, action: () => this.setStyle(id) })),
        { sep: true },
        { label: 'Perspective', checked: this.viewer.projection === 'perspective', action: () => this.setProjection('perspective') },
        { label: 'Orthographic', checked: this.viewer.projection === 'orthographic', action: () => this.setProjection('orthographic') },
        { sep: true },
        { head: 'Mouse tool' },
        ...[['select', 'Select', 'V'], ['pull', 'Pull (drag a residue to switch)', 'G'], ['paint', 'Paint colours', 'B'], ['label', 'Labels and annotations', 'T'], ['measure', 'Measure', 'M']]
          .map(([id, l, k]) => ({ label: l, kbd: k, checked: this.tool === id, action: () => this.tools.set(id), disabled: !has && id !== 'select' })),
        { sep: true },
        { label: 'Ensemble overlay', kbd: 'E', checked: this.state.ghosts, action: () => this.toggleGhosts() },
        { label: 'Flexibility putty', kbd: 'P', checked: this.state.putty, action: () => this.togglePutty() },
        { label: 'Molecular surface', kbd: 'U', checked: !!this.surface?.enabled, action: () => this.surface.set(!this.surface.enabled), disabled: !has },
        { label: 'Surface and transparency…', action: () => transparencyDialog(this), disabled: !has },
        { label: 'Measure (click two atoms)', kbd: 'M', checked: !!this.measureMode, action: () => this.action('measure'), disabled: !this.ens },
        { label: 'Auto-rotate', kbd: 'A', checked: !!this.viewer.controls.spin, action: () => this.action('spin') },
        { label: 'Reset view', kbd: 'R', action: () => this.action('reset-view') },
        { sep: true },
        { label: 'Interface language…', action: () => this.languageMenu($('#btn-lang')) },
        { label: 'Preferences…', action: () => preferencesDialog(this) },
      ] },
      { label: 'Color', items: () => [
        ...SCHEMES.map((s, i) => ({ label: s.label, checked: this.state.colorScheme === s.id, action: () => this.setColorScheme(s.id) })),
        { sep: true }, { head: 'Background' },
        ...BACKGROUNDS.map(b => ({ label: b.label, checked: vs.background === b.color, action: () => this.setBackground(b.color) })),
        { sep: true }, { head: 'Palette' },
        ...this.palettes.all().map(p => ({ label: p.name, checked: this.palettes.activeId === p.id, action: () => this.palettes.apply(p.id) })),
        { label: 'Palettes and custom colours…', action: () => this.openPalettes() },
      ] },
      { label: 'Conformation', items: () => [
        { label: 'Previous', kbd: '←', action: () => this.prev(), disabled: !has },
        { label: 'Next', kbd: '→', action: () => this.next(), disabled: !has },
        { label: this.playing ? 'Pause tour' : 'Play tour', kbd: 'Space', action: () => this.togglePlay(), disabled: !has },
        { sep: true }, { head: 'Order' },
        ...[['rank', 'Rank'], ['metric', 'By energy'], ['path', 'Similarity path'], ['input', 'Input order']].map(([id, l]) => ({ label: l, checked: this.state.orderMode === id, action: () => this.setOrder(id) })),
        { sep: true },
        { label: 'Superpose on selection', action: () => this.superposeOnSelection(), disabled: !has || this.state.selection.size < 3 },
        { label: 'Restore original superposition', action: () => this.resetSuperposition(), disabled: !has },
      ] },
      { label: 'Export', items: () => [
        { label: 'Image (PNG, TIFF, PDF)…', kbd: 'Ctrl+E', action: () => exportImageDialog(this), disabled: !has },
        { label: 'Figure panels (several conformations)…', action: () => panelsDialog(this), disabled: !has },
        { label: 'Plot studio (RMSD matrix, RMSF, energy)…', action: () => this.openPlotStudio(), disabled: !has },
        { label: 'Quick PNG snapshot', kbd: 'S', action: () => this.quickSnapshot(), disabled: !has },
        { label: 'Video (MP4, WebM, frames)…', action: () => exportVideoDialog(this), disabled: !has },
        { label: 'Transition video A → B…', action: () => exportVideoDialog(this, { content: 'transition' }), disabled: !has || this.ens.n < 2 },
        { sep: true },
        { label: 'Current conformation (PDB)…', action: () => this.exportPDB(false), disabled: !has },
        { label: 'All conformations (PDB)…', action: () => this.exportPDB(true), disabled: !has },
        { label: 'Energy table (CSV)…', action: () => this.exportCSV(), disabled: !has },
        { label: 'Interactive HTML page (share, supplement)…', action: () => this.exportHTML(), disabled: !has || !this.api.online },
        { sep: true },
        { label: 'Methods text & citation…', action: () => methodsDialog(this), disabled: !has || !this.api.online },
      ] },
      { label: 'Tools', items: () => toolsMenuItems(this) },
      { label: 'Plugins', items: () => [
        ...(this.plugins?.menuItems || []).map(m => ({ label: m.label, kbd: m.kbd, action: m.action })),
        ...((this.plugins?.menuItems || []).length ? [{ sep: true }] : []),
        { label: 'Manage plugins…', action: () => pluginsDialog(this, this.info) },
        { label: 'New plugin from template…', action: () => newPluginDialog(this), disabled: !this.api.online },
        { label: 'Open plugin folder', action: () => { const d = this.info?.plugin_dirs?.[0]; if (this.native?.available() && d) this.native.reveal(d + '/x'); else toast(d || 'unknown', 'info', 'Plugin folder', 8000); } },
        { label: 'Plugin developer guide', action: () => pluginGuideDialog(this), disabled: !this.api.online },
        { sep: true },
        { label: 'Reload window (apply plugin changes)', action: () => location.reload() },
      ] },
      { label: 'Window', items: () => this.wm.menuItems() },
      { label: 'Help', items: () => [
        { label: 'User guide', kbd: 'F1', action: () => userGuideDialog(this) },
        { label: 'Keyboard and mouse', kbd: '?', action: () => shortcutsDialog() },
        { label: 'Console commands', action: () => { this.console.setVisible(true); this.console.run('help'); } },
        { label: 'Plugin developer guide', action: () => pluginGuideDialog(this), disabled: !this.api.online },
        { label: 'Methods text & citation…', action: () => methodsDialog(this), disabled: !has || !this.api.online },
        { label: 'Open the demo ensembles', action: () => this.openDemo(), disabled: !this.api.online },
        { label: 'What is new in this version', action: () => userGuideDialog(this, 'changelog') },
        { label: 'RCV website', action: () => openWebsite(this) },
        { sep: true },
        { label: 'About', action: () => aboutDialog(this, this.info) },
      ] },
    ];
    menuBar($('#menus'), menus);
  }

  // Side panels: hidden by default below 820 px (shown as overlays on demand), collapsible above.
  _panelShown(side) { return this.wm.isShown(side); }
  _togglePanel(side) { this.wm.toggle(side); }

  // ---------------------------------------------------------------- interface language
  languageMenu(anchor) {
    popupMenu(anchor || { x: window.innerWidth - 260, y: 44 }, LANGUAGES.map(l => ({ label: l.name, checked: currentLanguage() === l.code, action: () => this.setLanguage(l.code) })));
    document.querySelector('.menu-pop')?.setAttribute('translate', 'no');
  }
  async setLanguage(code) {
    await setLanguage(code);
    this._syncLanguage();
    this.rebuildMenus();
    this._status();
    this.emit('language', code);
  }
  _syncLanguage() {
    for (const b of $$('.lang-btn')) {
      b.replaceChildren(icon('globe'), el('span', { class: 'lang-name', text: languageName(currentLanguage()) }));
      b.setAttribute('translate', 'no');
    }
  }

  refreshSchemes() {
    const sel = $('#sel-color');
    sel.replaceChildren(...SCHEMES.map(s => el('option', { value: s.id, text: s.label })));
    sel.value = this.state.colorScheme;
  }

  _syncToolbar() {
    const vs = this.viewer.settings;
    $('#sel-rep').value = vs.representation;
    $('#sel-color').value = this.state.colorScheme;
    $('#sel-style').value = vs.style;
    $('#btn-bg span').style.background = vs.background;
    $('#tg-ghosts').classList.toggle('on', this.state.ghosts);
    $('#tg-putty').classList.toggle('on', this.state.putty);
    $('#tg-surface')?.classList.toggle('on', !!this.surface?.enabled);
    $('#tg-pull')?.classList.toggle('on', !!this.pullCtl?.mode);
    $('#tg-measure')?.classList.toggle('on', !!this.measureMode);
    $('#tg-ss2d')?.classList.toggle('on', !!this.ss2d?.visible);
    $('#tg-console')?.classList.toggle('on', !!this.console?.visible);
    $('#tg-focus')?.classList.toggle('on', !!this.wm?.focus);
    $('#tg-fullscreen')?.classList.toggle('on', !!this.wm?.fullscreen);
    $('#sel-order').value = this.state.orderMode;
    const spin = $('[data-act="spin"]');
    spin?.classList.toggle('on', !!this.viewer.controls.spin);
    const play = $('[data-act="play"]');
    play.replaceChildren(icon(this.playing ? 'pause' : 'play'));
    play.title = this.playing ? 'Pause (Space)' : 'Play conformation tour (Space)';
    const dark = document.documentElement.dataset.theme === 'dark';
    $('#btn-theme').replaceChildren(icon(dark ? 'sun' : 'moon'));
  }

  _status() {
    const sb = $('#statusbar');
    const ens = this.ens;
    sb.replaceChildren(...[
      el('span', {}, el('span', { class: `dot ${this.api.online ? '' : 'off'}` }), this.api.online ? `Running locally · ${this.info?.system?.cpus || '?'} CPU cores` : 'Engine not reachable'),
      ens ? el('span', { text: `${ens.name} · ${ens.nRes} nt · ${ens.nAtoms} atoms · ${ens.n} conformation${ens.n === 1 ? '' : 's'}` }) : null,
      ens && this.state.selection.size ? el('span', { text: `${this.state.selection.size} selected` }) : null,
      el('span', { class: 'grow' }),
      el('span', { text: 'Drag rotate · right-drag pan · wheel zoom · click select · ←/→ morph' })].filter(Boolean));
  }

  // With no ensemble open, the start page is the place to be.
  _emptyState() {
    $('#empty').hidden = true;
    if (!this.ens && this.welcome && this.page === 'viewer') this.showPage('welcome');
  }

  // ---------------------------------------------------------------- documents
  async loadServerDocument(id, { select = true, session = null } = {}) {
    const doc = await this.api.get(`/api/documents/${id}`);
    const ens = this.addDocument(doc, { serverId: id, select });
    if (select && session?.view) this.applyViewState(session.view);
    return ens;
  }

  addDocument(doc, { serverId = null, select = true } = {}) {
    const ens = new Ensemble(doc, { serverId });
    for (const m of this._clientMetrics || []) this._computeClientMetric(ens, m);
    this.ensembles.set(ens.id, ens);
    this.emit('documents');
    if (select) { this.selectEnsemble(ens.id); this.showPage('viewer'); }
    return ens;
  }

  selectEnsemble(id, { keepView = false, conf = null, selection = null } = {}) {
    const ens = this.ensembles.get(id);
    if (!ens) return;
    this.stop();
    this.morph = null;
    this.ens = ens;
    const ms = ens.metricByKey;
    // each ensemble remembers its metric; first visit uses the document's default (engine score for predictions)
    const want = ens.lastMetric || ens.doc.default_metric;
    this.state.metric = want in ms ? want : ens.metrics[0]?.key;
    this.state.selection = new Set(selection || []);
    this._computeOrder();
    this.state.conf = conf ?? this.order[0] ?? 0;
    this.state.prev = null;
    this.display = ens.coords[this.state.conf].slice();
    this._colorCache.clear();
    this.viewer.setEnsemble(ens, this.display, { keepView });
    const disp = this.displayOf(ens);
    this.viewer.setDisplaySets({ hidden: disp.hidden, sticks: disp.sticks, spheres: disp.spheres, labels: disp.labels, atomColors: disp.atomColors });
    if (selection) this.viewer.setSelection(selection);
    this.viewer.setRepresentation(this.viewer.settings.representation);
    this.viewer.setStyle(this.viewer.settings.style);
    this._applyGhosts();
    this._applyPutty();
    this.refreshColors();
    this._emptyState();
    this._hud();
    this._status();
    this.rebuildMenus();
    this.emit('ensemble', ens);
    this.emit('documents');
  }

  async closeEnsemble(id) {
    const ens = this.ensembles.get(id);
    if (!ens) return;
    this.ensembles.delete(id);
    if (ens.serverId && this.api.online) this.api.del(`/api/documents/${ens.serverId}`).catch(() => {});
    if (this.ens === ens) {
      const next = [...this.ensembles.keys()][0];
      if (next) this.selectEnsemble(next);
      else { this.ens = null; this.viewer.setEnsemble(null); this._emptyState(); this._hud(); this._status(); this.emit('ensemble'); }
    }
    this.emit('documents');
  }

  async openServerPaths(paths, { name, chain } = {}) {
    const rcv = paths.length === 1 && /\.rcv$/i.test(paths[0]);
    const done = busy(rcv ? 'Opening the .rcv file…' : 'Reading structures and computing energies…');
    try {
      const s = await this.api.post('/api/open', { paths, name, chain });
      await this.loadServerDocument(s.id, { session: s.session });
      if (!rcv) this._warnings();
    } catch (err) { toast(err.message, 'error', 'Could not open'); } finally { done(); }
  }

  async openRcvFile(file) {
    if (!this.api.online) { toast('Opening .rcv files needs the local RCV engine.', 'error', 'Offline'); return; }
    const done = busy(`Opening ${file.name}…`);
    try {
      const r = await this.api.raw('POST', '/api/rcv/open', await file.arrayBuffer(), { 'Content-Type': 'application/octet-stream', 'X-Name': encodeURIComponent(file.name) });
      const s = await r.json();
      await this.loadServerDocument(s.id, { session: s.session });
    } catch (err) { toast(`${file.name}: ${err.message}`, 'error', 'Could not open'); } finally { done(); }
  }

  async openLocalFiles(files, { name, chain } = {}) {
    if (!files.length) return;
    const sessions = files.filter(f => /\.json$/i.test(f.name) && !/manifest\.json$/i.test(f.name));
    for (const f of sessions) await this.loadSessionText(await readFileText(f), f.name);
    const rcvs = files.filter(f => /\.rcv$/i.test(f.name));
    for (const f of rcvs) await this.openRcvFile(f);
    const rest = files.filter(f => !sessions.includes(f) && !rcvs.includes(f));
    if (!rest.length) return;
    if (!this.api.online) { toast('Reading PDB/mmCIF files needs the local viewer server. Sessions open offline.', 'error', 'Offline'); return; }
    const done = busy(`Reading ${rest.length} file(s) and computing energies…`);
    try {
      const payload = [];
      for (const f of rest) payload.push({ name: f.webkitRelativePath || f.relativePath || f.name, text: await readFileText(f) });
      const s = await this.api.post('/api/upload', { files: payload, name, chain });
      await this.loadServerDocument(s.id);
      this._warnings();
    } catch (err) { toast(err.message, 'error', 'Could not open'); } finally { done(); }
  }

  async openDropped(dt, opts) {
    const files = [];
    const walk = async (entry, prefix = '') => {
      if (entry.isFile) {
        const f = await new Promise((res, rej) => entry.file(res, rej));
        Object.defineProperty(f, 'relativePath', { value: prefix + f.name });
        files.push(f);
      } else if (entry.isDirectory) {
        const reader = entry.createReader();
        let batch;
        do {
          batch = await new Promise((res, rej) => reader.readEntries(res, rej));
          for (const e of batch) await walk(e, prefix + entry.name + '/');
        } while (batch.length);
      }
    };
    const items = [...(dt.items || [])].map(i => i.webkitGetAsEntry?.()).filter(Boolean);
    if (items.length) for (const e of items) await walk(e);
    else files.push(...dt.files);
    const keep = files.filter(f => /\.(pdb|ent|cif|mmcif|json|rcv)$/i.test(f.name));
    const fasta = files.filter(f => /\.(fa|fasta|fna)$/i.test(f.name));
    if (fasta.length && !keep.length) { this.showPage('welcome'); await this.welcome.loadFastaFile(fasta[0]); return; }
    if (!keep.length) { toast('No .rcv, PDB, mmCIF or workspace files found in the drop.', 'error'); return; }
    await this.openLocalFiles(keep, opts);
  }

  _warnings() {
    for (const w of this.ens?.doc.warnings || []) toast(w, 'info', 'Note');
  }

  // ---------------------------------------------------------------- order, metric, colours
  _computeOrder() {
    this.order = this.ens ? this.ens.order(this.state.orderMode, this.state.metric) : [];
  }
  setOrder(mode) {
    this.state.orderMode = mode;
    this._computeOrder();
    this.savePrefs();
    this._syncToolbar();
    this.emit('order');
  }
  setMetric(key) {
    if (!this.ens || !this.ens.metricByKey[key]) return;
    this.state.metric = key;
    this.ens.lastMetric = key;
    if (this.state.orderMode === 'metric') this._computeOrder();
    this._hud();
    this.emit('metric');
  }

  // Per-ensemble display state set from the console / selections panel: hidden residues, extra sticks and
  // spheres, custom labels, custom colours and named selections.
  displayOf(ens = this.ens) {
    const fresh = () => ({ hidden: new Set(), sticks: new Set(), spheres: new Set(), labels: new Map(), colors: new Map(), atomColors: new Map(), named: new Map(), version: 0 });
    if (!ens) return fresh();
    const d = (ens.display ||= fresh());
    d.atomColors ||= new Map();
    return d;
  }
  applyDisplay() {
    if (!this.ens) return;
    const d = this.displayOf();
    d.version++;
    this.viewer.setDisplaySets({ hidden: d.hidden, sticks: d.sticks, spheres: d.spheres, labels: d.labels, atomColors: d.atomColors });
    this._colorCache.clear();
    this.refreshColors();
    this.emit('display');
  }

  colorsFor(k, prevK = null) {
    const ens = this.ens, scheme = this.state.colorScheme;
    const ov = ens.display?.colors;
    const key = `${ens.id}|${scheme}|${k}|${scheme === 'displacement' ? prevK : ''}|${scheme === 'uniform' ? this.state.uniformColor : ''}|${ov && ov.size ? ens.display.version : ''}`;
    if (this._colorCache.has(key)) return this._colorCache.get(key);
    const ctx = { conf: k, coords: ens.coords[k], prevCoords: prevK !== null ? ens.coords[prevK] : null, uniform: this.state.uniformColor };
    const hex = residueColors(ens, scheme, ctx);
    if (ov && ov.size) for (const [r, h] of ov) if (r < hex.length) hex[r] = h;
    const lin = new Float32Array(ens.nRes * 3);
    hex.forEach((h, i) => { const c = hexToLin(h); lin[3 * i] = c[0]; lin[3 * i + 1] = c[1]; lin[3 * i + 2] = c[2]; });
    const out = { hex, lin, ctx };
    if (this._colorCache.size > 400) this._colorCache.clear();
    this._colorCache.set(key, out);
    return out;
  }

  refreshColors() {
    if (!this.ens) return;
    const c = this.colorsFor(this.state.conf, this.state.prev);
    this.curColors = c.lin.slice();
    this.viewer.setResidueColorsLinear(this.curColors);
    this._legend();
    this.emit('colors');
  }

  legendData() {
    if (!this.ens) return null;
    const c = this.colorsFor(this.state.conf, this.state.prev);
    return schemeLegend(this.ens, this.state.colorScheme, c.ctx);
  }

  _legend() {
    const box = $('#hud-legend');
    box.replaceChildren();
    const lg = this.legendData();
    if (!lg) return;
    box.append(el('div', { class: 'lg-title', text: lg.title }));
    if (lg.ramp) {
      box.append(el('span', { class: 'lg', text: lg.min }), el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${lg.ramp.join(',')})` } }), el('span', { class: 'lg', text: lg.max }));
    } else for (const it of lg.items) box.append(el('span', { class: 'lg' }, el('span', { class: 'sw', style: { background: it.color } }), it.label));
  }

  setColorScheme(id) {
    this.state.colorScheme = id;
    this._colorCache.clear();
    this.refreshColors();
    this.savePrefs();
    this._syncToolbar();
    this.emit('settings');
  }

  setRepresentation(r) {
    this.viewer.setRepresentation(r);
    this.savePrefs();
    this._syncToolbar();
  }
  setStyle(s) {
    this.viewer.setStyle(s);
    this.savePrefs();
    this._syncToolbar();
  }
  setBackground(color) {
    this.viewer.settings.background = color;
    this.viewer.needsRender = true;
    this._applyBackground();
    this.savePrefs();
    this._syncToolbar();
  }
  setProjection(p) { this.viewer.setProjection(p); this.prefs.projection = p; this.savePrefs(); }
  setLabelMode(m) { this.viewer.labelMode = m; this.viewer.needsRender = true; this.savePrefs(); }
  isDarkColor(hex) { return luminance(hex) < 0.25; }
  isDarkBackground() { return this.isDarkColor(this.viewer.settings.background); }
  _applyBackground() {
    const vp = $('#viewport');
    const dark = this.isDarkBackground();
    vp.style.background = this.viewer.settings.background;
    vp.style.color = dark ? '#f0efec' : '#141414';
    vp.style.setProperty('--hud-sub', dark ? '#c3c2b7' : '#52514e');
  }
  _backgroundMenu(btn) {
    popupMenu(btn, [
      ...BACKGROUNDS.map(b => ({ label: b.label, checked: this.viewer.settings.background === b.color, action: () => this.setBackground(b.color) })),
      { sep: true },
      { label: 'Custom colour…', action: () => {
        const inp = el('input', { type: 'color', value: this.viewer.settings.background });
        inp.addEventListener('input', () => this.setBackground(inp.value));
        inp.click();
      } },
    ]);
  }

  setTheme(mode, save = true) {
    const sys = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
    document.documentElement.dataset.theme = mode === 'system' ? sys : mode;
    if (save) { this.prefs.theme = mode; savePref('prefs', this.prefs); }
    if (this.viewer) { this._syncToolbar(); this.emit('settings'); this.emit('metric'); }
  }

  toggleGhosts() { this.state.ghosts = !this.state.ghosts; this._applyGhosts(); this._syncToolbar(); }
  togglePull() { this.tools.set(this.tool === 'pull' ? 'select' : 'pull'); this._syncToolbar(); if (this.pullCtl.mode) toast('Drag any residue: the structure stretches, and on release it switches to the nearest stable conformation. Press G again to rotate. Alt+drag pulls at any time.', 'info', 'Pull tool', 6000); }
  _applyGhosts() {
    if (!this.ens) return;
    this.viewer.fitAll = this.state.ghosts;
    this.viewer.setGhosts(this.state.ghosts ? this.ens.coords.filter((_, k) => k !== this.state.conf) : null);
  }
  togglePutty() { this.state.putty = !this.state.putty; this._applyPutty(); this._syncToolbar(); }
  _applyPutty() {
    if (!this.ens) return;
    if (!this.state.putty) { this.viewer.setPutty(null); return; }
    const r = this.ens.rmsf, max = Math.max(1e-6, ...r);
    this.viewer.setPutty(r.map(v => 0.35 + 2.2 * v / max));
  }

  // ---------------------------------------------------------------- conformation navigation & morph
  // target: optional coordinates to land on instead of ens.coords[k] (same conformation, e.g. re-fitted by the pull tool)
  goTo(k, { animate = true, duration = null, force = false, target: landing = null } = {}) {
    const ens = this.ens;
    if (!ens || k === null || k === undefined || k < 0 || k >= ens.n) return;
    if (k === this.state.conf && !this.morph && !force) return;
    const from = this.state.conf;
    const startCoords = this.display.slice();
    const startColors = (this.curColors || this.colorsFor(from).lin).slice();
    if (from !== k) this.state.prev = from;
    this.state.conf = k;
    const target = this.colorsFor(k, this.state.prev);
    const dur = (duration ?? this.morphDuration) * 1000;
    this.viewer.setPreview(null);
    this.emit('conformation', { from, to: k });
    this._hud();
    if (!animate || dur <= 0) {
      this.morph = null;
      this.display.set(landing || ens.coords[k]);
      this.viewer.setCoords(this.display);
      this.curColors = target.lin.slice();
      this.viewer.setResidueColorsLinear(this.curColors);
      this._settled();
      return;
    }
    this.morph = { from, to: k, t0: performance.now(), dur, start: startCoords, target: landing || ens.coords[k], cs: startColors, ct: target.lin };
    if (!this._morphAnimator) {
      this._morphAnimator = true;
      this.viewer._addAnimator(() => { const go = this._stepMorph(); if (!go) this._morphAnimator = false; return go; });
    }
  }

  _stepMorph() {
    const m = this.morph;
    if (!m) return false;
    const t = clamp((performance.now() - m.t0) / m.dur, 0, 1);
    const e = easeInOut(t);
    const d = this.display, a = m.start, b = m.target;
    for (let i = 0; i < d.length; i++) d[i] = a[i] + (b[i] - a[i]) * e;
    this.viewer.setCoords(d);
    const cc = this.curColors = this.curColors && this.curColors.length === m.cs.length ? this.curColors : new Float32Array(m.cs.length);
    for (let i = 0; i < cc.length; i++) cc[i] = m.cs[i] + (m.ct[i] - m.cs[i]) * e;
    this.viewer.setResidueColorsLinear(cc);
    this.emit('morph', { from: m.from, to: m.to, t: e });
    $('#hud-morph').hidden = false;
    $('#hud-morph').textContent = `${this.ens.shortId(m.from)} → ${this.ens.shortId(m.to)}`;
    if (t >= 1) { this.morph = null; this._settled(); return false; }
    return true;
  }

  _settled() {
    $('#hud-morph').hidden = true;
    this.emit('morph', { done: true });
    if (this.state.colorScheme === 'displacement') this.refreshColors(); else this._legend();
    if (this.state.ghosts) this._applyGhosts();
    this._hud();
    if (this.playing) {
      clearTimeout(this._playTimer);
      this._playTimer = setTimeout(() => this._playStep(), this.holdTime * 1000);
    }
  }

  next() { this._step(1); }
  prev() { this._step(-1); }
  _step(dir) {
    if (!this.ens) return;
    const i = this.order.indexOf(this.state.conf);
    let j = i + dir;
    if (j >= this.order.length) j = this.loop ? 0 : this.order.length - 1;
    if (j < 0) j = this.loop ? this.order.length - 1 : 0;
    if (j !== i) this.goTo(this.order[j]);
  }
  first() { if (this.ens) this.goTo(this.order[0]); }
  last() { if (this.ens) this.goTo(this.order[this.order.length - 1]); }

  togglePlay() { if (this.playing) this.stop(); else this.play(); }
  play() {
    if (!this.ens || this.ens.n < 2) return;
    this.playing = true;
    this._syncToolbar();
    this._playStep();
  }
  _playStep() {
    if (!this.playing) return;
    const i = this.order.indexOf(this.state.conf);
    if (i === this.order.length - 1 && !this.loop) { this.stop(); return; }
    this.next();
  }
  stop() {
    this.playing = false;
    clearTimeout(this._playTimer);
    if (this.viewer) this._syncToolbar();
  }

  // Continuous scrubbing along the order (timeline drag): pos in [0, n-1]
  scrub(pos) {
    const ens = this.ens;
    if (!ens) return;
    this.stop();
    this.morph = null;
    const n = this.order.length;
    pos = clamp(pos, 0, n - 1);
    const i = Math.floor(pos), j = Math.min(n - 1, i + 1), t = pos - i;
    const a = this.order[i], b = this.order[j];
    this.applyFrame(a, b, t);
    const near = this.order[Math.round(pos)];
    if (near !== this.state.conf) {
      this.state.prev = this.state.conf;
      this.state.conf = near;
      this.emit('conformation', { from: this.state.prev, to: near, scrub: true });
      this._hud();
    }
    this.emit('morph', { from: a, to: b, t });
  }
  endScrub(pos) {
    if (!this.ens) return;
    const k = this.order[Math.round(clamp(pos, 0, this.order.length - 1))];
    if (k !== this.state.conf) { this.state.prev = this.state.conf; this.state.conf = k; }
    this.goTo(k, { duration: 0.25, force: true });     // settle from the interpolated frame onto k
  }

  // Set the displayed state to an interpolation between two conformations (used by scrub and video).
  applyFrame(a, b, t) {
    const ens = this.ens;
    const A = ens.coords[a], B = ens.coords[b], d = this.display;
    for (let i = 0; i < d.length; i++) d[i] = A[i] + (B[i] - A[i]) * t;
    this.viewer.setCoords(d);
    if (this.capturing && this.surface?.enabled) this.surface.update(d, { spacing: 1.1 });
    const ca = this.colorsFor(a).lin, cb = this.colorsFor(b, a).lin;
    const cc = this.curColors = new Float32Array(ca.length);
    for (let i = 0; i < cc.length; i++) cc[i] = ca[i] + (cb[i] - ca[i]) * t;
    this.viewer.setResidueColorsLinear(cc);
  }

  snapshotView() {
    return { conf: this.state.conf, prev: this.state.prev, display: this.display.slice(), colors: this.curColors?.slice(), q: this.viewer.controls.quaternion.clone(), playing: this.playing };
  }
  restoreView(s) {
    this.display.set(s.display);
    this.viewer.setCoords(this.display);
    if (s.colors) { this.curColors = s.colors; this.viewer.setResidueColorsLinear(this.curColors); }
    this.viewer.controls.quaternion.copy(s.q);
    this.viewer.needsRender = true;
  }

  previewConformation(k) {
    clearTimeout(this._prevTimer);
    if (k === null || k === undefined || !this.ens || k === this.state.conf || this.morph) { this.viewer.setPreview(null); return; }
    this._prevTimer = setTimeout(() => this.viewer.setPreview(this.ens.coords[k]), 60);
  }

  // ---------------------------------------------------------------- selection & picking
  select(residues, mode = 'set') {
    const s = this.state.selection;
    if (mode === 'set') { s.clear(); residues.forEach(r => s.add(r)); }
    else if (mode === 'add') residues.forEach(r => s.add(r));
    else if (mode === 'toggle') residues.forEach(r => (s.has(r) ? s.delete(r) : s.add(r)));
    this.viewer.setSelection([...s]);
    this._status();
    this.rebuildMenus();
    this.emit('selection', [...s]);
  }
  hoverResidues(list, opts = {}) {
    this.viewer.setHighlight(list);
    if (!opts.silent) this.emit('hover', list || []);
  }
  centerOnSelection() {
    const s = [...this.state.selection];
    if (!s.length || !this.ens) return;
    const pts = s.map(r => this.viewer.labelAnchor(r));
    const c = [0, 1, 2].map(i => pts.reduce((a, p) => a + p[i], 0) / pts.length);
    this.viewer.centerOn(c);
  }
  superposeOnSelection() {
    try {
      this.ens.superposeOn([...this.state.selection], this.state.conf);
      this.ens.superposedOn = [...this.state.selection];
      this.display.set(this.ens.coords[this.state.conf]);
      this.viewer.setCoords(this.display);
      this._applyGhosts();
      toast(`All conformations fitted on ${this.state.selection.size} selected residues.`, 'ok', 'Superposition');
    } catch (err) { toast(err.message, 'error', 'Superposition'); }
  }
  resetSuperposition() {
    const ens = this.ens;
    const fresh = new Ensemble(ens.doc, { serverId: ens.serverId });
    ens.coords = fresh.coords;
    ens.superposedOn = null;
    ens._bounds();
    this.display.set(ens.coords[this.state.conf]);
    this.viewer.setCoords(this.display);
    this._applyGhosts();
  }

  toggleMeasureMode(on = !this.measureMode, { quiet = false } = {}) {
    if (!quiet && this.tools && (on ? this.tool !== 'measure' : this.tool === 'measure')) { this.tools.set(on ? 'measure' : 'select'); return; }
    this.measureMode = on;
    this._measureAtoms = [];
    this.viewer.removeOverlay('measure-pending');
    this._syncToolbar();
    this.emit('measurements');
    if (on && !quiet) toast('Click two atoms to measure their distance (more in the console: angle, dihedral). Press M or Esc to stop.', 'info', 'Measure', 5000);
  }

  _onPick(p, e) {
    if (!this.ens) return;
    if (this.measureMode) {
      if (!p) return;
      const atoms = this._measureAtoms = [...(this._measureAtoms || []), p.atom];
      this.viewer.addOverlay('measure-pending', ({ coords }) => atoms.map(a => ({ type: 'point', p: [coords[3 * a], coords[3 * a + 1], coords[3 * a + 2]], color: '#eb6834', r: 5 })));
      if (atoms.length === 2) {
        try {
          const m = this.measure.add('distance', atoms);
          toast(`${this.measure.describe(m)} = ${fmt(this.measure.value(m, this.display), 2)} Å (plotted per conformation as "${m.key}")`, 'ok', 'Distance', 5000);
        } catch (err) { toast(err.message, 'error'); }
        this._measureAtoms = [];
        this.viewer.removeOverlay('measure-pending');
      }
      return;
    }
    if (!p) { if (!e.shiftKey && this.state.selection.size) this.select([], 'set'); return; }
    this.select([p.residue], e.shiftKey || e.ctrlKey || e.metaKey ? 'toggle' : (this.state.selection.size === 1 && this.state.selection.has(p.residue) ? 'toggle' : 'set'));
  }
  // Right-click (without dragging; right-drag pans) on the structure: quick actions for the residue or the
  // selection it belongs to. Every action is a console command, so it is logged with the ensemble.
  _bindContextMenu() {
    const cv = this.viewer.canvas;
    let down = null;
    cv.addEventListener('pointerdown', e => { if (e.button === 2) down = { x: e.clientX, y: e.clientY }; });
    cv.addEventListener('pointerup', e => {
      if (e.button !== 2 || !down) return;
      const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
      down = null;
      if (moved > 4 || !this.ens || this.pullCtl?.dragging) return;
      const p = this.viewer.pick(e.clientX, e.clientY);
      if (p) this._residueMenu(p, e.clientX, e.clientY);
    });
  }

  _residueMenu(p, x, y) {
    const ens = this.ens, r = p.residue;
    const sel = this.state.selection;
    const onSel = sel.has(r) && sel.size > 1;
    const target = onSel ? 'sele' : `pos ${r + 1}`;
    const what = onSel ? `${sel.size} selected residues` : ens.resLabel(r);
    const run = cmd => this.console.run(cmd, { echo: this.console.visible });
    const pi = ens.pairInfo(this.state.conf, r);
    const hel = (ens.conformations[this.state.conf].helix_of_residue || [])[r];
    const fam = hel >= 0 ? ens.doc.helix_families?.[hel] : null;
    popupMenu({ x, y }, [
      { head: `${ens.resLabel(r)}${pi ? ' · paired with ' + ens.resLabel(pi.partner) : ' · unpaired'} · ${ens.atomNames[p.atom]}` },
      { label: `Select ${ens.resLabel(r)}`, action: () => run(`select pos ${r + 1}`) },
      { label: 'Add to selection', action: () => this.select([r], 'add') },
      pi ? { label: `Select the pair (${ens.resLabel(r)}–${ens.resLabel(pi.partner)})`, action: () => run(`select pos ${r + 1}+${pi.partner + 1}`) } : null,
      fam ? { label: `Select helix ${fam.name}`, action: () => run(`select helix ${fam.name}`) } : null,
      { label: 'Select residues within 6 Å', action: () => run(`select within 6 of pos ${r + 1}`) },
      { sep: true }, { head: `Display — ${what}` },
      { label: 'Show as sticks', action: () => run(`show sticks, ${target}`) },
      { label: 'Show as spheres', action: () => run(`show spheres, ${target}`) },
      { label: 'Hide cartoon', action: () => run(`hide cartoon, ${target}`) },
      { label: 'Label', action: () => run(`label ${target}`) },
      { label: 'Show everything again', action: () => run('show everything') },
      { sep: true }, { head: `Colour — ${what}` },
      ...['orange', 'purple', 'green', 'red', 'gold'].map(c => ({ label: c[0].toUpperCase() + c.slice(1), action: () => run(`color ${c}, ${target}`) })),
      { label: 'Scheme colour (remove custom)', action: () => run(`color default, ${target}`) },
      { sep: true },
      { label: 'Zoom here', action: () => run(`zoom ${target}`) },
      { label: `Measure from ${ens.resLabel(r)}/${ens.atomNames[p.atom]}…`, action: () => { this.toggleMeasureMode(true); this._onPick(p, { shiftKey: false }); } },
    ].filter(Boolean));
  }

  _onHover(p, e) {
    const res = p && this.ens ? [p.residue] : [];
    if (res.join() !== (this._hoverRes || []).join()) { this._hoverRes = res; this.emit('hover', res); }
    if (!p || !this.ens || !e) { tooltip.hide(); return; }
    const ens = this.ens, r = p.residue, k = this.state.conf;
    const pi = ens.pairInfo(k, r);
    const atomRep = this.viewer.settings.representation !== 'cartoon';
    tooltip.show([
      tipRow(ens.resLabel(r), atomRep ? `atom ${ens.atomNames[p.atom]}` : (ens.residues[r].comp !== ens.residues[r].name ? ens.residues[r].comp : '')),
      el('div', { class: 'tl', text: pi ? `paired with ${ens.resLabel(pi.partner)}${pi.pk ? ' (pseudoknot)' : ''}` : 'unpaired in this conformation' }),
      el('div', { class: 'tl', text: `RMSF ${fmt(ens.rmsf[r], 1)} Å across ${ens.n} conformations` }),
    ], e.clientX, e.clientY);
  }

  // ---------------------------------------------------------------- HUD
  _hud() {
    const t = $('#hud-title');
    const ens = this.ens;
    if (!ens) { t.replaceChildren(); $('#hud-legend').replaceChildren(); return; }
    const k = this.state.conf, m = ens.metricByKey[this.state.metric] || {};
    const v = ens.value(k, this.state.metric);
    t.replaceChildren(el('div', { class: 't1', text: ens.name }),
      el('div', { class: 't2' }, `${ens.shortId(k)} · `, el('span', { text: `${this.order.indexOf(k) + 1} of ${ens.n}` }),
        v !== null ? [' · ', el('span', { text: m.label }), ` ${fmtVal(v, m.unit)}${m.unit ? ' ' + m.unit : ''}`] : null));
  }

  // ---------------------------------------------------------------- client metrics (plugins)
  addClientMetric(def) {
    if (typeof def.compute !== 'function') throw new Error('metric needs compute(ensemble, k)');
    this._clientMetrics = this._clientMetrics || [];
    this._clientMetrics.push(def);
    for (const ens of this.ensembles.values()) this._computeClientMetric(ens, def);
    this.emit('values');
  }
  _computeClientMetric(ens, def) {
    for (let k = 0; k < ens.n; k++) {
      let v = null;
      try { v = def.compute(ens, k); } catch (err) { console.error(err); }
      ens.conformations[k].values[def.key] = typeof v === 'number' && Number.isFinite(v) ? v : null;
    }
    if (!ens.metricByKey[def.key]) {
      const m = { key: def.key, label: def.label || def.key, unit: def.unit || '', group: 'plugin', lower_is_better: def.lower_is_better ?? null, description: def.description || `From plugin ${def.plugin}` };
      ens.metrics.push(m);
      ens.metricByKey[def.key] = m;
    }
  }

  // ---------------------------------------------------------------- actions
  action(name) {
    switch (name) {
      case 'open': return openDialog(this);
      case 'predict': this.showPage('welcome'); setTimeout(() => this.welcome.seqIn.focus(), 50); return;
      case 'home': return this.showPage('welcome');
      case 'pull': return this.togglePull();
      case 'undo': return this.history.undo();
      case 'redo': return this.history.redo();
      case 'clear': return popupMenu($('#btn-clear'), clearItems(this));
      case 'measure': return this.ens && this.tools.set(this.tool === 'measure' ? 'select' : 'measure');
      case 'ss2d': this.ss2d.setVisible(!this.ss2d.visible); this._syncToolbar(); return;
      case 'console': this.console.toggle(); this._syncToolbar(); return;
      case 'fullscreen': return this.wm.toggleFullscreen();
      case 'focus': return this.wm.setFocus(!this.wm.focus);
      case 'reset-view': this.viewer.resetView(); return;
      case 'spin': this.viewer.controls.spin = this.viewer.controls.spin ? 0 : 0.5; this.viewer.needsRender = true; this._syncToolbar(); return;
      case 'export-image': return this.ens && exportImageDialog(this);
      case 'export-video': return this.ens && exportVideoDialog(this);
      case 'prev': return this.prev();
      case 'next': return this.next();
      case 'play': return this.togglePlay();
    }
  }

  async quickSnapshot() {
    if (!this.ens) return;
    const [w, h] = this.viewer.cssSize;
    const r = await exportFigure(this, { format: 'png', width: w * 2, height: h * 2, dpi: 300, transparent: false, background: null, legend: true, labels: true, panel: false });
    const saved = await saveBlob(r.blob, `${slug(this.ens.name)}_${this.ens.shortId(this.state.conf)}.png`);
    if (saved) toast(saved, 'ok', 'Snapshot saved');
  }

  exportTransition(a, b) { if (this.ens && this.ens.n > 1) exportVideoDialog(this, { content: 'transition', from: a, to: b }); }
  openVideoDialog(preset = {}) { if (this.ens) exportVideoDialog(this, preset); }
  openPanels() { if (this.ens) panelsDialog(this); }
  openPalettes() { paletteDialog(this); }
  openPluginGuide() { return pluginGuideDialog(this); }
  openUserGuide(which = 'guide') { return userGuideDialog(this, which); }
  setTransparency(key, v, { quiet = false } = {}) {
    const keys = key === 'all' ? ['cartoon', 'bases', 'atoms', 'surface'] : [key];
    for (const k of keys) this.viewer.setTransparency(k, v);
    if (!quiet) this.emit('surface');
  }
  openPlotStudio(type = null) { if (this.ens) plotStudio(this, type); }
  fetchDialog() { fetchDialog(this); }

  // One self-contained web page with this ensemble (viewer, morphs, energies, 2D) — no install needed to open it.
  async exportHTML() {
    const ens = this.ens;
    if (!ens) return;
    const done = busy('Writing the interactive page…');
    let blob;
    try {
      const view = this.viewState();
      delete view.parts;
      const body = { view, title: ens.name, ...(ens.serverId ? { document: ens.serverId } : { doc: ens.doc }) };
      const r = await this.api.raw('POST', '/api/export-html', JSON.stringify(body), { 'Content-Type': 'application/json' });
      blob = new Blob([await r.arrayBuffer()], { type: 'text/html' });
    } catch (err) { toast(err.message, 'error', 'HTML export'); return; } finally { done(); }
    const saved = await saveBlob(blob, `${slug(ens.name)}.html`, [{ description: 'Web page', accept: { 'text/html': ['.html'] } }]);
    if (saved) toast(`${saved} · ${(blob.size / 1e6).toFixed(1)} MB · opens in any current browser, no installation (supplementary material, reviewers).`, 'ok', 'Interactive page saved', 7000);
  }

  // Prepare a prediction of another sequence (e.g. a mutant from the switch analysis) on the start page.
  predictSequence(seq, name) {
    this.showPage('welcome');
    this.welcome.seqIn.value = `>${name}\n${seq}`;
    this.welcome.nameIn.value = name;
    this.welcome._estimate();
    this.welcome.seqIn.scrollIntoView({ block: 'center' });
    toast(`${name} is ready to predict: check the settings and click Predict.`, 'info', 'Mutant', 6000);
  }

  async exportPDB(all, which = null) {
    const ens = this.ens;
    const remarks = [`RCV EXPORT ${new Date().toISOString().slice(0, 19)}`, `ENSEMBLE ${ens.name}`, 'COORDINATES SUPERPOSED AS DISPLAYED'];
    let text = '';
    if (all) {
      text = remarks.map(r => `REMARK 250 ${r}`).join('\n') + '\n';
      this.order.forEach((k, i) => {
        text += `REMARK 250 MODEL ${i + 1} = ${ens.shortId(k)} ${this.state.metric} ${fmt(ens.value(k, this.state.metric))}\n`;
      });
      this.order.forEach((k, i) => { text += toPDB(ens, ens.coords[k], { model: i + 1 }); });
      text += 'END\n';
    } else {
      const k = which ?? this.state.conf;
      text = toPDB(ens, this.morph && which === null ? this.display : ens.coords[k], { remarks: [...remarks, `CONFORMATION ${ens.shortId(k)}`] }) + 'END\n';
    }
    const name = `${slug(ens.name)}_${all ? 'all' : ens.shortId(which ?? this.state.conf)}.pdb`;
    const saved = await saveBlob(new Blob([text], { type: 'chemical/x-pdb' }), name, [{ description: 'PDB file', accept: { 'chemical/x-pdb': ['.pdb'] } }]);
    if (saved) toast(saved, 'ok', 'PDB saved');
  }

  async exportCSV() {
    const ens = this.ens;
    const ms = ens.metrics.filter(m => ens.conformations.some(c => typeof c.values[m.key] === 'number'));
    const q = s => `"${String(s).replace(/"/g, '""')}"`;
    const rows = [['order', 'id', 'label', 'rank', 'cluster', 'dbn', ...ms.map(m => `${m.label}${m.unit ? ' (' + m.unit + ')' : ''}`)].map(q).join(',')];
    this.order.forEach((k, i) => {
      const c = ens.conformations[k];
      rows.push([i + 1, q(c.id), q(c.label), c.rank ?? '', c.cluster ?? '', q(c.dbn), ...ms.map(m => (typeof c.values[m.key] === 'number' ? c.values[m.key] : ''))].join(','));
    });
    const saved = await saveBlob(new Blob([rows.join('\n') + '\n'], { type: 'text/csv' }), `${slug(ens.name)}_energies.csv`, [{ description: 'CSV table', accept: { 'text/csv': ['.csv'] } }]);
    if (saved) toast(saved, 'ok', 'Table saved');
  }

  // ---------------------------------------------------------------- view state, .rcv files and workspaces
  // Features that keep per-ensemble state (named selections, measurements, scenes, colour overrides, notes)
  // register a session part: { save(ens) -> json, load(ens, json) }.
  registerSessionPart(id, part) { (this._sessionParts = this._sessionParts || new Map()).set(id, part); }

  viewState() {
    const parts = {};
    for (const [id, p] of this._sessionParts || []) { try { const v = p.save(this.ens); if (v !== undefined && v !== null) parts[id] = v; } catch (err) { console.error(err); } }
    return {
      conf: this.state.conf, prev: this.state.prev, metric: this.state.metric, orderMode: this.state.orderMode, colorScheme: this.state.colorScheme,
      uniformColor: this.state.uniformColor, selection: [...this.state.selection], ghosts: this.state.ghosts, putty: this.state.putty,
      camera: this.viewer.controls.state, settings: { ...this.viewer.settings }, projection: this.viewer.projection, labelMode: this.viewer.labelMode,
      superposedOn: this.ens?.superposedOn || null, parts,
    };
  }

  applyViewState(v) {
    if (!v || !this.ens) return;
    if (v.settings) {
      const s = { ...v.settings };
      Object.assign(this.viewer.settings, s);
      if (s.style) this.viewer.setStyle(s.style);
      if (s.representation) this.viewer.setRepresentation(s.representation);
    }
    if (v.orderMode) this.state.orderMode = v.orderMode;
    if (v.colorScheme) this.state.colorScheme = v.colorScheme;
    if (v.uniformColor) this.state.uniformColor = v.uniformColor;
    // feature state first: measurements are metrics the view may refer to
    for (const [id, p] of this._sessionParts || []) if (v.parts && id in v.parts) { try { p.load(this.ens, v.parts[id]); } catch (err) { console.error(err); } }
    if (v.metric && this.ens.metricByKey[v.metric]) this.setMetric(v.metric);
    if (v.projection) this.viewer.setProjection(v.projection);
    if (v.labelMode) this.viewer.labelMode = v.labelMode;
    if (Array.isArray(v.superposedOn) && v.superposedOn.length >= 3) {
      try { this.ens.superposeOn(v.superposedOn, Number.isInteger(v.conf) ? v.conf : 0); this.ens.superposedOn = v.superposedOn; } catch { /* ignore */ }
    }
    this._computeOrder();
    this._colorCache.clear();
    if (Number.isInteger(v.conf) && v.conf < this.ens.n) this.goTo(v.conf, { animate: false, force: true });
    if (Number.isInteger(v.prev) && v.prev < this.ens.n) { this.state.prev = v.prev; this.emit('conformation', { from: v.prev, to: this.state.conf }); }
    this.state.ghosts = !!v.ghosts; this._applyGhosts();
    this.state.putty = !!v.putty; this._applyPutty();
    if (v.selection) this.select(v.selection.filter(r => r < this.ens.nRes), 'set');
    if (v.camera) { this.viewer.controls.state = v.camera; }
    this.refreshColors();
    this._applyBackground();
    this._syncToolbar();
    this.viewer.needsRender = true;
    this.emit('order');
    this.emit('settings');
  }

  // The current ensemble with this view as .rcv bytes (Blob) plus size statistics.
  async packRcv(ens = this.ens) {
    let thumbnail = null;
    try { thumbnail = (await composeFigure(this, { width: 400, height: 250, transparent: false, background: null, legend: false, labels: false })).toDataURL('image/jpeg', 0.8); } catch { /* optional */ }
    const session = { format: 'rcv.view/1', app: 'RCV', version: this.info?.version, view: this.viewState() };
    const body = ens.serverId ? { document: ens.serverId, session, thumbnail } : { doc: ens.doc, session, thumbnail };
    const r = await this.api.raw('POST', '/api/rcv/pack', JSON.stringify(body), { 'Content-Type': 'application/json' });
    let stats = {};
    try { stats = JSON.parse(r.headers.get('X-RCV-Stats') || '{}'); } catch { /* optional */ }
    return { blob: new Blob([await r.arrayBuffer()], { type: 'application/vnd.rcv+zip' }), stats };
  }

  // Save the current ensemble as one .rcv file: all conformations, energies, 2D structures, provenance and this view.
  async saveRcv() {
    const ens = this.ens;
    if (!ens) return;
    if (!this.api.online) { toast('Saving .rcv files needs the local RCV engine; use File ▸ Save workspace instead.', 'error', 'Offline'); return; }
    const done = busy('Packing the .rcv file…');
    let blob, stats = {};
    try { ({ blob, stats } = await this.packRcv()); } catch (err) { toast(err.message, 'error', 'Could not save'); return; } finally { done(); }
    const saved = await saveBlob(blob, `${slug(ens.name)}.rcv`, [{ description: 'RCV ensemble', accept: { 'application/vnd.rcv+zip': ['.rcv'] } }]);
    if (saved) {
      const ratio = stats.pdb_bytes ? stats.pdb_bytes / blob.size : null;
      toast(`${saved} · ${(blob.size / 1024).toFixed(0)} kB${ratio ? ` · ${ratio.toFixed(0)}× smaller than the ${ens.n} models as PDB` : ''} · conformations, energies, 2D structures, provenance and this view in one file`, 'ok', 'Saved', 7000);
      this.emit('saved', { name: saved, bytes: blob.size });
    }
  }

  // Workspace: every open ensemble in one JSON file (larger than .rcv; opens without the engine).
  async saveWorkspace() {
    if (!this.ens) return;
    const session = {
      format: 'rcv.session/1', created: new Date().toISOString(), app: 'RCV', current: this.ens.id,
      view: this.viewState(),
      ensembles: [...this.ensembles.values()].map(e => e.doc),
    };
    const blob = new Blob([JSON.stringify(session)], { type: 'application/json' });
    const saved = await saveBlob(blob, `${slug(this.ens.name)}.workspace.json`, [{ description: 'RCV workspace', accept: { 'application/json': ['.json'] } }]);
    if (saved) toast(`${saved} · ${(blob.size / 1e6).toFixed(1)} MB`, 'ok', 'Workspace saved');
  }

  openSessionFile() {
    const inp = el('input', { type: 'file', accept: '.rcv,.json' });
    inp.addEventListener('change', async () => { if (inp.files[0]) await this.openLocalFiles([inp.files[0]]); });
    inp.click();
  }

  async loadSessionText(text, name) {
    let obj;
    try { obj = JSON.parse(text); } catch { toast(`${name} is not valid JSON`, 'error'); return; }
    if (DOC_SCHEMAS.includes(obj.schema)) { this.addDocument(obj); return; }
    if (!['rcv.session/1', 'rnaconfswitch.session/1', 'rnafold3d.session/1'].includes(obj.format)) { toast(`${name} is not an RCV workspace`, 'error'); return; }
    for (const d of obj.ensembles || []) this.addDocument(d, { select: false });
    const cur = obj.current && this.ensembles.has(obj.current) ? obj.current : [...this.ensembles.keys()][0];
    if (!cur) return;
    this.selectEnsemble(cur);
    this.showPage('viewer');
    this.applyViewState(obj.view || {});
    toast(name, 'ok', 'Workspace opened');
  }

  // ---------------------------------------------------------------- input
  _bindKeys() {
    window.addEventListener('keydown', e => {
      if (e.target.closest('input, textarea, select, dialog')) return;
      const k = e.key;
      const ctrl = e.ctrlKey || e.metaKey;
      if (ctrl && k.toLowerCase() === 'o') { e.preventDefault(); openDialog(this); return; }
      if (ctrl && k.toLowerCase() === 'n') { e.preventDefault(); this.action('predict'); return; }
      if (ctrl && k.toLowerCase() === 'h') { e.preventDefault(); this.showPage('welcome'); return; }
      if (ctrl && e.shiftKey && k.toLowerCase() === 'f') { e.preventDefault(); this.fetchDialog(); return; }
      if (k === 'F11') { e.preventDefault(); this.wm.toggleFullscreen(); return; }
      if (k === 'F1') { e.preventDefault(); userGuideDialog(this); return; }
      if (this.page !== 'viewer') return;
      if (ctrl && k.toLowerCase() === 's') { e.preventDefault(); if (e.shiftKey) this.saveWorkspace(); else this.saveRcv(); return; }
      if (ctrl && k.toLowerCase() === 'e') { e.preventDefault(); if (this.ens) exportImageDialog(this); return; }
      if (ctrl && (k.toLowerCase() === 'y' || (k.toLowerCase() === 'z' && e.shiftKey))) { e.preventDefault(); this.history.redo(); return; }
      if (ctrl && k.toLowerCase() === 'z') { e.preventDefault(); this.history.undo(); return; }
      if (ctrl && k.toLowerCase() === 'a') { e.preventDefault(); if (this.ens) this.select([...Array(this.ens.nRes).keys()], 'set'); return; }
      if (ctrl || e.altKey) return;
      const reps = ['cartoon', 'cartoon+sticks', 'sticks', 'ballstick', 'spacefill'];
      switch (k) {
        case 'ArrowRight': this.next(); break;
        case 'ArrowLeft': this.prev(); break;
        case 'Home': this.first(); break;
        case 'End': this.last(); break;
        case ' ': this.togglePlay(); break;
        case 'Escape': if (this.wm.escape()) break; if (this.tool !== 'select') { this.tools.set('select'); break; } this.stop(); this.select([], 'set'); break;
        case 'v': case 'V': this.tools.set('select'); break;
        case 'b': case 'B': this.tools.set(this.tool === 'paint' ? 'select' : 'paint'); break;
        case 't': case 'T': this.tools.set(this.tool === 'label' ? 'select' : 'label'); break;
        case 'r': case 'R': this.viewer.resetView(); break;
        case 'a': case 'A': this.action('spin'); break;
        case 'e': case 'E': this.toggleGhosts(); break;
        case 'g': case 'G': this.togglePull(); break;
        case 'p': case 'P': this.togglePutty(); break;
        case 'd': case 'D': this.action('ss2d'); break;
        case 'u': case 'U': if (this.ens) this.surface.set(!this.surface.enabled); break;
        case 'm': case 'M': this.action('measure'); break;
        case '`': this.action('console'); break;
        case 'f': case 'F': this.centerOnSelection(); break;
        case 'h': case 'H': this.wm.setFocus(!this.wm.focus); break;
        case 's': case 'S': this.quickSnapshot(); break;
        case 'c': case 'C': {
          const i = SCHEMES.findIndex(s => s.id === this.state.colorScheme);
          this.setColorScheme(SCHEMES[(i + 1) % SCHEMES.length].id);
          break;
        }
        case 'l': case 'L': {
          const modes = ['selection', 'every10', 'all', 'none'];
          this.setLabelMode(modes[(modes.indexOf(this.viewer.labelMode) + 1) % modes.length]);
          toast(`Labels: ${{ selection: 'selected residues', every10: 'every 10th residue', all: 'all residues', none: 'none' }[this.viewer.labelMode]}`, 'info', '', 1500);
          break;
        }
        case '?': shortcutsDialog(); break;
        default:
          if (k >= '1' && k <= '5') this.setRepresentation(reps[+k - 1]);
          else return;
      }
      e.preventDefault();
    });
  }

  _bindDrop() {
    const dz = $('#dropzone');
    let depth = 0;
    window.addEventListener('dragenter', e => { if ([...(e.dataTransfer?.types || [])].includes('Files')) { depth++; dz.hidden = false; } });
    window.addEventListener('dragleave', () => { depth = Math.max(0, depth - 1); if (!depth) dz.hidden = true; });
    window.addEventListener('dragover', e => e.preventDefault());
    window.addEventListener('drop', async e => {
      e.preventDefault();
      depth = 0; dz.hidden = true;
      if (e.target.closest && e.target.closest('dialog')) return;
      await this.openDropped(e.dataTransfer, {});
    });
  }
}

const app = new App();
window.rcv = app;                  // for debugging, the command console and automated tests
app.init().catch(err => { console.error(err); toast(err.message || String(err), 'error', 'Start-up failed'); });
