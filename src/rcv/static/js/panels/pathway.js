// Pathway tab: how conformation A turns into conformation B.
//   * secondary-structure refolding paths (ViennaRNA findpath) - direct and through intermediates - with barriers,
//     flux-ranked routes (transition path theory) and an order-of-magnitude rate;
//   * the bottleneck: the cooperative unit of base pairs that must open together before the highest saddle, the
//     nucleus of B that forms after it, and the rigid module whose motion carries them;
//   * rigid modules and hinges between the two 3D models, and a mechanism morph that plays the motion in order;
//   * an in-silico phi-value scan of every point mutation (which mutations speed the switch up or slow it down).
// Computed by the local engine (rcv/transition.py); energies under the conditions of the Conditions dialog.
import { el, svg, fmt, fmtSigned } from '../util.js';
import { icon } from '../icons.js';
import { tooltip, tipRow, toast, select, popupMenu } from '../ui.js';
import { linear, extent, yAxis, xAxis } from './charts.js';
import { chartMenuItems } from '../export/chartexport.js';
import { registerScheme, sampleRamp, ORANGE_RAMP, VIRIDIS, SLOTS } from '../colors.js';
import { divergingColor } from './switch.js';
import { drawMap, landscapeOf, loadLandscape } from '../landscape.js';
import { conditionsDialog, conditionsParts, conditionsOf, requestConditions } from '../conditions.js';
import { buildMechanism } from '../mechmorph.js';

const K0S = [1e4, 1e5, 1e6, 1e7, 1e8];
const ROUTE_COLORS = ['var(--series-2)', 'var(--series-1)', '#1baf7a', '#6250d6'];

export function fmtTime(s) {
  if (s === null || s === undefined || !Number.isFinite(s)) return '–';
  if (s < 1e-6) return `${fmt(s * 1e9, 1)} ns`;
  if (s < 1e-3) return `${fmt(s * 1e6, 1)} µs`;
  if (s < 1) return `${fmt(s * 1e3, 1)} ms`;
  if (s < 120) return `${fmt(s, 1)} s`;
  if (s < 7200) return `${fmt(s / 60, 1)} min`;
  if (s < 172800) return `${fmt(s / 3600, 1)} h`;
  if (s < 3.15e7) return `${fmt(s / 86400, 1)} days`;
  return `${(s / 3.15e7).toExponential(1)} years`;
}

export function verdict(halfLife) {
  if (halfLife === null || halfLife === undefined) return { kind: 'same', text: 'A and B have the same secondary structure: the change is a tertiary rearrangement without base-pair changes, so there is no 2D barrier.' };
  if (halfLife < 1e-3) return { kind: 'good', text: 'Fast: well below a millisecond.' };
  if (halfLife < 1) return { kind: 'good', text: 'Fast: within milliseconds.' };
  if (halfLife < 60) return { kind: 'good', text: 'Feasible within seconds.' };
  if (halfLife < 3600) return { kind: 'good', text: 'Feasible within minutes.' };
  if (halfLife < 86400) return { kind: 'warn', text: 'Slow: hours.' };
  if (halfLife < 3.15e7) return { kind: 'warn', text: 'Very slow: days to months.' };
  return { kind: 'bad', text: 'Not feasible by refolding alone under these conditions (half-life over a year): the switch needs unfolding, a helper such as a helicase or chaperone, or happens during transcription.' };
}

const sup = x => String(x).replace(/[0-9-]/g, c => '⁰¹²³⁴⁵⁶⁷⁸⁹'['0123456789'.indexOf(c)] || '⁻');
const k0Label = k => `10${sup(Math.round(Math.log10(k)))} s⁻¹`;

export class PathwayPanel {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    this.a = null; this.b = null;
    this.k0 = 1e6;
    this.result = null; this.scan = null;
    this.loading = false; this.scanning = false; this.error = null;
    app.on('ensemble', () => { this.result = this.scan = null; this.a = this.b = null; this.error = null; if (this.visible()) this.load(); else this.render(); });
    app.on('tab', id => { if (id === 'pathway' && !this.result && !this.loading && app.ens) this.load(); });
    app.on('conditions', () => { this.result = this.scan = null; if (this.visible()) this.load(); });
    app.on('conformation', () => this._markCurrent());
    new ResizeObserver(() => { const w = root.clientWidth; if (this.result && Math.abs(w - (this._w || 0)) > 8) { this._w = w; this.render(); } }).observe(root);
    this.render();
  }

  visible() { return !this.root.hidden; }
  width() { return Math.max(240, this.root.clientWidth - 28); }

  _defaultPair() {
    const app = this.app, ens = app.ens, cur = app.state.conf, prev = app.state.prev;
    const key = k => ens.conformations[k].layout_key || ens.conformations[k].dbn;
    if (prev !== null && prev !== cur && key(prev) !== key(cur)) return [prev, cur];
    let best = null;
    for (let k = 0; k < ens.n; k++) {
      if (k === cur) continue;
      const other = key(k) !== key(cur);
      const score = (other ? 1000 : 0) + ens.rmsd[cur][k];
      if (!best || score > best[0]) best = [score, k];
    }
    return [cur, best ? best[1] : cur];
  }

  async load() {
    const app = this.app, ens = app.ens;
    if (!ens) { this.render(); return; }
    if (ens.n < 2) { this.error = 'A pathway needs two conformations; this ensemble has one.'; this.render(); return; }
    if (!app.api.online) { this.error = 'The pathway analysis needs the local RCV engine.'; this.render(); return; }
    if (!app.info?.capabilities?.vienna) { this.error = 'The pathway analysis needs ViennaRNA, which is not installed.'; this.render(); return; }
    if (this.a === null || this.b === null) [this.a, this.b] = this._defaultPair();
    this.error = null;
    this.loading = true;
    this.render();
    loadLandscape(app, ens);
    try {
      const body = { a: this.a, b: this.b, k0: this.k0, conditions: requestConditions(app, ens, { probing: false }) };
      if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
      const r = await app.api.post('/api/transition', body);
      if (app.ens !== ens) return;
      this.result = r;
      this.scan = null;
      ens._pathway = r;                     // for figures, videos and the HTML page
    } catch (err) { this.error = err.message; } finally { this.loading = false; }
    this.render();
    app.emit('pathway', this.result);
  }

  async runScan() {
    const app = this.app, ens = app.ens, r = this.result;
    if (!r || !ens) return;
    this.scanning = true;
    this.render();
    try {
      const route = r.best.route.map(k => r.nodes[k].dbn);
      const body = { route, conditions: requestConditions(app, ens, { probing: false }) };
      if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
      const sc = await app.api.post('/api/transition/scan', body);
      if (app.ens !== ens || this.result !== r) return;
      this.scan = sc;
      ens._phiScan = sc;
    } catch (err) { toast(err.message, 'error', 'Kinetic scan'); } finally { this.scanning = false; }
    this.render();
  }

  // representative 3D conformation of a network node (A, B, or an ensemble state), or null for a 2D-only minimum
  nodeConf(k) {
    const r = this.result;
    if (k === 0) return r.a;
    if (k === 1) return r.b;
    const n = r.nodes[k];
    return n?.conformations?.length ? n.conformations[0] : null;
  }

  nodeName(k) {
    const r = this.result, n = r.nodes[k];
    if (k === 0) return 'A';
    if (k === 1) return 'B';
    if (n.kind === 'ensemble') return n.label || `S${k}`;
    const idx = r.nodes.slice(0, k + 1).filter(x => x.kind === 'minimum').length;
    return `I${idx}`;
  }

  render() {
    const app = this.app, ens = app.ens, root = this.root;
    root.replaceChildren();
    if (!ens) { root.append(el('div', { class: 'muted', text: 'No ensemble loaded.' })); return; }
    const head = el('div', { class: 'section' });
    root.append(head);
    head.append(el('div', { class: 'section-title' }, 'Transition A → B', el('span', { class: 'grow' }),
      el('button', { class: 'linkbtn', title: 'Temperature, Na⁺ and ligand used for the energies', onclick: () => conditionsDialog(app) }, ...conditionsParts(conditionsOf(app, ens), ens))));
    if (ens.n >= 2) {
      const opts = [...Array(ens.n).keys()].map(k => ({ value: String(k), label: ens.shortId(k) }));
      const selA = select(opts, String(this.a ?? app.state.conf), { 'aria-label': 'Conformation A' });
      const selB = select(opts, String(this.b ?? app.state.conf), { 'aria-label': 'Conformation B' });
      const go = () => { this.a = +selA.value; this.b = +selB.value; if (this.a === this.b) { toast('Choose two different conformations.', 'info'); return; } this.result = this.scan = null; this.load(); };
      selA.addEventListener('change', go); selB.addEventListener('change', go);
      head.append(el('div', { class: 'sw-ab pw-ab' },
        el('span', { class: 'sw-key a', text: 'A' }), selA,
        el('button', { class: 'btn small ghost icon-only', title: 'Swap A and B', 'aria-label': 'Swap A and B', onclick: () => { [this.a, this.b] = [this.b, this.a]; this.result = this.scan = null; this.load(); } }, icon('swap')),
        el('span', { class: 'sw-key b', text: 'B' }), selB),
        el('div', { class: 'pw-prev' }, el('button', { class: 'linkbtn', title: 'Use the previous and the shown conformation', onclick: () => { const p = app.state.prev, c = app.state.conf; if (p === null || p === c) { toast('Show another conformation first: the previous one becomes A.', 'info'); return; } this.a = p; this.b = c; this.result = this.scan = null; this.load(); } }, 'previous → shown')));
    }
    if (this.loading) { head.append(el('div', { class: 'row' }, el('span', { class: 'spinner small' }), ' computing refolding paths…')); return; }
    if (this.error) { head.append(el('div', { class: 'warnbox', text: this.error })); return; }
    const r = this.result;
    if (!r) { head.append(el('button', { class: 'btn', onclick: () => this.load() }, 'Analyse')); return; }
    this._verdict(head, r);
    this._mapSection(r);
    if (r.A.dbn !== r.B.dbn) {
      this._profileSection(r);
      this._pathwaysSection(r);
      this._bottleneckSection(r);
    }
    this._motionSection(r);
    if (r.A.dbn !== r.B.dbn && !r.same_basin) this._scanSection(r);
    root.append(el('div', { class: 'section' },
      el('div', { class: 'desc', text: `Secondary-structure level: ViennaRNA ${r.method?.vienna || ''} Turner energies; refolding paths by findpath (Flamm et al. 2001), search width ${r.method?.width ?? '–'}; routes through the ensemble's 2D states and local minima from ${r.method?.samples ?? 0} Boltzmann samples; flux by transition path theory with rates k = k₀·exp(−(saddle − E)/RT).` }),
      el('div', { class: 'desc', text: 'Barriers are those of the paths found (upper bounds); rates scale with k₀ and are order-of-magnitude estimates. 3D modules come from the two models and are only as reliable as the models.' })));
  }

  _verdict(head, r) {
    const app = this.app, ens = app.ens;
    const same = r.A.dbn === r.B.dbn || r.same_basin;
    const v = r.same_basin && r.A.dbn !== r.B.dbn
      ? { kind: 'same', text: 'A and B relax into the same basin of the 2D landscape (the same local minimum): switching between them needs no barrier crossing.' }
      : verdict(same ? null : r.half_life_s);
    const card = el('div', { class: `pw-verdict ${v.kind}` });
    if (!same) {
      const faster = r.rate && r.rate_direct ? r.rate / r.rate_direct : null;
      card.append(
        el('div', { class: 'pw-big' }, el('span', { class: 'pw-num', text: `${fmt(r.effective_barrier, 1)}` }), el('span', { class: 'pw-unit', text: ' kcal/mol' }),
          el('span', { class: 'pw-sub', text: ' effective barrier A → B' })),
        el('div', { class: 'pw-line', text: `Direct path ${fmt(r.direct.barrier, 1)} · best single route ${fmt(r.best.route_barrier ?? r.best.barrier, 1)} kcal/mol · ΔG(A→B) ${fmtSigned(r.dG, 1)} kcal/mol · ${r.bp_distance} base-pair changes` }),
        el('div', { class: 'pw-line' }, `Half-life of A ≈ `, el('b', { text: fmtTime(r.half_life_s) }), ` with k₀ = `,
          (() => { const s = select(K0S.map(k => ({ value: String(k), label: k0Label(k) })), String(this.k0), { class: 'select small', 'aria-label': 'Prefactor k0' });
            s.addEventListener('change', () => { this.k0 = +s.value; this.result = null; this.load(); }); return s; })(),
          faster && faster > 1.5 ? ` · ${fmt(faster, faster < 10 ? 1 : 0)}× faster than the direct path through ${r.nodes.length - 2} intermediate${r.nodes.length - 2 === 1 ? '' : 's'}` : ''),
        el('div', { class: 'pw-verdict-text', text: v.text }));
      if ((r.A.above_basin ?? 0) > 0.05 || (r.B.above_basin ?? 0) > 0.05) card.append(el('div', { class: 'pw-line', text: `A lies ${fmt(r.A.above_basin, 1)} kcal/mol above the minimum of its basin and B ${fmt(r.B.above_basin, 1)} kcal/mol above its own; barriers and the time to switch are measured between the basins.` }));
      if (r.trap) card.append(el('div', { class: 'pw-line', text: r.trap.kind === 'deep'
        ? `Trap: on the way the RNA can fall into ${this.nodeName(r.trap.node)}, ${fmt(r.trap.depth, 1)} kcal/mol below A and deeper than B, and must climb out of it to reach B; the effective barrier counts that wait.`
        : r.trap.kind === 'intermediate'
          ? `Slow step: on the way the RNA falls into ${this.nodeName(r.trap.node)}, ${fmt(r.trap.depth, 1)} kcal/mol below A, and spends most of the switching time there before it reaches B; the effective barrier counts that wait.`
          : `Trap: ${this.nodeName(r.trap.node)} lies ${fmt(r.trap.depth, 1)} kcal/mol below A on the A side (committor < 0.5), so the RNA mostly waits there before it reaches B; the effective barrier counts that wait.` }));
    } else {
      card.append(el('div', { class: 'pw-verdict-text', text: v.text }));
      if (r.same_basin && r.A.dbn !== r.B.dbn) card.append(el('div', { class: 'pw-line', text: `A lies ${fmt(r.A.above_basin, 1)} and B ${fmt(r.B.above_basin, 1)} kcal/mol above the minimum of their shared basin (ΔG A → B ${fmtSigned(r.dG, 1)} kcal/mol): both are fluctuations of one fold rather than separate states.` }));
      card.append(el('div', { class: 'pw-line', text: `C4′ RMSD between A and B ${fmt(r.rmsd_ab, 1)} Å` }));
    }
    head.append(card);
    void ens;
  }

  _mapSection(r) {
    const app = this.app, ens = app.ens;
    const L = landscapeOf(ens);
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Routes on the ensemble map', el('span', { class: 'grow' }), this._exportBtn(() => this.mapSvg, 'Pathway map', 'landscape')));
    const W = this.width(), H = Math.min(300, Math.max(200, W * 0.62));
    const s = svg('svg', { class: 'chart pw-map', viewBox: `0 0 ${W} ${H}`, width: '100%', height: H, role: 'img', 'aria-label': 'Ensemble map with the transition routes' });
    this.mapSvg = s;
    sec.append(s);
    if (!L) {
      sec.append(el('div', { class: 'desc', text: 'Computing the ensemble map…' }));
      loadLandscape(app, ens).then(x => { if (x && this.result === r) this.render(); });
      this.root.append(sec);
      return;
    }
    const toConfs = nodes => nodes.map(k => this.nodeConf(k)).filter(k => k !== null);
    // routes through the same 3D conformations are drawn once (their 2D-only intermediates differ)
    const merged = new Map();
    for (const p of (r.same_basin ? [] : r.tpt?.pathways?.length ? r.tpt.pathways : [{ nodes: r.best.route, flux_share: 1 }])) {
      const confs = toConfs(p.nodes);
      const key = confs.join('>');
      const m = merged.get(key) || { confs, share: 0, via: [] };
      m.share += p.flux_share;
      m.via.push(p.nodes.slice(1, -1).map(k => this.nodeName(k)).join('·') || 'A → B');
      merged.set(key, m);
    }
    const drawn = [...merged.values()].sort((x, y) => y.share - x.share).slice(0, 3);
    const routes = drawn.map((m, i) => ({ nodes: m.confs, width: 1.2 + 4 * Math.min(1, m.share), color: ROUTE_COLORS[i % ROUTE_COLORS.length], opacity: i ? 0.75 : 1 }));
    const paths = drawn;
    const inRoute = new Set(routes.flatMap(x => x.nodes));
    const g = svg('g');
    const res = drawMap(g, ens, L, { x0: 6, y0: 4, x1: W - 6, y1: H - 6 }, {
      routes, current: app.state.conf, color: k => (k === r.a ? 'var(--series-1)' : k === r.b ? 'var(--series-2)' : inRoute.has(k) ? 'var(--text-2)' : 'var(--axis)'),
      radius: k => (k === r.a || k === r.b ? 5.5 : 3.5), marks: [{ k: r.a, text: 'A' }, { k: r.b, text: 'B' }], fontSize: 9.5 });
    s.append(g);
    this._mapPos = res;
    const hit = svg('rect', { class: 'hit', x: 0, y: 0, width: W, height: H });
    const nearest = e => {
      const b = s.getBoundingClientRect();
      const x = (e.clientX - b.left) * W / b.width, y = (e.clientY - b.top) * H / b.height;
      let bk = -1, bd = 20;
      for (let k = 0; k < ens.n; k++) { const d = Math.hypot(res.px(k) - x, res.py(k) - y); if (d < bd) { bd = d; bk = k; } }
      return bk;
    };
    hit.addEventListener('pointermove', e => {
      const k = nearest(e);
      if (k < 0) { tooltip.hide(); app.previewConformation(null); return; }
      app.previewConformation(k);
      tooltip.show([tipRow(ens.shortId(k), k === r.a ? 'A' : k === r.b ? 'B' : 'conformation'), el('div', { class: 'tl', text: 'click: show · Shift+click: make it B' })], e.clientX, e.clientY);
    });
    hit.addEventListener('pointerleave', () => { tooltip.hide(); app.previewConformation(null); });
    hit.addEventListener('click', e => {
      const k = nearest(e);
      if (k < 0) return;
      if (e.shiftKey && k !== this.a) { this.b = k; this.result = this.scan = null; this.load(); } else app.goTo(k);
    });
    s.append(hit);
    const legend = el('div', { class: 'pw-legend' });
    paths.forEach((m, i) => legend.append(el('span', {}, el('i', { style: { background: ROUTE_COLORS[i % ROUTE_COLORS.length] } }),
      `${m.confs.map(k => ens.shortId(k)).join(' → ')}: ${Math.round(m.share * 100)} % of the A → B flux${m.via.some(v => v !== 'A → B') ? ` (${m.via.length} route${m.via.length === 1 ? '' : 's'}: ${[...new Set(m.via)].slice(0, 3).join(', ')}${m.via.length > 3 ? ', …' : ''})` : ''}`)));
    sec.append(legend, el('div', { class: 'desc', text: 'Disks are secondary-structure topologies; arrows run through the conformations of the routes (intermediates without a 3D model are not drawn). Shift+click a point to make it B.' }));
    this.root.append(sec);
  }

  _profileSection(r) {
    const app = this.app, ens = app.ens;
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Energy along the refolding path', el('span', { class: 'grow' }), this._exportBtn(() => this.profSvg, 'Refolding path', 'pathway')));
    const W = this.width(), H = 190, ml = 50, mr = 10, mt = 14, mb = 34;
    const s = svg('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}`, width: '100%', height: H, role: 'img', 'aria-label': 'Free energy along the refolding path' });
    this.profSvg = s;
    const best = r.best.profile, direct = r.direct.profile;
    const nx = Math.max(best.length, direct.length) - 1 || 1;
    const x = linear(0, nx, ml, W - mr);
    const [lo, hi] = extent([...best.map(p => p.E), ...direct.map(p => p.E)], 0.12);
    const y = linear(lo, hi, H - mb, mt);
    const g = svg('g');
    yAxis(g, y, ml, W - mr, { ticks: 4, title: 'ΔG (kcal/mol)', grid: true });
    const ticks = [0, Math.round(nx / 2), nx].filter((v, i, a) => a.indexOf(v) === i);
    xAxis(g, ml, W - mr, H - mb, ticks.map(v => [x(v), String(v)]), { title: 'Base-pair moves along the path' });
    const line = (pts, i0, attrs) => svg('path', { d: pts.map((p, i) => `${i ? 'L' : 'M'} ${x(i0 + i).toFixed(1)} ${y(p.E).toFixed(1)}`).join(' '), fill: 'none', ...attrs });
    const twoRoutes = r.best.route.length > 2 || best.some(p => p.phase && p.phase !== 'route');
    if (twoRoutes) g.append(line(direct, 0, { stroke: 'var(--series-1)', 'stroke-width': 1.4, 'stroke-dasharray': '4 3' }));
    const ra = best.filter(p => p.phase === 'relax_a').length;
    const rb = best.filter(p => p.phase === 'relax_b').length;
    const core = best.slice(ra, best.length - rb);
    if (ra) g.append(line(best.slice(0, ra + 1), 0, { stroke: 'var(--text-3)', 'stroke-width': 1.6, 'stroke-dasharray': '2 3' }));
    g.append(line(core, ra, { stroke: 'var(--series-2)', 'stroke-width': 2 }));
    if (rb) g.append(line(best.slice(best.length - rb - 1), best.length - rb - 1, { stroke: 'var(--text-3)', 'stroke-width': 1.6, 'stroke-dasharray': '2 3' }));
    const mark = (i, label, hollow = true) => {
      g.append(svg('circle', { cx: x(i), cy: y(best[i].E), r: 3.5, fill: hollow ? 'var(--surface)' : 'var(--text-2)', stroke: 'var(--text-2)', 'stroke-width': 1.5 }));
      g.append(svg('text', { x: x(i), y: y(best[i].E) + 13, 'text-anchor': 'middle', 'font-size': 9.5, fill: 'var(--text-2)', text: label }));
    };
    // route nodes: the start of every leg (the basin of A first), then the basin of B; A and B themselves are filled
    // when they are not minima and relax into their basins first
    if (r.same_basin) mark(ra, 'basin');
    else {
      core.forEach((p, i) => { if (i === 0 || p.leg !== core[i - 1].leg) mark(ra + (i === 0 ? 0 : i - 1), i === 0 && ra ? 'basin A' : this.nodeName(r.best.route[p.leg])); });
      mark(ra + core.length - 1, rb ? 'basin B' : 'B');
    }
    if (ra) mark(0, 'A', false);
    if (rb) mark(best.length - 1, 'B', false);
    const sd = r.bottleneck?.saddle_step;
    if (sd !== undefined && sd !== null && best[sd]) {
      const sx = x(sd), sy = y(best[sd].E);
      g.append(svg('path', { d: `M ${sx} ${sy - 4} L ${sx + 5} ${sy - 12} L ${sx - 5} ${sy - 12} Z`, fill: 'var(--text)' }));
      g.append(svg('text', { x: sx, y: sy - 15, 'text-anchor': 'middle', 'font-size': 9.5, fill: 'var(--text)', text: `saddle ${fmt(best[sd].E - core[0].E, 1)}` }));
    }
    // hover: structure and events at each step
    const evAt = new Map();
    for (const e of r.events) (evAt.get(e.step) || evAt.set(e.step, []).get(e.step)).push(e);
    const band = (W - ml - mr) / nx;
    best.forEach((p, i) => {
      const hit = svg('rect', { class: 'hit', x: x(i) - band / 2, y: mt, width: band, height: H - mt - mb });
      hit.addEventListener('pointermove', e => {
        const evs = evAt.get(i) || [];
        const res = evs.flatMap(v => [v.i, v.j]);
        app.hoverResidues(res);
        tooltip.show([tipRow(`${fmt(p.E, 1)} kcal/mol`, `step ${i}`), el('div', { class: 'tl mono', text: p.s }),
          ...evs.map(v => el('div', { class: 'tl', text: `${v.type === 'open' ? 'opens' : 'forms'} ${ens.resLabel(v.i)}–${ens.resLabel(v.j)}` }))], e.clientX, e.clientY);
      });
      hit.addEventListener('pointerleave', () => { tooltip.hide(); app.hoverResidues([]); });
      g.append(hit);
    });
    s.append(g);
    sec.append(s, el('div', { class: 'pw-legend' },
      r.same_basin ? null : el('span', {}, el('i', { style: { background: 'var(--series-2)' } }), r.best.route.length > 2 ? 'best route through intermediates' : 'refolding path'),
      twoRoutes ? el('span', {}, el('i', { class: 'dash', style: { background: 'var(--series-1)' } }), 'direct path') : null,
      ra || rb ? el('span', {}, el('i', { class: 'dash', style: { background: 'var(--text-3)' } }), 'relaxation into the basins') : null),
      el('div', { class: 'desc', text: 'Hover a step to see the structure and the base pairs that open or form (highlighted in 3D).' }));
    this.root.append(sec);
  }

  _pathwaysSection(r) {
    const paths = r.tpt?.pathways || [];
    if (!paths.length) return;
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Routes ranked by reactive flux'));
    const sad = new Map(r.edges.map(e => [`${Math.min(e.a, e.b)}-${Math.max(e.a, e.b)}`, e.saddle]));
    const tb = el('tbody');
    paths.forEach((p, i) => {
      let top = -Infinity;
      for (let j = 0; j + 1 < p.nodes.length; j++) top = Math.max(top, sad.get(`${Math.min(p.nodes[j], p.nodes[j + 1])}-${Math.max(p.nodes[j], p.nodes[j + 1])}`) ?? -Infinity);
      tb.append(el('tr', {}, el('td', {}, el('span', { class: 'pw-swatch', style: { background: ROUTE_COLORS[i % ROUTE_COLORS.length] } }), String(i + 1)),
        el('td', { class: 'num', text: `${Math.round(p.flux_share * 100)} %` }),
        el('td', { text: p.nodes.map(k => this.nodeName(k)).join(' → ') }),
        el('td', { class: 'num', text: Number.isFinite(top) ? fmt(top - (r.A.basin_energy ?? r.A.energy), 1) : '–' })));
    });
    const covered = paths.reduce((a, p) => a + p.flux_share, 0);
    sec.append(el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: '#' }), el('th', { class: 'num', text: 'Flux' }), el('th', { text: 'Route' }), el('th', { class: 'num', text: 'Barrier' }))), tb),
      el('div', { class: 'desc', text: `${paths.length === 1 ? '1 route carries' : `${paths.length} routes carry`} ${Math.round(covered * 100)} % of the A → B flux.` }),
      el('div', { class: 'desc', text: 'I = local minima of the 2D landscape, S = 2D states of the ensemble; barrier = highest saddle above A (kcal/mol).' }));
    this.root.append(sec);
  }

  _bottleneckSection(r) {
    const app = this.app, ens = app.ens, b = r.bottleneck;
    if (!b) return;
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Bottleneck'));
    const helixText = hs => hs.map(h => {
      const i0 = h[0][0], i1 = h[h.length - 1][0], j0 = h[h.length - 1][1], j1 = h[0][1];
      return `${ens.resLabel(i0)}${h.length > 1 ? '–' + ens.resLabel(i1) : ''} / ${ens.resLabel(j0)}${h.length > 1 ? '–' + ens.resLabel(j1) : ''} (${h.length} bp)`;
    }).join(', ');
    const unitN = b.unit.length, nucN = b.nucleus.length;
    sec.append(el('div', { class: 'pw-bneck' },
      el('div', {}, el('b', { text: 'Cooperative unit: ' }), unitN ? `${unitN} base pair${unitN === 1 ? '' : 's'} must open together on the climb of ${fmt(b.climb_energy, 1)} kcal/mol to the highest saddle — ${helixText(b.unit_helices)}.` : 'no pair has to open before the saddle.'),
      el('div', {}, el('b', { text: 'Nucleus of B: ' }), nucN ? `${nucN === 1 ? '1 base pair forms' : `${nucN} base pairs form`} right after the saddle — ${helixText(b.nucleus_helices)}.` : '–'),
      r.coupling && r.motion ? (() => {
        const m = r.motion.modules[r.coupling.module];
        return el('div', {}, el('b', { text: 'Rate-limiting coordinated motion: ' }), `module M${m.id + 1} (${m.size} nt, ${rangesText(ens, m.residues)}) turns ${fmt(m.rotation, 0)}° and shifts ${fmt(m.translation, 1)} Å relative to M1; it carries ${r.coupling.unit_residues.length} nucleotide${r.coupling.unit_residues.length === 1 ? '' : 's'} of the cooperative unit.`);
      })() : null),
      el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '8px' } },
        el('button', { class: 'btn small', onclick: () => { app.select([...new Set(b.unit.flat())], 'set'); app.centerOnSelection?.(); } }, icon('target'), 'Select the unit'),
        el('button', { class: 'btn small', onclick: () => this.colorBy('role') }, icon('layers'), 'Colour by mechanism'),
        el('button', { class: 'btn small ghost', onclick: () => this.colorBy('order') }, 'Colour by event order')));
    this.root.append(sec);
  }

  _motionSection(r) {
    const app = this.app, ens = app.ens, mo = r.motion;
    if (!mo) return;
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Coordinated motions (rigid modules)'));
    if (!mo.modules.length) {
      sec.append(el('div', { class: 'desc', text: 'No group of residues stays rigid between A and B (fit within 2 Å).' }));
      this.root.append(sec);
      return;
    }
    const tb = el('tbody');
    for (const m of mo.modules) {
      const tr = el('tr', { class: 'clickable' },
        el('td', {}, el('span', { class: 'pw-swatch', style: { background: SLOTS[m.id % SLOTS.length] } }), `M${m.id + 1}`),
        el('td', { text: rangesText(ens, m.residues) }), el('td', { class: 'num', text: String(m.size) }),
        el('td', { class: 'num', text: m.id === 0 ? 'reference' : `${fmt(m.rotation, 0)}°` }),
        el('td', { class: 'num', text: m.id === 0 ? '–' : `${fmt(m.translation, 1)} Å` }),
        el('td', { class: 'num', text: m.window ? `${fmt(m.window[0], 2)}–${fmt(m.window[1], 2)}` : '–' }));
      tr.addEventListener('click', () => app.select(m.residues, 'set'));
      tr.addEventListener('pointerenter', () => app.hoverResidues(m.residues));
      tr.addEventListener('pointerleave', () => app.hoverResidues([]));
      tb.append(tr);
    }
    const hinges = mo.hinges.slice(0, 6).map(k => ens.resLabel(k)).join(', ');
    sec.append(el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: '' }), el('th', { text: 'Residues' }), el('th', { class: 'num', text: 'nt' }),
      el('th', { class: 'num', text: 'Turn' }), el('th', { class: 'num', text: 'Shift' }), el('th', { class: 'num', text: 'When', title: 'Time window of the base-pair events of the module along the route (0 = A, 1 = B)' }))), tb),
      el('div', { class: 'desc', text: mo.modules.length === 1 ? `1 module moves as a rigid body. Hinges (largest backbone change): ${hinges || '–'}.` : `${mo.modules.length} modules move as rigid bodies (turn and shift relative to M1, the largest). Hinges (largest backbone change): ${hinges || '–'}.` }),
      el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '6px' } },
        el('button', { class: 'btn small primary', onclick: () => this.play(true) }, icon('play'), 'Play the mechanism'),
        el('button', { class: 'btn small', onclick: () => this.play(false) }, 'Rigid-module morph'),
        el('button', { class: 'btn small ghost', onclick: () => this.colorBy('modules') }, icon('layers'), 'Colour by module'),
        el('button', { class: 'btn small ghost', onclick: () => app.transitionVideo?.(r.a, r.b, { mechanism: { motion: r.motion, timing: r.timing } }) }, icon('film'), 'Mechanism video…')),
      el('div', { class: 'desc', text: 'The mechanism morph turns each module about its centre in the order of its base-pair events on the route and keeps the chain connected; it illustrates the analysed mechanism and is not a simulation.' }));
    this.root.append(sec);
  }

  _scanSection(r) {
    const app = this.app, ens = app.ens;
    const sec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Mutations that change the switching rate (φ-value scan)'));
    this.root.append(sec);
    if (this.scanning) { sec.append(el('div', { class: 'row' }, el('span', { class: 'spinner small' }), ` evaluating ${3 * ens.nRes} point mutations…`)); return; }
    if (!this.scan) {
      sec.append(el('div', { class: 'desc', text: 'For every point mutation: the change of the barrier ΔΔG‡ (blue = faster switching, orange = slower) and the φ value ΔΔG‡/ΔΔG — near 1 where the nucleotide already looks as in B at the saddle (its A pairs open or its B pairs formed), near 0 where it still looks as in A.' }),
        el('button', { class: 'btn small primary', onclick: () => this.runScan() }, icon('flask'), `Run the scan (${3 * ens.nRes} mutants)`));
      return;
    }
    const sc = this.scan;
    sec.append(el('div', { class: 'desc', text: `Wild type: barrier ${fmt(sc.barrier_wt, 1)} kcal/mol along the best route, ΔG(A→B) ${fmtSigned(sc.dg_wt, 1)} kcal/mol. φ is shown where |ΔΔG| ≥ ${sc.phi_min_ddg} kcal/mol.` }));
    sec.append(this._scanHeat(sc));
    const list = (title, items) => el('div', { class: 'sw-list' }, el('div', { class: 'section-title', style: { marginTop: '10px' } }, title),
      items.length ? el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'Mutation' }), el('th', { class: 'num', text: 'ΔΔG‡' }), el('th', { class: 'num', text: 'ΔΔG' }), el('th', { class: 'num', text: 'φ' }), el('th', { class: 'num', text: 'Rate' }))),
        el('tbody', {}, ...items.map(x => { const tr = el('tr', { class: 'clickable' }, el('td', {}, el('b', { text: `${x.from}${ens.residues[x.pos].num}${x.to}` })),
          el('td', { class: 'num', text: fmtSigned(x.ddg_ts, 1) }), el('td', { class: 'num', text: fmtSigned(x.ddg_eq, 1) }), el('td', { class: 'num', text: x.phi === null ? '–' : fmt(x.phi, 2) }),
          el('td', { class: 'num', text: `×${rateFactor(x.ddg_ts, r.kT)}` })); tr.addEventListener('click', () => { app.select([x.pos], 'set'); }); return tr; })))
        : el('div', { class: 'desc', text: 'none' }));
    sec.append(list('Speed up switching', sc.faster), list('Slow down switching', sc.slower),
      el('div', { class: 'row', style: { marginTop: '10px', gap: '8px', flexWrap: 'wrap' } },
        el('button', { class: 'btn small', onclick: () => this.colorBy('kinetic') }, icon('layers'), 'Colour 3D by kinetic sensitivity'),
        el('button', { class: 'btn small ghost', onclick: () => this.copyScan() }, icon('copy'), 'Copy scan as table')),
      el('div', { class: 'desc', text: 'Route kept and its energies re-evaluated for each mutant; a mutation that makes a pair non-canonical removes it. Rate = factor on the A → B rate (exp(−ΔΔG‡/RT)).' }));
  }

  _scanHeat(sc) {
    const app = this.app, ens = app.ens, n = ens.nRes;
    const dark = document.documentElement.dataset.theme === 'dark';
    const W = Math.max(240, this.root.clientWidth - 28);
    const lab = 16, cw = Math.max(4, Math.min(14, (W - lab) / n)), gap = cw >= 7 ? 1 : 0, ch = 13, phiH = 22;
    const H = 4 * ch + phiH + 34;
    const width = lab + cw * n;
    const s = svg('svg', { class: 'chart sw-heat', width, height: H, viewBox: `0 0 ${width} ${H}`, role: 'img', 'aria-label': 'Kinetic mutation scan heat map' });
    'ACGU'.split('').forEach((b, row) => s.append(svg('text', { x: 2, y: 10 + row * ch + ch / 2, 'font-size': 10, fill: 'var(--text-2)', text: b })));
    s.append(svg('text', { x: 2, y: 4 * ch + 20, 'font-size': 10, fill: 'var(--text-2)', text: 'φ' }));
    for (let i = 0; i < n; i++) {
      const wt = ens.sequence[i];
      'ACGU'.split('').forEach((b, row) => {
        const x = lab + i * cw, y = 4 + row * ch;
        if (b === wt) { s.append(svg('rect', { x: x + gap / 2, y: y + gap / 2, width: cw - gap, height: ch - gap, fill: 'var(--surface)', stroke: 'var(--axis)', 'stroke-width': 0.5 })); return; }
        const v = sc.matrix[i][b];
        s.append(svg('rect', { x: x + gap / 2, y: y + gap / 2, width: cw - gap, height: ch - gap, fill: divergingColor(v.ddg_ts, dark), class: 'sw-cell', 'data-i': i, 'data-b': b }));
      });
      const phi = sc.positions[i].phi;
      if (phi !== null) {
        const h = Math.max(1, Math.min(1.2, Math.abs(phi)) / 1.2 * (phiH - 4));
        s.append(svg('rect', { x: lab + i * cw + gap / 2, y: 4 * ch + 8 + (phiH - 4 - h), width: cw - gap, height: h, fill: phi < 0 ? 'var(--bad)' : 'var(--text-2)' }));
      }
      if (i === 0 || (i + 1) % 10 === 0) s.append(svg('text', { x: lab + i * cw + cw / 2, y: H - 4, 'text-anchor': 'middle', 'font-size': 9, fill: 'var(--text-3)', text: ens.residues[i].num }));
    }
    s.addEventListener('pointermove', e => {
      const c = e.target.closest('.sw-cell');
      if (!c) { tooltip.hide(); app.hoverResidues([]); return; }
      const i = +c.dataset.i, b = c.dataset.b, v = sc.matrix[i][b];
      app.hoverResidues([i]);
      tooltip.show([tipRow(`${ens.sequence[i]}${ens.residues[i].num}${b}`, `ΔΔG‡ ${fmtSigned(v.ddg_ts, 1)} kcal/mol`),
        el('div', { class: 'tl', text: `ΔΔG ${fmtSigned(v.ddg_eq, 1)} kcal/mol · φ ${v.phi === null ? '–' : fmt(v.phi, 2)} · rate ×${rateFactor(v.ddg_ts, this.result.kT)}${v.breaks ? ' · breaks a pair on the route' : ''}` })], e.clientX, e.clientY);
    });
    s.addEventListener('pointerleave', () => { tooltip.hide(); app.hoverResidues([]); });
    s.addEventListener('click', e => { const c = e.target.closest('.sw-cell'); if (c) app.select([+c.dataset.i], 'set'); });
    const legend = el('div', { class: 'sw-legend' }, el('span', { text: 'faster' }),
      el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${[-4, -2, 0, 2, 4].map(v => divergingColor(v, dark)).join(',')})` } }),
      el('span', { text: 'slower' }), el('span', { class: 'muted', text: '±4 kcal/mol · bars: median φ per position' }));
    return el('div', {}, el('div', { class: 'sw-heatbox' }, s), legend);
  }

  colorBy(kind) {
    const app = this.app, ens = app.ens, r = this.result;
    if (!r) return;
    const n = ens.nRes;
    if (kind === 'role') {
      const role = r.timing.map(t => t.role);
      const mod = r.motion?.module_of || [];
      const coupled = r.coupling ? new Set(r.motion.modules[r.coupling.module].residues) : new Set();
      registerScheme({
        id: 'pathway-role', label: 'Mechanism (Pathway)',
        residueColors: () => [...Array(n).keys()].map(i => role[i] === 'unit' ? '#eb6834' : role[i] === 'nucleus' ? '#2a78d6' : coupled.has(i) ? '#f5c6a8' : mod[i] >= 0 ? '#d6d4cc' : '#9b998f'),
        legend: () => ({ title: 'Transition mechanism', items: [{ color: '#eb6834', label: 'cooperative unit (opens before the saddle)' }, { color: '#2a78d6', label: 'nucleus of B (forms after it)' },
          { color: '#f5c6a8', label: 'rate-limiting module' }, { color: '#d6d4cc', label: 'other rigid modules' }, { color: '#9b998f', label: 'flexible (hinges, loops)' }] }),
      });
      this._apply('pathway-role');
    } else if (kind === 'order') {
      const t = r.timing.map(x => x.first);
      registerScheme({
        id: 'pathway-order', label: 'Event order (Pathway)',
        residueColors: () => t.map(v => (v === null ? '#c4c2ba' : sampleRamp(VIRIDIS, v))),
        legend: () => ({ title: 'When the pairing changes on the route A → B', ramp: VIRIDIS, min: 'early', max: 'late' }),
      });
      this._apply('pathway-order');
    } else if (kind === 'modules') {
      const mod = r.motion?.module_of || [];
      registerScheme({
        id: 'pathway-modules', label: 'Rigid modules (Pathway)',
        residueColors: () => mod.map(m => (m >= 0 ? SLOTS[m % SLOTS.length] : '#9b998f')),
        legend: () => ({ title: 'Rigid modules between A and B', items: [...(r.motion?.modules || []).map(m => ({ color: SLOTS[m.id % SLOTS.length], label: `M${m.id + 1} (${m.size} nt)` })), { color: '#9b998f', label: 'flexible' }] }),
      });
      this._apply('pathway-modules');
    } else if (kind === 'kinetic' && this.scan) {
      const v = this.scan.positions.map(p => p.max_abs_ddg_ts);
      const max = Math.max(1e-6, ...v);
      registerScheme({
        id: 'pathway-kinetic', label: 'Kinetic sensitivity (Pathway)',
        residueColors: () => v.map(x => sampleRamp(ORANGE_RAMP, x / max)),
        legend: () => ({ title: 'Largest |ΔΔG‡| of a point mutation', ramp: ORANGE_RAMP, min: '0', max: `${fmt(max, 1)} kcal/mol` }),
      });
      this._apply('pathway-kinetic');
    }
  }

  _apply(id) {
    this.app.refreshSchemes();
    this.app.setColorScheme(id);
  }

  play(staged) {
    const app = this.app, ens = app.ens, r = this.result;
    if (!r?.motion) return;
    const mech = buildMechanism(ens.coords[r.a], ens.coords[r.b], ens.atomRes, r.motion, r.timing, { staged, names: ens.atomNames });
    app.playMechanism(r.a, r.b, mech, { duration: staged ? 6 : 3.5 });
  }

  copyScan() {
    const ens = this.app.ens, sc = this.scan;
    if (!sc) return;
    const rows = ['position\twild type\tmutant\tddG_TS (kcal/mol)\tddG_eq (kcal/mol)\tphi\tbreaks route pair'];
    sc.matrix.forEach((row, i) => { for (const [b, v] of Object.entries(row)) rows.push([ens.residues[i].num, ens.sequence[i], b, v.ddg_ts, v.ddg_eq, v.phi ?? '', v.breaks].join('\t')); });
    navigator.clipboard?.writeText(rows.join('\n')).then(() => toast(`${rows.length - 1} mutations copied as a tab-separated table.`, 'ok'));
  }

  _markCurrent() {
    if (!this.result || !this._mapPos || !this.mapSvg) return;
    // the map is cheap to redraw; keep the current-conformation ring up to date
    if (this.visible()) this.render();
  }

  _exportBtn(getSvg, name, studio = null) {
    return el('button', { class: 'ss-btn chart-export', title: 'Export this chart (PDF, SVG, PNG, TIFF) or open it in the plot studio', 'aria-label': `Export ${name}`,
      onclick: e => popupMenu(e.currentTarget, [...chartMenuItems(getSvg(), `${this.app.ens?.name || 'rcv'} ${name}`, toast),
        ...(studio ? [{ sep: true }, { label: 'Open in the plot studio (journal sizes)…', action: () => this.app.openPlotStudio(studio) }] : [])]) }, icon('download'));
  }
}

export function rangesText(ens, residues) {
  const rs = [...residues].sort((a, b) => a - b);
  const out = [];
  let s = rs[0], p = rs[0];
  for (let i = 1; i <= rs.length; i++) {
    if (i < rs.length && rs[i] === p + 1) { p = rs[i]; continue; }
    out.push(s === p ? ens.resLabel(s) : `${ens.resLabel(s)}–${ens.resLabel(p)}`);
    s = p = rs[i];
  }
  return out.join(', ');
}

function rateFactor(ddg, kT) {
  const f = Math.exp(-ddg / kT);
  if (f >= 100 || f < 0.01) return f.toExponential(0);
  return fmt(f, f < 1 ? 2 : 1);
}
