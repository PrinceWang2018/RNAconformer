// Energy side panel: current value and change, energy profile over the conformation order,
// energy landscape (funnel vs RMSD, or similarity map), per-term breakdown and base-pair changes.
// Markers follow conformation transitions frame by frame.
import { el, svg, fmt, fmtVal, fmtSigned, loadPref, savePref } from '../util.js';
import { icon } from '../icons.js';
import { linear, extent, yAxis, xAxis } from './charts.js';
import { chartMenuItems } from '../export/chartexport.js';
import { tooltip, tipRow, select, popupMenu, toast } from '../ui.js';
import { ArcDiagram } from './arcs.js';
import { BLUE_RAMP, sampleRamp } from '../colors.js';
import { drawMap, landscapeOf, loadLandscape, qualityText } from '../landscape.js';

const GROUPS = [['engine', 'Prediction engine'], ['confidence', 'Confidence model'], ['evaluator', 'Evaluator (any model)'], ['descriptor', 'Structure descriptors'], ['data', 'Probing data'], ['measurement', 'Measurements'], ['plugin', 'Plugins']];

export function direction(metric, delta) {
  if (delta === null || !metric || metric.lower_is_better === null || metric.lower_is_better === undefined || Math.abs(delta) < 1e-9) return 'neutral';
  const better = metric.lower_is_better ? delta < 0 : delta > 0;
  return better ? 'good' : 'bad';
}

export class EnergyPanel {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    this.landMode = loadPref('landMode', 'funnel');
    this.metricSel = el('select', { class: 'select', 'aria-label': 'Energy metric' });
    this.metricSel.addEventListener('change', () => app.setMetric(this.metricSel.value));
    this.hero = el('div', { class: 'hero' });
    this.desc = el('div', { class: 'desc' });
    this.profile = svg('svg', { class: 'chart', role: 'img' });
    this.land = svg('svg', { class: 'chart', role: 'img' });
    this.landCaption = el('div', { class: 'desc' });
    this.seg = el('div', { class: 'seg' },
      el('button', { 'data-m': 'funnel', text: 'Funnel', title: 'Energy versus C4′ RMSD to the best conformation' }),
      el('button', { 'data-m': 'map', text: 'Map', title: 'Topology-aware map: one disk per secondary-structure topology, distances by C4′ RMSD, coloured by energy' }),
      el('button', { 'data-m': 'mds', text: 'MDS', title: 'Similarity map (classical MDS of pairwise RMSD), coloured by energy' }));
    this.seg.addEventListener('click', e => {
      const m = e.target.closest('button')?.dataset.m;
      if (m) { this.landMode = m; savePref('landMode', m); this.render(); }
    });
    this.terms = el('div', { class: 'terms' });
    this.termsSection = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Energy terms', el('span', { class: 'grow' }),
      el('span', { class: 'note', text: 'bar = current · tick = previous' })), this.terms);
    const arcHost = el('div');
    this.arcs = new ArcDiagram(app, arcHost);
    this.metricRow = el('div', { class: 'metric-row' }, this.metricSel);
    this.charts = el('div');
    root.append(this.metricRow, el('div', { class: 'section' }, this.hero, this.desc), this.charts);
    this.charts.append(
      el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Energy profile', el('span', { class: 'grow' }),
        el('span', { class: 'note', id: 'profile-note' }), this._exportBtn(() => this.profile, 'Energy profile', 'profile')), this.profile),
      el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Energy landscape', el('span', { class: 'grow' }), this.seg, this._exportBtn(() => this.land, 'Energy landscape', () => (this.landMode === 'map' ? 'landscape' : 'funnel'))), this.land, this.landCaption),
      this.termsSection,
      el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Base pairs', el('span', { class: 'grow' }),
        el('span', { class: 'note', text: 'current above · previous below' }), this._exportBtn(() => this.arcs.svg, 'Base pairs', 'ss2d')), arcHost),
    );
    for (const ev of ['ensemble', 'conformation', 'metric', 'order', 'values']) app.on(ev, () => this.render());
    app.on('morph', m => this.updateMarkers(m));
    new ResizeObserver(() => { const w = this.root.clientWidth; if (Math.abs(w - (this._w || 0)) > 4) { this._w = w; this.render(); } }).observe(root);
  }

  width() { return Math.max(240, this.root.clientWidth - 28); }

  // small export button beside a chart title: PDF / SVG / PNG / TIFF, or open the plot studio
  _exportBtn(getSvg, name, studioType) {
    return el('button', { class: 'ss-btn chart-export', title: 'Export this chart (PDF, SVG, PNG, TIFF) or open it in the plot studio', 'aria-label': `Export ${name}`,
      onclick: e => popupMenu(e.currentTarget, [...chartMenuItems(getSvg(), `${this.app.ens?.name || 'rcv'} ${name}`, toast), { sep: true },
        { label: 'Open in the plot studio (journal sizes)…', action: () => this.app.openPlotStudio(typeof studioType === 'function' ? studioType() : studioType) }]) }, icon('download'));
  }

  render() {
    const app = this.app, ens = app.ens;
    this.charts.hidden = this.metricRow.hidden = !ens;
    if (!ens) {
      this.metricSel.replaceChildren();
      this.hero.replaceChildren(el('div', { class: 'label', text: 'No ensemble loaded' }), el('div', { class: 'ctx', text: 'Open structures or run a prediction to see energies.' }));
      for (const s of [this.profile, this.land]) s.replaceChildren();
      this.terms.replaceChildren();
      this.desc.textContent = '';
      this.arcs.render(this.width());
      return;
    }
    this._metricOptions();
    this._hero();
    this._profile();
    this._landscape();
    this._terms();
    this.arcs.render(this.width());
    for (const b of this.seg.children) b.classList.toggle('on', b.dataset.m === this.landMode);
  }

  _metricOptions() {
    const ens = this.app.ens, key = this.app.state.metric;
    const sig = ens.id + ':' + ens.metrics.length;
    if (this._sig !== sig) {
      this._sig = sig;
      const opts = [];
      for (const [g, label] of GROUPS) {
        const ms = ens.metrics.filter(m => (m.group || 'descriptor') === g && ens.conformations.some(c => typeof c.values[m.key] === 'number'));
        if (ms.length) opts.push({ group: label, options: ms.map(m => ({ value: m.key, label: m.unit ? `${m.label} (${m.unit})` : m.label })) });
      }
      const s = select(opts, key);
      this.metricSel.replaceChildren(...s.children);
    }
    this.metricSel.value = key;
  }

  _hero() {
    const app = this.app, ens = app.ens, key = app.state.metric, m = ens.metricByKey[key] || { label: key };
    const k = app.state.conf, p = app.state.prev;
    const v = ens.value(k, key);
    const pv = p !== null && p !== undefined && p !== k ? ens.value(p, key) : null;
    const d = v !== null && pv !== null ? v - pv : null;
    const dir = direction(m, d);
    const sorted = ens.order('metric', key);
    const rank = sorted.indexOf(k) + 1;
    this.hero.replaceChildren(
      el('div', { class: 'label' }, el('span', { text: m.label }), ` · ${ens.shortId(k)}`),
      el('div', { class: 'value' }, fmtVal(v, m.unit), m.unit ? el('small', { text: m.unit }) : null),
      d === null ? el('div', { class: 'ctx', text: p === null || p === k ? 'Switch conformation to see how the energy changes' : 'no value for the previous conformation' })
        : el('div', { class: `delta ${dir}` }, icon(d < 0 ? 'down' : 'up'), `${fmtSigned(d)} ${m.unit || ''} vs ${ens.shortId(p)}`,
          el('span', { class: 'muted', style: { fontWeight: 400 }, text: dir === 'good' ? ' · better' : dir === 'bad' ? ' · worse' : '' })),
      el('div', { class: 'ctx', text: m.lower_is_better === null || m.lower_is_better === undefined ? `${ens.n} conformation${ens.n === 1 ? '' : 's'}` : `Rank ${rank} of ${ens.n} by this metric (${m.lower_is_better ? 'lower is better' : 'higher is better'})` }),
    );
    this.desc.textContent = m.description || '';
  }

  // ------------------------------------------------------------------ profile
  _profile() {
    const app = this.app, ens = app.ens, key = app.state.metric;
    const order = app.order;
    const W = this.width(), H = 172, ml = 56, mr = 10, mt = 10, mb = 38;
    const s = this.profile;
    s.replaceChildren();
    s.setAttribute('viewBox', `0 0 ${W} ${H}`);
    s.setAttribute('height', H);
    const vals = order.map(k => ens.value(k, key));
    const n = order.length;
    const x = linear(0, Math.max(1, n - 1), ml + 6, W - mr - 6);
    const [lo, hi] = extent(vals, 0.12);
    const y = linear(lo, hi, H - mb, mt);
    const m = ens.metricByKey[key] || {};
    s.setAttribute('aria-label', `${m.label || key} for each conformation in ${app.state.orderMode} order`);
    document.getElementById('profile-note').textContent = { rank: 'rank order', metric: 'sorted by energy', path: 'similarity path', input: 'input order' }[app.state.orderMode] || '';
    const g = svg('g');
    yAxis(g, y, ml, W - mr, { ticks: 4, title: `${m.label || key}${m.unit ? ` (${m.unit})` : ''}`, grid: true });
    // x ticks: positions in the order
    const every = Math.max(1, Math.ceil(n / 8));
    const labs = [];
    for (let i = 0; i < n; i += every) {
      if (i > 0 && n - 1 - i < every * 0.6 && (n - 1) % every) continue;   // keep room for the final label
      labs.push([x(i), String(i + 1)]);
    }
    if (!labs.length || labs[labs.length - 1][1] !== String(n)) labs.push([x(n - 1), String(n)]);
    xAxis(g, ml, W - mr, H - mb, labs, { title: `Conformation (${{ rank: 'rank', metric: 'energy', path: 'similarity path', input: 'input' }[app.state.orderMode] || 'order'} order)` });
    const pts = vals.map((v, i) => (v === null ? null : [x(i), y(v)]));
    const path = pts.filter(Boolean).map(p => p.join(',')).join(' ');
    if (path) g.append(svg('polyline', { points: path, fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 1.6, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
    pts.forEach(p => { if (p) g.append(svg('circle', { cx: p[0], cy: p[1], r: 3.5, fill: 'var(--series-1)', stroke: 'var(--surface)', 'stroke-width': 1.5 })); });
    const cur = order.indexOf(app.state.conf), prev = app.state.prev !== null ? order.indexOf(app.state.prev) : -1;
    if (prev >= 0 && prev !== cur && pts[prev]) g.append(svg('circle', { cx: pts[prev][0], cy: pts[prev][1], r: 7.5, fill: 'none', stroke: 'var(--text-3)', 'stroke-width': 1.5, 'stroke-dasharray': '2.5 2' }));
    if (cur >= 0 && pts[cur]) g.append(svg('circle', { cx: pts[cur][0], cy: pts[cur][1], r: 7.5, fill: 'none', stroke: 'var(--text)', 'stroke-width': 2 }));
    this.pMarker = svg('circle', { r: 5.5, fill: 'var(--series-2)', stroke: 'var(--surface)', 'stroke-width': 2, visibility: 'hidden' });
    g.append(this.pMarker);
    this.pScale = { x, y, order, key };
    // hover / click bands
    const band = n > 1 ? (x(1) - x(0)) : W;
    order.forEach((k, i) => {
      const r = svg('rect', { class: 'hit', x: x(i) - band / 2, y: mt, width: band, height: H - mt - mb + 6 });
      r.addEventListener('pointermove', e => {
        const v = vals[i], dv = v !== null && ens.value(app.state.conf, key) !== null ? v - ens.value(app.state.conf, key) : null;
        tooltip.show([tipRow(`${fmtVal(v, m.unit)} ${m.unit || ''}`, ens.shortId(k), 'var(--series-1)'), el('div', { class: 'tl', text: k === app.state.conf ? 'current conformation' : `${fmtSigned(dv)} vs current · click to morph` })], e.clientX, e.clientY);
      });
      r.addEventListener('pointerenter', () => app.previewConformation(k));
      r.addEventListener('pointerleave', () => { tooltip.hide(); app.previewConformation(null); });
      r.addEventListener('click', () => app.goTo(k));
      g.append(r);
    });
    s.append(g);
  }

  // ------------------------------------------------------------------ landscape
  _landscape() {
    const app = this.app, ens = app.ens, key = app.state.metric;
    const m = ens.metricByKey[key] || {};
    const W = this.width(), H = 212, ml = 56, mr = 12, mt = 10, mb = 38;
    const s = this.land;
    s.replaceChildren();
    s.setAttribute('viewBox', `0 0 ${W} ${H}`);
    s.setAttribute('height', H);
    const idx = [...Array(ens.n).keys()];
    const vals = idx.map(k => ens.value(k, key));
    const finite = vals.filter(v => v !== null);
    let best = 0;
    if (finite.length) {
      const sorted = ens.order('metric', key);
      best = sorted[0];
    }
    const g = svg('g');
    let px, py;
    const cur = app.state.conf, prev = app.state.prev;
    if (this.landMode === 'map' || this.landMode === 'mds') {
      const L = landscapeOf(ens);
      if (!L && this.landMode === 'map') loadLandscape(app, ens).then(x => { if (x && app.ens === ens) this.render(); });
      const data = L || { positions: ens.embedding, mds: ens.embedding, classes: [], edges: [] };
      const mode = this.landMode === 'map' && L ? 'map' : 'mds';
      g.append(svg('rect', { x: ml - 6, y: mt - 2, width: W - ml - mr + 12, height: H - mt - mb + 4, rx: 6, fill: 'none', stroke: 'var(--grid)' }));
      const [lo, hi] = extent(finite, 0);
      const col = v => {
        if (v === null) return 'var(--axis)';
        let t = (v - lo) / ((hi - lo) || 1);
        if (m.lower_is_better) t = 1 - t;
        return sampleRamp(BLUE_RAMP.slice(2), t);
      };
      this._pointColor = k => col(vals[k]);
      const weights = ens.probing?.fit?.reweighting?.weights;
      const radius = weights && this.app.state.metric === 'probe_weight' ? k => 2.5 + 7 * Math.sqrt(weights[k] / Math.max(...weights)) : () => 4.5;
      const r = drawMap(g, ens, data, { x0: ml - 2, y0: mt + 2, x1: W - mr - 2, y1: H - mb - 2 }, { mode, color: k => col(vals[k]), radius, current: cur, prev, fontSize: 9.5 });
      px = r.px; py = r.py;
      g.append(svg('text', { x: ml - 4, y: H - mb + 26, text: m.lower_is_better === false ? 'darker = higher' : 'darker = lower (better)' }));
      if (mode === 'map') this.landCaption.replaceChildren(el('span', { text: 'Topology-aware map: each disk is one secondary-structure topology (abstract shape, label T1, T2, …); inside a disk distances follow C4′ RMSD and 2D states stay together. Click a point to morph to it.' }), ' ', el('span', { text: `${qualityText(L)}.` }));
      else this.landCaption.textContent = 'Each point is a conformation; distances approximate pairwise C4′ RMSD (classical MDS). Click a point to morph to it.';
    } else {
      const rm = idx.map(k => ens.rmsd[best][k]);
      const x = linear(0, Math.max(1, ...rm) * 1.08, ml + 4, W - mr);
      const [lo, hi] = extent(vals, 0.12);
      const y = linear(lo, hi, H - mb, mt);
      yAxis(g, y, ml, W - mr, { ticks: 4, title: `${m.label || key}${m.unit ? ` (${m.unit})` : ''}`, grid: true });
      const xt = [0, ...[2, 5, 10, 15, 20, 30, 40, 60].filter(v => v < x.domain[1])].slice(0, 7);
      xAxis(g, ml, W - mr, H - mb, xt.map(v => [x(v), String(v)]), { title: `C4′ RMSD to ${ens.shortId(best)} (Å)` });
      px = k => x(rm[k]);
      py = k => (vals[k] === null ? H - mb : y(vals[k]));
      this._pointColor = () => 'var(--series-1)';
      this.landCaption.textContent = `Funnel plot: ${m.label || key} against structural distance from the best conformation by this metric (${ens.shortId(best)}).`;
      if (prev !== null && prev !== cur) {
        g.append(svg('path', { d: `M ${px(prev)} ${py(prev)} L ${px(cur)} ${py(cur)}`, fill: 'none', stroke: 'var(--text-3)', 'stroke-width': 1.5, 'stroke-dasharray': '4 3' }));
      }
      for (const k of idx) g.append(svg('circle', { cx: px(k), cy: py(k), r: 4.5, fill: this._pointColor(k), stroke: 'var(--surface)', 'stroke-width': 2 }));
      g.append(svg('path', { d: diamond(px(best), py(best), 6), fill: 'none', stroke: 'var(--text-2)', 'stroke-width': 1.5 }));
      if (prev !== null && prev !== cur) g.append(svg('circle', { cx: px(prev), cy: py(prev), r: 8, fill: 'none', stroke: 'var(--text-3)', 'stroke-width': 1.5, 'stroke-dasharray': '2.5 2' }));
      g.append(svg('circle', { cx: px(cur), cy: py(cur), r: 8, fill: 'none', stroke: 'var(--text)', 'stroke-width': 2 }));
    }
    this.lMarker = svg('circle', { r: 5.5, fill: 'var(--series-2)', stroke: 'var(--surface)', 'stroke-width': 2, visibility: 'hidden' });
    g.append(this.lMarker);
    this.lScale = { px, py };
    // nearest-point hover layer
    const hit = svg('rect', { class: 'hit', x: 0, y: 0, width: W, height: H });
    const nearest = e => {
      const r = s.getBoundingClientRect();
      const sx = (e.clientX - r.left) * W / r.width, sy = (e.clientY - r.top) * H / r.height;
      let bk = -1, bd = 24;
      for (const k of idx) { const d = Math.hypot(px(k) - sx, py(k) - sy); if (d < bd) { bd = d; bk = k; } }
      return bk;
    };
    hit.addEventListener('pointermove', e => {
      const k = nearest(e);
      if (k < 0) { tooltip.hide(); app.previewConformation(null); hit.style.cursor = 'default'; return; }
      hit.style.cursor = 'pointer';
      app.previewConformation(k);
      tooltip.show([tipRow(`${fmtVal(vals[k], m.unit)} ${m.unit || ''}`, ens.shortId(k)),
        el('div', { class: 'tl', text: `RMSD ${fmt(ens.rmsd[cur][k])} Å to current${ens.conformations[k].cluster ? ' · cluster ' + ens.conformations[k].cluster : ''}` })], e.clientX, e.clientY);
    });
    hit.addEventListener('pointerleave', () => { tooltip.hide(); app.previewConformation(null); });
    hit.addEventListener('click', e => { const k = nearest(e); if (k >= 0) app.goTo(k); });
    g.append(hit);
    s.setAttribute('aria-label', this.landMode === 'map' ? 'Similarity map of conformations coloured by energy' : 'Energy versus RMSD funnel plot');
    s.append(g);
  }

  // ------------------------------------------------------------------ term breakdown
  _terms() {
    const app = this.app, ens = app.ens, key = app.state.metric;
    const m = ens.metricByKey[key];
    const terms = (m && m.terms) || [];
    this.termsSection.hidden = !terms.length;
    this.terms.replaceChildren();
    if (!terms.length) return;
    const k = app.state.conf, p = app.state.prev;
    const rows = terms.map(t => ({ t, v: ens.value(k, t), pv: p !== null && p !== k ? ens.value(p, t) : null }))
      .filter(r => r.v !== null && (Math.abs(r.v) > 1e-6 || (r.pv !== null && Math.abs(r.pv) > 1e-6)));
    let lo = 0, hi = 0;
    for (const r of rows) for (const v of [r.v, r.pv]) if (v !== null) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    if (hi - lo < 1e-9) hi = lo + 1;
    for (const r of rows) {
      const info = ens.terms[r.t] || { label: r.t };
      const BW = 100;
      const x = v => (v - lo) / (hi - lo) * BW;
      const bar = svg('svg', { viewBox: `0 0 ${BW} 14`, preserveAspectRatio: 'none' });
      const x0 = x(0), x1 = x(r.v);
      bar.append(svg('rect', { x: Math.min(x0, x1), y: 3, width: Math.max(0.8, Math.abs(x1 - x0)), height: 8, rx: 1.5, fill: 'var(--series-1)' }));
      if (lo < 0) bar.append(svg('line', { x1: x0, x2: x0, y1: 0, y2: 14, stroke: 'var(--axis)', 'stroke-width': 0.8, 'vector-effect': 'non-scaling-stroke' }));
      if (r.pv !== null) bar.append(svg('line', { x1: x(r.pv), x2: x(r.pv), y1: 0, y2: 14, stroke: 'var(--text)', 'stroke-width': 2, 'vector-effect': 'non-scaling-stroke' }));
      const d = r.pv !== null ? r.v - r.pv : null;
      const dir = direction({ lower_is_better: m.lower_is_better }, d);
      const name = el('span', { class: 'tname', text: info.label, title: info.description || '' });
      this.terms.append(name, bar, el('span', { class: 'tval' }, fmtVal(r.v), d !== null && Math.abs(d) > 5e-4 ? el('span', { class: `d ${dir}`, text: fmtSigned(d) }) : null));
    }
    const total = ens.value(k, key), ptotal = p !== null && p !== k ? ens.value(p, key) : null;
    const dt = total !== null && ptotal !== null ? total - ptotal : null;
    this.terms.append(el('span', { class: 'tname', style: { fontWeight: 600, color: 'var(--text)' }, text: 'Total' }), el('span'),
      el('span', { class: 'tval', style: { fontWeight: 600 } }, fmtVal(total), dt !== null ? el('span', { class: `d ${direction(m, dt)}`, text: fmtSigned(dt) }) : null));
  }

  // ------------------------------------------------------------------ live markers during transitions
  updateMarkers(mo) {
    const app = this.app, ens = app.ens;
    if (!ens || !this.pScale || !this.lScale) return;
    if (!mo || mo.done) {
      this.pMarker?.setAttribute('visibility', 'hidden');
      this.lMarker?.setAttribute('visibility', 'hidden');
      return;
    }
    const { from, to, t } = mo;
    const { x, y, order, key } = this.pScale;
    const ia = order.indexOf(from), ib = order.indexOf(to);
    const va = ens.value(from, key), vb = ens.value(to, key);
    if (ia >= 0 && ib >= 0 && va !== null && vb !== null) {
      this.pMarker.setAttribute('cx', x(ia + (ib - ia) * t));
      this.pMarker.setAttribute('cy', y(va + (vb - va) * t));
      this.pMarker.setAttribute('visibility', 'visible');
    }
    const { px, py } = this.lScale;
    this.lMarker.setAttribute('cx', px(from) + (px(to) - px(from)) * t);
    this.lMarker.setAttribute('cy', py(from) + (py(to) - py(from)) * t);
    this.lMarker.setAttribute('visibility', 'visible');
  }
}

function diamond(cx, cy, r) {
  return `M ${cx} ${cy - r} L ${cx + r} ${cy} L ${cx} ${cy + r} L ${cx - r} ${cy} Z`;
}
