// Plot studio: publication plots of an ensemble, drawn as SVG in points (sizes in millimetres, journal column widths)
// with the figure font, exported as vector PDF (embedded font) / SVG, or PNG / TIFF at 300–1200 dpi.
//   RMSD matrix (optional average-linkage clustering with dendrogram) · RMSD to a reference · per-residue RMSF ·
//   energy vs RMSD (funnel) · energy profile · 2D structure
import { el, svg, fmt, niceTicks, slug, FIG_FONT, ensureFigureFont } from './util.js';
import { dialog, field, select, toast } from './ui.js';
import { BLUE_RAMP, ORANGE_RAMP, sampleRamp, SLOTS, pairedFraction } from './colors.js';
import { kabsch, applyTransform } from './ensemble.js';
import { layoutFor, drawStructure2D, orientForBox } from './ss2dlayout.js';
import { exportChart, copyChart } from './export/chartexport.js';
import { upgma } from './cluster.js';
export { upgma };

const MM = 72 / 25.4;
const INK = '#1a1a1a', SUB = '#4d4d4d', AXIS = '#1a1a1a', GRID = '#e2e2e2';
const RAMPS = {
  blues: BLUE_RAMP, oranges: ORANGE_RAMP,
  grays: ['#f2f2f2', '#d9d9d9', '#bdbdbd', '#969696', '#737373', '#525252', '#252525'],
  purples: ['#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#4a1486'],
  greens: ['#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#238b45', '#005a32'],
};
export const PLOTS = [
  ['matrix', 'Pairwise RMSD matrix'], ['toref', 'RMSD to a reference'], ['rmsf', 'Per-residue flexibility (RMSF)'],
  ['funnel', 'Energy vs RMSD (funnel)'], ['profile', 'Energy profile'], ['ss2d', '2D structure (vector)'],
];

let _mctx = null;
const textW = (s, size, bold = false) => { _mctx = _mctx || document.createElement('canvas').getContext('2d'); _mctx.font = `${bold ? 'bold' : 'normal'} ${size}px ${FIG_FONT}`; return _mctx.measureText(String(s)).width; };
const T = (x, y, s, size, o = {}) => svg('text', { x, y, 'font-size': size, 'font-family': FIG_FONT, fill: o.fill || INK, 'text-anchor': o.anchor || 'start', 'font-weight': o.bold ? 'bold' : null, transform: o.rotate ? `rotate(${o.rotate} ${x} ${y})` : null, text: s });
// Legend items ({ name, color, kind: 'line' | 'dot' }) laid out in rows that wrap inside [x0, x1]: never clipped.
function legendLayout(items, x0, x1, fs) {
  const size = fs * 0.85, rows = [[]];
  let x = x0;
  for (const it of items) {
    const w = 15 + textW(it.name, size);
    if (x > x0 && x + w > x1) { rows.push([]); x = x0; }
    rows[rows.length - 1].push({ ...it, x });
    x += w + 12;
  }
  return { rows, size, height: items.length ? rows.length * fs * 1.3 + fs * 0.35 : 0 };
}
function drawLegend(g, lay, top) {
  lay.rows.forEach((row, r) => {
    const cy = top + lay.size * 0.75 + r * lay.size / 0.85 * 1.3;
    for (const it of row) {
      if (it.kind === 'dot') g.append(svg('circle', { cx: it.x + 5, cy, r: 2.6, fill: it.color }));
      else g.append(svg('line', { x1: it.x, y1: cy, x2: it.x + 11, y2: cy, stroke: it.color, 'stroke-width': 1.4 }));
      g.append(svg('text', { x: it.x + 15, y: cy + lay.size * 0.35, 'font-size': lay.size, 'font-family': FIG_FONT, fill: INK, text: it.name }));
    }
  });
}
const L = (x1, y1, x2, y2, o = {}) => svg('line', { x1, y1, x2, y2, stroke: o.stroke || AXIS, 'stroke-width': o.width ?? 0.6, 'stroke-linecap': o.cap || 'butt', 'stroke-dasharray': o.dash || null });

// ------------------------------------------------------------------ data helpers
function c4Atoms(ens) { return [...ens.trace].filter(a => a >= 0); }
function rmsdFit(ens, A, B, atoms) {
  const Bc = Float32Array.from(B);
  const f = kabsch(Bc, A, atoms);
  applyTransform(Bc, f.R, f.t);
  let s = 0;
  for (const i of atoms) { const dx = A[3 * i] - Bc[3 * i], dy = A[3 * i + 1] - Bc[3 * i + 1], dz = A[3 * i + 2] - Bc[3 * i + 2]; s += dx * dx + dy * dy + dz * dz; }
  return Math.sqrt(s / atoms.length);
}
function meanCoords(ens) {
  const m = new Float32Array(ens.coords[0].length);
  for (const c of ens.coords) for (let i = 0; i < m.length; i++) m[i] += c[i] / ens.n;
  return m;
}

// ------------------------------------------------------------------ axes
function axes(g, { x0, x1, y0, y1, xs, ys, xticks, yticks, xlabel, ylabel, fs, grid, xfmt = v => fmt(v, 1), yfmt = v => fmt(v, 1), xTickLabels = null }) {
  if (grid) for (const v of yticks) g.append(L(x0, ys(v), x1, ys(v), { stroke: GRID, width: 0.4 }));
  g.append(L(x0, y1, x1, y1), L(x0, y0, x0, y1));
  for (const v of yticks) { const y = ys(v); g.append(L(x0 - 3, y, x0, y), T(x0 - 4.5, y + fs * 0.35, yfmt(v), fs, { anchor: 'end' })); }
  if (xTickLabels) for (const [v, lab] of xTickLabels) { const x = xs(v); g.append(L(x, y1, x, y1 + 3), T(x, y1 + 3 + fs * 1.05, lab, fs, { anchor: 'middle' })); }
  else for (const v of xticks) { const x = xs(v); g.append(L(x, y1, x, y1 + 3), T(x, y1 + 3 + fs * 1.05, xfmt(v), fs, { anchor: 'middle' })); }
  if (xlabel) g.append(T((x0 + x1) / 2, y1 + 3 + fs * 2.45, xlabel, fs, { anchor: 'middle' }));
  if (ylabel) { const yw = Math.max(...yticks.map(v => textW(yfmt(v), fs))); const xx = x0 - 7 - yw - fs * 0.6; g.append(T(xx, (y0 + y1) / 2, ylabel, fs, { anchor: 'middle', rotate: -90 })); }
}
const ticksFor = (lo, hi, n = 5) => niceTicks(lo, hi, n).filter(v => v >= lo - 1e-9 && v <= hi + 1e-9);

// ------------------------------------------------------------------ plots
export function drawPlot(app, o) {
  const ens = app.ens;
  const W = o.w * MM, H = o.h * MM, fs = o.font;
  const root = svg('svg', { viewBox: `0 0 ${W.toFixed(2)} ${H.toFixed(2)}`, class: 'plot-svg', 'font-family': FIG_FONT });
  const g = svg('g');
  root.append(svg('rect', { x: 0, y: 0, width: W, height: H, fill: '#ffffff' }), g);
  const title = o.title?.trim();
  let top = 6 + (title ? fs * 1.6 : 0);
  if (title) g.append(T(6, 6 + fs * 1.1, title, fs * 1.1, { bold: true }));
  const key = o.metric || app.state.metric, m = ens.metricByKey[key] || { label: key, unit: '' };
  const mlab = `${m.label}${m.unit ? ` (${m.unit})` : ''}`;
  const order = o.order === 'input' ? [...Array(ens.n).keys()] : app.order;
  const pal = SLOTS;
  const hl = o.highlight ? app.state.conf : -1;
  const ringed = (x, y, r, fill) => svg('circle', { cx: x, cy: y, r, fill, stroke: '#ffffff', 'stroke-width': 0.6 });

  if (o.type === 'matrix') {
    const n = ens.n;
    let idx = order.slice();
    let tree = null;
    if (o.cluster && n >= 3 && n <= 300) { tree = upgma(ens.rmsd); idx = tree.order; }
    const ramp = RAMPS[o.ramp] || BLUE_RAMP;
    const max = Math.max(1e-6, ...ens.rmsd.flat());
    const labels = n <= 60;
    const cbW = 34;
    const dendH = tree ? Math.min(40, H * 0.16) : 0;
    // Row/column labels never overlap: the label size follows the cell size; when full names would get too small,
    // the common prefix ("Model ", "conf_0") moves into an axis title and only the numbers are written; when even
    // those are too small (<5 pt), every k-th label is written.
    const ids = idx.map(k => ens.shortId(k));
    let prefix = ids.reduce((a, b) => { let q = 0; while (q < a.length && a[q] === b[q]) q++; return a.slice(0, q); });
    prefix = prefix.replace(/\d+$/, '');
    const shortIds = ids.map(t => t.slice(prefix.length).replace(/^0+(?=\d)/, ''));
    const canShort = prefix.length >= 2 && shortIds.every(t => t.length);
    const maxLfs = fs * 0.85, minLfs = Math.min(5, maxLfs);
    const layout = (texts, lfs, titled) => {
      const tw = labels ? Math.max(...texts.map(t => textW(t, lfs))) : 0;
      const labW = labels ? tw + 5 + (titled ? fs * 1.4 : 0) : fs * 2;
      const x0 = 6 + labW, y0 = top + dendH + 2;
      const size = Math.min(W - x0 - cbW - 10, H - y0 - labW - 6);
      return { tw, labW, x0, y0, size, cell: size / n, lfs, texts, titled };
    };
    const fit = (texts, titled) => {
      let lfs = maxLfs, lay = layout(texts, lfs, titled);
      for (let it = 0; it < 4; it++) { lfs = Math.max(minLfs, Math.min(maxLfs, lay.cell * 0.9)); lay = layout(texts, lfs, titled); }
      return lay;
    };
    let lay = fit(ids, false);
    if (labels && canShort && lay.cell * 0.9 < Math.min(maxLfs, 6.5)) lay = fit(shortIds, true);
    const { x0, y0, size, cell, lfs } = lay;
    const every = labels ? Math.max(1, Math.ceil(lfs / (cell * 0.9) - 1e-9)) : 1;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      const v = ens.rmsd[idx[i]][idx[j]];
      g.append(svg('rect', { x: x0 + j * cell, y: y0 + i * cell, width: cell + 0.05, height: cell + 0.05, fill: sampleRamp(ramp, v / max) }));
      if (o.values && n <= 12) {
        const c = sampleRamp(ramp, v / max), lum = parseInt(c.slice(1, 3), 16) * 0.299 + parseInt(c.slice(3, 5), 16) * 0.587 + parseInt(c.slice(5, 7), 16) * 0.114;
        g.append(T(x0 + (j + 0.5) * cell, y0 + (i + 0.5) * cell + fs * 0.3, fmt(v, 1), Math.min(fs * 0.85, cell * 0.32), { anchor: 'middle', fill: lum > 150 ? INK : '#ffffff' }));
      }
    }
    g.append(svg('rect', { x: x0, y: y0, width: size, height: size, fill: 'none', stroke: AXIS, 'stroke-width': 0.6 }));
    if (hl >= 0) { const p = idx.indexOf(hl); if (p >= 0) g.append(svg('rect', { x: x0 + p * cell, y: y0 + p * cell, width: cell, height: cell, fill: 'none', stroke: '#eb6834', 'stroke-width': 1 })); }
    if (labels) {
      lay.texts.forEach((t, i) => {
        if (i % every) return;
        g.append(T(x0 - 3, y0 + (i + 0.5) * cell + lfs * 0.35, t, lfs, { anchor: 'end' }));
        const xx = x0 + (i + 0.5) * cell, yy = y0 + size + 3;
        g.append(T(xx + lfs * 0.35, yy, t, lfs, { anchor: 'end', rotate: -90 }));
      });
      if (lay.titled) {
        const name = /^conf/i.test(prefix) ? 'Conformation' : prefix.replace(/[\s_:#-]+$/, '');
        const ty = x0 - 3 - lay.tw - 4, cy = y0 + size / 2;
        g.append(T(ty, cy, name, fs * 0.9, { anchor: 'middle', rotate: -90 }));
        g.append(T(x0 + size / 2, y0 + size + 3 + lay.tw + 4 + fs * 0.8, name, fs * 0.9, { anchor: 'middle' }));
      }
    }
    if (tree) {
      // dendrogram above the matrix (heights scaled to dendH)
      const pos = new Map(idx.map((k, i) => [k, x0 + (i + 0.5) * cell]));
      const hmax = tree.root.h || 1;
      const yOf = h => y0 - 2 - (h / hmax) * (dendH - 2);
      const place = c => {
        if (!c.left) return { x: pos.get(c.id), y: y0 - 2 };
        const a = place(c.left), b = place(c.right), y = yOf(c.h);
        g.append(svg('path', { d: `M ${a.x} ${a.y} L ${a.x} ${y} L ${b.x} ${y} L ${b.x} ${b.y}`, fill: 'none', stroke: SUB, 'stroke-width': 0.6 }));
        return { x: (a.x + b.x) / 2, y };
      };
      place(tree.root);
    }
    // colour bar
    const cbx = x0 + size + 10, cbh = size * 0.7, cby = y0 + (size - cbh) / 2, steps = 40;
    for (let s = 0; s < steps; s++) g.append(svg('rect', { x: cbx, y: cby + cbh * (1 - (s + 1) / steps), width: 8, height: cbh / steps + 0.05, fill: sampleRamp(ramp, (s + 0.5) / steps) }));
    g.append(svg('rect', { x: cbx, y: cby, width: 8, height: cbh, fill: 'none', stroke: AXIS, 'stroke-width': 0.5 }));
    for (const v of ticksFor(0, max, 4)) { const y = cby + cbh * (1 - v / max); g.append(L(cbx + 8, y, cbx + 10.5, y, { width: 0.5 }), T(cbx + 12, y + fs * 0.32, fmt(v, 1), fs * 0.85)); }
    g.append(T(cbx, cby - fs * 0.7, 'RMSD (Å)', fs * 0.85));
    return root;
  }

  if (o.type === 'ss2d') {
    const k = app.state.conf;
    const target = svgBuilder(g);
    const colors = app.colorsFor(k, app.state.prev).hex;
    drawStructure2D(target, ens, orientForBox(ens, layoutFor(ens, k), W - 12, H - top - 6), { k, prev: o.highlight ? app.state.prev : null, colors, box: [6, top, W - 12, H - top - 6], dark: false, unit: fs / 9 });
    return root;
  }

  // xy plots share the frame
  const yTickProbe = v => textW(fmt(v, 1), fs);
  const left = 6 + fs * 1.6 + 30, right = W - 8, bottom = H - fs * 3.4 - (o.type === 'rmsf' && o.track ? 8 : 0);
  if (o.type === 'toref' || o.type === 'profile') {
    let vals, ylabel, refName = '';
    if (o.type === 'toref') {
      const atoms = c4Atoms(ens);
      let ref;
      if (o.ref === 'mean') { ref = meanCoords(ens); refName = 'ensemble mean'; }
      else {
        const k = o.ref === 'current' ? app.state.conf : o.ref === 'previous' ? (app.state.prev ?? order[0]) : o.ref === 'best'
          ? [...Array(ens.n).keys()].sort((a, b) => (ens.value(a, key) ?? Infinity) - (ens.value(b, key) ?? Infinity))[0] : order[0];
        ref = ens.coords[k]; refName = ens.shortId(k);
      }
      vals = order.map(k => rmsdFit(ens, ref, ens.coords[k], atoms));
      ylabel = o.ylabel || `C4′ RMSD to ${refName} (Å)`;
    } else { vals = order.map(k => ens.value(k, key)); ylabel = o.ylabel || mlab; }
    const fin = vals.filter(v => v !== null && Number.isFinite(v));
    let lo = o.type === 'toref' ? 0 : Math.min(...fin), hi = Math.max(...fin);
    if (o.type !== 'toref') { const p = (hi - lo) * 0.08 || 1; lo -= p; hi += p; } else hi = hi * 1.08 || 1;
    const yt = ticksFor(lo, hi, 5);
    const x0 = left + Math.max(0, Math.max(...yt.map(yTickProbe)) - 30), x1 = right, y0 = top + 4, y1 = bottom;
    const n = order.length, band = (x1 - x0) / n;
    const xs = i => x0 + (i + 0.5) * band, ys = v => y1 - (v - lo) / (hi - lo) * (y1 - y0);
    const labelEvery = Math.max(1, Math.ceil(n / Math.floor((x1 - x0) / (fs * 2.2))));
    const xl = order.map((k, i) => [i, n <= 14 ? ens.shortId(k).replace(/^conf_0*/, '') : String(i + 1)]).filter(([i]) => i % labelEvery === 0);
    axes(g, { x0, x1, y0, y1, xs, ys, xticks: [], yticks: yt, xTickLabels: xl, xlabel: o.xlabel || (n <= 14 ? 'Conformation' : 'Conformation (order)'), ylabel, fs, grid: o.grid });
    if (o.style === 'bars' || o.type === 'toref') {
      const bw = Math.min(band * 0.7, 14);
      const base = lo <= 0 && hi >= 0 ? ys(0) : y1;            // bars grow from zero (or the axis)
      vals.forEach((v, i) => {
        if (v === null) return;
        const yy = ys(v);
        g.append(svg('rect', { x: xs(i) - bw / 2, y: Math.min(yy, base), width: bw, height: Math.abs(base - yy), fill: order[i] === hl ? '#eb6834' : pal[0] }));
      });
    } else {
      const pts = vals.map((v, i) => (v === null ? null : [xs(i), ys(v)])).filter(Boolean);
      g.append(svg('polyline', { points: pts.map(p => p.join(',')).join(' '), fill: 'none', stroke: pal[0], 'stroke-width': 1.2, 'stroke-linejoin': 'round' }));
      vals.forEach((v, i) => { if (v !== null) g.append(ringed(xs(i), ys(v), order[i] === hl ? 3.2 : 2.2, order[i] === hl ? '#eb6834' : pal[0])); });
    }
    return root;
  }

  if (o.type === 'rmsf') {
    const n = ens.nRes;
    const series = [{ name: 'RMSF (ensemble)', vals: ens.rmsf, color: pal[0] }];
    if (o.moved && app.state.prev !== null && app.state.prev !== app.state.conf) {
      series.push({ name: `Moved since ${ens.shortId(app.state.prev)}`, vals: ens.residueDisplacement(ens.coords[app.state.conf], ens.coords[app.state.prev]), color: pal[1] });
    }
    const hi = Math.max(1e-6, ...series.flatMap(s => s.vals)) * 1.08;
    const yt = ticksFor(0, hi, 5);
    const x0 = left + Math.max(0, Math.max(...yt.map(yTickProbe)) - 30), x1 = right;
    const leg = series.length > 1 ? legendLayout(series.map(s => ({ name: s.name, color: s.color, kind: 'line' })), x0 + 4, x1, fs) : null;
    const y0 = top + 4 + (leg ? leg.height : 0), y1 = bottom;
    const xs = i => x0 + (n > 1 ? i / (n - 1) : 0.5) * (x1 - x0), ys = v => y1 - v / hi * (y1 - y0);
    const step = n > 120 ? 20 : 10;
    const xl = [[0, ens.residues[0].num]];
    for (let i = step - 1; i < n; i += step) xl.push([i, ens.residues[i].num]);
    axes(g, { x0, x1, y0, y1, xs, ys, xticks: [], yticks: yt, xTickLabels: xl, xlabel: o.xlabel || 'Residue', ylabel: o.ylabel || (series.length > 1 ? 'RMSF, displacement (Å)' : 'RMSF (Å)'), fs, grid: o.grid });
    for (const s of series) {
      const pts = s.vals.map((v, i) => [xs(i), ys(v)]);
      if (s === series[0]) g.append(svg('path', { d: `M ${x0} ${y1} ` + pts.map(p => `L ${p[0]} ${p[1]}`).join(' ') + ` L ${x1} ${y1} Z`, fill: s.color, 'fill-opacity': 0.12 }));
      g.append(svg('polyline', { points: pts.map(p => p.join(',')).join(' '), fill: 'none', stroke: s.color, 'stroke-width': 1.2, 'stroke-linejoin': 'round' }));
    }
    if (leg) drawLegend(g, leg, top);
    if (o.track) {
      // pairing track of the shown conformation (paired = filled) and the pairing frequency across the ensemble
      const f = pairedFraction(ens), part = ens.partner(app.state.conf), ty = y1 + fs * 3.25;
      for (let i = 0; i < n; i++) {
        const xa = xs(i) - (x1 - x0) / Math.max(1, n - 1) / 2, wcell = (x1 - x0) / Math.max(1, n - 1);
        g.append(svg('rect', { x: xa, y: ty, width: wcell + 0.05, height: 4, fill: sampleRamp(['#f2f2f2', ...BLUE_RAMP.slice(3)], f[i]) }));
        if (part[i] >= 0) g.append(svg('rect', { x: xa, y: ty + 5, width: wcell + 0.05, height: 2, fill: INK }));
      }
      g.append(T(x0 - (x1 - x0) / Math.max(1, n - 1) / 2 - 3, ty + 5, 'paired', fs * 0.75, { anchor: 'end', fill: SUB }));
    }
    return root;
  }

  if (o.type === 'funnel') {
    const best = [...Array(ens.n).keys()].sort((a, b) => (ens.value(a, key) ?? Infinity) - (ens.value(b, key) ?? Infinity))[0];
    const xsV = ens.rmsd[best].slice(), ysV = [...Array(ens.n).keys()].map(k => ens.value(k, key));
    const states = [...new Set(ens.conformations.map(c => c.dbn))];
    const fin = ysV.filter(v => v !== null);
    let lo = Math.min(...fin), hi = Math.max(...fin); const p = (hi - lo) * 0.08 || 1; lo -= p; hi += p;
    const xhi = Math.max(1, ...xsV) * 1.08;
    const yt = ticksFor(lo, hi, 5), xt = ticksFor(0, xhi, 5);
    const x0 = left + Math.max(0, Math.max(...yt.map(yTickProbe)) - 30), x1 = right;
    const leg = states.length > 1 && states.length <= 8 ? legendLayout(states.map((_, i) => ({ name: `2D state ${i + 1}`, color: pal[i % pal.length], kind: 'dot' })), x0 + 4, x1, fs) : null;
    const y0 = top + 4 + (leg ? leg.height : 0), y1 = bottom;
    const xs = v => x0 + v / xhi * (x1 - x0), ys = v => y1 - (v - lo) / (hi - lo) * (y1 - y0);
    axes(g, { x0, x1, y0, y1, xs, ys, xticks: xt, yticks: yt, xlabel: o.xlabel || `C4′ RMSD to ${ens.shortId(best)} (Å)`, ylabel: o.ylabel || mlab, fs, grid: o.grid });
    for (let k = 0; k < ens.n; k++) {
      if (ysV[k] === null) continue;
      const si = states.indexOf(ens.conformations[k].dbn);
      g.append(ringed(xs(xsV[k]), ys(ysV[k]), k === hl ? 3.4 : 2.4, states.length > 1 ? pal[si % pal.length] : pal[0]));
      if (k === hl) g.append(svg('circle', { cx: xs(xsV[k]), cy: ys(ysV[k]), r: 5, fill: 'none', stroke: INK, 'stroke-width': 0.7 }));
    }
    if (leg) drawLegend(g, leg, top);
    return root;
  }
  return root;
}

// SVG builder with explicit attributes (for the 2D-structure plot).
function svgBuilder(root) {
  const st = (o) => ({ stroke: o.stroke || 'none', 'stroke-width': o.width || 1, 'stroke-dasharray': o.dash ? o.dash.join(' ') : null, 'stroke-linecap': 'round', 'stroke-linejoin': 'round', opacity: o.alpha ?? null });
  return {
    line(x1, y1, x2, y2, o) { root.append(svg('line', { x1, y1, x2, y2, ...st(o) })); },
    path(cmds, o) { root.append(svg('path', { d: cmds.map(c => c[0] + ' ' + c.slice(1).map(v => +v.toFixed(2)).join(' ')).join(' '), fill: o.fill || 'none', ...st(o) })); },
    circle(cx, cy, r, o) { root.append(svg('circle', { cx, cy, r, fill: o.fill || 'none', stroke: o.stroke || 'none', 'stroke-width': o.width || 1 })); },
    text(x, y, s, o) { root.append(svg('text', { x, y, 'font-size': o.size, 'font-family': FIG_FONT, 'font-weight': (o.weight || 400) >= 600 ? 'bold' : null, fill: o.fill, 'text-anchor': o.anchor || 'start', text: s })); },
  };
}

// ------------------------------------------------------------------ dialog
const SIZES = [['89x70', 'Single column 89 × 70 mm'], ['89x89', 'Single column 89 × 89 mm'], ['120x85', '1.5 columns 120 × 85 mm'], ['183x90', 'Double column 183 × 90 mm'], ['183x120', 'Double column 183 × 120 mm'], ['custom', 'Custom']];

export function plotStudio(app, initial = null) {
  const ens = app.ens;
  if (!ens) return;
  const prefs = app.prefs.plots || {};
  const o = { type: initial || prefs.type || 'matrix', size: prefs.size || '89x70', w: prefs.w || 89, h: prefs.h || 70, font: prefs.font || 8, title: '', xlabel: '', ylabel: '',
    ramp: prefs.ramp || 'blues', cluster: prefs.cluster ?? true, values: prefs.values ?? true, grid: prefs.grid ?? false, highlight: prefs.highlight ?? true,
    ref: prefs.ref || 'first', style: prefs.style || 'line', moved: prefs.moved ?? true, track: prefs.track ?? true, metric: app.state.metric, order: 'current' };
  const list = el('div', { class: 'ps-list' });
  const preview = el('div', { class: 'ps-preview' });
  const opts = el('div', { class: 'ps-opts' });
  let current = null;
  const redraw = () => {
    const [w, h] = o.size === 'custom' ? [+o.w, +o.h] : o.size.split('x').map(Number);
    o.w = w; o.h = h; o.font = +o.font || 8;
    current = drawPlot(app, o);
    preview.replaceChildren(current);
    for (const b of list.children) b.classList.toggle('on', b.dataset.type === o.type);
    app.prefs.plots = { type: o.type, size: o.size, w: o.w, h: o.h, font: o.font, ramp: o.ramp, cluster: o.cluster, values: o.values, grid: o.grid, highlight: o.highlight, ref: o.ref, style: o.style, moved: o.moved, track: o.track };
    app.savePrefs();
  };
  for (const [id, label] of PLOTS) list.append(el('button', { class: 'ps-item', 'data-type': id, onclick: () => { o.type = id; o.title = ''; o.xlabel = ''; o.ylabel = ''; buildOpts(); redraw(); } }, label));
  const inp = (key, attrs = {}) => { const i = el('input', { class: 'input', value: o[key] ?? '', ...attrs }); i.addEventListener('input', () => { o[key] = attrs.type === 'number' ? +i.value : i.value; redraw(); }); return i; };
  const chk = (key, label) => { const c = el('input', { type: 'checkbox', checked: o[key] ? true : null }); c.addEventListener('change', () => { o[key] = c.checked; redraw(); }); return el('label', { class: 'check' }, c, label); };
  const sel = (key, options) => { const s = select(options.map(([v, l]) => ({ value: v, label: l })), String(o[key])); s.addEventListener('change', () => { o[key] = s.value; buildOpts(); redraw(); }); return s; };
  function buildOpts() {
    const t = o.type;
    const metricOpts = ens.metrics.filter(m => ens.conformations.some(c => typeof c.values[m.key] === 'number')).map(m => [m.key, m.label]);
    opts.replaceChildren(...[
      field('Size', sel('size', SIZES)),
      o.size === 'custom' ? el('div', { class: 'grid2' }, field('Width (mm)', inp('w', { type: 'number', min: 30, max: 400 })), field('Height (mm)', inp('h', { type: 'number', min: 20, max: 400 }))) : null,
      field('Font size (pt)', sel('font', [['6', '6 pt'], ['7', '7 pt'], ['8', '8 pt'], ['9', '9 pt'], ['10', '10 pt'], ['12', '12 pt']])),
      field('Title (optional)', inp('title', { placeholder: 'none' })),
      t !== 'matrix' && t !== 'ss2d' ? el('div', { class: 'grid2' }, field('x-axis label', inp('xlabel', { placeholder: 'automatic' })), field('y-axis label', inp('ylabel', { placeholder: 'automatic' }))) : null,
      t === 'matrix' ? field('Colour ramp', sel('ramp', [['blues', 'Blues'], ['grays', 'Grays (print)'], ['purples', 'Purples'], ['greens', 'Greens'], ['oranges', 'Oranges']])) : null,
      t === 'matrix' ? chk('cluster', 'Cluster (average linkage) with dendrogram') : null,
      t === 'matrix' ? chk('values', 'Values in cells (up to 12 conformations)') : null,
      t === 'toref' ? field('Reference', sel('ref', [['first', 'First in the current order'], ['best', 'Best by the energy metric'], ['current', 'Shown conformation'], ['previous', 'Previous conformation'], ['mean', 'Ensemble mean structure']])) : null,
      t === 'profile' ? field('Style', sel('style', [['line', 'Line with markers'], ['bars', 'Bars']])) : null,
      t === 'profile' || t === 'funnel' ? field('Energy / metric', sel('metric', metricOpts)) : null,
      t === 'rmsf' ? chk('moved', 'Also: displacement since the previous conformation') : null,
      t === 'rmsf' ? chk('track', 'Pairing track under the axis') : null,
      t !== 'matrix' && t !== 'ss2d' ? chk('grid', 'Light gridlines') : null,
      t !== 'ss2d' ? chk('highlight', 'Highlight the shown conformation') : chk('highlight', 'Mark pairs new since the previous conformation'),
      el('div', { class: 'desc', text: 'Text uses Liberation Sans (open licence, metric-compatible with Helvetica/Arial); PDF embeds it. Colours follow the active palette (Color ▸ Palettes).' })].filter(Boolean));
  }
  buildOpts();
  const dpiSel = select([{ value: '300', label: '300 dpi' }, { value: '600', label: '600 dpi' }, { value: '1200', label: '1200 dpi' }], '600');
  const name = () => `${slug(ens.name)}_${o.type}`;
  const run = fmtId => async btn => {
    btn.disabled = true;
    try { await ensureFigureFont(); const saved = await exportChart(current, fmtId, name(), { dpi: +dpiSel.value, title: o.title || PLOTS.find(p => p[0] === o.type)[1] }); if (saved) toast(`${saved} · ${o.w} × ${o.h} mm`, 'ok', 'Plot exported'); }
    catch (err) { toast(err.message, 'error', 'Export failed'); } finally { btn.disabled = false; }
  };
  dialog({
    title: 'Plot studio', width: 1080,
    body: b => b.append(el('div', { class: 'ps-grid' }, list, el('div', { class: 'ps-center' }, preview, el('div', { class: 'desc', text: 'Preview at the chosen physical size (scaled to fit). Right-click a chart in the side panel for quick exports too.' })), opts)),
    actions: [{ node: el('span', { class: 'muted', text: 'Resolution for PNG/TIFF' }) }, { node: dpiSel }, { spacer: true },
      { label: 'Copy image', action: async () => { try { await copyChart(current, +dpiSel.value); toast('Plot copied.', 'ok'); } catch { toast('Copying images is not available here.', 'error'); } } },
      { label: 'PNG', action: run('png') }, { label: 'TIFF', action: run('tiff') }, { label: 'SVG', action: run('svg') }, { label: 'PDF', primary: true, action: run('pdf') }],
  });
  ensureFigureFont().then(redraw);
  redraw();
}
