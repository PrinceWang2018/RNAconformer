// RCV Lite: the interactive page written by "Export ▸ Interactive HTML". One self-contained file (no server, no
// install) for supplementary material, reviewers and collaborators: rotate, morph between conformations, play the
// tour, change representation and colours, read energies and the 2D structure of every conformation.
import { Viewer } from './viewer/viewer.js';
import { Ensemble } from './ensemble.js';
import { residueColors, legend as schemeLegend } from './colors.js';
import { layoutFor, blendLayouts, drawStructure2D, canvas2DTarget } from './ss2dlayout.js';
import { el, svg, easeInOut, fmt, clamp } from './util.js';
import { drawMap } from './landscape.js';
import { buildMechanism } from './mechmorph.js';

const data = JSON.parse(document.getElementById('rcv-data').textContent);
const doc = data.doc;
const ens = new Ensemble(doc);
const view = data.view || {};
const $ = s => document.querySelector(s);
const hexToLin = hex => {
  const n = parseInt(hex.slice(1), 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255].map(v => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); });
};
const SCHEMES = [['nucleotide', 'Nucleotide'], ['helix', 'Helices'], ['elements', 'Structure elements'], ['pairing', 'Base pairing'], ['pairfreq', 'Pairing frequency'], ['switching', 'Switching residues'], ['flexibility', 'Flexibility (RMSF)'], ['position', 'Position 5′→3′']];
const REPS = [['cartoon', 'Cartoon'], ['cartoon+sticks', 'Cartoon + atoms'], ['sticks', 'Sticks'], ['ballstick', 'Ball & stick'], ['spacefill', 'Spheres']];

const state = { conf: view.conf ?? 0, prev: null, scheme: SCHEMES.some(s => s[0] === view.colorScheme) ? view.colorScheme : 'nucleotide', playing: false, morph: null, metric: doc.metrics.some(m => m.key === view.metric) ? view.metric : doc.default_metric };
const order = [...Array(ens.n).keys()].sort((a, b) => (ens.conformations[a].rank ?? a) - (ens.conformations[b].rank ?? b));
if (state.conf >= ens.n) state.conf = order[0];

document.title = `${data.title || ens.name} · RCV`;
$('#title').textContent = data.title || ens.name;
$('#subtitle').textContent = `${ens.nRes} nt · ${ens.n} conformation${ens.n === 1 ? '' : 's'} · ${doc.origin === 'rnafold3d-engine' ? 'predicted with the rnafold3d engine' : 'experimental or imported structures'}`;
$('#sequence').textContent = ens.sequence;

const viewer = new Viewer({ canvas: $('#gl'), container: $('#vp'), overlay: $('#ov'), labels: $('#lb'), onPick: () => {}, onHover: (p, e) => hover(p, e), onDblPick: p => { if (p) viewer.centerOn(viewer.labelAnchor(p.residue)); } });
viewer.settings.background = view.background || '#ffffff';
viewer.settings.style = view.style && viewer.styles[view.style] ? view.style : 'studio';
viewer.labelMode = 'none';
const display = ens.coords[state.conf].slice();
viewer.setEnsemble(ens, display);
viewer.setRepresentation(REPS.some(r => r[0] === view.representation) ? view.representation : 'cartoon');
viewer.setStyle(viewer.settings.style);
if (view.camera) viewer.controls.state = view.camera;
$('#vp').style.background = viewer.settings.background;

let colorsLin = null;
const colorsFor = k => {
  const hex = residueColors(ens, state.scheme, { conf: k, coords: ens.coords[k], prevCoords: null, uniform: '#5b7a99' });
  const lin = new Float32Array(ens.nRes * 3);
  hex.forEach((h, i) => { const c = hexToLin(h); lin.set(c, 3 * i); });
  return { hex, lin };
};
function applyColors() {
  const c = colorsFor(state.conf);
  colorsLin = c.lin;
  viewer.setResidueColorsLinear(colorsLin);
  const lg = schemeLegend(ens, state.scheme, { conf: state.conf });
  const box = $('#legend');
  box.replaceChildren(el('div', { class: 'lg-title', text: lg?.title || '' }));
  if (lg?.ramp) box.append(el('span', { text: lg.min }), el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${lg.ramp.join(',')})` } }), el('span', { text: lg.max }));
  else for (const it of lg?.items || []) box.append(el('span', { class: 'lg' }, el('span', { class: 'sw', style: { background: it.color } }), it.label));
  draw2D();
}

// ---------------------------------------------------------------- morph
function goTo(k, dur = 1.1) {
  if (k === state.conf && !state.morph) return;
  const from = state.conf;
  state.prev = from;
  state.conf = k;
  const start = display.slice(), target = ens.coords[k];
  const c0 = colorsLin.slice(), c1 = colorsFor(k).lin;
  const L0 = shown2D ? shown2D.slice() : layoutFor(ens, from), L1 = layoutFor(ens, k);
  const t0 = performance.now();
  state.morph = true;
  viewer._addAnimator(() => {
    const t = clamp((performance.now() - t0) / (dur * 1000), 0, 1), e = easeInOut(t);
    for (let i = 0; i < display.length; i++) display[i] = start[i] + (target[i] - start[i]) * e;
    viewer.setCoords(display);
    const cc = new Float32Array(c0.length);
    for (let i = 0; i < cc.length; i++) cc[i] = c0[i] + (c1[i] - c0[i]) * e;
    colorsLin = cc;
    viewer.setResidueColorsLinear(cc);
    shown2D = blendLayouts(L0, L1, e);
    draw2D({ from, t: e });
    if (t >= 1) { state.morph = null; applyColors(); renderList(); renderMap(); if (state.playing) playTimer = setTimeout(step, 900); }
    return t < 1;
  });
  renderList();
}
let playTimer = 0;
function step(dir = 1) {
  const i = order.indexOf(state.conf);
  goTo(order[(i + dir + order.length) % order.length]);
}
function togglePlay() {
  state.playing = !state.playing;
  $('#play').textContent = state.playing ? 'Pause' : 'Play';
  clearTimeout(playTimer);
  if (state.playing) step();
}

// ---------------------------------------------------------------- 2D panel
let shown2D = null;
function draw2D(m = null) {
  const cv = $('#ss'), dpr = Math.min(2, window.devicePixelRatio || 1);
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!w || !h) return;
  cv.width = w * dpr; cv.height = h * dpr;
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, h);
  const coords = shown2D && !m ? layoutFor(ens, state.conf) : (shown2D || layoutFor(ens, state.conf));
  const dark = matchMedia('(prefers-color-scheme: dark)').matches;
  drawStructure2D(canvas2DTarget(ctx), ens, coords, { k: state.conf, prev: state.prev, from: m ? m.from : null, t: m ? m.t : 1, colors: colorsFor(m && m.t < 0.5 ? m.from : state.conf).hex, box: [6, 6, w - 12, h - 12], dark });
  $('#ss-cap').textContent = `2D of ${ens.shortId(state.conf)} · ${ens.conformations[state.conf].pairs.length} base pairs`;
}

// ---------------------------------------------------------------- lists & charts
const metricSel = $('#metric');
for (const m of doc.metrics.filter(m => ens.conformations.some(c => typeof c.values[m.key] === 'number'))) metricSel.append(el('option', { value: m.key, text: `${m.label}${m.unit ? ' (' + m.unit + ')' : ''}` }));
metricSel.value = state.metric;
metricSel.addEventListener('change', () => { state.metric = metricSel.value; renderList(); });
function renderList() {
  const m = doc.metrics.find(x => x.key === state.metric) || {};
  const list = $('#confs');
  list.replaceChildren(...order.map((k, i) => {
    const c = ens.conformations[k], v = ens.value(k, state.metric);
    const row = el('button', { class: `conf${k === state.conf ? ' on' : ''}`, onclick: () => { state.playing && togglePlay(); goTo(k); } },
      el('span', { class: 'n', text: String(i + 1) }), el('b', { text: ens.shortId(k) }), el('span', { class: 'v', text: v === null ? '' : fmt(v) }),
      el('span', { class: 'dbn', text: c.dbn }));
    return row;
  }));
  // energy profile
  const s = $('#profile'), W = s.clientWidth || 280, H = 90;
  const vals = order.map(k => ens.value(k, state.metric));
  const fin = vals.filter(v => v !== null);
  s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  s.replaceChildren();
  if (fin.length < 2) return;
  let lo = Math.min(...fin), hi = Math.max(...fin);
  if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
  const X = i => 10 + i * (W - 20) / Math.max(1, order.length - 1), Y = v => H - 12 - (v - lo) / (hi - lo) * (H - 24);
  const ns = 'http://www.w3.org/2000/svg';
  const mk = (tag, a) => { const n = document.createElementNS(ns, tag); for (const [k, v] of Object.entries(a)) n.setAttribute(k, v); s.append(n); return n; };
  mk('polyline', { points: vals.map((v, i) => (v === null ? null : `${X(i)},${Y(v)}`)).filter(Boolean).join(' '), fill: 'none', stroke: 'var(--accent)', 'stroke-width': 2 });
  vals.forEach((v, i) => { if (v !== null) mk('circle', { cx: X(i), cy: Y(v), r: order[i] === state.conf ? 5.5 : 3, fill: order[i] === state.conf ? '#eb6834' : 'var(--accent)' }); });
  $('#profile-cap').textContent = `${m.label || state.metric}${m.unit ? ' (' + m.unit + ')' : ''} per conformation (rank order)`;
}

function hover(p, e) {
  const tip = $('#tip');
  if (!p || !e) { tip.hidden = true; return; }
  const r = p.residue, pi = ens.pairInfo(state.conf, r);
  tip.textContent = `${ens.resLabel(r)} · ${pi ? 'paired with ' + ens.resLabel(pi.partner) : 'unpaired'} · RMSF ${fmt(ens.rmsf[r], 1)} Å`;
  tip.style.left = `${e.clientX + 12}px`; tip.style.top = `${e.clientY + 12}px`;
  tip.hidden = false;
}

// ---------------------------------------------------------------- controls
const repSel = $('#rep'), colSel = $('#scheme');
for (const [v, l] of REPS) repSel.append(el('option', { value: v, text: l }));
for (const [v, l] of SCHEMES) colSel.append(el('option', { value: v, text: l }));
repSel.value = viewer.settings.representation;
colSel.value = state.scheme;
repSel.addEventListener('change', () => viewer.setRepresentation(repSel.value));
colSel.addEventListener('change', () => { state.scheme = colSel.value; applyColors(); });
$('#prev').addEventListener('click', () => step(-1));
$('#next').addEventListener('click', () => step(1));
$('#play').addEventListener('click', togglePlay);
$('#overlay').addEventListener('change', e => { viewer.fitAll = e.target.checked; viewer.setGhosts(e.target.checked ? ens.coords.filter((_, k) => k !== state.conf) : null); });
$('#reset').addEventListener('click', () => viewer.resetView());
window.addEventListener('keydown', e => {
  if (e.target.closest('select, input')) return;
  if (e.key === 'ArrowRight') step(1); else if (e.key === 'ArrowLeft') step(-1); else if (e.key === ' ') { togglePlay(); e.preventDefault(); }
});
new ResizeObserver(() => draw2D()).observe($('#ss'));
if (data.methods) $('#methods').textContent = data.methods;
else $('#methods-box').hidden = true;
$('#made').textContent = `Made with RCV ${data.version || ''} — RNA Conformation Viewer by Zixiang Wang (open source, GPL-3.0)`;
// ---------------------------------------------------------------- RCV 0.7: ensemble map and pathway
const analysis = data.analysis || {};
const fmtTime = x => {
  if (x === null || x === undefined || !Number.isFinite(x)) return '–';
  if (x < 1e-3) return `${fmt(x * 1e6, 1)} µs`; if (x < 1) return `${fmt(x * 1e3, 1)} ms`; if (x < 120) return `${fmt(x, 1)} s`;
  if (x < 7200) return `${fmt(x / 60, 1)} min`; if (x < 172800) return `${fmt(x / 3600, 1)} h`; if (x < 3.15e7) return `${fmt(x / 86400, 1)} days`;
  return `${(x / 3.15e7).toExponential(1)} years`;
};
function mapRoutes(r) {
  if (!r) return [];
  const conf = k => (k === 0 ? r.a : k === 1 ? r.b : r.nodes[k]?.conformations?.[0] ?? null);
  const merged = new Map();
  for (const p of (r.same_basin ? [] : r.tpt?.pathways?.length ? r.tpt.pathways : [{ nodes: r.best.route, flux_share: 1 }])) {
    const confs = p.nodes.map(conf).filter(k => k !== null), key = confs.join('>');
    const m = merged.get(key) || { confs, share: 0 };
    m.share += p.flux_share;
    merged.set(key, m);
  }
  return [...merged.values()].sort((a, b) => b.share - a.share).slice(0, 3)
    .map((m, i) => ({ nodes: m.confs, width: 1.2 + 4 * Math.min(1, m.share), color: ['var(--series-2)', 'var(--series-1)', '#1baf7a'][i] }));
}
let mapPos = null;
function renderMap() {
  const L = doc.ensemble?.landscape;
  if (!L) { $('#map-box').hidden = true; return; }
  const s = $('#map'), W = s.clientWidth || 300, H = 220;
  s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  s.replaceChildren();
  const g = svg('g');
  const r = analysis.pathway;
  mapPos = drawMap(g, ens, L, { x0: 4, y0: 2, x1: W - 4, y1: H - 4 }, { current: state.conf, routes: mapRoutes(r), fontSize: 9,
    marks: r ? [{ k: r.a, text: 'A' }, { k: r.b, text: 'B' }] : [], color: k => (k === state.conf ? 'var(--series-2)' : 'var(--series-1)'), radius: () => 3.5 });
  s.append(g);
  const q = L.quality?.map;
  $('#map-cap').textContent = `Each disk is one secondary-structure topology (abstract shape); inside, distances follow C4′ RMSD. Click a point to morph.${q ? ` Trustworthiness ${fmt(q.trustworthiness, 2)}, topology purity ${fmt(q.topology_purity, 2)}.` : ''}`;
}
$('#map').addEventListener('click', e => {
  if (!mapPos) return;
  const s = $('#map'), b = s.getBoundingClientRect(), W = s.clientWidth || 300;
  const x = (e.clientX - b.left) * W / b.width, y = (e.clientY - b.top) * 220 / b.height;
  let best = -1, bd = 18;
  for (let k = 0; k < ens.n; k++) { const d = Math.hypot(mapPos.px(k) - x, mapPos.py(k) - y); if (d < bd) { bd = d; best = k; } }
  if (best >= 0) { state.playing && togglePlay(); goTo(best); }
});
function renderPathway() {
  const r = analysis.pathway;
  if (!r) return;
  $('#pw-box').hidden = false;
  const box = $('#pw-text');
  const same = r.A.dbn === r.B.dbn;
  const lines = [el('div', {}, el('b', { text: `${ens.shortId(r.a)} → ${ens.shortId(r.b)}` }), same ? ' · same secondary structure (tertiary rearrangement)' :
    r.same_basin ? ` · same 2D basin: A lies ${fmt(r.A.above_basin, 1)} and B ${fmt(r.B.above_basin, 1)} kcal/mol above its minimum, so no barrier separates them` :
    ` · effective barrier ${fmt(r.effective_barrier, 1)} kcal/mol (direct path ${fmt(r.direct.barrier, 1)}) · ΔG ${fmt(r.dG, 1)} kcal/mol · half-life ≈ ${fmtTime(r.half_life_s)} (k₀ = ${r.k0.toExponential(0)} s⁻¹)`)];
  const b = r.bottleneck;
  if (!same && !r.same_basin && ((r.A.above_basin ?? 0) > 0.05 || (r.B.above_basin ?? 0) > 0.05)) lines.push(el('div', { text: `Kinetics between basins: A lies ${fmt(r.A.above_basin, 1)} and B ${fmt(r.B.above_basin, 1)} kcal/mol above the minimum it relaxes into.` }));
  if (b && b.unit.length) lines.push(el('div', { text: `Cooperative unit: ${b.unit.length} base pairs open together before the highest saddle (${b.unit_helices.map(h => `${ens.resLabel(h[0][0])}–${ens.resLabel(h[h.length - 1][0])}/${ens.resLabel(h[h.length - 1][1])}–${ens.resLabel(h[0][1])}`).join(', ')}).` }));
  if (r.motion?.modules?.length) lines.push(el('div', { text: `Rigid modules: ${r.motion.modules.map(m => `M${m.id + 1} ${m.size} nt${m.id ? ` turns ${fmt(m.rotation, 0)}°` : ' (reference)'}`).join(' · ')}.` }));
  lines.push(el('div', { text: `${r.conditions_text}. Rates are order-of-magnitude estimates.` }));
  box.replaceChildren(...lines);
  // energy along the best route
  const s = $('#pw-prof'), W = s.clientWidth || 300, H = 110;
  s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  s.replaceChildren();
  const P = r.best.profile;
  if (P.length > 1) {
    const E = P.map(p => p.E), lo = Math.min(...E), hi = Math.max(...E);
    const X = i => 8 + i * (W - 16) / (P.length - 1), Y = v => H - 14 - (v - lo) / ((hi - lo) || 1) * (H - 30);
    s.append(svg('polyline', { points: P.map((p, i) => `${X(i)},${Y(p.E)}`).join(' '), fill: 'none', stroke: 'var(--series-2)', 'stroke-width': 2 }));
    const sd = b?.saddle_step;
    const e0 = P.find(p => p.phase !== 'relax_a')?.E ?? P[0].E;      // the basin of A
    if (sd !== undefined && sd !== null && P[sd]) s.append(svg('text', { x: X(sd), y: Y(P[sd].E) - 5, 'text-anchor': 'middle', 'font-size': 10, fill: 'var(--text-2)', text: `saddle +${fmt(P[sd].E - e0, 1)}` }));
    s.append(svg('text', { x: 8, y: H - 2, 'font-size': 10, fill: 'var(--text-3)', text: 'A' }), svg('text', { x: W - 8, y: H - 2, 'text-anchor': 'end', 'font-size': 10, fill: 'var(--text-3)', text: 'B' }));
  }
  $('#pw-play').hidden = !r.motion?.modules?.length;
  // phi scan heat map
  const sc = analysis.phi, ph = $('#pw-phi');
  if (sc) {
    const n = ens.nRes, W2 = ph.clientWidth || 300, cw = (W2 - 14) / n, ch = Math.max(4, Math.min(10, cw * 1.6)), H2 = 4 * ch + 4;
    ph.setAttribute('viewBox', `0 0 ${W2} ${H2}`);
    ph.setAttribute('height', H2);
    ph.replaceChildren();
    const col = v => { const t = Math.max(-1, Math.min(1, v / 4)); const a = Math.abs(t); return t < 0 ? `rgba(28,92,171,${0.1 + 0.9 * a})` : `rgba(207,84,34,${0.1 + 0.9 * a})`; };
    'ACGU'.split('').forEach((nt, row) => {
      ph.append(svg('text', { x: 0, y: row * ch + ch - 1, 'font-size': Math.min(9, ch), fill: 'var(--text-3)', text: nt }));
      for (let i = 0; i < n; i++) if (nt !== ens.sequence[i]) ph.append(svg('rect', { x: 12 + i * cw, y: row * ch, width: cw + 0.2, height: ch, fill: col(sc.matrix[i][nt].ddg_ts) }));
    });
    $('#pw-cap').textContent = 'φ-value scan: change of the barrier for every point mutation (blue = faster switching, orange = slower).';
  } else ph.hidden = true;
}
$('#pw-play').addEventListener('click', () => {
  const r = analysis.pathway;
  if (!r?.motion) return;
  state.playing && togglePlay();
  const mech = buildMechanism(ens.coords[r.a], ens.coords[r.b], ens.atomRes, r.motion, r.timing, { staged: true, names: ens.atomNames });
  state.prev = state.conf; state.conf = r.a;
  display.set(ens.coords[r.a]);
  viewer.setCoords(display);
  applyColors();
  const c0 = colorsFor(r.a).lin, c1 = colorsFor(r.b).lin, t0 = performance.now();
  state.morph = true;
  viewer._addAnimator(() => {
    const t = clamp((performance.now() - t0) / 6000, 0, 1);
    mech.at(t, display);
    viewer.setCoords(display);
    const cc = new Float32Array(c0.length);
    for (let i = 0; i < cc.length; i++) cc[i] = c0[i] + (c1[i] - c0[i]) * t;
    colorsLin = cc;
    viewer.setResidueColorsLinear(cc);
    if (t >= 1) { state.morph = null; state.prev = r.a; state.conf = r.b; applyColors(); renderList(); renderMap(); }
    return t < 1;
  });
});
new ResizeObserver(() => { renderMap(); }).observe($('#map'));

applyColors();
renderList();
renderMap();
renderPathway();
document.body.dataset.ready = '1';
