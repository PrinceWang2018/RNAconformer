// RCV Lite: the interactive page written by "Export ▸ Interactive HTML". One self-contained file (no server, no
// install) for supplementary material, reviewers and collaborators: rotate, morph between conformations, play the
// tour, change representation and colours, read energies and the 2D structure of every conformation.
import { Viewer } from './viewer/viewer.js';
import { Ensemble } from './ensemble.js';
import { residueColors, legend as schemeLegend } from './colors.js';
import { layoutFor, blendLayouts, drawStructure2D, canvas2DTarget } from './ss2dlayout.js';
import { el, easeInOut, fmt, clamp } from './util.js';

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
    if (t >= 1) { state.morph = null; applyColors(); renderList(); if (state.playing) playTimer = setTimeout(step, 900); }
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
applyColors();
renderList();
document.body.dataset.ready = '1';
