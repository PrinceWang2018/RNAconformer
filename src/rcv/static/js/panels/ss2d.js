// Secondary-structure corner window: the 2D drawing of the shown conformation, with nucleotides coloured
// like the 3D view and linked to it (hover, click-select, labels). When the conformation changes, the
// drawing morphs to the new structure in step with the 3D morph; gained pairs are emphasised.
import { el, svg, loadPref, savePref } from '../util.js';
import { icon } from '../icons.js';
import { tooltip, tipRow } from '../ui.js';
import { layoutFor, fit2D, blendLayouts, drawStructure2D, canvas2DTarget, orientForBox } from '../ss2dlayout.js';

export { layoutFor, fit2D, blendLayouts, drawStructure2D, canvas2DTarget } from '../ss2dlayout.js';

const SIZES = { s: [230, 190], m: [320, 260], l: [460, 380] };
const CORNERS = ['tr', 'br', 'bl', 'tl'];

// Backend that builds SVG elements (static exports, PDF through the chart walker).
export function svgTarget(root) {
  const st = (o, extra = {}) => ({ stroke: o.stroke || 'none', 'stroke-width': o.width || 1, 'stroke-dasharray': o.dash ? o.dash.join(' ') : null, 'stroke-linecap': 'round', 'stroke-linejoin': 'round', opacity: o.alpha ?? null, ...extra });
  return {
    line(x1, y1, x2, y2, o) { root.append(svg('line', { x1, y1, x2, y2, ...st(o) })); },
    path(cmds, o) { root.append(svg('path', { d: cmds.map(c => c[0] + ' ' + c.slice(1).map(v => +v.toFixed(2)).join(' ')).join(' '), fill: o.fill || 'none', ...st(o) })); },
    circle(cx, cy, r, o) { root.append(svg('circle', { cx, cy, r, fill: o.fill || 'none', stroke: o.stroke || 'none', 'stroke-width': o.width || 1, 'data-r': o.residue })); },
    text(x, y, s, o) { root.append(svg('text', { x, y, 'font-size': o.size, 'font-weight': o.weight || 500, fill: o.fill, 'text-anchor': o.anchor || 'start', 'data-r': o.residue, text: s })); },
  };
}

// ------------------------------------------------------------------ the corner window
export class SecondaryInset {
  constructor(app, viewport) {
    this.app = app;
    const p = loadPref('ss2d', {});
    this.visible = p.visible ?? true;
    this.size = SIZES[p.size] ? p.size : 'm';
    this.corner = CORNERS.includes(p.corner) ? p.corner : 'tr';
    this.hover = new Set();
    this.shown = null;              // Float64Array currently drawn
    this.anim = null;
    this.title = el('span', { class: 'ss-title' });
    this.sub = el('span', { class: 'ss-sub' });
    const btn = (ic, title, fn) => el('button', { class: 'ss-btn', title, 'aria-label': title, onclick: e => { e.stopPropagation(); fn(); } }, icon(ic));
    this.sizeBtn = btn('expand', 'Window size (small / medium / large)', () => this.cycleSize());
    this.cornerBtn = btn('corner', 'Move to another corner', () => this.cycleCorner());
    this.head = el('div', { class: 'ss-head', 'data-win-head': '' }, el('span', { class: 'ss-tt' }, this.title, this.sub), el('span', { class: 'grow' }),
      this.cornerBtn, this.sizeBtn, el('span', { class: 'win-btns', 'data-win': 'ss2d' }));
    this.maxed = false;
    this.svg = svg('svg', { class: 'ss-svg', role: 'img' });
    this.legend = el('div', { class: 'ss-legend' });
    this.root = el('div', { class: 'ss-inset', 'aria-label': 'Secondary structure of the shown conformation' }, this.head, this.svg, this.legend);
    viewport.append(this.root);
    new ResizeObserver(() => { if (this.maxed && !this.root.hidden) this.render(); }).observe(this.root);
    this.svg.addEventListener('pointermove', e => this._pointer(e));
    this.svg.addEventListener('pointerleave', () => { if (this.hover.size) { this.hover.clear(); app.hoverResidues([]); tooltip.hide(); this.draw(); } });
    this.svg.addEventListener('click', e => {
      const r = this._residueAt(e);
      if (r === null) return;
      app.select([r], e.shiftKey || e.ctrlKey || e.metaKey ? 'toggle' : 'set');
    });
    this.svg.addEventListener('dblclick', e => { const r = this._residueAt(e); if (r !== null && app.ens) app.viewer.centerOn(app.viewer.labelAnchor(r)); });
    app.on('ensemble', () => { this.shown = null; this.anim = null; this.render(); });
    app.on('conformation', ev => this._onConformation(ev));
    app.on('morph', m => this._onMorph(m));
    app.on('selection', () => this.draw());
    app.on('settings', () => this.draw());
    app.on('hover', list => {
      const next = new Set(list || []);
      if (next.size === this.hover.size && [...next].every(r => this.hover.has(r))) return;
      this.hover = next;
      this.draw();
    });
    app.on('colors', () => this.draw());
    this._apply();
  }

  _save() { savePref('ss2d', { visible: this.visible, size: this.size, corner: this.corner }); }
  _apply() {
    this.root.hidden = !this.visible || !this.app.ens;
    this.root.dataset.corner = this.corner;
    this.root.classList.toggle('maxed', this.maxed);
    this.cornerBtn.hidden = this.sizeBtn.hidden = this.maxed;
    const [w, h] = SIZES[this.size];
    this.root.style.width = this.maxed ? '' : `${w}px`;
    this.svg.style.height = this.maxed ? '' : `${h - 54}px`;
    const vp = this.root.parentElement;
    if (vp) {
      vp.classList.toggle('has-ss2d', !this.root.hidden && !this.maxed);
      vp.style.setProperty('--ss2d-w', `${w}px`);
      vp.style.setProperty('--ss2d-h', `${h}px`);
      vp.dataset.ss2dCorner = this.corner;
    }
  }
  setVisible(on) { this.visible = on; this._save(); this._apply(); this.render(); this.app.rebuildMenus?.(); this.app._syncToolbar?.(); }
  // Window ▸ Maximise: the drawing fills the window under the menu bar (Esc restores)
  setMaximized(on) { this.maxed = !!on; this._apply(); requestAnimationFrame(() => this.render()); }
  resetPlace() { this.size = 'm'; this.corner = 'tr'; this.visible = true; this.maxed = false; this._save(); this._apply(); this.render(); }
  cycleSize() { const ks = Object.keys(SIZES); this.size = ks[(ks.indexOf(this.size) + 1) % ks.length]; this._save(); this._apply(); this.draw(); }
  cycleCorner() { this.corner = CORNERS[(CORNERS.indexOf(this.corner) + 1) % CORNERS.length]; this._save(); this._apply(); }

  // conformation change: animate the drawing in step with the 3D morph
  _onConformation(ev) {
    const ens = this.app.ens;
    if (!ens || this.root.hidden) return;
    if (ev.scrub) { this.anim = null; return; }
    const target = layoutFor(ens, ev.to);
    if (!this.shown || this.shown.length !== target.length) { this.shown = target.slice(); this.render(); return; }
    this.anim = { from: this.shown.slice(), to: target, fit: fit2D(this.shown, target), src: ev.from, dst: ev.to };
  }
  _onMorph(m) {
    if (!this.app.ens || this.root.hidden) return;
    if (m.done) { this._finish(); return; }
    const a = this.anim;
    if (!a) {
      // scrubbing or video frames: blend the two layouts directly
      if (m.from === undefined) return;
      const A = layoutFor(this.app.ens, m.from), B = layoutFor(this.app.ens, m.to);
      this.shown = blendLayouts(A, B, m.t);
      this._frame = { from: m.from, t: m.t, to: m.to };
      this.draw();
      return;
    }
    this.shown = blendLayouts(a.from, a.to, m.t, a.fit, this.shown);
    this._frame = { from: a.src, t: m.t, to: a.dst };
    this.draw();
  }
  _finish() {
    const ens = this.app.ens;
    if (!ens) return;
    this.anim = null;
    this._frame = null;
    this.shown = layoutFor(ens, this.app.state.conf).slice();
    this.render();
  }

  render() {
    this._apply();
    const ens = this.app.ens;
    if (!ens || this.root.hidden) return;
    if (!this.shown || this.shown.length !== 2 * ens.nRes) this.shown = layoutFor(ens, this.app.state.conf).slice();
    this.draw();
  }

  draw() {
    const app = this.app, ens = app.ens;
    if (!ens || this.root.hidden) return;
    if (!this.shown || this.shown.length !== 2 * ens.nRes) { this.anim = null; this._frame = null; this.shown = layoutFor(ens, app.state.conf).slice(); }
    const k = app.state.conf, prev = app.state.prev;
    const fr = this._frame;
    const shownK = fr ? (fr.t >= 0.5 ? fr.to : fr.from) : k;
    const colors = app.colorsFor(shownK, prev).hex;
    const w = this.svg.clientWidth || SIZES[this.size][0] - 2, h = this.svg.clientHeight || SIZES[this.size][1] - 54;
    this.svg.setAttribute('viewBox', `0 0 ${w} ${h}`);
    const g = svg('g');
    const dark = document.documentElement.dataset.theme === 'dark';
    this._geom = drawStructure2D(svgTarget(g), ens, orientForBox(ens, this.shown, w, h), { k: fr ? fr.to : k, prev, from: fr ? fr.from : null, t: fr ? fr.t : 1, colors, sel: app.state.selection, hover: this.hover, box: [4, 4, w - 8, h - 8], dark });
    this.svg.replaceChildren(g);
    const c = ens.conformations[k];
    const npk = c.pairs.filter(p => p[3]).length;
    this.title.textContent = '2D';
    this.sub.textContent = ` · ${ens.shortId(fr ? fr.to : k)} · ${c.pairs.length} bp${npk ? ` (${npk} pseudoknotted)` : ''}`;
    const gained = prev !== null && prev !== k ? c.pairs.filter(p => !ens.conformations[prev].pairs.some(q => q[0] === p[0] && q[1] === p[1])).length : 0;
    const key = (color, label, dash) => el('span', { class: 'lg' }, el('span', { class: 'ln', style: { borderTop: `2px ${dash ? 'dashed' : 'solid'} ${color}` } }), label);
    this.legend.replaceChildren(...[key('var(--text-2)', 'pair'), gained ? key('#eb6834', `${gained} new vs ${ens.shortId(prev)}`) : null, npk ? key('#6250d6', 'pseudoknot', true) : null].filter(Boolean));
  }

  _residueAt(e) {
    const t = e.target.closest?.('[data-r]');
    if (t) return +t.getAttribute('data-r');
    if (!this._geom || !this.app.ens) return null;
    const rect = this.svg.getBoundingClientRect();
    const vb = this.svg.viewBox.baseVal;
    const x = (e.clientX - rect.left) * (vb.width / rect.width), y = (e.clientY - rect.top) * (vb.height / rect.height);
    let best = null, bd = Infinity;
    for (let i = 0; i < this.app.ens.nRes; i++) { const d = Math.hypot(this._geom.X(i) - x, this._geom.Y(i) - y); if (d < bd) { bd = d; best = i; } }
    return bd <= Math.max(6, this._geom.r * 1.6) ? best : null;
  }
  _pointer(e) {
    const app = this.app, ens = app.ens;
    const r = this._residueAt(e);
    if (r === null) { if (this.hover.size) { this.hover.clear(); app.hoverResidues([]); tooltip.hide(); this.draw(); } return; }
    if (!this.hover.has(r) || this.hover.size !== 1) { this.hover = new Set([r]); app.hoverResidues([r], { from2d: true }); this.draw(); }
    const pi = ens.pairInfo(app.state.conf, r);
    tooltip.show([tipRow(ens.resLabel(r), 'click to select · double-click to centre'),
      el('div', { class: 'tl', text: pi ? `paired with ${ens.resLabel(pi.partner)}${pi.pk ? ' (pseudoknot)' : ''}` : 'unpaired in this conformation' })], e.clientX, e.clientY);
  }

  // Static rendering for figures and video frames (canvas context, box in pixels).
  drawTo(ctx, box, { k = this.app.state.conf, from = null, t = 1, dark = false, frame = true, unit = 1 } = {}) {
    const app = this.app, ens = app.ens;
    if (!ens) return;
    const [bw, bh] = [box[2], box[3]];
    const coords = orientForBox(ens, from !== null && from !== k && t < 1 ? blendLayouts(layoutFor(ens, from), layoutFor(ens, k), t) : layoutFor(ens, k), bw, bh);
    const colors = app.colorsFor(from !== null && t < 0.5 ? from : k).hex;
    const [x, y, w, h] = box;
    if (frame) {
      ctx.save();
      ctx.fillStyle = dark ? 'rgba(26,26,25,0.82)' : 'rgba(255,255,255,0.86)';
      ctx.strokeStyle = dark ? 'rgba(255,255,255,0.14)' : 'rgba(0,0,0,0.12)';
      ctx.lineWidth = Math.max(1, h / 300);
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(x, y, w, h, Math.min(w, h) * 0.04); else ctx.rect(x, y, w, h);
      ctx.fill(); ctx.stroke();
      ctx.restore();
    }
    const pad = Math.min(w, h) * 0.04;
    drawStructure2D(canvas2DTarget(ctx), ens, coords, { k, prev: from ?? app.state.prev, from, t, colors, box: [x + pad, y + pad, w - 2 * pad, h - 2 * pad], dark, sel: new Set(), unit });
  }
}
