// Annotations for figures: text labels with optional arrows anchored to a residue or an atom, region labels that mark
// a group of residues, and free notes placed on the view (titles, captions). Anchored annotations follow every morph
// and conformation; all of them are drawn in image, panel and video exports and saved in .rcv files.
// Label tool (T): click the structure to add a label, drag a label to move it, double-click to edit, Delete to remove.
import { el, svg, figFont, FIG_FONT, clamp } from './util.js';
import { icon } from './icons.js';
import { SLOTS } from './colors.js';

let counter = 0;
const SIZES = [9, 11, 13, 16, 20, 26];

export class Annotations {
  constructor(app) {
    this.app = app;
    this.selected = null;
    this.svg = svg('svg', { class: 'ann-layer' });
    app.viewer.container.append(this.svg);
    app.viewer.renderHooks.add((proj, w, h) => this._draw(proj, w, h));
    app.registerSessionPart('annotations', {
      save: ens => (ens?.annotations || []).map(a => ({ ...a })),
      load: (ens, list) => { ens.annotations = (list || []).map(a => ({ ...a })); this.app.viewer.needsRender = true; this.app.emit('annotations'); },
    });
    window.addEventListener('keydown', e => {
      if (e.target.closest?.('input, textarea, select') || !this.selected) return;
      if (e.key === 'Delete' || e.key === 'Backspace') { this.remove(this.selected.id); e.preventDefault(); }
      if (e.key === 'Enter') { this.edit(this.selected.id); e.preventDefault(); }
    });
    app.on('ensemble', () => { this.selected = null; this.closeEditor(); });
  }

  list(ens = this.app.ens) { return ens ? (ens.annotations ||= []) : []; }
  get interactive() { return this.app.tool === 'label'; }

  // ---------------------------------------------------------------- model
  add(a) {
    const ens = this.app.ens;
    if (!ens) return null;
    this.app.history?.commit();                 // a new annotation is its own undo step
    const item = { id: `a${Date.now().toString(36)}${counter++}`, kind: 'label', text: '', color: '#141414', size: 13, bold: true, arrow: true, box: true, dx: 34, dy: -38, ...a };
    this.list().push(item);
    this.changed();
    return item;
  }
  update(id, patch) { const a = this.list().find(x => x.id === id); if (a) { Object.assign(a, patch); this.changed(); } return a; }
  remove(id) {
    this.app.history?.commit();
    const L = this.list(), i = L.findIndex(x => x.id === id);
    if (i >= 0) L.splice(i, 1);
    if (this.selected?.id === id) this.selected = null;
    this.closeEditor();
    this.changed();
  }
  clear() { if (this.app.ens) this.app.ens.annotations = []; this.selected = null; this.closeEditor(); this.changed(); }
  changed() { this.app.viewer.needsRender = true; this.app.emit('annotations'); }

  defaultText(r, atom) {
    const ens = this.app.ens;
    return atom !== undefined && atom !== null ? `${ens.resLabel(r)}:${ens.atomNames[atom]}` : ens.resLabel(r);
  }

  // anchor in model coordinates for the shown coordinates
  anchor(a, coords = this.app.display) {
    const ens = this.app.ens;
    if (a.kind === 'note') return null;
    if (a.atom !== undefined && a.atom !== null && a.atom >= 0) return [coords[3 * a.atom], coords[3 * a.atom + 1], coords[3 * a.atom + 2]];
    const res = a.kind === 'region' ? a.residues : [a.residue];
    const v = this.app.viewer;
    const pts = res.map(r => (coords === v.coords ? v.labelAnchor(r) : ens.residueCenter(coords, r)));
    return [0, 1, 2].map(i => pts.reduce((s, p) => s + p[i], 0) / pts.length);
  }

  // ---------------------------------------------------------------- live layer (SVG over the 3D view)
  _draw(proj, w, h) {
    const ens = this.app.ens, root = this.svg;
    root.setAttribute('viewBox', `0 0 ${w} ${h}`);
    root.classList.toggle('interactive', this.interactive);
    const list = this.list(ens);
    while (root.children.length > list.length) root.lastChild.remove();
    list.forEach((a, i) => {
      let g = root.children[i];
      if (!g || g.dataset.id !== a.id) {
        g = svg('g', { class: 'ann', 'data-id': a.id });
        if (root.children[i]) root.children[i].replaceWith(g); else root.append(g);
        this._bind(g);
      }
      g.replaceChildren();
      this.paint(svgTarget(g), a, proj, w, h, 1, ens, this.app.display);
      g.classList.toggle('sel', this.selected?.id === a.id);
    });
  }

  // Draw one annotation into a target (live SVG layer, or a canvas for exports). s = pixel scale.
  paint(T, a, proj, w, h, s, ens, coords) {
    const size = (a.size || 13) * s;
    const font = figFont(size, a.bold);
    let tx, ty, ax = null, ay = null;
    if (a.kind === 'note') { tx = (a.x ?? 0.5) * w; ty = (a.y ?? 0.08) * h; }
    else {
      const p = this.anchorFor(a, ens, coords);
      if (!p) return;
      const [px, py, pz] = proj(...p);
      if (pz < -1 || pz > 1) return;
      ax = px; ay = py;
      tx = px + (a.dx ?? 30) * s; ty = py + (a.dy ?? -30) * s;
      if (a.kind === 'region') {
        for (const r of a.residues) {
          const q = proj(...ens.residueCenter(coords, r));
          T.circle(q[0], q[1], 5.5 * s, { fill: a.color, alpha: 0.28 });
        }
      }
    }
    const lines = String(a.text || '').split('\n');
    const lw = Math.max(...lines.map(l => T.measure(l, font)), size * 0.6);
    const lh = size * 1.18, bw = lw + size * 0.7, bh = lh * lines.length + size * 0.35;
    const bx = tx - bw / 2, by = ty - bh / 2;
    if (ax !== null && a.arrow) {
      // arrow from the label box edge to the anchor
      const cx = clamp(ax, bx, bx + bw), cy = clamp(ay, by, by + bh);
      const dx = ax - cx, dy = ay - cy, L = Math.hypot(dx, dy);
      if (L > 6 * s) {
        const ux = dx / L, uy = dy / L, head = Math.min(9 * s, L * 0.4), end = [ax - ux * 3 * s, ay - uy * 3 * s];
        T.line(cx, cy, end[0] - ux * head * 0.6, end[1] - uy * head * 0.6, { stroke: a.color, width: 1.6 * s });
        T.poly([[end[0], end[1]], [end[0] - ux * head - uy * head * 0.45, end[1] - uy * head + ux * head * 0.45], [end[0] - ux * head + uy * head * 0.45, end[1] - uy * head - ux * head * 0.45]], { fill: a.color });
      }
    }
    if (a.box) T.rect(bx, by, bw, bh, 4 * s, { fill: 'rgba(255,255,255,0.88)', stroke: a.color, width: 1 * s });
    lines.forEach((l, k) => T.text(tx, by + size * 0.18 + lh * (k + 0.5) + size * 0.36, l, { font, fill: a.color, anchor: 'middle' }));
    if (!a.box) {
      // invisible hit area so a bare text label can be grabbed
      T.rect(bx, by, bw, bh, 0, { fill: 'rgba(0,0,0,0)', hit: true });
    }
  }

  anchorFor(a, ens, coords) {
    if (a.atom !== undefined && a.atom !== null && a.atom >= 0) return [coords[3 * a.atom], coords[3 * a.atom + 1], coords[3 * a.atom + 2]];
    const res = a.kind === 'region' ? a.residues : [a.residue];
    if (!res || !res.length || res.some(r => r === undefined || r >= ens.nRes)) return null;
    const pts = res.map(r => {
      const g = ens.glyco[r] >= 0 ? ens.glyco[r] : ens.trace[r];
      if (ens.rings[r]) { const o = ens.rings[r].outline; return [0, 1, 2].map(i => o.reduce((s2, k) => s2 + coords[3 * k + i], 0) / o.length); }
      return g >= 0 ? [coords[3 * g], coords[3 * g + 1], coords[3 * g + 2]] : ens.residueCenter(coords, r);
    });
    return [0, 1, 2].map(i => pts.reduce((s2, p) => s2 + p[i], 0) / pts.length);
  }

  // Exports: draw every annotation on a canvas rendered with `project` (pixel scale s relative to the screen).
  drawTo(ctx, project, w, h, s, coords = this.app.display) {
    const ens = this.app.ens;
    for (const a of this.list(ens)) this.paint(canvasTarget(ctx), a, project, w, h, s, ens, coords);
  }

  // ---------------------------------------------------------------- interaction
  _bind(g) {
    let drag = null;
    g.addEventListener('pointerdown', e => {
      if (!this.interactive || e.button !== 0) return;
      e.stopPropagation(); e.preventDefault();
      const a = this.list().find(x => x.id === g.dataset.id);
      if (!a) return;
      this.selected = a;
      drag = { x: e.clientX, y: e.clientY, dx: a.dx, dy: a.dy, nx: a.x, ny: a.y, moved: false };
      g.setPointerCapture(e.pointerId);
      this.app.viewer.needsRender = true;
    });
    g.addEventListener('pointermove', e => {
      if (!drag) return;
      const a = this.selected;
      const ddx = e.clientX - drag.x, ddy = e.clientY - drag.y;
      if (Math.abs(ddx) + Math.abs(ddy) > 2) drag.moved = true;
      if (a.kind === 'note') { const [w, h] = this.app.viewer.cssSize; a.x = clamp(drag.nx + ddx / w, 0.02, 0.98); a.y = clamp(drag.ny + ddy / h, 0.02, 0.98); }
      else { a.dx = drag.dx + ddx; a.dy = drag.dy + ddy; }
      this.app.viewer.needsRender = true;
    });
    g.addEventListener('pointerup', () => { if (drag?.moved) this.changed(); drag = null; });
    g.addEventListener('dblclick', e => { if (!this.interactive) return; e.stopPropagation(); this.edit(g.dataset.id); });
    g.addEventListener('contextmenu', e => {
      e.preventDefault(); e.stopPropagation();
      this.selected = this.list().find(x => x.id === g.dataset.id);
      this.edit(g.dataset.id);
    });
  }

  // Label tool: click on the structure (atom labels with Alt or in atom representations).
  onPick(p, e) {
    const app = this.app, ens = app.ens;
    const atomLevel = e.altKey || app.viewer.settings.representation !== 'cartoon';
    const sel = app.state.selection;
    let a;
    if (sel.size > 1 && sel.has(p.residue) && !atomLevel) {
      const res = [...sel].sort((x, y) => x - y);
      a = this.add({ kind: 'region', residues: res, text: `${ens.resLabel(res[0])}–${ens.resLabel(res[res.length - 1])}`, color: SLOTS[1], arrow: false, dx: 0, dy: -46 });
    } else {
      a = this.add({ kind: 'label', residue: p.residue, atom: atomLevel ? p.atom : null, text: this.defaultText(p.residue, atomLevel ? p.atom : null), color: '#141414' });
    }
    this.selected = a;
    this.edit(a.id, true);
    return a;
  }

  addNote(text = 'Title') {
    const a = this.add({ kind: 'note', text, x: 0.5, y: 0.07, size: 20, bold: true, arrow: false, box: false, color: '#141414' });
    this.selected = a;
    this.edit(a.id, true);
    return a;
  }

  // ---------------------------------------------------------------- editor popover
  closeEditor() { this.editor?.remove(); this.editor = null; }

  edit(id, fresh = false) {
    const a = this.list().find(x => x.id === id);
    if (!a) return;
    this.closeEditor();
    const app = this.app;
    const ta = el('textarea', { class: 'input', rows: 2, 'aria-label': 'Annotation text' });
    ta.value = a.text;
    ta.addEventListener('input', () => { a.text = ta.value; app.viewer.needsRender = true; });
    const swatches = el('div', { class: 'ann-sw' }, ...['#141414', '#ffffff', ...SLOTS.slice(0, 8)].map(c => el('button', { class: `ann-c ${a.color === c ? 'on' : ''}`, style: { background: c }, title: c,
      onclick: () => { a.color = c; app.viewer.needsRender = true; for (const b of swatches.children) b.classList.toggle('on', b.title === c); } })));
    const custom = el('input', { type: 'color', value: /^#[0-9a-f]{6}$/i.test(a.color) ? a.color : '#141414', title: 'Other colour' });
    custom.addEventListener('input', () => { a.color = custom.value; app.viewer.needsRender = true; });
    swatches.append(custom);
    const size = el('input', { type: 'range', min: 0, max: SIZES.length - 1, step: 1, value: Math.max(0, SIZES.findIndex(s => s >= a.size)), 'aria-label': 'Text size' });
    size.addEventListener('input', () => { a.size = SIZES[+size.value]; app.viewer.needsRender = true; });
    const tog = (key, label) => { const b = el('button', { class: `btn small ${a[key] ? 'on' : ''}`, text: label, onclick: () => { a[key] = !a[key]; b.classList.toggle('on', a[key]); app.viewer.needsRender = true; } }); return b; };
    const done = () => { this.closeEditor(); if (!a.text.trim()) this.remove(a.id); else this.changed(); };
    const box = el('div', { class: 'ann-editor', role: 'dialog', 'aria-label': 'Edit annotation' },
      el('div', { class: 'ann-head' }, el('b', { text: a.kind === 'note' ? 'Note' : a.kind === 'region' ? 'Region label' : 'Label' }), el('span', { class: 'grow' }),
        el('button', { class: 'ss-btn', title: 'Delete', onclick: () => this.remove(a.id) }, icon('trash')), el('button', { class: 'ss-btn', title: 'Done', onclick: done }, icon('check'))),
      ta, swatches,
      el('div', { class: 'ann-row' }, el('span', { class: 'muted', text: 'Size' }), size, tog('bold', 'Bold'), a.kind !== 'note' ? tog('arrow', 'Arrow') : null, tog('box', 'Box')));
    ta.addEventListener('keydown', e => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); done(); } if (e.key === 'Escape') { e.preventDefault(); done(); } e.stopPropagation(); });
    app.viewer.container.append(box);
    this.editor = box;
    // place near the label
    const [w, h] = app.viewer.cssSize;
    const g = this.svg.querySelector(`g[data-id="${a.id}"]`);
    const r = g?.getBBox ? g.getBBox() : { x: w / 2, y: h / 2, width: 0, height: 0 };
    box.style.left = `${clamp(r.x + r.width + 12, 8, w - 280)}px`;
    box.style.top = `${clamp(r.y - 10, 8, h - 190)}px`;
    setTimeout(() => { ta.focus(); if (fresh) ta.select(); }, 0);
  }
}

// ------------------------------------------------------------------ drawing targets
function svgTarget(g) {
  let measureCtx = null;
  return {
    measure(text, font) { measureCtx = measureCtx || document.createElement('canvas').getContext('2d'); measureCtx.font = font; return measureCtx.measureText(text).width; },
    line(x1, y1, x2, y2, o) { g.append(svg('line', { x1, y1, x2, y2, stroke: o.stroke, 'stroke-width': o.width, 'stroke-linecap': 'round' })); },
    poly(pts, o) { g.append(svg('path', { d: 'M' + pts.map(p => p.join(' ')).join(' L') + ' Z', fill: o.fill })); },
    rect(x, y, w, h, r, o) { g.append(svg('rect', { x, y, width: w, height: h, rx: r, fill: o.fill, stroke: o.stroke || null, 'stroke-width': o.width || null, class: o.hit ? 'ann-hit' : null })); },
    circle(cx, cy, r, o) { g.append(svg('circle', { cx, cy, r, fill: o.fill, opacity: o.alpha ?? null, 'pointer-events': 'none' })); },
    text(x, y, s, o) { g.append(svg('text', { x, y, 'text-anchor': o.anchor, fill: o.fill, style: `font: ${o.font}`, text: s })); },
  };
}

export function canvasTarget(ctx) {
  return {
    measure(text, font) { ctx.font = font; return ctx.measureText(text).width; },
    line(x1, y1, x2, y2, o) { ctx.save(); ctx.strokeStyle = o.stroke; ctx.lineWidth = o.width; ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); ctx.restore(); },
    poly(pts, o) { ctx.save(); ctx.fillStyle = o.fill; ctx.beginPath(); pts.forEach((p, i) => (i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]))); ctx.closePath(); ctx.fill(); ctx.restore(); },
    rect(x, y, w, h, r, o) {
      if (o.hit) return;
      ctx.save(); ctx.beginPath(); if (ctx.roundRect) ctx.roundRect(x, y, w, h, r); else ctx.rect(x, y, w, h);
      if (o.fill) { ctx.fillStyle = o.fill; ctx.fill(); } if (o.stroke) { ctx.strokeStyle = o.stroke; ctx.lineWidth = o.width; ctx.stroke(); } ctx.restore();
    },
    circle(cx, cy, r, o) { ctx.save(); ctx.globalAlpha = o.alpha ?? 1; ctx.fillStyle = o.fill; ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2); ctx.fill(); ctx.restore(); },
    text(x, y, s, o) { ctx.save(); ctx.font = o.font; ctx.fillStyle = o.fill; ctx.textAlign = o.anchor === 'middle' ? 'center' : o.anchor; ctx.textBaseline = 'alphabetic'; ctx.fillText(s, x, y); ctx.restore(); },
  };
}

export { FIG_FONT };
