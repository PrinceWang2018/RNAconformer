// Mouse tools of the 3D view. The left button does what the active tool says; dragging on empty space always rotates.
//   Select (V)   click residues to select them (Shift adds)
//   Pull (G)     drag a residue: the structure stretches and switches to the nearest stable conformation
//   Paint (B)    click or drag over residues (or atoms) to colour them with the paint colour; Shift = scheme colour
//   Label (T)    click to add a label (Alt: atom label; inside a multi-residue selection: region label); drag labels
//   Measure (M)  click two atoms for a distance
import { el, $ } from './util.js';
import { icon } from './icons.js';
import { popupMenu, toast } from './ui.js';
import { SLOTS } from './colors.js';

export const TOOLS = [
  ['select', 'Select', 'cursor', 'V', 'Select residues (click; Shift adds)'],
  ['pull', 'Pull', 'hand', 'G', 'Pull a residue to switch conformation'],
  ['paint', 'Paint', 'brush', 'B', 'Paint residues or atoms with a colour (drag; Shift restores the scheme colour)'],
  ['label', 'Label', 'tag', 'T', 'Add and edit labels (click the structure; Alt = atom label; drag labels to move)'],
  ['measure', 'Measure', 'ruler', 'M', 'Measure: click two atoms'],
];

export class MouseTools {
  constructor(app) {
    this.app = app;
    this.paintColor = app.prefs.paintColor || '#eb6834';
    this.paintAtoms = false;
    this.stroke = null;
    const seg = el('div', { class: 'seg tool-seg', role: 'radiogroup', 'aria-label': 'Mouse tool' });
    for (const [id, label, ic, key, tip] of TOOLS) seg.append(el('button', { 'data-tool': id, title: `${tip} (${key})`, 'aria-label': label }, icon(ic), el('span', { text: label })));
    seg.addEventListener('click', e => { const t = e.target.closest('button')?.dataset.tool; if (t) this.set(t === app.tool && t !== 'select' ? 'select' : t); });
    this.chip = el('button', { class: 'paint-chip', title: 'Paint colour', 'aria-label': 'Paint colour', onclick: e => this.colorMenu(e.currentTarget) }, el('span'));
    this.root = el('div', { class: 'tb-group tool-group' }, seg, this.chip);
    this.seg = seg;
    app.tool = 'select';
    this._sync();
  }

  set(tool) {
    const app = this.app;
    if (app.tool === tool) return;
    if (app.tool === 'pull') app.pullCtl.setMode(false);
    if (app.tool === 'measure') app.toggleMeasureMode(false, { quiet: true });
    app.tool = tool;
    if (tool === 'pull') app.pullCtl.setMode(true);
    if (tool === 'measure') app.toggleMeasureMode(true, { quiet: true });
    if (tool !== 'label') app.annotations.closeEditor();
    app.viewer.container.dataset.tool = tool;
    this._sync();
    app.viewer.needsRender = true;
    app.rebuildMenus?.();
    app.emit('tool', tool);
  }

  _sync() {
    for (const b of this.seg.children) { const on = b.dataset.tool === this.app.tool; b.classList.toggle('on', on); b.setAttribute('aria-checked', on ? 'true' : 'false'); }
    this.chip.hidden = this.app.tool !== 'paint';
    this.chip.firstChild.style.background = this.paintColor;
    this.chip.title = `Paint colour ${this.paintColor} · ${this.paintAtoms ? 'atoms' : 'residues'}`;
  }

  colorMenu(anchor) {
    const app = this.app;
    const pick = c => { this.paintColor = c; app.prefs.paintColor = c; app.savePrefs(); this._sync(); };
    popupMenu(anchor, [
      { head: 'Paint colour (active palette)' },
      ...SLOTS.slice(0, 10).map((c, i) => ({ label: `● ${c}`, checked: this.paintColor === c, action: () => pick(c), color: c })),
      { label: 'Black', checked: this.paintColor === '#141414', action: () => pick('#141414') },
      { label: 'White', checked: this.paintColor === '#ffffff', action: () => pick('#ffffff') },
      { label: 'Light gray', checked: this.paintColor === '#c4c2ba', action: () => pick('#c4c2ba') },
      { label: 'Other colour…', action: () => { const inp = el('input', { type: 'color', value: this.paintColor }); inp.addEventListener('input', () => pick(inp.value)); inp.click(); } },
      { sep: true },
      { label: 'Paint residues', checked: !this.paintAtoms, action: () => { this.paintAtoms = false; this._sync(); } },
      { label: 'Paint single atoms', checked: this.paintAtoms, action: () => { this.paintAtoms = true; this._sync(); if (app.viewer.settings.representation === 'cartoon') toast('Atoms are visible in the Cartoon + atoms, Sticks, Ball & stick and Spheres representations.', 'info', 'Paint atoms'); } },
      { sep: true },
      { label: 'Palettes…', action: () => app.openPalettes() },
    ]);
    // colour dots in the menu
    for (const b of document.querySelectorAll('.menu-pop .menu-item')) {
      const m = /● (#[0-9a-f]{6})/i.exec(b.textContent);
      if (m) { const lab = b.children[1]; lab.textContent = m[1]; lab.prepend(el('i', { class: 'menu-dot', style: { background: m[1] } })); }
    }
  }

  // Called for every left-button press on the canvas (after the Pull tool had its chance). true = gesture taken.
  onPointerDown(e) {
    const app = this.app;
    if (e.button !== 0 || !app.ens) return false;
    if (app.tool === 'paint') {
      const p = app.viewer.pick(e.clientX, e.clientY);
      if (!p) return false;
      app.history?.commit();
      this.stroke = { atoms: new Set(), residues: new Set(), erase: e.shiftKey, pointerId: e.pointerId };
      this._paintAt(p);
      const cv = app.viewer.canvas;
      try { cv.setPointerCapture(e.pointerId); } catch { /* synthetic */ }
      const move = ev => { const q = app.viewer.pick(ev.clientX, ev.clientY); if (q) this._paintAt(q); };
      const up = () => { cv.removeEventListener('pointermove', move); cv.removeEventListener('pointerup', up); cv.removeEventListener('pointercancel', up); this._commit(); };
      cv.addEventListener('pointermove', move);
      cv.addEventListener('pointerup', up);
      cv.addEventListener('pointercancel', up);
      return true;
    }
    if (app.tool === 'label') {
      const p = app.viewer.pick(e.clientX, e.clientY);
      if (!p) { app.annotations.selected = null; app.annotations.closeEditor(); app.viewer.needsRender = true; return false; }
      app.annotations.onPick(p, e);
      return true;
    }
    return false;
  }

  _paintAt(p) {
    const S = this.stroke, app = this.app, d = app.displayOf();
    if (this.paintAtoms) {
      if (S.atoms.has(p.atom)) return;
      S.atoms.add(p.atom);
      if (S.erase) d.atomColors.delete(p.atom); else d.atomColors.set(p.atom, this.paintColor);
      app.viewer.setDisplaySets({ atomColors: d.atomColors });
    } else {
      if (S.residues.has(p.residue)) return;
      S.residues.add(p.residue);
      if (S.erase) d.colors.delete(p.residue); else d.colors.set(p.residue, this.paintColor);
      app._colorCache.clear();
      app.refreshColors();
    }
  }

  // one stroke = one logged command (reproducible, undoable)
  _commit() {
    const S = this.stroke, app = this.app;
    this.stroke = null;
    if (!S) return;
    const ens = app.ens;
    if (S.residues.size) {
      const sel = `pos ${[...S.residues].sort((a, b) => a - b).map(r => r + 1).join('+')}`;
      (ens.commandLog ||= []).push(`color ${S.erase ? 'default' : this.paintColor}, ${sel}`);
    }
    if (S.atoms.size) (ens.commandLog ||= []).push(`# painted ${S.atoms.size} atom(s) ${S.erase ? 'back to the scheme colour' : this.paintColor}`);
    app.applyDisplay();
  }
}
