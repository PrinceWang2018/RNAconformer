// Sequence bar above the 3D view: nucleotides coloured like the structure, the dot-bracket of the shown
// conformation underneath (positions whose pairing changed since the previous conformation are marked),
// click / Shift+click / drag to select, hover to highlight in 3D and in the 2D window.
import { el, loadPref, savePref } from '../util.js';
import { tooltip, tipRow } from '../ui.js';

export class SequenceBar {
  constructor(app, host) {
    this.app = app;
    this.visible = loadPref('seqbar', true);
    this.track = el('div', { class: 'sq-track', translate: 'no' });
    this.root = el('div', { class: 'seqbar', role: 'listbox', 'aria-label': 'Sequence (click to select residues)' }, this.track);
    host.before(this.root);
    let drag = null;
    const resAt = e => { const t = e.target.closest?.('[data-r]'); return t ? +t.dataset.r : null; };
    this.track.addEventListener('pointerdown', e => {
      const r = resAt(e);
      if (r === null || !app.ens) return;
      e.preventDefault();
      if (e.shiftKey && this.anchor !== undefined && this.anchor !== null) { this._range(this.anchor, r, 'add'); return; }
      drag = { start: r, mode: e.ctrlKey || e.metaKey ? 'toggle' : 'set' };
      this.anchor = r;
      app.select([r], drag.mode);
      this.track.setPointerCapture(e.pointerId);
    });
    this.track.addEventListener('pointermove', e => {
      const r = drag ? this._resAtPoint(e) : resAt(e);
      if (drag && r !== null) { this._range(drag.start, r, drag.mode === 'toggle' ? 'add' : 'set'); return; }
      if (r === null) { if (this._hover !== null) { this._hover = null; app.hoverResidues([]); tooltip.hide(); } return; }
      if (r !== this._hover) { this._hover = r; app.hoverResidues([r]); }
      const ens = app.ens, pi = ens.pairInfo(app.state.conf, r);
      tooltip.show([tipRow(ens.resLabel(r), `position ${r + 1}`), el('div', { class: 'tl', text: pi ? `paired with ${ens.resLabel(pi.partner)}` : 'unpaired' })], e.clientX, e.clientY);
    });
    const end = () => { drag = null; };
    this.track.addEventListener('pointerup', end);
    this.track.addEventListener('pointercancel', end);
    this.track.addEventListener('pointerleave', () => { if (!drag && this._hover !== null) { this._hover = null; app.hoverResidues([]); tooltip.hide(); } });
    for (const ev of ['ensemble', 'conformation', 'colors', 'display']) app.on(ev, () => this.render());
    app.on('selection', () => this._marks());
    app.on('hover', list => this._hoverMarks(list));
    this._apply();
  }

  _resAtPoint(e) {
    const t = document.elementFromPoint(e.clientX, Math.min(e.clientY, this.track.getBoundingClientRect().bottom - 4));
    const c = t?.closest?.('[data-r]');
    return c ? +c.dataset.r : null;
  }

  _range(a, b, mode) {
    const lo = Math.min(a, b), hi = Math.max(a, b);
    this.app.select(Array.from({ length: hi - lo + 1 }, (_, i) => lo + i), mode);
  }

  setVisible(on) { this.visible = on; savePref('seqbar', on); this._apply(); this.render(); this.app.viewer.resize(); this.app.rebuildMenus?.(); }
  _apply() { this.root.hidden = !this.visible || !this.app.ens; }

  render() {
    this._apply();
    const app = this.app, ens = app.ens;
    if (!ens || this.root.hidden) return;
    const k = app.state.conf, prev = app.state.prev;
    const hex = app.colorsFor(k, prev).hex;
    const dbn = ens.conformations[k].dbn || '';
    const pdbn = prev !== null && prev !== k ? ens.conformations[prev].dbn || '' : null;
    const partner = ens.partner(k), pp = prev !== null && prev !== k ? ens.partner(prev) : null;
    const cells = [];
    for (let r = 0; r < ens.nRes; r++) {
      const num = ens.residues[r].num;
      const tick = r === 0 || (r + 1) % 10 === 0;
      const changed = pp && pp[r] !== partner[r];
      cells.push(el('span', { class: `sq-c${changed ? ' chg' : ''}`, 'data-r': r },
        el('i', { class: 'sq-n', text: tick ? String(num) : '' }),
        el('b', { class: 'sq-b', style: { background: hex[r], color: textOn(hex[r]) }, text: ens.residues[r].name }),
        el('i', { class: 'sq-d', text: dbn[r] || '.' })));
    }
    this.track.replaceChildren(el('span', { class: 'sq-lab' }, el('i', { text: '' }), el('b', { text: "5′" }), el('i', { text: pdbn !== null ? 'Δ' : '' })), ...cells, el('span', { class: 'sq-lab' }, el('i', { text: '' }), el('b', { text: "3′" }), el('i', { text: '' })));
    this.root.title = pdbn !== null ? `Dot-bracket of ${ens.shortId(k)}; underlined positions changed pairing since ${ens.shortId(prev)}` : `Dot-bracket of ${ens.shortId(k)}`;
    this._marks();
  }

  _marks() {
    const sel = this.app.state.selection;
    for (const c of this.track.querySelectorAll('.sq-c')) c.classList.toggle('sel', sel.has(+c.dataset.r));
  }
  _hoverMarks(list) {
    const set = new Set(list || []);
    for (const c of this.track.querySelectorAll('.sq-c.hov')) c.classList.remove('hov');
    for (const r of set) this.track.querySelector(`.sq-c[data-r="${r}"]`)?.classList.add('hov');
    if (set.size === 1) {
      const c = this.track.querySelector(`.sq-c[data-r="${[...set][0]}"]`);
      if (c) {
        const tr = this.track.getBoundingClientRect(), cr = c.getBoundingClientRect();
        if (cr.left < tr.left || cr.right > tr.right) c.scrollIntoView({ block: 'nearest', inline: 'center' });
      }
    }
  }
}

function textOn(hex) {
  const h = (hex || '#888888').replace('#', '');
  const v = parseInt(h, 16);
  const r = (v >> 16) & 255, g = (v >> 8) & 255, b = v & 255;
  return 0.299 * r + 0.587 * g + 0.114 * b > 150 ? '#141414' : '#ffffff';
}
