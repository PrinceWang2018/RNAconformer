// Colour palettes for figures: built-in palettes (RCV, colour-blind-safe sets published by Okabe & Ito and by Paul
// Tol, pastel, grayscale for print) and the user's own. The active palette drives helix, structure-element and
// nucleotide colours, the paint tool, annotations, measurements and the plot studio. A colour-vision check flags
// colours that become hard to tell apart for colour-blind readers.
import { el, loadPref, savePref, rgbToHex } from './util.js';
import { icon } from './icons.js';
import { dialog, toast } from './ui.js';
import { applyPalette, SLOTS } from './colors.js';

export { BUILTIN, checkPalette } from './palettecheck.js';
import { BUILTIN, checkPalette } from './palettecheck.js';

// ------------------------------------------------------------------ manager
export class Palettes {
  constructor(app) {
    this.app = app;
    this.user = loadPref('palettes', []);
    this.activeId = loadPref('palette', 'rcv');
    this.apply(this.activeId, { quiet: true });
  }
  all() { return [...BUILTIN, ...this.user]; }
  get active() { return this.all().find(p => p.id === this.activeId) || BUILTIN[0]; }
  colors() { return [...SLOTS]; }
  apply(id, { quiet = false } = {}) {
    const p = this.all().find(x => x.id === id) || BUILTIN[0];
    this.activeId = p.id;
    savePref('palette', p.id);
    applyPalette(p);
    if (!quiet) { this.app._colorCache?.clear(); this.app.refreshColors?.(); this.app.emit?.('settings'); this.app.emit?.('palette'); }
  }
  saveUser() { savePref('palettes', this.user); }
}

export function paletteDialog(app) {
  const P = app.palettes;
  let selId = P.activeId;
  const listBox = el('div', { class: 'pal-list' });
  const editor = el('div', { class: 'pal-editor' });
  const swatches = cols => el('span', { class: 'pal-sw' }, ...cols.slice(0, 10).map(c => el('i', { style: { background: c } })));
  const render = () => {
    listBox.replaceChildren(...P.all().map(p => {
      const row = el('button', { class: `pal-item ${p.id === selId ? 'on' : ''}`, onclick: () => { selId = p.id; render(); } },
        el('span', { class: 'pal-name', text: p.name }), swatches(p.colors), p.id === P.activeId ? el('span', { class: 'badge', text: 'in use' }) : null);
      return row;
    }));
    const p = P.all().find(x => x.id === selId) || BUILTIN[0];
    const userPal = !BUILTIN.includes(p);
    const chk = checkPalette(p.colors);
    const fmtChk = (k, lab) => {
      const w = chk[k];
      if (!Number.isFinite(w.d)) return null;
      const bad = w.d < (k === 'normal' ? 15 : 8);
      return el('div', { class: bad ? 'state-bad' : 'state-ok' }, `${lab}: closest pair ${p.colors[w.i]} / ${p.colors[w.j]} — ΔE ${w.d.toFixed(1)}${bad ? ' (hard to tell apart)' : ''}`);
    };
    const colorsRow = el('div', { class: 'pal-colors' });
    p.colors.forEach((c, i) => {
      const inp = el('input', { type: 'color', value: c, disabled: userPal ? null : true, title: `Colour ${i + 1}` });
      inp.addEventListener('input', () => { p.colors[i] = inp.value; P.saveUser(); if (P.activeId === p.id) P.apply(p.id); });
      const rm = userPal && p.colors.length > 4 ? el('button', { class: 'ss-btn', title: 'Remove colour', onclick: () => { p.colors.splice(i, 1); P.saveUser(); render(); } }, icon('x')) : null;
      colorsRow.append(el('span', { class: 'pal-c' }, inp, rm));
    });
    if (userPal && p.colors.length < 12) colorsRow.append(el('button', { class: 'btn small', onclick: () => { p.colors.push(p.colors[p.colors.length - 1]); P.saveUser(); render(); } }, icon('plus'), 'Colour'));
    const nt = p.nucleotides || {};
    const ntRow = el('div', { class: 'pal-colors' }, ...['A', 'C', 'G', 'U'].map(b => {
      const inp = el('input', { type: 'color', value: nt[b] || '#888888', disabled: userPal ? null : true });
      inp.addEventListener('input', () => { p.nucleotides = { ...(p.nucleotides || {}), [b]: inp.value }; P.saveUser(); if (P.activeId === p.id) P.apply(p.id); });
      return el('label', { class: 'pal-nt' }, el('b', { text: b }), inp);
    }));
    const nameIn = el('input', { class: 'input', value: p.name, disabled: userPal ? null : true });
    nameIn.addEventListener('change', () => { p.name = nameIn.value.trim() || p.name; P.saveUser(); render(); });
    editor.replaceChildren(
      el('div', { class: 'field' }, el('label', { text: 'Name' }), nameIn),
      el('div', { class: 'flabel', text: 'Categorical colours (helices, structure elements, annotations, plots)' }), colorsRow,
      el('div', { class: 'flabel', style: { marginTop: '8px' }, text: 'Nucleotide colours' }), ntRow,
      el('div', { class: 'pal-check' }, fmtChk('normal', 'Normal vision'), fmtChk('deutan', 'Deuteranopia'), fmtChk('protan', 'Protanopia'), fmtChk('tritan', 'Tritanopia')),
      el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '10px' } },
        el('button', { class: 'btn small primary', onclick: () => { P.apply(p.id); render(); toast(`${p.name} is now used for all colour schemes.`, 'ok', 'Palette'); } }, icon('check'), 'Use this palette'),
        el('button', { class: 'btn small', onclick: () => { const c = { id: `user-${Date.now()}`, name: `${p.name} (copy)`, colors: [...p.colors], nucleotides: { ...(p.nucleotides || {}) } }; P.user.push(c); P.saveUser(); selId = c.id; render(); } }, icon('copy'), 'Duplicate to edit'),
        userPal ? el('button', { class: 'btn small ghost', onclick: () => { P.user = P.user.filter(x => x !== p); P.saveUser(); if (P.activeId === p.id) P.apply('rcv'); selId = P.activeId; render(); } }, icon('trash'), 'Delete') : null,
        el('button', { class: 'btn small ghost', onclick: () => navigator.clipboard?.writeText(JSON.stringify({ name: p.name, colors: p.colors, nucleotides: p.nucleotides }, null, 1)).then(() => toast('Palette copied as JSON.', 'ok')) }, 'Copy as JSON'),
        el('button', { class: 'btn small ghost', onclick: async () => {
          try {
            const o = JSON.parse(await navigator.clipboard.readText());
            if (!Array.isArray(o.colors) || o.colors.length < 4 || !o.colors.every(c => /^#[0-9a-f]{6}$/i.test(c))) throw new Error('needs at least 4 colours as #rrggbb');
            const c = { id: `user-${Date.now()}`, name: String(o.name || 'Imported palette').slice(0, 40), colors: o.colors.slice(0, 12), nucleotides: o.nucleotides };
            P.user.push(c); P.saveUser(); selId = c.id; render();
          } catch (err) { toast(`Copy a palette (JSON) first: ${err.message}`, 'error', 'Import'); }
        } }, 'Paste palette')),
      userPal ? null : el('div', { class: 'desc', text: 'Built-in palettes cannot be changed: duplicate one to make your own.' }));
  };
  render();
  dialog({ title: 'Colour palettes', width: 760, body: b => b.append(el('div', { class: 'pal-grid' }, listBox, editor)), actions: [{ label: 'Close', primary: true, close: true }] });
}

export { rgbToHex };
