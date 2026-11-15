// Window layout of the main window: every part can be shown or hidden (Window menu), the side panels, the 2D
// structure window and the console can be maximised to fill the window (button in their title bar, double-click
// on the title bar, Esc to restore), H hides every panel around the 3D view, and F11 switches to full screen.
// The choice of shown parts is remembered between sessions.
//
// 0.7.5: the side panels are resized by dragging their inner edge (double-click: default width); the widths are
// remembered for each screen size. A hidden side panel folds to a narrow rail of icons that opens it again with
// one click. Presentation mode (F5) shows only the 3D view, without the menu bar, with large labels, full screen.
import { el, $, loadPref, savePref } from './util.js';
import { icon } from './icons.js';
import { toast } from './ui.js';

export const WINDOWS = [
  { id: 'left', label: 'Ensembles and conformations panel', max: true },
  { id: 'right', label: 'Analysis panel (energy, details, switch)', max: true },
  { id: 'ss2d', label: '2D structure window', kbd: 'D', max: true },
  { id: 'console', label: 'Command console', kbd: '`', max: true },
  { id: 'seqbar', label: 'Sequence bar' },
  { id: 'toolbar', label: 'Toolbar' },
  { id: 'timeline', label: 'Playback bar' },
  { id: 'statusbar', label: 'Status bar' },
  { id: 'hud', label: 'Title and legend on the 3D view' },
];
const SIMPLE = ['toolbar', 'timeline', 'statusbar', 'hud'];

// side panel widths (px): smallest, largest, and room always kept for the 3D view
export const WIDTH_LIMITS = { left: [190, 520], right: [312, 760], stage: 360 };
// icons of the right panel's tabs on its rail (plugin tabs: puzzle)
const TAB_ICONS = { energy: 'chart', switch: 'swap', pathway: 'target', details: 'info' };
const screenKey = () => `${window.screen?.width || 0}x${window.screen?.height || 0}`;

export class WindowManager {
  constructor(app) {
    this.app = app;
    const p = loadPref('windows', {});
    this.shown = { left: true, right: true, toolbar: true, timeline: true, statusbar: true, hud: true, ...(p.shown || {}) };
    this.widths = p.widths || {};  // screen size ("2560x1440") -> { left, right } in px
    this.maxed = null;            // id of the maximised window, or null
    this.focus = false;           // H: only the 3D view (and the menu bar)
    this.present = false;         // F5: presentation mode
    this.root = $('#app');
    for (const id of ['left', 'right', ...SIMPLE]) this.root.classList.toggle(`no-${id}`, !this.shown[id]);
    this._buildSplitters();
    this._buildRails();
    this._applyWidths();
    document.addEventListener('fullscreenchange', () => {
      // Esc in a browser leaves full screen before RCV sees the key: leave presentation mode with it
      if (this.present && this._presentFull && !this.fullscreen) { this._presentFull = false; this.setPresent(false); }
      this._syncFullscreen();
    });
    window.addEventListener('resize', () => { this._syncFullscreen(); this._applyWidths(); });
  }

  _save() { savePref('windows', { shown: this.shown, widths: this.widths }); }
  _narrow() { return window.innerWidth <= 820; }

  isShown(id) {
    const a = this.app;
    switch (id) {
      case 'left': case 'right': return this._narrow() ? this.root.classList.contains(`show-${id}`) : !!this.shown[id];
      case 'ss2d': return !!a.ss2d?.visible;
      case 'console': return !!a.console?.visible;
      case 'seqbar': return !!a.seqbar?.visible;
      default: return !!this.shown[id];
    }
  }

  set(id, on) {
    const a = this.app;
    if (!on && this.maxed === id) this.restore();
    switch (id) {
      case 'left': case 'right':
        if (this._narrow()) this.root.classList.toggle(`show-${id}`, on);
        else { this.shown[id] = on; this.root.classList.toggle(`no-${id}`, !on); this._save(); }
        break;
      case 'ss2d': a.ss2d?.setVisible(on); break;
      case 'console': a.console?.setVisible(on); break;
      case 'seqbar': a.seqbar?.setVisible(on); break;
      default: this.shown[id] = on; this.root.classList.toggle(`no-${id}`, !on); this._save();
    }
    if (on && this.focus) this.setFocus(false);
    if (on && this.present) this.setPresent(false);
    this._changed();
  }
  toggle(id) { this.set(id, !this.isShown(id)); }

  showAll() {
    for (const w of WINDOWS) if (w.id !== 'console' && !this.isShown(w.id)) this.set(w.id, true);
    this.setFocus(false);
  }
  resetLayout() {
    this.restore();
    this.setFocus(false);
    for (const w of WINDOWS) this.set(w.id, w.id !== 'console');
    delete this.widths[screenKey()];
    this._applyWidths();
    this._save();
    this.app.ss2d?.resetPlace?.();
    toast('Window layout reset.', 'ok', '', 2000);
  }

  // ---------------------------------------------------------------- maximise one window
  canMaximize(id) { return !!WINDOWS.find(w => w.id === id)?.max; }
  maximize(id) {
    if (this.maxed === id) return;
    this.restore(false);
    if (!this.isShown(id)) this.set(id, true);
    this.maxed = id;
    if (id === 'ss2d') this.app.ss2d?.setMaximized(true);
    else if (id === 'console') this.app.console?.root.classList.add('maxed');
    else this.root.dataset.max = id;
    this._changed();
  }
  restore(notify = true) {
    const id = this.maxed;
    if (!id) return false;
    this.maxed = null;
    if (id === 'ss2d') this.app.ss2d?.setMaximized(false);
    else if (id === 'console') this.app.console?.root.classList.remove('maxed');
    delete this.root.dataset.max;
    if (notify) this._changed();
    return true;
  }
  toggleMax(id) { if (this.maxed === id) this.restore(); else this.maximize(id); }

  // ---------------------------------------------------------------- 3D view only (H)
  setFocus(on) {
    if (on === this.focus || (on && this.present)) return;
    this.focus = on;
    if (on) this.restore(false);
    this.root.classList.toggle('focus3d', on);
    if (on) toast('Panels hidden: only the 3D view is shown. Press H or Esc to show them again.', 'info', '', 3000);
    this._changed();
  }

  // ---------------------------------------------------------------- presentation mode (F5)
  // Only the 3D view: no menu bar, panels or bars, large labels, full screen. Keys keep working (arrows, Space,
  // R, S ...); Esc or F5 leaves it, and so does the button that appears when the mouse moves.
  async setPresent(on) {
    if (on === this.present) return;
    this.present = on;
    if (on) { this.restore(false); this.setFocus(false); }
    this.root.classList.toggle('present', on);
    this._presentBtn ||= this._buildPresentButton();
    this._presentBtn.hidden = true;
    if (on && !this.fullscreen) { this._presentFull = true; await this.toggleFullscreen(); }
    else if (!on && this._presentFull) { this._presentFull = false; if (this.fullscreen) await this.toggleFullscreen(); }
    if (on) toast('Presentation mode: press Esc or F5 to leave.', 'info', '', 2500);
    this._changed();
  }
  _buildPresentButton() {
    const b = el('button', { class: 'btn present-exit', type: 'button', hidden: true, onclick: () => this.setPresent(false) },
      icon('restore'), el('span', { text: 'Exit presentation (Esc)' }));
    document.body.append(b);
    let timer = null;
    window.addEventListener('pointermove', () => {
      if (!this.present) return;
      b.hidden = false;
      clearTimeout(timer);
      timer = setTimeout(() => { if (!b.matches(':hover')) b.hidden = true; }, 2200);
    });
    return b;
  }

  // ---------------------------------------------------------------- side panel widths (drag the inner edge)
  _buildSplitters() {
    for (const side of ['left', 'right']) {
      const bar = this.root.querySelector(`.sidebar.${side}`);
      if (!bar) continue;
      const h = el('div', { class: `splitter ${side}`, role: 'separator', 'aria-orientation': 'vertical', tabindex: '-1',
        title: 'Drag to resize the panel; double-click for the default width' });
      h.addEventListener('pointerdown', e => this._dragStart(side, bar, h, e));
      h.addEventListener('dblclick', () => this.setWidth(side, null));
      bar.append(h);
    }
  }
  _dragStart(side, bar, h, e) {
    if (e.button !== 0 || this._narrow()) return;
    e.preventDefault();
    try { h.setPointerCapture(e.pointerId); } catch { /* synthetic event */ }
    this.root.classList.add('resizing');
    const r = bar.getBoundingClientRect();
    let raf = 0, last = null;
    const move = ev => {
      last = side === 'left' ? ev.clientX - r.left : r.right - ev.clientX;
      if (!raf) raf = requestAnimationFrame(() => { raf = 0; this._setVar(side, this._clamp(side, last)); });
    };
    const up = () => {
      h.removeEventListener('pointermove', move);
      h.removeEventListener('pointerup', up);
      h.removeEventListener('pointercancel', up);
      this.root.classList.remove('resizing');
      if (last !== null) this.setWidth(side, last);
    };
    h.addEventListener('pointermove', move);
    h.addEventListener('pointerup', up);
    h.addEventListener('pointercancel', up);
  }
  _clamp(side, w) {
    const [lo, hi] = WIDTH_LIMITS[side];
    const other = this.root.querySelector(`.sidebar.${side === 'left' ? 'right' : 'left'}`);
    const room = window.innerWidth - (other?.offsetWidth || 0) - WIDTH_LIMITS.stage;
    return Math.round(Math.max(lo, Math.min(hi, room, w)));
  }
  _setVar(side, w) {
    if (w == null) this.root.style.removeProperty(`--${side}-user`);
    else this.root.style.setProperty(`--${side}-user`, `${w}px`);
  }
  // width in px for this screen size, or null for the default width
  setWidth(side, w) {
    const key = screenKey(), cur = { ...(this.widths[key] || {}) };
    if (w == null) delete cur[side]; else cur[side] = this._clamp(side, w);
    if (Object.keys(cur).length) this.widths[key] = cur; else delete this.widths[key];
    this._applyWidths();
    this._save();
    this._changed();
  }
  widthOf(side) { return this.widths[screenKey()]?.[side] ?? null; }
  _applyWidths() {
    const cur = this.widths[screenKey()] || {};
    for (const side of ['left', 'right']) this._setVar(side, cur[side] == null ? null : this._clamp(side, cur[side]));
  }

  // ---------------------------------------------------------------- rails: a hidden side panel folds to icons
  _buildRails() {
    this.rails = {};
    for (const side of ['left', 'right']) {
      const rail = el('nav', { class: `rail ${side}`, 'aria-label': side === 'left' ? 'Ensembles and conformations panel' : 'Analysis panel (energy, details, switch)' });
      this.root.append(rail);
      this.rails[side] = rail;
    }
    this._fillRails();
  }
  _fillRails() {
    if (!this.rails) return;
    const btn = (ic, title, fn) => el('button', { class: 'rail-btn', type: 'button', title, 'aria-label': title, onclick: fn }, icon(ic));
    const open = (side, then) => () => { this.set(side, true); then?.(); };
    this.rails.left.replaceChildren(
      btn('next', 'Show the ensembles and conformations panel', open('left')),
      el('span', { class: 'rail-sep' }),
      btn('layers', 'Ensembles', open('left')),
      btn('grid', 'Conformations', open('left')));
    const tabs = this.app._tabs || [];
    this.rails.right.replaceChildren(
      btn('prev', 'Show the analysis panel', open('right')),
      el('span', { class: 'rail-sep' }),
      ...tabs.map(t => btn(TAB_ICONS[t.id] || 'puzzle', t.btn.textContent, open('right', () => this.app.showTab(t.id)))));
  }

  // ---------------------------------------------------------------- full screen (F11)
  get fullscreen() { return !!document.fullscreenElement || !!this._nativeFull; }
  async toggleFullscreen() {
    const bridge = window.pywebview?.api;
    try {
      if (bridge?.toggle_fullscreen) { await bridge.toggle_fullscreen(); this._nativeFull = !this._nativeFull; }
      else if (document.fullscreenElement) await document.exitFullscreen();
      else await document.documentElement.requestFullscreen({ navigationUI: 'hide' });
    } catch (err) { toast(err.message || String(err), 'error', 'Full screen'); }
    this._changed();
  }
  _syncFullscreen() { this.root.classList.toggle('is-fullscreen', this.fullscreen); this.app.rebuildMenus?.(); }

  // Esc: restore a maximised window, then leave the 3D-only view. true when something was restored.
  escape() {
    if (this.present) { this.setPresent(false); return true; }
    if (this.restore()) return true;
    if (this.focus) { this.setFocus(false); return true; }
    return false;
  }

  // ---------------------------------------------------------------- title-bar buttons
  // Fills every <span class="win-btns" data-win="id"> under root with maximise/restore and hide buttons.
  decorate(root = document) {
    for (const host of root.querySelectorAll('.win-btns[data-win]')) {
      if (host.childElementCount) continue;
      const id = host.dataset.win;
      const maxBtn = this.canMaximize(id) ? el('button', { class: 'ss-btn win-max', type: 'button', onclick: e => { e.stopPropagation(); this.toggleMax(id); } }) : null;
      const hide = el('button', { class: 'ss-btn win-hide', type: 'button', title: 'Hide (Window menu shows it again)', 'aria-label': 'Hide',
        onclick: e => { e.stopPropagation(); this.set(id, false); } }, icon('x'));
      host.append(...[maxBtn, hide].filter(Boolean));
      const head = host.closest('[data-win-head]') || host.parentElement;
      if (this.canMaximize(id) && head && !head._winDbl) {
        head._winDbl = true;
        head.addEventListener('dblclick', e => { if (!e.target.closest('button, input, select, a, .tab')) this.toggleMax(id); });
      }
    }
    this._syncButtons();
  }
  _syncButtons() {
    for (const b of document.querySelectorAll('.win-btns[data-win] .win-max')) {
      const id = b.parentElement.dataset.win, on = this.maxed === id;
      b.replaceChildren(icon(on ? 'restore' : 'maximize'));
      b.title = on ? 'Restore (Esc)' : 'Maximise (double-click the title bar)';
      b.setAttribute('aria-label', on ? 'Restore' : 'Maximise');
    }
  }

  _changed() {
    this._syncButtons();
    this._fillRails();
    this.app.viewer?.resize();
    this.app.rebuildMenus?.();
    this.app._syncToolbar?.();
    this.app.emit?.('layout', { maxed: this.maxed, focus: this.focus });
  }

  // ---------------------------------------------------------------- Window menu
  menuItems() {
    const has = !!this.app.ens;
    return [
      { head: 'Show' },
      ...WINDOWS.map(w => ({ label: w.label, kbd: w.kbd, checked: this.isShown(w.id), action: () => this.toggle(w.id),
        disabled: (w.id === 'ss2d' || w.id === 'seqbar') && !has })),
      { sep: true },
      { label: 'Show all panels', action: () => this.showAll() },
      { label: 'Only the 3D view (hide panels)', kbd: 'H', checked: this.focus, action: () => this.setFocus(!this.focus) },
      { label: 'Full screen', kbd: 'F11', checked: this.fullscreen, action: () => this.toggleFullscreen() },
      { label: 'Presentation mode', kbd: 'F5', checked: this.present, action: () => this.setPresent(!this.present), disabled: !has },
      { sep: true }, { head: 'Maximise' },
      ...WINDOWS.filter(w => w.max).map(w => ({ label: w.label, checked: this.maxed === w.id, action: () => this.toggleMax(w.id),
        disabled: w.id === 'ss2d' && !has })),
      { label: 'Restore', kbd: 'Esc', action: () => this.restore(), disabled: !this.maxed },
      { sep: true },
      { label: 'Default panel widths', action: () => { this.setWidth('left', null); this.setWidth('right', null); }, disabled: !this.widths[screenKey()] },
      { label: 'Reset window layout', action: () => this.resetLayout() },
    ];
  }
}
