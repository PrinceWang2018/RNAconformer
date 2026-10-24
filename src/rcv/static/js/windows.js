// Window layout of the main window: every part can be shown or hidden (Window menu), the side panels, the 2D
// structure window and the console can be maximised to fill the window (button in their title bar, double-click
// on the title bar, Esc to restore), H hides every panel around the 3D view, and F11 switches to full screen.
// The choice of shown parts is remembered between sessions.
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

export class WindowManager {
  constructor(app) {
    this.app = app;
    const p = loadPref('windows', {});
    this.shown = { left: true, right: true, toolbar: true, timeline: true, statusbar: true, hud: true, ...(p.shown || {}) };
    this.maxed = null;            // id of the maximised window, or null
    this.focus = false;           // H: only the 3D view (and the menu bar)
    this.root = $('#app');
    for (const id of ['left', 'right', ...SIMPLE]) this.root.classList.toggle(`no-${id}`, !this.shown[id]);
    document.addEventListener('fullscreenchange', () => this._syncFullscreen());
    window.addEventListener('resize', () => this._syncFullscreen());
  }

  _save() { savePref('windows', { shown: this.shown }); }
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
    if (on === this.focus) return;
    this.focus = on;
    if (on) this.restore(false);
    this.root.classList.toggle('focus3d', on);
    if (on) toast('Panels hidden: only the 3D view is shown. Press H or Esc to show them again.', 'info', '', 3000);
    this._changed();
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
      { sep: true }, { head: 'Maximise' },
      ...WINDOWS.filter(w => w.max).map(w => ({ label: w.label, checked: this.maxed === w.id, action: () => this.toggleMax(w.id),
        disabled: w.id === 'ss2d' && !has })),
      { label: 'Restore', kbd: 'Esc', action: () => this.restore(), disabled: !this.maxed },
      { sep: true },
      { label: 'Reset window layout', action: () => this.resetLayout() },
    ];
  }
}
