// Undo / redo of everything the user edits on an ensemble: custom colours (residues and atoms), hidden or highlighted
// residues, labels, annotations, measurements, scenes, named selections, the selection, colour scheme,
// representation, style, background, overlay, putty, surface, transparency and plugin state (e.g. tracked distances).
// Viewing (camera, conformation) is not undone. One history per ensemble; changes are grouped (a drag or a script is
// one step).
import { toast } from './ui.js';

const MAX = 150;
const EVENTS = ['display', 'selection', 'settings', 'measurements', 'annotations', 'scenes', 'named', 'surface', 'plugin-state'];
const NOT_UNDONE = new Set(['commands', 'conditions', 'probing']);   // settings and loaded data, not edits

export class History {
  constructor(app) {
    this.app = app;
    this.stacks = new Map();          // ensemble id -> { undo: [], redo: [], current }
    this.restoring = false;
    this.timer = 0;
    for (const ev of EVENTS) app.on(ev, () => this.schedule());
    app.on('ensemble', ens => { if (ens && !this.stacks.has(ens.id)) this.stacks.set(ens.id, { undo: [], redo: [], current: this.capture() }); this.sync(); });
  }

  get stack() { const ens = this.app.ens; return ens ? this.stacks.get(ens.id) : null; }

  capture() {
    const app = this.app, ens = app.ens;
    if (!ens) return null;
    const parts = {};
    for (const [id, p] of app._sessionParts || []) {
      if (NOT_UNDONE.has(id)) continue;
      try { parts[id] = p.save(ens); } catch { /* ignore */ }
    }
    const vs = app.viewer.settings;
    return JSON.stringify({
      parts, selection: [...app.state.selection].sort((a, b) => a - b), scheme: app.state.colorScheme, uniform: app.state.uniformColor,
      rep: vs.representation, style: vs.style, bg: vs.background, ghosts: app.state.ghosts, putty: app.state.putty,
      alpha: vs.alpha ? { ...vs.alpha } : null, surface: app.surface?.state?.() || null, labelMode: app.viewer.labelMode,
    });
  }

  schedule() {
    if (this.restoring || !this.app.ens) return;
    clearTimeout(this.timer);
    this.timer = setTimeout(() => this.commit(), 260);
  }

  commit() {
    clearTimeout(this.timer);
    const S = this.stack;
    if (!S || this.restoring) return;
    const now = this.capture();
    if (now === S.current) return;
    if (S.current !== null) S.undo.push(S.current);
    if (S.undo.length > MAX) S.undo.shift();
    S.current = now;
    S.redo = [];
    this.sync();
  }

  undo() {
    this.commit();
    const S = this.stack;
    if (!S || !S.undo.length) { toast('Nothing to undo.', 'info', '', 1500); return; }
    S.redo.push(S.current);
    S.current = S.undo.pop();
    this.restore(S.current);
  }

  redo() {
    this.commit();
    const S = this.stack;
    if (!S || !S.redo.length) { toast('Nothing to redo.', 'info', '', 1500); return; }
    S.undo.push(S.current);
    S.current = S.redo.pop();
    this.restore(S.current);
  }

  restore(json) {
    const app = this.app, ens = app.ens;
    if (!ens || !json) return;
    const s = JSON.parse(json);
    this.restoring = true;
    try {
      for (const [id, p] of app._sessionParts || []) {
        if (NOT_UNDONE.has(id)) continue;
        try { p.load(ens, s.parts?.[id] ?? null); } catch (err) { console.error(err); }
      }
      if (s.rep !== app.viewer.settings.representation) app.setRepresentation(s.rep);
      if (s.style !== app.viewer.settings.style) app.setStyle(s.style);
      if (s.bg !== app.viewer.settings.background) app.setBackground(s.bg);
      if (s.uniform) app.state.uniformColor = s.uniform;
      if (s.scheme !== app.state.colorScheme) app.setColorScheme(s.scheme);
      if (!!s.ghosts !== app.state.ghosts) app.toggleGhosts();
      if (!!s.putty !== app.state.putty) app.togglePutty();
      if (s.alpha && app.setTransparency) for (const [k, v] of Object.entries(s.alpha)) app.setTransparency(k, v, { quiet: true });
      if (app.surface) app.surface.restore(s.surface);
      if (s.labelMode) app.viewer.labelMode = s.labelMode;
      app.select(s.selection || [], 'set');
      app.applyDisplay();
      app.viewer.needsRender = true;
    } finally {
      clearTimeout(this.timer);
      this.restoring = false;
    }
    this.sync();
  }

  sync() {
    const S = this.stack;
    const u = document.getElementById('btn-undo'), r = document.getElementById('btn-redo');
    if (u) u.disabled = !S || !S.undo.length;
    if (r) r.disabled = !S || !S.redo.length;
    this.app.rebuildMenus?.();
  }

  get canUndo() { return !!this.stack?.undo.length; }
  get canRedo() { return !!this.stack?.redo.length; }
}

// Clear menu actions (each one is a single undoable step).
export function clearItems(app) {
  const has = !!app.ens;
  const d = () => app.displayOf();
  const clearDisplay = () => { const x = d(); x.colors.clear(); x.atomColors.clear(); x.hidden.clear(); x.sticks.clear(); x.spheres.clear(); x.labels.clear(); app.applyDisplay(); };
  const all = () => { clearDisplay(); app.annotations.clear(); app.measure.clear(); app.select([], 'set'); };
  return [
    { label: 'Clear all edits (colours, display, labels, annotations, measurements, selection)', action: () => { all(); toast('Cleared. Undo (Ctrl+Z) brings it back.', 'info', 'Clear', 3000); }, disabled: !has },
    { sep: true },
    { label: 'Custom colours (residues and atoms)', action: () => { d().colors.clear(); d().atomColors.clear(); app.applyDisplay(); }, disabled: !has },
    { label: 'Display changes (hidden, sticks, spheres)', action: () => { const x = d(); x.hidden.clear(); x.sticks.clear(); x.spheres.clear(); app.applyDisplay(); }, disabled: !has },
    { label: 'Residue labels', action: () => { d().labels.clear(); app.applyDisplay(); }, disabled: !has },
    { label: 'Annotations', action: () => app.annotations.clear(), disabled: !has },
    { label: 'Measurements', action: () => app.measure.clear(), disabled: !has },
    { label: 'Selection', kbd: 'Esc', action: () => app.select([], 'set'), disabled: !has },
  ];
}
