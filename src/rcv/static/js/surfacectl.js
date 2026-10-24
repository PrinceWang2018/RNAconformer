// Surface and transparency controls: the molecular surface (whole molecule or a selection) follows the shown
// conformation (rebuilt after every morph, every frame in video export) and takes the residue colours; each kind of
// object (cartoon, bases, atoms, surface, ensemble overlay) has its own transparency.
import { el } from './util.js';
import { dialog } from './ui.js';
import { gaussianSurface } from './viewer/surface.js';

export class SurfaceController {
  constructor(app) {
    this.app = app;
    this.enabled = false;
    this.residues = null;            // null = whole molecule
    this.spacing = 0.9;
    app.on('morph', m => {
      if (!this.enabled) return;
      if (m.done) this.update(); else if (!app.capturing) app.viewer.setSurfaceVisible(false);
    });
    app.on('ensemble', () => { this.residues = null; this.update(); });
    app.registerSessionPart('surface', {
      save: () => this.state(),
      load: (ens, s) => this.restore(s),
    });
  }

  state() { return { enabled: this.enabled, residues: this.residues ? [...this.residues] : null, alpha: { ...(this.app.viewer.settings.alpha || {}) } }; }
  restore(s) {
    if (!s) { this.set(false); return; }
    if (s.alpha) for (const [k, v] of Object.entries(s.alpha)) this.app.viewer.setTransparency(k, v);
    this.residues = s.residues ? new Set(s.residues) : null;
    this.set(!!s.enabled, { quiet: true });
  }

  set(on, { residues, quiet = false } = {}) {
    this.enabled = on;
    if (residues !== undefined) this.residues = residues ? new Set(residues) : null;
    this.update();
    this.app._syncToolbar?.();
    if (!quiet) this.app.emit('surface');
  }

  atoms() {
    const ens = this.app.ens, out = [];
    const hidden = this.app.displayOf().hidden;
    for (let r = 0; r < ens.nRes; r++) {
      if (this.residues && !this.residues.has(r)) continue;
      if (!this.residues && hidden.has(r)) continue;
      const res = ens.residues[r];
      for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) out.push(a);
    }
    return out;
  }

  update(coords = this.app.display, { spacing = this.spacing } = {}) {
    const app = this.app, v = app.viewer;
    if (!this.enabled || !app.ens || !coords) { v.setSurface(null); return; }
    const t0 = performance.now();
    const data = gaussianSurface(app.ens, coords, this.atoms(), { spacing });
    v.setSurface(data);
    this.lastMs = performance.now() - t0;
  }
}

export function transparencyDialog(app) {
  const v = app.viewer;
  const A = { cartoon: 1, bases: 1, atoms: 1, surface: 1, overlay: 0.2, ...(v.settings.alpha || {}) };
  const row = (key, label) => {
    const r = el('input', { type: 'range', min: 5, max: 100, step: 5, value: Math.round(A[key] * 100), 'aria-label': `${label} opacity` });
    const val = el('span', { class: 'tr-val', text: `${Math.round(A[key] * 100)} %` });
    r.addEventListener('input', () => { app.setTransparency(key, +r.value / 100); val.textContent = `${r.value} %`; });
    return el('div', { class: 'tr-row' }, el('span', { class: 'tr-lab', text: label }), r, val);
  };
  const surf = el('input', { type: 'checkbox', checked: app.surface.enabled ? true : null });
  surf.addEventListener('change', () => app.surface.set(surf.checked));
  const selOnly = el('input', { type: 'checkbox', checked: app.surface.residues ? true : null, disabled: app.state.selection.size ? null : true });
  selOnly.addEventListener('change', () => app.surface.set(true, { residues: selOnly.checked ? [...app.state.selection] : null }));
  dialog({
    title: 'Surface and transparency', width: 460,
    body: b => b.append(
      el('label', { class: 'check' }, surf, 'Molecular surface (Gaussian, coloured like the residues)'),
      el('label', { class: 'check' }, selOnly, `Only around the selection (${app.state.selection.size} residues)`),
      el('div', { class: 'section-title', style: { marginTop: '10px' }, text: 'Opacity' }),
      row('cartoon', 'Cartoon'), row('bases', 'Bases'), row('atoms', 'Atoms and bonds'), row('surface', 'Surface'), row('overlay', 'Ensemble overlay'),
      el('div', { class: 'desc', text: 'Console: show surface[, selection] · hide surface · set transparency, 0.5, surface. A see-through surface over the cartoon is a classic figure style: surface 40–60 %.' })),
    actions: [{ label: 'Close', primary: true, close: true }],
  });
}
