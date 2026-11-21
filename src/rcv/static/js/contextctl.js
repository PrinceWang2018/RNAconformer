// Bound molecules: what the opened structure holds besides the RNA chain that is shown (ligands, ions, protein, DNA
// and other RNA chains; document.context). Drawn by viewer/context.js; this controller keeps the drawing on the
// shown conformation, morphs it with the RNA, and holds the one switch of the View menu.
import { el, loadPref, savePref } from './util.js';
import { tipRow } from './ui.js';

export class ContextController {
  constructor(app) {
    this.app = app;
    this.enabled = loadPref('context', true) !== false;
    this.start = null;                 // shown coordinates when the running morph started
    app.viewer.setContextEnabled(this.enabled);
    app.on('ensemble', () => this.at());
    app.on('conformation', ev => { if (!ev?.scrub) this.start = app.viewer.context.snapshot(); });
    app.on('morph', m => {
      if (!this.available) return;
      if (m.done) this.at();
      else if (app.morph) app.viewer.setContextBlend(this.start, m.to, m.t);      // a morph from wherever the view was
      else app.viewer.setContextFrame(m.from, m.to, m.t);                         // scrubbing, mechanism playback
    });
  }

  get available() { return !!this.app.ens?.context; }
  get summary() { return this.app.ens?.context?.summary || ''; }

  at(k = this.app.state.conf) { if (this.available) this.app.viewer.setContextFrame(k); }
  frame(a, b, t) { if (this.available) this.app.viewer.setContextFrame(a, b, t); }

  set(on) {
    this.enabled = !!on;
    savePref('context', this.enabled);
    this.app.viewer.setContextEnabled(this.enabled);
    this.at();
    this.app.rebuildMenus?.();
    this.app.emit('context');
  }

  // Tooltip rows for what is under the pointer, when it is in front of the RNA atom found there (hit: Viewer.pick).
  tipAt(e, hit) {
    if (!this.available || !this.enabled || !e || e.buttons) return null;
    const c = this.app.viewer.pickContext(e.clientX, e.clientY);
    if (!c || (hit && hit.score <= c.score)) return null;
    return [tipRow(c.name, c.what, c.color),
      el('div', { class: 'tl', text: c.atom ? `chain ${c.chain} · atom ${c.atom}` : `chain ${c.chain}` }),
      el('div', { class: 'tl', text: 'Bound to the RNA in this structure; not part of the analyses.' })];
  }
}
