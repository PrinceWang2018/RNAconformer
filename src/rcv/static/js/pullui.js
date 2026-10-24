// Pull-to-switch interaction: drag a residue (Pull tool, or Alt+drag), the structure stretches elastically,
// the nearest stable conformation is previewed live, and on release the structure switches to it. If no
// computed conformation matches the pull well, the engine can sample new stable conformations near it.
import { el, fmtVal, fmtSigned } from './util.js';
import { icon } from './icons.js';
import { toast } from './ui.js';
import { ElasticPull, rankConformations, chooseStable, fittedCoords } from './pull.js';
import { direction } from './panels/energy.js';

const GOOD_MATCH = 3.0;          // Å: a computed conformation realises the pull this well -> just switch

export class PullController {
  constructor(app) {
    this.app = app;
    this.mode = false;
    this.active = null;
    this.hud = el('div', { class: 'pull-hud', hidden: true });
    this.card = el('div', { class: 'pull-card', hidden: true, role: 'status' });
    app.viewer.container.append(this.hud, this.card);
    app.viewer.controls.interceptor = e => this.onPointerDown(e);
  }

  setMode(on) {
    this.mode = on;
    this.app.viewer.container.classList.toggle('pull-mode', on);
    this.app.emit('pullmode', on);
  }

  onPointerDown(e) {
    const app = this.app;
    if (e.button !== 0 || !app.ens || !(this.mode || e.altKey)) return false;
    const p = app.viewer.pick(e.clientX, e.clientY);
    if (!p) return false;
    this.begin(p.residue, e);
    return true;
  }

  begin(residue, e) {
    const app = this.app, v = app.viewer, ens = app.ens;
    app.stop();
    if (app.morph) { app.morph = null; }
    this.hideCard();
    const k = app.state.conf;
    const pull = new ElasticPull(ens, app.display, residue, { pairs: ens.conformations[k].pairs });
    const a = ens.trace[residue];
    this.active = { pull, residue, start: k, anchor: [app.display[3 * a], app.display[3 * a + 1], app.display[3 * a + 2]],
      startDisplay: app.display.slice(), last: 0, choice: null, pointerId: e.pointerId, viewDir: v.viewDirModel(), previewK: null };
    try { v.canvas.setPointerCapture(e.pointerId); } catch { /* synthetic events (tests, plugins) */ }
    v.container.classList.add('pulling');
    this._move = ev => this.move(ev);
    this._up = ev => this.end(ev, false);
    this._key = ev => { if (ev.key === 'Escape') this.end(ev, true); };
    v.canvas.addEventListener('pointermove', this._move);
    v.canvas.addEventListener('pointerup', this._up);
    v.canvas.addEventListener('pointercancel', this._up);
    window.addEventListener('keydown', this._key);
    this.hud.hidden = false;
    this._hud(e, `Pulling ${ens.resLabel(residue)}`, 'move the mouse; release to switch · Esc cancels');
  }

  move(e) {
    const A = this.active;
    if (!A) return;
    const app = this.app, v = app.viewer;
    const T = v.screenToModel(e.clientX, e.clientY, A.anchor);
    A.pull.setTarget(T);
    A.pull.solve(25);
    A.pull.apply(app.display);
    v.setCoords(app.display);
    const now = performance.now();
    if (now - A.last > 90) {
      A.last = now;
      this._evaluate();
    }
    this._hudPos(e);
  }

  _evaluate() {
    const A = this.active, app = this.app, ens = app.ens;
    const ranked = rankConformations(ens, A.pull, { viewDir: A.viewDir });
    const sel = chooseStable(ens, ranked, app.state.metric, A.start);
    A.ranked = ranked;
    A.choice = sel.choice;
    const c = sel.choice;
    const pulled = Math.max(...A.pull.displacement());
    if (c.k !== A.start && pulled > 1.0) {
      if (A.previewK !== c.k) { A.previewK = c.k; app.viewer.setPreview(fittedCoords(ens, c.k, c.fit)); }
      app.emit('morph', { from: A.start, to: c.k, t: 1 });
    } else {
      if (A.previewK !== null) { A.previewK = null; app.viewer.setPreview(null); }
      app.emit('morph', { done: true });
    }
    const m = ens.metricByKey[app.state.metric] || {};
    const dv = ens.value(c.k, app.state.metric) !== null && ens.value(A.start, app.state.metric) !== null ? ens.value(c.k, app.state.metric) - ens.value(A.start, app.state.metric) : null;
    const target = c.k === A.start ? 'stays in the current conformation' : `→ ${ens.shortId(c.k)}`;
    this.hudText = [`${ens.resLabel(A.residue)} pulled ${pulled.toFixed(1)} Å ${target}`,
      `match ${c.local.toFixed(1)} Å${dv !== null && c.k !== A.start ? ` · ${m.label || ''} ${fmtSigned(dv)}` : ''}${sel.near.length > 1 ? ` · ${sel.near.length} close candidates` : ''}`];
    this.hudDir = direction(m, dv);
    this._hud(null, ...this.hudText);
  }

  end(e, cancel) {
    const A = this.active;
    if (!A) return;
    const app = this.app, v = app.viewer, ens = app.ens;
    try { v.canvas.releasePointerCapture(A.pointerId); } catch { /* released */ }
    v.canvas.removeEventListener('pointermove', this._move);
    v.canvas.removeEventListener('pointerup', this._up);
    v.canvas.removeEventListener('pointercancel', this._up);
    window.removeEventListener('keydown', this._key);
    v.container.classList.remove('pulling');
    this.hud.hidden = true;
    v.setPreview(null);
    this.active = null;
    if (cancel) { app.goTo(A.start, { force: true, duration: 0.4, target: A.startDisplay }); return; }
    this.active = A;
    this._evaluate();
    this.active = null;
    const c = A.choice;
    const pulled = Math.max(...A.pull.displacement());
    if (pulled < 1.0) { app.goTo(A.start, { force: true, duration: 0.25, target: A.startDisplay }); return; }
    const deformed = A.pull.trace();
    if (c.k === A.start) app.goTo(A.start, { force: true, duration: 0.6, target: A.startDisplay });
    else app.goTo(c.k, { target: fittedCoords(ens, c.k, c.fit) });
    this.last = { residue: A.residue, start: A.start, choice: c, deformed, pull: A.pull, doc: ens.serverId, ens, viewDir: A.viewDir };
    this._card(c, A, pulled);
  }

  // ---------------------------------------------------------------- feedback
  _hud(e, line1, line2) {
    this.hud.replaceChildren(el('b', { text: line1 }), line2 ? el('span', { class: `pull-sub ${this.hudDir || ''}`, text: line2 }) : null);
    if (e) this._hudPos(e);
  }
  _hudPos(e) {
    const r = this.app.viewer.container.getBoundingClientRect();
    this.hud.style.left = `${Math.min(r.width - 280, e.clientX - r.left + 18)}px`;
    this.hud.style.top = `${Math.max(8, e.clientY - r.top - 54)}px`;
  }

  hideCard() { this.card.hidden = true; clearTimeout(this._cardTimer); }

  _card(c, A, pulled) {
    const app = this.app, ens = app.ens;
    const good = c.local <= GOOD_MATCH;
    const stayed = c.k === A.start;
    const m = ens.metricByKey[app.state.metric] || {};
    const dv = !stayed && ens.value(c.k, app.state.metric) !== null && ens.value(A.start, app.state.metric) !== null
      ? ens.value(c.k, app.state.metric) - ens.value(A.start, app.state.metric) : null;
    const title = stayed ? 'No other computed conformation fits this pull better'
      : good ? `Switched to ${ens.shortId(c.k)}, the nearest stable conformation`
        : `Switched to ${ens.shortId(c.k)}, the closest computed conformation`;
    const detail = stayed
      ? `Your pull of ${ens.resLabel(A.residue)} (${pulled.toFixed(1)} Å) is best matched by the current conformation.`
      : `${ens.resLabel(A.residue)} and its neighbours lie within ${c.local.toFixed(1)} Å of where you pulled them` +
        (dv !== null ? `; ${m.label} ${fmtSigned(dv)}${m.unit ? ' ' + m.unit : ''}.` : '.') +
        (good ? '' : ' None of the computed conformations reaches the pull closely.');
    const explore = el('button', { class: `btn small ${good && !stayed ? '' : 'primary'}`, title: 'Run the prediction engine with your pull as a temporary restraint, then relax without it',
      onclick: () => this.explore() }, icon('flask'), 'Explore new conformations');
    const undo = el('button', { class: 'btn small ghost', onclick: () => { this.hideCard(); app.goTo(A.start); } }, 'Undo');
    const close = el('button', { class: 'btn ghost icon-only small', 'aria-label': 'Close', onclick: () => this.hideCard() }, icon('x'));
    this.card.replaceChildren(el('div', { class: 'pc-icon' }, icon(stayed || !good ? 'info' : 'check')),
      el('div', { class: 'pc-body' }, el('b', { text: title }), el('span', { text: detail })),
      el('div', { class: 'pc-actions' }, stayed ? null : undo, app.api.online ? explore : null, close));
    this.card.hidden = false;
    clearTimeout(this._cardTimer);
    if (good && !stayed) this._cardTimer = setTimeout(() => this.hideCard(), 12000);
  }

  async explore() {
    const app = this.app, L = this.last;
    if (!L || !L.doc) { toast('Exploration needs an ensemble opened through the engine.', 'error'); return; }
    this.hideCard();
    let job;
    try {
      job = await app.api.post('/api/explore', { document: L.doc, conformation: L.start, residue: L.residue, deformed: L.deformed,
        cores: app.info?.system?.cpus ? Math.max(1, app.info.system.cpus - 1) : undefined });
    } catch (err) { toast(err.message, 'error', 'Cannot explore'); return; }
    toast(`Sampling new conformations near your pull of ${L.ens.resLabel(L.residue)} on this computer (about ${Math.round((job.eta_s || 120) / 60)} min). You can keep working.`, 'info', 'Exploring', 8000);
    app.jobsDock.track(job, { title: `Explore pull · ${L.ens.resLabel(L.residue)}`, onDone: async j => {
      const ens = await app.replaceDocument(L.ens.id, j.document.id);
      const added = new Set(j.result.added || []);
      const ranked = rankConformations(ens, L.pull, { viewDir: L.viewDir }).filter(r => added.has(ens.conformations[r.k].id));
      if (ranked.length) {
        const best = chooseStable(ens, ranked, app.state.metric).choice;
        app.goTo(best.k);
        toast(`${added.size} new conformation${added.size === 1 ? '' : 's'} added. Showing ${ens.shortId(best.k)} (match ${best.local.toFixed(1)} Å to your pull).`, 'ok', 'Exploration finished', 9000);
      } else toast('The exploration finished but produced no valid conformation.', 'error');
    } });
  }
}
