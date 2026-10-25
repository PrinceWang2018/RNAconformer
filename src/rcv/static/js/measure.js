// Measurements (distance, angle, dihedral between atoms) that follow the structure while it morphs.
// Each measurement is also an ensemble metric, so the energy panel can plot it per conformation
// (profile, funnel, map) and it is exported with the energy table.
import { el, fmt } from './util.js';
import { icon } from './icons.js';

import { SLOTS } from './colors.js';
const COLORS = SLOTS;              // the active palette
const UNIT = { distance: 'Å', angle: '°', dihedral: '°' };
const N_ATOMS = { distance: 2, angle: 3, dihedral: 4 };

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(a[0], a[1], a[2]);

export function measureValue(kind, P) {
  if (kind === 'distance') return norm(sub(P[0], P[1]));
  if (kind === 'angle') {
    const u = sub(P[0], P[1]), v = sub(P[2], P[1]);
    return Math.acos(Math.max(-1, Math.min(1, dot(u, v) / ((norm(u) * norm(v)) || 1)))) * 180 / Math.PI;
  }
  const b0 = sub(P[0], P[1]), b1 = sub(P[2], P[1]), b2 = sub(P[3], P[2]);
  const l = norm(b1) || 1;
  const n1 = b1.map(x => x / l);
  const v = sub(b0, n1.map(x => x * dot(b0, n1))), w = sub(b2, n1.map(x => x * dot(b2, n1)));
  return Math.atan2(dot(cross(n1, v), w), dot(v, w)) * 180 / Math.PI;
}

export class Measurements {
  constructor(app) {
    this.app = app;
    this.counter = 0;
    app.viewer.addOverlay('measurements', ({ coords, ens }) => this._overlay(coords, ens));
    app.registerSessionPart('measurements', {
      save: ens => (ens?.measurements || []).map(m => ({ name: m.name, kind: m.kind, atoms: m.atoms, color: m.color })),
      load: (ens, list) => { if (ens.measurements?.length) this.clear(); ens.measurements = []; for (const m of list || []) this.add(m.kind, m.atoms, { name: m.name, color: m.color, quiet: true }); },
    });
    app.on('ensemble', () => this._syncMetrics());
  }

  list(ens = this.app.ens) { return ens ? (ens.measurements ||= []) : []; }

  point(ens, coords, a) { return [coords[3 * a], coords[3 * a + 1], coords[3 * a + 2]]; }

  atomLabel(ens, a) { return `${ens.resLabel(ens.atomRes[a])}/${ens.atomNames[a]}`; }

  add(kind, atoms, { name = null, color = null, quiet = false } = {}) {
    const app = this.app, ens = app.ens;
    if (!ens) throw new Error('no ensemble open');
    if (atoms.length !== N_ATOMS[kind]) throw new Error(`${kind} needs ${N_ATOMS[kind]} atoms`);
    const list = this.list(ens);
    const id = ++this.counter;
    const m = { id, kind, atoms: atoms.slice(), name: name || `${kind[0]}${list.length + 1}`, color: color || COLORS[(list.length + 1) % COLORS.length], ens: ens.id };
    m.key = `measure.${m.name}`;
    const existing = list.findIndex(x => x.name === m.name);
    if (existing >= 0) list.splice(existing, 1);
    list.push(m);
    this._register(ens, m);
    app.viewer.needsRender = true;
    app.emit('measurements');
    if (!quiet) app.emit('values');
    return m;
  }

  remove(name) {
    const list = this.list();
    const i = list.findIndex(m => m.name === name);
    if (i < 0) return false;
    const [m] = list.splice(i, 1);
    this.app.dropMetric(this.app.ens, m.key);
    this.app.viewer.needsRender = true;
    this.app.emit('measurements');
    this.app.emit('values');
    return true;
  }

  clear() { for (const m of [...this.list()]) this.remove(m.name); }

  value(m, coords, ens = this.app.ens) { return measureValue(m.kind, m.atoms.map(a => this.point(ens, coords, a))); }

  describe(m, ens = this.app.ens) { return `${m.name}: ${m.kind} ${m.atoms.map(a => this.atomLabel(ens, a)).join(' – ')}`; }

  _register(ens, m) {
    const label = `${m.kind === 'distance' ? 'Distance' : m.kind === 'angle' ? 'Angle' : 'Dihedral'} ${m.name} (${m.atoms.map(a => this.atomLabel(ens, a)).join('–')})`;
    for (let k = 0; k < ens.n; k++) ens.conformations[k].values[m.key] = this.value(m, ens.coords[k], ens);
    const def = { key: m.key, label, unit: UNIT[m.kind], group: 'measurement', lower_is_better: null, description: `Measured in each conformation (${m.kind} between atoms).` };
    if (!ens.metricByKey[m.key]) ens.metrics.push(def); else Object.assign(ens.metricByKey[m.key], def);
    ens.metricByKey[m.key] = ens.metrics.find(x => x.key === m.key);
  }

  _syncMetrics() {
    const ens = this.app.ens;
    if (!ens) return;
    for (const m of this.list(ens)) this._register(ens, m);
  }

  _overlay(coords, ens) {
    const out = [];
    for (const m of this.list(ens)) {
      const P = m.atoms.map(a => this.point(ens, coords, a));
      const v = measureValue(m.kind, P);
      const lab = `${m.name} ${fmt(v, 1)} ${UNIT[m.kind]}`;
      if (m.kind === 'distance') out.push({ type: 'line', a: P[0], b: P[1], color: m.color, dash: true, width: 2, label: lab });
      else {
        for (let i = 0; i + 1 < P.length; i++) out.push({ type: 'line', a: P[i], b: P[i + 1], color: m.color, dash: true, width: 1.6, label: i === Math.floor((P.length - 2) / 2) ? lab : null });
      }
      for (const p of P) out.push({ type: 'point', p, color: m.color, r: 3.5 });
    }
    return out;
  }

  // Panel section (Details tab): every measurement with its value now and its range over the ensemble.
  renderSection() {
    const app = this.app, ens = app.ens;
    const list = this.list(ens);
    const box = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Measurements', el('span', { class: 'grow' }),
      el('button', { class: 'btn small ghost', title: 'Click two atoms in the 3D view to measure a distance (M)', onclick: () => app.toggleMeasureMode() }, icon('ruler'), app.measureMode ? 'Measuring…' : 'Measure')));
    if (!list.length) {
      box.append(el('div', { class: 'desc', text: 'Measure mode (M): click two atoms. Or type in the console, e.g. distance 12/N1, 34/N3 · angle 5/P, 6/P, 7/P · dihedral … Values follow every morph and can be plotted per conformation.' }));
      return box;
    }
    const tb = el('tbody');
    for (const m of list) {
      const vals = ens.conformations.map(c => c.values[m.key]).filter(v => typeof v === 'number');
      const cur = this.value(m, app.display);
      tb.append(el('tr', {},
        el('td', {}, el('span', { class: 'sw', style: { background: m.color, display: 'inline-block', width: '10px', height: '10px', borderRadius: '2px', marginRight: '6px' } }), m.name),
        el('td', { class: 'muted', text: m.atoms.map(a => this.atomLabel(ens, a)).join(' – ') }),
        el('td', { class: 'num', text: `${fmt(cur, 1)} ${UNIT[m.kind]}` }),
        el('td', { class: 'num muted', text: vals.length ? `${fmt(Math.min(...vals), 1)}–${fmt(Math.max(...vals), 1)}` : '' }),
        el('td', {}, el('button', { class: 'btn small ghost icon-only', title: 'Plot per conformation', onclick: () => { app.setMetric(m.key); app.showTab('energy'); } }, icon('image')),
          el('button', { class: 'btn small ghost icon-only', title: 'Delete', onclick: () => this.remove(m.name) }, icon('trash')))));
    }
    box.append(el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'Name' }), el('th', { text: 'Atoms' }), el('th', { class: 'num', text: 'Now' }), el('th', { class: 'num', text: 'Ensemble' }), el('th', { text: '' }))), tb));
    return box;
  }
}
