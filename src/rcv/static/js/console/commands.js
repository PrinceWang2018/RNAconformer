// RCV command language: short typed commands for the things molecular-graphics users do all day
// (select, show/hide, colour, label, zoom, measure, superpose, save views, export), plus RNA-ensemble
// commands (conformations, morphs, transitions, switching residues). Written for RCV; the verbs are the
// common vocabulary of molecular viewers. One command per line; arguments are separated by commas.
import { evaluateSelection, resolveAtom, describeResidues, SelectionError, SELECTION_KEYWORDS } from './selection.js';
import { kabsch, applyTransform } from '../ensemble.js';
import { SCHEMES } from '../colors.js';
import { fmt, slug, saveBlob } from '../util.js';
import { conditionsOf, conditionsText, setConditions, LIGANDS, STANDARD_SALT } from '../conditions.js';
import { fmtTime } from '../panels/pathway.js';

export const COLOR_NAMES = {
  red: '#e34948', green: '#1baf7a', darkgreen: '#008300', blue: '#2a78d6', navy: '#184f95', orange: '#eb6834', yellow: '#eda100',
  gold: '#e8b800', purple: '#6250d6', violet: '#8e6bd8', magenta: '#d43d9b', pink: '#e87ba4', cyan: '#19a6c4', teal: '#14857a',
  brown: '#a1522c', white: '#ffffff', black: '#141414', gray: '#8f8d86', grey: '#8f8d86', lightgray: '#c4c2ba', lightgrey: '#c4c2ba', salmon: '#f08a6c', wheat: '#e8d3a5',
};
const REPS = { cartoon: 'cartoon', 'cartoon+atoms': 'cartoon+sticks', sticks: 'sticks', ballstick: 'ballstick', 'ball+stick': 'ballstick', spheres: 'spacefill', spacefill: 'spacefill' };

export class CommandError extends Error {}

// Split "a, b(c, d), 'e, f'" into top-level comma-separated arguments.
export function splitArgs(s) {
  const out = [];
  let depth = 0, quote = null, cur = '';
  for (const ch of String(s)) {
    if (quote) { if (ch === quote) quote = null; else cur += ch; continue; }
    if (ch === '"' || ch === "'") {
      // an apostrophe right after a letter/digit is part of a name (C4', 5')
      if (ch === "'" && /[A-Za-z0-9]$/.test(cur)) { cur += ch; continue; }
      quote = ch; continue;
    }
    if (ch === '(') depth++;
    if (ch === ')') depth--;
    if (ch === ',' && depth === 0) { out.push(cur.trim()); cur = ''; continue; }
    cur += ch;
  }
  if (cur.trim() || out.length) out.push(cur.trim());
  return out;
}

export function parseColor(s) {
  const t = String(s || '').trim().toLowerCase();
  if (COLOR_NAMES[t]) return COLOR_NAMES[t];
  if (/^#?[0-9a-f]{6}$/.test(t)) return t.startsWith('#') ? t : '#' + t;
  if (/^#?[0-9a-f]{3}$/.test(t)) { const h = t.replace('#', ''); return '#' + h.split('').map(c => c + c).join(''); }
  return null;
}

export class CommandSet {
  constructor(app) {
    this.app = app;
    this.commands = new Map();
    this._define();
  }

  // ---------------------------------------------------------------- helpers
  ctx() {
    const app = this.app;
    return { ens: app.ens, conf: app.state.conf, prev: app.state.prev, selection: app.state.selection, named: app.displayOf().named, coords: app.display };
  }
  ens() { if (!this.app.ens) throw new CommandError('no ensemble is open'); return this.app.ens; }
  sel(expr, fallback = 'all') {
    this.ens();
    try { return evaluateSelection(expr && expr.trim() ? expr : fallback, this.ctx()); } catch (e) { throw new CommandError(e.message); }
  }
  atom(spec) { this.ens(); try { return resolveAtom(spec, this.ctx()); } catch (e) { throw new CommandError(e.message); } }
  conf(spec) {
    const ens = this.ens(), app = this.app;
    const s = String(spec ?? '').trim();
    if (!s) throw new CommandError('name a conformation (id, label, or its number in the current order)');
    const low = s.toLowerCase();
    if (low === 'current' || low === 'this') return app.state.conf;
    if (low === 'previous' || low === 'prev') { if (app.state.prev === null) throw new CommandError('there is no previous conformation yet'); return app.state.prev; }
    if (low === 'first') return app.order[0];
    if (low === 'last') return app.order[app.order.length - 1];
    if (/^\d+$/.test(s)) { const i = +s - 1; if (i < 0 || i >= app.order.length) throw new CommandError(`there are ${ens.n} conformations`); return app.order[i]; }
    const k = ens.conformations.findIndex(c => c.id.toLowerCase() === low || (c.label || '').toLowerCase() === low || ens.shortId(ens.conformations.indexOf(c)).toLowerCase() === low);
    if (k < 0) throw new CommandError(`no conformation "${s}"`);
    return k;
  }
  setNamed(name, set) { this.app.displayOf().named.set(name, new Set(set)); this.app.emit('named'); }

  // ---------------------------------------------------------------- registry
  def(names, { args = '', help, run, complete = null }) {
    const list = [].concat(names);
    const cmd = { name: list[0], aliases: list.slice(1), args, help, run, complete };
    for (const n of list) this.commands.set(n, cmd);
  }
  names() { return [...new Set([...this.commands.values()].map(c => c.name))].sort(); }

  // Execute one line; returns an array of output lines ({ text, kind }).
  async execute(line) {
    const src = String(line).trim();
    if (!src || src.startsWith('#')) return [];
    const m = /^([A-Za-z_][\w+-]*)\s*(.*)$/s.exec(src);
    if (!m) throw new CommandError(`cannot read "${src}"`);
    const name = m[1].toLowerCase();
    const cmd = this.commands.get(name);
    if (!cmd) {
      const near = this.names().filter(n => n.startsWith(name.slice(0, 2)));
      throw new CommandError(`unknown command "${name}"${near.length ? ` — did you mean ${near.slice(0, 4).join(', ')}?` : ''} (type help)`);
    }
    const out = await cmd.run(splitArgs(m[2]), m[2], name);
    if (out === undefined || out === null) return [];
    return [].concat(out).map(x => (typeof x === 'string' ? { text: x, kind: 'out' } : x));
  }

  // ---------------------------------------------------------------- commands
  _define() {
    const app = this.app;
    const say = t => ({ text: t, kind: 'out' });
    const need = (args, n, usage) => { if (args.filter(a => a !== '').length < n) throw new CommandError(`usage: ${usage}`); };

    this.def(['help', '?'], { args: '[command]', help: 'List the commands, or explain one.', run: args => {
      if (args[0]) {
        const c = this.commands.get(args[0].toLowerCase());
        if (!c) throw new CommandError(`no command "${args[0]}"`);
        return [say(`${c.name} ${c.args}`), say(`  ${c.help}`), ...(c.aliases.length ? [say(`  also: ${c.aliases.join(', ')}`)] : [])];
      }
      const rows = this.names().map(n => { const c = this.commands.get(n); return say(`${(c.name + ' ' + c.args).padEnd(38)} ${c.help}`); });
      return [say('Commands (arguments separated by commas; selections like: resi 5-12, resn G and paired, helix H2, switching, rmsf > 3, within 5 of resi 20):'), ...rows,
        say('Keys: ↑/↓ history · Tab completes · Esc closes · several lines can be pasted at once (a script).')];
    } });

    // -------- selections
    this.def(['select', 'sel'], { args: '[name,] selection', help: 'Select residues (shown as atoms and labelled); with a name, also store it as a named selection.', run: args => {
      need(args, 1, 'select [name,] selection');
      let name = null, expr = args[0];
      if (args.length >= 2) { name = args[0]; expr = args.slice(1).join(','); }
      if (name && !/^[A-Za-z_][\w-]*$/.test(name)) throw new CommandError(`"${name}" is not a valid name (letters, digits, _ and -)`);
      if (name && SELECTION_KEYWORDS.includes(name)) throw new CommandError(`"${name}" is a selection keyword; choose another name`);
      const s = this.sel(expr);
      app.select([...s], 'set');
      if (name) this.setNamed(name, s);
      return say(`${name ? name + ': ' : ''}${describeResidues(app.ens, s)}`);
    } });
    this.def(['deselect'], { help: 'Clear the selection.', run: () => { app.select([], 'set'); return say('selection cleared'); } });
    this.def(['count'], { args: 'selection', help: 'Count residues matching a selection without selecting them.', run: (args, raw) => {
      const s = this.sel(raw);
      return say(describeResidues(app.ens, s));
    } });

    // -------- representation
    const showHide = show => (args, raw, verb) => {
      const ens = this.ens();
      const d = app.displayOf();
      const what = (args[0] || '').toLowerCase();
      const s = args.length > 1 ? this.sel(args.slice(1).join(',')) : null;
      const all = () => new Set([...Array(ens.nRes).keys()]);
      const target = s || all();
      if (!what || what === 'everything' || what === 'all') {
        if (show) { for (const r of target) d.hidden.delete(r); }
        else { for (const r of target) d.hidden.add(r); for (const r of target) { d.sticks.delete(r); d.spheres.delete(r); } }
      } else if (['cartoon', 'tube', 'backbone'].includes(what)) {
        if (!s && show) { app.setRepresentation(REPS.cartoon); d.hidden.clear(); }
        else for (const r of target) (show ? d.hidden.delete(r) : d.hidden.add(r));
      } else if (['sticks', 'atoms', 'lines'].includes(what)) {
        if (!s) app.setRepresentation(show ? 'sticks' : 'cartoon');
        else for (const r of target) (show ? d.sticks.add(r) : d.sticks.delete(r));
      } else if (['spheres', 'spacefill', 'vdw'].includes(what)) {
        if (!s) app.setRepresentation(show ? 'spacefill' : 'cartoon');
        else for (const r of target) (show ? d.spheres.add(r) : d.spheres.delete(r));
      } else if (what === 'ballstick' || what === 'ball+stick') {
        if (show) app.setRepresentation('ballstick'); else app.setRepresentation('cartoon');
      } else if (what === 'labels' || what === 'label') {
        if (show) for (const r of target) { if (!d.labels.has(r)) d.labels.set(r, null); }
        else { if (s) for (const r of target) d.labels.delete(r); else { d.labels.clear(); app.setLabelMode('none'); } }
      } else if (what === 'ensemble' || what === 'ghosts' || what === 'overlay') {
        if (show !== app.state.ghosts) app.toggleGhosts();
        return say(`ensemble overlay ${show ? 'on' : 'off'}`);
      } else if (what === 'putty' || what === 'flexibility') {
        if (show !== app.state.putty) app.togglePutty();
        return say(`flexibility putty ${show ? 'on' : 'off'}`);
      } else if (what === 'surface' || what === 'surf') {
        app.surface.set(show, { residues: show ? (s ? [...s] : null) : undefined });
        return say(`surface ${show ? 'shown' + (s ? ' around ' + describeResidues(ens, s) : '') : 'hidden'}`);
      } else if (what === '2d' || what === 'ss2d') {
        app.ss2d.setVisible(show);
        return say(`2D window ${show ? 'shown' : 'hidden'}`);
      } else if (what === 'sequence' || what === 'seq') {
        app.seqbar?.setVisible(show);
        return say(`sequence bar ${show ? 'shown' : 'hidden'}`);
      } else throw new CommandError(`unknown representation "${what}" (cartoon, sticks, spheres, ballstick, surface, labels, everything, ensemble, putty, 2d, sequence)`);
      app.applyDisplay();
      return say(`${verb} ${what || 'everything'}${s ? ' for ' + describeResidues(ens, s) : ''}`);
    };
    this.def(['show'], { args: 'what[, selection]', help: 'Show cartoon / sticks / spheres / ballstick / labels / everything (for a selection, or globally); also ensemble, putty, 2d, sequence.', run: showHide(true) });
    this.def(['hide'], { args: 'what[, selection]', help: 'Hide a representation (for a selection, or globally).', run: showHide(false) });
    this.def(['rep', 'representation'], { args: 'cartoon | cartoon+atoms | sticks | ballstick | spheres', help: 'Set the representation of the whole molecule.', run: args => {
      const r = REPS[(args[0] || '').toLowerCase()];
      if (!r) throw new CommandError('representations: cartoon, cartoon+atoms, sticks, ballstick, spheres');
      app.setRepresentation(r);
      return say(`representation ${args[0]}`);
    } });
    this.def(['style'], { args: 'studio | illustrative | glossy | depth | plain', help: 'Rendering style.', run: args => {
      const s = (args[0] || '').toLowerCase();
      if (!app.viewer.styles[s]) throw new CommandError(`styles: ${Object.keys(app.viewer.styles).join(', ')}`);
      app.setStyle(s);
      return say(`style ${s}`);
    } });

    // -------- colour
    this.def(['color', 'colour'], { args: 'color | scheme[, selection]', help: 'Colour residues (name or #hex), reset them (default), or switch the colour scheme.', run: args => {
      need(args, 1, 'color red, resi 5-10   ·   color helix   ·   color default, all');
      const ens = this.ens();
      const what = args[0].trim(), low = what.toLowerCase();
      const scheme = SCHEMES.find(s => s.id === low || s.label.toLowerCase() === low);
      const d = app.displayOf();
      if (scheme && args.length < 2) { app.setColorScheme(scheme.id); return say(`colour scheme: ${scheme.label}`); }
      const s = this.sel(args.slice(1).join(','), 'all');
      if (low === 'default' || low === 'scheme' || low === 'reset') { for (const r of s) d.colors.delete(r); app.applyDisplay(); return say(`custom colours removed from ${describeResidues(ens, s)}`); }
      const hex = parseColor(what);
      if (!hex) throw new CommandError(`unknown colour "${what}" (names: ${Object.keys(COLOR_NAMES).join(', ')}, #rrggbb, or a scheme: ${SCHEMES.map(x => x.id).join(', ')})`);
      for (const r of s) d.colors.set(r, hex);
      app.applyDisplay();
      return say(`${hex} → ${describeResidues(ens, s)}`);
    } });
    this.def(['bg', 'background'], { args: 'color', help: 'Viewport background colour.', run: args => {
      const hex = parseColor(args[0]);
      if (!hex) throw new CommandError('usage: bg white | black | #f7f6f2 …');
      app.setBackground(hex);
      return say(`background ${hex}`);
    } });

    // -------- labels
    this.def(['label'], { args: 'selection[, text] | off', help: 'Label residues (optionally with your own text); "label off" removes custom labels.', run: args => {
      need(args, 1, 'label resi 20-22   ·   label resi 34, anticodon   ·   label off');
      const d = app.displayOf();
      if (['off', 'none', 'clear'].includes(args[0].toLowerCase())) { d.labels.clear(); app.applyDisplay(); return say('custom labels removed'); }
      const s = this.sel(args[0]);
      const text = args.length > 1 ? args.slice(1).join(', ') : null;
      for (const r of s) d.labels.set(r, text);
      app.applyDisplay();
      return say(`labelled ${describeResidues(app.ens, s)}`);
    } });

    // -------- view
    const atomsOf = s => { const ens = this.ens(), out = []; for (const r of s) { const res = ens.residues[r]; for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) out.push(a); } return out; };
    const centroid = atoms => { const c = app.display; const m = [0, 0, 0]; for (const a of atoms) { m[0] += c[3 * a]; m[1] += c[3 * a + 1]; m[2] += c[3 * a + 2]; } return m.map(v => v / Math.max(1, atoms.length)); };
    this.def(['zoom'], { args: '[selection]', help: 'Centre and zoom on a selection (whole molecule without one).', run: args => {
      const ens = this.ens();
      if (!args[0]) { app.viewer.resetView(); return say('zoomed to the whole molecule'); }
      const s = this.sel(args.join(','));
      if (!s.size) throw new CommandError('the selection is empty');
      const atoms = atomsOf(s), c = centroid(atoms), disp = app.display;
      let r = 0;
      for (const a of atoms) r = Math.max(r, Math.hypot(disp[3 * a] - c[0], disp[3 * a + 1] - c[1], disp[3 * a + 2] - c[2]));
      app.viewer.centerOn(c);
      app.viewer.controls.zoom = Math.max(0.3, Math.min(12, ens.fitRadius / Math.max(4, r + 3)));
      app.viewer.needsRender = true;
      return say(`zoom on ${describeResidues(ens, s)}`);
    } });
    this.def(['center', 'centre'], { args: '[selection]', help: 'Centre the view on a selection.', run: args => {
      const s = this.sel(args.join(','), 'all');
      app.viewer.centerOn(centroid(atomsOf(s)));
      return say(`centred on ${describeResidues(app.ens, s)}`);
    } });
    this.def(['orient'], { args: '[selection]', help: 'Turn the view so the selection’s long axis is horizontal, then centre on it.', run: args => {
      const s = this.sel(args.join(','), 'all');
      const ens = this.ens(), disp = app.display;
      const pts = [...s].map(r => ens.trace[r]).filter(a => a >= 0).map(a => [disp[3 * a], disp[3 * a + 1], disp[3 * a + 2]]);
      if (pts.length < 3) throw new CommandError('select at least 3 residues to orient on');
      const q = principalQuaternion(pts);
      app.viewer.controls.quaternion.copy(q);
      app.viewer.centerOn(centroid(atomsOf(s)), false);
      app.viewer.needsRender = true;
      return say(`oriented on ${describeResidues(ens, s)}`);
    } });
    this.def(['reset'], { help: 'Reset the view.', run: () => { app.viewer.resetView(); return say('view reset'); } });
    this.def(['turn', 'rotate'], { args: 'x | y | z, degrees', help: 'Rotate the view about a screen axis.', run: args => {
      need(args, 2, 'turn y, 90');
      const ax = { x: [1, 0, 0], y: [0, 1, 0], z: [0, 0, 1] }[args[0].toLowerCase()];
      const deg = parseFloat(args[1]);
      if (!ax || !Number.isFinite(deg)) throw new CommandError('usage: turn x|y|z, degrees');
      const c = app.viewer.controls;
      const q = c.quaternion.clone().setFromAxisAngle({ x: ax[0], y: ax[1], z: ax[2] }, deg * Math.PI / 180);
      c.quaternion.premultiply(q);
      app.viewer.needsRender = true;
      return say(`turned ${deg}° about ${args[0]}`);
    } });
    this.def(['spin'], { args: 'on | off', help: 'Auto-rotate.', run: args => {
      const on = !['off', '0', 'no'].includes((args[0] || 'on').toLowerCase());
      app.viewer.controls.spin = on ? 0.5 : 0; app.viewer.needsRender = true; app._syncToolbar();
      return say(`spin ${on ? 'on' : 'off'}`);
    } });
    this.def(['projection'], { args: 'perspective | orthographic', help: 'Camera projection.', run: args => {
      const p = (args[0] || '').toLowerCase().startsWith('o') ? 'orthographic' : 'perspective';
      app.setProjection(p); return say(`projection ${p}`);
    } });

    // -------- conformations
    this.def(['conf', 'goto'], { args: 'id | label | number | first | last | previous', help: 'Morph to a conformation (number = position in the current order).', run: args => {
      const k = this.conf(args[0]); app.goTo(k); return say(`→ ${app.ens.shortId(k)}`);
    } });
    this.def(['next'], { help: 'Next conformation.', run: () => { app.next(); return say(`→ ${app.ens?.shortId(app.state.conf)}`); } });
    this.def(['prev', 'previous'], { help: 'Previous conformation.', run: () => { app.prev(); return say(`→ ${app.ens?.shortId(app.state.conf)}`); } });
    this.def(['play'], { help: 'Play the conformation tour.', run: () => { app.play(); return say('playing'); } });
    this.def(['stop', 'pause'], { help: 'Stop playback.', run: () => { app.stop(); return say('stopped'); } });
    this.def(['morph'], { args: 'from, to[, seconds]', help: 'Show the morph from one conformation to another.', run: args => {
      need(args, 2, 'morph conf_001, conf_005, 2');
      const a = this.conf(args[0]), b = this.conf(args[1]);
      const sec = args[2] ? parseFloat(args[2]) : app.morphDuration;
      app.goTo(a, { animate: false, force: true });
      setTimeout(() => app.goTo(b, { duration: sec }), 50);
      return say(`${app.ens.shortId(a)} → ${app.ens.shortId(b)} in ${sec} s`);
    } });
    this.def(['order'], { args: 'rank | energy | path | input', help: 'Order of conformations for browsing and playback.', run: args => {
      const m = { rank: 'rank', energy: 'metric', metric: 'metric', path: 'path', similarity: 'path', input: 'input' }[(args[0] || '').toLowerCase()];
      if (!m) throw new CommandError('orders: rank, energy, path, input');
      app.setOrder(m); return say(`order: ${args[0]}`);
    } });
    this.def(['metric'], { args: 'key or label', help: 'Energy/metric shown in the charts.', run: (args, raw) => {
      const ens = this.ens(), q = raw.trim().toLowerCase();
      const m = ens.metrics.find(x => x.key.toLowerCase() === q) || ens.metrics.find(x => x.label.toLowerCase() === q) || ens.metrics.find(x => x.label.toLowerCase().includes(q));
      if (!m) throw new CommandError(`metrics: ${ens.metrics.map(x => x.key).join(', ')}`);
      app.setMetric(m.key); return say(`metric: ${m.label}`);
    } });

    // -------- measurements
    const measure = kind => (args) => {
      const n = { distance: 2, angle: 3, dihedral: 4 }[kind];
      let name = null;
      let specs = args.filter(a => a !== '');
      if (specs.length === n + 1) { name = specs[0]; specs = specs.slice(1); }
      if (specs.length !== n) throw new CommandError(`usage: ${kind} [name,] ${Array.from({ length: n }, (_, i) => `atom${i + 1}`).join(', ')}   (atom = residue/atom, e.g. 12/N1, 34/C4', pos 5/P)`);
      const atoms = specs.map(s => this.atom(s));
      const m = app.measure.add(kind, atoms, { name });
      const v = app.measure.value(m, app.display);
      const vals = app.ens.conformations.map(c => c.values[m.key]);
      return [say(`${app.measure.describe(m)} = ${fmt(v, 2)} ${kind === 'distance' ? 'Å' : '°'}`), say(`  across the ensemble: ${fmt(Math.min(...vals), 2)} – ${fmt(Math.max(...vals), 2)} (metric "${m.key}", plot it with: metric ${m.key})`)];
    };
    this.def(['distance', 'dist'], { args: '[name,] atom1, atom2', help: 'Measure a distance (follows morphs; plotted per conformation).', run: measure('distance') });
    this.def(['angle'], { args: '[name,] atom1, atom2, atom3', help: 'Measure an angle.', run: measure('angle') });
    this.def(['dihedral', 'torsion'], { args: '[name,] atom1, atom2, atom3, atom4', help: 'Measure a dihedral angle.', run: measure('dihedral') });
    this.def(['measurements'], { args: '[clear]', help: 'List (or clear) measurements.', run: args => {
      if ((args[0] || '').toLowerCase() === 'clear') { app.measure.clear(); return say('measurements removed'); }
      const list = app.measure.list();
      return list.length ? list.map(m => say(`${app.measure.describe(m)} = ${fmt(app.measure.value(m, app.display), 2)}`)) : say('no measurements');
    } });

    // -------- comparison
    const rmsdOf = (a, b, s, fitIt) => {
      const ens = this.ens();
      const atoms = [...s].map(r => ens.trace[r]).filter(x => x >= 0);
      if (atoms.length < 3) throw new CommandError('need at least 3 residues with C4′ atoms');
      const A = ens.coords[a], B = Float32Array.from(ens.coords[b]);
      if (fitIt) { const f = kabsch(B, A, atoms); applyTransform(B, f.R, f.t); }
      let sum = 0;
      const per = [];
      for (const r of s) {
        const i = ens.trace[r];
        if (i < 0) continue;
        const d = Math.hypot(A[3 * i] - B[3 * i], A[3 * i + 1] - B[3 * i + 1], A[3 * i + 2] - B[3 * i + 2]);
        per.push([r, d]); sum += d * d;
      }
      return { rmsd: Math.sqrt(sum / per.length), per: per.sort((x, y) => y[1] - x[1]), n: per.length };
    };
    this.def(['align', 'super'], { args: 'confA, confB[, selection]', help: 'C4′ RMSD between two conformations after optimal superposition on the selection; lists the most displaced residues.', run: args => {
      need(args, 2, 'align conf_001, conf_004[, selection]');
      const a = this.conf(args[0]), b = this.conf(args[1]);
      const s = this.sel(args.slice(2).join(','), 'all');
      const r = rmsdOf(a, b, s, true), ens = app.ens;
      return [say(`RMSD ${fmt(r.rmsd, 2)} Å over ${r.n} C4′ atoms (${ens.shortId(a)} vs ${ens.shortId(b)}, fitted on ${describeResidues(ens, s)})`),
        say(`  most displaced: ${r.per.slice(0, 6).map(([x, d]) => `${ens.resLabel(x)} ${fmt(d, 1)} Å`).join(', ')}`)];
    } });
    this.def(['rms', 'rms_cur'], { args: 'confA, confB[, selection]', help: 'C4′ RMSD in the current superposition (no fitting).', run: args => {
      need(args, 2, 'rms conf_001, conf_004[, selection]');
      const a = this.conf(args[0]), b = this.conf(args[1]);
      const s = this.sel(args.slice(2).join(','), 'all');
      const r = rmsdOf(a, b, s, false);
      return say(`RMSD ${fmt(r.rmsd, 2)} Å over ${r.n} C4′ atoms without refitting`);
    } });
    this.def(['fit'], { args: 'selection | reset', help: 'Superpose every conformation on the selection (C4′); "fit reset" restores the original superposition.', run: (args, raw) => {
      if ((args[0] || '').toLowerCase() === 'reset') { app.resetSuperposition(); return say('original superposition restored'); }
      const s = this.sel(raw);
      if (s.size < 3) throw new CommandError('select at least 3 residues to fit on');
      app.select([...s], 'set');
      app.superposeOnSelection();
      return say(`all conformations fitted on ${describeResidues(app.ens, s)}`);
    } });
    this.def(['pairs'], { args: '[conformation]', help: 'List the base pairs of a conformation (default: the shown one).', run: args => {
      const ens = this.ens(), k = args[0] ? this.conf(args[0]) : app.state.conf;
      const ps = ens.conformations[k].pairs;
      return [say(`${ens.shortId(k)}: ${ps.length} base pairs`), say('  ' + ps.map(([i, j, t, pk]) => `${ens.resLabel(i)}-${ens.resLabel(j)}${pk ? '(pk)' : ''}`).join(' ')), say(`  ${ens.sequence}`), say(`  ${ens.conformations[k].dbn}`)];
    } });
    this.def(['info'], { help: 'Summary of the open ensemble and the shown conformation.', run: () => {
      const ens = this.ens(), k = app.state.conf, c = ens.conformations[k], m = ens.metricByKey[app.state.metric];
      const v = ens.value(k, app.state.metric);
      return [say(`${ens.name}: ${ens.nRes} nt, ${ens.nAtoms} atoms, ${ens.n} conformations (${ens.doc.origin})`), say(`sequence  ${ens.sequence}`),
        say(`shown     ${c.dbn}  ${ens.shortId(k)}${v !== null ? ` · ${m.label} ${fmt(v)} ${m.unit || ''}` : ''}`),
        say(`named selections: ${[...app.displayOf().named.keys()].join(', ') || 'none'} · measurements: ${app.measure.list().length} · scenes: ${app.scenes.list().length}`)];
    } });

    // -------- scenes
    this.def(['view', 'scene'], { args: 'save name | name | list | delete name', help: 'Save and recall views (camera, conformation, colours, representation, selection).', run: args => {
      this.ens();
      const a = (args[0] || 'list').trim();
      const [verb, ...rest] = a.split(/\s+/);
      const name = (rest.join(' ') || args[1] || '').trim();
      if (verb === 'save' || verb === 'store') { const s = app.scenes.save(name || undefined); return say(`scene "${s.name}" saved`); }
      if (verb === 'list' || verb === 'ls') { const l = app.scenes.list(); return say(l.length ? l.map(s => s.name).join(', ') : 'no scenes yet (view save <name>)'); }
      if (verb === 'delete' || verb === 'del') { if (!app.scenes.remove(name)) throw new CommandError(`no scene "${name}"`); return say(`scene "${name}" deleted`); }
      const target = verb === 'recall' ? name : a;
      try { app.scenes.recall(target); } catch (e) { throw new CommandError(e.message); }
      return say(`scene "${target}"`);
    } });

    // -------- objects
    this.def(['ls', 'list'], { help: 'Named selections, measurements and scenes.', run: () => {
      const d = app.displayOf(), ens = this.ens();
      const out = [...d.named.entries()].map(([n, s]) => say(`selection ${n}: ${describeResidues(ens, s)}`));
      out.push(...app.measure.list().map(m => say(`measurement ${app.measure.describe(m)}`)));
      out.push(...app.scenes.list().map(s => say(`scene ${s.name} (${ens.conformations.find(c => c.id === s.confId)?.id || ''})`)));
      return out.length ? out : say('nothing stored yet');
    } });
    this.def(['delete', 'remove'], { args: 'name', help: 'Delete a named selection, measurement or scene.', run: args => {
      need(args, 1, 'delete name');
      const n = args[0], d = app.displayOf();
      if (d.named.delete(n)) { app.emit('named'); return say(`selection ${n} deleted`); }
      if (app.measure.remove(n)) return say(`measurement ${n} deleted`);
      if (app.scenes.remove(n)) return say(`scene ${n} deleted`);
      throw new CommandError(`nothing called "${n}"`);
    } });

    // -------- settings
    this.def(['set'], { args: 'setting, value', help: 'morph (s), hold (s), loop on|off, labels …, quality high|low, outline #hex, elements on|off, transparency 0–1[, object].', run: args => {
      need(args, 2, 'set morph, 1.5');
      const k = args[0].toLowerCase(), v = args[1].trim(), on = !['off', 'no', '0', 'false'].includes(v.toLowerCase());
      switch (k) {
        case 'morph': app.morphDuration = Math.max(0.1, parseFloat(v) || 1.2); break;
        case 'hold': app.holdTime = Math.max(0, parseFloat(v) || 0); break;
        case 'loop': app.loop = on; break;
        case 'labels': app.setLabelMode(v); break;
        case 'quality': app.viewer.settings.quality = v === 'low' ? 'low' : 'high'; break;
        case 'outline': { const h = parseColor(v); if (!h) throw new CommandError('outline needs a colour'); app.viewer.settings.outlineColor = h; break; }
        case 'elements': app.viewer.settings.elementColors = on; app.refreshColors(); break;
        case 'transparency': case 'opacity': {
          let x = parseFloat(v); if (!Number.isFinite(x)) throw new CommandError('set transparency, 0.5[, cartoon|bases|atoms|surface|overlay|all]');
          if (x > 1) x /= 100;
          const obj = (args[2] || 'all').trim().toLowerCase();
          if (!['cartoon', 'bases', 'atoms', 'surface', 'overlay', 'all'].includes(obj)) throw new CommandError('objects: cartoon, bases, atoms, surface, overlay, all');
          app.setTransparency(obj, k === 'transparency' ? 1 - x : x);
          return say(`${obj} ${k} ${Math.round(x * 100)} %`);
        }
        default: throw new CommandError('settings: morph, hold, loop, labels, quality, outline, elements, transparency, opacity');
      }
      app.savePrefs(); app.viewer.needsRender = true;
      return say(`${k} = ${v}`);
    } });

    // -------- files & export
    this.def(['png', 'snapshot'], { args: '[width[, height[, dpi]]]', help: 'Save a PNG image of the view (default: 2× the window, 300 dpi).', run: async args => {
      const ens = this.ens();
      const [cw, ch] = app.viewer.cssSize;
      const w = parseInt(args[0], 10) || cw * 2, h = parseInt(args[1], 10) || Math.round(w * ch / cw), dpi = parseInt(args[2], 10) || 300;
      const { exportFigure } = await import('../export/figure.js');
      const r = await exportFigure(app, { format: 'png', width: w, height: h, dpi, transparent: false, background: null, legend: true, labels: true, panel: false });
      const saved = await saveBlob(r.blob, `${slug(ens.name)}_${slug(ens.shortId(app.state.conf))}.png`);
      return say(saved ? `saved ${saved} (${w}×${h} px, ${dpi} dpi)` : 'not saved');
    } });
    this.def(['save'], { args: 'rcv | pdb | pdb all | csv | workspace', help: 'Save the ensemble (.rcv), the shown or all conformations (PDB), the energy table or the workspace.', run: async args => {
      this.ens();
      const what = (args[0] || 'rcv').toLowerCase();
      if (what.endsWith('rcv')) await app.saveRcv();
      else if (what.startsWith('pdb') && /all/.test(what)) await app.exportPDB(true);
      else if (what.startsWith('pdb')) await app.exportPDB(false);
      else if (what.endsWith('csv')) await app.exportCSV();
      else if (what.startsWith('work') || what.endsWith('json')) await app.saveWorkspace();
      else throw new CommandError('save rcv | pdb | pdb all | csv | workspace');
      return say(`save ${what}`);
    } });
    this.def(['video'], { args: 'transition A, B | tour | rotate', help: 'Open the video dialog (transition from A to B, tour or rotation).', run: args => {
      const ens = this.ens();
      const [kind, ...rest] = (args[0] || 'tour').split(/\s+/);
      if (kind.toLowerCase() === 'transition') {
        const a = this.conf(rest.join(' ') || 'previous'), b = this.conf(args[1] || 'current');
        app.exportTransition(a, b);
        return say(`transition video ${ens.shortId(a)} → ${ens.shortId(b)}`);
      }
      app.openVideoDialog({ content: kind.toLowerCase() });
      return say(`video: ${kind}`);
    } });
    // -------- analyses (RCV 0.7): pathways, conditions, probing data
    this.def(['pathway', 'transition'], { args: 'A, B', help: 'Refolding paths, barrier, bottleneck and rigid modules from conformation A to B (Pathway tab).', run: async args => {
      need(args, 2, 'pathway conf_001, conf_005');
      const ens = this.ens();
      const a = this.conf(args[0]), b = this.conf(args[1]);
      if (a === b) throw new CommandError('A and B must be different conformations');
      const P = app.pathway;
      P.a = a; P.b = b; P.result = P.scan = null;
      app.showTab('pathway');
      await P.load();
      if (P.error) throw new CommandError(P.error);
      const r = P.result;
      if (r.A.dbn === r.B.dbn) return say(`${ens.shortId(a)} → ${ens.shortId(b)}: same 2D structure (tertiary rearrangement), C4′ RMSD ${fmt(r.rmsd_ab, 1)} Å`);
      if (r.same_basin) return say(`${ens.shortId(a)} → ${ens.shortId(b)}: same 2D basin (A ${fmt(r.A.above_basin, 1)}, B ${fmt(r.B.above_basin, 1)} kcal/mol above its minimum), no barrier between them; C4′ RMSD ${fmt(r.rmsd_ab, 1)} Å`);
      return say(`${ens.shortId(a)} → ${ens.shortId(b)}: effective barrier ${fmt(r.effective_barrier, 1)} kcal/mol (direct path ${fmt(r.direct.barrier, 1)}), ΔG ${fmt(r.dG, 1)} kcal/mol, half-life ≈ ${fmtTime(r.half_life_s)}; cooperative unit ${r.bottleneck?.unit.length ?? 0} bp; ${r.motion?.modules.length ?? 0} rigid modules`);
    } });
    this.def(['phiscan'], { help: 'φ-value scan of every point mutation along the current pathway (run "pathway A, B" first).', run: async () => {
      const P = app.pathway;
      if (!P.result) throw new CommandError('run "pathway A, B" first');
      if (P.result.same_basin || P.result.A.dbn === P.result.B.dbn) throw new CommandError('A and B are in the same 2D basin: there is no barrier to scan');
      await P.runScan();
      if (!P.scan) throw new CommandError('the scan failed');
      const ens = this.ens(), name = x => `${x.from}${ens.residues[x.pos].num}${x.to}`;
      return [say(`barrier ${fmt(P.scan.barrier_wt, 1)} kcal/mol; ${3 * ens.nRes} mutants`),
        say(`fastest: ${P.scan.faster.slice(0, 5).map(x => `${name(x)} (${fmt(x.ddg_ts, 1)})`).join(', ') || '–'}`),
        say(`slowest: ${P.scan.slower.slice(0, 5).map(x => `${name(x)} (+${fmt(x.ddg_ts, 1)})`).join(', ') || '–'}`)];
    } });
    this.def(['mechanism'], { args: '[staged | rigid]', help: 'Play the mechanism morph of the current pathway (rigid modules in the order of the base-pair events).', run: args => {
      const P = app.pathway;
      if (!P.result?.motion) throw new CommandError('run "pathway A, B" first');
      const staged = (args[0] || 'staged').toLowerCase() !== 'rigid';
      P.play(staged);
      return say(staged ? 'playing the mechanism' : 'playing the rigid-module morph');
    } });
    this.def(['conditions', 'cond'], { args: '[°C][, Na mM][, ligand µM | no ligand]', help: 'Show or set the conditions of the 2D analyses, e.g. "conditions 25, 150 mM, theophylline 10 µM" or "conditions standard".', run: args => {
      const ens = this.ens();
      if (!args.length) return say(conditionsText(conditionsOf(app, ens), ens));
      const patch = {};
      for (const tok of args.map(t => t.trim()).filter(Boolean)) {
        const t = tok.toLowerCase();
        if (t === 'standard') { Object.assign(patch, { temperature: 37, sodium: STANDARD_SALT, ligand: null }); continue; }
        if (t === 'no ligand' || t === 'none' || t === 'ligand none') { patch.ligand = null; continue; }
        let m;
        if ((m = t.match(/^(-?[\d.]+)\s*(°c|c)?$/))) { patch.temperature = Math.min(95, Math.max(0, parseFloat(m[1]))); continue; }
        if ((m = t.match(/^(?:na\+?\s*)?([\d.]+)\s*(mm|m)$/))) { const v = parseFloat(m[1]) / (m[2] === 'mm' ? 1000 : 1); patch.sodium = Math.min(3, Math.max(0.005, v)); continue; }
        if ((m = t.match(/^([a-z][a-z-]*)\s+([\d.]+)\s*(µm|um|nm|mm)?$/)) && LIGANDS[m[1]]) {
          const f = { 'nm': 1e-3, 'mm': 1e3 }[m[3]] || 1;
          patch.ligand = { ...LIGANDS[m[1]], preset: m[1], conc_uM: parseFloat(m[2]) * f };
          continue;
        }
        throw new CommandError(`cannot read "${tok}": use a temperature (25), a salt concentration (150 mM) or a ligand (theophylline 10 µM)`);
      }
      setConditions(app, patch, ens);
      return say(conditionsText(conditionsOf(app, ens), ens));
    } });
    this.def(['probing', 'shape'], { args: 'path[, SHAPE | DMS] | clear', help: 'Load a SHAPE or DMS reactivity file (normalised by the 2–8 % rule) and fit the ensemble to it; "probing clear" removes it.', run: async (args, raw) => {
      need(args, 1, 'probing C:\\data\\tar.shape, SHAPE');
      const ens = this.ens();
      if (args[0].trim().toLowerCase() === 'clear') { app.probing.clear(ens); return say('probing data removed'); }
      const kind = (args[1] || 'SHAPE').trim().toUpperCase();
      if (!['SHAPE', 'DMS'].includes(kind)) throw new CommandError('the reagent is SHAPE or DMS');
      const path = args[0].trim().replace(/^["']|["']$/g, '');
      const r = await app.api.post('/api/read-text', { path });
      const parsed = await app.api.post('/api/probing/parse', { text: r.text, sequence: ens.sequence, kind, normalise: true });
      await app.probing.set(ens, { name: path.split(/[\\/]/).pop(), kind, values: parsed.normalised, raw: parsed.values, factor: parsed.factor, format: parsed.format });
      const f = ens.probing?.fit;
      return say(`${kind}: ${parsed.n_values} of ${ens.nRes} nucleotides${f ? `; agreement (AUC) of ${ens.shortId(app.state.conf)} ${fmt(f.agreement[app.state.conf], 2)}${f.reweighting ? `; reweighting χ² ${fmt(f.reweighting.chi2_prior, 2)} → ${fmt(f.reweighting.chi2, 2)}` : ''}` : ''}`);
    } });
    this.def(['load', 'open'], { args: 'path', help: 'Open an .rcv file, structure files or a prediction folder on this computer.', run: async (args, raw) => {
      need(args, 1, 'load C:\\data\\tar.rcv');
      await app.openServerPaths([raw.trim().replace(/^["']|["']$/g, '')]);
      return say(`opened ${app.ens?.name || ''}`);
    } });
    this.def(['fetch'], { args: 'PDB id[, chain]', help: 'Download an RNA structure from the RCSB PDB (needs internet) and open it, e.g. fetch 1ANR.', run: async args => {
      need(args, 1, 'fetch 1ANR');
      const id = args[0].trim();
      if (!/^[0-9][A-Za-z0-9]{3}$/.test(id)) throw new CommandError('a PDB id has 4 characters, e.g. 1ANR');
      const s = await app.api.post('/api/fetch', { id, chain: args[1] || null });
      await app.loadServerDocument(s.id);
      return say(`${id.toUpperCase()}: ${s.n_conformations} model(s), ${s.length} nt`);
    } });
    this.def(['run', '@'], { args: 'path', help: 'Run a command script (a text file with one command per line).', run: async (args, raw) => {
      need(args, 1, 'run C:\\data\\figure.rcvs');
      const r = await app.api.post('/api/read-text', { path: raw.trim().replace(/^["']|["']$/g, '') });
      const lines = r.text.split(/\r?\n/);
      const out = [say(`running ${r.path} (${lines.filter(l => l.trim() && !l.trim().startsWith('#')).length} commands)`)];
      for (const l of lines) {
        if (!l.trim() || l.trim().startsWith('#')) continue;
        out.push({ text: `RCV> ${l.trim()}`, kind: 'cmd' });
        try { out.push(...await this.execute(l)); } catch (e) { out.push({ text: e.message, kind: 'err' }); break; }
      }
      return out;
    } });
    this.def(['echo', 'print'], { args: 'text', help: 'Print text (useful in scripts).', run: (a, raw) => say(raw) });
  }
}

// Rotation that lays the principal axis of the points along screen X (and the second along Y).
function principalQuaternion(pts) {
  const m = [0, 0, 0];
  for (const p of pts) { m[0] += p[0]; m[1] += p[1]; m[2] += p[2]; }
  m[0] /= pts.length; m[1] /= pts.length; m[2] /= pts.length;
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] += (p[i] - m[i]) * (p[j] - m[j]);
  const lead = (M, orth) => {
    let v = [1, 0.5, 0.25];
    for (let it = 0; it < 80; it++) {
      let w = [0, 1, 2].map(i => M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2]);
      if (orth) { const d = w[0] * orth[0] + w[1] * orth[1] + w[2] * orth[2]; w = w.map((x, i) => x - d * orth[i]); }
      const l = Math.hypot(...w) || 1;
      v = w.map(x => x / l);
    }
    return v;
  };
  const e1 = lead(C, null), e2 = lead(C, e1);
  const e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
  // rows e1,e2,e3 form the rotation model -> view; convert to a quaternion
  const R = [e1, e2, e3];
  const tr = R[0][0] + R[1][1] + R[2][2];
  let w, x, y, z;
  if (tr > 0) { const s = Math.sqrt(tr + 1) * 2; w = 0.25 * s; x = (R[2][1] - R[1][2]) / s; y = (R[0][2] - R[2][0]) / s; z = (R[1][0] - R[0][1]) / s; }
  else if (R[0][0] > R[1][1] && R[0][0] > R[2][2]) { const s = Math.sqrt(1 + R[0][0] - R[1][1] - R[2][2]) * 2; w = (R[2][1] - R[1][2]) / s; x = 0.25 * s; y = (R[0][1] + R[1][0]) / s; z = (R[0][2] + R[2][0]) / s; }
  else if (R[1][1] > R[2][2]) { const s = Math.sqrt(1 + R[1][1] - R[0][0] - R[2][2]) * 2; w = (R[0][2] - R[2][0]) / s; x = (R[0][1] + R[1][0]) / s; y = 0.25 * s; z = (R[1][2] + R[2][1]) / s; }
  else { const s = Math.sqrt(1 + R[2][2] - R[0][0] - R[1][1]) * 2; w = (R[1][0] - R[0][1]) / s; x = (R[0][2] + R[2][0]) / s; y = (R[1][2] + R[2][1]) / s; z = 0.25 * s; }
  return { x, y, z, w, isQuaternion: true, _q: true };
}
