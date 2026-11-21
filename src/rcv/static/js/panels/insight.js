// Advanced analyses (0.7.4), shown as closed sections inside the existing tabs so that the interface stays as it
// was: nothing is computed until a section is opened and its button pressed.
//   Switch tab   — experiments that tell state A from B (rescue mutants, probing positions, label pairs); variants
//   Pathway tab  — folding during transcription
//   Details tab  — structural motifs, target sites for antisense oligonucleotides, comparison of two ensembles
// The computations are in rcv/insight.py (POST /api/insight/<name>).
import { el, svg, fmt, fmtSigned } from '../util.js';
import { icon } from '../icons.js';
import { toast, select } from '../ui.js';
import { registerScheme, sampleRamp, ORANGE_RAMP, BLUE_RAMP } from '../colors.js';
import { requestConditions } from '../conditions.js';

const store = new Map();          // `${section}|${ensemble id}` -> { open, result, loading, error, key, ... }
const stateOf = (id, ens) => { const k = `${id}|${ens.id}`; if (!store.has(k)) store.set(k, { open: false }); return store.get(k); };
const pct = x => (x === null || x === undefined ? '–' : `${fmt(100 * x, x > 0.995 || x < 0.005 ? 1 : 0)}%`);
const slug = s => String(s || 'rna').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 40) || 'rna';
const range = (a, b) => Array.from({ length: b - a + 1 }, (_, i) => a + i);
const span = (ens, a, b) => (a === b ? `${ens.residues[a].num}` : `${ens.residues[a].num}–${ens.residues[b].num}`);

// A closed section; `body(st, draw)` returns its content and is called only while the section is open.
export function fold(ens, id, title, hint, body) {
  const st = stateOf(id, ens);
  const det = el('details', { class: 'section adv', open: st.open ? true : null },
    el('summary', { class: 'section-title adv-title' }, el('span', { text: title }), el('span', { class: 'adv-tag', text: 'advanced' })));
  const box = el('div', { class: 'adv-body' });
  det.append(box);
  const draw = () => {
    if (!det.open) return;
    box.replaceChildren(el('div', { class: 'desc', text: hint }), ...[body(st, draw)].flat(3).filter(Boolean));
  };
  st.draw = draw;                                   // the newest drawing of this section receives the results
  st.det = det;
  det.addEventListener('toggle', () => {
    st.open = det.open;
    // the tab was drawn again between the click and this event: pass the click on to the section now shown
    if (st.det !== det && !det.isConnected && st.det.open !== det.open) st.det.open = det.open;
    draw();
  });
  draw();
  return det;
}

async function run(app, ens, st, name, extra = {}) {
  if (!app.api.online) { st.error = 'This analysis needs the local RCV engine.'; st.draw(); return; }
  st.loading = true; st.error = null; st.draw();
  try {
    const body = { conditions: requestConditions(app, ens, { probing: false }), ...extra };
    if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
    st.result = await app.api.post(`/api/insight/${name}`, body);
  } catch (err) { st.error = err.message; st.result = null; }
  st.loading = false;
  st.draw();
}

export const status = st => (st.loading ? el('div', { class: 'row' }, el('span', { class: 'spinner small' }), ' computing…')
  : st.error ? el('div', { class: 'warnbox', text: st.error }) : null);

function pick(app, residues) { app.select(residues, 'set'); app.centerOnSelection(); }

export function table(head, rows) {
  return el('div', { class: 'wl-table adv-table' }, el('table', { class: 'table' },
    el('thead', {}, el('tr', {}, ...head.map(h => el('th', { class: h.num ? 'num' : null, title: h.title || null, text: h.label })))),
    el('tbody', {}, ...rows)));
}

function copyButton(label, text) {
  return el('button', { class: 'btn small ghost', onclick: () => navigator.clipboard?.writeText(text()).then(() => toast('Copied as a tab-separated table.', 'ok')) }, icon('copy'), label);
}

function colourBy(app, id, label, values, { ramp = ORANGE_RAMP, title, min = '0', max, scale = null } = {}) {
  const top = scale ?? Math.max(1e-9, ...values.filter(v => v !== null && Number.isFinite(v)));
  registerScheme({
    id, label,
    residueColors: () => values.map(v => (v === null || !Number.isFinite(v) ? '#b9b7ae' : sampleRamp(ramp, Math.min(1, Math.max(0, v / top))))),
    legend: () => ({ title: title || label, ramp, min, max: max ?? fmt(top, 1) }),
  });
  app.refreshSchemes();
  app.setColorScheme(id);
}

// ------------------------------------------------------------------------------------------- Switch tab
export function experimentsSection(app, ens, r) {
  const sc = r?.scan;
  if (!sc) return null;
  const num = i => ens.residues[i].num;
  return fold(ens, 'experiments', 'Experiments that tell A from B',
    'Suggested experiments for the two states chosen above: double mutants that break a helix of one state and restore it (the classic test that a helix exists and matters), nucleotides whose SHAPE/DMS signal differs, and label positions whose distance differs (FRET, PELDOR). 2D energies from ViennaRNA; distances from the 3D models.',
    st => {
      const key = `${sc.a}-${sc.b}-${r.conditions_text}`;
      if (st.key !== key) { st.key = key; st.result = null; st.error = null; }
      const x = st.result;
      if (!x) return [status(st), st.loading ? null : el('button', { class: 'btn small', onclick: () => run(app, ens, st, 'experiments', { a: sc.a, b: sc.b }) }, icon('flask'), 'Suggest experiments')];
      const mut = m => `${m.from}${num(m.pos)}${m.to}`;
      const rescue = x.rescue.length ? table(
        [{ label: 'Helix of' }, { label: 'First mutation' }, { label: 'A share', num: true, title: 'Share of A in the two-state A/B equilibrium with the first mutation' },
          { label: 'Rescue' }, { label: 'A share', num: true, title: 'Share of A with both mutations' }, { label: '' }],
        x.rescue.map(d => {
          const tr = el('tr', { class: 'clickable', title: `Pair ${ens.sequence[d.pair[0]]}${num(d.pair[0])}–${ens.sequence[d.pair[1]]}${num(d.pair[1])}: ΔΔG ${fmtSigned(d.ddg_first, 1)} kcal/mol, with rescue ${fmtSigned(d.ddg_double, 1)}` },
            el('td', {}, el('span', { class: `sw-key ${d.state === 'A' ? 'a' : 'b'}`, text: d.state })),
            el('td', {}, el('b', { text: mut(d.first) })), el('td', { class: 'num', text: pct(d.frac_a_first) }),
            el('td', {}, el('b', { text: `+ ${mut(d.second)}` })), el('td', { class: 'num', text: pct(d.frac_a_double) }),
            el('td', {}, el('button', { class: 'linkbtn', title: 'Open the start page with the double mutant, ready to predict', text: '3D',
              onclick: e => { e.stopPropagation(); app.predictSequence(d.sequence_double, `${slug(ens.name)}_${mut(d.first)}_${mut(d.second)}`); } })));
          tr.addEventListener('click', () => pick(app, d.pair));
          return tr;
        })) : el('div', { class: 'desc', text: 'No clean disrupt-and-rescue pair found: the two states share their helices, or every disrupting mutation also affects the other state.' });
      const seg = x.probing.segments;
      const probing = seg.length ? el('div', { class: 'chips' }, ...seg.map(s => el('button', { class: 'chip-btn', title: 'Select these nucleotides',
        onclick: () => pick(app, range(s.start, s.end)) }, el('b', { text: span(ens, s.start, s.end) }), ` reactive in ${s.reactive_in}`)))
        : el('div', { class: 'desc', text: 'The two states pair the same nucleotides: chemical probing cannot tell them apart.' });
      const labels = x.labels.length ? table(
        [{ label: 'Label pair' }, { label: 'A (Å)', num: true }, { label: 'B (Å)', num: true }, { label: 'FRET A', num: true }, { label: 'FRET B', num: true }],
        x.labels.map(l => {
          const tr = el('tr', { class: 'clickable' }, el('td', {}, el('b', { text: `${ens.sequence[l.i]}${num(l.i)} – ${ens.sequence[l.j]}${num(l.j)}` })),
            el('td', { class: 'num', text: fmt(l.d_a, 0) }), el('td', { class: 'num', text: fmt(l.d_b, 0) }), el('td', { class: 'num', text: fmt(l.fret_a, 2) }), el('td', { class: 'num', text: fmt(l.fret_b, 2) }));
          tr.addEventListener('click', () => pick(app, [l.i, l.j]));
          return tr;
        })) : el('div', { class: 'desc', text: 'No label pair separates the 3D models of the two states clearly (or a state has no 3D model).' });
      const tsv = () => ['experiment\tdetail\tA share wild type\tA share first\tA share double',
        ...x.rescue.map(d => `rescue (helix of ${d.state})\t${mut(d.first)} + ${mut(d.second)}\t${d.frac_a_wt}\t${d.frac_a_first}\t${d.frac_a_double}`),
        ...seg.map(s => `probing\t${num(s.start)}-${num(s.end)} reactive in ${s.reactive_in}\t\t\t`),
        ...x.labels.map(l => `label pair\t${num(l.i)}-${num(l.j)}: ${l.d_a} Å in A, ${l.d_b} Å in B\t\t\t`)].join('\n');
      return [
        el('div', { class: 'adv-sub', text: 'Disrupt and rescue' }),
        el('div', { class: 'desc', text: `Wild type: ${pct(sc.frac_a_wt)} A. The first mutation should move the balance, the second bring it back.` }), rescue,
        el('div', { class: 'adv-sub', text: 'Chemical probing (SHAPE, DMS)' }),
        el('div', { class: 'desc', text: `${x.probing.informative} of ${x.probing.length} nucleotides are paired in one state and unpaired in the other.` }), probing,
        el('div', { class: 'adv-sub', text: 'Distance labels (FRET, PELDOR)' }),
        el('div', { class: 'desc', text: `C4′–C4′ distances in the conformations of each state; FRET efficiency for R₀ = ${fmt(x.r0, 0)} Å. Linkers and dyes add their own length.` }), labels,
        el('div', { class: 'row adv-actions' }, copyButton('Copy as table', tsv))];
    });
}

export function variantsSection(app, ens, r) {
  if (!r) return null;
  const sc = r.scan;
  const num = i => ens.residues[i].num;
  return fold(ens, 'variants', 'Variants: what does a mutation do?',
    'Paste variants (SNVs from a patient, a design, a screen). For each one: how much the pairing of the whole RNA changes (1 − correlation of the per-nucleotide unpaired probabilities, the riboSNitch measure), how many base pairs of the most stable structure change and, when two states are chosen, the shift of the A/B balance.',
    st => {
      const input = el('input', { class: 'input adv-input', type: 'text', value: st.text || '', placeholder: 'G34U, A12C+G30C', spellcheck: 'false', translate: 'no', 'aria-label': 'Variants' });
      input.addEventListener('input', () => { st.text = input.value; });
      const go = () => { st.text = input.value; run(app, ens, st, 'variants', { text: input.value, a: sc ? sc.a : null, b: sc ? sc.b : null }); };
      input.addEventListener('keydown', e => { if (e.key === 'Enter') go(); });
      const out = [el('div', { class: 'row adv-actions' }, input, el('button', { class: 'btn small', onclick: go }, 'Analyse')), status(st)];
      const x = st.result;
      if (x && !st.loading) {
        const head = [{ label: 'Variant' }, { label: 'Effect' }, { label: 'Change', num: true, title: '1 − correlation of the unpaired probabilities, wild type vs variant (0 = unchanged)' },
          { label: 'Pairs', num: true, title: 'Base pairs of the most stable structure gained or lost' }];
        if (x.two_state) head.push({ label: 'ΔΔG', num: true, title: 'Change of ΔG(A) − ΔG(B), kcal/mol' }, { label: 'A share', num: true });
        head.push({ label: '' });
        out.push(table(head, x.rows.map(v => {
          const tr = el('tr', { class: 'clickable', title: `Largest change around nucleotide ${num(v.hotspot)}` },
            el('td', {}, el('b', { text: v.positions.map(p => `${ens.sequence[p]}${num(p)}${v.sequence[p]}`).join('+'), translate: 'no' })),
            el('td', {}, el('span', { class: `adv-effect ${v.effect}`, text: v.effect })),
            el('td', { class: 'num', text: fmt(v.change, 2) }), el('td', { class: 'num', text: String(v.pairs_changed) }),
            x.two_state ? el('td', { class: 'num', text: fmtSigned(v.ddg, 1) }) : null, x.two_state ? el('td', { class: 'num', text: pct(v.frac_a) }) : null,
            el('td', {}, el('button', { class: 'linkbtn', title: 'Open the start page with this variant, ready to predict', text: '3D',
              onclick: e => { e.stopPropagation(); app.predictSequence(v.sequence, `${slug(ens.name)}_${v.variant}`); } })));
          tr.addEventListener('click', () => pick(app, [...new Set([...v.positions, v.hotspot])]));
          return tr;
        })), el('div', { class: 'desc', text: `Effect: small, moderate (change ≥ 0.1 or ≥ 3 pairs) or large (change ≥ 0.3 or ≥ 8 pairs). ${x.two_state ? `Wild type: ${pct(x.wt.frac_a)} A. ` : ''}A click selects the mutated nucleotide and the place that changes most. Positions are the residue numbers shown in RCV. Substitutions only; modified nucleotides are not modelled.` }),
        // "Not analysed:" and every reason are separate texts: each is translated on its own
        (x.skipped || []).length ? el('div', { class: 'desc', style: 'color: var(--warn)' }, el('span', { text: 'Not analysed:' }), ' ',
          ...x.skipped.flatMap((k, i) => (i ? ['; ', el('span', { text: k.reason })] : [el('span', { text: k.reason })])), '.') : null,
        el('div', { class: 'row adv-actions' }, copyButton('Copy as table', () => ['variant\teffect\tchange\tpairs changed\tddG\tA share\tmost stable structure',
          ...x.rows.map(v => [v.variant, v.effect, v.change, v.pairs_changed, v.ddg ?? '', v.frac_a ?? '', v.mfe_dbn].join('\t'))].join('\n'))));
      }
      return out;
    });
}

// ------------------------------------------------------------------------------------------ Pathway tab
function growthChart(x, width) {
  const n = x.length, rows = x.rows, m = rows.length;
  const left = 34, bottom = 16, cw = Math.max(1, Math.min(6, Math.floor((width - left - 4) / n))), rh = Math.max(1, Math.min(4, Math.floor(260 / m)));
  const W = left + cw * n + 14, H = rh * m + bottom + 4;
  const dpr = window.devicePixelRatio || 1;
  const cv = el('canvas', { class: 'adv-growth', width: Math.round(W * dpr), height: Math.round(H * dpr), style: { width: `${W}px`, height: `${H}px` }, role: 'img',
    'aria-label': 'Structure of the growing transcript: one row per transcript length' });
  const g = cv.getContext('2d');
  g.scale(dpr, dpr);
  const css = getComputedStyle(document.documentElement);
  const text = css.getPropertyValue('--text-3').trim() || '#888', faint = css.getPropertyValue('--border').trim() || '#ddd';
  rows.forEach((r, y) => {
    g.fillStyle = faint;
    g.fillRect(left, 2 + y * rh, cw * r.length, Math.max(1, rh - (rh > 2 ? 1 : 0)));
    for (let i = 0; i < r.code.length; i++) {
      const c = r.code[i];
      if (c === '.') continue;
      g.fillStyle = c === 'f' ? '#2f6fbd' : '#e08a2b';
      g.fillRect(left + i * cw, 2 + y * rh, cw, Math.max(1, rh - (rh > 2 ? 1 : 0)));
    }
  });
  g.fillStyle = text; g.font = '10px sans-serif'; g.textBaseline = 'middle'; g.textAlign = 'right';
  for (const y of [0, Math.floor(m / 2), m - 1]) g.fillText(String(rows[y].length), left - 4, 2 + y * rh + rh / 2);
  g.textAlign = 'center'; g.textBaseline = 'top';
  for (const p of [1, Math.round(n / 2), n]) g.fillText(String(p), left + (p - 0.5) * cw, H - bottom + 4);
  return cv;
}

function formationChart(x, width) {
  const rows = x.rows, W = Math.max(220, Math.min(width, 460)), H = 90, l = 30, b = 16, t = 6, r = 8;
  const X = L => l + (W - l - r) * (L - rows[0].length) / Math.max(1, rows[rows.length - 1].length - rows[0].length), Y = v => t + (H - t - b) * (1 - v);
  const s = svg('svg', { class: 'chart', width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Fraction of the base pairs of states A and B formed at each transcript length' });
  for (const v of [0, 0.5, 1]) {
    s.append(svg('line', { x1: l, x2: W - r, y1: Y(v), y2: Y(v), stroke: 'var(--border)', 'stroke-width': 1 }),
      svg('text', { x: l - 4, y: Y(v) + 3, 'text-anchor': 'end', 'font-size': 10, fill: 'var(--text-3)', text: `${v * 100}%` }));
  }
  for (const [key, col] of [['a', '#1c5cab'], ['b', '#cf5422']]) {
    if (rows[0][key] === undefined) continue;
    s.append(svg('polyline', { points: rows.map(q => `${X(q.length).toFixed(1)},${Y(q[key]).toFixed(1)}`).join(' '), fill: 'none', stroke: col, 'stroke-width': 1.8 }),
      svg('text', { x: W - r, y: Y(rows[rows.length - 1][key]) - 4, 'text-anchor': 'end', 'font-size': 10, fill: col, text: key.toUpperCase() }));
  }
  s.append(svg('text', { x: (l + W - r) / 2, y: H - 3, 'text-anchor': 'middle', 'font-size': 10, fill: 'var(--text-3)', text: 'transcript length (nt)' }));
  return s;
}

export function cotranscriptionSection(app, ens, a = null, b = null) {
  if (!ens) return null;
  const num = i => ens.residues[i].num;
  return fold(ens, 'cotranscription', 'Folding during transcription',
    'The most stable structure of the growing transcript, one nucleotide at a time from the 5′ end. Shows which helices form first, which are only temporary, where the structure rearranges, and which of the two conformations chosen above the nascent RNA resembles: a riboswitch decides while it is being made.',
    st => {
      const key = `${a}-${b}`;
      if (st.key !== key) { st.key = key; st.result = null; st.error = null; }
      const x = st.result;
      if (!x) return [status(st), st.loading ? null : el('button', { class: 'btn small', onclick: () => run(app, ens, st, 'cotranscription', { conf_a: a, conf_b: b }) }, icon('play'), 'Fold the growing transcript')];
      const width = Math.max(240, (app.tabBody ? app.tabBody('pathway').clientWidth : 360) - 28);
      const ev = x.events.slice().sort((p, q) => q.broken - p.broken).slice(0, 6).sort((p, q) => p.length - q.length);
      const line = k => (x[`first_${k}`] ? `${k.toUpperCase()}: 90% of its base pairs are formed at ${x[`first_${k}`]} nt.`
        : x[`peak_${k}`] ? `${k.toUpperCase()}: never more than ${pct(x[`peak_${k}`].fraction)} formed while the RNA is made (at ${x[`peak_${k}`].length} nt).` : null);
      return [
        el('div', { class: 'adv-scroll' }, growthChart(x, width)),
        el('div', { class: 'adv-legend' }, el('i', { style: { background: '#2f6fbd' } }), 'pair kept in the full-length structure', el('i', { style: { background: '#e08a2b' } }), 'temporary pair',
          el('span', { class: 'muted', text: ' · rows: transcript length, columns: nucleotide' })),
        x.rows[0].a !== undefined || x.rows[0].b !== undefined ? formationChart(x, width) : null,
        ...['a', 'b'].map(line).filter(Boolean).map(tx => el('div', { class: 'desc', text: tx })),
        el('div', { class: 'adv-sub', text: 'Rearrangements' }),
        ev.length ? el('div', { class: 'chips' }, ...ev.map(e => el('button', { class: 'chip-btn', title: 'Select the nucleotide whose synthesis triggers the rearrangement',
          onclick: () => pick(app, [e.length - 1]) }, el('b', { text: `${num(e.length - 1)}` }), ` ${e.broken} pairs break, ${e.formed} form`)))
          : el('div', { class: 'desc', text: 'None: every helix that forms during transcription is kept.' }),
        el('div', { class: 'desc', text: `Up to ${x.max_transient} temporary base pairs at one length. Each length is at equilibrium (${x.conditions_text}): this shows what the nascent RNA can form, not how fast; polymerase speed and pausing are not modelled.` })];
    });
}

// ------------------------------------------------------------------------------------------ Details tab
export function motifsSection(app, ens) {
  const k = app.state.conf;
  return fold(ens, 'motifs', 'Structural motifs',
    'Named building blocks of the shown conformation: hairpin loops by type (GNRA, UNCG …), bulges and internal loops with the bend they make, junctions with their coaxial stacks, pseudoknots, non-canonical pairs, base triples and A-minor contacts. A click selects the motif in 3D and 2D.',
    st => {
      if (st.key !== k) { const had = st.key !== undefined && st.result; st.key = k; st.result = null; st.error = null; if (had && !st.loading) queueMicrotask(() => run(app, ens, st, 'motifs', { conformation: k })); }
      const x = st.result;
      if (!x) return [status(st), st.loading ? null : el('button', { class: 'btn small', onclick: () => run(app, ens, st, 'motifs', { conformation: k }) }, icon('layers'), `Find motifs in ${ens.shortId(k)}`)];
      const counts = Object.entries(x.counts).map(([t, c]) => `${c} ${t}${c > 1 ? 's' : ''}`).join(', ');
      return [
        el('div', { class: 'desc', text: x.motifs.length ? `${ens.shortId(k)}: ${x.helices} helices; ${counts}.` : `${ens.shortId(k)}: ${x.helices} helices, no loops or contacts found.` }),
        x.motifs.length ? table([{ label: 'Motif' }, { label: 'Where' }, { label: 'Detail' }], x.motifs.map(m => {
          const res = m.residues, first = res[0], last = res[res.length - 1];
          const tr = el('tr', { class: 'clickable' }, el('td', {}, el('b', { text: m.label, translate: m.type.includes('pair') || m.type.includes('triple') || m.type.includes('A-minor') ? 'no' : null })),
            el('td', { class: 'mono', text: span(ens, first, last) }), el('td', { class: 'adv-detail', text: m.detail }));
          tr.addEventListener('click', () => pick(app, res));
          return tr;
        })) : null,
        el('div', { class: 'row adv-actions' }, x.motifs.length ? copyButton('Copy as table', () => ['type\tmotif\tresidues\tdetail',
          ...x.motifs.map(m => [m.type, m.label, m.residues.map(r => ens.residues[r].num).join(','), m.detail].join('\t'))].join('\n')) : null),
        el('div', { class: 'desc', text: 'Contacts are found from the atom positions (two or more hydrogen bonds between coplanar bases); base edges are named, cis/trans is not. Predicted models have few tertiary contacts; experimental structures show them.' })];
    });
}

export function pocketsSection(app, ens) {
  const k = app.state.conf;
  const others = ens.conformations.map((_c, i) => i).filter(i => i !== k);
  return fold(ens, 'pockets', 'Pockets: where a compound could bind',
    'Open space on the surface of the shown conformation that is enclosed enough to hold a small molecule: its volume, how enclosed it is, and how much of the lining lies outside a Watson-Crick helix — compounds bind loops, bulges and junctions, not the grooves of a helix. On thirty experimental RNA-compound complexes, the site the compound occupies is the top-ranked pocket in 21 and among the pockets found in 24. Choosing a second conformation measures the same sites in both, so a pocket that opens in one and closes in the other is visible. A click selects the pocket in 3D and 2D.',
    (st, draw) => {
      if (st.key !== k) { const had = st.key !== undefined && st.result; st.key = k; st.result = null; st.error = null; if (had) st.rerun = !st.loading; }
      if (st.other !== null && st.other !== undefined && !others.includes(st.other)) st.other = null;
      const body = () => ({ conformation: k, other: st.other === null || st.other === undefined ? null : st.other, burial_min: st.deep ? 6 : 5 });
      const go = () => run(app, ens, st, 'pockets', body());
      if (st.rerun) { st.rerun = false; queueMicrotask(go); }
      const other = select([{ value: '', label: 'nothing' }, ...others.map(i => ({ value: String(i), label: ens.shortId(i) }))],
        st.other === null || st.other === undefined ? '' : String(st.other), { 'aria-label': 'Compare the pockets with this conformation' });
      other.addEventListener('change', () => { st.other = other.value === '' ? null : +other.value; st.result = null; draw(); });
      const depth = select([{ value: 'sites', label: 'sites (enclosed on 5 of 7 sides)' }, { value: 'deep', label: 'only deep pockets (6 of 7)' }],
        st.deep ? 'deep' : 'sites', { 'aria-label': 'How enclosed a site must be' });
      depth.addEventListener('change', () => { st.deep = depth.value === 'deep'; st.result = null; draw(); });
      const out = [el('div', { class: 'row adv-actions' }, el('span', { class: 'muted', text: 'Compare with' }), other,
        el('span', { class: 'muted', text: 'Find' }), depth,
        el('button', { class: 'btn small', onclick: go }, icon('target'), `Find pockets in ${ens.shortId(k)}`)), status(st)];
      const x = st.result;
      if (!x || st.loading) return out;
      const cmp = x.other !== null && x.other !== undefined;
      const nameOther = cmp ? ens.shortId(x.other) : '';
      const vol = v => `${fmt(v, 0)} Å³`;
      const where = p => (p.segments || []).map(s => span(ens, s.start, s.end)).join(', ');
      const openIn = r => el('span', { class: `badge ${r.open_in === 'both' ? '' : r.open_in === 'A' ? 'open-a' : 'open-b'}`,
        title: r.open_in === 'both' ? 'about as open in both conformations' : `more open in ${r.open_in === 'A' ? ens.shortId(k) : nameOther}`,
        text: r.open_in === 'both' ? 'both' : (r.open_in === 'A' ? ens.shortId(k) : nameOther) });
      // while comparing, the panel has no room for every column: 'enclosed' and 'outside a helix' move into the tooltip
      const head = [{ label: 'Pocket', num: true }, { label: 'Volume', num: true, title: 'Grid volume of the pocket, ±10% at this grid spacing' },
        ...(cmp ? [] : [{ label: 'Enclosed', num: true, title: 'Mean share of the seven scan directions in which the structure closes the site' },
          { label: 'Outside a helix', num: true, title: 'Share of the lining nucleotides that are unpaired or non-canonically paired, weighted by contact' }]),
        { label: 'Score', num: true, title: 'Enclosed × outside a helix × size, 0 to 1: a ranking of the sites, not an affinity' },
        ...(cmp ? [{ label: `In ${nameOther}`, num: true, title: `Open volume at the same nucleotides in ${nameOther}` }, { label: 'More open in' }]
          : [{ label: 'Nucleotides' }])];                   // while comparing, the nucleotides are in the table below
      const rows = x.pockets.map(p => {
        const tr = el('tr', { class: 'clickable', title: `${where(p)}; ${p.residues.length} nucleotides; enclosed ${pct(p.burial)}, outside a helix ${pct(p.loop_fraction)}, largest width ${fmt(p.extent, 0)} Å${p.wide ? ' — a groove or cleft rather than a pocket' : ''}` },
          el('td', { class: 'num' }, el('b', { text: String(p.pocket) })),
          el('td', { class: 'num', text: vol(p.volume) }),
          ...(cmp ? [] : [el('td', { class: 'num', text: pct(p.burial) }),
            el('td', { class: 'num', text: p.loop_fraction === null ? '–' : pct(p.loop_fraction) })]),
          el('td', { class: 'num' }, el('b', { text: fmt(p.score, 2) })),
          ...(cmp ? [el('td', { class: 'num', title: `${vol(p.volume_site)} at these nucleotides in ${ens.shortId(k)}`, text: vol(p.volume_other) }), el('td', {}, openIn(p))]
            : [el('td', { class: 'mono adv-detail', text: where(p) })]));
        tr.addEventListener('click', () => pick(app, p.residues));
        return tr;
      });
      const colours = () => {
        const v = new Array(ens.nRes).fill(null);
        for (const p of x.pockets) for (const r of p.residues) v[r] = Math.max(v[r] ?? 0, p.score);
        colourBy(app, 'pockets', 'Pocket score', v, { scale: 1, max: '1', title: 'Score of the best pocket a nucleotide lines' });
      };
      out.push(el('div', { class: 'desc', text: x.n_pockets
        ? `${ens.shortId(k)}: ${x.n_pockets} site${x.n_pockets > 1 ? 's' : ''} of at least ${fmt(x.params.min_volume, 0)} Å³${cmp ? `, ${nameOther}: ${x.n_pockets_other}` : ''}. Grid ${fmt(x.spacing, 2)} Å, probe ${fmt(x.params.probe, 1)} Å.`
        : `No site of at least ${fmt(x.params.min_volume, 0)} Å³ is enclosed enough in ${ens.shortId(k)}. Try a conformation with a formed tertiary structure.` }));
      if (x.n_pockets) out.push(table(head, rows));
      if (cmp && x.changes?.length) {
        const ch = x.changes.filter(c => Math.abs(c.difference) >= 0.25);
        out.push(el('div', { class: 'adv-sub', text: 'Sites that differ between the two conformations' }),
          el('div', { class: 'desc', text: ch.length
            ? `The open volume at the same nucleotides, in ${ens.shortId(k)} and in ${nameOther}. A site that is open in one conformation and closed in the other is what a conformation-selective compound would use.`
            : `Every site holds about the same volume in ${ens.shortId(k)} and in ${nameOther}: nothing here would tell the two apart.` }));
        if (ch.length) out.push(table([{ label: 'Site' }, { label: ens.shortId(k), num: true }, { label: nameOther, num: true },
          { label: 'More open in', title: 'Difference of the two volumes over the larger one, in the tooltip: +100% is open here only, −100% open there only' },
          { label: 'Nucleotides' }],
        ch.map(c => {
          const tr = el('tr', { class: 'clickable', title: `difference ${fmtSigned(100 * c.difference, 0)}% of the larger volume` },
            el('td', {}, el('b', { text: c.pocket && c.other_pocket ? `${c.pocket} / ${c.other_pocket}` : c.pocket ? `${c.pocket} / –` : `– / ${c.other_pocket}` })),
            el('td', { class: 'num', text: vol(c.volume) }), el('td', { class: 'num', text: vol(c.volume_other) }),
            el('td', {}, openIn(c)), el('td', { class: 'mono adv-detail', text: where(c) }));
          tr.addEventListener('click', () => pick(app, c.residues));
          return tr;
        })));
      }
      out.push(el('div', { class: 'row adv-actions' },
        x.n_pockets ? el('button', { class: 'btn small', onclick: colours }, icon('layers'), 'Colour 3D by pocket') : null,
        x.n_pockets ? copyButton('Copy as table', () => [`pocket\tvolume (A^3)\tenclosed\toutside a helix\tscore${cmp ? `\tvolume in ${ens.shortId(x.other)} (A^3)\tmore open in` : ''}\tnucleotides`,
          ...x.pockets.map(p => [p.pocket, p.volume, p.burial, p.loop_fraction, p.score, ...(cmp ? [p.volume_other, p.open_in === 'A' ? ens.shortId(k) : p.open_in === 'B' ? ens.shortId(x.other) : 'both'] : []),
            p.residues.map(r => ens.residues[r].num).join(',')].join('\t'))].join('\n')) : null));
      out.push(el('div', { class: 'desc', text: 'A buriedness scan of the free space around the atoms (the LIGSITE approach), clustered into sites. This is geometry: the score ranks sites, it is not a predicted affinity, and a pocket in a predicted model is only as reliable as the model. Docking of compounds into these pockets is not part of this release.' }));
      return out;
    });
}

export function targetsSection(app, ens) {
  return fold(ens, 'targets', 'Target sites (antisense, probes)',
    'Where an antisense oligonucleotide, siRNA or hybridisation probe can bind: the free energy needed to open every window of the RNA (lower = more accessible) and the share of the conformations in which the window is unpaired.',
    st => {
      const width = select([14, 16, 18, 20, 22, 25].map(w => ({ value: String(w), label: `${w} nt` })), String(st.width || 18), { 'aria-label': 'Window width' });
      width.addEventListener('change', () => { st.width = +width.value; });
      const out = [el('div', { class: 'row adv-actions' }, el('span', { class: 'muted', text: 'Window' }), width,
        el('button', { class: 'btn small', onclick: () => run(app, ens, st, 'targets', { width: +width.value }) }, icon('target'), 'Find accessible sites')), status(st)];
      const x = st.result;
      if (x && !st.loading) {
        out.push(table([{ label: 'Window' }, { label: 'ΔG open', num: true, title: 'Free energy to make the window single-stranded, kcal/mol' },
          { label: 'Open in', num: true, title: 'Share of the conformations in which at least 80% of the window is unpaired' }, { label: 'GC', num: true }, { label: 'Antisense DNA (5′→3′)' }],
        x.best.map(w => {
          const tr = el('tr', { class: 'clickable' }, el('td', {}, el('b', { text: span(ens, w.start, w.end) })),
            el('td', { class: 'num', text: fmt(w.dg_open, 1) }), el('td', { class: 'num', text: pct(w.open_conformations) }), el('td', { class: 'num', text: pct(w.gc) }),
            el('td', { class: 'mono adv-seq', translate: 'no', title: 'Click to copy', text: w.antisense_dna,
              onclick: e => { e.stopPropagation(); navigator.clipboard?.writeText(w.antisense_dna).then(() => toast('Sequence copied.', 'ok')); } }));
          tr.addEventListener('click', () => pick(app, range(w.start, w.end)));
          return tr;
        })),
        el('div', { class: 'row adv-actions' },
          el('button', { class: 'btn small', onclick: () => colourBy(app, 'accessibility', 'Accessibility (unpaired probability)', x.unpaired, { ramp: BLUE_RAMP, scale: 1, max: '100%', title: 'Probability that a nucleotide is unpaired' }) }, icon('layers'), 'Colour 3D by accessibility'),
          copyButton('Copy as table', () => ['start\tend\tdG open (kcal/mol)\topen in conformations\tGC\ttarget (RNA)\tantisense DNA',
            ...x.best.map(w => [ens.residues[w.start].num, ens.residues[w.end].num, w.dg_open, w.open_conformations, w.gc, w.sequence, w.antisense_dna].join('\t'))].join('\n'))),
        el('div', { class: 'desc', text: `Opening energies from the ViennaRNA partition function (${x.conditions_text}). Binding strength, off-targets and chemistry of the oligonucleotide are not part of this ranking.` }));
      }
      return out;
    });
}

function profileChart(values, { color = 'var(--accent)', label, unit = '' } = {}) {
  const n = values.length, W = Math.max(260, Math.min(460, n * 4 + 40)), H = 70, l = 30, b = 14, t = 4;
  const top = Math.max(1e-9, ...values.filter(v => v !== null));
  const s = svg('svg', { class: 'chart', width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': label });
  const bw = (W - l - 4) / n;
  values.forEach((v, i) => { if (v !== null) { const h = (H - t - b) * v / top; s.append(svg('rect', { x: l + i * bw, y: H - b - h, width: Math.max(0.6, bw - 0.4), height: Math.max(0.5, h), fill: color })); } });
  s.append(svg('text', { x: l - 4, y: t + 8, 'text-anchor': 'end', 'font-size': 10, fill: 'var(--text-3)', text: `${fmt(top, top < 2 ? 2 : 0)}${unit}` }),
    svg('text', { x: l - 4, y: H - b, 'text-anchor': 'end', 'font-size': 10, fill: 'var(--text-3)', text: '0' }),
    svg('text', { x: (l + W) / 2, y: H - 2, 'text-anchor': 'middle', 'font-size': 10, fill: 'var(--text-3)', text: label }));
  return s;
}

export function compareSection(app, ens) {
  const others = [...app.ensembles.values()].filter(e => e.id !== ens.id);
  return fold(ens, 'compare', 'Compare with another ensemble',
    'Two ensembles of the same RNA side by side: wild type and mutant, with and without ligand, two prediction methods, a prediction and the experiment. Shows where the shown conformations differ in 3D, where two predictions agree locally (a sign that the region can be trusted), and where the base pairing differs.',
    st => {
      if (!others.length) return el('div', { class: 'desc', text: 'Open a second ensemble of this RNA (File menu, or the start page) to compare.' });
      const sel = select(others.map(e => ({ value: e.id, label: `${e.name} (${e.n})` })), others.some(e => e.id === st.other) ? st.other : others[0].id, { 'aria-label': 'Other ensemble' });
      sel.addEventListener('change', () => { st.other = sel.value; });
      const go = () => {
        const o = app.ensembles.get(sel.value);
        st.other = sel.value; st.otherName = o.name;
        run(app, ens, st, 'compare', { ...(o.serverId ? { other: o.serverId } : { other_doc: o.doc }), ka: app.state.conf, kb: 0 });
      };
      const out = [el('div', { class: 'row adv-actions' }, sel, el('button', { class: 'btn small', onclick: go }, 'Compare')), status(st)];
      const x = st.result;
      if (x && !st.loading) {
        const num = i => ens.residues[i].num;
        const dif = x.differences.length ? `${x.differences.length} sequence difference${x.differences.length > 1 ? 's' : ''} (${x.differences.slice(0, 4).map(d => `${d.a}${num(d.pos_a)}${d.b}`).join(', ')}${x.differences.length > 4 ? ' …' : ''})` : 'same sequence';
        out.push(
          el('dl', { class: 'kv' },
            el('dt', { text: 'Compared' }), el('dd', { text: `${ens.shortId(x.ka)} with the first conformation of ${st.otherName}; ${x.matched} matched nucleotides, ${dif}` }),
            el('dt', { text: 'RMSD (C4′)' }), el('dd', { text: `${fmt(x.rmsd, 1)} Å · TM-score ${fmt(x.tm, 2)}` }),
            el('dt', { text: 'Local agreement' }), el('dd', { title: 'Share of the distances to neighbours within 15 Å that agree within 0.5, 1, 2 and 4 Å (insensitive to domain movements)', text: pct(x.local_agreement_mean) }),
            el('dt', { text: 'Closest pair' }), el('dd', { text: `${fmt(x.cross_rmsd_min, 1)} Å between any two conformations of the ensembles (mean ${fmt(x.cross_rmsd_mean, 1)} Å)` })),
          profileChart(x.displacement, { label: 'displacement per nucleotide', unit: ' Å' }),
          x.moved.length ? el('div', { class: 'chips' }, ...x.moved.map(m => el('button', { class: 'chip-btn', title: 'Select the region that moved', onclick: () => pick(app, range(m.start, m.end)) },
            el('b', { text: span(ens, m.start, m.end) }), ` moved up to ${fmt(m.max, 0)} Å`))) : el('div', { class: 'desc', text: 'No region stands out: the two conformations differ evenly or not at all.' }),
          x.pairing_changed.length ? el('div', { class: 'chips' }, ...x.pairing_changed.map(m => el('button', { class: 'chip-btn', title: 'Select the region whose base pairing differs', onclick: () => pick(app, range(m.start, m.end)) },
            el('b', { text: span(ens, m.start, m.end) }), ' pairing differs'))) : null,
          el('div', { class: 'row adv-actions' },
            el('button', { class: 'btn small', onclick: () => colourBy(app, 'compare-displacement', 'Displacement (comparison)', x.displacement, { title: `Displacement against ${st.otherName}`, max: `${fmt(Math.max(...x.displacement.filter(v => v !== null)), 0)} Å` }) }, icon('layers'), 'Colour 3D by displacement'),
            el('button', { class: 'btn small ghost', onclick: () => colourBy(app, 'compare-agreement', 'Local agreement (comparison)', x.local_agreement, { ramp: BLUE_RAMP, scale: 1, max: '100%', title: `Local agreement with ${st.otherName}` }) }, 'Colour by local agreement')),
          el('div', { class: 'desc', text: 'The superposition follows the part that stays the same, so a moving domain shows as displaced. High local agreement between two independent predictions marks regions both methods build the same way.' }));
      }
      return out;
    });
}
