// Switch tab: the ensemble's secondary-structure states with Boltzmann populations (ViennaRNA, adjustable
// temperature) and a single-mutation scan showing which mutations push the equilibrium towards state A or B.
// Links to 3D: colour by mutation sensitivity, select the mutated residue, predict the mutant in 3D.
import { el, svg, fmt, fmtSigned, clamp } from '../util.js';
import { icon } from '../icons.js';
import { tooltip, tipRow, toast, select } from '../ui.js';
import { registerScheme, sampleRamp, BLUE_RAMP, ORANGE_RAMP } from '../colors.js';

const SAT = 4.0;            // kcal/mol at which the diverging colour saturates
const NEUTRAL_LIGHT = '#e7e5de', NEUTRAL_DARK = '#3a3a37';

export function divergingColor(ddg, dark = false) {
  if (ddg === null || ddg === undefined || !Number.isFinite(ddg)) return dark ? '#2a2a28' : '#ffffff';
  const t = clamp(Math.abs(ddg) / SAT, 0, 1);
  const neutral = dark ? NEUTRAL_DARK : NEUTRAL_LIGHT;
  if (t < 0.04) return neutral;
  const ramp = ddg < 0 ? BLUE_RAMP.slice(2) : ORANGE_RAMP.slice(1);
  return sampleRamp([neutral, ...ramp], t);
}

export class SwitchPanel {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    this.temperature = 37;
    this.result = null;
    this.a = null; this.b = null;
    this.picked = null;
    app.on('ensemble', () => { this.result = null; this.a = this.b = null; this.picked = null; if (this.visible()) this.load(); else this.render(); });
    app.on('conformation', () => this._markCurrent());
    app.on('tab', id => { if (id === 'switch' && !this.result && app.ens) this.load(); });
    new ResizeObserver(() => { const w = root.clientWidth; if (this.result && Math.abs(w - (this._w || 0)) > 8) { this._w = w; this.render(); } }).observe(root);
    this.render();
  }

  visible() { return !this.root.hidden; }

  stateOf(k) { return this.result ? this.result.states.findIndex(s => s.conformations.includes(k)) : -1; }

  async load() {
    const app = this.app, ens = app.ens;
    if (!ens) { this.render(); return; }
    if (!app.api.online) { this.error = 'The switch analysis needs the local RCV engine.'; this.render(); return; }
    if (!app.info?.capabilities?.vienna) { this.error = 'The switch analysis needs ViennaRNA, which is not installed.'; this.render(); return; }
    this.error = null;
    this.loading = true;
    this.render();
    try {
      const body = { temperature: this.temperature, a: this.a, b: this.b };
      if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
      let r = await app.api.post('/api/switch', body);
      if ((this.a === null || this.b === null) && r.states.length > 1) {
        // default comparison: previous vs shown conformation if they differ in 2D, else the two most populated states
        const sk = r.states.findIndex(s => s.conformations.includes(app.state.conf));
        const sp = app.state.prev !== null ? r.states.findIndex(s => s.conformations.includes(app.state.prev)) : -1;
        const byPop = r.states.map((s, i) => [s.relative, i]).sort((x, y) => y[0] - x[0]).map(x => x[1]);
        [this.a, this.b] = sp >= 0 && sk >= 0 && sp !== sk ? [sp, sk] : [byPop[0], byPop[1]];
        r = await app.api.post('/api/switch', { ...body, a: this.a, b: this.b });
      }
      this.result = r;
      this.picked = null;
    } catch (err) { this.error = err.message; } finally { this.loading = false; }
    this.render();
  }

  render() {
    const app = this.app, ens = app.ens, root = this.root;
    root.replaceChildren();
    if (!ens) { root.append(el('div', { class: 'muted', text: 'No ensemble loaded.' })); return; }
    const head = el('div', { class: 'section' },
      el('div', { class: 'section-title' }, 'Secondary-structure states', el('span', { class: 'grow' }),
        el('label', { class: 'sw-temp', title: 'Temperature for the Turner energies and populations' }, 'T',
          (() => { const i = el('input', { class: 'input', type: 'number', min: 0, max: 95, step: 1, value: this.temperature, 'aria-label': 'Temperature in °C' });
            i.addEventListener('change', () => { this.temperature = clamp(+i.value || 37, 0, 95); this.load(); }); return i; })(), '°C')),
      el('div', { class: 'desc', text: 'Distinct 2D structures of the ensemble’s conformations with ViennaRNA free energies and populations (Boltzmann). A 2D-level estimate: the 3D model scores are not used.' }));
    root.append(head);
    if (this.loading) { head.append(el('div', { class: 'row' }, el('span', { class: 'spinner small' }), ' computing…')); return; }
    if (this.error) { head.append(el('div', { class: 'warnbox', text: this.error })); return; }
    const r = this.result;
    if (!r) { head.append(el('button', { class: 'btn', onclick: () => this.load() }, 'Analyse')); return; }
    // ---- states table
    const tb = el('tbody');
    const maxRel = Math.max(...r.states.map(s => s.relative));
    r.states.forEach((s, i) => {
      const confs = s.conformations;
      tb.append(el('tr', { class: `sw-state ${confs.includes(app.state.conf) ? 'on' : ''}`, 'data-s': i },
        el('td', {}, el('b', { text: `S${i + 1}` })),
        el('td', { class: 'mono sw-dbn', title: s.dbn, text: s.dbn }),
        el('td', { class: 'num', text: fmt(s.energy, 1) }),
        el('td', { class: 'num', title: 'Probability of exactly this structure in the full ViennaRNA ensemble', text: pct(s.p_ensemble) }),
        el('td', {}, el('div', { class: 'sw-bar', title: `${pct(s.relative)} of the population among these states` }, el('i', { style: { width: `${(s.relative / maxRel) * 100}%` } }), el('span', { text: pct(s.relative) }))),
        el('td', { class: 'num' }, confs.length ? el('button', { class: 'linkbtn sw-conf', title: `${confs.map(k => ens.shortId(k)).join(', ')} — click to show ${ens.shortId(confs[0])}`, text: String(confs.length), onclick: () => app.goTo(confs.includes(app.state.conf) ? confs[(confs.indexOf(app.state.conf) + 1) % confs.length] : confs[0]) }) : el('span', { class: 'muted', text: '–' }))));
    });
    head.append(el('div', { class: 'wl-table sw-tablebox' }, el('table', { class: 'table sw-table' },
      el('thead', {}, el('tr', {}, el('th', { text: '' }), el('th', { text: '2D' }), el('th', { class: 'num', text: 'ΔG kcal/mol' }), el('th', { class: 'num', text: 'P ensemble' }), el('th', { text: 'Share' }), el('th', { class: 'num', text: 'Models', title: 'Conformations with this 2D structure (click to step through them)' }))), tb)),
      el('div', { class: 'desc', text: `MFE ${fmt(r.mfe.energy, 1)} kcal/mol (P ${pct(r.mfe_probability)})${r.states.some(s => s.dbn === r.mfe.dbn) ? ' — present in the ensemble' : ' — not among the ensemble’s states'} · ensemble free energy ${fmt(r.ensemble_energy, 1)} kcal/mol · ${r.temperature} °C · ViennaRNA ${r.vienna || ''}` }));
    // ---- mutation scan
    const scanSec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Mutations that shift the switch'));
    root.append(scanSec);
    if (r.states.length < 2) {
      scanSec.append(el('div', { class: 'desc', text: 'All conformations share one 2D state, so there is no 2D switch to shift. Explore new conformations with the Pull tool, or open an ensemble with several 2D states.' }));
      return;
    }
    const opts = r.states.map((s, i) => ({ value: String(i), label: `S${i + 1} · ${fmt(s.energy, 1)} kcal/mol` }));
    const selA = select(opts, String(this.a)), selB = select(opts, String(this.b));
    const reload = () => { this.a = +selA.value; this.b = +selB.value; if (this.a !== this.b) this.load(); else toast('Choose two different states.', 'info'); };
    selA.addEventListener('change', reload); selB.addEventListener('change', reload);
    scanSec.append(el('div', { class: 'sw-ab' }, el('span', { class: 'sw-key a', text: 'A' }), selA, el('span', { class: 'sw-key b', text: 'B' }), selB));
    const sc = r.scan;
    if (!sc) return;
    scanSec.append(el('div', { class: 'desc', text: `Wild type: ΔG(A) − ΔG(B) = ${fmtSigned(sc.ddg_wt, 1)} kcal/mol → ${pct(sc.frac_a_wt)} A in the two-state A/B equilibrium. Each cell shows ΔΔG of a point mutation: blue favours A, orange favours B.` }));
    scanSec.append(this._heatmap(sc));
    this.detail = el('div', { class: 'sw-detail' });
    scanSec.append(this.detail);
    this._detail();
    const list = (title, items) => el('div', { class: 'sw-list' }, el('div', { class: 'section-title', style: { marginTop: '10px' } }, title),
      items.length ? el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'Mutation' }), el('th', { class: 'num', text: 'ΔΔG' }), el('th', { class: 'num', text: 'A share' }), el('th', { class: 'num', text: 'P(A)' }), el('th', { class: 'num', text: 'P(B)' }))),
        el('tbody', {}, ...items.map(x => { const tr = el('tr', { class: 'clickable' }, el('td', {}, el('b', { text: `${x.from}${ens.residues[x.pos].num}${x.to}` })), el('td', { class: 'num', text: fmtSigned(x.ddg, 1) }),
          el('td', { class: 'num', text: pct(x.frac_a) }), el('td', { class: 'num', text: pct(x.p_a) }), el('td', { class: 'num', text: pct(x.p_b) })); tr.addEventListener('click', () => this.pick(x.pos, x.to)); return tr; })))
        : el('div', { class: 'desc', text: 'none' }));
    scanSec.append(list('Strongest shifts towards A', sc.towards_a), list('Strongest shifts towards B', sc.towards_b),
      el('div', { class: 'row', style: { marginTop: '10px', gap: '8px', flexWrap: 'wrap' } },
        el('button', { class: 'btn small', onclick: () => this.colorBySensitivity() }, icon('layers'), 'Colour 3D by mutation sensitivity'),
        el('button', { class: 'btn small ghost', onclick: () => this.copyTable() }, icon('copy'), 'Copy scan as table')),
      el('div', { class: 'desc', text: 'P(A), P(B): probabilities of the two structures in the full ViennaRNA ensemble of the mutant (other structures can dominate). A share: two-state A/B equilibrium only. A mutation that breaks a pair removes it from that state.' }));
  }

  _heatmap(sc) {
    const app = this.app, ens = app.ens, n = ens.nRes;
    const dark = document.documentElement.dataset.theme === 'dark';
    const W = Math.max(240, this.root.clientWidth - 28);
    const lab = 16, cw = Math.max(4, Math.min(14, (W - lab) / n)), gap = cw >= 7 ? 1 : 0, ch = 14, H = 4 * ch + 30;
    const width = lab + cw * n;
    const s = svg('svg', { class: 'chart sw-heat', width, height: H, viewBox: `0 0 ${width} ${H}`, role: 'img', 'aria-label': 'Mutation scan heat map' });
    'ACGU'.split('').forEach((b, row) => s.append(svg('text', { x: 2, y: 12 + row * ch + ch / 2, 'font-size': 10, fill: 'var(--text-2)', text: b })));
    for (let i = 0; i < n; i++) {
      const wt = ens.sequence[i];
      'ACGU'.split('').forEach((b, row) => {
        const x = lab + i * cw, y = 6 + row * ch;
        if (b === wt) {
          s.append(svg('rect', { x: x + gap / 2, y: y + gap / 2, width: cw - gap, height: ch - gap, fill: 'var(--surface)', stroke: 'var(--axis)', 'stroke-width': 0.5 }));
          if (cw >= 6) s.append(svg('circle', { cx: x + cw / 2, cy: y + ch / 2, r: 1.4, fill: 'var(--text-3)' }));
          return;
        }
        const v = sc.matrix[i][b];
        const rect = svg('rect', { x: x + gap / 2, y: y + gap / 2, width: cw - gap, height: ch - gap, rx: cw >= 8 ? 1.5 : 0, fill: divergingColor(v.ddg, dark), class: 'sw-cell', 'data-i': i, 'data-b': b });
        if (this.picked && this.picked.pos === i && this.picked.to === b) rect.setAttribute('stroke', 'var(--text)'), rect.setAttribute('stroke-width', 1.6);
        s.append(rect);
      });
      if (i === 0 || (i + 1) % 10 === 0) s.append(svg('text', { x: lab + i * cw + cw / 2, y: H - 4, 'text-anchor': 'middle', 'font-size': 9, fill: 'var(--text-3)', text: ens.residues[i].num }));
    }
    s.addEventListener('pointermove', e => {
      const c = e.target.closest('.sw-cell');
      if (!c) { tooltip.hide(); app.hoverResidues([]); return; }
      const i = +c.dataset.i, b = c.dataset.b, v = sc.matrix[i][b];
      app.hoverResidues([i]);
      tooltip.show([tipRow(`${ens.sequence[i]}${ens.residues[i].num}${b}`, `ΔΔG ${fmtSigned(v.ddg, 1)} kcal/mol`),
        el('div', { class: 'tl', text: `${v.ddg < 0 ? 'favours A' : v.ddg > 0 ? 'favours B' : 'neutral'} · A share ${pct(sc.frac_a_wt)} → ${pct(v.frac_a)}${v.breaks_a ? ' · breaks a pair of A' : ''}${v.breaks_b ? ' · breaks a pair of B' : ''}` })], e.clientX, e.clientY);
    });
    s.addEventListener('pointerleave', () => { tooltip.hide(); app.hoverResidues([]); });
    s.addEventListener('click', e => { const c = e.target.closest('.sw-cell'); if (c) this.pick(+c.dataset.i, c.dataset.b); });
    const legend = el('div', { class: 'sw-legend' }, el('span', { text: 'favours A' }),
      el('span', { class: 'ramp', style: { background: `linear-gradient(90deg, ${[-SAT, -SAT / 2, 0, SAT / 2, SAT].map(v => divergingColor(v, dark)).join(',')})` } }),
      el('span', { text: 'favours B' }), el('span', { class: 'muted', text: `±${SAT} kcal/mol` }));
    return el('div', {}, el('div', { class: 'sw-heatbox' }, s), legend);
  }

  pick(pos, to) {
    this.picked = { pos, to };
    this.app.select([pos], 'set');
    this.render();
  }

  _detail() {
    const app = this.app, ens = app.ens, box = this.detail, sc = this.result?.scan;
    if (!box || !sc) return;
    if (!this.picked) { box.replaceChildren(el('div', { class: 'desc', text: 'Click a cell for details and to predict the mutant in 3D.' })); return; }
    const { pos, to } = this.picked;
    const v = sc.matrix[pos][to];
    const seq = ens.sequence.slice(0, pos) + to + ens.sequence.slice(pos + 1);
    const name = `${ens.sequence[pos]}${ens.residues[pos].num}${to}`;
    box.replaceChildren(el('div', { class: 'sw-card' },
      el('div', {}, el('b', { text: name }), ` · ΔΔG ${fmtSigned(v.ddg, 1)} kcal/mol · A share ${pct(sc.frac_a_wt)} → ${pct(v.frac_a)}`),
      el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '6px' } },
        el('button', { class: 'btn small primary', onclick: () => app.predictSequence(seq, `${slugName(ens.name)}_${name}`) , title: 'Open the start page with the mutant sequence, ready to predict' }, icon('flask'), 'Predict in 3D'),
        el('button', { class: 'btn small', onclick: () => { app.select([pos], 'set'); app.centerOnSelection(); } }, icon('target'), 'Show residue'),
        el('button', { class: 'btn small ghost', onclick: () => navigator.clipboard?.writeText(`>${name}\n${seq}\n`).then(() => toast('Mutant sequence copied (FASTA).', 'ok')) }, icon('copy'), 'Copy FASTA'))));
  }

  _markCurrent() {
    const k = this.app.state.conf;
    for (const tr of this.root.querySelectorAll('.sw-state')) tr.classList.toggle('on', !!this.result?.states[+tr.dataset.s]?.conformations.includes(k));
  }

  colorBySensitivity() {
    const sc = this.result?.scan;
    if (!sc) return;
    const sens = sc.matrix.map(row => Math.max(0, ...Object.values(row).map(v => Math.abs(v.ddg))));
    const max = Math.max(1e-6, ...sens);
    const a = this.result.states[sc.a], b = this.result.states[sc.b];
    registerScheme({
      id: 'switch-sensitivity', label: 'Mutation sensitivity (switch)',
      residueColors: () => sens.map(v => sampleRamp(ORANGE_RAMP, v / max)),
      legend: () => ({ title: `Largest |ΔΔG| of a point mutation between S${sc.a + 1} and S${sc.b + 1}`, ramp: ORANGE_RAMP, min: '0', max: `${fmt(max, 1)} kcal/mol` }),
    });
    this.app.refreshSchemes();
    this.app.setColorScheme('switch-sensitivity');
    toast(`Residues coloured by how strongly a point mutation shifts the S${sc.a + 1}/S${sc.b + 1} balance (${a.conformations.length + b.conformations.length} conformations).`, 'info', 'Colour scheme');
  }

  copyTable() {
    const ens = this.app.ens, sc = this.result?.scan;
    if (!sc) return;
    const rows = ['position\twild type\tmutant\tddG (kcal/mol)\tA share\tbreaks A\tbreaks B'];
    sc.matrix.forEach((row, i) => { for (const [b, v] of Object.entries(row)) rows.push([ens.residues[i].num, ens.sequence[i], b, v.ddg, v.frac_a, v.breaks_a, v.breaks_b].join('\t')); });
    navigator.clipboard?.writeText(rows.join('\n')).then(() => toast(`${rows.length - 1} mutations copied as a tab-separated table.`, 'ok'));
  }
}

function pct(x) {
  if (x === null || x === undefined || !Number.isFinite(x)) return '–';
  if (x >= 0.995) return '>99%';
  if (x < 0.001) return '<0.1%';
  return x < 0.1 ? `${(x * 100).toFixed(1)}%` : `${Math.round(x * 100)}%`;
}
function slugName(s) { return String(s || 'rna').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 40) || 'rna'; }
