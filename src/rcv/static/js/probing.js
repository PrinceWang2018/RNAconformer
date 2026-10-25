// Chemical probing data (SHAPE or DMS reactivities) for the open ensemble: load a profile, colour the structure by
// reactivity, score how well each conformation agrees with the data (AUC of reactivity vs unpaired nucleotides),
// and reweight the ensemble to fit the data (maximum entropy, computed by rcv/probing.py). The profile, the fit and
// the weights are saved in .rcv files. Data-informed 2D populations: Conditions dialog (Switch tab).
import { el, svg, fmt } from './util.js';
import { icon } from './icons.js';
import { dialog, field, select, toast, tooltip, tipRow } from './ui.js';
import { registerScheme, sampleRamp, ORANGE_RAMP } from './colors.js';

const MISSING = '#c4c2ba';
export const reactivityColor = v => (v === null || v === undefined ? MISSING : sampleRamp(['#f4f2ec', ...ORANGE_RAMP], Math.min(1, Math.max(0, v / 1.5))));

export class Probing {
  constructor(app) {
    this.app = app;
    app.registerSessionPart('probing', {
      save: ens => (ens?.probing ? { ...ens.probing } : null),
      load: (ens, v) => { if (v && Array.isArray(v.values) && v.values.length === ens.nRes) this._install(ens, v, { quiet: true }); else if (!v && ens.probing) this.clear(ens, { quiet: true }); },
    });
    app.on('ensemble', ens => { if (ens?.probing) this._metrics(ens); });
  }

  // ---------------------------------------------------------------- loading
  dialog() {
    const app = this.app, ens = app.ens;
    if (!ens) return;
    if (!app.api.online) { toast('Reading probing data needs the local RCV engine.', 'error', 'Offline'); return; }
    let text = '', name = '';
    const kind = select([{ value: 'SHAPE', label: 'SHAPE (1M7, NMIA, NAI, …: all four nucleotides)' }, { value: 'DMS', label: 'DMS (A and C)' }], ens.probing?.kind || 'SHAPE');
    const norm = el('input', { type: 'checkbox', checked: true });
    const area = el('textarea', { class: 'input mono', rows: 6, placeholder: '1\t0.12\n2\t0.85\n3\t-999\n…  (position and reactivity; .shape, .map, CSV with a reactivity column, or one value per line)', 'aria-label': 'Reactivity data' });
    const fileBtn = el('button', { class: 'btn small', type: 'button' }, icon('open'), 'Choose a file…');
    const preview = el('div', { class: 'desc' });
    const chart = el('div');
    const inp = el('input', { type: 'file', accept: '.txt,.shape,.map,.csv,.tsv,.dat,.rdat,text/plain', hidden: true });
    inp.addEventListener('change', async () => { const f = inp.files[0]; if (!f) return; name = f.name; area.value = text = await f.text(); check(); });
    fileBtn.addEventListener('click', () => inp.click());
    let parsed = null;
    const check = async () => {
      text = area.value;
      parsed = null;
      if (!text.trim()) { preview.replaceChildren(); chart.replaceChildren(); return; }
      try {
        parsed = await app.api.post('/api/probing/parse', { text, sequence: ens.sequence, kind: kind.value, normalise: norm.checked });
        const vals = norm.checked ? parsed.normalised : parsed.values;
        preview.replaceChildren(el('span', { text: `${parsed.n_values} of ${ens.nRes} nucleotides with data` }), el('span', { text: ' · ' }), el('span', { text: parsed.format }),
          parsed.factor ? el('span', { text: ' · ' }) : '', parsed.factor ? el('span', { text: `normalised by ${fmt(parsed.factor, 3)} (2–8 % rule)` }) : '',
          ...(parsed.notes || []).map(x => el('span', { text: ` · ${x}` })));
        chart.replaceChildren(profileChart(ens, vals, Math.max(240, 520)));
      } catch (err) { preview.replaceChildren(err.message); chart.replaceChildren(); }
    };
    area.addEventListener('input', () => { clearTimeout(this._t); this._t = setTimeout(check, 300); });
    kind.addEventListener('change', check); norm.addEventListener('change', check);
    const d = dialog({
      title: 'Load probing data', width: 600,
      body: b => b.append(
        el('div', { class: 'desc', text: `Reactivities for ${ens.name} (${ens.nRes} nt), position 1 = the first nucleotide of the sequence.` }),
        el('div', { class: 'row', style: { gap: '10px', alignItems: 'end' } }, field('Reagent', kind), el('label', { class: 'check' }, norm, 'Normalise (2–8 % rule)'), fileBtn, inp),
        area, preview, chart),
      actions: [{ label: 'Cancel', close: true }, { label: 'Load', primary: true, action: async () => {
        if (!parsed) await check();
        if (!parsed || !parsed.n_values) { toast('No reactivity values found.', 'error', 'Probing data'); return; }
        const values = norm.checked ? parsed.normalised : parsed.values;
        await this.set(ens, { name: name || 'pasted data', kind: kind.value, values, raw: parsed.values, factor: parsed.factor || null, format: parsed.format });
        d.close();
      } }],
    });
  }

  async set(ens, data) {
    const app = this.app;
    const p = { ...data, theta: data.theta ?? null, fit: null };
    try {
      const body = { values: p.values, kind: p.kind, theta: p.theta };
      if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
      p.fit = await app.api.post('/api/probing/fit', body);
    } catch (err) { toast(err.message, 'error', 'Probing data'); }
    this._install(ens, p);
    if (p.fit) {
      const w = p.fit.reweighting;
      toast(`${p.kind} data for ${p.values.filter(v => v !== null).length} nucleotides. Agreement (AUC) ${fmt(Math.min(...p.fit.agreement.filter(x => x !== null)), 2)}–${fmt(Math.max(...p.fit.agreement.filter(x => x !== null)), 2)}${w ? `; reweighting χ² ${fmt(w.chi2_prior, 2)} → ${fmt(w.chi2, 2)}` : ''}.`, 'ok', 'Probing data', 7000);
    }
  }

  async refit(theta) {
    const ens = this.app.ens;
    if (!ens?.probing) return;
    await this.set(ens, { ...ens.probing, theta });
  }

  _install(ens, p, { quiet = false } = {}) {
    ens.probing = p;
    this._metrics(ens);
    const app = this.app;
    registerScheme({
      id: 'reactivity', label: 'Reactivity (probing data)',
      residueColors: e => (e.probing ? e.probing.values.map(reactivityColor) : e.residues.map(() => MISSING)),
      legend: e => ({ title: `${e.probing?.kind || 'Probing'} reactivity (normalised)`, ramp: ['#f4f2ec', ...ORANGE_RAMP], min: '0', max: '≥ 1.5' }),
    });
    app.refreshSchemes();
    if (!quiet) app.emit('values');
    app.emit('probing', ens);
  }

  // per-conformation metrics: agreement (AUC) and data weight
  _metrics(ens) {
    const p = ens.probing;
    const add = (key, def, vals) => {
      for (let k = 0; k < ens.n; k++) ens.conformations[k].values[key] = vals?.[k] ?? null;
      const m = { key, group: 'data', ...def };
      if (ens.metricByKey[key]) Object.assign(ens.metricByKey[key], m); else { ens.metrics.push(m); ens.metricByKey[key] = m; }
    };
    if (!p?.fit) return;
    add('probe_auc', { label: `${p.kind} agreement (AUC)`, unit: '', lower_is_better: false,
      description: `Area under the ROC curve of the ${p.kind} reactivities as a predictor of unpaired nucleotides in this conformation (1 = perfect, 0.5 = no information).` }, p.fit.agreement);
    if (p.fit.reweighting) {
      add('probe_weight', { label: `${p.kind} data weight`, unit: '%', lower_is_better: false,
        description: `Weight of this conformation after maximum-entropy reweighting of the ensemble to fit the ${p.kind} data (prior: uniform).` }, p.fit.reweighting.weights.map(w => 100 * w));
    }
  }

  clear(ens = this.app.ens, { quiet = false } = {}) {
    if (!ens?.probing) return;
    ens.probing = null;
    for (const key of ['probe_auc', 'probe_weight']) this.app.dropMetric(ens, key);
    if (ens.conditions) ens.conditions.useProbing = false;
    if (this.app.state.colorScheme === 'reactivity' && this.app.ens === ens) this.app.setColorScheme('nucleotide');
    if (!quiet) this.app.emit('values');
    this.app.emit('probing', ens);
  }

  // ---------------------------------------------------------------- Details tab section
  renderSection() {
    const app = this.app, ens = app.ens, p = ens?.probing;
    const box = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Probing data', el('span', { class: 'grow' }),
      el('button', { class: 'btn small ghost', title: 'Load a SHAPE or DMS reactivity profile', onclick: () => this.dialog() }, icon('upload'), p ? 'Replace' : 'Load')));
    if (!p) {
      box.append(el('div', { class: 'desc', text: 'Load SHAPE or DMS reactivities to colour the structure by reactivity, see how well each conformation agrees with the data and reweight the ensemble to fit them.' }));
      return box;
    }
    const k = app.state.conf;
    const W = Math.max(240, (app.details?.root?.clientWidth || 300) - 28);
    box.append(el('div', { class: 'desc' }, el('span', { text: `${p.kind} · ${p.name} · ` }), el('span', { text: `${p.values.filter(v => v !== null).length} of ${ens.nRes} nucleotides` }),
      p.factor ? el('span', { text: ' · ' }) : null, p.factor ? el('span', { text: `normalised (factor ${fmt(p.factor, 3)})` }) : null),
      profileChart(ens, p.values, W, ens.conformations[k]));
    const fit = p.fit;
    if (fit) {
      const auc = fit.agreement[k];
      const best = fit.agreement.reduce((b, v, j) => (v !== null && (b === null || v > fit.agreement[b]) ? j : b), null);
      const rw = fit.reweighting;
      box.append(el('dl', { class: 'kv' },
        el('dt', { text: 'Agreement (AUC)' }), el('dd', { text: `${auc === null ? '–' : fmt(auc, 2)} for ${ens.shortId(k)} · best ${best !== null ? `${ens.shortId(best)} (${fmt(fit.agreement[best], 2)})` : '–'}` }),
        rw ? el('dt', { text: 'Reweighting' }) : null,
        rw ? el('dd', { text: `χ²/n ${fmt(rw.chi2_prior, 2)} → ${fmt(rw.chi2, 2)} · effective conformations ${fmt(rw.n_eff_conformations, 1)} of ${ens.n} · weight of ${ens.shortId(k)} ${fmt(100 * rw.weights[k], 1)} %` }) : null));
      if (rw) {
        const th = el('input', { class: 'input', type: 'number', min: 0.1, step: 'any', value: rw.theta, style: { width: '90px' }, 'aria-label': 'theta' });
        box.append(el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '6px', alignItems: 'center' } },
          el('span', { class: 'muted', text: 'θ (how far the data may move the weights; smaller = more):' }), th,
          el('button', { class: 'btn small ghost', onclick: () => this.refit(+th.value || null) }, 'Refit')));
      }
    }
    box.append(el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginTop: '8px' } },
      el('button', { class: 'btn small', onclick: () => app.setColorScheme('reactivity') }, icon('palette'), 'Colour by reactivity'),
      fit ? el('button', { class: 'btn small', onclick: () => app.setMetric('probe_auc') }, icon('chart'), 'Plot agreement') : null,
      fit?.reweighting ? el('button', { class: 'btn small', onclick: () => app.setMetric('probe_weight') }, icon('chart'), 'Plot data weights') : null,
      el('button', { class: 'btn small ghost', onclick: () => this.clear() }, icon('trash'), 'Remove')),
      el('div', { class: 'desc', text: 'Agreement: reactivity as a predictor of unpaired nucleotides (base pairs annotated in 3D). Reweighting: maximum-entropy fit of the ensemble average (unpaired = reactive) to the data; weights near the prior mean the data do not distinguish the conformations.' }));
    return box;
  }
}

// Bar chart of a reactivity profile; with a conformation, paired nucleotides are marked under the axis.
export function profileChart(ens, vals, W = 300, conf = null) {
  const n = ens.nRes, H = 86, ml = 26, mr = 4, mt = 6, mb = 18;
  const s = svg('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}`, width: '100%', height: H, role: 'img', 'aria-label': 'Reactivity profile' });
  const bw = (W - ml - mr) / n;
  const max = Math.max(1.5, ...vals.filter(v => v !== null && v !== undefined).map(v => Math.min(v, 3)));
  const y = v => H - mb - Math.min(v, max) / max * (H - mt - mb);
  s.append(svg('line', { x1: ml, x2: W - mr, y1: H - mb + 0.5, y2: H - mb + 0.5, stroke: 'var(--axis)' }));
  for (const t of [0, 1]) s.append(svg('text', { x: ml - 4, y: y(t) + 3, 'text-anchor': 'end', 'font-size': 9, fill: 'var(--text-3)', text: String(t) }));
  const paired = new Set();
  if (conf) for (const p of conf.pairs || []) { paired.add(p[0]); paired.add(p[1]); }
  vals.forEach((v, i) => {
    const x = ml + i * bw;
    if (v === null || v === undefined) s.append(svg('rect', { x: x + bw * 0.15, y: H - mb - 2, width: Math.max(0.6, bw * 0.7), height: 2, fill: MISSING }));
    else s.append(svg('rect', { x: x + bw * 0.1, y: y(Math.max(0, v)), width: Math.max(0.6, bw * 0.8), height: Math.max(0.5, H - mb - y(Math.max(0, v))), fill: reactivityColor(v) }));
    if (conf) s.append(svg('rect', { x: x + bw * 0.1, y: H - mb + 3, width: Math.max(0.6, bw * 0.8), height: 3, fill: paired.has(i) ? 'var(--series-1)' : 'var(--grid)' }));
  });
  s.append(svg('text', { x: ml, y: H - 2, 'font-size': 9, fill: 'var(--text-3)', text: conf ? 'blue = paired in this conformation' : `1–${n}` }));
  s.addEventListener('pointermove', e => {
    const r = s.getBoundingClientRect();
    const i = Math.floor(((e.clientX - r.left) * W / r.width - ml) / bw);
    if (i < 0 || i >= n) { tooltip.hide(); return; }
    tooltip.show([tipRow(vals[i] === null || vals[i] === undefined ? 'no data' : fmt(vals[i], 2), ens.resLabel(i)), conf ? el('div', { class: 'tl', text: paired.has(i) ? 'paired' : 'unpaired' }) : null], e.clientX, e.clientY);
  });
  s.addEventListener('pointerleave', () => tooltip.hide());
  return s;
}
