// SAXS profiles, FRET distances and cryo-EM maps for the open ensemble (Details tab ▸ closed section "Experimental
// data", Tools ▸ Load SAXS, FRET or cryo-EM data…). The same path as the probing data (probing.js): every
// conformation's agreement with the data, and a maximum-entropy reweighting of the ensemble to fit all loaded data
// together (rcv/expdata.py, POST /api/expdata/parse and /api/expdata/fit). Data, agreement and weights are saved in
// .rcv files; a map is kept by the local engine only while RCV runs (load it again to refit after reopening).
import { el, svg, fmt } from './util.js';
import { icon } from './icons.js';
import { dialog, field, select, toast } from './ui.js';
import { fold, status, table } from './panels/insight.js';

const KINDS = [{ value: 'saxs', label: 'SAXS profile' }, { value: 'fret', label: 'FRET distances' }, { value: 'map', label: 'Cryo-EM map' }];
const KIND_LABEL = Object.fromEntries(KINDS.map(k => [k.value, k.label]));
const METRIC = {
  saxs: { key: 'saxs_chi2', label: 'SAXS χ²', lower_is_better: true, description: 'Reduced χ² of the SAXS profile computed from this conformation (Debye formula, scale and background fitted) against the measured one.' },
  fret: { key: 'fret_chi2', label: 'FRET χ²', lower_is_better: true, description: 'Reduced χ² of the FRET distances or efficiencies of this conformation against the measured ones.' },
  map: { key: 'map_cc', label: 'Map correlation', lower_is_better: false, description: 'Correlation of the density simulated from this conformation, placed in the map, with the cryo-EM map above its contour level.' },
};
const HINT = {
  saxs: 'q and I(q) with errors, one point per line (.dat from ATSAS, BioXTAS RAW, …). q in 1/Å or 1/nm.',
  fret: 'One pair per line: "i j distance [sd]" in Å, or "i j E=0.45 R0=54" for an efficiency. Residue numbers 1…n; the C4′ atoms unless atom=… is given.',
  map: 'An MRC/CCP4 density map of this RNA (.mrc, .map). Each conformation is placed in the map and its simulated density compared above the contour level.',
};

export class ExpData {
  constructor(app) {
    this.app = app;
    app.registerSessionPart('expdata', {
      save: ens => (ens?.expdata?.sets?.length ? { sets: ens.expdata.sets, fit: ens.expdata.fit, theta: ens.expdata.theta ?? null } : null),
      load: (ens, v) => { if (v?.sets?.length) { ens.expdata = { sets: v.sets, fit: v.fit || null, theta: v.theta ?? null }; this._metrics(ens); } else if (!v && ens.expdata) this.clear(ens, true); },
    });
    app.on('ensemble', ens => { if (ens?.expdata) this._metrics(ens); });
  }

  dialog(kind0 = 'saxs') {
    const app = this.app, ens = app.ens;
    if (!ens) return;
    if (!app.api.online) { toast('Reading experimental data needs the local RCV engine.', 'error', 'Offline'); return; }
    const kind = select(KINDS, kind0);
    const area = el('textarea', { class: 'input mono', rows: 6, 'aria-label': 'Data', spellcheck: 'false', style: { width: '100%', boxSizing: 'border-box' } });
    const res = el('input', { class: 'input', type: 'number', min: 1, step: 'any', placeholder: 'auto', style: { width: '90px' }, 'aria-label': 'Resolution (Å)' });
    const cont = el('input', { class: 'input', type: 'number', step: 'any', placeholder: 'mean + 2 sd', style: { width: '120px' }, 'aria-label': 'Contour level' });
    const mapRow = el('div', { class: 'row', style: { gap: '10px', alignItems: 'end' } }, field('Resolution (Å)', res), field('Contour level', cont));
    const inp = el('input', { type: 'file', hidden: true });
    const fileBtn = el('button', { class: 'btn small', type: 'button', onclick: () => inp.click() }, icon('open'), 'Choose a file…');
    const hint = el('div', { class: 'desc' });
    const preview = el('div', { class: 'desc' });
    let name = '', mapData = null, parsed = null;
    const sync = () => {
      const k = kind.value;
      hint.textContent = HINT[k];
      area.hidden = k === 'map';
      mapRow.hidden = k !== 'map';
      inp.accept = k === 'map' ? '.mrc,.map,.ccp4,.mrcs' : '.dat,.txt,.csv,.tsv,.out,text/plain';
      area.placeholder = k === 'saxs' ? '0.0100  1.52e+03  12.1\n0.0125  1.50e+03  11.8\n…' : '12  48  52.5  4\n3  60  E=0.42  R0=54\n…';
      parsed = null;
      preview.replaceChildren();
    };
    inp.addEventListener('change', async () => {
      const f = inp.files[0];
      if (!f) return;
      name = f.name;
      if (kind.value === 'map') {
        mapData = await new Promise((ok, bad) => { const r = new FileReader(); r.onload = () => ok(r.result); r.onerror = bad; r.readAsDataURL(f); });
        preview.textContent = `${f.name} (${fmt(f.size / 1048576, 1)} MB)`;
      } else area.value = await f.text();
      parsed = null;
    });
    kind.addEventListener('change', sync);
    sync();
    const parse = async () => {
      const body = { kind: kind.value, name: name || (kind.value === 'map' ? 'map' : 'pasted data'), n: ens.nRes };
      if (kind.value === 'map') { if (!mapData) throw new Error('Choose a map file first.'); Object.assign(body, { data: mapData, resolution: res.value || null, contour: cont.value === '' ? null : +cont.value }); }
      else body.text = area.value;
      return app.api.post('/api/expdata/parse', body);
    };
    const d = dialog({
      title: 'Load SAXS, FRET or cryo-EM data', width: 600,
      body: b => b.append(
        el('div', { class: 'desc', text: `Data for ${ens.name} (${ens.nRes} nt, ${ens.n} conformations). Each data type replaces the earlier one of its kind; all loaded data are fitted together.` }),
        el('div', { class: 'row', style: { gap: '10px', alignItems: 'end' } }, field('Data type', kind), fileBtn, inp),
        hint, area, mapRow, preview),
      actions: [{ label: 'Cancel', close: true }, { label: 'Load', primary: true, action: async () => {
        try { parsed = parsed || await parse(); } catch (err) { preview.textContent = err.message; return; }
        preview.textContent = `${parsed.n_points} points${parsed.q_range ? ` · q ${parsed.q_range[0]}–${parsed.q_range[1]} 1/Å` : ''}${(parsed.notes || []).map(x => ` · ${x}`).join('')}`;
        d.close();
        await this.add(ens, parsed);
      } }],
    });
  }

  async add(ens, data) {
    const x = ens.expdata || (ens.expdata = { sets: [], fit: null, theta: null });
    x.sets = x.sets.filter(s => s.kind !== data.kind).concat([data]);
    await this.refit(ens);
  }

  async refit(ens = this.app.ens, theta) {
    const x = ens?.expdata;
    if (!x?.sets?.length) return;
    if (theta !== undefined) x.theta = theta;
    const app = this.app;
    x.loading = true; x.error = null; app.emit('expdata', ens);
    try {
      const body = { datasets: x.sets, theta: x.theta };
      if (ens.serverId) body.document = ens.serverId; else body.doc = ens.doc;
      x.fit = await app.api.post('/api/expdata/fit', body);
      const rw = x.fit.reweighting;
      toast(`${x.sets.map(s => KIND_LABEL[s.kind]).join(', ')}: ${rw ? `reweighting χ² ${rw.chi2_prior.map(v => fmt(v, 2)).join(', ')} → ${rw.chi2.map(v => fmt(v, 2)).join(', ')}` : 'agreement computed'}.`, 'ok', 'Experimental data', 7000);
    } catch (err) { x.error = err.message; toast(err.message, 'error', 'Experimental data'); }
    x.loading = false;
    this._metrics(ens);
    app.emit('values');
    app.emit('expdata', ens);
  }

  remove(ens, kind) {
    const x = ens?.expdata;
    if (!x) return;
    x.sets = x.sets.filter(s => s.kind !== kind);
    this.app.dropMetric(ens, METRIC[kind].key);
    if (!x.sets.length) { this.clear(ens); return; }
    this.refit(ens);
  }

  clear(ens = this.app.ens, quiet = false) {
    if (!ens?.expdata) return;
    ens.expdata = null;
    for (const m of Object.values(METRIC)) this.app.dropMetric(ens, m.key);
    this.app.dropMetric(ens, 'exp_weight');
    if (!quiet) this.app.emit('values');
    this.app.emit('expdata', ens);
  }

  _metrics(ens) {
    const x = ens.expdata;
    if (!x?.fit) return;
    const add = (key, def, vals) => {
      for (let k = 0; k < ens.n; k++) ens.conformations[k].values[key] = vals?.[k] ?? null;
      const m = { key, group: 'data', unit: '', ...def };
      if (ens.metricByKey[key]) Object.assign(ens.metricByKey[key], m); else { ens.metrics.push(m); ens.metricByKey[key] = m; }
    };
    for (const d of x.fit.datasets) add(METRIC[d.kind].key, METRIC[d.kind], d.agreement);
    if (x.fit.reweighting) add('exp_weight', { label: 'Data weight (SAXS/FRET/EM)', unit: '%', lower_is_better: false, description: 'Weight of this conformation after maximum-entropy reweighting of the ensemble to fit the loaded SAXS, FRET and cryo-EM data together (prior: uniform).' }, x.fit.reweighting.weights.map(w => 100 * w));
  }
}

// log I(q): measured (points) and computed from the shown conformation (line)
function saxsChart(set, curve, W) {
  const H = 110, ml = 30, mr = 6, mt = 6, mb = 18;
  const s = svg('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}`, width: '100%', height: H, role: 'img', 'aria-label': 'SAXS profile' });
  const q = set.q, I = set.I;
  const pos = I.filter(v => v > 0);
  const lo = Math.log10(Math.min(...pos)), hi = Math.log10(Math.max(...pos));
  const x = v => ml + (v - q[0]) / Math.max(1e-9, q[q.length - 1] - q[0]) * (W - ml - mr);
  const y = v => Math.min(H - mb, Math.max(mt, mt + (hi - Math.log10(Math.max(v, 1e-30))) / Math.max(1e-9, hi - lo) * (H - mt - mb)));
  s.append(svg('line', { x1: ml, x2: W - mr, y1: H - mb + 0.5, y2: H - mb + 0.5, stroke: 'var(--axis)' }));
  q.forEach((v, i) => s.append(svg('circle', { cx: x(v), cy: y(I[i]), r: 1.4, fill: 'var(--text-3)' })));
  if (curve) s.append(svg('polyline', { points: q.flatMap((v, i) => (curve[i] > 0 ? [`${x(v)},${y(curve[i])}`] : [])).join(' '), fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 1.6 }));
  s.append(svg('text', { x: ml, y: H - 3, 'font-size': 9, fill: 'var(--text-3)', text: `q ${fmt(q[0], 3)}–${fmt(q[q.length - 1], 3)} 1/Å · log I · line: this conformation` }));
  return s;
}

export function expdataSection(app, ens) {
  const k = app.state.conf;
  return fold(ens, 'expdata', 'Experimental data: SAXS, FRET, cryo-EM',
    'Load a SAXS profile, FRET distances or a cryo-EM map to see how well each conformation agrees with them and to reweight the ensemble so that its average fits all the data together (maximum entropy).',
    () => {
      const x = ens.expdata;
      const load = el('div', { class: 'row adv-actions' }, ...KINDS.map(kd => el('button', { class: 'btn small', onclick: () => app.expdata.dialog(kd.value) }, icon('upload'), `${kd.label}…`)));
      if (!x?.sets?.length) return [load];
      const st = { loading: x.loading, error: x.error };
      const fit = x.fit;
      const rows = x.sets.map(set => {
        const d = fit?.datasets?.find(z => z.kind === set.kind);
        const ag = d?.agreement || [];
        const best = ag.reduce((b, v, j) => (v !== null && (b === null || (d.lower_is_better ? v < ag[b] : v > ag[b])) ? j : b), null);
        return el('tr', {},
          el('td', {}, el('b', { text: KIND_LABEL[set.kind] })), el('td', { class: 'adv-detail', text: set.name, translate: 'no' }),
          el('td', { class: 'num', text: d ? fmt(ag[k], set.kind === 'map' ? 3 : 2) : '–' }),
          el('td', { class: 'num', text: best !== null ? ens.shortId(best) : '–' }),
          el('td', {}, el('button', { class: 'ss-btn', title: 'Remove these data', 'aria-label': 'Remove these data', onclick: () => app.expdata.remove(ens, set.kind) }, icon('x'))));
      });
      const out = [status(st),
        table([{ label: 'Data' }, { label: 'File' }, { label: `${ens.shortId(k)}`, num: true, title: 'Agreement of the shown conformation: χ² (lower is better) or map correlation (higher is better)' }, { label: 'Best', num: true }, { label: '' }], rows)];
      const saxs = x.sets.find(s => s.kind === 'saxs');
      const sd = fit?.datasets?.find(z => z.kind === 'saxs');
      if (saxs) out.push(saxsChart(saxs, sd?.curves?.[k], Math.max(240, (app.details?.root?.clientWidth || 300) - 40)));
      const rw = fit?.reweighting;
      if (rw) {
        const th = el('input', { class: 'input', type: 'number', min: 0.1, step: 'any', value: rw.theta, style: { width: '90px' }, 'aria-label': 'theta' });
        out.push(el('dl', { class: 'kv' },
          el('dt', { text: 'Reweighting' }), el('dd', { text: `χ²/n ${rw.chi2_prior.map(v => fmt(v, 2)).join(' · ')} → ${rw.chi2.map(v => fmt(v, 2)).join(' · ')}` }),
          el('dt', { text: 'Effective conformations' }), el('dd', { text: `${fmt(rw.n_eff_conformations, 1)} of ${ens.n} · weight of ${ens.shortId(k)} ${fmt(100 * rw.weights[k], 1)} %` })),
        el('div', { class: 'desc', text: 'θ (how far the data may move the weights; smaller = more):' }),
        el('div', { class: 'row', style: { gap: '6px', alignItems: 'center' } }, th,
          el('button', { class: 'btn small ghost', onclick: () => app.expdata.refit(ens, +th.value || null) }, 'Refit')));
      }
      out.push(el('div', { class: 'row adv-actions' },
        ...x.sets.map(s => (fit ? el('button', { class: 'btn small', onclick: () => app.setMetric(METRIC[s.kind].key) }, icon('chart'), `Plot ${METRIC[s.kind].label}`) : null)),
        rw ? el('button', { class: 'btn small', onclick: () => app.setMetric('exp_weight') }, icon('chart'), 'Plot data weights') : null,
        el('button', { class: 'btn small ghost', onclick: () => app.expdata.clear(ens) }, icon('trash'), 'Remove all')),
        load,
        el('div', { class: 'desc', text: 'SAXS: Debye formula over the heavy atoms with excluded volume, no hydration shell (judge the fit at q < 0.25 1/Å). FRET: distances between the named atoms; dye linkers are not modelled. Maps: rigid placement of each conformation, density simulated at the map resolution.' }));
      return out;
    });
}
