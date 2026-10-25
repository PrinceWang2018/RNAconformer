// Details tab: conformation facts, full metric table (table view of the charts), selection, ensemble info.
import { el, fmt, fmtVal, fmtSigned } from '../util.js';
import { icon } from '../icons.js';
import { direction } from './energy.js';
import { askText } from '../ui.js';

export class DetailsPanel {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    for (const ev of ['ensemble', 'conformation', 'metric', 'selection', 'values', 'order', 'measurements', 'named', 'scenes', 'probing']) app.on(ev, () => this.render());
  }

  render() {
    const app = this.app, ens = app.ens, root = this.root;
    root.replaceChildren();
    if (!ens) { root.append(el('div', { class: 'muted', text: 'No ensemble loaded.' })); return; }
    const k = app.state.conf, p = app.state.prev;
    const c = ens.conformations[k];
    // ---- selection
    const sel = [...app.state.selection].sort((a, b) => a - b);
    const selBox = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Selection', el('span', { class: 'grow' }),
      sel.length ? el('button', { class: 'btn small ghost', text: 'Clear', onclick: () => app.select([], 'set') }) : null));
    if (!sel.length) selBox.append(el('div', { class: 'desc', text: 'Click residues in the 3D view to select them (Shift+click adds). Selected residues are shown as atoms and labelled.' }));
    else {
      const chips = el('div', { class: 'chips' });
      for (const r of sel) chips.append(el('span', { class: 'chip' }, ens.resLabel(r), el('button', { title: 'Remove', 'aria-label': `Remove ${ens.resLabel(r)}`, onclick: () => app.select([r], 'toggle') }, icon('x'))));
      selBox.append(chips);
      const t = el('table', { class: 'table', style: { marginTop: '8px' } },
        el('thead', {}, el('tr', {}, el('th', { text: 'Residue' }), el('th', { text: 'Paired with (current)' }), el('th', { class: 'num', text: 'RMSF Å' }))));
      const tb = el('tbody');
      for (const r of sel.slice(0, 30)) {
        const pi = ens.pairInfo(k, r);
        tb.append(el('tr', {}, el('td', { text: ens.resLabel(r) }), el('td', { text: pi ? `${ens.resLabel(pi.partner)}${pi.type ? ' (' + pi.type[0] + '–' + pi.type[1] + ')' : ''}${pi.pk ? ' PK' : ''}` : 'unpaired' }),
          el('td', { class: 'num', text: fmt(ens.rmsf[r], 1) })));
      }
      t.append(tb);
      selBox.append(t, el('div', { class: 'row', style: { marginTop: '8px' } },
        el('button', { class: 'btn small', text: 'Center on selection', onclick: () => app.centerOnSelection() }),
        el('button', { class: 'btn small', text: 'Superpose on selection', title: 'Re-fit every conformation on the selected residues (C4′ atoms)', onclick: () => app.superposeOnSelection() })));
    }
    // ---- conformation
    const val = c.validation;
    const conf = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Conformation'),
      el('dl', { class: 'kv' },
        el('dt', { text: 'ID' }), el('dd', { text: c.id }),
        c.label !== c.id ? el('dt', { text: 'Label' }) : null, c.label !== c.id ? el('dd', { text: c.label }) : null,
        el('dt', { text: 'Source' }), el('dd', { text: [c.source?.file, c.source?.model ? `model ${c.source.model}` : null, c.source?.chain ? `chain ${c.source.chain}` : null].filter(Boolean).join(' · ') || '–' }),
        el('dt', { text: 'Rank' }), el('dd', { text: String(c.rank ?? '–') }),
        c.cluster ? el('dt', { text: 'Cluster' }) : null, c.cluster ? el('dd', { text: String(c.cluster) }) : null,
        el('dt', { text: 'Validation' }), el('dd', {}, val ? el('span', { class: val.passed ? 'state-ok' : 'state-bad' }, val.passed ? 'passed' : 'issues: ' + Object.entries(val.issue_counts || {}).map(([a, b]) => `${a} ${b}`).join(', ')) : '–',
          val ? el('span', { class: 'muted', text: ` · ${val.severe_clashes ?? 0} severe clashes · ${fmt(val.overlaps_per_1000 ?? null, 1)} overlaps/1000 atoms` }) : null),
        val?.sugar_puckers ? el('dt', { text: 'Sugar puckers' }) : null,
        val?.sugar_puckers ? el('dd', { text: `C3′-endo ${val.sugar_puckers.C3endo}, C2′-endo ${val.sugar_puckers.C2endo}, other ${val.sugar_puckers.other}` }) : null),
      el('div', { class: 'section-title', style: { marginTop: '10px' } }, '2D structure found in 3D'),
      el('div', { class: 'seq' }, alignDbn(ens.sequence, c.dbn)),
      c.input_dbn ? el('div', { class: 'section-title', style: { marginTop: '8px' } }, 'Input 2D hypothesis' + (c.secondary_structure_id ? ` (${c.secondary_structure_id})` : '')) : null,
      c.input_dbn ? el('div', { class: 'seq', text: c.input_dbn }) : null);
    // ---- values table (current vs previous)
    const tbl = el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'Metric' }), el('th', { class: 'num', text: ens.shortId(k) }),
      el('th', { class: 'num', text: p !== null && p !== k ? ens.shortId(p) : 'previous' }), el('th', { class: 'num', text: 'Δ' }))));
    const tb = el('tbody');
    for (const m of ens.metrics) {
      const v = ens.value(k, m.key);
      if (v === null && ens.conformations.every((_, j) => ens.value(j, m.key) === null)) continue;
      const pv = p !== null && p !== k ? ens.value(p, m.key) : null;
      const d = v !== null && pv !== null ? v - pv : null;
      const dir = direction(m, d);
      tb.append(el('tr', { class: m.key === app.state.metric ? 'on' : '' },
        el('td', { text: m.label + (m.unit ? ` (${m.unit})` : ''), title: m.description || '' }),
        el('td', { class: 'num', text: fmtVal(v, m.unit) }), el('td', { class: 'num', text: fmtVal(pv, m.unit) }),
        el('td', { class: 'num' }, el('span', { class: dir === 'good' ? 'state-ok' : dir === 'bad' ? 'state-bad' : 'muted', text: d === null ? '–' : fmtSigned(d) }))));
      tb.lastChild.style.cursor = 'pointer';
      tb.lastChild.addEventListener('click', () => app.setMetric(m.key));
    }
    tbl.append(tb);
    const values = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'All metrics', el('span', { class: 'grow' }), el('span', { class: 'note', text: 'click a row to plot it' })), tbl);
    // ---- all conformations (table view of the charts)
    const key = app.state.metric, mk = ens.metricByKey[key] || {};
    const all = el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: '#' }), el('th', { text: 'Conformation' }),
      el('th', { class: 'num', text: mk.label || key }), el('th', { class: 'num', text: 'RMSD to current' }), el('th', { class: 'num', text: 'Pairs' }))));
    const tb2 = el('tbody');
    app.order.forEach((j, i) => {
      const tr = el('tr', { class: j === k ? 'on' : '' }, el('td', { text: String(i + 1) }), el('td', { text: ens.shortId(j) }),
        el('td', { class: 'num', text: fmtVal(ens.value(j, key), mk.unit) }), el('td', { class: 'num', text: fmt(ens.rmsd[k][j]) }),
        el('td', { class: 'num', text: String(ens.conformations[j].pairs.length) }));
      tr.style.cursor = 'pointer';
      tr.addEventListener('click', () => app.goTo(j));
      tb2.append(tr);
    });
    all.append(tb2);
    const allSec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'All conformations'), all);
    // ---- ensemble
    const d = ens.doc;
    const rms = ens.rmsd.flat().filter(v => v > 0);
    const prov = d.provenance;
    const ensSec = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Ensemble'),
      ...(d.warnings || []).map(w => el('div', { class: 'warnbox', text: w })),
      el('dl', { class: 'kv' },
        el('dt', { text: 'Name' }), el('dd', { text: ens.name }),
        el('dt', { text: 'Origin' }), el('dd', { text: d.origin === 'rnafold3d-engine' ? `rnafold3d prediction (${prov?.software?.engine_version ? 'engine ' + prov.software.engine_version : 'engine'})` : d.origin === 'rnafold3d-package' ? 'rnafold3d ensemble package' : 'structure files' }),
        el('dt', { text: 'Length' }), el('dd', { text: `${ens.nRes} nt, ${ens.nAtoms} heavy atoms` }),
        el('dt', { text: 'Conformations' }), el('dd', { text: String(ens.n) }),
        rms.length ? el('dt', { text: 'Pairwise RMSD' }) : null, rms.length ? el('dd', { text: `${fmt(Math.min(...rms))} – ${fmt(Math.max(...rms))} Å (C4′)` }) : null,
        el('dt', { text: 'Superposition' }), el('dd', { text: d.ensemble.superposition }),
        prov?.secondary_structure_mode ? el('dt', { text: '2D source' }) : null, prov?.secondary_structure_mode ? el('dd', { text: prov.secondary_structure_mode.replaceAll('_', ' ') }) : null,
        prov?.parameters?.trajectories ? el('dt', { text: 'Sampling' }) : null, prov?.parameters?.trajectories ? el('dd', { text: `${prov.parameters.trajectories} trajectories, seed ${prov.parameters.seed}` }) : null),
      el('div', { class: 'section-title', style: { marginTop: '10px' } }, 'Sequence'),
      el('div', { class: 'seq', text: ens.sequence }),
      ...(prov?.notes || []).map(n => el('div', { class: 'desc', text: n })));
    root.append(selBox, this._objects(), app.measure ? app.measure.renderSection() : null, app.probing ? app.probing.renderSection() : null, conf, values, allSec, ensSec);
  }
}

// Named selections and saved scenes, with quick actions (the console offers the same and more).
DetailsPanel.prototype._objects = function () {
  const app = this.app, ens = app.ens;
  const d = app.displayOf();
  const box = el('div', { class: 'section' }, el('div', { class: 'section-title' }, 'Named selections & scenes', el('span', { class: 'grow' }),
    el('button', { class: 'btn small ghost', title: 'Store the current selection under a name', disabled: app.state.selection.size ? null : true, onclick: () => {
      askText('Name the selection', 'Name', `sel${d.named.size + 1}`, { pattern: /^[A-Za-z_][\w-]*$/, hint: 'Letters, digits, _ and - (start with a letter).' })
        .then(name => { if (name) app.console.run(`select ${name}, sele`); });
    } }, icon('plus'), 'Selection'),
    el('button', { class: 'btn small ghost', title: 'Save the current view (camera, conformation, colours, selection) as a scene', onclick: () => app.console.run(`view save scene${app.scenes.list().length + 1}`) }, icon('bookmark'), 'Scene')));
  const named = [...d.named.entries()];
  const scenes = app.scenes.list();
  if (!named.length && !scenes.length) {
    box.append(el('div', { class: 'desc', text: 'Store selections and views to come back to them, e.g. in the console: select loop, resi 23-27 · color orange, loop · view save closeup. Both are saved in .rcv files.' }));
    return box;
  }
  const chips = el('div', { class: 'chips obj-chips' });
  for (const [name, set] of named) {
    chips.append(el('span', { class: 'chip', title: `${set.size} residues` },
      el('button', { class: 'chip-main', text: name, onclick: () => app.select([...set], 'set') }),
      el('button', { title: 'Zoom', 'aria-label': `Zoom on ${name}`, onclick: () => app.console.run(`zoom ${name}`, { echo: false }) }, icon('target')),
      el('button', { title: 'Delete', 'aria-label': `Delete ${name}`, onclick: () => app.console.run(`delete ${name}`, { echo: false }) }, icon('x'))));
  }
  for (const s of scenes) {
    chips.append(el('span', { class: 'chip scene' },
      el('button', { class: 'chip-main', onclick: () => app.scenes.recall(s.name) }, icon('bookmark'), s.name),
      el('button', { title: 'Delete', 'aria-label': `Delete scene ${s.name}`, onclick: () => app.scenes.remove(s.name) }, icon('x'))));
  }
  box.append(chips);
  return box;
};

function alignDbn(seq, dbn) {
  const wrap = el('div');
  for (let i = 0; i < seq.length; i += 40) {
    wrap.append(el('div', { text: seq.slice(i, i + 40) }), el('div', { style: { color: 'var(--text)', marginBottom: '4px' }, text: (dbn || '').slice(i, i + 40) }));
  }
  return wrap;
}
