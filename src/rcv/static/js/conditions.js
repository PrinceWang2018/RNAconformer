// Solution conditions of the secondary-structure analyses (Switch and Pathway tabs): temperature, Na+, a ligand
// that binds an aptamer motif, and (Switch tab) reactivity pseudo-energies from probing data. Kept per ensemble and
// saved in .rcv files; changing them emits 'conditions' so the analyses are computed again.
import { el, fmt } from './util.js';
import { dialog, field, select, toast } from './ui.js';

export const STANDARD_SALT = 1.021;            // mol/L, standard state of the Turner parameters (ViennaRNA default)
export const LIGANDS = {
  theophylline: { name: 'Theophylline aptamer', motif_seq: 'GAUACCAG&CCCUUGGCAGC', motif_struct: '(...((((&)...)))...)', kd_uM: 0.32 },
};

export function defaultConditions() {
  return { temperature: 37, sodium: STANDARD_SALT, ligand: null, useProbing: false };
}

export function conditionsOf(app, ens = app.ens) {
  if (!ens) return defaultConditions();
  if (!ens.conditions) ens.conditions = defaultConditions();
  return ens.conditions;
}

// JSON for the analysis requests.
export function requestConditions(app, ens = app.ens, { probing = true } = {}) {
  const c = conditionsOf(app, ens);
  const out = { temperature: c.temperature, sodium: c.sodium, ligand: c.ligand || null };
  if (probing && c.useProbing && ens?.probing?.values) out.probing = { kind: ens.probing.kind, values: ens.probing.values };
  return out;
}

export function conditionsText(c, ens = null) {
  const parts = [`${fmt(c.temperature, 0)} °C`, Math.abs(c.sodium - STANDARD_SALT) < 1e-6 ? '1.021 M Na⁺ (standard)' : `${fmt(c.sodium * 1000, 0)} mM Na⁺`];
  if (c.ligand) parts.push(`${c.ligand.name || 'ligand'} ${fmt(c.ligand.conc_uM ?? 0, 2)} µM`);
  if (c.useProbing && ens?.probing) parts.push(`${ens.probing.kind} data`);
  return parts.join(' · ');
}

// The same summary as separate text pieces (translatable one by one), e.g. for buttons and labels.
export function conditionsParts(c, ens = null) {
  const out = [el('span', { text: `${fmt(c.temperature, 0)} °C` })];
  const sep = () => el('span', { text: ' · ' });
  if (Math.abs(c.sodium - STANDARD_SALT) < 1e-6) out.push(sep(), el('span', { text: '1.021 M Na⁺' }), el('span', { text: ' ' }), el('span', { text: '(standard)' }));
  else out.push(sep(), el('span', { text: `${fmt(c.sodium * 1000, 0)} mM Na⁺` }));
  if (c.ligand) out.push(sep(), el('span', { text: c.ligand.name || 'ligand' }), el('span', { text: ` ${fmt(c.ligand.conc_uM ?? 0, 2)} µM` }));
  if (c.useProbing && ens?.probing) out.push(sep(), el('span', { text: `${ens.probing.kind}` }), el('span', { text: ' ' }), el('span', { text: 'data' }));
  return out;
}

export function setConditions(app, patch, ens = app.ens) {
  if (!ens) return;
  ens.conditions = { ...conditionsOf(app, ens), ...patch };
  app.emit('conditions', ens.conditions);
}

export function conditionsDialog(app) {
  const ens = app.ens;
  if (!ens) return;
  const c = { ...conditionsOf(app, ens) };
  const temp = el('input', { class: 'input', type: 'number', min: 0, max: 95, step: 1, value: c.temperature, 'aria-label': 'Temperature in °C' });
  const na = el('input', { class: 'input', type: 'number', min: 5, max: 3000, step: 1, value: Math.round(c.sodium * 1000), 'aria-label': 'Na+ concentration in mM' });
  const ligSel = select([{ value: '', label: 'No ligand' }, ...Object.entries(LIGANDS).map(([k, v]) => ({ value: k, label: v.name })), { value: 'custom', label: 'Custom aptamer…' }],
    c.ligand ? (c.ligand.preset || 'custom') : '');
  const lig = c.ligand || {};
  const mseq = el('input', { class: 'input mono', value: lig.motif_seq || '', placeholder: 'GAUACCAG&CCCUUGGCAGC', 'aria-label': 'Aptamer motif sequence' });
  const mst = el('input', { class: 'input mono', value: lig.motif_struct || '', placeholder: '(...((((&)...)))...)', 'aria-label': 'Aptamer motif structure' });
  const kd = el('input', { class: 'input', type: 'number', min: 0.000001, step: 'any', value: lig.kd_uM ?? 1, 'aria-label': 'Kd in µM' });
  const conc = el('input', { class: 'input', type: 'number', min: 0, step: 'any', value: lig.conc_uM ?? 10, 'aria-label': 'Ligand concentration in µM' });
  const probe = el('input', { type: 'checkbox', checked: c.useProbing ? true : null, disabled: ens.probing ? null : true });
  const custom = el('div', { class: 'cond-custom' },
    field('Aptamer motif (sequence; & splits an interior loop)', mseq), field('Motif structure (dot-bracket, same length)', mst));
  const ligBox = el('div', {}, custom, el('div', { class: 'row', style: { gap: '10px' } }, field('Kd (µM)', kd), field('Concentration (µM)', conc)));
  const sync = () => {
    const v = ligSel.value;
    ligBox.hidden = !v;
    custom.hidden = v !== 'custom';
    if (v && v !== 'custom') { const p = LIGANDS[v]; kd.value = c.ligand?.preset === v ? (c.ligand.kd_uM ?? p.kd_uM) : p.kd_uM; }
  };
  ligSel.addEventListener('change', sync);
  sync();
  const body = el('div', { class: 'cond-dialog' },
    el('div', { class: 'desc', text: 'Conditions of the secondary-structure energies (ViennaRNA, Turner 2004) used by the Switch and Pathway tabs.' }),
    el('div', { class: 'row', style: { gap: '10px' } }, field('Temperature (°C)', temp), field('Na⁺ (mM; 1021 = standard)', na)),
    field('Ligand', ligSel, 'A ligand stabilises every structure that forms its aptamer (a hairpin or interior loop) by −RT ln(1 + c/Kd).'),
    ligBox,
    el('label', { class: 'row', style: { gap: '6px', marginTop: '8px' } }, probe, el('span', { text: ens.probing ? `Use the ${ens.probing.kind} reactivities as pseudo-energies for the populations (Switch tab)` : 'Probing data: load a SHAPE or DMS profile first (Tools ▸ Load probing data…)' })));
  const d = dialog({
    title: 'Conditions', width: 560, body: b => b.append(body),
    actions: [
      { label: 'Standard conditions', action: () => { setConditions(app, defaultConditions()); d.close(); } },
      { label: 'Cancel', close: true },
      { label: 'Apply', primary: true, action: () => {
        const t = Math.min(95, Math.max(0, +temp.value || 37));
        const s = Math.min(3, Math.max(0.005, (+na.value || 1021) / 1000));
        let ligand = null;
        const v = ligSel.value;
        if (v) {
          const base = v === 'custom' ? { name: 'Custom aptamer', motif_seq: mseq.value.trim().toUpperCase().replace(/T/g, 'U'), motif_struct: mst.value.trim() } : { ...LIGANDS[v], preset: v };
          if (!base.motif_seq || base.motif_seq.length !== base.motif_struct.length) { toast('The motif sequence and structure must have the same length.', 'error', 'Conditions'); return; }
          ligand = { ...base, kd_uM: Math.max(1e-9, +kd.value || 1), conc_uM: Math.max(0, +conc.value || 0) };
        }
        setConditions(app, { temperature: t, sodium: s, ligand, useProbing: !!probe.checked && !!ens.probing });
        d.close();
      } },
    ],
  });
}
