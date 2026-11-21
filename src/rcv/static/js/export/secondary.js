// Export ▸ Secondary structure: the sequence and the 2D structure of conformations as text in the formats other
// programs read (mirrors rcv/ss2dexport.py, used by `rcv ss`):
//   dot-bracket (.dbn; the ViennaRNA / FASTA-like form read by RNAfold, RNAplot, forna, VARNA, RNAstructure dot2ct ...)
//   CT, the connect table of mfold and RNAstructure (.ct)
//   BPSEQ (.bpseq)
// The 2D structure of a conformation is the one RCV shows for it: the base pairs found in its 3D coordinates.
import { el, saveBlob, slug } from '../util.js';
import { dialog, select, toast, field } from '../ui.js';

const OPEN = '([{<', CLOSE = ')]}>';

// Partner of every position (0-based, -1 = unpaired) from a dot-bracket with ( ) [ ] { } < >.
export function partners(dbn) {
  const out = new Array(dbn.length).fill(-1);
  const stacks = [[], [], [], []];
  for (let k = 0; k < dbn.length; k++) {
    const o = OPEN.indexOf(dbn[k]), c = CLOSE.indexOf(dbn[k]);
    if (o >= 0) stacks[o].push(k);
    else if (c >= 0 && stacks[c].length) { const i = stacks[c].pop(); out[i] = k; out[k] = i; }
  }
  return out;
}

const energyText = e => (typeof e === 'number' && Number.isFinite(e) ? e.toFixed(2) : null);
const titleOf = r => [r.name, r.id].filter(Boolean).join(' ').replace(/\s+/g, ' ').trim() || 'RNA';

// records: [{ name, id, seq, dbn, energy (kcal/mol or null), also (ids of conformations with the same 2D) }]
export function toDotBracket(records) {
  return records.map(r => {
    const e = energyText(r.energy);
    const extra = [e !== null ? `dG=${e} kcal/mol` : '', r.also?.length ? `also ${r.also.join(',')}` : ''].filter(Boolean).join(' ');
    return `>${titleOf(r)}${extra ? ' ' + extra : ''}\n${r.seq}\n${r.dbn}\n`;
  }).join('');
}

export function toCT(records) {
  return records.map(r => {
    const p = partners(r.dbn), n = r.seq.length, e = energyText(r.energy);
    const head = `${String(n).padStart(5)}  ${e !== null ? `ENERGY = ${e}  ` : ''}${titleOf(r)}`;
    const rows = [...r.seq].map((b, i) => `${String(i + 1).padStart(5)} ${b} ${String(i).padStart(7)} ${String(i + 1 < n ? i + 2 : 0).padStart(4)} ${String(p[i] + 1).padStart(4)} ${String(r.numbers?.[i] ?? i + 1).padStart(4)}`);
    return [head, ...rows].join('\n') + '\n';
  }).join('');
}

export function toBpseq(records) {
  return records.map(r => {
    const p = partners(r.dbn), e = energyText(r.energy);
    return [`# ${titleOf(r)}${e !== null ? ` dG=${e} kcal/mol` : ''}`, ...[...r.seq].map((b, i) => `${i + 1} ${b} ${p[i] + 1}`)].join('\n') + '\n';
  }).join('\n');
}

export const FORMATS = {
  dbn: { label: 'Dot-bracket (.dbn) — ViennaRNA, forna, VARNA', ext: 'dbn', write: toDotBracket, mime: 'text/plain' },
  ct: { label: 'CT, connect table (.ct) — mfold, RNAstructure', ext: 'ct', write: toCT, mime: 'text/plain' },
  bpseq: { label: 'BPSEQ (.bpseq)', ext: 'bpseq', write: toBpseq, mime: 'text/plain' },
};

// Conformations of an ensemble as records. which: 'current' | 'all' | 'distinct' (every different 2D once, with the
// conformations that share it); order: conformation indices as listed.
export function recordsOf(ens, which, order, current) {
  const energy = k => ens.value(k, 'ss_dg') ?? ens.value(k, 'dg2d');
  const numbers = ens.residues.map(r => (/^-?\d+$/.test(String(r.num).trim()) ? parseInt(r.num, 10) : null));   // numbering of the file, when it is plain
  const plain = numbers.length > 0 && numbers.every((v, i) => v !== null && (i === 0 || v > numbers[i - 1]));
  const rec = k => ({ name: ens.name, id: ens.shortId(k), seq: ens.sequence, dbn: ens.conformations[k].dbn, energy: energy(k), numbers: plain ? numbers : null });
  if (which === 'current') return [rec(current)];
  const all = (order?.length ? order : ens.conformations.map((_, k) => k)).map(rec);
  if (which !== 'distinct') return all;
  const seen = new Map();
  for (const r of all) {
    const first = seen.get(r.dbn);
    if (first) first.also.push(r.id); else seen.set(r.dbn, { ...r, also: [] });
  }
  return [...seen.values()];
}

export function exportSecondaryDialog(app) {
  const ens = app.ens;
  if (!ens) return;
  const nDistinct = new Set(ens.conformations.map(c => c.dbn)).size;
  const fmtSel = select(Object.entries(FORMATS).map(([value, f]) => ({ value, label: f.label })), 'dbn', { 'aria-label': 'Format' });
  const whichSel = select([
    { value: 'current', label: `Shown conformation (${ens.shortId(app.state.conf)})` },
    { value: 'distinct', label: `Every different 2D structure once (${nDistinct})` },
    { value: 'all', label: `All conformations (${ens.n})` }], ens.n > 1 ? 'distinct' : 'current', { 'aria-label': 'Conformations' });
  const preview = el('textarea', { class: 'input mono ss-export', rows: 12, readonly: true, spellcheck: 'false', wrap: 'off', 'aria-label': 'Preview', translate: 'no' });
  const text = () => FORMATS[fmtSel.value].write(recordsOf(ens, whichSel.value, app.order, app.state.conf));
  const refresh = () => { preview.value = text(); };
  fmtSel.addEventListener('change', refresh);
  whichSel.addEventListener('change', refresh);
  const mods = (ens.doc.modifications || []).length;
  const d = dialog({ title: 'Export the secondary structure', width: 640,
    body: el('div', {},
      el('p', { class: 'desc', text: 'The sequence and the 2D structure (base pairs found in the 3D coordinates of each conformation), as text for other programs.' }),
      el('div', { class: 'grid2' }, field('Format', fmtSel), field('Conformations', whichSel)),
      preview,
      mods ? el('p', { class: 'desc', text: 'Modified nucleotides are written as their parent base (A, C, G, U): these formats have no place for them.' }) : null),
    actions: [
      { label: 'Close', close: true },
      { label: 'Copy', action: async () => { try { await navigator.clipboard.writeText(text()); toast('Copied to the clipboard.', 'ok'); } catch { preview.select(); document.execCommand?.('copy'); toast('Copied to the clipboard.', 'ok'); } } },
      { label: 'Save…', primary: true, action: async () => {
        const f = FORMATS[fmtSel.value];
        const name = `${slug(ens.name)}_${whichSel.value === 'current' ? slug(ens.shortId(app.state.conf)) : whichSel.value === 'all' ? 'all' : '2d'}.${f.ext}`;
        const saved = await saveBlob(new Blob([text()], { type: f.mime }), name, [{ description: `${f.ext.toUpperCase()} file`, accept: { 'text/plain': [`.${f.ext}`] } }]);
        if (saved) { toast(saved, 'ok', 'Secondary structure saved'); d.close(); }
      } }] });
  refresh();
  return d;
}

