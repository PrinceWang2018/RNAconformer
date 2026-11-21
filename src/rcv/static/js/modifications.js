// Modified nucleotides (0.7.6): a sub-window of the start page to put modifications on a sequence. Click nucleotides
// (Shift: range, Ctrl/⌘: add), pick a modification from the list for that base, see what it does to pairing; the
// result goes back into the sequence as [code] tokens (modseq.js), which the engine reads (rnafold3d/modifications.py).
import { el } from './util.js';
import { icon } from './icons.js';
import { dialog, toast } from './ui.js';
import * as modseq from './modseq.js';

const PAIRING = {
  measured: ['measured', 'Pairing from measured nearest-neighbour parameters (ViennaRNA)'],
  estimated: ['estimated', 'Pairing from Turner parameters of a reference pair with an estimated correction'],
  blocked: ['blocks pairing', 'Watson-Crick edge blocked: kept unpaired in the 2D structures'],
  unchanged: ['no change', 'No known effect on base pairing'],
};

function splitFasta(text) {
  const header = [], seqLines = [], structure = [];
  for (const raw of String(text || '').split(/\r?\n/)) {
    const l = raw.trim();
    if (!l) continue;
    if (l.startsWith('>')) { if (!header.length) header.push(l); continue; }
    if (/^[().[\]{}<>]+$/.test(l) && seqLines.length) structure.push(l); else seqLines.push(l);
  }
  return { header: header[0] || '', seq: modseq.cleanLine(seqLines.join('')), structure: structure.join('') };
}

// The pair whose stacking energies a correction is added to, when it is not the parent base with this partner
// (8-oxoG with A: the energies of a U-A pair).
function refPair(e, partner) {
  const ref = e.ddg_ref?.[partner];
  return ref && ref !== e.parent + partner ? ` (on ${ref[0]}–${ref[1]})` : '';
}

export async function modificationsDialog(app, { text, onApply }) {
  const cat = await modseq.loadCatalogue(app.api);
  if (!cat) { toast('The modification list comes from the engine, which is not reachable.', 'error', 'Modified nucleotides'); return; }
  const parts = splitFasta(text);
  const parsed = modseq.parse(parts.seq);
  if (parsed.problems.length || parsed.bad.length) {
    toast(parsed.problems[0] || `Unsupported letters: ${parsed.bad.join(' ')}`, 'error', 'Modified nucleotides');
    return;
  }
  const n = parsed.seq.length;
  // state: the unmodified letter of every position, its modification, phosphorothioate linkages
  const mod = new Map(), ps = new Set();
  for (const m of parsed.mods) (m.code === 'ps' ? ps.add(m.pos - 1) : mod.set(m.pos - 1, m.code));
  const orig = [...parsed.seq].map((b, i) => (mod.has(i) ? modseq.entry(mod.get(i))?.parent || b : b));
  let sel = new Set(), anchor = null, chosen = null;

  // ---- sequence map
  const map = el('div', { class: 'mod-map', translate: 'no', role: 'listbox', 'aria-multiselectable': 'true', 'aria-label': 'Sequence' });
  const cells = orig.map((b, i) => {
    const c = el('button', { type: 'button', class: 'mod-nt', 'data-i': i });
    c.addEventListener('click', e => {
      if (e.shiftKey && anchor !== null) { const [a, z] = [Math.min(anchor, i), Math.max(anchor, i)]; if (!(e.ctrlKey || e.metaKey)) sel.clear(); for (let k = a; k <= z; k++) sel.add(k); }
      else if (e.ctrlKey || e.metaKey) { sel.has(i) ? sel.delete(i) : sel.add(i); anchor = i; }
      else { sel = new Set([i]); anchor = i; }
      // the details follow the selection: its modification, or nothing when the chosen one is for another base
      const one = sel.size === 1 ? [...sel][0] : null;
      if (one !== null && mod.has(one)) chosen = mod.get(one);
      else if (chosen && ![...sel].every(k => orig[k] === modseq.entry(chosen)?.parent)) chosen = null;
      refresh();
    });
    return c;
  });
  for (let i = 0; i < n; i += 10) {
    map.append(el('div', { class: 'mod-row' }, el('span', { class: 'mod-pos', text: String(i + 1) }), ...cells.slice(i, i + 10)));
  }
  const selInfo = el('div', { class: 'mod-selinfo' });
  const quick = el('div', { class: 'mod-quick' }, el('span', { class: 'muted', text: 'Select:' }),
    ...['A', 'C', 'G', 'U'].map(b => el('button', { class: 'linkbtn', type: 'button', translate: 'no', text: `all ${b}`,
      onclick: () => { sel = new Set(orig.map((x, i) => (x === b ? i : -1)).filter(i => i >= 0)); refresh(); } })),
    el('button', { class: 'linkbtn', type: 'button', text: 'none', onclick: () => { sel.clear(); refresh(); } }));

  // ---- the list of modifications for the selected base(s)
  const search = el('input', { class: 'input', type: 'search', placeholder: 'Search: m6A, pseudouridine, inosine, 8-oxo, 2\'-O-methyl, LNA…', 'aria-label': 'Search modifications' });
  const list = el('div', { class: 'mod-list', role: 'listbox', 'aria-label': 'Modifications' });
  const detail = el('div', { class: 'mod-detail' });
  const applyBtn = el('button', { class: 'btn primary small', type: 'button', disabled: true, onclick: () => apply(chosen) });
  const removeBtn = el('button', { class: 'btn small', type: 'button', onclick: () => { for (const i of sel) { mod.delete(i); } refresh(); } }, 'Remove modification');
  const psBox = el('input', { type: 'checkbox', 'aria-label': 'Phosphorothioate' });
  psBox.addEventListener('change', () => { for (const i of sel) if (i > 0) (psBox.checked ? ps.add(i) : ps.delete(i)); refresh(); });
  const psLabel = el('label', { class: 'check', title: 'One non-bridging oxygen of the phosphate before this nucleotide replaced by sulfur (written as * in the sequence: A*G). Nuclease-resistant antisense and siRNA backbones; no effect on the 2D.' },
    psBox, el('span', { text: 'Phosphorothioate linkage before the selected nucleotides' }));
  search.addEventListener('input', () => fillList());

  function selectedBases() { return new Set([...sel].map(i => orig[i])); }

  function fillList() {
    const bases = selectedBases();
    const q = search.value.trim().toLowerCase();
    const rows = cat.entries.filter(e => e.kind !== 'backbone' && (!bases.size || (bases.size === 1 && bases.has(e.parent)))
      && (!q || [e.code, e.name, e.ccd, ...(e.aliases || []), e.note].some(t => String(t || '').toLowerCase().includes(q))));
    list.replaceChildren();
    if (bases.size > 1) { list.append(el('p', { class: 'help', text: 'The selection holds different bases: select nucleotides of one base (A, C, G or U) to modify them together.' })); return; }
    if (!rows.length) { list.append(el('p', { class: 'help', text: 'No modification matches.' })); return; }
    let lastCat = null;
    for (const e of rows) {
      if (e.category !== lastCat) { lastCat = e.category; list.append(el('div', { class: 'mod-cat', text: cat.categories[e.category] || e.category })); }
      const [badge] = PAIRING[e.pairing] || ['', ''];
      const row = el('button', { type: 'button', class: `mod-item${chosen === e.code ? ' on' : ''}`, role: 'option', 'aria-selected': chosen === e.code ? 'true' : 'false',
        onclick: () => { chosen = e.code; fillList(); showDetail(); },
        ondblclick: () => { chosen = e.code; apply(e.code); } },
        el('b', { class: 'mod-code', translate: 'no', text: modseq.display(e.code) }),
        el('span', { class: 'mod-name', translate: 'no', text: e.name }),
        el('span', { class: 'grow' }),
        el('span', { class: `mod-badge p-${e.pairing}`, text: badge }));
      list.append(row);
    }
  }

  function showDetail() {
    const e = modseq.entry(chosen);
    detail.replaceChildren();
    if (!e) { detail.append(el('p', { class: 'help', text: 'Choose a modification to see what it does.' })); applyBtn.disabled = true; applyBtn.textContent = 'Apply'; return; }
    const [badge, what] = PAIRING[e.pairing] || ['', ''];
    const partners = e.pairing === 'blocked' ? 'none (kept unpaired)' : e.partners ? e.partners.join(', ') : 'as the unmodified base';
    detail.append(
      el('div', { class: 'mod-dhead' }, el('b', { translate: 'no', text: `${modseq.display(e.code)} · ${e.name}` }), el('span', { class: `mod-badge p-${e.pairing}`, title: what, text: badge })),
      el('div', { class: 'mod-facts' },
        el('span', {}, 'Pairs with: ', el('b', { translate: 'no', text: partners })),
        e.ddg ? el('span', { translate: 'no', text: `ΔΔG per stack: ${Object.entries(e.ddg).map(([p, v]) => `${p} ${v > 0 ? '+' : ''}${v.toFixed(2)}${refPair(e, p)}`).join(', ')} kcal/mol` }) : null,
        e.ccd ? el('span', {}, 'PDB component: ', el('b', { translate: 'no', text: e.ccd })) : null,
        e.build !== e.parent ? el('span', { text: `${e.parent} → built on ${e.build}` }) : null,
        e.syn ? el('span', { text: 'syn when unpaired' }) : null,
        e.pucker ? el('span', { text: `Sugar: ${e.pucker}` }) : null,
        e.core_only ? el('span', { text: '3D: base core only' }) : null),
      el('p', { class: 'mod-note', translate: 'no', text: e.note }));
    const k = [...sel].filter(i => orig[i] === e.parent).length;
    applyBtn.disabled = !k;
    applyBtn.textContent = k ? `Apply to ${k} selected` : 'Select nucleotides first';
  }

  function apply(code) {
    const e = modseq.entry(code);
    if (!e) return;
    let k = 0;
    for (const i of sel) if (orig[i] === e.parent) { mod.set(i, code); k++; }
    if (!k) { toast(`${modseq.display(code)} modifies ${e.parent}: select ${e.parent} nucleotides first.`, 'info', 'Modified nucleotides'); return; }
    refresh();
  }

  // ---- presets
  const presets = el('div', { class: 'mod-presets' }, el('span', { class: 'muted', text: 'Presets:' }),
    el('button', { class: 'linkbtn', type: 'button', title: 'N6-methyladenosine at the central A of every DRACH motif (D = A/G/U, R = A/G, H = A/C/U), the consensus of the METTL3-METTL14 writer. Not every DRACH site is methylated in a cell.',
      onclick: () => { const s = modseq.drachSites(orig.join('')).filter(i => !mod.has(i)); for (const i of s) mod.set(i, 'm6A'); refresh(); toast(s.length ? `m6A placed at ${s.length} DRACH site(s).` : 'No unmodified DRACH motif in this sequence.', 'info', 'Modified nucleotides'); } }, 'm6A at DRACH motifs'),
    el('button', { class: 'linkbtn', type: 'button', title: 'Every U becomes N1-methylpseudouridine, as in mRNA vaccines.',
      onclick: () => { orig.forEach((b, i) => { if (b === 'U') mod.set(i, 'm1Y'); }); refresh(); } }, 'every U → m1Ψ'),
    el('button', { class: 'linkbtn', type: 'button', onclick: () => { orig.forEach((b, i) => { if (b === 'U') mod.set(i, 'Y'); }); refresh(); } }, 'every U → Ψ'),
    el('button', { class: 'linkbtn', type: 'button', onclick: () => { mod.clear(); ps.clear(); refresh(); } }, 'clear all'));

  // ---- summary of the modifications and the resulting notation
  const table = el('div', { class: 'mod-summary' });
  const preview = el('div', { class: 'mod-preview mono', translate: 'no' });

  function currentMods() {
    const out = [];
    for (let i = 0; i < n; i++) {
      if (mod.has(i)) out.push({ pos: i + 1, code: mod.get(i) });
      if (ps.has(i)) out.push({ pos: i + 1, code: 'ps' });
    }
    return out;
  }
  function currentSeq() { return orig.map((b, i) => (mod.has(i) ? modseq.entry(mod.get(i))?.build || b : b)).join(''); }

  function refresh() {
    cells.forEach((c, i) => {
      const code = mod.get(i);
      const e = code ? modseq.entry(code) : null;
      c.className = `mod-nt${sel.has(i) ? ' sel' : ''}${e ? ` mod c-${e.category}` : ''}${ps.has(i) ? ' ps' : ''}`;
      c.replaceChildren(el('span', { class: 'b', text: orig[i] }), e ? el('span', { class: 'm', text: modseq.display(code) }) : '');
      c.title = `${i + 1} ${orig[i]}${e ? ` → ${e.name}` : ''}${ps.has(i) ? ' · phosphorothioate before it' : ''}`;
      c.setAttribute('aria-selected', sel.has(i) ? 'true' : 'false');
    });
    const picked = [...sel].sort((a, b) => a - b);
    selInfo.textContent = !picked.length ? 'Click a nucleotide (Shift: range, Ctrl/⌘: add more).'
      : picked.length > 1 ? `Selected: ${picked.length} nucleotides`
        : mod.has(picked[0]) ? `Selected: ${orig[picked[0]]}${picked[0] + 1} · now ${modseq.display(mod.get(picked[0]))}` : `Selected: ${orig[picked[0]]}${picked[0] + 1}`;
    removeBtn.disabled = !picked.some(i => mod.has(i));
    psBox.disabled = !picked.some(i => i > 0);
    psBox.checked = picked.length > 0 && picked.filter(i => i > 0).every(i => ps.has(i));
    if (picked.length === 1 && mod.has(picked[0]) && !chosen) chosen = mod.get(picked[0]);
    fillList();
    showDetail();
    const ms = currentMods();
    table.replaceChildren();
    if (!ms.length) table.append(el('p', { class: 'help', text: 'No modified nucleotides yet.' }));
    else {
      const tb = el('tbody');
      for (const m of ms) {
        const e = modseq.entry(m.code);
        tb.append(el('tr', {},
          el('td', { class: 'num', text: String(m.pos) }),
          el('td', { translate: 'no' }, el('b', { text: modseq.display(m.code) })),
          el('td', { translate: 'no', text: e?.name || m.code }),
          el('td', {}, el('span', { class: `mod-badge p-${e?.pairing || 'unchanged'}`, text: (PAIRING[e?.pairing] || PAIRING.unchanged)[0] })),
          el('td', {}, el('button', { class: 'btn ghost icon-only small', type: 'button', 'aria-label': 'Remove', title: 'Remove',
            onclick: () => { m.code === 'ps' ? ps.delete(m.pos - 1) : mod.delete(m.pos - 1); refresh(); } }, icon('x')))));
      }
      table.append(el('div', { class: 'mod-count', text: ms.length === 1 ? `1 modification · ${modseq.summary(ms)}` : `${ms.length} modifications · ${modseq.summary(ms)}` }),
        el('div', { class: 'mod-table' }, el('table', { class: 'table' },
          el('thead', {}, el('tr', {}, el('th', { class: 'num', text: '#' }), el('th', { text: 'Code' }), el('th', { text: 'Modification' }), el('th', { text: 'Pairing' }), el('th', { text: '' }))), tb)));
    }
    preview.textContent = modseq.format(currentSeq(), ms);
  }

  const d = dialog({
    title: 'Modified nucleotides', width: 940,
    body: b => b.append(
      el('p', { class: 'help', text: 'Click nucleotides, then choose a modification for them. Measured parameters (ViennaRNA) or estimated corrections change base pairing in the 2D structures; Watson-Crick-blocking modifications stay unpaired; the modification\'s atoms are built into every 3D model.' }),
      el('div', { class: 'mod-grid' },
        el('div', { class: 'mod-left' }, map, selInfo, quick, presets),
        el('div', { class: 'mod-right' }, search, list, detail,
          el('div', { class: 'row', style: { gap: '8px', flexWrap: 'wrap' } }, applyBtn, removeBtn), psLabel)),
      el('div', { class: 'section-title', text: 'In this sequence' }), table,
      el('div', { class: 'section-title', text: 'Written in the sequence as' }), preview),
    actions: [{ label: 'Cancel', close: true },
      { label: 'Use in the sequence', primary: true, action: () => {
        const seqText = modseq.format(currentSeq(), currentMods());
        onApply([parts.header, seqText, parts.structure].filter(Boolean).join('\n'));
        d.close();
      } }],
  });
  refresh();
  setTimeout(() => search.focus(), 0);
}
