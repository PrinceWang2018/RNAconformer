// Tools menu: the command-line tools of RCV as dialogs (convert, inspect, extract, batch, fetch, HTML page).
import { el, fmt, slug, saveBlob, readFileText } from './util.js';
import { icon } from './icons.js';
import { dialog, field, toast, busy } from './ui.js';
import { ZipWriter } from './export/encoders.js';
import { toPDB } from './ensemble.js';

const STRUCT_RE = /\.(pdb|ent|cif|mmcif)$/i;

// Pick files: native Open dialog (paths) in the desktop window, else the browser chooser (File objects).
function pickFiles(app, { accept = '', multiple = true, kind = 'structures' } = {}) {
  if (app.native?.available()) return app.native.chooseFiles(kind).then(paths => ({ paths: paths || [], files: [] }));
  return new Promise(resolve => {
    const inp = el('input', { type: 'file', accept, multiple: multiple ? true : null, hidden: true });
    inp.addEventListener('change', () => resolve({ paths: [], files: [...inp.files] }));
    document.body.append(inp);
    inp.click();
    setTimeout(() => inp.remove(), 60000);
  });
}

async function rcvRequest(app, route, pick, headers = {}) {
  if (pick.paths.length) return app.api.raw('POST', route, null, { 'X-Path': encodeURIComponent(pick.paths[0]), ...headers });
  return app.api.raw('POST', route, await pick.files[0].arrayBuffer(), { 'Content-Type': 'application/octet-stream', ...headers });
}

// ------------------------------------------------------------------ fetch from the PDB
export function fetchDialog(app) {
  const idIn = el('input', { class: 'input mono', placeholder: 'e.g. 1ANR', maxlength: 4, 'aria-label': 'PDB identifier', style: { textTransform: 'uppercase' } });
  const chainIn = el('input', { class: 'input', placeholder: 'auto (longest RNA chain)', 'aria-label': 'Chain' });
  const examples = el('div', { class: 'chips' }, ...[['1ANR', 'HIV-1 TAR, NMR (20 models)'], ['1EHZ', 'tRNA-Phe, X-ray'], ['2KOC', 'UUCG hairpin, NMR'], ['1Y26', 'adenine riboswitch']]
    .map(([id, d]) => el('button', { class: 'chip-btn', type: 'button', title: d, onclick: () => { idIn.value = id; idIn.focus(); } }, el('b', { text: id }), ` ${d}`)));
  const go = async () => {
    const id = idIn.value.trim().toUpperCase();
    if (!/^[0-9][A-Z0-9]{3}$/.test(id)) { toast('A PDB identifier has 4 characters and starts with a digit, e.g. 1ANR.', 'error', 'Fetch'); return; }
    d.close();
    const done = busy(`Downloading ${id} from the RCSB PDB…`);
    try {
      const s = await app.api.post('/api/fetch', { id, chain: chainIn.value.trim() || null });
      await app.loadServerDocument(s.id);
      (app.ens.commandLog ||= []).push(`fetch ${id}${chainIn.value.trim() ? ', ' + chainIn.value.trim() : ''}`);
      toast(`${id}: ${s.n_conformations} model${s.n_conformations === 1 ? '' : 's'}, ${s.length} nt (saved in the projects folder).`, 'ok', 'Fetched from the PDB');
    } catch (err) { toast(err.message, 'error', `Could not fetch ${id}`); } finally { done(); }
  };
  idIn.addEventListener('keydown', e => { if (e.key === 'Enter') go(); });
  const d = dialog({
    title: 'Fetch a structure from the PDB', width: 520,
    body: b => b.append(el('div', { class: 'grid2' }, field('PDB identifier', idIn), field('RNA chain', chainIn)),
      el('div', { class: 'flabel', text: 'Examples' }), examples,
      el('div', { class: 'desc', text: 'Downloads the mmCIF file from files.rcsb.org (internet needed) into the projects folder and opens its RNA chain; multi-model NMR entries give one conformation per model.' })),
    actions: [{ label: 'Cancel', close: true }, { label: 'Fetch and open', primary: true, action: go }],
  });
  setTimeout(() => idIn.focus(), 0);
}

// ------------------------------------------------------------------ convert to .rcv
export function convertDialog(app) {
  let pick = null;
  const list = el('div', { class: 'filelist tl-files' }, el('div', { class: 'fl-item muted', text: 'No files chosen yet.' }));
  const nameIn = el('input', { class: 'input', placeholder: 'from the files', 'aria-label': 'Ensemble name' });
  const chainIn = el('input', { class: 'input', placeholder: 'auto', 'aria-label': 'Chain' });
  const separate = el('input', { type: 'checkbox' });
  const keep = el('input', { type: 'checkbox', checked: true });
  const show = () => {
    const names = pick ? (pick.paths.length ? pick.paths.map(p => p.split(/[\\/]/).pop()) : pick.files.map(f => f.name)) : [];
    list.replaceChildren(...(names.length ? names.map(n => el('div', { class: 'fl-item' }, icon('file'), el('span', { class: 'nm', text: n }))) : [el('div', { class: 'fl-item muted', text: 'No files chosen yet.' })]));
  };
  const choose = async () => { pick = await pickFiles(app, { accept: '.pdb,.ent,.cif,.mmcif' }); show(); };
  const run = async () => {
    if (!pick || !(pick.paths.length || pick.files.length)) { toast('Choose structure files first.', 'info'); return; }
    d.close();
    const groups = separate.checked ? (pick.paths.length ? pick.paths.map(p => ({ paths: [p], files: [] })) : pick.files.map(f => ({ paths: [], files: [f] }))) : [pick];
    const done = busy(`Converting ${groups.length} ensemble${groups.length === 1 ? '' : 's'}…`);
    const outputs = [];
    try {
      for (const g of groups) {
        let s;
        const opts = { name: groups.length === 1 ? nameIn.value.trim() || undefined : undefined, chain: chainIn.value.trim() || undefined };
        if (g.paths.length) s = await app.api.post('/api/open', { paths: g.paths, ...opts });
        else {
          const payload = [];
          for (const f of g.files) payload.push({ name: f.name, text: await readFileText(f) });
          s = await app.api.post('/api/upload', { files: payload, ...opts });
        }
        const r = await app.api.raw('POST', '/api/rcv/pack', JSON.stringify({ document: s.id }), { 'Content-Type': 'application/json' });
        const stats = JSON.parse(r.headers.get('X-RCV-Stats') || '{}');
        outputs.push({ name: `${slug(s.name)}.rcv`, data: new Uint8Array(await r.arrayBuffer()), stats, summary: s });
        if (keep.checked) await app.loadServerDocument(s.id, { select: groups.length === 1 });
        else await app.api.del(`/api/documents/${s.id}`).catch(() => {});
      }
    } catch (err) { toast(err.message, 'error', 'Conversion failed'); return; } finally { done(); }
    let blob, name;
    if (outputs.length === 1) { blob = new Blob([outputs[0].data], { type: 'application/vnd.rcv+zip' }); name = outputs[0].name; }
    else { const z = new ZipWriter(); for (const o of outputs) z.add(o.name, o.data); blob = z.blob(); name = 'rcv_files.zip'; }
    const saved = await saveBlob(blob, name, outputs.length === 1 ? [{ description: 'RCV ensemble', accept: { 'application/vnd.rcv+zip': ['.rcv'] } }] : undefined);
    const pdb = outputs.reduce((s, o) => s + (o.stats.pdb_bytes || 0), 0);
    if (saved) toast(`${saved} · ${outputs.length} ensemble${outputs.length === 1 ? '' : 's'} · ${(blob.size / 1024).toFixed(0)} kB (${pdb ? (pdb / blob.size).toFixed(0) + '× smaller than PDB' : ''})`, 'ok', 'Converted to .rcv');
  };
  const d = dialog({
    title: 'Convert structures to .rcv', width: 600,
    body: b => b.append(
      el('div', { class: 'desc', text: 'Several PDB/mmCIF files of one RNA (or one multi-model file) become one .rcv file with every conformation, energies, 2D structures and provenance.' }),
      el('div', { class: 'row', style: { gap: '8px', margin: '8px 0' } }, el('button', { class: 'btn', onclick: choose }, icon('file'), 'Choose files…')),
      list,
      el('div', { class: 'grid2', style: { marginTop: '10px' } }, field('Ensemble name', nameIn), field('RNA chain', chainIn)),
      el('label', { class: 'check' }, separate, 'Each file is a different RNA (one .rcv per file, saved together as a ZIP)'),
      el('label', { class: 'check' }, keep, 'Also open the result in the viewer')),
    actions: [{ label: 'Cancel', close: true }, { label: 'Convert and save', primary: true, action: run }],
  });
}

// ------------------------------------------------------------------ inspect an .rcv file
export async function inspectDialog(app) {
  const pick = await pickFiles(app, { accept: '.rcv', multiple: false, kind: 'rcv' });
  if (!pick.paths.length && !pick.files.length) return;
  const done = busy('Reading the .rcv file…');
  let info;
  try { info = await (await rcvRequest(app, '/api/rcv/info', pick)).json(); } catch (err) { toast(err.message, 'error', 'Not a valid .rcv file'); return; } finally { done(); }
  const prov = info.provenance || {};
  const sw = prov.software || {};
  const tb = el('tbody', {}, ...info.conformations.slice(0, 200).map(c => el('tr', {}, el('td', { text: c.id }), el('td', { class: 'mono', text: c.dbn || '' }), el('td', { class: 'num', text: c.value === null || c.value === undefined ? '–' : fmt(c.value) }))));
  const openIt = async () => {
    d.close();
    if (pick.paths.length) await app.openServerPaths([pick.paths[0]]); else await app.openRcvFile(pick.files[0]);
  };
  const d = dialog({
    title: `${info.name} — .rcv file`, width: 720,
    body: b => b.append(el('div', { class: 'tl-info' },
      info.thumbnail ? el('img', { class: 'tl-thumb', src: `data:image/${info.contents.includes('thumbnail.jpg') ? 'jpeg' : 'png'};base64,${info.thumbnail}`, alt: 'Thumbnail', draggable: 'false' }) : null,
      el('dl', { class: 'kv' },
        el('dt', { text: 'Sequence' }), el('dd', { class: 'mono', text: `${info.sequence} (${info.sequence.length} nt)` }),
        el('dt', { text: 'Conformations' }), el('dd', { text: `${info.n_conformations} · ${info.n_atoms} atoms each` }),
        el('dt', { text: 'Origin' }), el('dd', { text: info.origin === 'rnafold3d-engine' ? `rnafold3d prediction${sw.engine_version ? ' (engine ' + sw.engine_version + ')' : ''}` : 'structure files' }),
        el('dt', { text: 'File' }), el('dd', { text: `${(info.bytes / 1024).toFixed(0)} kB — about ${(info.pdb_bytes_estimate / info.bytes).toFixed(0)}× smaller than the same models as PDB (${(info.pdb_bytes_estimate / 1024).toFixed(0)} kB)` }),
        el('dt', { text: 'Saved' }), el('dd', { text: `${(info.saved || '').replace('T', ' ').slice(0, 19)} by RCV ${info.app_version || '?'} · format ${info.format}${info.has_session ? ' · includes a saved view' : ''}` }),
        el('dt', { text: 'Contents' }), el('dd', { class: 'mono', text: info.contents.join('  ') }))),
      el('div', { class: 'section-title', style: { marginTop: '10px' } }, `Conformations · ${info.metric.label}${info.metric.unit ? ' (' + info.metric.unit + ')' : ''}`),
      el('div', { class: 'wl-table', style: { maxHeight: '220px' } }, el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'ID' }), el('th', { text: '2D structure' }), el('th', { class: 'num', text: info.metric.label }))), tb)),
      info.methods ? el('details', { class: 'tl-methods' }, el('summary', { text: 'Methods text stored in the file' }), el('p', { class: 'desc', text: info.methods })) : null),
    actions: [{ label: 'Close', close: true }, { label: 'Open in the viewer', primary: true, action: openIt }],
  });
}

// ------------------------------------------------------------------ extract to PDB
export function extractDialog(app) {
  let mode = 'zip';
  const seg = el('div', { class: 'seg' }, el('button', { 'data-m': 'zip', text: 'One PDB per conformation (ZIP)' }), el('button', { 'data-m': 'multi', text: 'One multi-model PDB' }));
  const sync = () => { for (const b of seg.children) b.classList.toggle('on', b.dataset.m === mode); };
  seg.addEventListener('click', e => { const m = e.target.closest('button')?.dataset.m; if (m) { mode = m; sync(); } });
  sync();
  const fromFile = async () => {
    const pick = await pickFiles(app, { accept: '.rcv', multiple: false, kind: 'rcv' });
    if (!pick.paths.length && !pick.files.length) return;
    d.close();
    const done = busy('Extracting…');
    try {
      const r = await rcvRequest(app, '/api/rcv/extract', pick, { 'X-Mode': mode });
      const base = (pick.paths[0] || pick.files[0].name).split(/[\\/]/).pop().replace(/\.rcv$/i, '');
      const blob = new Blob([await r.arrayBuffer()], { type: mode === 'zip' ? 'application/zip' : 'chemical/x-pdb' });
      const saved = await saveBlob(blob, `${slug(base)}${mode === 'zip' ? '_pdb.zip' : '_models.pdb'}`);
      if (saved) toast(`${saved} · ${(blob.size / 1024).toFixed(0)} kB`, 'ok', 'Extracted');
    } catch (err) { toast(err.message, 'error', 'Extraction failed'); } finally { done(); }
  };
  const fromOpen = async () => {
    const ens = app.ens;
    if (!ens) return;
    d.close();
    if (mode === 'multi') { await app.exportPDB(true); return; }
    const z = new ZipWriter();
    const enc = new TextEncoder();
    for (const k of app.order) z.add(`${slug(ens.shortId(k))}.pdb`, enc.encode(toPDB(ens, ens.coords[k], { remarks: [`RCV EXPORT ${ens.name}`, `CONFORMATION ${ens.shortId(k)}`] }) + 'END\n'));
    const blob = z.blob();
    const saved = await saveBlob(blob, `${slug(ens.name)}_pdb.zip`);
    if (saved) toast(`${saved} · ${ens.n} PDB files`, 'ok', 'Extracted');
  };
  const d = dialog({
    title: 'Extract conformations as PDB', width: 560,
    body: b => b.append(el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Output' }), seg),
      el('div', { class: 'desc', text: 'From an .rcv file on this computer, or from the ensemble open in the viewer (coordinates as displayed).' })),
    actions: [{ label: 'Cancel', close: true }, { label: 'From the open ensemble', action: fromOpen }, { label: 'Choose an .rcv file…', primary: true, action: fromFile }],
  });
  if (!app.ens) d.buttons['From the open ensemble'].disabled = true;
}

// ------------------------------------------------------------------ batch from FASTA
export function batchFromFasta(app) {
  app.showPage('welcome');
  app.welcome.adv.open = true;
  app.welcome.chooseFasta();
}

export function toolsMenuItems(app) {
  const has = !!app.ens, online = app.api.online;
  return [
    { label: 'Convert structures to .rcv…', action: () => convertDialog(app), disabled: !online },
    { label: 'Inspect an .rcv file…', action: () => inspectDialog(app), disabled: !online },
    { label: 'Extract conformations as PDB…', action: () => extractDialog(app), disabled: !online && !has },
    { sep: true },
    { label: 'Batch prediction from a FASTA file…', action: () => batchFromFasta(app), disabled: !online },
    { label: 'Fetch a structure from the PDB…', kbd: 'Ctrl+Shift+F', action: () => fetchDialog(app), disabled: !online },
    { sep: true },
    { label: 'Interactive HTML page of this ensemble…', action: () => app.exportHTML(), disabled: !has || !online },
    { label: 'Figure panels…', action: () => app.openPanels(), disabled: !has },
    { label: 'Plot studio (RMSD, RMSF, energy)…', action: () => app.openPlotStudio?.(), disabled: !has },
    { sep: true },
    { label: 'Command console', kbd: '`', action: () => app.action('console') },
  ];
}
