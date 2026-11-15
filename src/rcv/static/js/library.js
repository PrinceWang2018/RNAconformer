// Project library (File > Project library…, Ctrl+L): every prediction run, .rcv file and structure opened in RCV,
// plus the runs and .rcv files in the watched folders, with a thumbnail and the parameters — searchable by name,
// sequence, folder or parameter, and opened with one click. The local server keeps the index (rcv/library.py);
// the viewer draws a thumbnail the first time an item is shown.
import { el } from './util.js';
import { icon } from './icons.js';
import { dialog, toast, askText } from './ui.js';
import { composeFigure } from './export/figure.js';

const KINDS = [['all', 'All'], ['run', 'Prediction runs'], ['rcv', '.rcv files'], ['structure', 'Structures']];
const SORTS = [['when', 'Last used'], ['name', 'Name'], ['length', 'Length'], ['opens', 'Most opened']];
const KIND_ICON = { run: 'flask', rcv: 'cube', structure: 'file', pdb: 'globe' };

// Draws a thumbnail for an ensemble the first time it is shown in this session (if it came from a file).
export function watchThumbnails(app) {
  const done = new Set();
  let timer = null;
  app.on('ensemble', () => {
    clearTimeout(timer);
    timer = setTimeout(async () => {
      const ens = app.ens;
      if (!ens?.serverId || done.has(ens.serverId) || !app.api.online || app.page !== 'viewer') return;
      done.add(ens.serverId);
      try {
        const c = await composeFigure(app, { width: 320, height: 200, transparent: false, background: null, legend: false, labels: false });
        await app.api.post('/api/library/thumb', { document: ens.serverId, image: c.toDataURL('image/jpeg', 0.78) });
      } catch { /* a thumbnail is optional */ }
    }, 1800);
  });
}

const paramText = p => {
  if (!p) return '';
  const engine = String(p.engine || '').startsWith('rnafold3d') ? 'RCV engine' : /alphafold/i.test(p.engine || '') ? 'AlphaFold 3' : p.engine;
  return [engine, p.models ? `${p.models} models` : '', p.trajectories ? `${p.trajectories} trajectories` : '',
    p.seed !== undefined ? `seed ${p.seed}` : '', p.templates ? `${p.templates} template${p.templates === 1 ? '' : 's'}` : (p.templates === 0 ? 'no template' : ''),
    p['2d'] ? `2D ${p['2d']}` : '', p.temperature ? `${p.temperature} °C` : ''].filter(Boolean).join(' · ');
};
const dateText = ts => ts ? new Date(ts * 1000).toISOString().replace('T', ' ').slice(0, 16) : '';
const kindOf = v => (v.kind === 'pdb' ? 'structure' : v.kind);

export async function projectLibrary(app) {
  if (!app.api.online) { toast('The project library is kept by the local RCV engine, which is not reachable.', 'error', 'Project library'); return; }
  let data;
  try { data = await app.api.get('/api/library/index'); } catch (err) { toast(err.message, 'error', 'Project library'); return; }
  let kind = 'all', sort = 'when';
  const search = el('input', { class: 'input', type: 'search', placeholder: 'Search by name, sequence, folder or parameter…', 'aria-label': 'Search the library', spellcheck: 'false' });
  const kindSeg = el('div', { class: 'seg', role: 'radiogroup', 'aria-label': 'Kind' }, ...KINDS.map(([k, t]) => el('button', { type: 'button', 'data-k': k, text: t })));
  const sortSel = el('select', { class: 'select', 'aria-label': 'Sort' }, ...SORTS.map(([k, t]) => el('option', { value: k, text: t })));
  const grid = el('div', { class: 'lib-grid', role: 'list' });
  const count = el('span', { class: 'muted lib-count' });
  const foldersBox = el('div', { class: 'lib-folders' });

  const open = async v => {
    if (!v.exists) { toast(`${v.path} was moved or deleted.`, 'error', v.name); return; }
    d.close();
    try { await app.openServerPaths([v.path]); app.showPage('viewer'); } catch (err) { toast(err.message, 'error', v.name); }
  };
  const forget = async v => {
    await app.api.post('/api/library/forget', { path: v.path });
    data.items = data.items.filter(x => x.key !== v.key);
    render();
  };
  const card = v => {
    const thumb = v.thumb
      ? el('img', { class: 'lib-thumb', alt: '', loading: 'lazy', src: `/api/library/thumb/${v.key}?token=${encodeURIComponent(app.api.token || '')}&v=${Math.round(v.when)}` })
      : el('div', { class: 'lib-thumb lib-noimg' }, icon(KIND_ICON[v.kind] || 'cube'));
    const meta = [v.length ? `${v.length} nt` : '', v.n_conformations ? `${v.n_conformations} conformation${v.n_conformations === 1 ? '' : 's'}` : ''].filter(Boolean).join(' · ');
    return el('div', { class: `lib-card${v.exists ? '' : ' missing'}`, role: 'listitem' },
      el('button', { class: 'lib-open', type: 'button', title: v.path, onclick: () => open(v) }, thumb,
        el('div', { class: 'lib-body' },
          el('b', { text: v.name || v.path.split(/[\\/]/).pop(), translate: 'no' }),
          el('span', { class: 'lib-meta' }, el('span', { class: 'lib-kind', text: KINDS.find(k => k[0] === kindOf(v))?.[1] || v.kind }), meta ? ' · ' : '', meta ? el('span', { text: meta }) : ''),
          paramText(v.params) ? el('span', { class: 'lib-params', text: paramText(v.params), translate: 'no' }) : null,
          el('span', { class: 'lib-when' }, v.exists ? dateText(v.when) : el('span', { class: 'state-bad', text: 'File not found' }), v.opens ? ' · ' : '', v.opens ? el('span', { text: `opened ${v.opens}×` }) : ''))),
      el('button', { class: 'ss-btn lib-forget', type: 'button', title: 'Remove from the library (the file stays where it is)', 'aria-label': 'Remove from the library', onclick: () => forget(v) }, icon('x')));
  };
  const render = () => {
    for (const b of kindSeg.children) b.classList.toggle('on', b.dataset.k === kind);
    const words = search.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    let list = data.items.filter(v => kind === 'all' || kindOf(v) === kind);
    if (words.length) list = list.filter(v => {
      const hay = [v.name, v.sequence, v.path, v.kind, v.origin, paramText(v.params), v.length ? `${v.length} nt` : ''].join(' ').toLowerCase();
      return words.every(w => hay.includes(w));
    });
    const by = { name: (a, b) => (a.name || '').localeCompare(b.name || ''), length: (a, b) => (b.length || 0) - (a.length || 0),
      opens: (a, b) => (b.opens || 0) - (a.opens || 0), when: (a, b) => b.when - a.when }[sort];
    list = list.slice().sort(by);
    grid.replaceChildren(...list.slice(0, 300).map(card));
    if (!list.length) grid.append(el('div', { class: 'lib-empty', text: data.items.length ? 'Nothing matches the search.' : 'The library fills as you open prediction runs, .rcv files and structures; runs in the watched folders appear here too.' }));
    count.textContent = `${list.length} of ${data.items.length}`;
  };
  const drawFolders = () => {
    const row = (f, removable) => el('div', { class: 'lib-folder' }, icon('folder'), el('span', { class: 'mono', text: f, title: f, translate: 'no' }), el('span', { class: 'grow' }),
      removable ? el('button', { class: 'ss-btn', type: 'button', title: 'Stop watching this folder', 'aria-label': 'Stop watching this folder', onclick: async () => {
        data.folders = (await app.api.post('/api/library/folders', { path: f, add: false })).folders; await reload();
      } }, icon('x')) : el('span', { class: 'muted', text: 'projects folder' }));
    foldersBox.replaceChildren(
      ...data.default_folders.map(f => row(f, false)), ...data.folders.map(f => row(f, true)),
      el('button', { class: 'btn small', type: 'button', onclick: async () => {
        const f = app.native?.available() ? await app.native.chooseFolder() : await askText('Watch a folder', 'Folder on this computer', '');
        if (!f) return;
        try { data.folders = (await app.api.post('/api/library/folders', { path: f, add: true })).folders; await reload(); }
        catch (err) { toast(err.message, 'error', 'Watched folders'); }
      } }, icon('plus'), 'Watch a folder…'));
  };
  const reload = async () => { data = await app.api.get('/api/library/index'); drawFolders(); render(); };
  search.addEventListener('input', render);
  kindSeg.addEventListener('click', e => { const k = e.target.closest('button')?.dataset.k; if (k) { kind = k; render(); } });
  sortSel.addEventListener('change', () => { sort = sortSel.value; render(); });
  drawFolders();
  render();
  const d = dialog({
    title: 'Project library', width: 1060,
    body: b => b.append(
      el('div', { class: 'lib-bar' }, search, kindSeg, sortSel, count),
      grid,
      el('details', { class: 'adv lib-watch' }, el('summary', { text: 'Watched folders: prediction runs and .rcv files in them are listed even if not opened yet' }), foldersBox)),
    actions: [{ label: 'Close', close: true }],
  });
  setTimeout(() => search.focus(), 50);
  return d;
}
