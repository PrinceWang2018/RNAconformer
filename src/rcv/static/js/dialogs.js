// Application dialogs: open, prediction, image/video export, preferences, plugins, shortcuts, about.
import { el, fmt, slug, saveBlob, readFileText, debounce } from './util.js';
import { icon } from './icons.js';
import { dialog, field, select, toast } from './ui.js';
import { exportFigure, composeFigure } from './export/figure.js';
import { recordVideo, videoFormats, framePlan, transitionPath } from './export/video.js';
import { SCHEMES } from './colors.js';
import { composePanels, exportPanels, defaultPanelConfs } from './export/panels.js';
import { LANGUAGES, currentLanguage } from './i18n.js';

const STRUCT_RE = /\.(pdb|ent|cif|mmcif)$/i;

// ------------------------------------------------------------------ open
export function openDialog(app, initialTab = 'local') {
  let tab = initialTab;
  const nameIn = el('input', { class: 'input', placeholder: 'optional', 'aria-label': 'Ensemble name' });
  const chainIn = el('input', { class: 'input', placeholder: 'auto (longest RNA chain)', 'aria-label': 'Chain', style: { width: '100%' } });
  const seg = el('div', { class: 'seg', style: { marginBottom: '12px' } },
    el('button', { 'data-t': 'local', text: 'This computer' }),
    el('button', { 'data-t': 'server', text: 'Folders' }),
    el('button', { 'data-t': 'library', text: 'Prediction runs' }));
  const pane = el('div');
  const opts = el('div', { class: 'grid2', style: { marginTop: '12px' } }, field('Ensemble name', nameIn), field('RNA chain', chainIn));
  const d = dialog({
    title: 'Open structures', width: 640,
    body: b => b.append(seg, pane, opts),
    actions: [{ node: el('span', { class: 'muted', style: { fontSize: '12px' }, id: 'open-hint' }) }, { spacer: true }, { label: 'Cancel', close: true },
      { id: 'open', label: 'Open', primary: true, action: () => openServerSelection() }],
  });
  const openBtn = d.buttons.open;
  const opt = () => ({ name: nameIn.value.trim() || undefined, chain: chainIn.value.trim() || undefined });

  // --- local files
  const fileIn = el('input', { type: 'file', multiple: true, accept: '.rcv,.pdb,.ent,.cif,.mmcif,.json', hidden: true });
  const dirIn = el('input', { type: 'file', hidden: true });
  dirIn.setAttribute('webkitdirectory', '');
  const handleFiles = async files => {
    d.close();
    await app.openLocalFiles([...files], opt());
  };
  fileIn.addEventListener('change', () => fileIn.files.length && handleFiles(fileIn.files));
  dirIn.addEventListener('change', () => dirIn.files.length && handleFiles([...dirIn.files].filter(f => STRUCT_RE.test(f.name) || f.name === 'manifest.json')));

  // --- server folders
  let cwd = null, selected = new Set();
  const pathIn = el('input', { class: 'input', 'aria-label': 'Folder path' });
  const list = el('div', { class: 'filelist' });
  const places = el('div', { class: 'places' });
  const browseTo = async p => {
    try {
      const r = await app.api.get(`/api/browse?path=${encodeURIComponent(p || '')}`);
      cwd = r; selected = new Set();
      pathIn.value = r.path;
      places.replaceChildren(...r.places.map(pl => el('button', { class: 'btn small', text: pl.name, title: pl.path, onclick: () => browseTo(pl.path) })));
      list.replaceChildren();
      if (r.parent) list.append(el('div', { class: 'fl-item', ondblclick: () => browseTo(r.parent), onclick: () => browseTo(r.parent) }, icon('folder'), el('span', { class: 'nm', text: '..' })));
      for (const e of r.entries) {
        const item = el('div', { class: 'fl-item', title: e.path },
          icon(e.type === 'dir' ? 'folder' : 'file'), el('span', { class: 'nm', text: e.name }),
          e.package ? el('span', { class: 'badge engine', text: 'ensemble' }) : null,
          e.rcv ? el('span', { class: 'badge engine', text: 'RCV' }) : null,
          e.size ? el('span', { class: 'sz', text: `${Math.max(1, Math.round(e.size / 1024))} KB` }) : null);
        if (e.type === 'dir') item.addEventListener('click', () => browseTo(e.path));
        else if (e.rcv) item.addEventListener('click', async () => { d.close(); await app.openServerPaths([e.path]); });
        else item.addEventListener('click', () => { item.classList.toggle('sel'); if (selected.has(e.path)) selected.delete(e.path); else selected.add(e.path); hint(); });
        list.append(item);
      }
      hint();
    } catch (err) { toast(err.message, 'error', 'Cannot open folder'); }
  };
  const hint = () => {
    const h = document.getElementById('open-hint');
    if (!h || tab !== 'server' || !cwd) return;
    const nFiles = cwd.entries.filter(e => e.type === 'file').length;
    h.textContent = cwd.is_package ? 'This folder is an ensemble package' : selected.size ? `${selected.size} file(s) selected` : nFiles ? `Open all ${nFiles} structure files in this folder` : 'No structure files here';
  };
  const openServerSelection = async () => {
    if (tab !== 'server' || !cwd) return;
    const paths = selected.size ? [...selected] : [cwd.path];
    d.close();
    await app.openServerPaths(paths, opt());
  };

  const render = () => {
    for (const b of seg.children) b.classList.toggle('on', b.dataset.t === tab);
    openBtn.hidden = tab !== 'server';
    const h = document.getElementById('open-hint'); if (h) h.textContent = '';
    pane.replaceChildren();
    if (tab === 'local') {
      const drop = el('div', { class: 'preview-box', style: { height: '170px', backgroundImage: 'none', flexDirection: 'column', gap: '10px', borderStyle: 'dashed' } },
        el('div', { style: { fontWeight: 600 }, text: 'Drop an .rcv file, PDB / mmCIF files or an ensemble folder here' }),
        el('div', { class: 'row', style: { flex: 'none' } },
          el('button', { class: 'btn primary', onclick: () => fileIn.click() }, icon('file'), 'Choose files…'),
          el('button', { class: 'btn', onclick: () => dirIn.click() }, icon('folder'), 'Choose folder…')));
      drop.addEventListener('dragover', e => { e.preventDefault(); drop.style.borderColor = 'var(--accent)'; });
      drop.addEventListener('dragleave', () => { drop.style.borderColor = ''; });
      drop.addEventListener('drop', async e => { e.preventDefault(); e.stopPropagation(); d.close(); await app.openDropped(e.dataTransfer, opt()); });
      pane.append(drop, fileIn, dirIn,
        el('div', { class: 'desc', text: 'All structures must be the same RNA. Multi-model files (NMR ensembles) give one conformation per model; several files give one conformation each. rnafold3d run folders (with metadata/manifest.json) keep their scores. RCV files (.rcv) and saved sessions open here too.' }));
    } else if (tab === 'server') {
      pane.append(el('div', { class: 'pathbar' }, pathIn, el('button', { class: 'btn', text: 'Go', onclick: () => browseTo(pathIn.value) })), places, list,
        el('div', { class: 'desc', text: 'Browse folders on the machine running the viewer (for example WSL paths or /mnt/c/...). Click files to select them, or open the whole folder.' }));
      pathIn.addEventListener('keydown', e => { if (e.key === 'Enter') browseTo(pathIn.value); });
      browseTo(cwd?.path || '');
    } else {
      const box = el('div', { class: 'filelist', style: { height: '360px' } }, el('div', { class: 'fl-item muted', text: 'Searching…' }));
      pane.append(box);
      app.api.get('/api/library').then(r => {
        box.replaceChildren();
        if (!r.packages.length) box.append(el('div', { class: 'fl-item muted', text: 'No rnafold3d runs found under ' + r.roots.join(', ') }));
        for (const p of r.packages) {
          const item = el('div', { class: 'fl-item', title: p.path },
            icon('cube'), el('div', { class: 'nm' }, el('div', { text: p.name, style: { fontWeight: 560 } }), el('div', { class: 'sz', text: p.path })),
            el('span', { class: 'sz', text: `${p.n_conformations} conf · ${p.length ?? '?'} nt` }));
          item.addEventListener('click', async () => { d.close(); await app.openServerPaths([p.path], opt()); });
          box.append(item);
        }
      }).catch(err => box.replaceChildren(el('div', { class: 'fl-item', text: err.message })));
    }
  };
  seg.addEventListener('click', e => { const t = e.target.closest('button')?.dataset.t; if (t) { tab = t; render(); } });
  if (!app.api.online && tab !== 'local') tab = 'local';
  for (const b of seg.children) if (b.dataset.t !== 'local') b.disabled = !app.api.online;
  render();
}

// ------------------------------------------------------------------ image export
const SIZES = [
  { id: 'view2', label: 'Current view × 2' },
  { id: 'col1', label: 'Journal single column (89 mm)', mm: 89 },
  { id: 'col15', label: 'Journal 1.5 columns (120 mm)', mm: 120 },
  { id: 'col2', label: 'Journal double column (183 mm)', mm: 183 },
  { id: 'slide', label: 'Slide 1920 × 1080', px: [1920, 1080] },
  { id: 'custom', label: 'Custom pixels' },
];

export function exportImageDialog(app) {
  const v = app.viewer;
  const prefs = app.prefs.imageExport || {};
  let format = prefs.format || 'png';
  const fmtSeg = el('div', { class: 'seg' }, ...['png', 'tiff', 'pdf'].map(f => el('button', { 'data-f': f, text: f.toUpperCase() })));
  const size = select(SIZES.map(s => ({ value: s.id, label: s.label })), prefs.size || 'col1');
  const aspect = select([{ value: 'view', label: 'Same as viewport' }, { value: '1', label: 'Square 1:1' }, { value: '1.3333', label: '4:3' }, { value: '1.7778', label: '16:9' }, { value: '0.75', label: '3:4 (portrait)' }], prefs.aspect || 'view');
  const dpi = select([{ value: '300', label: '300 dpi' }, { value: '600', label: '600 dpi' }, { value: '1200', label: '1200 dpi (line art)' }], String(prefs.dpi || 600));
  const wIn = el('input', { class: 'input', type: 'number', min: 64, max: 12000, value: prefs.w || 2400 });
  const hIn = el('input', { class: 'input', type: 'number', min: 64, max: 12000, value: prefs.h || 1800 });
  const bg = select([{ value: 'view', label: 'Viewport colour' }, { value: '#ffffff', label: 'White' }, { value: 'transparent', label: 'Transparent' }], prefs.bg || '#ffffff');
  const legend = el('input', { type: 'checkbox', checked: prefs.legend !== false ? true : null });
  const labels = el('input', { type: 'checkbox', checked: prefs.labels !== false ? true : null });
  const panel = el('input', { type: 'checkbox', checked: prefs.panel ? true : null });
  const info = el('div', { class: 'desc', style: { marginTop: 0 } });
  const prevImg = el('img', { alt: 'Export preview' });
  const preview = el('div', { class: 'preview-box' }, prevImg);
  const dims = () => {
    const [cw, ch] = v.cssSize;
    const ar = aspect.value === 'view' ? cw / ch : parseFloat(aspect.value);
    const s = SIZES.find(x => x.id === size.value);
    const D = +dpi.value;
    let w, h;
    if (s.id === 'view2') { w = Math.round(cw * 2); h = Math.round(w / ar); }
    else if (s.mm) { w = Math.round(s.mm / 25.4 * D); h = Math.round(w / ar); }
    else if (s.px) { [w, h] = s.px; }
    else { w = +wIn.value; h = +hIn.value; }
    return { w, h, D };
  };
  const refresh = debounce(() => {
    for (const b of fmtSeg.children) b.classList.toggle('on', b.dataset.f === format);
    const custom = size.value === 'custom';
    wIn.disabled = hIn.disabled = !custom;
    aspect.disabled = custom || size.value === 'slide';
    const { w, h, D } = dims();
    if (!custom) { wIn.value = w; hIn.value = h; }
    info.textContent = `${w.toLocaleString()} × ${h.toLocaleString()} px · ${(w / D * 25.4).toFixed(0)} × ${(h / D * 25.4).toFixed(0)} mm at ${D} dpi` +
      (panel.checked ? ' · plus energy panel' : '') + (format === 'pdf' ? ' · structure raster, panel & text vector' : '');
    try {
      const pw = 520, ph = Math.max(40, Math.round(pw * h / w));
      composeFigure(app, { width: pw, height: ph, transparent: bg.value === 'transparent', background: bg.value === 'view' || bg.value === 'transparent' ? null : bg.value, legend: legend.checked, labels: labels.checked })
        .then(c => { prevImg.src = c.toDataURL('image/png'); });
    } catch (err) { info.textContent = err.message; }
  }, 120);
  fmtSeg.addEventListener('click', e => { const f = e.target.closest('button')?.dataset.f; if (f) { format = f; refresh(); } });
  for (const c of [size, aspect, dpi, wIn, hIn, bg, legend, labels, panel]) c.addEventListener('input', refresh);
  dialog({
    title: 'Export image', width: 600,
    body: b => b.append(preview, info,
      el('div', { class: 'field', style: { marginTop: '10px' } }, el('span', { class: 'flabel', text: 'Format' }), fmtSeg),
      el('div', { class: 'grid2' }, field('Size', size), field('Aspect ratio', aspect)),
      el('div', { class: 'grid3' }, field('Resolution', dpi), field('Width (px)', wIn), field('Height (px)', hIn)),
      el('div', { class: 'grid2' }, field('Background', bg), el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Include' }),
        el('label', { class: 'check' }, legend, 'Colour legend'), el('label', { class: 'check' }, labels, 'Residue labels'),
        el('label', { class: 'check' }, panel, 'Energy panel beside the structure')))),
    actions: [{ label: 'Cancel', close: true }, { label: 'Export', primary: true, action: async btn => {
      const { w, h, D } = dims();
      app.prefs.imageExport = { format, size: size.value, aspect: aspect.value, dpi: D, w, h, bg: bg.value, legend: legend.checked, labels: labels.checked, panel: panel.checked };
      app.savePrefs();
      btn.disabled = true; btn.textContent = 'Rendering…';
      try {
        const r = await exportFigure(app, { format, width: w, height: h, dpi: D, transparent: bg.value === 'transparent', background: bg.value === 'view' || bg.value === 'transparent' ? null : bg.value, legend: legend.checked, labels: labels.checked, panel: panel.checked });
        const name = `${slug(app.ens.name)}_${app.ens.shortId(app.state.conf)}.${r.ext}`;
        const saved = await saveBlob(r.blob, name);
        if (saved) toast(`${saved} · ${(r.blob.size / 1e6).toFixed(1)} MB`, 'ok', 'Image exported');
      } finally { btn.disabled = false; btn.textContent = 'Export'; }
    } }],
  });
  refresh();
}

// ------------------------------------------------------------------ video export
export function exportVideoDialog(app, preset = {}) {
  const prefs = app.prefs.videoExport || {};
  const ens = app.ens;
  let content = preset.content || prefs.content || (ens && ens.n > 1 ? 'tour' : 'rotate');
  if (content === 'transition' && (!ens || ens.n < 2)) content = 'rotate';
  const cSeg = el('div', { class: 'seg' }, ...[['transition', 'Transition A → B'], ['tour', 'Conformation tour'], ['rotate', 'Rotate 360°'], ['rock', 'Rock'], ['tour-rotate', 'Tour + rotation']].map(([v, l]) => el('button', { 'data-c': v, text: l })));
  const confOpts = ens ? app.order.map(k => ({ value: String(k), label: `${ens.shortId(k)}${ens.value(k, app.state.metric) !== null ? ' · ' + fmt(ens.value(k, app.state.metric)) : ''}` })) : [];
  const defA = preset.from ?? (app.state.prev !== null && app.state.prev !== app.state.conf ? app.state.prev : app.order[0]);
  const defB = preset.to ?? (app.state.conf !== defA ? app.state.conf : app.order[1] ?? app.order[0]);
  const fromSel = select(confOpts, String(defA ?? 0));
  const toSel = select(confOpts, String(defB ?? 0));
  let route = prefs.route || 'direct';
  const routeSeg = el('div', { class: 'seg' }, el('button', { 'data-r': 'direct', text: 'Direct morph', title: 'Linear interpolation from A to B' }),
    el('button', { 'data-r': 'path', text: 'Through intermediates', title: 'Pass through the ensemble conformations on the lowest-RMSD route from A to B' }));
  const boomerang = el('input', { type: 'checkbox', checked: prefs.boomerang ? true : null });
  const spin = el('input', { class: 'input', type: 'number', step: 5, min: 0, max: 360, value: prefs.spinDeg ?? 0 });
  const inset2d = el('input', { type: 'checkbox', checked: prefs.inset2d ?? true ? true : null });
  const routeInfo = el('div', { class: 'desc', style: { marginTop: '4px' } });
  const res = select([{ value: '1280x720', label: '1280 × 720 (HD)' }, { value: '1920x1080', label: '1920 × 1080 (Full HD)' }, { value: '2560x1440', label: '2560 × 1440' }, { value: '1080x1080', label: '1080 × 1080 (square)' }, { value: '3840x2160', label: '3840 × 2160 (4K)' }], prefs.res || '1920x1080');
  const fps = select([{ value: '24', label: '24 fps' }, { value: '30', label: '30 fps' }, { value: '60', label: '60 fps' }], String(prefs.fps || 30));
  const trans = el('input', { class: 'input', type: 'number', step: 0.1, min: 0.2, max: 10, value: prefs.transition ?? app.morphDuration });
  const hold = el('input', { class: 'input', type: 'number', step: 0.1, min: 0, max: 10, value: prefs.hold ?? 1.0 });
  const secs = el('input', { class: 'input', type: 'number', step: 1, min: 2, max: 120, value: prefs.seconds ?? 12 });
  const rock = el('input', { class: 'input', type: 'number', step: 5, min: 5, max: 90, value: prefs.rockDeg ?? 30 });
  const loop = el('input', { type: 'checkbox', checked: prefs.loop !== false ? true : null });
  const overlay = el('input', { type: 'checkbox', checked: prefs.overlay !== false ? true : null });
  const legend = el('input', { type: 'checkbox', checked: prefs.legend !== false ? true : null });
  const fmts = videoFormats();
  const format = select(fmts.map(f => ({ value: f.id, label: f.label })), prefs.format && fmts.find(f => f.id === prefs.format) ? prefs.format : fmts[0].id);
  const info = el('div', { class: 'desc' });
  const prog = el('div', { class: 'progress', hidden: true }, el('i'));
  const opts = () => {
    const [w, h] = res.value.split('x').map(Number);
    return { content, width: w, height: h, fps: +fps.value, transition: +trans.value, hold: +hold.value, seconds: +secs.value, rockDeg: +rock.value,
      loop: loop.checked, overlay: overlay.checked, legend: legend.checked, format: format.value, bitrate: Math.round(w * h * +fps.value * 0.15),
      from: +fromSel.value, to: +toSel.value, route, boomerang: boomerang.checked, spinDeg: +spin.value, inset2d: inset2d.checked };
  };
  const transBox = el('div', { class: 'vt-box' },
    el('div', { class: 'grid2' }, field('From (A)', fromSel), field('To (B)', toSel)),
    el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Route' }), routeSeg), routeInfo,
    el('div', { class: 'grid2' }, field('Slow turn during the clip (°)', spin), el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Return' }),
      el('label', { class: 'check' }, boomerang, 'Play back to A (boomerang)'))));
  const refresh = () => {
    for (const b of cSeg.children) b.classList.toggle('on', b.dataset.c === content);
    for (const b of routeSeg.children) b.classList.toggle('on', b.dataset.r === route);
    const tour = content === 'tour' || content === 'tour-rotate', tr = content === 'transition';
    transBox.hidden = !tr;
    trans.disabled = !(tour || tr);
    hold.disabled = !(tour || tr);
    loop.disabled = !tour;
    secs.disabled = tour || tr;
    rock.disabled = content !== 'rock';
    const o = opts();
    if (tr && app.ens) {
      const r = o.route === 'path' ? transitionPath(app.ens, o.from, o.to) : [o.from, o.to];
      routeInfo.textContent = o.from === o.to ? 'Choose two different conformations.' : `${r.map(k => app.ens.shortId(k)).join(' → ')} · C4′ RMSD A–B ${fmt(app.ens.rmsd[o.from][o.to])} Å`;
    }
    const n = app.ens ? framePlan(app, o).length : 0;
    info.textContent = `${n} frames · ${(n / o.fps).toFixed(1)} s · recording takes about as long as the video plays (MP4/WebM).` + ((tour || tr) && app.ens && app.ens.n < 2 ? ' This ensemble has a single conformation.' : '');
  };
  cSeg.addEventListener('click', e => { const c = e.target.closest('button')?.dataset.c; if (c) { content = c; refresh(); } });
  routeSeg.addEventListener('click', e => { const r = e.target.closest('button')?.dataset.r; if (r) { route = r; refresh(); } });
  for (const c of [res, fps, trans, hold, secs, rock, loop, fromSel, toSel, boomerang, spin]) c.addEventListener('input', refresh);
  let abort = null;
  dialog({
    title: 'Record video', width: 600,
    body: b => b.append(
      el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Content' }), cSeg),
      transBox,
      el('div', { class: 'grid2' }, field('Resolution', res), field('Frame rate', fps)),
      el('div', { class: 'grid3' }, field('Transition (s)', trans), field('Hold (s)', hold), field('Rotation length (s)', secs)),
      el('div', { class: 'grid3' }, field('Rock angle (°)', rock), field('Format', format), el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Options' }),
        el('label', { class: 'check' }, loop, 'Loop back to start'), el('label', { class: 'check' }, overlay, 'Title & energy overlay'), el('label', { class: 'check' }, legend, 'Colour legend'),
        el('label', { class: 'check' }, inset2d, '2D structure inset'))),
      info, prog,
      el('div', { class: 'desc', text: 'The tour follows the conformation order chosen in the left panel. A transition morphs from A to B, directly or through the ensemble conformations on the lowest-RMSD route. Morphs are interpolations between superposed conformations — a visual aid, not a simulated pathway.' })),
    actions: [{ label: 'Close', close: true }, { label: 'Record', primary: true, action: async btn => {
      if (abort) { abort.abort(); return; }
      const o = opts();
      if (content === 'transition' && o.from === o.to) { toast('Choose two different conformations for A and B.', 'error'); return; }
      app.prefs.videoExport = { content, res: res.value, fps: o.fps, transition: o.transition, hold: o.hold, seconds: o.seconds, rockDeg: o.rockDeg, loop: o.loop, overlay: o.overlay, legend: o.legend, format: o.format, route: o.route, boomerang: o.boomerang, spinDeg: o.spinDeg, inset2d: o.inset2d };
      app.savePrefs();
      abort = new AbortController();
      prog.hidden = false;
      btn.textContent = 'Stop';
      try {
        const r = await recordVideo(app, o, { signal: abort.signal, onProgress: f => { prog.firstChild.style.width = `${Math.round(f * 100)}%`; } });
        const name = content === 'transition' ? `${slug(app.ens.name)}_${slug(app.ens.shortId(o.from))}_to_${slug(app.ens.shortId(o.to))}.${r.ext}` : `${slug(app.ens.name)}_${content}.${r.ext}`;
        const saved = await saveBlob(r.blob, name);
        if (saved) toast(`${saved} · ${r.frames} frames · ${(r.blob.size / 1e6).toFixed(1)} MB`, 'ok', 'Video saved');
      } catch (err) {
        if (err.name !== 'AbortError') throw err;
        toast('Recording cancelled', 'info');
      } finally { abort = null; btn.textContent = 'Record'; prog.hidden = true; }
    } }],
  });
  refresh();
}

// ------------------------------------------------------------------ preferences
export function preferencesDialog(app) {
  const s = app.state;
  const theme = select([{ value: 'system', label: 'Follow system' }, { value: 'light', label: 'Light' }, { value: 'dark', label: 'Dark' }], app.prefs.theme || 'system');
  const lang = select(LANGUAGES.map(l => ({ value: l.code, label: l.code === 'en' ? l.name : `${l.name} (${l.english})` })), currentLanguage(), { translate: 'no', 'aria-label': 'Interface language' });
  const quality = select([{ value: 'high', label: 'High (ambient occlusion, outlines)' }, { value: 'low', label: 'Fast (no ambient occlusion)' }], app.viewer.settings.quality);
  const proj = select([{ value: 'perspective', label: 'Perspective' }, { value: 'orthographic', label: 'Orthographic' }], app.viewer.projection);
  const labels = select([{ value: 'selection', label: 'Selected residues' }, { value: 'every10', label: 'Every 10th residue' }, { value: 'all', label: 'All residues' }, { value: 'none', label: 'None' }], app.viewer.labelMode);
  const outline = el('input', { type: 'color', value: app.viewer.settings.outlineColor, class: 'swatch-btn', style: { width: '60px' } });
  const uniform = el('input', { type: 'color', value: s.uniformColor, class: 'swatch-btn', style: { width: '60px' } });
  const elemCol = el('input', { type: 'checkbox', checked: app.viewer.settings.elementColors ? true : null });
  const devMode = el('input', { type: 'checkbox', checked: app.prefs.devMode ? true : null });
  const projDir = el('input', { class: 'input', value: app.info?.system?.projects_dir || '', 'aria-label': 'Projects folder', disabled: app.api.online ? null : true });
  const browseDir = el('button', { class: 'btn', type: 'button', hidden: !app.native?.available(), onclick: async () => { const d = await app.native.chooseFolder(); if (d) projDir.value = d; } }, 'Choose…');
  dialog({
    title: 'Preferences', width: 520,
    body: b => b.append(
      el('div', { class: 'grid2' }, field('Interface language', lang), field('Interface theme', theme)),
      el('div', { class: 'grid2' }, field('Rendering quality', quality), el('div')),
      el('div', { class: 'grid2' }, field('Projection', proj), field('Residue labels', labels)),
      el('div', { class: 'grid2' }, field('Outline colour', outline), field('Uniform colour scheme', uniform)),
      el('label', { class: 'check' }, elemCol, 'Colour N, O, P atoms by element in atom representations'),
      el('label', { class: 'check' }, devMode, 'Developer mode (for plugin authors): F12 / Ctrl+Shift+I open the browser developer tools, Shift+right-click shows the browser menu'),
      el('div', { class: 'field', style: { marginTop: '12px' } }, el('label', { text: 'Projects folder (new predictions are saved here)' }),
        el('div', { class: 'row', style: { gap: '8px' } }, projDir, browseDir))),
    actions: [{ label: 'Close', close: true }, { label: 'Apply', primary: true, action: async () => {
      if (app.api.online && projDir.value.trim() && projDir.value.trim() !== app.info?.system?.projects_dir) {
        app.info.system = await app.api.post('/api/settings', { projects_dir: projDir.value.trim() });
        toast(`New predictions are saved in ${app.info.system.projects_dir}`, 'ok', 'Projects folder');
      }
      app.setTheme(theme.value);
      if (lang.value !== currentLanguage()) await app.setLanguage(lang.value);
      app.viewer.settings.quality = quality.value;
      app.viewer.setProjection(proj.value);
      app.setLabelMode(labels.value);
      app.viewer.settings.outlineColor = outline.value;
      app.viewer.settings.elementColors = elemCol.checked;
      app.state.uniformColor = uniform.value;
      app.prefs.devMode = devMode.checked;
      app.prefs.quality = quality.value; app.prefs.projection = proj.value; app.prefs.outlineColor = outline.value; app.prefs.elementColors = elemCol.checked; app.prefs.uniformColor = uniform.value;
      app.savePrefs();
      app.refreshColors();
      app.viewer.needsRender = true;
    } }],
  });
}

// ------------------------------------------------------------------ plugins
export function pluginsDialog(app, info) {
  const host = app.plugins;
  const body = el('div');
  const render = () => {
    body.replaceChildren();
    const t = el('table', { class: 'table' }, el('thead', {}, el('tr', {}, el('th', { text: 'Plugin' }), el('th', { text: 'Status' }), el('th', { text: 'Enabled' }))));
    const tb = el('tbody');
    for (const p of host.plugins) {
      const cb = el('input', { type: 'checkbox', checked: p.enabled ? true : null, 'aria-label': `Enable ${p.name || p.file}` });
      cb.addEventListener('change', () => { host.setEnabled(p.url, cb.checked); note.hidden = false; });
      tb.append(el('tr', {}, el('td', {}, el('div', { style: { fontWeight: 560 }, text: `${p.name || p.file}${p.version ? ' ' + p.version : ''}` }), el('div', { class: 'muted', style: { fontSize: '11.5px' }, text: p.description || p.file || p.url }), p.error ? el('div', { class: 'state-bad', style: { fontSize: '11.5px' }, text: p.error }) : null),
        el('td', { text: `${p.status}${p.builtin ? ' · built-in' : ' · user'}` }), el('td', {}, cb)));
    }
    t.append(tb);
    body.append(t);
    const pyp = info?.python_plugins || [];
    body.append(el('div', { class: 'section-title', style: { marginTop: '14px' } }, 'Analysis plugins (Python)'),
      el('div', { class: 'desc', text: pyp.length ? pyp.map(p => `${p.name} (${p.path})`).join(' · ') : 'None installed.' }),
      ...(info?.plugin_errors || []).map(e => el('div', { class: 'warnbox', text: e })),
      el('div', { class: 'desc', text: `Plugin folder: ${(info?.plugin_dirs || ['the plugins folder'])[0]} — put *.js (interface) or *.py (analysis metrics) files there and restart RCV. See docs/RCV.md for the plugin API.` }));
  };
  const note = el('div', { class: 'warnbox', hidden: true, text: 'Reload the window to apply plugin changes.' });
  render();
  dialog({ title: 'Plugins', width: 620, body: b => b.append(note, body), actions: [{ label: 'Reload window', action: () => location.reload() }, { label: 'Close', primary: true, close: true }] });
}

// ------------------------------------------------------------------ help
export function shortcutsDialog() {
  const rows = [
    ['Drag', 'Rotate'], ['Right-drag / Shift+drag', 'Pan'], ['Wheel / pinch', 'Zoom'], ['Ctrl+drag', 'Roll'],
    ['Click', 'Select residue (Shift: add)'], ['Double-click', 'Centre on residue'], ['Right-click', 'Residue menu: select, show, colour, label, zoom, measure'], ['Esc', 'Clear selection / stop playback'],
    ['← →', 'Previous / next conformation (morph)'], ['Home / End', 'First / last conformation'], ['Space', 'Play / pause conformation tour'],
    ['1 – 5', 'Cartoon, cartoon + atoms, sticks, ball & stick, spheres'], ['C', 'Next colour scheme'], ['E', 'Ensemble overlay'], ['P', 'Flexibility putty'],
    ['L', 'Cycle residue labels'], ['F', 'Centre on selection'], ['R', 'Reset view'], ['A', 'Auto-rotate'], ['S', 'Quick PNG snapshot'],
    ['V', 'Select tool'], ['G / Alt+drag', 'Pull tool: drag a residue to switch conformation'], ['B', 'Paint tool: click or drag over residues to colour them'],
    ['T', 'Label tool: click a residue to label it; drag labels; double-click to edit'], ['M', 'Measure: click two atoms'], ['Delete', 'Remove the selected label or note'],
    ['U', 'Molecular surface'], ['D', '2D structure window'], ['`', 'Command console'],
    ['H', 'Only the 3D view (hide panels)'], ['F11', 'Full screen'], ['Esc', 'Restore a maximised window'], ['F1', 'User guide'],
    ['Ctrl+Z / Ctrl+Y', 'Undo / redo'], ['Ctrl+A', 'Select all residues'], ['Ctrl+F', 'Select residues (console)'], ['Ctrl+P', 'Export image'],
    ['Ctrl++ / Ctrl+− / Ctrl+0', 'Zoom the 3D view in / out / reset'], ['Ctrl+Shift+F', 'Fetch a structure from the PDB'],
    ['Ctrl+O', 'Open .rcv file or structures'], ['Ctrl+S', 'Save as .rcv file'], ['Ctrl+Shift+S', 'Save workspace'], ['Ctrl+E', 'Export image'], ['Ctrl+N', 'New prediction'], ['?', 'This list'],
  ];
  const t = el('table', { class: 'table kbd-table' }, el('tbody', {}, ...rows.map(([k, v]) => el('tr', {}, el('td', {}, ...k.split(' / ').flatMap((x, i) => [i ? ' / ' : '', el('kbd', { text: x })])), el('td', { text: v })))));
  dialog({ title: 'Keyboard and mouse', width: 520, body: t, actions: [{ label: 'Close', primary: true, close: true }] });
}

export function aboutDialog(app, info) {
  dialog({
    title: 'About RCV', width: 560,
    body: b => b.append(
      el('div', { class: 'about-brand' }, el('img', { src: 'img/icon.svg', alt: '' }), el('div', {}, el('b', { text: 'RCV' }), el('div', { class: 'muted', text: `RNA Conformation Viewer · version ${info?.version || '–'}` }))),
      el('div', { class: 'about-author' }, el('span', { class: 'muted', text: 'Author' }), el('b', { text: 'Zixiang Wang' }),
        el('span', { class: 'about-mail', text: 'princewang2018@gmail.com' }),
        el('button', { class: 'btn small ghost', onclick: () => navigator.clipboard?.writeText('princewang2018@gmail.com').then(() => toast('E-mail address copied.', 'ok')) }, icon('copy'), 'Copy')),
      el('p', { text: 'Predict RNA conformational ensembles from sequence on your own computer, and explore, compare and present them: morph and pull between conformations with the energy and the 2D structure following, analyse switches and mutations, and export publication figures, videos, .rcv files and interactive web pages.' }),
      el('p', { class: 'desc', text: 'Predictions come from the physics-inspired rnafold3d engine. Engine scores and evaluator energies are model scores for comparing conformations of one sequence; they are not experimental free energies. Morph animations are visual interpolations, not simulated pathways.' }),
      el('p', { class: 'desc', text: 'Copyright © 2026 Zixiang Wang. Free and open-source software (GPL-3.0). Rendering uses three.js (MIT License, © three.js authors). Secondary-structure energies use ViennaRNA (separate license) when installed. The command console uses RCV’s own command language; its verbs follow the common vocabulary of molecular viewers.' })),
    actions: [{ label: 'Close', primary: true, close: true }],
  });
}

// ------------------------------------------------------------------ methods text
export async function methodsDialog(app) {
  const ens = app.ens;
  if (!ens) return;
  let r;
  try {
    r = ens.serverId ? await app.api.get(`/api/documents/${ens.serverId}/methods`) : await app.api.post('/api/methods', { doc: ens.doc });
  } catch (err) { toast(err.message, 'error', 'Methods text'); return; }
  const extra = [];
  const st = app.state;
  const metric = ens.metricByKey[st.metric];
  const schemeLabel = (SCHEMES.find(s => s.id === st.colorScheme)?.label || st.colorScheme).toLowerCase();
  extra.push(`Figures show ${ens.shortId(st.conf)}${st.prev !== null && st.prev !== st.conf ? ` and ${ens.shortId(st.prev)}` : ''} rendered with RCV (${app.viewer.settings.style} style); residues are coloured by ${schemeLabel}${metric ? `, and energies are reported as ${metric.label} (${metric.unit || 'dimensionless'})` : ''}.`);
  const log = ens.commandLog || [];
  const box = (label, text, rows) => {
    const ta = el('textarea', { class: 'input mono methods-ta', rows, readonly: true });
    ta.value = text;
    return el('div', { class: 'field' }, el('div', { class: 'row', style: { alignItems: 'baseline' } }, el('label', { text: label }), el('span', { class: 'grow' }),
      el('button', { class: 'btn small', onclick: () => navigator.clipboard?.writeText(ta.value).then(() => toast(`${label} copied.`, 'ok')) }, icon('copy'), 'Copy')), ta);
  };
  dialog({
    title: 'Methods text & citation', width: 720,
    body: b => b.append(
      el('div', { class: 'desc', text: 'Built from what this ensemble records (engine version, parameters, seeds, validation, superposition). Check it against your own protocol before use.' }),
      box('Methods', `${r.methods} ${extra.join(' ')}`, 9),
      box('Software citation', r.citation, 3),
      log.length ? box('Commands used on this ensemble (reproducible figure recipe)', log.join('\n'), Math.min(8, log.length + 1)) : null),
    actions: [{ label: 'Close', primary: true, close: true }],
  });
}

// ------------------------------------------------------------------ multi-panel figure
export function panelsDialog(app) {
  const ens = app.ens;
  if (!ens) return;
  const prefs = app.prefs.panelsExport || {};
  let chosen = defaultPanelConfs(app, 4);
  let format = prefs.format || 'png';
  const listBox = el('div', { class: 'pn-list' });
  const renderList = () => {
    listBox.replaceChildren(...app.order.map(k => {
      const cb = el('input', { type: 'checkbox', checked: chosen.includes(k) ? true : null });
      cb.addEventListener('change', () => { chosen = cb.checked ? [...chosen, k] : chosen.filter(x => x !== k); chosen.sort((a, b) => app.order.indexOf(a) - app.order.indexOf(b)); refresh(); });
      const v = ens.value(k, app.state.metric);
      return el('label', { class: 'check pn-item' }, cb, el('b', { text: ens.shortId(k) }), el('span', { class: 'muted', text: `${v !== null ? fmt(v) : ''} · ${ens.conformations[k].dbn}` }));
    }));
  };
  const quick = el('div', { class: 'row', style: { gap: '6px', flexWrap: 'wrap', marginBottom: '6px' } },
    el('button', { class: 'btn small', onclick: () => { chosen = defaultPanelConfs(app, 6); renderList(); refresh(); } }, 'One per 2D state'),
    el('button', { class: 'btn small', onclick: () => { chosen = app.order.slice(0, 4); renderList(); refresh(); } }, 'First 4'),
    el('button', { class: 'btn small', onclick: () => { chosen = [app.state.prev ?? app.order[0], app.state.conf].filter((x, i, a) => a.indexOf(x) === i); renderList(); refresh(); } }, 'Previous + current'),
    el('button', { class: 'btn small ghost', onclick: () => { chosen = []; renderList(); refresh(); } }, 'None'));
  const fmtSeg = el('div', { class: 'seg' }, ...['png', 'tiff', 'pdf'].map(f => el('button', { 'data-f': f, text: f.toUpperCase() })));
  const width = select([{ value: '89', label: 'Single column (89 mm)' }, { value: '120', label: '1.5 columns (120 mm)' }, { value: '183', label: 'Double column (183 mm)' }], String(prefs.mm || 183));
  const aspect = select([{ value: '1', label: 'Square panels' }, { value: '1.3333', label: '4:3 panels' }, { value: '0.8', label: 'Portrait 4:5' }], String(prefs.aspect || 1));
  const cols = select([{ value: '0', label: 'Auto' }, { value: '1', label: '1' }, { value: '2', label: '2' }, { value: '3', label: '3' }, { value: '4', label: '4' }, { value: '6', label: '6' }], String(prefs.cols || 0));
  const dpi = select([{ value: '300', label: '300 dpi' }, { value: '600', label: '600 dpi' }], String(prefs.dpi || 600));
  const bg = select([{ value: '#ffffff', label: 'White' }, { value: 'view', label: 'Viewport colour' }, { value: 'transparent', label: 'Transparent (PNG/TIFF)' }], prefs.bg || '#ffffff');
  const cbox = (label, on) => { const c = el('input', { type: 'checkbox', checked: on ? true : null }); return [c, el('label', { class: 'check' }, c, label)]; };
  const [letters, lettersL] = cbox('Panel letters (A, B, …)', prefs.letters !== false);
  const [captions, captionsL] = cbox('Name and energy', prefs.captions !== false);
  const [inset, insetL] = cbox('2D structure inset', prefs.inset2d !== false);
  const [legend, legendL] = cbox('Colour legend', prefs.legend !== false);
  const [frame, frameL] = cbox('Thin frames', !!prefs.frame);
  const info = el('div', { class: 'desc' });
  const prevImg = el('img', { alt: 'Figure preview' });
  const opts = (scale = 1) => {
    const D = +dpi.value, mm = +width.value;
    return { confs: chosen.slice(0, 12), width: Math.round(mm / 25.4 * D * scale), aspect: +aspect.value, cols: +cols.value, dpi: D, format,
      background: bg.value === 'view' || bg.value === 'transparent' ? null : bg.value, transparent: bg.value === 'transparent' && format !== 'pdf',
      letters: letters.checked, captions: captions.checked, inset2d: inset.checked, legend: legend.checked, frame: frame.checked };
  };
  const refresh = debounce(async () => {
    for (const b of fmtSeg.children) b.classList.toggle('on', b.dataset.f === format);
    const o = opts();
    info.textContent = chosen.length ? `${Math.min(12, chosen.length)} panel${chosen.length === 1 ? '' : 's'} · ${o.width} px wide (${width.value} mm at ${o.dpi} dpi)${chosen.length > 12 ? ' · at most 12 panels' : ''}` : 'Choose conformations to show.';
    if (!chosen.length) { prevImg.removeAttribute('src'); return; }
    try {
      const c = await composePanels(app, { ...opts(), width: 560 });
      prevImg.src = c.toDataURL('image/png');
    } catch (err) { info.textContent = err.message; }
  }, 150);
  fmtSeg.addEventListener('click', e => { const f = e.target.closest('button')?.dataset.f; if (f) { format = f; refresh(); } });
  for (const c of [width, aspect, cols, dpi, bg, letters, captions, inset, legend, frame]) c.addEventListener('input', refresh);
  renderList();
  dialog({
    title: 'Figure panels (several conformations)', width: 760,
    body: b => b.append(el('div', { class: 'pn-grid' },
      el('div', {}, el('div', { class: 'flabel', text: 'Conformations (same orientation, current colours)' }), quick, listBox),
      el('div', {}, el('div', { class: 'preview-box pn-preview' }, prevImg), info,
        el('div', { class: 'field', style: { marginTop: '8px' } }, el('span', { class: 'flabel', text: 'Format' }), fmtSeg),
        el('div', { class: 'grid2' }, field('Figure width', width), field('Panels', aspect)),
        el('div', { class: 'grid3' }, field('Columns', cols), field('Resolution', dpi), field('Background', bg)),
        el('div', { class: 'pn-checks' }, lettersL, captionsL, insetL, legendL, frameL)))),
    actions: [{ label: 'Cancel', close: true }, { label: 'Export', primary: true, action: async btn => {
      if (!chosen.length) { toast('Choose at least one conformation.', 'info'); return; }
      const o = opts();
      app.prefs.panelsExport = { format, mm: +width.value, aspect: +aspect.value, cols: +cols.value, dpi: o.dpi, bg: bg.value, letters: o.letters, captions: o.captions, inset2d: o.inset2d, legend: o.legend, frame: o.frame };
      app.savePrefs();
      btn.disabled = true; btn.textContent = 'Rendering…';
      try {
        const r = await exportPanels(app, o);
        const saved = await saveBlob(r.blob, `${slug(ens.name)}_panels.${r.ext}`);
        if (saved) toast(`${saved} · ${r.canvas.width} × ${r.canvas.height} px · ${(r.blob.size / 1e6).toFixed(1)} MB`, 'ok', 'Figure exported');
      } finally { btn.disabled = false; btn.textContent = 'Export'; }
    } }],
  });
  refresh();
}

// ------------------------------------------------------------------ plugin development
export async function pluginGuideDialog(app) { return userGuideDialog(app, 'plugins'); }

// ------------------------------------------------------------------ help: user guide, what is new, website
const DOC_TITLES = { guide: 'User guide', changelog: 'What is new in this version', plugins: 'Plugin developer guide' };
const DOC_LINKS = { 'RCV.md': 'guide', 'PLUGINS.md': 'plugins', 'CHANGELOG.md': 'changelog' };

// Web pages open in the system browser (desktop window) or a new browser window.
export async function openExternal(app, url) {
  try {
    if (window.pywebview?.api?.open_url) { await window.pywebview.api.open_url(url); return; }
  } catch (err) { toast(err.message || String(err), 'error', 'Link'); return; }
  window.open(url, '_blank', 'noopener');
}
export function openWebsite(app) { return openExternal(app, app.info?.website || 'https://www.rnaconfviewer.com/'); }

// The documentation installed with RCV (English), with a table of contents. which: guide | changelog | plugins.
export async function userGuideDialog(app, which = 'guide') {
  const title = DOC_TITLES[which] || 'Help';
  let md;
  try { md = await (await app.api.raw('GET', `/api/docs/${which}`)).text(); }
  catch (err) { toast(app.api.online ? err.message : 'The help files are read from the local RCV engine, which is not reachable.', 'error', title); return; }
  const { renderMarkdown, slugify } = await import('./markdown.js');
  const urls = [];
  let d = null;
  const body = renderMarkdown(md, {
    image: async src => {
      const name = src.split(/[\\/]/).pop();
      const r = await app.api.raw('GET', `/api/docs/images/${encodeURIComponent(name)}`);
      const u = URL.createObjectURL(await r.blob());
      urls.push(u);
      return u;
    },
    link: href => {
      if (/^https?:/i.test(href)) return openExternal(app, href);
      if (/^mailto:/i.test(href)) return openExternal(app, href);
      const [file, anchor] = href.split('#');
      const doc = DOC_LINKS[file.split('/').pop()];
      if (file && doc && doc !== which) { d?.close(); userGuideDialog(app, doc); return; }
      if (anchor) body.querySelector(`#${CSS.escape(slugify(decodeURIComponent(anchor)))}`)?.scrollIntoView({ block: 'start' });
    },
  });
  body.setAttribute('translate', 'no');
  const heads = [...body.querySelectorAll('h3, h4')];
  const toc = el('nav', { class: 'guide-toc', 'aria-label': 'Contents', translate: 'no' },
    ...heads.map(h => el('button', { type: 'button', class: h.tagName === 'H4' ? 'l3' : '', text: h.textContent, onclick: () => h.scrollIntoView({ block: 'start' }) })));
  const online = { guide: 'guide/', changelog: 'guide/changelog.html', plugins: 'guide/plugins.html' }[which];
  d = dialog({
    title, width: 1120,
    body: b => b.append(el('div', { class: 'guide-grid' }, toc, el('div', { class: 'guide-body' }, body))),
    actions: [
      which === 'plugins' ? { label: 'New plugin from template…', action: () => newPluginDialog(app) } : null,
      { label: 'Online version', action: () => openExternal(app, (app.info?.website || 'https://www.rnaconfviewer.com/') + online) },
      { spacer: true }, { label: 'Close', primary: true, close: true }].filter(Boolean),
    onClose: () => urls.forEach(u => URL.revokeObjectURL(u)),
  });
  return d;
}

export function newPluginDialog(app) {
  const nameIn = el('input', { class: 'input', value: 'My plugin', maxlength: 40, 'aria-label': 'Plugin name' });
  let kind = 'js';
  const seg = el('div', { class: 'seg' }, el('button', { 'data-k': 'js', text: 'Interface plugin (JavaScript)' }), el('button', { 'data-k': 'py', text: 'Analysis plugin (Python)' }));
  const sync = () => { for (const b of seg.children) b.classList.toggle('on', b.dataset.k === kind); };
  seg.addEventListener('click', e => { const k = e.target.closest('button')?.dataset.k; if (k) { kind = k; sync(); } });
  sync();
  dialog({
    title: 'New plugin from template', width: 560,
    body: b => b.append(field('Name', nameIn), el('div', { class: 'field' }, el('span', { class: 'flabel', text: 'Kind' }), seg),
      el('div', { class: 'desc', text: 'Writes a commented, working template into your plugin folder. Edit it in any text editor, then reload the window (JavaScript) or restart RCV (Python).' })),
    actions: [{ label: 'Cancel', close: true }, { label: 'Create', primary: true, action: async () => {
      try {
        const r = await app.api.post('/api/plugins/new', { name: nameIn.value, kind });
        toast(`${r.path}`, 'ok', 'Plugin template created', 9000);
        if (app.native?.available()) app.native.reveal(r.path);
        document.querySelector('dialog[open]')?.close();
      } catch (err) { toast(err.message, 'error', 'New plugin'); }
    } }],
  });
}
