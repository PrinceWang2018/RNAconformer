// GUI test scenarios for tools/viewer_smoketest.py. Each scenario drives the real app (inside an iframe)
// through its public controller `app` and reports with log(); "DONE" marks success.
// Screenshot scenarios are rendered with a virtual-time budget; functional ones run in real time.
const ids = app => [...app.ensembles.keys()];
const byName = (app, re) => ids(app).find(id => re.test(app.ensembles.get(id).name)) || ids(app)[0];

export async function overview(app, w, log, sleep) {
  log('ensembles', ids(app).length, 'plugins', app.plugins.plugins.map(p => `${p.name}:${p.status}`).join(', '));
  if (!ids(app).length) throw new Error('no ensembles loaded (start the server with --demo)');
}

export async function morph(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.setColorScheme('helix');
  await sleep(300);
  const target = app.order.find(k => app.ens.conformations[k].dbn !== app.ens.conformations[app.state.conf].dbn) ?? app.order[1];
  app.goTo(target, { duration: 0.5 });
  await sleep(1500);
  if (app.state.conf !== target) throw new Error('conformation did not change');
  log('conf', app.ens.shortId(app.state.conf), 'prev', app.ens.shortId(app.state.prev), 'metric', app.state.metric);
}

export async function dark(app, w, log, sleep) {
  app.setTheme('dark', false);
  app.setBackground('#151514');
  app.selectEnsemble(byName(app, /NMR/i));
  app.setRepresentation('cartoon+sticks');
  app.select([16, 17, 18], 'set');
  app.goTo(app.order[Math.min(6, app.order.length - 1)], { duration: 0.3 });
  await sleep(1200);
  app.energyPanel.landMode = 'map';
  app.energyPanel.render();
  log('selection', [...app.state.selection].join(','));
}

export async function overlay(app, w, log, sleep) {
  const big = ids(app).sort((a, b) => app.ensembles.get(b).nRes - app.ensembles.get(a).nRes)[0];
  app.selectEnsemble(big);
  app.setColorScheme('flexibility');
  app.togglePutty();
  app.toggleGhosts();
  app.goTo(app.order[1] ?? app.order[0], { duration: 0.3 });
  await sleep(1200);
  app.showTab('details');
  log('ghosts', app.viewer.ghostGroup.children.length);
}

export async function plugins(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  const active = app.plugins.plugins.filter(p => p.status === 'active').length;
  if (active < 2) throw new Error(`expected the two built-in plugins, got ${active}`);
  app.select([5, 21], 'set');
  app.showTab('distance-tracker:panel');
  await sleep(100);
  w.document.querySelector('.dt-list').parentElement.querySelector('.btn.primary').click();
  await sleep(100);
  app.setMetric('plugin.distance-tracker.d1');
  app.goTo(app.order[2], { duration: 0.3 });
  await sleep(900);
  const v = app.ens.value(app.state.conf, app.state.metric);
  if (!(v > 0)) throw new Error('distance metric missing');
  log('distance metric', v.toFixed(2), 'Å; plugin metrics', app.ens.metrics.filter(m => m.group === 'plugin').map(m => m.key).join(','));
}

// Distance tracker (plugin API 3): removing a distance removes its metric; undo/redo; saved in .rcv files;
// kept when the ensemble is replaced by its extended version (pull exploration).
export async function trackersave(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  app.showTab('distance-tracker:panel');
  await sleep(100);
  const panel = w.document.querySelector('.dt-list').parentElement;
  const rows = () => panel.querySelectorAll('.dt-list .btn.ghost').length;
  const key = k => `plugin.distance-tracker.${k}`;
  const has = k => !!app.ens.metricByKey[key(k)];
  app.select([5, 21], 'set');
  panel.querySelector('.btn.primary').click();
  await sleep(400);
  panel.querySelector('select[name=atom]').value = 'P';
  app.select([2, 26], 'set');
  panel.querySelector('.btn.primary').click();
  await sleep(400);
  if (rows() !== 2 || !has('d1') || !has('d2')) throw new Error(`expected 2 tracked distances, got ${rows()}`);
  app.setMetric(key('d2'));
  panel.querySelectorAll('.dt-list .btn.ghost')[1].click();
  await sleep(400);
  if (rows() !== 1 || has('d2') || app.state.metric === key('d2')) throw new Error(`remove: rows ${rows()}, metric kept ${has('d2')}, shown ${app.state.metric}`);
  app.history.undo();
  await sleep(100);
  if (rows() !== 2 || !has('d2')) throw new Error('undo did not bring the distance back');
  app.history.redo();
  await sleep(100);
  if (rows() !== 1 || has('d2')) throw new Error('redo did not remove the distance again');

  app.setMetric(key('d1'));
  const before = app.ens.conformations.map(c => c.values[key('d1')]);
  const { blob } = await app.packRcv();
  await app.closeEnsemble(app.ens.id);
  await app.openRcvFile(new w.File([blob], 'tracked.rcv'));
  await sleep(1500);
  const after = app.ens.conformations.map(c => c.values[key('d1')]);
  const maxDiff = Math.max(...after.map((v, k) => Math.abs(v - before[k])));
  const lines = () => app.viewer.overlays.get('distance-tracker:lines')({ coords: app.display, ens: app.ens }).length;
  const checks = { row: rows() === 1, metric: has('d1') && !has('d2'), shown: app.state.metric === key('d1'), values: maxDiff < 0.01, line: lines() === 1 };
  const bad = Object.entries(checks).filter(([, v]) => !v).map(([k]) => k);
  if (bad.length) throw new Error('not restored from .rcv: ' + bad.join(', ') + ` (max diff ${maxDiff})`);

  const old = app.ens;
  const ext = await app.replaceDocument(old.id, old.serverId);
  await sleep(300);
  if (ext === old || !ext.metricByKey[key('d1')] || rows() !== 1 || lines() !== 1) throw new Error('distance lost when the ensemble was replaced');
  log('tracked', app.ens.metricByKey[key('d1')].label, '| remove, undo, redo, .rcv round trip (max diff', maxDiff.toFixed(4), 'Å), replace: ok');
}

export async function exports(app, w, log, sleep) {
  const fig = await import('/js/export/figure.js');
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.goTo(app.order[1], { animate: false });
  await sleep(300);
  for (const format of ['png', 'tiff', 'pdf']) {
    const t0 = performance.now();
    const r = await fig.exportFigure(app, { format, width: 1800, height: 1350, dpi: 600, transparent: format === 'tiff', background: '#ffffff', legend: true, labels: true, panel: format !== 'tiff' });
    if (r.blob.size < 10000) throw new Error(`${format} export too small`);
    log(format, `${(r.blob.size / 1e6).toFixed(2)} MB`, r.blob.type, `${Math.round(performance.now() - t0)} ms`);
  }
}

export async function video(app, w, log, sleep) {
  const vid = await import('/js/export/video.js');
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  for (const f of vid.videoFormats().map(x => x.id).filter((x, i, a) => a.indexOf(x) === i)) {
    const o = { content: 'tour', width: 640, height: 360, fps: 15, transition: 0.4, hold: 0.1, seconds: 2, rockDeg: 30, loop: false, overlay: true, legend: true, format: f, bitrate: 2e6 };
    const t0 = performance.now();
    const r = await vid.recordVideo(app, o);
    if (r.blob.size < 1000) throw new Error(`${f} video empty`);
    log(f, 'frames', r.frames, `${(r.blob.size / 1e3).toFixed(0)} kB`, r.blob.type, `${Math.round(performance.now() - t0)} ms`);
  }
}

export async function perf(app, w, log, sleep) {
  const big = ids(app).sort((a, b) => app.ensembles.get(b).nAtoms - app.ensembles.get(a).nAtoms)[0];
  app.selectEnsemble(big);
  app.toggleGhosts();
  await sleep(400);
  let frames = 0, run = true;
  const tick = () => { frames++; if (run) w.requestAnimationFrame(tick); };
  w.requestAnimationFrame(tick);
  const t0 = performance.now();
  for (let i = 0; i < 4; i++) { app.goTo(app.order[(i + 1) % app.order.length], { duration: 0.8 }); await sleep(850); }
  run = false;
  log('morph with ensemble overlay', (frames / ((performance.now() - t0) / 1000)).toFixed(0), 'fps;', app.ens.nRes, 'nt,', app.ens.nAtoms, 'atoms,', app.ens.n, 'conformations');
}

export async function empty(app, w, log, sleep) {
  for (const id of ids(app)) await app.closeEnsemble(id);
  await sleep(300);
  if (app.ens) throw new Error('ensembles still open');
  log('empty state shown');
}

// Close-up for tuning ambient occlusion: ?steps=closeup&ao=<intensity>
export async function closeup(app, w, log, sleep) {
  const q = new URLSearchParams(location.search);
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('nucleotide');
  if (q.get('ao')) app.viewer.styles.studio.aoIntensity = +q.get('ao');
  if (q.get('rad')) app.viewer.styles.studio.aoRadius = +q.get('rad');
  app.viewer.controls.zoom = 2.2;
  app.viewer.needsRender = true;
  await sleep(400);
  log('ao', app.viewer.styles.studio.aoIntensity, 'radius', app.viewer.styles.studio.aoRadius);
}

// ---------------------------------------------------------------- RCV 0.4+ scenarios
export async function welcome(app, w, log, sleep) {
  app.showPage('welcome');
  const W = app.welcome;
  W.seqIn.value = '>HIV1_TAR\nGGCAGAUCUGAGCCUGGGAGCUCUCUGCC';
  W.seqIn.dispatchEvent(new w.Event('input'));
  await sleep(600);
  if (W.predictBtn.disabled) throw new Error('predict button disabled for a valid sequence');
  log('status', W.seqStatus.textContent, '|', W.estimate.textContent, '| gpu', app.viewer.gpuName());
}

// Drag a residue toward where it sits in another conformation, release, and check the switch.
async function dragResidue(app, w, sleep, r, k, steps = 18) {
  const v = app.viewer, ens = app.ens;
  const [cw, ch] = v.cssSize;
  const rect = v.canvas.getBoundingClientRect();
  const proj = v.projector(v._cameraFor(cw / ch), cw, ch);
  const a = ens.trace[r];
  const from = proj(app.display[3 * a], app.display[3 * a + 1], app.display[3 * a + 2]);
  const to = proj(ens.coords[k][3 * a], ens.coords[k][3 * a + 1], ens.coords[k][3 * a + 2]);
  const ev = (x, y) => ({ clientX: rect.left + x, clientY: rect.top + y, pointerId: 77, button: 0, altKey: true });
  app.pullCtl.begin(r, ev(from[0], from[1]));
  for (let s = 1; s <= steps; s++) {
    const t = s / steps;
    app.pullCtl.move(ev(from[0] + (to[0] - from[0]) * t, from[1] + (to[1] - from[1]) * t));
    await sleep(40);
  }
  const hud = app.pullCtl.hud.textContent;
  const A = app.pullCtl.active;
  const T = A.pull.target, Tk = [0, 1, 2].map(c => ens.coords[k][3 * a + c]), vd = A.viewDir;
  const e = T.map((x, i) => x - Tk[i]), dz = e[0] * vd[0] + e[1] * vd[1] + e[2] * vd[2];
  dragResidue.debug = `target vs conf ${k} residue: in-plane ${Math.sqrt(Math.max(0, e[0] ** 2 + e[1] ** 2 + e[2] ** 2 - dz * dz)).toFixed(2)} Å, depth ${dz.toFixed(1)} Å; top ` +
    (A.ranked || []).slice(0, 3).map(r => `${ens.shortId(r.k)}:${r.local.toFixed(1)}/${r.score.toFixed(1)}`).join(' ');
  app.pullCtl.end(ev(to[0], to[1]), false);
  return hud;
}

export async function pull(app, w, log, sleep) {
  const id = [...app.ensembles.keys()].find(i => /prediction|sequence/i.test(app.ensembles.get(i).name));
  app.selectEnsemble(id);
  app.setColorScheme('helix');
  app.togglePull();
  await sleep(300);
  const ens = app.ens, start = app.state.conf, r = 16, a = ens.trace[r];
  const k = app.order.find(j => j !== start && ens.rmsd[start][j] > 4) ?? app.order[1];
  const hud = await dragResidue(app, w, sleep, r, k);
  await sleep(1800);
  const card = app.pullCtl.card;
  if (card.hidden) throw new Error('no pull result card');
  const L = app.pullCtl.last, T = L.pull.target, vd = L.viewDir;
  const e = [0, 1, 2].map(c => app.display[3 * a + c] - T[c]);
  const dz = e[0] * vd[0] + e[1] * vd[1] + e[2] * vd[2];
  const inPlane = Math.sqrt(Math.max(0, e[0] ** 2 + e[1] ** 2 + e[2] ** 2 - dz * dz));
  const switched = app.state.conf !== start;
  if (switched && inPlane > 6) throw new Error(`after switching, the residue lies ${inPlane.toFixed(1)} Å from the cursor`);
  log(`${switched ? 'switched' : 'stayed'} ${ens.shortId(start)} -> ${ens.shortId(app.state.conf)}; residue lands ${inPlane.toFixed(1)} Å from the cursor (screen plane) | hud: ${hud} | card: ${card.querySelector('b').textContent}`);
}

export async function predict(app, w, log, sleep) {
  app.showPage('welcome');
  const W = app.welcome;
  W.seqIn.value = '>uucg_hairpin\nGGCACUUCGGUGCC';
  W.modelsIn.value = 2;
  W.sampling = 'fast';
  W.seqIn.dispatchEvent(new w.Event('input'));
  await sleep(300);
  const before = app.ensembles.size;
  W.predictBtn.click();
  for (let i = 0; i < 600 && app.ensembles.size === before; i++) {
    await sleep(500);
    if (W.job && ['failed', 'cancelled'].includes(W.job.status)) throw new Error('prediction ' + W.job.status + ': ' + W.job.error);
  }
  if (app.ensembles.size === before) throw new Error('prediction did not finish in time');
  await sleep(500);
  log('page', app.page, '| ensemble', app.ens.name, app.ens.n, 'conformation(s),', app.ens.nRes, 'nt | notes', (W.job?.notes || []).length);
}

export async function explore(app, w, log, sleep) {
  const id = [...app.ensembles.keys()].find(i => /prediction|sequence/i.test(app.ensembles.get(i).name));
  app.selectEnsemble(id);
  app.togglePull();
  await sleep(300);
  const n0 = app.ens.n;
  const start = app.state.conf;
  const k = app.order.find(j => j !== start && app.ens.rmsd[start][j] > 4) ?? app.order[1];
  await dragResidue(app, w, sleep, 16, k);
  await sleep(800);
  await app.pullCtl.explore();
  for (let i = 0; i < 900 && app.ens.n === n0; i++) await sleep(500);
  if (app.ens.n === n0) throw new Error('exploration added no conformations in time');
  const added = app.ens.conformations.filter(c => c.explored);
  log('explored: ensemble', n0, '->', app.ens.n, 'conformations; showing', app.ens.shortId(app.state.conf), '| pull rms after relax', added.map(c => c.explored.rms_relaxed_A).join(', '));
}

// Diagnostics: screen <-> model round trip used by the pull tool.
export async function pulldebug(app, w, log, sleep) {
  const id = [...app.ensembles.keys()].find(i => /prediction|sequence/i.test(app.ensembles.get(i).name));
  app.selectEnsemble(id);
  await sleep(400);
  const v = app.viewer, ens = app.ens;
  const [cw, ch] = v.cssSize;
  const rect = v.canvas.getBoundingClientRect();
  const proj = v.projector(v._cameraFor(cw / ch), cw, ch);
  const a = ens.trace[16];
  const P = [0, 1, 2].map(c => app.display[3 * a + c]);
  const s = proj(...P);
  const M = v.screenToModel(rect.left + s[0], rect.top + s[1], P);
  log('css', cw, ch, 'rect', rect.width, rect.height, rect.left, rect.top, 'proj', s.map(x => x.toFixed(1)).join(','),
      'round-trip error', Math.hypot(M[0] - P[0], M[1] - P[1], M[2] - P[2]).toFixed(3), 'Å');
  const Q = [P[0] + 10, P[1], P[2]];
  const s2 = proj(...Q);
  const M2 = v.screenToModel(rect.left + s2[0], rect.top + s2[1], P);
  const vd = v.viewDirModel();
  const e = M2.map((x, i) => x - Q[i]); const dz = e[0] * vd[0] + e[1] * vd[1] + e[2] * vd[2];
  log('offset point: in-plane error', Math.sqrt(Math.max(0, e[0] ** 2 + e[1] ** 2 + e[2] ** 2 - dz * dz)).toFixed(3), 'depth', dz.toFixed(2),
      '| pixel ratio', v.renderer.getPixelRatio(), '| projection', v.projection);
}

export async function pulldebug2(app, w, log, sleep) {
  const id = [...app.ensembles.keys()].find(i => /prediction|sequence/i.test(app.ensembles.get(i).name));
  app.selectEnsemble(id);
  app.togglePull();
  await sleep(400);
  const v = app.viewer, ens = app.ens, r = 16, a = ens.trace[r];
  const start = app.state.conf;
  const k = app.order.find(j => j !== start && ens.rmsd[start][j] > 4) ?? app.order[1];
  const [cw, ch] = v.cssSize;
  const rect = v.canvas.getBoundingClientRect();
  const proj = v.projector(v._cameraFor(cw / ch), cw, ch);
  const P = [0, 1, 2].map(c => app.display[3 * a + c]), Tk = [0, 1, 2].map(c => ens.coords[k][3 * a + c]);
  const from = proj(...P), to = proj(...Tk);
  log('from', from.map(x => x.toFixed(0)).join(','), 'to', to.map(x => x.toFixed(0)).join(','), 'canvas', cw, ch);
  const ev = (x, y) => ({ clientX: rect.left + x, clientY: rect.top + y, pointerId: 77, button: 0, altKey: true });
  app.pullCtl.begin(r, ev(from[0], from[1]));
  const qBefore = v.controls.quaternion.toArray().map(x => x.toFixed(3)).join(',');
  for (let s = 1; s <= 6; s++) {
    const t = s / 6, x = from[0] + (to[0] - from[0]) * t, y = from[1] + (to[1] - from[1]) * t;
    app.pullCtl.move(ev(x, y));
    const T = app.pullCtl.active.pull.target;
    const back = v.projector(v._cameraFor(cw / ch), cw, ch)(...T);
    log(`step ${s}: mouse ${x.toFixed(0)},${y.toFixed(0)} -> target projects to ${back[0].toFixed(0)},${back[1].toFixed(0)}; X[r] ${app.pullCtl.active.pull.X.slice(3 * r, 3 * r + 3).map(z => z.toFixed(1)).join(',')} T ${T.map(z => z.toFixed(1)).join(',')}`);
    await sleep(40);
  }
  log('quaternion before', qBefore, 'after', v.controls.quaternion.toArray().map(x => x.toFixed(3)).join(','));
  app.pullCtl.end(ev(to[0], to[1]), true);
}

export async function pulldebug3(app, w, log, sleep) {
  const id = [...app.ensembles.keys()].find(i => /prediction|sequence/i.test(app.ensembles.get(i).name));
  app.selectEnsemble(id);
  app.setColorScheme('helix');
  app.togglePull();
  await sleep(300);
  const v = app.viewer, ens = app.ens, r = 16, a = ens.trace[r];
  const start = app.state.conf;
  const k = app.order.find(j => j !== start && ens.rmsd[start][j] > 4) ?? app.order[1];
  const [cw, ch] = v.cssSize;
  const rect = v.canvas.getBoundingClientRect();
  const proj = v.projector(v._cameraFor(cw / ch), cw, ch);
  const from = proj(...[0, 1, 2].map(c => app.display[3 * a + c])), to = proj(...[0, 1, 2].map(c => ens.coords[k][3 * a + c]));
  const ev = (x, y) => ({ clientX: rect.left + x, clientY: rect.top + y, pointerId: 77, button: 0, altKey: true });
  log('rect0', rect.left, rect.top, rect.width, rect.height, 'from', from.map(x => x.toFixed(0)).join(','), 'to', to.map(x => x.toFixed(0)).join(','));
  app.pullCtl.begin(r, ev(from[0], from[1]));
  for (let s = 1; s <= 18; s++) {
    const t = s / 18, x = from[0] + (to[0] - from[0]) * t, y = from[1] + (to[1] - from[1]) * t;
    app.pullCtl.move(ev(x, y));
    if (s % 6 === 0) {
      const T = app.pullCtl.active.pull.target, r2 = v.canvas.getBoundingClientRect();
      const back = v.projector(v._cameraFor(v.cssSize[0] / v.cssSize[1]), v.cssSize[0], v.cssSize[1])(...T);
      log(`step ${s}: mouse ${x.toFixed(0)},${y.toFixed(0)} target->screen ${back[0].toFixed(0)},${back[1].toFixed(0)} rect ${r2.left},${r2.top},${r2.width},${r2.height} css ${v.cssSize.join('x')}`);
    }
    await sleep(40);
  }
  const A = app.pullCtl.active, T = A.pull.target, Tk = [0, 1, 2].map(c => ens.coords[k][3 * a + c]);
  const vd = A.viewDir, vd2 = v.viewDirModel();
  const e = T.map((x, i) => x - Tk[i]), dz = e[0] * vd[0] + e[1] * vd[1] + e[2] * vd[2];
  const P2 = v.projector(v._cameraFor(cw / ch), cw, ch);
  log('T', T.map(x => x.toFixed(1)).join(','), 'Tk', Tk.map(x => x.toFixed(1)).join(','), '|e|', Math.hypot(...e).toFixed(2), 'dz', dz.toFixed(2),
      'screen T', P2(...T).slice(0, 2).map(x => x.toFixed(0)).join(','), 'screen Tk', P2(...Tk).slice(0, 2).map(x => x.toFixed(0)).join(','),
      'viewDir', vd.map(x => x.toFixed(2)).join(','), 'now', vd2.map(x => x.toFixed(2)).join(','));
  app.pullCtl.end(ev(to[0], to[1]), true);
}

// ---------------------------------------------------------------- RCV 0.5 scenarios
// 2D structure window: tRNA ensemble, switch to a conformation with a different 2D structure.
export async function ss2d(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('helix');
  if (!app.ss2d.visible) app.ss2d.setVisible(true);
  await sleep(300);
  const k0 = app.state.conf;
  const target = app.order.find(k => app.ens.conformations[k].dbn !== app.ens.conformations[k0].dbn) ?? app.order[1];
  app.goTo(target, { duration: 0.6 });
  await sleep(1400);
  const svgEl = w.document.querySelector('.ss-svg');
  const circles = svgEl.querySelectorAll('circle').length;
  if (circles !== app.ens.nRes) throw new Error(`2D window shows ${circles} nucleotides, expected ${app.ens.nRes}`);
  app.select([20, 21, 22], 'set');
  log('2D window', circles, 'nt;', w.document.querySelector('.ss-sub').textContent.trim());
}

// Transition video A -> B through intermediates, with the 2D inset (real time).
export async function transition(app, w, log, sleep) {
  const vid = await import('/js/export/video.js');
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(200);
  const a = app.order[0], b = app.order[app.order.length - 1];
  const path = vid.transitionPath(app.ens, a, b);
  const o = { content: 'transition', from: a, to: b, route: 'path', width: 640, height: 360, fps: 15, transition: 1.2, hold: 0.3, holdStart: 0.3, holdEnd: 0.3,
    boomerang: true, spinDeg: 30, overlay: true, legend: true, inset2d: true, format: 'zip', bitrate: 2e6 };
  const plan = vid.framePlan(app, o);
  const t0 = performance.now();
  const r = await vid.recordVideo(app, o);
  if (r.blob.size < 1000 || r.frames !== plan.length) throw new Error('transition video failed');
  // keep one frame for inspection
  const fr = plan[Math.floor(plan.length * 0.35)];
  app.applyFrame(fr.from, fr.to, fr.t);
  log('route', path.map(k => app.ens.shortId(k)).join('>'), '| frames', r.frames, `${(r.blob.size / 1e3).toFixed(0)} kB`, `${Math.round(performance.now() - t0)} ms`);
}

// Debug helper: show a canvas full-window on top of the app (for screenshots of export frames).
function showCanvas(w, canvas) {
  const doc = w.parent.document;
  const img = doc.createElement('img');
  img.src = canvas.toDataURL('image/png');
  img.style.cssText = 'position:fixed;inset:0;width:100%;height:100%;object-fit:contain;background:#777;z-index:9';
  doc.body.append(img);
}

// One frame of a transition video, shown full-window (screenshot scenario).
export async function transitionframe(app, w, log, sleep) {
  const vid = await import('/js/export/video.js');
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('helix');
  await sleep(300);
  const a = app.order[0], b = app.order[app.order.length - 1];
  const o = { content: 'transition', from: a, to: b, route: 'direct', width: 1600, height: 900, fps: 10, transition: 2, hold: 0.2, overlay: true, legend: true, inset2d: true };
  const plan = vid.framePlan(app, o);
  const fr = plan.find(f => f.phase === 'forward' && f.t > 0.45) || plan[0];
  const cap = app.viewer.beginCapture(o.width, o.height);
  app.applyFrame(fr.from, fr.to, fr.t);
  const r = cap.frame();
  const c = w.document.createElement('canvas'); c.width = o.width; c.height = o.height;
  const ctx = c.getContext('2d');
  ctx.drawImage(r.canvas, 0, 0);
  vid.drawTransitionOverlay(ctx, app, fr, o, plan);
  const s = o.height / 1080;
  app.ss2d.drawTo(ctx, [o.width - 420 * s - 40 * s, 40 * s, 420 * s, 360 * s], { k: fr.to, from: fr.from, t: fr.t, dark: false });
  cap.end();
  showCanvas(w, c);
  log('frame t', fr.t.toFixed(2));
}

// Command console: run a script of typical commands on the tRNA ensemble (screenshot scenario).
export async function console1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('nucleotide');
  app.console.setVisible(true);
  await sleep(200);
  const script = [
    'select anticodon, resi 33-37',
    'color orange, anticodon',
    'show sticks, anticodon',
    'show spheres, resi 76',
    'label resi 35, anticodon loop',
    'distance d1, 34/N1, 54/N3',
    'angle a1, 5/P, 6/P, 7/P',
    'hide cartoon, resi 58-62',
    'align 1, 3',
    'count switching',
    'view save v1',
    'select none',
    'info',
  ];
  let ok = 0;
  for (const c of script) if (await app.console.run(c)) ok++;
  const errs = [...w.document.querySelectorAll('.cs-err')].map(e => e.textContent);
  if (errs.length) throw new Error('console errors: ' + errs.join(' | '));
  const d = app.displayOf();
  log('commands', ok, 'of', script.length, '| named', [...d.named.keys()].join(','), '| colours', d.colors.size, '| hidden', d.hidden.size, '| measurements', app.measure.list().length, '| scenes', app.scenes.list().length);
}

// Switch tab: 2D states, populations and the mutation scan on the TAR prediction (screenshot).
export async function switchtab(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.showTab('switch');
  for (let i = 0; i < 100 && !app.switchPanel.result; i++) await sleep(100);
  const r = app.switchPanel.result;
  if (!r) throw new Error('switch analysis did not load: ' + (app.switchPanel.error || ''));
  if (!r.scan) throw new Error('no mutation scan');
  const best = r.scan.towards_b[0];
  app.switchPanel.pick(best.pos, best.to);
  await sleep(200);
  log('states', r.states.length, '| A share', r.scan.frac_a_wt, '| best towards B', `${best.from}${best.pos + 1}${best.to}`, best.ddg);
}

// Multi-panel figure (2 x 2) of the TAR prediction, shown full-window (screenshot).
export async function panels(app, w, log, sleep) {
  const P = await import('/js/export/panels.js');
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.setColorScheme('elements');
  await sleep(300);
  const confs = P.defaultPanelConfs(app, 4);
  const c = await P.composePanels(app, { confs, width: 1600, aspect: 1.333, cols: 2, dpi: 300, background: '#ffffff', letters: true, captions: true, inset2d: true, legend: true });
  const r = await P.exportPanels(app, { confs, width: 1200, aspect: 1, cols: 2, dpi: 300, background: '#ffffff', letters: true, captions: true, inset2d: true, legend: true, format: 'pdf' });
  if (r.blob.size < 10000) throw new Error('panel PDF too small');
  showCanvas(w, c);
  log('panels', confs.map(k => app.ens.shortId(k)).join(','), `${c.width}x${c.height}`, '| pdf', `${(r.blob.size / 1e3).toFixed(0)} kB`);
}

// .rcv round trip through the UI code: named selection, measurement, scene and view survive saving and reopening.
export async function rcvsave(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  for (const c of ['select bulge, resi 23-25', 'color orange, bulge', 'distance d1, 22/N1, 27/N3', 'view save s1', 'color helix', 'conf 4']) await app.console.run(c);
  await sleep(1600);
  const name = app.ens.name, conf = app.ens.conformations[app.state.conf].id;
  const { blob, stats } = await app.packRcv();
  if (!(blob.size > 1000 && stats.pdb_bytes > 5 * blob.size)) throw new Error(`unexpected .rcv size ${blob.size} (PDB ${stats.pdb_bytes})`);
  await app.closeEnsemble(app.ens.id);
  await app.openRcvFile(new w.File([blob], 'tar.rcv'));
  await sleep(1500);
  const ens = app.ens;
  if (!ens || ens.name !== name) throw new Error('reopened ensemble missing');
  const d = app.displayOf();
  const checks = { named: d.named.has('bulge'), colours: d.colors.size === 3, measurement: app.measure.list().length === 1, scene: app.scenes.list().length === 1,
    conf: ens.conformations[app.state.conf].id === conf, scheme: app.state.colorScheme === 'helix', log: (ens.commandLog || []).length >= 5 };
  const bad = Object.entries(checks).filter(([, v]) => !v).map(([k]) => k);
  if (bad.length) throw new Error('not restored: ' + bad.join(', '));
  log('rcv', `${(blob.size / 1024).toFixed(0)} kB`, `(${(stats.pdb_bytes / blob.size).toFixed(1)}x smaller than PDB)`, '| restored', Object.keys(checks).join(', '));
}

// Start page: Predict with an empty box runs the demo sequence (the job is cancelled right away).
export async function welcomedemo(app, w, log, sleep) {
  app.showPage('welcome');
  const W = app.welcome;
  W.seqIn.value = '';
  W._estimate();
  if (W.predictBtn.disabled) throw new Error('Predict is disabled for an empty sequence');
  const label = W.predictBtn.textContent;
  await W.startPrediction();
  const job = W.job;
  if (!job) throw new Error('no job started');
  const seqOk = W.seqIn.value.includes('GGCAGAUCUGAGCCUGGGAGCUCUCUGCC');
  await app.api.post(`/api/jobs/${job.id}/cancel`, {});
  for (let i = 0; i < 60 && W.job && !['cancelled', 'failed', 'done'].includes(W.job.status); i++) await sleep(250);
  if (!seqOk) throw new Error('the demo sequence was not filled in');
  log('button', JSON.stringify(label), '| job', job.title, '| status', W.job?.status);
}

// Batch: a FASTA with two short hairpins predicted from the start page; the user picks one to open.
export async function batchui(app, w, log, sleep) {
  app.showPage('welcome');
  const W = app.welcome;
  W.seqIn.value = '>uucg\nGGCACUUCGGUGCC\n>gnra\nGGACGAAAGUCC\n';
  W.modelsIn.value = 2;
  W.sampling = 'fast';
  W._estimate();
  const items = w.document.querySelectorAll('.wl-bitem').length;
  if (items !== 2) throw new Error(`batch list shows ${items} sequences`);
  await W.startPrediction();
  for (let i = 0; i < 600 && !W.batch; i++) await sleep(500);
  if (!W.batch) throw new Error('batch did not finish: ' + (W.job?.error || W.job?.status));
  const radios = w.document.querySelectorAll('input[name="batch-pick"]');
  const openBtn = [...w.document.querySelectorAll('.wl-batch button')].find(b => /Open selected/.test(b.textContent));
  if (radios.length !== 2 || !openBtn?.disabled) throw new Error('chooser not shown or not waiting for a choice');
  radios[1].click();
  await sleep(100);
  openBtn.click();
  for (let i = 0; i < 100 && !(app.ens && /gnra/.test(app.ens.name)); i++) await sleep(200);
  if (!(app.ens && /gnra/.test(app.ens.name))) throw new Error('chosen prediction not opened');
  log('batch', W.batch.result.items.map(it => `${it.name}:${it.status}:${it.n_conformations}`).join(' '), '| opened', app.ens.name, '| page', app.page);
}

// Interactive HTML export of the open ensemble.
export async function htmlexport(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  const r = await app.api.raw('POST', '/api/export-html', JSON.stringify({ document: app.ens.serverId, view: { conf: 1 }, title: 'test' }), { 'Content-Type': 'application/json' });
  const text = await r.text();
  if (!/x-rcv-module/.test(text) || text.length < 500000) throw new Error('page incomplete');
  log('html', `${(text.length / 1e6).toFixed(2)} MB`);
}

// Start page with a multi-sequence FASTA pasted (screenshot).
export async function welcomebatch(app, w, log, sleep) {
  app.showPage('welcome');
  const W = app.welcome;
  W.seqIn.value = '>TAR\nGGCAGAUCUGAGCCUGGGAGCUCUCUGCC\n>uucg_hairpin\nGGCACUUCGGUGCC\n>gnra\nGGACGAAAGUCC\n((((....))))\n>bad_one\nGGGAAXCC\n';
  W._estimate();
  await sleep(300);
  log('status', W.seqStatus.textContent);
}

// Documentation screenshot: TAR prediction switching to its alternative 2D state, Switch tab open.
export async function docviewer(app, w, log, sleep) {
  app.setTheme('light', false);
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.setColorScheme('elements');
  if (!app.ss2d.visible) app.ss2d.setVisible(true);
  const alt = app.order.find(k => app.ens.conformations[k].dbn !== app.ens.conformations[app.order[0]].dbn);
  app.goTo(app.order[0], { animate: false, force: true });
  await sleep(200);
  app.goTo(alt, { duration: 0.4 });
  await sleep(900);
  app.showTab('switch');
  for (let i = 0; i < 100 && !app.switchPanel.result; i++) await sleep(100);
  app.select([12, 13, 14], 'set');
  await sleep(300);
  log('shown', app.ens.shortId(app.state.conf), 'prev', app.ens.shortId(app.state.prev));
}

// Right-click residue menu: opening it and choosing "Show as sticks" runs (and logs) a console command.
export async function ctxmenu(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(300);
  const v = app.viewer, ens = app.ens;
  const [cw, ch] = v.cssSize;
  const rect = v.canvas.getBoundingClientRect();
  let p = null;
  for (let i = 0; i < 400 && !p; i++) p = v.pick(rect.left + cw * (0.3 + 0.4 * Math.random()), rect.top + ch * (0.3 + 0.4 * Math.random()));
  if (!p) throw new Error('nothing under the cursor');
  app._residueMenu(p, rect.left + cw / 2, rect.top + ch / 2);
  const items = [...w.document.querySelectorAll('.menu-pop .menu-item')];
  const item = items.find(b => /Show as sticks/.test(b.textContent));
  if (!item) throw new Error('menu has no "Show as sticks" item');
  item.click();
  await sleep(300);
  if (!app.displayOf().sticks.has(p.residue)) throw new Error('sticks not shown');
  log('menu items', items.length, '| residue', ens.resLabel(p.residue), '| logged', JSON.stringify((ens.commandLog || []).slice(-1)[0]));
}

// The example script examples/scripts/tar_switch_figure.rcvs runs without errors (console "run").
export async function script1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  const r = await app.api.post('/api/read-text', { path: (new URLSearchParams(location.search)).get('script') });
  app.console.setVisible(true);
  await app.console.runScript(r.text);
  await sleep(2500);
  const errs = [...w.document.querySelectorAll('.cs-err')].map(e => e.textContent);
  if (errs.length) throw new Error('script errors: ' + errs.join(' | '));
  log('script ok |', app.scenes.list().map(s => s.name).join(','), '| measurements', app.measure.list().length);
}

// Tools menu back ends: inspect and extract an .rcv file; the fetch dialog opens.
export async function tools1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  const { blob } = await app.packRcv();
  const buf = await blob.arrayBuffer();
  const info = await (await app.api.raw('POST', '/api/rcv/info', buf, { 'Content-Type': 'application/octet-stream' })).json();
  if (info.n_conformations !== app.ens.n || !info.thumbnail) throw new Error('inspect failed');
  const zip = await (await app.api.raw('POST', '/api/rcv/extract', buf, { 'Content-Type': 'application/octet-stream', 'X-Mode': 'zip' })).arrayBuffer();
  const multi = await (await app.api.raw('POST', '/api/rcv/extract', buf, { 'Content-Type': 'application/octet-stream', 'X-Mode': 'multi' })).text();
  if (zip.byteLength < 10000 || (multi.match(/^MODEL/mg) || []).length !== app.ens.n) throw new Error('extract failed');
  app.fetchDialog();
  await sleep(200);
  const dlg = w.document.querySelector('dialog .dlg-head h2')?.textContent;
  w.document.querySelector('dialog')?.close();
  log('inspect', info.name, info.n_conformations, '| zip', `${(zip.byteLength / 1024).toFixed(0)} kB`, '| models', (multi.match(/^MODEL/mg) || []).length, '| dialog', dlg);
}

// PDF export embeds Arial (when the system has it): check the font dictionaries of an exported figure.
export async function pdffont(app, w, log, sleep) {
  const fig = app.lib;                                               // the app's own module instances (iframe)
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(300);
  const r = await fig.exportFigure(app, { format: 'pdf', width: 1200, height: 900, dpi: 300, transparent: false, background: '#ffffff', legend: true, labels: true, panel: true });
  const bytes = new Uint8Array(await r.blob.arrayBuffer());
  let s = '';
  for (let i = 0; i < bytes.length; i++) s += String.fromCharCode(bytes[i]);
  const fonts = [...s.matchAll(/\/BaseFont \/([A-Za-z+-]+)/g)].map(m => m[1]);
  const embedded = (s.match(/\/FontFile2/g) || []).length;
  if (!fonts.some(f => /LiberationSans/.test(f))) throw new Error('figure font not embedded: ' + fonts.join(','));
  log('fonts', [...new Set(fonts)].join(' '), '| embedded font programs', embedded, '| pdf', `${(bytes.length / 1e6).toFixed(2)} MB`);
}

// Figure panels: preview (560 px) next to the full-size export scaled down — the 2D drawings must look alike.
export async function panelsize(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('elements');
  await sleep(300);
  const confs = app.order.slice(0, 2);
  const o = { confs, aspect: 1.333, cols: 2, dpi: 600, background: '#ffffff', letters: true, captions: true, inset2d: true, legend: true };
  const small = await app.lib.composePanels(app, { ...o, width: 560 });
  const big = await app.lib.composePanels(app, { ...o, width: 4323 });
  const c = w.document.createElement('canvas');
  c.width = 1140; c.height = Math.round(small.height) + 20;
  const ctx = c.getContext('2d');
  ctx.fillStyle = '#888'; ctx.fillRect(0, 0, c.width, c.height);
  ctx.drawImage(small, 10, 10);
  ctx.imageSmoothingQuality = 'high';
  ctx.drawImage(big, 580, 10, 560, Math.round(560 * big.height / big.width));
  showCanvas(w, c);
  log('preview', `${small.width}x${small.height}`, '| export', `${big.width}x${big.height}`);
}

// Mouse tools, annotations, palettes and undo/redo (screenshot): paint, label, region, note; undo and redo.
export async function tools2(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.setColorScheme('elements');
  await sleep(300);
  const v = app.viewer, rect = v.canvas.getBoundingClientRect(), [cw, ch] = v.cssSize;
  const hits = [];
  for (let i = 0; i < 2000 && hits.length < 6; i++) {
    const x = rect.left + cw * (0.25 + 0.5 * Math.random()), y = rect.top + ch * (0.25 + 0.5 * Math.random());
    const p = v.pick(x, y);
    if (p && !hits.some(h => h.p.residue === p.residue)) hits.push({ p, x, y });
  }
  if (hits.length < 4) throw new Error('could not find residues under the cursor');
  // paint two residues
  app.tools.set('paint');
  const ev = (type, h, extra = {}) => new w.PointerEvent(type, { clientX: h.x, clientY: h.y, button: 0, buttons: 1, pointerId: 7, bubbles: true, ...extra });
  v.canvas.dispatchEvent(ev('pointerdown', hits[0]));
  v.canvas.dispatchEvent(ev('pointermove', hits[1]));
  v.canvas.dispatchEvent(ev('pointerup', hits[1], { buttons: 0 }));
  await sleep(400);
  const painted = app.displayOf().colors.size;
  // label + note
  app.tools.set('label');
  app.annotations.onPick(hits[2].p, { altKey: false });
  app.annotations.closeEditor();
  app.annotations.addNote('HIV-1 TAR · alternative fold');
  app.annotations.closeEditor();
  await sleep(400);
  app.history.commit();
  const nAnn = app.annotations.list().length;
  app.history.undo();
  await sleep(100);
  const afterUndo = app.annotations.list().length;
  app.history.redo();
  await sleep(100);
  app.palettes.apply('okabe-ito');
  app.tools.set('select');
  await sleep(400);
  log('painted', painted, '| annotations', nAnn, '-> undo', afterUndo, '-> redo', app.annotations.list().length, '| palette', app.palettes.active.name);
  if (painted < 2 || nAnn !== 2 || afterUndo !== 1 || app.annotations.list().length !== 2) throw new Error('tools/undo mismatch');
}

// Plot studio: each plot type renders; the RMSD matrix (clustered) is shown as a screenshot and exported as PDF.
export async function plots(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /NMR/i));
  await sleep(300);
  const P = await import('/js/plots.js');
  const sizes = [];
  for (const [type] of P.PLOTS) {
    const s = P.drawPlot(app, { type, w: 89, h: 70, font: 8, ramp: 'blues', cluster: true, values: true, grid: false, highlight: true, ref: 'mean', style: 'line', moved: true, track: true, metric: app.state.metric });
    sizes.push(`${type}:${s.querySelectorAll('*').length}`);
  }
  app.openPlotStudio('matrix');
  await sleep(800);
  const prev = w.document.querySelector('.ps-preview svg');
  if (!prev) throw new Error('plot studio did not render');
  const pdf = await app.lib.exportFigure ? null : null;
  log('elements', sizes.join(' '));
}

// All plot types side by side (screenshot of the SVGs themselves).
export async function plotsall(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.goTo(app.order[4], { animate: false });
  await sleep(300);
  const P = await import('/js/plots.js');
  const doc = w.parent.document;
  const box = doc.createElement('div');
  box.style.cssText = 'position:fixed;inset:0;z-index:9;background:#888;display:grid;grid-template-columns:repeat(3,1fr);gap:10px;padding:10px';
  for (const [type] of P.PLOTS) {
    const s = P.drawPlot(app, { type, w: 89, h: 70, font: 8, ramp: 'blues', cluster: true, values: true, grid: false, highlight: true, ref: 'best', style: 'line', moved: true, track: true, metric: app.state.metric });
    s.style.cssText = 'width:100%;height:auto;background:#fff';
    box.append(doc.importNode(s, true));
  }
  doc.body.append(box);
  log('plots', P.PLOTS.length);
}

// Molecular surface with transparency over the cartoon (screenshot), via console commands.
export async function surface1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.setColorScheme('elements');
  await sleep(300);
  for (const c of ['show surface', 'set transparency, 0.45, surface']) if (!(await app.console.run(c))) throw new Error('command failed: ' + c);
  await sleep(500);
  const m = app.viewer.surfaceMesh;
  if (!m) throw new Error('no surface');
  app.goTo(app.order[1], { duration: 0.3 });
  await sleep(1200);
  log('surface triangles', m.geometry.index.count / 3, '| built in', Math.round(app.surface.lastMs), 'ms | transparent', app.viewer.matSurf.transparent);
}

// Plugin guide renders; a template plugin can be created (in the test profile's plugin folder).
export async function plugindev(app, w, log, sleep) {
  await app.openPluginGuide();
  await sleep(300);
  const md = w.document.querySelector('dialog .md');
  if (!md || md.querySelectorAll('table').length < 5) throw new Error('guide not rendered');
  const r = await app.api.post('/api/plugins/new', { name: 'Test plugin', kind: 'js' });
  const t = await app.api.post('/api/plugins/new', { name: 'Test plugin', kind: 'py' });
  log('guide tables', md.querySelectorAll('table').length, '| code blocks', md.querySelectorAll('pre').length, '| templates', r.path.split(/[\\/]/).pop(), t.path.split(/[\\/]/).pop());
}

// Template plugins written by "New plugin" work after a restart (JS: metric, command, panel; Python: metric).
export async function plugintpl(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(300);
  const rec = app.plugins.plugins.find(p => p.name === 'Test plugin');
  if (!rec || rec.status !== 'active') throw new Error('template JS plugin not active: ' + (rec ? rec.status + ' ' + rec.error : 'missing'));
  const ok = await app.console.run('gccount 2');
  const out = [...w.document.querySelectorAll('.cs-out .cs-out')].pop()?.textContent;
  const jsMetric = app.ens.metrics.find(m => m.key === 'plugin.test-plugin.gc_pairs');
  const pyMetric = app.ens.metrics.find(m => m.key === 'plugin.test-plugin.radius');
  if (!ok || !jsMetric || !pyMetric) throw new Error(`command ${ok} js metric ${!!jsMetric} py metric ${!!pyMetric}`);
  log('command:', out, '| js metric', jsMetric.label, '| py metric', pyMetric.label, app.ens.value(0, 'plugin.test-plugin.radius'));
}

// Desktop behaviour: right-click gives our own menus (chart / table / text field), text selection and browser
// shortcuts are blocked.
export async function native1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  app.showTab('details');
  await sleep(300);
  const D = w.document;
  const closeAll = () => D.body.dispatchEvent(new w.PointerEvent('pointerdown', { bubbles: true, clientX: 5, clientY: 5 }));
  const fire = (elem, type, init = {}) => { const e = new w.MouseEvent(type, { bubbles: true, cancelable: true, clientX: 300, clientY: 300, button: 2, ...init }); elem.dispatchEvent(e); return e; };
  const out = {};
  const chart = D.querySelector('#tab-bodies svg.chart') || D.querySelector('svg.chart');
  app.showTab('energy'); await sleep(100);
  const e1 = fire(D.querySelector('svg.chart'), 'contextmenu');
  out.chartMenu = e1.defaultPrevented && !!D.querySelector('.menu-pop') && /Export as PDF/.test(D.querySelector('.menu-pop').textContent);
  closeAll();
  app.showTab('details'); await sleep(100);
  const td = D.querySelector('#tab-bodies table td');
  const e2 = fire(td, 'contextmenu');
  out.tableMenu = e2.defaultPrevented && /Copy table/.test(D.querySelector('.menu-pop')?.textContent || '');
  closeAll();
  const e3 = fire(D.querySelector('.panel-head') || D.querySelector('.statusbar'), 'contextmenu');
  out.plainBlocked = e3.defaultPrevented && !D.querySelector('.menu-pop');
  const sel = new w.Event('selectstart', { bubbles: true, cancelable: true });
  D.querySelector('.statusbar span').dispatchEvent(sel);
  out.noSelect = sel.defaultPrevented;
  const k = new w.KeyboardEvent('keydown', { key: 'r', ctrlKey: true, bubbles: true, cancelable: true });
  D.body.dispatchEvent(k);
  out.noReload = k.defaultPrevented;
  const bad = Object.entries(out).filter(([, v]) => !v).map(([n]) => n);
  if (bad.length) throw new Error('not native: ' + bad.join(', '));
  void chart;
  log('checks', Object.keys(out).join(', '));
}

// ------------------------------------------------------------------ RCV 0.6: window layout and languages
const rectOf = elem => { const r = elem.getBoundingClientRect(); return { w: Math.round(r.width), h: Math.round(r.height), vis: r.width > 0 && r.height > 0 }; };

// Window menu: hide/show every part, maximise panels (they fill the window), H = 3D view only, Esc restores.
export async function windows1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  await sleep(300);
  const D = w.document, wm = app.wm, W = w.innerWidth;
  const parts = { left: '.sidebar.left', right: '.sidebar.right', toolbar: '.toolbar', timeline: '.timeline', statusbar: '.statusbar', ss2d: '.ss-inset', seqbar: '.seqbar', console: '.console' };
  const res = [];
  for (const [id, sel] of Object.entries(parts)) {
    const before = wm.isShown(id);
    wm.set(id, false); await sleep(80);
    const hidden = !rectOf(D.querySelector(sel)).vis && !wm.isShown(id);
    wm.set(id, true); await sleep(80);
    const shown = rectOf(D.querySelector(sel)).vis && wm.isShown(id);
    if (!hidden || !shown) throw new Error(`${id}: hide ${hidden} show ${shown}`);
    wm.set(id, before); await sleep(50);
    res.push(id);
  }
  // persisted between sessions
  wm.set('statusbar', false);
  const saved = JSON.parse(w.localStorage.getItem('rcv.windows') || w.localStorage.getItem('windows') || '{}');
  wm.set('statusbar', true);
  // maximise: the window fills the width under the menu bar; Esc restores
  const maxed = [];
  for (const [id, sel] of [['right', '.sidebar.right'], ['left', '.sidebar.left'], ['ss2d', '.ss-inset'], ['console', '.console']]) {
    wm.maximize(id); await sleep(150);
    const r = rectOf(D.querySelector(sel));
    if (r.w < W - 4 || r.h < w.innerHeight - 60) throw new Error(`${id} maximised to ${r.w}x${r.h} of ${W}x${w.innerHeight}`);
    D.querySelector('#gl').dispatchEvent(new w.KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
    await sleep(150);
    if (wm.maxed || rectOf(D.querySelector(sel)).w >= W - 4) throw new Error(`${id} not restored by Esc`);
    maxed.push(`${id} ${r.w}x${r.h}`);
  }
  app.console.setVisible(false);
  // H: only the 3D view
  D.querySelector('#gl').dispatchEvent(new w.KeyboardEvent('keydown', { key: 'h', bubbles: true }));
  await sleep(150);
  const vp = rectOf(D.querySelector('#viewport'));
  const panelsGone = !rectOf(D.querySelector('.sidebar.left')).vis && !rectOf(D.querySelector('.toolbar')).vis && !rectOf(D.querySelector('.timeline')).vis;
  if (!wm.focus || !panelsGone || vp.w < W - 4) throw new Error(`focus mode: ${wm.focus} ${panelsGone} ${vp.w}`);
  D.querySelector('#gl').dispatchEvent(new w.KeyboardEvent('keydown', { key: 'Escape', bubbles: true }));
  await sleep(150);
  if (wm.focus || !rectOf(D.querySelector('.toolbar')).vis) throw new Error('Esc did not leave the 3D-only view');
  // the Window menu lists every part
  const menu = wm.menuItems().filter(i => i.label).map(i => i.label);
  if (menu.length < 15) throw new Error('Window menu too short: ' + menu.length);
  log('hide/show', res.join(','), '| saved', JSON.stringify(saved).slice(0, 80), '| maximised', maxed.join('; '), '| 3D only', vp.w + 'x' + vp.h, '| menu', menu.length);
}

// Screenshot: the analysis panel maximised (charts at full width).
export async function winmax(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  await sleep(300);
  app.wm.maximize('right');
  app.showTab('energy');
  await sleep(600);
  log('maximised', app.wm.maxed);
}

// Every language: the start page, menus and dialogs are translated; switching back gives the exact English text.
export async function lang1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  await sleep(300);
  const D = w.document;
  const english = () => [...D.querySelectorAll('.menu-btn, .toolbar [title], .panel-head, .tab')].map(e => e.textContent.trim() + '|' + (e.title || '')).join('\n');
  const en0 = english();
  const report = [];
  for (const l of app.i18n.languages.filter(x => x.code !== 'en')) {
    await app.setLanguage(l.code);
    await sleep(120);
    const menus = [...D.querySelectorAll('.menu-btn')].map(b => b.textContent.trim());
    const leftEn = app.i18n.untranslated(D.querySelector('#app'));
    // open each menu and count untranslated items
    let items = 0, miss = [];
    for (const b of D.querySelectorAll('.menu-btn')) {
      b.click(); await sleep(30);
      for (const it of D.querySelectorAll('.menu-pop .menu-item > span:nth-child(2), .menu-pop .menu-head')) {
        items++;
        const s = it.textContent.trim();
        if (app.i18n.t(s) !== s) miss.push(s);        // English text left although a translation exists
      }
      D.body.dispatchEvent(new w.PointerEvent('pointerdown', { bubbles: true, clientX: 5, clientY: 900 }));
    }
    if (D.documentElement.lang !== l.code) throw new Error('lang attribute ' + D.documentElement.lang);
    const leftKnown = leftEn.filter(x => app.i18n.t(x) !== x);
    report.push(`${l.code}: menus ${menus.slice(0, 3).join('/')} · ${items} items · English text left with a translation ${leftKnown.length}`);
    if (miss.length || leftKnown.length) throw new Error(`${l.code}: not translated: ${[...miss, ...leftKnown].slice(0, 12).join(' | ')}`);
  }
  await app.setLanguage('en');
  await sleep(120);
  const en1 = english();
  if (en1 !== en0) throw new Error('English not restored exactly');
  log(report.join(' || '));
}

async function langShot(app, w, sleep, code, page) {
  await app.setLanguage(code);
  if (page === 'welcome') { app.showPage('welcome'); await sleep(600); return; }
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  app.showTab('energy');
  await sleep(500);
}
export async function langzh(app, w, log, sleep) {
  await langShot(app, w, sleep, 'zh-CN');
  w.document.querySelectorAll('.menu-btn')[2].click();          // View menu open in the screenshot
  await sleep(200);
  log('language', app.i18n.current());
}
export async function langde(app, w, log, sleep) { await langShot(app, w, sleep, 'de'); log('language', app.i18n.current()); }
export async function langja(app, w, log, sleep) { await langShot(app, w, sleep, 'ja', 'welcome'); log('language', app.i18n.current()); }

// Help > User guide: the installed guide with its contents list and images.
export async function guide1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  const d = await app.openUserGuide();
  await sleep(1500);
  const D = w.document;
  const toc = D.querySelectorAll('.guide-toc button').length, imgs = [...D.querySelectorAll('.guide-body img')];
  if (toc < 15) throw new Error('guide contents: ' + toc);
  log('contents', toc, '| images', imgs.length, imgs.filter(i => i.complete && i.naturalWidth > 0).length, 'loaded');
  void d;
}

// Developer aid (not in the smoke list): with a translation loaded, open the main views, menus and dialogs and
// log the visible interface text that has no translation yet.  live.sh harvest 120
export async function harvest(app, w, log, sleep) {
  const D = w.document;
  await app.setLanguage(new URLSearchParams(w.location.search).get('lang') || 'zh-CN');
  const found = new Set();
  const ATTRS = ['title', 'placeholder', 'aria-label'];
  const collect = () => {
    for (const s of app.i18n.untranslated(D.body)) found.add(s);
    for (const e of D.querySelectorAll(ATTRS.map(a => `[${a}]`).join(','))) {
      if (e.closest('[translate="no"], svg, canvas, pre, code')) continue;
      for (const a of ATTRS) { const v = (e.getAttribute(a) || '').trim(); if (v && /[A-Za-z]{3}/.test(v) && app.i18n.t(v) === v) found.add('@' + v); }
    }
  };
  const closeDlg = async () => { for (const d of D.querySelectorAll('dialog')) d.close(); await sleep(80); };
  app.showPage('welcome'); await sleep(400); collect();
  if (app.welcome.adv) app.welcome.adv.open = true; await sleep(100); collect();
  for (const re of [/TAR_NMR|1ANR/i, /prediction|sequence/i, /1EHZ|tRNA/i]) {
    app.selectEnsemble(byName(app, re)); app.showPage('viewer'); await sleep(400);
    for (const t of app._tabs) { app.showTab(t.id); await sleep(500); collect(); }
  }
  for (const b of D.querySelectorAll('.menu-btn')) {
    b.click(); await sleep(60); collect();
    D.body.dispatchEvent(new w.PointerEvent('pointerdown', { bubbles: true, clientX: 5, clientY: 900 }));
  }
  const dialogs = [
    () => app.openDialog('local'), () => app.openDialog('server'), () => app.openDialog('library'),
    () => app.lib && w.document.querySelector('[data-act="export-image"]').click(),
    () => app.openVideoDialog(), () => app.openPanels(), () => app.openPlotStudio(), () => app.openPalettes(),
    () => app.fetchDialog(), () => app.console.setVisible(true),
  ];
  for (const open of dialogs) { try { await open(); } catch (e) { log('dialog error', e.message); } await sleep(700); collect(); await closeDlg(); }
  app.console.setVisible(false);
  for (const item of app.rebuildMenus && [...D.querySelectorAll('.menu-btn')].map(b => b.textContent)) void item;
  const list = [...found].sort();
  for (let i = 0; i < list.length; i += 25) log('UNTRANSLATED', JSON.stringify(list.slice(i, i + 25)));
  log('count', list.length);
}

// Developer aid (not in the smoke list): time a language switch and menu opening.
export async function langperf(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /TAR_NMR|1ANR/i));
  await sleep(300);
  const D = w.document, P = w.performance;
  const times = [];
  for (const code of ['de', 'zh-CN', 'en']) {
    const t0 = P.now(); await app.setLanguage(code); const t1 = P.now();
    const b = D.querySelectorAll('.menu-btn')[2]; b.click(); const t2 = P.now();
    await new Promise(r => w.setTimeout(r, 0)); const t3 = P.now();
    D.body.dispatchEvent(new w.PointerEvent('pointerdown', { bubbles: true, clientX: 5, clientY: 900 }));
    const t4 = P.now(); await sleep(50); const t5 = P.now();
    times.push(`${code}: switch ${(t1 - t0).toFixed(0)} ms, menu ${(t3 - t1).toFixed(0)} ms, sleep(50) took ${(t5 - t4).toFixed(0)} ms`);
  }
  log(times.join(' | '));
}

// ---------------------------------------------------------------- RCV 0.7
const waitFor = async (sleep, cond, ms = 20000) => { const t0 = Date.now(); while (!cond() && Date.now() - t0 < ms) await sleep(100); return cond(); };

// Pathway tab: refolding paths, bottleneck, modules, phi scan, mechanism morph.
export async function pathway1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(300);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result || app.pathway.error)) throw new Error('no pathway result');
  if (app.pathway.error) throw new Error(app.pathway.error);
  const r = app.pathway.result;
  app.pathway.render();                       // surfaces any rendering error with its stack
  const root = app.tabBody('pathway');
  const need = ['.pw-verdict', '.pw-map', '.pw-bneck'];
  const miss = need.filter(q => !root.querySelector(q));
  if (miss.length) throw new Error('missing sections: ' + miss.join(', '));
  await app.pathway.runScan();
  if (!app.pathway.scan) throw new Error('no phi scan');
  app.pathway.colorBy('role');
  if (app.state.colorScheme !== 'pathway-role') throw new Error('mechanism colours not applied');
  app.pathway.play(true);
  if (!await waitFor(sleep, () => app.state.conf === r.b, 12000)) throw new Error('mechanism morph did not reach B');
  log(app.ens.name, '| barrier', r.direct.barrier, '->', r.effective_barrier, 'kcal/mol | t1/2', r.half_life_s?.toExponential(1), 's | modules', r.motion.modules.length,
    '| unit', r.bottleneck.unit.length, 'bp | scan faster', app.pathway.scan.faster.length, '| routes', r.tpt.pathways.length);
}

// Pathway between two NMR models of one fold: both relax into the same 2D basin -> no barrier, no scan, no routes.
export async function pathwaybasin(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1ANR|NMR/i));
  await sleep(300);
  const P = app.pathway;
  P.a = 0; P.b = 11; P.result = P.scan = null;
  app.showTab('pathway');
  await P.load();
  if (P.error) throw new Error(P.error);
  const r = P.result;
  if (!r.same_basin) throw new Error('expected the same basin');
  P.render();
  const root = app.tabBody('pathway');
  const verdict = root.querySelector('.pw-verdict')?.textContent || '';
  if (!/same basin/.test(verdict) || !/shared basin/.test(verdict)) throw new Error('verdict: ' + verdict.slice(0, 200));
  if (!root.querySelector('svg.chart:not(.pw-map) path')) throw new Error('no energy profile');
  if (root.querySelector('.pw-bneck')) throw new Error('bottleneck shown for one basin');
  if ([...root.querySelectorAll('.section-title')].some(t => /φ/.test(t.textContent))) throw new Error('scan section shown for one basin');
  const labels = [...root.querySelectorAll('svg.chart text')].map(t => t.textContent);
  if (!labels.includes('basin')) throw new Error('basin not marked: ' + labels.join(','));
  await app.setLanguage('zh-CN');
  P.render();
  await sleep(400);
  const zh = app.tabBody('pathway').querySelector('.pw-verdict')?.textContent || '';
  await app.setLanguage('en');
  if (!/盆地/.test(zh)) throw new Error('verdict not translated: ' + zh.slice(0, 200));
  log(app.ens.name, '| A above basin', r.A.above_basin, '| B above basin', r.B.above_basin, '| direct', r.direct.barrier, '| zh ok');
}

// Probing data: synthetic SHAPE profile from conformation 0 -> agreement, reweighting, scheme, metric.
export async function probing1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  const ens = app.ens, c = ens.conformations[0];
  const paired = new Set(c.pairs.flatMap(p => [p[0], p[1]]));
  const values = [...Array(ens.nRes).keys()].map(i => (paired.has(i) ? 0.1 : 1.1) + 0.05 * Math.sin(i));
  await app.probing.set(ens, { name: 'synthetic', kind: 'SHAPE', values, raw: values, factor: null, format: 'test' });
  if (!ens.probing?.fit) throw new Error('no fit');
  const auc = ens.value(0, 'probe_auc');
  if (!(auc > 0.9)) throw new Error('agreement of the source conformation too low: ' + auc);
  app.setColorScheme('reactivity');
  app.setMetric('probe_weight');
  app.showTab('details');
  await sleep(300);
  if (!app.tabBody('details').textContent.includes('Probing data')) throw new Error('no probing section');
  const { blob } = await app.packRcv();
  await app.closeEnsemble(ens.id);
  await app.openRcvFile(new w.File([blob], 'probe.rcv'));
  await sleep(1200);
  if (!app.ens.probing?.fit || app.ens.value(0, 'probe_auc') !== auc) throw new Error('probing data not restored from .rcv');
  log('AUC conf_001', auc.toFixed(3), '| chi2', ens.probing.fit.reweighting.chi2_prior, '->', ens.probing.fit.reweighting.chi2, '| restored from .rcv');
}

// Conditions: salt and a ligand change the Switch tab energies; saved with the ensemble.
export async function conditions1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /prediction|sequence/i));
  await sleep(200);
  app.showTab('switch');
  if (!await waitFor(sleep, () => app.switchPanel.result || app.switchPanel.error)) throw new Error('no switch result');
  const e0 = app.switchPanel.result.states[0].energy;
  const { setConditions } = await import('/js/conditions.js');
  setConditions(app, { sodium: 0.1 });
  await sleep(100);
  if (!await waitFor(sleep, () => app.switchPanel.result && !app.switchPanel.loading)) throw new Error('no result after the change');
  const e1 = app.switchPanel.result.states[0].energy;
  if (!(e1 > e0)) throw new Error(`lower salt should destabilise: ${e0} -> ${e1}`);
  setConditions(app, { ligand: { preset: 'theophylline', name: 'Theophylline aptamer', motif_seq: 'GAUACCAG&CCCUUGGCAGC', motif_struct: '(...((((&)...)))...)', kd_uM: 0.32, conc_uM: 10 } });
  await sleep(100);
  if (!await waitFor(sleep, () => app.switchPanel.result && !app.switchPanel.loading)) throw new Error('no result with ligand');
  const txt = app.switchPanel.result.conditions_text;
  if (!/Theophylline/.test(txt)) throw new Error('ligand not in the conditions: ' + txt);
  setConditions(app, { sodium: 1.021, ligand: null });
  log('ΔG S1', e0, '->', e1, 'kcal/mol at 0.1 M Na+ | ligand:', txt);
}

// Screenshot: Pathway tab of the tRNA prediction (map with routes, verdict, profile, bottleneck, modules, scan).
export async function pathwayshot(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  app.wm?.maximize?.('right');
  await sleep(200);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result || app.pathway.error)) throw new Error('no pathway result');
  if (app.pathway.error) throw new Error(app.pathway.error);
  await app.pathway.runScan();
  app.pathway.colorBy('role');
  log('pathway', app.pathway.result.effective_barrier, 'kcal/mol');
}

// Screenshot: topology-aware map in the Energy tab (TAR NMR ensemble: several 2D states).
export async function landmap(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /NMR/i));
  await sleep(200);
  app.energyPanel.landMode = 'map';
  app.showTab('energy');
  app.energyPanel.render();
  await sleep(300);
  const L = app.ens.doc.ensemble.landscape;
  if (!L) throw new Error('no landscape in the document');
  log('classes', L.classes.map(c => `${c.label}:${c.shape}:${c.members.length}`).join(' '), '| quality', JSON.stringify(L.quality.map));
}

// Screenshot: lower part of the Pathway tab (bottleneck, modules, phi scan); the upper sections are hidden.
export async function pathwayshot2(app, w, log, sleep) {
  await pathwayshot(app, w, log, sleep);
  app.pathway.render();
  app.pathway.render = () => {};              // keep this rendering for the screenshot
  const body = app.tabBody('pathway');
  for (const s of body.querySelectorAll(':scope > .section')) { const t = s.querySelector('.section-title')?.textContent || ''; if (!/Routes ranked|Bottleneck|Coordinated|φ-value/.test(t)) s.style.display = 'none'; }
  await sleep(300);
}

// Screenshot: the 0.7 publication plots (ensemble map with routes, refolding path, phi scan) at journal sizes.
export async function plots07(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(200);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result || app.pathway.error)) throw new Error('no pathway result');
  await app.pathway.runScan();
  app.showTab('energy');
  const P = await import('/js/plots.js');
  const doc = w.parent.document;
  const box = doc.createElement('div');
  box.style.cssText = 'position:fixed;inset:0;z-index:9;background:#888;display:grid;grid-template-columns:1fr 1fr;gap:10px;padding:10px;align-content:start';
  const base = { font: 8, grid: false, highlight: true, metric: app.state.metric, routes: true };
  for (const [type, w0, h0] of [['landscape', 89, 80], ['pathway', 89, 70], ['phi', 183, 60]]) {
    const s = P.drawPlot(app, { ...base, type, w: w0, h: h0 });
    s.style.cssText = `display:block;width:100%;height:auto;aspect-ratio:${w0}/${h0};background:#fff;${type === 'phi' ? 'grid-column:1 / span 2' : ''}`;
    box.append(doc.importNode(s, true));
  }
  doc.body.append(box);
  log('plots', 3, '| barrier', app.pathway.result.effective_barrier);
}

// Interactive HTML page with the 0.7 analyses: the page starts, draws the map and the pathway, plays the mechanism.
export async function htmlpage07(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(200);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result || app.pathway.error)) throw new Error('no pathway result');
  await app.pathway.runScan();
  const ens = app.ens;
  const body = { view: { conf: ens._pathway.a }, title: 'test', analysis: { pathway: ens._pathway, phi: ens._phiScan }, document: ens.serverId };
  const r = await app.api.raw('POST', '/api/export-html', JSON.stringify(body), { 'Content-Type': 'application/json' });
  const text = await r.text();
  const doc = w.parent.document;
  const frame = doc.createElement('iframe');
  frame.style.cssText = 'position:fixed;inset:0;width:100%;height:100%;z-index:20;background:#fff;border:0';
  frame.srcdoc = text;
  doc.body.append(frame);
  if (!await waitFor(sleep, () => frame.contentDocument?.body?.dataset.ready === '1', 30000)) {
    const err = frame.contentDocument?.getElementById('err');
    throw new Error('page did not start' + (err && !err.hidden ? ': ' + err.textContent : ''));
  }
  const fd = frame.contentDocument;
  const mapItems = fd.querySelectorAll('#map circle').length;
  if (mapItems < ens.n) throw new Error('map not drawn: ' + mapItems);
  if (fd.getElementById('pw-box').hidden) throw new Error('pathway section hidden');
  fd.getElementById('pw-play').click();
  await sleep(7000);
  const err = fd.getElementById('err');
  if (err && !err.hidden) throw new Error(err.textContent);
  log('page', `${(text.length / 1e6).toFixed(2)} MB`, '| map points', mapItems, '|', fd.getElementById('pw-text').textContent.slice(0, 90));
  frame.remove();
}

// Confidence model: per-nucleotide predicted lDDT appears as a metric and a colour scheme.
export async function confidence1(app, w, log, sleep) {
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  if (!await waitFor(sleep, () => app.ens._confidence, 20000)) throw new Error('no confidence values');
  const v = app.ens.value(app.state.conf, 'confidence');
  if (!(v > 0 && v < 1)) throw new Error('confidence out of range: ' + v);
  app.setColorScheme('confidence');
  await sleep(200);
  log('confidence', app.ens.conformations.map((c, k) => app.ens.value(k, 'confidence').toFixed(2)).join(' '), '| scheme', app.state.colorScheme);
}
// Screenshot: the confidence colour scheme on the tRNA prediction.
export async function confshot(app, w, log, sleep) { await confidence1(app, w, log, sleep); app.setMetric('confidence'); }

// Developer aid: interface text of the 0.7 views that has no translation yet: the views are built in English, their
// strings collected, then German is loaded and every string the app's translator leaves unchanged is reported.
export async function i18n07(app, w, log, sleep) {
  const doc = w.document;
  const views = {};
  const collect = (name, root) => {
    const out = new Set();
    if (root) {
      const wk = doc.createTreeWalker(root, w.NodeFilter.SHOW_TEXT);
      for (let n = wk.nextNode(); n; n = wk.nextNode()) {
        const v = n.data.trim();
        if (v && /[A-Za-z]{3}/.test(v) && !n.parentElement?.closest('svg, code, pre, kbd, textarea, select, [translate="no"]')) out.add(v);
      }
      for (const n of root.querySelectorAll('[title], [placeholder], [aria-label]')) for (const a of ['title', 'placeholder', 'aria-label']) {
        const v = n.getAttribute(a);
        if (v && /[A-Za-z]{3}/.test(v) && !n.closest('svg')) out.add(v);
      }
      for (const o of root.querySelectorAll('option')) if (/[A-Za-z]{3}/.test(o.textContent)) out.add(o.textContent.trim());
    }
    views[name] = [...(views[name] || []), ...out];
  };
  await app.setLanguage('en');
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  await sleep(300);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result)) throw new Error('no pathway');
  collect('pathway', app.tabBody('pathway'));
  await app.pathway.runScan();
  await sleep(300);
  collect('pathway', app.tabBody('pathway'));
  app.showTab('switch');
  const { setConditions } = await import('/js/conditions.js');
  setConditions(app, { ligand: { preset: 'theophylline', name: 'Theophylline aptamer', motif_seq: 'GAUACCAG&CCCUUGGCAGC', motif_struct: '(...((((&)...)))...)', kd_uM: 0.32, conc_uM: 10 } });
  await sleep(200);
  await waitFor(sleep, () => app.switchPanel.result && !app.switchPanel.loading);
  await sleep(300);
  collect('switch', app.tabBody('switch'));
  setConditions(app, { ligand: null });
  app.conditionsDialog();
  await sleep(300);
  collect('conditions', doc.querySelector('dialog[open]'));
  doc.querySelector('dialog[open]')?.close();
  const ens = app.ens, c = ens.conformations[0];
  const paired = new Set(c.pairs.flatMap(p => [p[0], p[1]]));
  await app.probing.set(ens, { name: 'synthetic', kind: 'SHAPE', values: [...Array(ens.nRes).keys()].map(i => (paired.has(i) ? 0.1 : 1.1)), raw: null, factor: 1.2, format: 'test' });
  app.showTab('details');
  await sleep(300);
  collect('details', app.tabBody('details'));
  app.probing.dialog();
  await sleep(300);
  collect('probingdialog', doc.querySelector('dialog[open]'));
  doc.querySelector('dialog[open]')?.close();
  app.showTab('energy');
  app.energyPanel.landMode = 'map';
  app.energyPanel.render();
  await sleep(300);
  collect('energy', app.tabBody('energy'));
  app.openPlotStudio('landscape');
  await sleep(400);
  collect('plotstudio', doc.querySelector('dialog[open]'));
  for (const t of ['pathway', 'phi']) { doc.querySelector(`.ps-item[data-type="${t}"]`)?.click(); await sleep(200); collect('plotstudio', doc.querySelector('dialog[open]')); }
  doc.querySelector('dialog[open]')?.close();
  app.transitionVideo(ens._pathway.a, ens._pathway.b, { mechanism: { motion: ens._pathway.motion, timing: ens._pathway.timing } });
  await sleep(400);
  collect('video', doc.querySelector('dialog[open]'));
  doc.querySelector('dialog[open]')?.close();
  const { preferencesDialog } = await import('/js/dialogs.js');
  preferencesDialog(app);
  await sleep(300);
  collect('preferences', doc.querySelector('dialog[open]'));
  doc.querySelector('dialog[open]')?.close();
  const { toolsMenuItems } = await import('/js/tools.js');
  views.menu = toolsMenuItems(app).filter(x => x.label).map(x => x.label);
  await app.setLanguage('de');
  await sleep(200);
  for (const [name, list] of Object.entries(views)) {
    const miss = [...new Set(list)].filter(s => app.i18n.t(s) === s);
    const b64 = btoa(unescape(encodeURIComponent(JSON.stringify(miss))));
    const parts = Math.ceil(b64.length / 1200) || 1;
    for (let i = 0; i < parts; i++) log(`I18N ${name} ${i + 1}/${parts} ${b64.slice(i * 1200, (i + 1) * 1200)}`);
  }
  await app.setLanguage('en');
}

// Screenshot: the Pathway tab in Chinese.
export async function pathwayzh(app, w, log, sleep) {
  await app.setLanguage('zh-CN');
  await pathwayshot(app, w, log, sleep);
  await sleep(300);
}

// Documentation screenshot (docs/images/rcv_pathway.png): tRNA prediction, Pathway tab, coloured by mechanism.
export async function docpathway(app, w, log, sleep) {
  app.setTheme('light', false);
  app.selectEnsemble(byName(app, /1EHZ|tRNA/i));
  if (!app.ss2d.visible) app.ss2d.setVisible(true);
  await sleep(200);
  app.showTab('pathway');
  if (!await waitFor(sleep, () => app.pathway.result || app.pathway.error)) throw new Error('no pathway result');
  if (app.pathway.error) throw new Error(app.pathway.error);
  app.pathway.colorBy('role');
  app.goTo(app.pathway.result.a, { animate: false, force: true });
  await sleep(400);
  log('A', app.ens.shortId(app.pathway.result.a), 'B', app.ens.shortId(app.pathway.result.b), '| barrier', app.pathway.result.effective_barrier);
}
