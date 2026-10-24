// Video export: deterministic frame plans (conformation tour / rotation / rock / transition A → B) rendered
// offscreen, with an optional overlay (title, conformation, energy along the order or along the transition
// path with a moving marker) and an optional 2D-structure inset that morphs with the 3D structure.
import { easeInOut, sleep, fmt, FIG_FONT, ensureFigureFont } from '../util.js';
import { ZipWriter } from './encoders.js';

export function videoFormats() {
  const out = [];
  if (typeof MediaRecorder !== 'undefined') {
    for (const [id, label, mime] of [
      ['mp4', 'MP4 (H.264)', 'video/mp4;codecs=avc1.640028'],
      ['mp4', 'MP4 (H.264)', 'video/mp4'],
      ['webm', 'WebM (VP9)', 'video/webm;codecs=vp9'],
      ['webm', 'WebM (VP8)', 'video/webm;codecs=vp8'],
    ]) if (!out.find(o => o.id === id) && MediaRecorder.isTypeSupported(mime)) out.push({ id, label, mime });
  }
  out.push({ id: 'zip', label: 'PNG frames (ZIP, lossless)', mime: 'application/zip' });
  return out;
}

// Lowest-cost route from a to b through the ensemble (Dijkstra on squared C4' RMSD, so many small steps are
// preferred over one large jump): the visited conformations act as intermediates of the transition movie.
export function transitionPath(ens, a, b) {
  if (a === b) return [a];
  const n = ens.n, dist = new Array(n).fill(Infinity), prev = new Array(n).fill(-1), done = new Array(n).fill(false);
  dist[a] = 0;
  for (;;) {
    let u = -1;
    for (let i = 0; i < n; i++) if (!done[i] && (u < 0 || dist[i] < dist[u])) u = i;
    if (u < 0 || dist[u] === Infinity || u === b) break;
    done[u] = true;
    for (let v = 0; v < n; v++) {
      if (done[v] || v === u) continue;
      const w = (ens.rmsd[u][v] || 0) ** 2 + 0.25;          // small constant: avoid needless zig-zags
      if (dist[u] + w < dist[v]) { dist[v] = dist[u] + w; prev[v] = u; }
    }
  }
  const path = [b];
  while (path[0] !== a && prev[path[0]] >= 0) path.unshift(prev[path[0]]);
  return path[0] === a ? path : [a, b];
}

// Transition A -> B: hold on A, morph (directly or through intermediates, leg time ~ RMSD), hold on B, optional return.
function transitionFrames(app, o) {
  const ens = app.ens, fps = o.fps;
  const a = o.from ?? app.state.prev ?? app.order[0], b = o.to ?? app.state.conf;
  const route = o.route === 'path' ? transitionPath(ens, a, b) : [a, b];
  const legs = [];
  for (let i = 0; i + 1 < route.length; i++) legs.push([route[i], route[i + 1], Math.max(0.3, ens.rmsd[route[i]][route[i + 1]] || 1)]);
  const total = legs.reduce((s, l) => s + l[2], 0) || 1;
  const frames = [];
  const hold = (k, sec, phase) => { for (let f = 0; f < Math.round(sec * fps); f++) frames.push({ from: k, to: k, t: 0, phase }); };
  const morph = (list, dir) => {
    const nTot = Math.max(2, Math.round(o.transition * fps));
    let acc = 0;
    for (const [x, y, w] of list) {
      const nLeg = Math.max(1, Math.round(nTot * w / total));
      for (let f = 1; f <= nLeg; f++) {
        const s = f / nLeg;
        frames.push({ from: x, to: y, t: easeInOut(s), phase: dir, progress: (acc + w * s) / total });
      }
      acc += w;
    }
  };
  hold(a, o.holdStart ?? o.hold ?? 1, 'start');
  if (legs.length) morph(legs, 'forward');
  hold(b, o.holdEnd ?? o.hold ?? 1, 'end');
  if (o.boomerang && legs.length) {
    morph(legs.slice().reverse().map(([x, y, w]) => [y, x, w]), 'back');
    hold(a, Math.min(1, o.holdStart ?? o.hold ?? 1), 'start');
  }
  frames.route = route;
  return frames;
}

// A plan is a list of frame states: { from, to, t, angle } (angle = rotation about view Y in radians, cumulative)
export function framePlan(app, o) {
  if (o.content === 'transition') {
    const frames = transitionFrames(app, o);
    const N = frames.length, spin = (o.spinDeg || 0) * Math.PI / 180;
    frames.forEach((fr, i) => { fr.angle = spin ? spin * (N > 1 ? i / (N - 1) : 0) - spin / 2 : 0; });
    return frames;
  }
  const fps = o.fps, frames = [];
  const order = app.order.length ? app.order : [app.state.conf];
  const tour = o.content === 'tour' || o.content === 'tour-rotate';
  if (tour && order.length > 1) {
    const seq = [...order];
    if (o.loop) seq.push(order[0]);
    const nHold = Math.round(o.hold * fps), nTrans = Math.max(1, Math.round(o.transition * fps));
    for (let i = 0; i < seq.length - 1; i++) {
      for (let f = 0; f < nHold; f++) frames.push({ from: seq[i], to: seq[i], t: 0 });
      for (let f = 1; f <= nTrans; f++) frames.push({ from: seq[i], to: seq[i + 1], t: easeInOut(f / nTrans) });
    }
    if (!o.loop) for (let f = 0; f < nHold; f++) frames.push({ from: seq[seq.length - 1], to: seq[seq.length - 1], t: 0 });
  } else {
    const n = Math.max(2, Math.round(o.seconds * fps));
    const k = app.state.conf;
    for (let f = 0; f < n; f++) frames.push({ from: k, to: k, t: 0 });
  }
  const N = frames.length;
  frames.forEach((fr, i) => {
    if (o.content === 'rotate') fr.angle = (2 * Math.PI * i) / N;
    else if (o.content === 'rock') fr.angle = (o.rockDeg * Math.PI / 180) * Math.sin((2 * Math.PI * i) / N);
    else if (o.content === 'tour-rotate') fr.angle = (2 * Math.PI * i) / N;
    else fr.angle = 0;
  });
  return frames;
}

export async function recordVideo(app, o, { onProgress = () => {}, signal } = {}) {
  await ensureFigureFont();
  const viewer = app.viewer;
  const frames = framePlan(app, o);
  const saved = app.snapshotView();
  const baseQ = viewer.controls.quaternion.clone();
  const cap = viewer.beginCapture(o.width, o.height);
  const out = document.createElement('canvas');
  out.width = o.width; out.height = o.height;
  const ctx = out.getContext('2d');
  const fmtDef = videoFormats().find(f => f.id === o.format) || videoFormats()[0];
  let recorder = null, chunks = [], track = null, zip = null;
  if (fmtDef.id === 'zip') zip = new ZipWriter();
  else {
    const stream = out.captureStream(0);
    track = stream.getVideoTracks()[0];
    recorder = new MediaRecorder(stream, { mimeType: fmtDef.mime, videoBitsPerSecond: o.bitrate });
    recorder.ondataavailable = e => { if (e.data && e.data.size) chunks.push(e.data); };
    recorder.start(250);
  }
  const t0 = performance.now();
  app.capturing = true;
  try {
    for (let i = 0; i < frames.length; i++) {
      if (signal?.aborted) throw new DOMException('cancelled', 'AbortError');
      const fr = frames[i];
      app.applyFrame(fr.from, fr.to, fr.t);
      viewer.controls.quaternion.copy(baseQ);
      if (fr.angle) {
        const q = new baseQ.constructor().setFromAxisAngle({ x: 0, y: 1, z: 0 }, fr.angle);
        viewer.controls.quaternion.copy(q.multiply(baseQ));
      }
      const r = cap.frame();
      ctx.drawImage(r.canvas, 0, 0);
      if (o.annotations !== false) app.annotations?.drawTo(ctx, r.project, o.width, o.height, Math.max(0.25, r.scale), app.display);
      if (o.overlay) (o.content === 'transition' ? drawTransitionOverlay : drawOverlay)(ctx, app, fr, o, frames);
      if (o.inset2d && app.ss2d) {
        const s = o.height / 1080, w = 420 * s, h = 360 * s;
        app.ss2d.drawTo(ctx, [o.width - w - 40 * s, 40 * s, w, h], { k: fr.to, from: fr.from, t: fr.t, dark: app.isDarkBackground(), unit: 1.3 * s });
      }
      if (zip) {
        const blob = await new Promise(res => out.toBlob(res, 'image/png'));
        zip.add(`frame_${String(i + 1).padStart(5, '0')}.png`, new Uint8Array(await blob.arrayBuffer()));
      } else {
        track.requestFrame();
        const due = t0 + ((i + 1) * 1000) / o.fps;
        const wait = due - performance.now();
        if (wait > 0) await sleep(wait);
      }
      onProgress((i + 1) / frames.length, i + 1, frames.length);
    }
  } finally {
    app.capturing = false;
    cap.end();
    viewer.controls.quaternion.copy(baseQ);
    app.restoreView(saved);
    if (app.surface?.enabled) app.surface.update();
  }
  if (zip) return { blob: zip.blob(), ext: 'zip', frames: frames.length };
  await new Promise(res => { recorder.onstop = res; recorder.stop(); });
  return { blob: new Blob(chunks, { type: fmtDef.mime.split(';')[0] }), ext: fmtDef.id, frames: frames.length };
}

// Overlay: title + conformation (top-left), legend (bottom-left), energy profile with moving marker (bottom-right).
export function drawOverlay(ctx, app, fr, o) {
  const W = ctx.canvas.width, H = ctx.canvas.height, s = H / 1080;
  const ens = app.ens;
  const dark = app.isDarkBackground();
  const ink = dark ? 'rgba(240,239,236,0.96)' : 'rgba(20,20,20,0.92)';
  const sub = dark ? 'rgba(200,198,190,0.9)' : 'rgba(82,81,78,0.95)';
  const font = (px, w = 500) => `${w >= 600 ? 'bold' : 'normal'} ${Math.round(px * s)}px ${FIG_FONT}`;
  ctx.save();
  ctx.textBaseline = 'alphabetic';
  ctx.fillStyle = ink;
  ctx.font = font(34, 650);
  ctx.fillText(ens.name, 48 * s, 70 * s);
  const k = fr.t > 0.5 ? fr.to : fr.from;
  const key = app.state.metric, m = ens.metricByKey[key] || { label: key };
  const va = ens.value(fr.from, key), vb = ens.value(fr.to, key);
  const v = va !== null && vb !== null ? va + (vb - va) * fr.t : ens.value(k, key);
  ctx.fillStyle = sub;
  ctx.font = font(24);
  const morphing = fr.from !== fr.to && fr.t > 0 && fr.t < 1;
  ctx.fillText(morphing ? `${ens.shortId(fr.from)} → ${ens.shortId(fr.to)}` : `${ens.shortId(k)} · ${app.order.indexOf(k) + 1} / ${ens.n}`, 48 * s, 108 * s);
  // energy mini-chart
  const cw = 420 * s, ch = 150 * s, x0 = W - cw - 48 * s, y0 = H - ch - 56 * s;
  const order = app.order;
  const vals = order.map(j => ens.value(j, key));
  const fin = vals.filter(x => x !== null);
  if (fin.length > 1) {
    let lo = Math.min(...fin), hi = Math.max(...fin);
    if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
    const X = i => x0 + (order.length > 1 ? i / (order.length - 1) : 0.5) * cw;
    const Y = val => y0 + ch - ((val - lo) / (hi - lo)) * ch;
    ctx.fillStyle = dark ? 'rgba(26,26,25,0.72)' : 'rgba(255,255,255,0.78)';
    roundRect(ctx, x0 - 18 * s, y0 - 52 * s, cw + 36 * s, ch + 84 * s, 12 * s);
    ctx.fill();
    ctx.fillStyle = sub;
    ctx.font = font(19, 600);
    ctx.fillText(`${m.label}${m.unit ? ' (' + m.unit + ')' : ''}`, x0, y0 - 24 * s);
    ctx.fillStyle = ink;
    ctx.font = font(26, 650);
    ctx.textAlign = 'right';
    ctx.fillText(fmt(v), x0 + cw, y0 - 22 * s);
    ctx.textAlign = 'left';
    ctx.strokeStyle = dark ? '#3987e5' : '#2a78d6';
    ctx.lineWidth = 3 * s;
    ctx.lineJoin = 'round';
    ctx.beginPath();
    vals.forEach((val, i) => { if (val === null) return; const px = X(i), py = Y(val); if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py); });
    ctx.stroke();
    ctx.fillStyle = ctx.strokeStyle;
    vals.forEach((val, i) => { if (val === null) return; ctx.beginPath(); ctx.arc(X(i), Y(val), 4.5 * s, 0, Math.PI * 2); ctx.fill(); });
    const ia = order.indexOf(fr.from), ib = order.indexOf(fr.to);
    if (ia >= 0 && ib >= 0 && va !== null && vb !== null) {
      const px = X(ia + (ib - ia) * fr.t), py = Y(v);
      ctx.beginPath(); ctx.arc(px, py, 10 * s, 0, Math.PI * 2);
      ctx.fillStyle = '#eb6834'; ctx.fill();
      ctx.lineWidth = 3 * s; ctx.strokeStyle = dark ? '#1a1a19' : '#ffffff'; ctx.stroke();
    }
  }
  // legend
  const lg = app.legendData();
  if (lg && o.legend !== false) drawLegend(ctx, lg, 48 * s, H - 56 * s, s, ink, sub);
  ctx.restore();
}

// Transition overlay: title, "A → B" with both conformations' values, and the metric along the route.
export function drawTransitionOverlay(ctx, app, fr, o, frames) {
  const W = ctx.canvas.width, H = ctx.canvas.height, s = H / 1080;
  const ens = app.ens;
  const dark = app.isDarkBackground();
  const ink = dark ? 'rgba(240,239,236,0.96)' : 'rgba(20,20,20,0.92)';
  const sub = dark ? 'rgba(200,198,190,0.9)' : 'rgba(82,81,78,0.95)';
  const font = (px, w = 500) => `${w >= 600 ? 'bold' : 'normal'} ${Math.round(px * s)}px ${FIG_FONT}`;
  const route = frames.route || [fr.from, fr.to];
  const a = route[0], b = route[route.length - 1];
  const key = app.state.metric, m = ens.metricByKey[key] || { label: key, unit: '' };
  const val = k => ens.value(k, key);
  const va = val(fr.from), vb = val(fr.to);
  const v = va !== null && vb !== null ? va + (vb - va) * fr.t : null;
  ctx.save();
  ctx.fillStyle = ink;
  ctx.font = font(34, 650);
  ctx.fillText(ens.name, 48 * s, 70 * s);
  ctx.font = font(24);
  ctx.fillStyle = sub;
  const fmtV = x => (x === null ? '' : ` (${fmt(x)}${m.unit ? ' ' + m.unit : ''})`);
  ctx.fillText(`A ${ens.shortId(a)}${fmtV(val(a))}  →  B ${ens.shortId(b)}${fmtV(val(b))}`, 48 * s, 108 * s);
  if (route.length > 2) { ctx.font = font(19); ctx.fillText(`via ${route.slice(1, -1).map(k => ens.shortId(k)).join(', ')}`, 48 * s, 138 * s); }
  // status chip
  const state = fr.phase === 'start' ? `A · ${ens.shortId(a)}` : fr.phase === 'end' ? `B · ${ens.shortId(b)}` : `${ens.shortId(fr.from)} → ${ens.shortId(fr.to)}`;
  ctx.font = font(21, 600);
  const tw = ctx.measureText(state).width + 28 * s;
  ctx.fillStyle = fr.phase === 'start' || fr.phase === 'end' ? (dark ? 'rgba(57,135,229,0.9)' : 'rgba(42,120,214,0.92)') : 'rgba(235,104,52,0.92)';
  roundRect(ctx, 48 * s, 160 * s, tw, 36 * s, 18 * s);
  ctx.fill();
  ctx.fillStyle = '#ffffff';
  ctx.fillText(state, 62 * s, 185 * s);
  // metric along the route (x = path progress)
  const vals = route.map(val);
  if (vals.every(x => x !== null) && route.length > 1) {
    const cw = 420 * s, ch = 130 * s, x0 = W - cw - 48 * s, y0 = H - ch - 56 * s;
    let lo = Math.min(...vals), hi = Math.max(...vals);
    if (hi - lo < 1e-9) { lo -= 1; hi += 1; }
    const w = [0];
    for (let i = 0; i + 1 < route.length; i++) w.push(w[i] + Math.max(0.3, ens.rmsd[route[i]][route[i + 1]] || 1));
    const X = i => x0 + (w[i] / w[w.length - 1]) * cw, Y = x => y0 + ch - ((x - lo) / (hi - lo)) * ch;
    ctx.fillStyle = dark ? 'rgba(26,26,25,0.72)' : 'rgba(255,255,255,0.8)';
    roundRect(ctx, x0 - 18 * s, y0 - 52 * s, cw + 36 * s, ch + 92 * s, 12 * s);
    ctx.fill();
    ctx.fillStyle = sub; ctx.font = font(19, 600);
    ctx.fillText(`${m.label}${m.unit ? ' (' + m.unit + ')' : ''} along the transition`, x0, y0 - 24 * s);
    ctx.strokeStyle = dark ? '#3987e5' : '#2a78d6'; ctx.lineWidth = 3 * s; ctx.lineJoin = 'round';
    ctx.beginPath();
    vals.forEach((x, i) => (i ? ctx.lineTo(X(i), Y(x)) : ctx.moveTo(X(i), Y(x))));
    ctx.stroke();
    ctx.fillStyle = ctx.strokeStyle;
    vals.forEach((x, i) => { ctx.beginPath(); ctx.arc(X(i), Y(x), 5 * s, 0, Math.PI * 2); ctx.fill(); });
    ctx.fillStyle = sub; ctx.font = font(17, 600); ctx.textAlign = 'center';
    ctx.fillText('A', X(0), y0 + ch + 26 * s); ctx.fillText('B', X(route.length - 1), y0 + ch + 26 * s);
    ctx.textAlign = 'left';
    const i0 = route.indexOf(fr.from), i1 = route.indexOf(fr.to);
    if (v !== null && i0 >= 0 && i1 >= 0) {
      const px = X(i0) + (X(i1) - X(i0)) * fr.t;
      ctx.beginPath(); ctx.arc(px, Y(v), 10 * s, 0, Math.PI * 2);
      ctx.fillStyle = '#eb6834'; ctx.fill();
      ctx.lineWidth = 3 * s; ctx.strokeStyle = dark ? '#1a1a19' : '#ffffff'; ctx.stroke();
    }
  }
  const lg = app.legendData();
  if (lg && o.legend !== false) drawLegend(ctx, lg, 48 * s, H - 56 * s, s, ink, sub);
  ctx.restore();
}

export function drawLegend(ctx, lg, x, yBottom, s, ink, sub) {
  ctx.save();
  const font = (px, w = 500) => `${w >= 600 ? 'bold' : 'normal'} ${Math.round(px * s)}px ${FIG_FONT}`;
  ctx.textBaseline = 'middle';
  if (lg.ramp) {
    const w = 260 * s, h = 16 * s, y = yBottom - h;
    const g = ctx.createLinearGradient(x, 0, x + w, 0);
    lg.ramp.forEach((c, i) => g.addColorStop(i / (lg.ramp.length - 1), c));
    ctx.fillStyle = g;
    ctx.fillRect(x, y, w, h);
    ctx.fillStyle = sub;
    ctx.font = font(18);
    ctx.fillText(lg.min, x, y - 16 * s);
    ctx.textAlign = 'right';
    ctx.fillText(lg.max, x + w, y - 16 * s);
    ctx.textAlign = 'left';
    ctx.font = font(18, 600);
    ctx.fillText(lg.title, x, y - 44 * s);
  } else {
    const items = lg.items || [];
    let cx = x, cy = yBottom - 12 * s;
    const sw = 22 * s;
    ctx.font = font(21);
    const rows = [];
    let row = [];
    let wsum = 0;
    for (const it of items) {
      const tw = ctx.measureText(it.label).width + sw + 34 * s;
      if (wsum + tw > 900 * s && row.length) { rows.push(row); row = []; wsum = 0; }
      row.push({ ...it, tw }); wsum += tw;
    }
    if (row.length) rows.push(row);
    cy = yBottom - 12 * s - (rows.length - 1) * 34 * s;
    ctx.fillStyle = sub;
    ctx.font = font(18, 600);
    ctx.fillText(lg.title, x, cy - 36 * s);
    for (const r of rows) {
      cx = x;
      for (const it of r) {
        ctx.fillStyle = it.color;
        roundRect(ctx, cx, cy - sw / 2, sw, sw, 4 * s);
        ctx.fill();
        ctx.fillStyle = ink;
        ctx.font = font(21);
        ctx.fillText(it.label, cx + sw + 9 * s, cy + 1);
        cx += it.tw;
      }
      cy += 34 * s;
    }
  }
  ctx.restore();
}

function roundRect(ctx, x, y, w, h, r) {
  ctx.beginPath();
  if (ctx.roundRect) ctx.roundRect(x, y, w, h, r); else ctx.rect(x, y, w, h);
}
