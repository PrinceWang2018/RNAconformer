// Secondary-structure drawing without DOM dependencies: layouts (ViennaRNA NAView from the document, or a
// radial-loop fallback), rigid 2D fitting, layout blending for morphs, and a drawing routine that renders
// into a canvas or SVG backend. Used by the 2D window, figure/video export and the interactive HTML export.
const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
const FONT = '"RCV Figure Sans", "Liberation Sans", Arimo, "Helvetica Neue", Helvetica, sans-serif';

// Nested pair table of a conformation (pseudoknotted pairs are drawn separately as arcs).
export function pairTable(n, pairs) {
  const pt = new Int32Array(n).fill(-1);
  for (const [i, j, , pk] of pairs) if (!pk && pt[i] < 0 && pt[j] < 0) { pt[i] = j; pt[j] = i; }
  return pt;
}
export function nestedKey(n, pairs) {
  const s = new Array(n).fill('.');
  const pt = pairTable(n, pairs);
  for (let i = 0; i < n; i++) if (pt[i] > i) { s[i] = '('; s[pt[i]] = ')'; }
  return s.join('');
}

// Radial loop layout (fallback when the document carries no ViennaRNA NAView layout): every loop is a regular
// polygon with unit edges, helices grow outwards from it, the exterior loop runs along a line.
export function simpleLayout(n, pt) {
  const P = new Float64Array(2 * n);
  const set = (i, x, y) => { P[2 * i] = x; P[2 * i + 1] = y; };
  const stem = (i, j, dx, dy) => {
    while (i + 1 < j - 1 && pt[i + 1] === j - 1) {
      set(i + 1, P[2 * i] + dx, P[2 * i + 1] + dy);
      set(j - 1, P[2 * j] + dx, P[2 * j + 1] + dy);
      i++; j--;
    }
    const verts = [i];
    for (let k = i + 1; k < j;) {
      if (pt[k] > k) { verts.push(k, pt[k]); k = pt[k] + 1; } else { verts.push(k); k++; }
    }
    verts.push(j);
    const m = verts.length;
    const R = 1 / (2 * Math.sin(Math.PI / m)), a = R * Math.cos(Math.PI / m);
    const mx = (P[2 * i] + P[2 * j]) / 2, my = (P[2 * i + 1] + P[2 * j + 1]) / 2;
    const cx = mx + dx * a, cy = my + dy * a;
    const t0 = Math.atan2(P[2 * i + 1] - cy, P[2 * i] - cx), tj = Math.atan2(P[2 * j + 1] - cy, P[2 * j] - cx);
    let d = tj - t0;
    while (d <= -Math.PI) d += 2 * Math.PI;
    while (d > Math.PI) d -= 2 * Math.PI;
    const step = -Math.sign(d || 1) * (2 * Math.PI / m);
    for (let q = 1; q < m - 1; q++) set(verts[q], cx + R * Math.cos(t0 + step * q), cy + R * Math.sin(t0 + step * q));
    for (let q = 1; q < m - 1; q++) {
      const a1 = verts[q];
      if (pt[a1] > a1 && verts[q + 1] === pt[a1]) {
        const b1 = pt[a1];
        const ox = (P[2 * a1] + P[2 * b1]) / 2 - cx, oy = (P[2 * a1 + 1] + P[2 * b1 + 1]) / 2 - cy;
        const l = Math.hypot(ox, oy) || 1;
        stem(a1, b1, ox / l, oy / l);
      }
    }
  };
  let x = 0;
  for (let i = 0; i < n;) {
    if (pt[i] > i) { set(i, x, 0); set(pt[i], x + 1, 0); stem(i, pt[i], 0, -1); x += 2; i = pt[i] + 1; } else { set(i, x, 0); x += 1; i++; }
  }
  // centre
  let sx = 0, sy = 0;
  for (let i = 0; i < n; i++) { sx += P[2 * i]; sy += P[2 * i + 1]; }
  for (let i = 0; i < n; i++) { P[2 * i] -= sx / n; P[2 * i + 1] -= sy / n; }
  return P;
}

export function layoutFor(ens, k) {
  const c = ens.conformations[k];
  ens._layoutCache = ens._layoutCache || new Map();
  const lay = c.layout_key && ens.doc.layouts2d?.[c.layout_key];
  if (lay && lay.length === 2 * ens.nRes) {
    const key = `L:${c.layout_key}`;
    if (!ens._layoutCache.has(key)) ens._layoutCache.set(key, relaxLayout(Float64Array.from(lay), c.pairs.filter(p => !p[3])));
    return ens._layoutCache.get(key);
  }
  const key = nestedKey(ens.nRes, c.pairs);
  if (!ens._layoutCache.has(key)) ens._layoutCache.set(key, relaxLayout(simpleLayout(ens.nRes, pairTable(ens.nRes, c.pairs)), c.pairs.filter(p => !p[3])));
  return ens._layoutCache.get(key);
}

// NAView (and any layout) can squeeze bulges and small loops so that nucleotides overlap. Backbone steps are brought
// back to at least one step, helices are kept rigid (pairs and their diagonals), and nucleotides closer than 0.9
// steps are pushed apart; drawings without overlaps are returned unchanged. Coordinates are in backbone-step units.
export function relaxLayout(P, pairs = []) {
  const n = P.length / 2;
  if (n < 3 || n > 5000) return P;
  const dist = (A, i, j) => Math.hypot(A[2 * i] - A[2 * j], A[2 * i + 1] - A[2 * j + 1]);
  const steps = [];
  for (let i = 0; i + 1 < n; i++) steps.push(dist(P, i, i + 1));
  const d0 = [...steps].sort((a, b) => a - b)[steps.length >> 1] || 1;
  const min = 0.9 * d0;
  const partner = new Int32Array(n).fill(-1);
  for (const [i, j] of pairs) { partner[i] = j; partner[j] = i; }
  const links = [], linked = new Set();
  const link = (a, b, L) => { if (a > b) [a, b] = [b, a]; if (a === b || linked.has(a * n + b)) return; linked.add(a * n + b); links.push([a, b, L]); };
  steps.forEach((L, i) => link(i, i + 1, Math.max(L, 0.95 * d0)));
  for (let i = 0; i < n; i++) {
    const j = partner[i];
    if (j <= i) continue;
    link(i, j, dist(P, i, j));
    if (i + 1 < n && partner[i + 1] === j - 1) { link(i, j - 1, dist(P, i, j - 1)); link(i + 1, j, dist(P, i + 1, j)); }   // stacked pairs: rigid
  }
  const overlaps = A => {
    const cell = new Map(), key = (x, y) => x * 1000003 + y, out = [];
    for (let i = 0; i < n; i++) { const k = key(Math.floor(A[2 * i] / min), Math.floor(A[2 * i + 1] / min)); if (!cell.has(k)) cell.set(k, []); cell.get(k).push(i); }
    for (let i = 0; i < n; i++) {
      const gx = Math.floor(A[2 * i] / min), gy = Math.floor(A[2 * i + 1] / min);
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (const j of cell.get(key(gx + dx, gy + dy)) || []) {
        if (j > i && !linked.has(i * n + j) && dist(A, i, j) < min) out.push([i, j]);
      }
    }
    return out;
  };
  const shortSteps = steps.some(L => L < 0.9 * d0);
  if (!shortSteps && !overlaps(P).length) return P;
  const Q = Float64Array.from(P);
  for (let it = 0; it < 200; it++) {
    for (let rep = 0; rep < 4; rep++) for (const [i, j, L0] of links) {
      const vx = Q[2 * j] - Q[2 * i], vy = Q[2 * j + 1] - Q[2 * i + 1], L = Math.hypot(vx, vy) || 1e-6, c = (L - L0) / L * 0.5;
      Q[2 * i] += vx * c; Q[2 * i + 1] += vy * c; Q[2 * j] -= vx * c; Q[2 * j + 1] -= vy * c;
    }
    const ov = overlaps(Q);
    let worst = 0;
    for (const [i, j] of ov) {
      const vx = Q[2 * j] - Q[2 * i], vy = Q[2 * j + 1] - Q[2 * i + 1], L = Math.hypot(vx, vy);
      worst = Math.max(worst, min - L);
      const ux = L > 1e-6 ? vx / L : Math.cos(i + j), uy = L > 1e-6 ? vy / L : Math.sin(i + j), push = (min - L) * 0.5;
      Q[2 * i] -= ux * push; Q[2 * i + 1] -= uy * push; Q[2 * j] += ux * push; Q[2 * j + 1] += uy * push;
    }
    let stretch = 0;
    for (const [i, j, L0] of links) stretch = Math.max(stretch, Math.abs(dist(Q, i, j) - L0));
    if (worst < 0.03 * d0 && stretch < 0.05 * d0) break;
  }
  return Q;
}

// Rigid 2D fit (rotation angle + centroids) of B onto A.
export function fit2D(A, B) {
  const n = A.length / 2;
  let ax = 0, ay = 0, bx = 0, by = 0;
  for (let i = 0; i < n; i++) { ax += A[2 * i]; ay += A[2 * i + 1]; bx += B[2 * i]; by += B[2 * i + 1]; }
  ax /= n; ay /= n; bx /= n; by /= n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    const px = B[2 * i] - bx, py = B[2 * i + 1] - by, qx = A[2 * i] - ax, qy = A[2 * i + 1] - ay;
    sxx += px * qx + py * qy; sxy += px * qy - py * qx;
  }
  return { phi: Math.atan2(sxy, sxx), ca: [ax, ay], cb: [bx, by] };
}

// Frame at t between layout A (shown) and B (target): interpolate in A's frame, rotating towards B's orientation.
export function blendLayouts(A, B, t, f = fit2D(A, B), out = new Float64Array(A.length)) {
  const n = A.length / 2;
  const c = Math.cos(f.phi), s = Math.sin(f.phi);
  const cr = Math.cos(-f.phi * t), sr = Math.sin(-f.phi * t);
  const ox = f.ca[0] + (f.cb[0] - f.ca[0]) * t, oy = f.ca[1] + (f.cb[1] - f.ca[1]) * t;
  for (let i = 0; i < n; i++) {
    const bx = B[2 * i] - f.cb[0], by = B[2 * i + 1] - f.cb[1];
    const fx = c * bx - s * by, fy = s * bx + c * by;                  // B aligned on A (relative to A's centroid)
    const x = (A[2 * i] - f.ca[0]) * (1 - t) + fx * t, y = (A[2 * i + 1] - f.ca[1]) * (1 - t) + fy * t;
    out[2 * i] = cr * x - sr * y + ox; out[2 * i + 1] = sr * x + cr * y + oy;
  }
  return out;
}

const pairKey = p => `${p[0]}-${p[1]}`;

// ------------------------------------------------------------------ drawing (shared by the window and exports)
// Draws into a backend with line/circle/path/text (canvas-like or SVG builder). coords in layout units.
// unit: size of one screen pixel in the target (1 on screen; larger for high-resolution exports), so that a figure
// looks the same in a small preview and in the exported file.
export function drawStructure2D(target, ens, coords, { k, prev = null, colors, t = 1, from = null, sel = new Set(), hover = new Set(), box, dark = false, labels = true, mono = false, unit = 1 }) {
  const n = ens.nRes;
  const [bx, by, bw, bh] = box;
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (let i = 0; i < n; i++) { minx = Math.min(minx, coords[2 * i]); maxx = Math.max(maxx, coords[2 * i]); miny = Math.min(miny, coords[2 * i + 1]); maxy = Math.max(maxy, coords[2 * i + 1]); }
  const pad = 1.4;
  const sc = Math.min(bw / (maxx - minx + 2 * pad), bh / (maxy - miny + 2 * pad));
  const ox = bx + (bw - (maxx - minx) * sc) / 2 - minx * sc, oy = by + (bh - (maxy - miny) * sc) / 2 - miny * sc;
  const X = i => ox + coords[2 * i] * sc, Y = i => oy + coords[2 * i + 1] * sc;
  const r = clamp(sc * 0.44, 1.2 * unit, 12 * unit);
  const ink = dark ? '#f0efec' : '#141414', sub = dark ? '#8f8d86' : '#898781', back = dark ? '#5c5b56' : '#b9b7ae';
  const gainCol = '#eb6834', pairCol = dark ? '#c3c2b7' : '#52514e', pkCol = '#6250d6';
  // backbone
  const segs = [];
  for (let i = 0; i + 1 < n; i++) segs.push(['L', X(i + 1), Y(i + 1)]);
  target.path([['M', X(0), Y(0)], ...segs], { stroke: back, width: Math.max(0.8 * unit, r * 0.28), fill: null });
  // pairs: those of the target conformation; while morphing, pairs of the source fade out
  const cur = ens.conformations[k].pairs;
  const ref = prev !== null && prev !== k ? new Set(ens.conformations[prev].pairs.map(pairKey)) : null;
  const drawPair = (p, alpha, gained) => {
    const [i, j, , pk] = p;
    if (pk) {
      const mx = (X(i) + X(j)) / 2, my = (Y(i) + Y(j)) / 2, dx = X(j) - X(i), dy = Y(j) - Y(i);
      const bend = 0.35;
      target.path([['M', X(i), Y(i)], ['C', mx - dy * bend, my + dx * bend, mx - dy * bend, my + dx * bend, X(j), Y(j)]],
        { stroke: pkCol, width: Math.max(0.8 * unit, r * 0.3), dash: [r * 0.6, r * 0.5], fill: null, alpha });
      return;
    }
    target.line(X(i), Y(i), X(j), Y(j), { stroke: gained ? gainCol : pairCol, width: Math.max(0.8 * unit, r * (gained ? 0.42 : 0.3)), alpha });
  };
  if (from !== null && from !== k && t < 1) {
    const curSet = new Set(cur.map(pairKey));
    // pairs that break fade out during the first half, new pairs fade in during the second half
    for (const p of ens.conformations[from].pairs) if (!curSet.has(pairKey(p)) && t < 0.5) drawPair(p, 1 - 2 * t, false);
    const fromSet = new Set(ens.conformations[from].pairs.map(pairKey));
    for (const p of cur) { const kept = fromSet.has(pairKey(p)); if (kept || t > 0.5) drawPair(p, kept ? 1 : 2 * t - 1, !kept); }
  } else for (const p of cur) drawPair(p, 1, !!ref && !ref.has(pairKey(p)));
  // nucleotides
  for (let i = 0; i < n; i++) {
    const hi = hover.has(i), s = sel.has(i);
    target.circle(X(i), Y(i), r * (hi ? 1.25 : 1), { fill: mono ? (dark ? '#2b2b29' : '#ffffff') : colors[i], stroke: s ? (dark ? '#ffffff' : '#111111') : (hi ? gainCol : (dark ? '#1a1a19' : '#ffffff')), width: s || hi ? Math.max(1.2 * unit, r * 0.32) : Math.max(0.5 * unit, r * 0.14), residue: i });
    if (labels && r >= 5.5 * unit) target.text(X(i), Y(i) + r * 0.37, ens.residues[i].name, { size: r * 1.05, fill: mono ? ink : textOn(colors[i]), anchor: 'middle', weight: 700, residue: i });
  }
  // numbering every 10 nt and the 5'/3' ends
  const tick = n > 150 ? 20 : 10;
  const cx = (minx + maxx) / 2, cy = (miny + maxy) / 2;
  const numLabel = (i, label) => {
    let dx = coords[2 * i] - cx, dy = coords[2 * i + 1] - cy;
    if (i > 0 && i < n - 1) {                // outward normal of the backbone at i
      const tx = coords[2 * i + 2] - coords[2 * i - 2], ty = coords[2 * i + 3] - coords[2 * i - 1];
      let nx = -ty, ny = tx;
      if (nx * dx + ny * dy < 0) { nx = -nx; ny = -ny; }
      dx = nx; dy = ny;
    }
    const l = Math.hypot(dx, dy) || 1;
    const fs = clamp(r * 1.1, 7.5 * unit, 12 * unit);
    target.text(X(i) + dx / l * (r + fs * 0.9), Y(i) + dy / l * (r + fs * 0.9) + fs * 0.35, label, { size: fs, fill: sub, anchor: 'middle', weight: 400 });
  };
  numLabel(0, "5′");
  if (n > 1) numLabel(n - 1, "3′");
  for (let i = tick - 1; i < n - 1; i += tick) numLabel(i, String(ens.residues[i].num ?? i + 1));
  return { X, Y, r, sc };
}

function textOn(hex) {
  const h = (hex || '#888888').replace('#', '');
  const v = parseInt(h.length === 3 ? h.split('').map(c => c + c).join('') : h, 16);
  const r = (v >> 16) & 255, g = (v >> 8) & 255, b = v & 255;
  return 0.299 * r + 0.587 * g + 0.114 * b > 150 ? '#141414' : '#ffffff';
}

// Backend that draws on a 2D canvas context.
export function canvas2DTarget(ctx) {
  const stroke = (o) => { ctx.strokeStyle = o.stroke; ctx.lineWidth = o.width || 1; ctx.setLineDash(o.dash || []); };
  return {
    line(x1, y1, x2, y2, o) { ctx.save(); ctx.globalAlpha = o.alpha ?? 1; stroke(o); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); ctx.restore(); },
    path(cmds, o) {
      ctx.save(); ctx.globalAlpha = o.alpha ?? 1; ctx.lineJoin = 'round'; ctx.lineCap = 'round'; ctx.beginPath();
      for (const c of cmds) { if (c[0] === 'M') ctx.moveTo(c[1], c[2]); else if (c[0] === 'L') ctx.lineTo(c[1], c[2]); else if (c[0] === 'C') ctx.bezierCurveTo(c[1], c[2], c[3], c[4], c[5], c[6]); }
      if (o.fill) { ctx.fillStyle = o.fill; ctx.fill(); }
      if (o.stroke) { stroke(o); ctx.stroke(); }
      ctx.restore();
    },
    circle(x, y, r, o) { ctx.save(); ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); if (o.fill) { ctx.fillStyle = o.fill; ctx.fill(); } if (o.stroke) { ctx.setLineDash([]); ctx.strokeStyle = o.stroke; ctx.lineWidth = o.width || 1; ctx.stroke(); } ctx.restore(); },
    text(x, y, s, o) { ctx.save(); ctx.fillStyle = o.fill; ctx.font = `${(o.weight || 400) >= 600 ? 'bold' : 'normal'} ${o.size}px ${FONT}`; ctx.textAlign = o.anchor === 'middle' ? 'center' : o.anchor || 'left'; ctx.fillText(s, x, y); ctx.restore(); },
  };
}



// Tall drawings (long hairpins) are turned sideways in landscape boxes; decided once per ensemble so that morphs
// between conformations never flip the orientation.
export function ensembleAspect(ens) {
  if (ens._ss2dAspect) return ens._ss2dAspect;
  let a = 0;
  for (let k = 0; k < ens.n; k++) {
    const P = layoutFor(ens, k);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < P.length; i += 2) { x0 = Math.min(x0, P[i]); x1 = Math.max(x1, P[i]); y0 = Math.min(y0, P[i + 1]); y1 = Math.max(y1, P[i + 1]); }
    a = Math.max(a, (y1 - y0 + 2) / (x1 - x0 + 2));
  }
  return (ens._ss2dAspect = a);
}
export function orientForBox(ens, coords, w, h) {
  if (!(ensembleAspect(ens) > 1.4 && w / h > 1.15)) return coords;
  const out = new Float64Array(coords.length);
  for (let i = 0; i < coords.length; i += 2) { out[i] = -coords[i + 1]; out[i + 1] = coords[i]; }   // 5'/3' ends to the left
  return out;
}
