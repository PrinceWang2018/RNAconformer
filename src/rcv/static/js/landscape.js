// Topology-aware ensemble map (computed by rcv/landscape.py): every secondary-structure topology has its own disk,
// disks are arranged by 3D similarity, and conformations inside are placed by C4' RMSD. Drawn in the Energy tab
// (Map), the Pathway tab (with routes), the plot studio and the interactive HTML page. Only SVG primitives the
// chart exporter understands are used (circle, line, path, text, g with translate).
import { svg } from './util.js';

export function landscapeOf(ens) {
  return ens?.doc?.ensemble?.landscape || ens?._landscape || null;
}

// Fetch the map from the engine for documents saved before it existed (cached on the ensemble).
export async function loadLandscape(app, ens) {
  if (!ens) return null;
  const have = landscapeOf(ens);
  if (have) return have;
  if (!app?.api?.online) return null;
  if (!ens._landscapeLoading) {
    const body = ens.serverId ? { document: ens.serverId } : { doc: ens.doc };
    ens._landscapeLoading = app.api.post('/api/landscape', body).then(L => { ens._landscape = L; return L; }).catch(() => null);
  }
  return ens._landscapeLoading;
}

// Small 2D drawing (backbone and pair lines) of a secondary structure layout, fitted into a size x size box.
// layout: the document's 2D drawing, a flat list x0, y0, x1, y1, ... in screen orientation (y down), or [[x, y], ...]
export function glyph(layout, dbn, cx, cy, size, { stroke = 'var(--text-2)', pairStroke = 'var(--series-1)' } = {}) {
  const g = svg('g', { class: 'ls-glyph' });
  if (!layout || layout.length < 2) return g;
  if (typeof layout[0] === 'number') { const pts = []; for (let i = 0; i + 1 < layout.length; i += 2) pts.push([layout[i], layout[i + 1]]); layout = pts; }
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  for (const [x, y] of layout) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
  const s = size / Math.max(x1 - x0, y1 - y0, 1e-6);
  const X = x => cx + (x - (x0 + x1) / 2) * s, Y = y => cy + (y - (y0 + y1) / 2) * s;
  const d = layout.map(([x, y], i) => `${i ? 'L' : 'M'} ${X(x).toFixed(1)} ${Y(y).toFixed(1)}`).join(' ');
  g.append(svg('path', { d, fill: 'none', stroke, 'stroke-width': 0.9, 'stroke-linejoin': 'round' }));
  if (dbn) {
    const st = [];
    for (let i = 0; i < dbn.length && i < layout.length; i++) {
      if (dbn[i] === '(') st.push(i);
      else if (dbn[i] === ')' && st.length) {
        const j = st.pop();
        g.append(svg('line', { x1: X(layout[j][0]), y1: Y(layout[j][1]), x2: X(layout[i][0]), y2: Y(layout[i][1]), stroke: pairStroke, 'stroke-width': 0.7 }));
      }
    }
  }
  return g;
}

function arrowHead(x0, y0, x1, y1, size, color) {
  const a = Math.atan2(y1 - y0, x1 - x0);
  const p1 = [x1 - size * Math.cos(a - 0.45), y1 - size * Math.sin(a - 0.45)];
  const p2 = [x1 - size * Math.cos(a + 0.45), y1 - size * Math.sin(a + 0.45)];
  return svg('path', { d: `M ${x1.toFixed(1)} ${y1.toFixed(1)} L ${p1[0].toFixed(1)} ${p1[1].toFixed(1)} L ${p2[0].toFixed(1)} ${p2[1].toFixed(1)} Z`, fill: color, stroke: 'none' });
}

// Draw the map into SVG group g inside box {x0, y0, x1, y1}. Returns the screen position of every conformation.
// opts: mode 'map' | 'mds', color(k), radius(k), current, prev, routes [{nodes:[k], width, color, dash}],
//       glyphs, labels, edges, scaleBar, fontSize
export const SCREEN_INK = { diskFill: 'var(--surface-2)', diskStroke: 'var(--border-strong)', text: 'var(--text)', text2: 'var(--text-2)', text3: 'var(--text-3)',
  axis: 'var(--axis)', surface: 'var(--surface)', point: 'var(--series-1)', route: 'var(--series-2)', glyph: 'var(--text-2)', glyphPair: 'var(--series-1)' };
export const PRINT_INK = { diskFill: '#f3f2ee', diskStroke: '#b9b7b0', text: '#1a1a1a', text2: '#4d4d4d', text3: '#6b6b6b',
  axis: '#9a9a9a', surface: '#ffffff', point: '#2a78d6', route: '#eb6834', glyph: '#4d4d4d', glyphPair: '#2a78d6' };

export function drawMap(g, ens, L, box, opts = {}) {
  const C = { ...SCREEN_INK, ...(opts.ink || {}) };
  const mode = opts.mode || 'map';
  const pos = mode === 'mds' ? (L.mds || ens.embedding) : L.positions;
  const n = pos.length;
  const fs = opts.fontSize || 10;
  const showDisks = mode === 'map' && opts.disks !== false;
  // extent: points (and disks in map mode) with room for the labels
  let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
  const grow = (x, y, r = 0) => { x0 = Math.min(x0, x - r); x1 = Math.max(x1, x + r); y0 = Math.min(y0, y - r); y1 = Math.max(y1, y + r); };
  pos.forEach(([x, y]) => grow(x, y));
  if (showDisks) for (const c of L.classes) grow(c.center[0], c.center[1], c.radius);
  if (!Number.isFinite(x0)) { x0 = -1; x1 = 1; y0 = -1; y1 = 1; }
  const labelRoom = showDisks && opts.labels !== false ? fs * 2.6 : 0;
  const W = box.x1 - box.x0, H = box.y1 - box.y0 - labelRoom;
  const sc = Math.min(W / Math.max(x1 - x0, 1e-6), H / Math.max(y1 - y0, 1e-6)) * 0.94;
  const cx = (box.x0 + box.x1) / 2, cy = box.y0 + labelRoom + H / 2;
  const mx = (x0 + x1) / 2, my = (y0 + y1) / 2;
  const px = k => cx + (pos[k][0] - mx) * sc, py = k => cy - (pos[k][1] - my) * sc;
  const PX = x => cx + (x - mx) * sc, PY = y => cy - (y - my) * sc;
  const layer = svg('g', { class: 'ls-map' });
  if (showDisks) {
    // links between related topologies (one helix apart, or the nearest)
    if (opts.edges !== false) {
      for (const e of L.edges || []) {
        const A = L.classes[e.a], B = L.classes[e.b];
        const ax = PX(A.center[0]), ay = PY(A.center[1]), bx = PX(B.center[0]), by = PY(B.center[1]);
        const d = Math.hypot(bx - ax, by - ay) || 1;
        const ra = A.radius * sc, rb = B.radius * sc;
        if (d <= ra + rb) continue;
        const ux = (bx - ax) / d, uy = (by - ay) / d;
        layer.append(svg('line', { x1: ax + ux * ra, y1: ay + uy * ra, x2: bx - ux * rb, y2: by - uy * rb, stroke: C.axis, 'stroke-width': 1, 'stroke-dasharray': e.kind === 'one helix' ? '3 3' : '1 3' }));
      }
    }
    for (const c of L.classes) {
      const r = Math.max(6, c.radius * sc);
      layer.append(svg('circle', { cx: PX(c.center[0]), cy: PY(c.center[1]), r, fill: C.diskFill, stroke: C.diskStroke, 'stroke-width': 1, class: 'ls-disk' }));
    }
    if (opts.labels !== false) {
      for (const c of L.classes) {
        const r = Math.max(6, c.radius * sc);
        const lx = PX(c.center[0]), ly = PY(c.center[1]) - r - 4;
        const shape = c.shape.length > 14 ? c.shape.slice(0, 13) + '…' : c.shape;
        const txt = `${c.label}  ${shape}  n=${c.members.length}`;
        layer.append(svg('text', { x: lx, y: ly, 'text-anchor': 'middle', 'font-size': fs, fill: C.text2, text: txt }));
        if (opts.glyphs !== false && ens?.doc?.layouts2d) {
          const rep = c.members.find(k => ens.conformations?.[k]?.layout_key && ens.doc.layouts2d[ens.conformations[k].layout_key]);
          if (rep !== undefined) {
            const key = ens.conformations[rep].layout_key;
            const size = Math.min(fs * 2.4, Math.max(14, r * 0.9));
            const gx = lx + (txt.length * fs * 0.27) + size / 2 + 4;
            layer.append(glyph(ens.doc.layouts2d[key], key, gx, ly - size / 2 + fs * 0.3, size, { stroke: C.glyph, pairStroke: C.glyphPair }));
          }
        }
      }
    }
  }
  // routes (Pathway): polylines with arrowheads through the conformations of each route
  for (const rt of opts.routes || []) {
    const nodes = (rt.nodes || []).filter(k => k !== null && k !== undefined && k >= 0 && k < n);
    if (nodes.length < 2) continue;
    const col = rt.color || C.route;
    const w = rt.width || 2;
    for (let i = 0; i + 1 < nodes.length; i++) {
      const a = nodes[i], b = nodes[i + 1];
      const xa = px(a), ya = py(a), xb = px(b), yb = py(b);
      const d = Math.hypot(xb - xa, yb - ya);
      if (d < 1) continue;
      const shorten = Math.min(7, d / 3);
      const ex = xb - (xb - xa) / d * shorten, ey = yb - (yb - ya) / d * shorten;
      layer.append(svg('line', { x1: xa, y1: ya, x2: ex, y2: ey, stroke: col, 'stroke-width': w, 'stroke-dasharray': rt.dash || null, 'stroke-linecap': 'round', opacity: rt.opacity ?? 1 }));
      layer.append(arrowHead(xa, ya, ex, ey, 4 + w * 1.4, col));
    }
  }
  // points
  const color = opts.color || (() => C.point);
  const radius = opts.radius || (() => 4);
  for (let k = 0; k < n; k++) layer.append(svg('circle', { cx: px(k), cy: py(k), r: radius(k), fill: color(k), stroke: C.surface, 'stroke-width': opts.pointStroke ?? 1.5 }));
  if (opts.prev !== null && opts.prev !== undefined && opts.prev !== opts.current && opts.prev < n) {
    layer.append(svg('circle', { cx: px(opts.prev), cy: py(opts.prev), r: 8, fill: 'none', stroke: C.text3, 'stroke-width': 1.5, 'stroke-dasharray': '2.5 2' }));
  }
  if (opts.current !== null && opts.current !== undefined && opts.current < n) {
    layer.append(svg('circle', { cx: px(opts.current), cy: py(opts.current), r: opts.ringR ?? 8, fill: 'none', stroke: C.text, 'stroke-width': opts.ringW ?? 2 }));
  }
  for (const m of opts.marks || []) {           // labelled conformations (A, B)
    if (m.k === null || m.k === undefined || m.k >= n) continue;
    layer.append(svg('text', { x: px(m.k) + fs * 0.8, y: py(m.k) - fs * 0.6, 'font-size': fs + 1, 'font-weight': 700, fill: m.color || C.text, text: m.text }));
  }
  // scale bar
  if (opts.scaleBar !== false) {
    const target = (W / 7) / sc;               // about a seventh of the width, rounded to 1, 2 or 5 x 10^k
    const len = [0.5, 1, 2, 5, 10, 20, 50, 100].reduce((best, l) => (Math.abs(Math.log(l / target)) < Math.abs(Math.log(best / target)) ? l : best), 1);
    const bx1 = box.x1 - 2, bx0 = bx1 - len * sc, by = box.y1 - 2;
    if (bx0 > box.x0) {
      layer.append(svg('line', { x1: bx0, x2: bx1, y1: by, y2: by, stroke: C.text3, 'stroke-width': 1.5 }));
      layer.append(svg('text', { x: (bx0 + bx1) / 2, y: by - 4, 'text-anchor': 'middle', 'font-size': fs, fill: C.text3, text: `${len} Å` }));
    }
  }
  g.append(layer);
  return { px, py, scale: sc };
}

// Map caption with its quality numbers (this layout vs plain MDS).
export function qualityText(L) {
  const q = L?.quality;
  if (!q?.map) return '';
  const f = v => (v === null || v === undefined ? '–' : v.toFixed(2));
  return `Neighbourhoods kept: trustworthiness ${f(q.map.trustworthiness)} (MDS ${f(q.mds?.trustworthiness)}), topology purity ${f(q.map.topology_purity)} (MDS ${f(q.mds?.topology_purity)}), stress ${f(q.map.stress)} (MDS ${f(q.mds?.stress)})`;
}
