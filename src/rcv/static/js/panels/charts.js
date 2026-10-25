// Shared chart helpers (SVG) and a tiny SVG walker used to re-draw charts into PDF and canvas exports.
import { svg, niceTicks, fmt, FIG_FONT } from '../util.js';

export function linear(d0, d1, r0, r1) {
  const span = d1 - d0 || 1;
  const f = v => r0 + (v - d0) / span * (r1 - r0);
  f.invert = p => d0 + (p - r0) / (r1 - r0) * span;
  f.domain = [d0, d1];
  return f;
}

export function extent(values, pad = 0.08) {
  const v = values.filter(x => Number.isFinite(x));
  if (!v.length) return [0, 1];
  let lo = Math.min(...v), hi = Math.max(...v);
  if (lo === hi) { const p = Math.abs(lo) * 0.1 || 1; lo -= p; hi += p; }
  const p = (hi - lo) * pad;
  return [lo - p, hi + p];
}

// Y axis in publication style: axis line, outward ticks, hairline grid (optional), rotated title; returns the ticks.
export function yAxis(g, y, x0, x1, { ticks = 4, unit = '', title = '', grid = true } = {}) {
  const t = niceTicks(y.domain[0], y.domain[1], ticks).filter(v => v >= Math.min(...y.domain) - 1e-9 && v <= Math.max(...y.domain) + 1e-9);
  const gridG = svg('g', { class: 'grid' });
  const ax = svg('g', { class: 'axis' });
  const lab = svg('g', { class: 'ticklab' });
  const yTop = y(y.domain[1]), yBot = y(y.domain[0]);
  const step = t.length > 1 ? Math.abs(t[1] - t[0]) : 1;
  const stepDigits = step > 0 && step < 1 ? Math.ceil(-Math.log10(step) - 1e-9) : 0;     // 0.05 apart -> 2 decimals
  ax.append(svg('line', { x1: x0 + 0.5, x2: x0 + 0.5, y1: Math.min(yTop, yBot), y2: Math.max(yTop, yBot) }));
  for (const v of t) {
    const yy = Math.round(y(v)) + 0.5;
    if (grid) gridG.append(svg('line', { x1: x0, x2: x1, y1: yy, y2: yy }));
    ax.append(svg('line', { x1: x0 - 4, x2: x0, y1: yy, y2: yy }));
    lab.append(svg('text', { x: x0 - 6, y: y(v) + 3.5, 'text-anchor': 'end', text: fmt(v, Math.max(Math.abs(v) < 10 ? 1 : 0, stepDigits)) }));
  }
  g.append(gridG, ax, lab);
  const ttl = title || unit;
  if (ttl) {
    const yy = (yTop + yBot) / 2, xx = 10;
    g.append(svg('text', { class: 'atitle', x: xx, y: yy, 'text-anchor': 'middle', transform: `rotate(-90 ${xx} ${yy})`, text: ttl }));
  }
  return t;
}

// X axis line with outward ticks at the given positions; labels [[x, text]]; optional title.
export function xAxis(g, x0, x1, yb, labels, { title = '', titleY = null } = {}) {
  const ax = svg('g', { class: 'axis' });
  ax.append(svg('line', { x1: x0, x2: x1, y1: yb + 0.5, y2: yb + 0.5 }));
  const lab = svg('g', { class: 'ticklab' });
  for (const [xx, text] of labels) {
    ax.append(svg('line', { x1: Math.round(xx) + 0.5, x2: Math.round(xx) + 0.5, y1: yb, y2: yb + 4 }));
    lab.append(svg('text', { x: xx, y: yb + 15, 'text-anchor': 'middle', text }));
  }
  g.append(ax, lab);
  if (title) g.append(svg('text', { class: 'atitle', x: (x0 + x1) / 2, y: titleY ?? yb + 29, 'text-anchor': 'middle', text: title }));
}

// --------------------------------------------------------------------------- export walker
// Walks a rendered <svg> and emits primitives with resolved (computed) colours:
//   backend.line(x1,y1,x2,y2,style) / path(cmds,style) / circle(cx,cy,r,style) / rect(x,y,w,h,rx,style) / text(x,y,str,style)
// Supported subset (what the app's charts use): line, polyline, path (M L C Z, absolute), circle, rect, text, g[transform=translate].
export function walkSVG(root, backend, ox = 0, oy = 0) {
  const visit = (node, tx, ty) => {
    if (node.nodeType !== 1) return;
    const cs = getComputedStyle(node);
    if (cs.display === 'none' || cs.visibility === 'hidden') return;
    let dx = tx, dy = ty;
    const tr = node.getAttribute('transform');
    if (tr) {
      const m = /translate\(\s*([-\d.e]+)[ ,]*([-\d.e]*)\s*\)/.exec(tr);
      if (m) { dx += parseFloat(m[1]); dy += parseFloat(m[2] || '0'); }
    }
    const tag = node.tagName.toLowerCase();
    const style = {
      fill: cs.fill, stroke: cs.stroke, strokeWidth: parseFloat(cs.strokeWidth) || 1,
      opacity: parseFloat(cs.opacity) * parseFloat(cs.fillOpacity ?? 1), strokeOpacity: parseFloat(cs.opacity) * parseFloat(cs.strokeOpacity ?? 1),
      dash: cs.strokeDasharray && cs.strokeDasharray !== 'none' ? cs.strokeDasharray.split(/[ ,]+/).map(parseFloat).filter(Number.isFinite) : null,
      fontSize: parseFloat(cs.fontSize) || 10, fontWeight: cs.fontWeight, anchor: node.getAttribute('text-anchor') || cs.textAnchor || 'start',
    };
    const n = a => parseFloat(node.getAttribute(a) || '0');
    if (tag === 'line') backend.line(n('x1') + dx, n('y1') + dy, n('x2') + dx, n('y2') + dy, style);
    else if (tag === 'circle') backend.circle(n('cx') + dx, n('cy') + dy, n('r'), style);
    else if (tag === 'rect') { if (!node.classList.contains('hit')) backend.rect(n('x') + dx, n('y') + dy, n('width'), n('height'), n('rx'), style); }
    else if (tag === 'polyline') {
      const pts = (node.getAttribute('points') || '').trim().split(/\s+/).map(p => p.split(',').map(Number));
      if (pts.length > 1) backend.path([['M', pts[0][0] + dx, pts[0][1] + dy], ...pts.slice(1).map(p => ['L', p[0] + dx, p[1] + dy])], style);
    } else if (tag === 'path') {
      const d = node.getAttribute('d') || '';
      const toks = d.match(/[MLCZ]|-?\d*\.?\d+(?:e[-+]?\d+)?/gi) || [];
      const cmds = [];
      let i = 0;
      while (i < toks.length) {
        const c = toks[i++].toUpperCase();
        if (c === 'Z') { cmds.push(['Z']); continue; }
        const k = c === 'C' ? 6 : 2;
        const v = toks.slice(i, i + k).map(Number);
        i += k;
        cmds.push([c, ...v.map((x, j) => x + (j % 2 === 0 ? dx : dy))]);
      }
      if (cmds.length) backend.path(cmds, style);
    } else if (tag === 'text') {
      // rotated labels: transform="rotate(angle x y)" (axis titles); the angle is passed to the backend
      const rot = /rotate\(\s*([-\d.e]+)/.exec(tr || '');
      if (node.textContent) backend.text(n('x') + dx, n('y') + dy, node.textContent, { ...style, angle: rot ? parseFloat(rot[1]) : 0 });
    }
    if (tag === 'g' || tag === 'svg') for (const ch of node.children) visit(ch, dx, dy);
  };
  for (const ch of root.children) visit(ch, ox, oy);
}

export function parseColor(c) {
  if (!c || c === 'none' || c === 'transparent') return null;
  const m = /rgba?\(([^)]+)\)/.exec(c);
  if (!m) return null;
  const p = m[1].split(/[ ,/]+/).filter(Boolean).map(Number);
  if (p.length === 4 && p[3] === 0) return null;
  return { r: p[0] / 255, g: p[1] / 255, b: p[2] / 255, a: p.length === 4 ? p[3] : 1 };
}

// Canvas backend for the walker (used to place charts into PNG/TIFF figures and video frames).
export function canvasBackend(ctx, scale = 1) {
  const col = (c, a = 1) => { const p = parseColor(c); return p ? `rgba(${Math.round(p.r * 255)},${Math.round(p.g * 255)},${Math.round(p.b * 255)},${p.a * a})` : null; };
  const stroke = s => {
    const c = col(s.stroke, s.strokeOpacity);
    if (!c) return false;
    ctx.strokeStyle = c; ctx.lineWidth = s.strokeWidth * scale; ctx.setLineDash((s.dash || []).map(v => v * scale));
    ctx.lineCap = 'round'; ctx.lineJoin = 'round';
    return true;
  };
  return {
    line(x1, y1, x2, y2, s) { if (stroke(s)) { ctx.beginPath(); ctx.moveTo(x1 * scale, y1 * scale); ctx.lineTo(x2 * scale, y2 * scale); ctx.stroke(); } },
    circle(cx, cy, r, s) {
      ctx.beginPath(); ctx.arc(cx * scale, cy * scale, r * scale, 0, Math.PI * 2);
      const f = col(s.fill, s.opacity); if (f) { ctx.fillStyle = f; ctx.fill(); }
      if (stroke(s)) ctx.stroke();
    },
    rect(x, y, w, h, rx, s) {
      ctx.beginPath();
      if (ctx.roundRect && rx) ctx.roundRect(x * scale, y * scale, w * scale, h * scale, rx * scale); else ctx.rect(x * scale, y * scale, w * scale, h * scale);
      const f = col(s.fill, s.opacity); if (f) { ctx.fillStyle = f; ctx.fill(); }
      if (stroke(s)) ctx.stroke();
    },
    path(cmds, s) {
      ctx.beginPath();
      for (const c of cmds) {
        if (c[0] === 'M') ctx.moveTo(c[1] * scale, c[2] * scale);
        else if (c[0] === 'L') ctx.lineTo(c[1] * scale, c[2] * scale);
        else if (c[0] === 'C') ctx.bezierCurveTo(c[1] * scale, c[2] * scale, c[3] * scale, c[4] * scale, c[5] * scale, c[6] * scale);
        else if (c[0] === 'Z') ctx.closePath();
      }
      const f = col(s.fill, s.opacity); if (f) { ctx.fillStyle = f; ctx.fill(); }
      if (stroke(s)) ctx.stroke();
    },
    text(x, y, str, s) {
      const f = col(s.fill, s.opacity); if (!f) return;
      ctx.fillStyle = f;
      ctx.font = `${s.fontWeight >= 600 ? 'bold' : 'normal'} ${s.fontSize * scale}px ${FIG_FONT}`;
      ctx.textAlign = s.anchor === 'middle' ? 'center' : s.anchor === 'end' ? 'right' : 'left';
      ctx.textBaseline = 'alphabetic';
      ctx.setLineDash([]);
      if (s.angle) { ctx.save(); ctx.translate(x * scale, y * scale); ctx.rotate(s.angle * Math.PI / 180); ctx.fillText(str, 0, 0); ctx.restore(); }
      else ctx.fillText(str, x * scale, y * scale);
    },
  };
}
