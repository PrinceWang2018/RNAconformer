// Z-RNA regions (0.7.6): high-confidence left-handed double-helical segments, found again on every displayed frame
// (still conformations, morphs, scrubbing, video, pulling) and drawn with their own theme in 3D, in 2D and in exports.
// Same rules as rnafold3d/zrna.py (tests/test_zrna.py runs both on the same structures):
//   * at least MIN_PAIRS consecutive stacked Watson-Crick pairs (C1'-C1' 8.5-12.5 Å in the shown frame),
//   * every two-step window left-handed (sum of the two twists <= WINDOW_MAX), no step above STEP_MAX,
//   * the purine of nearly every pair syn and the pyrimidines anti.
// Default on; View ▸ Z-RNA regions turns it off.
import { loadPref, savePref } from './util.js';

export const Z = {
  MIN_PAIRS: 3, WINDOW_MAX: -30, STEP_MAX: 25, SYN_FRACTION: 0.75, ANTI_FRACTION: 0.9,
  SYN_LO: -20, SYN_HI: 120, C1_MIN: 8.5, C1_MAX: 12.5,
};
export const Z_COLOR = '#c2185b';
const RING = ['C2', 'C4', 'C6'];

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = a => Math.hypot(a[0], a[1], a[2]);

function dihedral(p0, p1, p2, p3) {
  const b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
  const l = norm(b1) || 1, u = [b1[0] / l, b1[1] / l, b1[2] / l];
  const v = sub(b0, u.map(x => x * dot(b0, u))), w = sub(b2, u.map(x => x * dot(b2, u)));
  return Math.atan2(dot(cross(u, v), w), dot(v, w)) * 180 / Math.PI;
}

// atoms: function(r, name) -> [x,y,z] or null
export function chiOf(atoms, r, purine) {
  const [g, ref] = purine ? ['N9', 'C4'] : ['N1', 'C2'];
  const p = ["O4'", "C1'", g, ref].map(n => atoms(r, n));
  return p.every(Boolean) ? dihedral(...p) : null;
}
const isSyn = chi => chi !== null && chi >= Z.SYN_LO && chi <= Z.SYN_HI;

function normal(atoms, r) {
  const p = RING.map(n => atoms(r, n));
  if (!p.every(Boolean)) return null;
  const n = cross(sub(p[1], p[0]), sub(p[2], p[0])), l = norm(n);
  return l > 1e-6 ? n.map(x => x / l) : null;
}

export function stepTwist(atoms, i, j, k, l) {
  const c = [i, j, k, l].map(r => atoms(r, "C1'"));
  if (!c.every(Boolean)) return null;
  const ns = [i, j, k, l].map(r => normal(atoms, r));
  if (!ns.every(Boolean)) return null;
  const down = [0, 1, 2].map(d => (c[2][d] - c[0][d]) + (c[3][d] - c[1][d]));
  const h = [0, 0, 0];
  for (const n of ns) { const s = dot(n, down) >= 0 ? 1 : -1; for (let d = 0; d < 3; d++) h[d] += s * n[d]; }
  const lh = norm(h);
  if (lh < 1e-6) return null;
  for (let d = 0; d < 3; d++) h[d] /= lh;
  let a1 = sub(c[1], c[0]), a2 = sub(c[3], c[2]);
  a1 = sub(a1, h.map(x => x * dot(a1, h))); a2 = sub(a2, h.map(x => x * dot(a2, h)));
  if (norm(a1) < 1e-6 || norm(a2) < 1e-6) return null;
  return Math.atan2(dot(cross(a1, a2), h), dot(a1, a2)) * 180 / Math.PI;
}

export function helixRuns(pairs) {
  const ps = pairs.map(([i, j]) => [Math.min(i, j), Math.max(i, j)]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const key = (i, j) => i * 100003 + j;
  const have = new Set(ps.map(([i, j]) => key(i, j))), used = new Set(), runs = [];
  for (const p of ps) {
    if (used.has(key(...p)) || have.has(key(p[0] - 1, p[1] + 1))) continue;
    const run = [p]; used.add(key(...p));
    for (;;) {
      const [a, b] = run[run.length - 1];
      if (!(have.has(key(a + 1, b - 1)) && b - a > 2)) break;
      run.push([a + 1, b - 1]); used.add(key(a + 1, b - 1));
    }
    runs.push(run);
  }
  return runs;
}

// seq: parent letters; atoms(r, name); pairs [[i, j], ...] 0-based. Returns regions with 0-based residues.
export function detectZ(seq, atoms, pairs, { checkDistance = true } = {}) {
  const chiCache = new Map();
  const chi = r => { if (!chiCache.has(r)) chiCache.set(r, chiOf(atoms, r, 'AG'.includes(seq[r]))); return chiCache.get(r); };
  const ok = pairs.filter(([i, j]) => {
    if (!checkDistance) return true;
    const a = atoms(i, "C1'"), b = atoms(j, "C1'");
    if (!a || !b) return false;
    const d = norm(sub(a, b));
    return d >= Z.C1_MIN && d <= Z.C1_MAX;
  });
  const out = [];
  for (const run of helixRuns(ok)) {
    if (run.length < Z.MIN_PAIRS) continue;
    const tw = [];
    for (let s = 0; s + 1 < run.length; s++) tw.push(stepTwist(atoms, run[s][0], run[s][1], run[s + 1][0], run[s + 1][1]));
    const zp = run.map(([i, j]) => {
      const pur = [i, j].filter(r => 'AG'.includes(seq[r])), pyr = [i, j].filter(r => 'CU'.includes(seq[r]));
      return [pur.length === 1 && isSyn(chi(pur[0])), pyr.every(r => chi(r) !== null && !isSyn(chi(r))), pyr.length];
    });
    let a = 0;
    while (a < run.length) {
      let b = a;
      while (b + 1 < run.length) {
        const t = tw[b];
        if (t === null || t > Z.STEP_MAX) break;
        if (b > a && (tw[b - 1] === null || tw[b - 1] + t > Z.WINDOW_MAX)) break;
        b++;
      }
      if (b - a + 1 >= Z.MIN_PAIRS) {
        let seg = []; for (let s = a; s <= b; s++) seg.push(s);
        while (seg.length && !zp[seg[0]][0]) seg = seg.slice(1);
        while (seg.length && !zp[seg[seg.length - 1]][0]) seg = seg.slice(0, -1);
        if (seg.length >= Z.MIN_PAIRS) {
          const steps = seg.slice(0, -1).map(s => tw[s]);
          const synF = seg.filter(s => zp[s][0]).length / seg.length;
          const nPyr = seg.reduce((s, x) => s + zp[x][2], 0);
          const antiF = nPyr ? seg.reduce((s, x) => s + (zp[x][1] ? zp[x][2] : 0), 0) / nPyr : 0;
          let windows = steps.length >= 2;
          for (let s = 0; s + 1 < steps.length; s++) if (steps[s] + steps[s + 1] > Z.WINDOW_MAX) windows = false;
          if (synF >= Z.SYN_FRACTION && antiF >= Z.ANTI_FRACTION && windows) {
            const prs = seg.map(s => run[s]);
            out.push({ pairs: prs, residues: [...new Set(prs.flat())].sort((x, y) => x - y),
              meanTwist: steps.reduce((s, x) => s + x, 0) / steps.length, twists: steps, synPurines: synF, antiPyrimidines: antiF });
          }
        }
      }
      a = b + 1;
    }
  }
  return out;
}

export function ensembleAtoms(ens, coords) {
  return (r, name) => {
    const a = ens.atomIndex[r]?.[name];
    return a === undefined ? null : [coords[3 * a], coords[3 * a + 1], coords[3 * a + 2]];
  };
}

function parentSeq(ens) {
  const s = ens.sequence || '';
  if (s.length === ens.nRes) return s.toUpperCase();
  return ens.residues.map((res, r) => ('N9' in ens.atomIndex[r] ? 'G' : 'C')).join('');
}

// ------------------------------------------------------------------ app controller
export class ZRegions {
  constructor(app) {
    this.app = app;
    this.enabled = loadPref('zrna', { enabled: true }).enabled !== false;
    this.regions = [];
    this.residues = new Set();
    this._morph = null;
    this._key = null;
    this.svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    this.svg.setAttribute('class', 'zrna-layer');
    app.viewer.container.append(this.svg);
    app.viewer.renderHooks.add((proj, w, h) => this._draw(proj, w, h));
    app.on('morph', m => { this._morph = m.done ? null : m; });
    app.on('ensemble', () => { this._morph = null; this._key = null; this.update(); });
  }

  set(on) {
    this.enabled = !!on;
    savePref('zrna', { enabled: this.enabled });
    this._key = null;
    this.update();
    this.app.viewer.needsRender = true;
    this.app.ss2d?.draw?.();
    this.app.rebuildMenus?.();
  }

  // pairs to test: those of the shown conformation, and of both ends while morphing or pulling
  _pairs() {
    const ens = this.app.ens, ks = new Set([this.app.state.conf]);
    if (this._morph && this._morph.from !== undefined) { ks.add(this._morph.from); ks.add(this._morph.to); }
    const seen = new Set(), out = [];
    for (const k of ks) for (const p of ens.conformations[k]?.pairs || []) {
      const key = `${p[0]},${p[1]}`;
      if (!seen.has(key)) { seen.add(key); out.push([p[0], p[1]]); }
    }
    return out;
  }

  detect(coords = this.app.display, k = null) {
    const ens = this.app.ens;
    if (!ens || !coords) return [];
    const pairs = k === null ? this._pairs() : ens.conformations[k].pairs.map(p => [p[0], p[1]]);
    return detectZ(parentSeq(ens), ensembleAtoms(ens, coords), pairs);
  }

  update() {
    const ens = this.app.ens, coords = this.app.display;
    if (!this.enabled || !ens || !coords) { this.regions = []; this.residues = new Set(); return; }
    // cheap change key: a few coordinates and the conformation(s)
    const n = coords.length, key = `${this.app.state.conf}|${this._morph ? `${this._morph.from}>${this._morph.to}` : ''}|${coords[0]},${coords[n >> 1]},${coords[n - 1]},${coords[(n * 3) >> 2]}`;
    if (key === this._key) return;
    this._key = key;
    const prev = this.residues;
    this.regions = this.detect(coords);
    this.residues = new Set(this.regions.flatMap(z => z.residues));
    if (prev.size !== this.residues.size || [...prev].some(r => !this.residues.has(r))) {
      this.app.emit('zrna', this.regions);
      this.app.ss2d?.draw?.();
    }
  }

  // ---------------------------------------------------------------- drawing (live SVG layer and exports)
  _draw(proj, w, h) {
    this.update();
    const root = this.svg;
    root.setAttribute('viewBox', `0 0 ${w} ${h}`);
    root.replaceChildren();
    if (!this.regions.length) return;
    this.paint(svgPainter(root), proj, 1, this.app.display);
  }

  drawTo(ctx, project, w, h, s, coords = this.app.display, k = null) {
    if (!this.enabled || !this.app.ens) return;
    const regions = coords === this.app.display && k === null ? this.regions : this.detect(coords, k);
    this.paint(canvasPainter(ctx), project, s, coords, regions);
  }

  paint(T, proj, s, coords, regions = this.regions) {
    const ens = this.app.ens;
    for (const z of regions) {
      const pts = z.residues.map(r => proj(...ens.residueCenter(coords, r))).filter(p => p[2] >= -1 && p[2] <= 1);
      if (!pts.length) continue;
      for (const p of pts) T.circle(p[0], p[1], 12 * s, { fill: 'rgba(194,24,91,0.24)' });
      // the helix axis: midpoints of the pairs, joined
      const mids = z.pairs.map(([i, j]) => {
        const a = ens.residueCenter(coords, i), b = ens.residueCenter(coords, j);
        return proj((a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2);
      });
      for (let q = 0; q + 1 < mids.length; q++) T.line(mids[q][0], mids[q][1], mids[q + 1][0], mids[q + 1][1], { stroke: Z_COLOR, width: 3 * s, dash: [6 * s, 4 * s] });
      const top = pts.reduce((m, p) => (p[1] < m[1] ? p : m), pts[0]);
      const cx = pts.reduce((a, p) => a + p[0], 0) / pts.length;
      const label = `Z-RNA ${ens.resLabel(z.residues[0])}–${ens.resLabel(z.residues[z.residues.length - 1])}`;
      T.badge(cx, top[1] - 18 * s, label, s);
    }
  }

  summary(regions = this.regions) {
    const ens = this.app.ens;
    return regions.map(z => ({
      residues: z.residues.map(r => ens.resLabel(r)), pairs: z.pairs.map(([i, j]) => `${ens.resLabel(i)}-${ens.resLabel(j)}`),
      meanTwist: Math.round(z.meanTwist * 10) / 10, synPurines: z.synPurines,
    }));
  }
}

const FONT = '600 12px system-ui, sans-serif';
function svgPainter(root) {
  const NS = 'http://www.w3.org/2000/svg';
  const mk = (tag, attrs) => { const e = document.createElementNS(NS, tag); for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) e.setAttribute(k, v); root.append(e); return e; };
  return {
    circle(cx, cy, r, o) { mk('circle', { cx, cy, r, fill: o.fill, 'pointer-events': 'none' }); },
    line(x1, y1, x2, y2, o) { mk('line', { x1, y1, x2, y2, stroke: o.stroke, 'stroke-width': o.width, 'stroke-dasharray': o.dash ? o.dash.join(' ') : null, 'stroke-linecap': 'round' }); },
    badge(x, y, text, s) {
      const w = text.length * 7 * s + 14 * s, h = 20 * s;
      mk('rect', { x: x - w / 2, y: y - h / 2, width: w, height: h, rx: 10 * s, fill: Z_COLOR, opacity: 0.92 });
      const t = mk('text', { x, y: y + 4 * s, 'text-anchor': 'middle', fill: '#ffffff', style: `font: ${FONT}` });
      t.textContent = text;
    },
  };
}
function canvasPainter(ctx) {
  return {
    circle(cx, cy, r, o) { ctx.save(); ctx.fillStyle = o.fill; ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2); ctx.fill(); ctx.restore(); },
    line(x1, y1, x2, y2, o) { ctx.save(); ctx.strokeStyle = o.stroke; ctx.lineWidth = o.width; ctx.setLineDash(o.dash || []); ctx.lineCap = 'round'; ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); ctx.restore(); },
    badge(x, y, text, s) {
      ctx.save(); ctx.font = `600 ${12 * s}px system-ui, sans-serif`;
      const w = ctx.measureText(text).width + 14 * s, h = 20 * s;
      ctx.globalAlpha = 0.92; ctx.fillStyle = Z_COLOR; ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(x - w / 2, y - h / 2, w, h, 10 * s); else ctx.rect(x - w / 2, y - h / 2, w, h);
      ctx.fill(); ctx.globalAlpha = 1; ctx.fillStyle = '#ffffff'; ctx.textAlign = 'center'; ctx.fillText(text, x, y + 4 * s); ctx.restore();
    },
  };
}
