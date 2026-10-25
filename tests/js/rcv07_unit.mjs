// Unit tests for RCV 0.7 browser modules (run by tests/test_rcv07.py when Node >= 18 is available):
//   node tests/js/rcv07_unit.mjs <case.json>   (case: { doc, result } - an ensemble document and a Pathway result)
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
let failures = 0;
const ok = (cond, msg) => { if (!cond) { failures++; console.error('FAIL', msg); } else console.log('ok', msg); };

// minimal DOM for the SVG builders
class FakeNode {
  constructor(tag) { this.tag = tag; this.attrs = {}; this.children = []; this.textContent = ''; }
  setAttribute(k, v) { this.attrs[k] = String(v); }
  append(...c) { this.children.push(...c); }
  addEventListener() {}
}
globalThis.Node = FakeNode;
globalThis.document = { createElementNS: (ns, tag) => new FakeNode(tag), createElement: tag => new FakeNode(tag), createTextNode: t => { const n = new FakeNode('#text'); n.textContent = t; return n; } };
const count = (n, tag) => (n.tag === tag ? 1 : 0) + n.children.reduce((s, c) => s + count(c, tag), 0);

const { buildMechanism, chainStrain, rotationFromCovariance, quatPow, quatRotate, smooth } = await import(path.join(STATIC, 'mechmorph.js'));
const { drawMap, glyph, qualityText } = await import(path.join(STATIC, 'landscape.js'));
const { Ensemble } = await import(path.join(STATIC, 'ensemble.js'));

const { doc, result: r } = JSON.parse(readFileSync(process.argv[2], 'utf8'));
const ens = new Ensemble(doc);

// ---- rotations
const rng = (s => () => (s = (s * 16807) % 2147483647) / 2147483647)(42);
let worst = 0;
for (let k = 0; k < 20; k++) {
  const ax = [rng() - 0.5, rng() - 0.5, rng() - 0.5], n = Math.hypot(...ax), u = ax.map(x => x / n), ang = 3 * rng();
  const c = Math.cos(ang), s = Math.sin(ang);
  const R = [c + u[0] * u[0] * (1 - c), u[0] * u[1] * (1 - c) - u[2] * s, u[0] * u[2] * (1 - c) + u[1] * s,
    u[1] * u[0] * (1 - c) + u[2] * s, c + u[1] * u[1] * (1 - c), u[1] * u[2] * (1 - c) - u[0] * s,
    u[2] * u[0] * (1 - c) - u[1] * s, u[2] * u[1] * (1 - c) + u[0] * s, c + u[2] * u[2] * (1 - c)];
  const H = Array(9).fill(0);
  for (let p = 0; p < 6; p++) {
    const v = [rng() * 4 - 2, rng() * 4 - 2, rng() * 4 - 2];
    const q = [R[0] * v[0] + R[1] * v[1] + R[2] * v[2], R[3] * v[0] + R[4] * v[1] + R[5] * v[2], R[6] * v[0] + R[7] * v[1] + R[8] * v[2]];
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) H[3 * a + b] += v[a] * q[b];
  }
  const Rr = rotationFromCovariance(H);
  worst = Math.max(worst, ...Rr.map((x, i) => Math.abs(x - R[i])));
}
ok(worst < 1e-9, `rotation from covariance (20 random rotations, max error ${worst.toExponential(1)})`);
const q = [Math.cos(0.6), 0, Math.sin(0.6), 0], h = quatPow(q, 0.5);
const v1 = quatRotate(q, [1, 2, 3]), v2 = quatRotate(h, quatRotate(h, [1, 2, 3]));
ok(Math.hypot(v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2]) < 1e-12, 'half rotation applied twice equals the rotation');
ok(smooth(0) === 0 && smooth(1) === 1 && Math.abs(smooth(0.5) - 0.5) < 1e-12, 'smooth step');

// ---- mechanism morph on the Pathway result
const A = ens.coords[r.a], B = ens.coords[r.b];
const trace = [...ens.trace];
for (const staged of [true, false]) {
  const t0 = Date.now();
  const m = buildMechanism(A, B, ens.atomRes, r.motion, r.timing, { staged, names: ens.atomNames });
  const frames = [], lin = [];
  for (let k = 0; k <= 24; k++) {
    const t = k / 24;
    frames.push(m.at(t));
    const L = new Float32Array(A.length);
    for (let i = 0; i < A.length; i++) L[i] = A[i] + (B[i] - A[i]) * t;
    lin.push(L);
  }
  let e0 = 0, e1 = 0;
  for (let i = 0; i < A.length; i++) { e0 = Math.max(e0, Math.abs(frames[0][i] - A[i])); e1 = Math.max(e1, Math.abs(frames[24][i] - B[i])); }
  ok(e0 < 1e-4 && e1 < 1e-4, `${staged ? 'staged' : 'rigid'} morph reproduces A and B exactly`);
  const sm = chainStrain(frames, trace), sl = chainStrain(lin, trace);
  ok(sm < sl, `${staged ? 'staged' : 'rigid'} morph keeps the chain better than linear interpolation (strain ${sm.toFixed(2)} vs ${sl.toFixed(2)}, ${Date.now() - t0} ms)`);
}

// ---- map drawing and 2D glyph
const L = doc.ensemble.landscape;
const g = new FakeNode('g');
const pos = drawMap(g, ens, L, { x0: 0, y0: 0, x1: 300, y1: 220 }, { current: 0, routes: [{ nodes: [r.a, r.b] }], marks: [{ k: r.a, text: 'A' }] });
ok(count(g, 'circle') >= ens.n + L.classes.length, `map draws ${ens.n} points and ${L.classes.length} disk(s)`);
ok([...Array(ens.n).keys()].every(k => pos.px(k) >= 0 && pos.px(k) <= 300 && pos.py(k) >= 0 && pos.py(k) <= 220), 'every point inside the box');
const key = ens.conformations[0].layout_key;
const gl = glyph(doc.layouts2d[key], key, 50, 50, 30);
ok(count(gl, 'path') === 1 && count(gl, 'line') === (key.match(/\(/g) || []).length, 'glyph: backbone path and one line per base pair (flat layout)');
ok(/trustworthiness/.test(qualityText(L)), 'map quality caption');

if (failures) { console.error(`${failures} failure(s)`); process.exit(1); }
console.log('all 0.7 browser tests passed');
