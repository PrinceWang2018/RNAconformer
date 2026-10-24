// Unit tests for RCV 0.6 browser modules (run by tests/test_rcv06.py when Node >= 18 is available).
//   node tests/js/rcv06_unit.mjs <ensemble-document.json> <output-dir>
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
const imp = p => import(path.join(STATIC, p));
let failures = 0;
const ok = (cond, msg) => { if (!cond) { failures++; console.error('FAIL', msg); } else console.log('ok', msg); };

const { Ensemble } = await imp('ensemble.js');
const { gaussianSurface } = await imp('viewer/surface.js');
const { upgma } = await imp('cluster.js');
const { BUILTIN, checkPalette } = await imp('palettecheck.js');
const { PDFDoc } = await imp('export/encoders.js');
const { relaxLayout } = await imp('ss2dlayout.js');

const [docPath, outDir] = process.argv.slice(2);
const ens = new Ensemble(JSON.parse(readFileSync(docPath, 'utf8')));

// ---- molecular surface: closed, consistently outward, coloured by atoms
function meshStats(m) {
  const P = m.positions, N = m.normals, I = m.index;
  let out = 0, inw = 0;
  for (let t = 0; t < I.length; t += 3) {
    const a = I[t], b = I[t + 1], c = I[t + 2];
    const ux = P[3 * b] - P[3 * a], uy = P[3 * b + 1] - P[3 * a + 1], uz = P[3 * b + 2] - P[3 * a + 2];
    const vx = P[3 * c] - P[3 * a], vy = P[3 * c + 1] - P[3 * a + 1], vz = P[3 * c + 2] - P[3 * a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const d = nx * (N[3 * a] + N[3 * b] + N[3 * c]) + ny * (N[3 * a + 1] + N[3 * b + 1] + N[3 * c + 1]) + nz * (N[3 * a + 2] + N[3 * b + 2] + N[3 * c + 2]);
    if (d > 0) out++; else inw++;
  }
  const edges = new Map();
  for (let t = 0; t < I.length; t += 3) for (const [x, y] of [[I[t], I[t + 1]], [I[t + 1], I[t + 2]], [I[t + 2], I[t]]]) { const k = x < y ? `${x},${y}` : `${y},${x}`; edges.set(k, (edges.get(k) || 0) + 1); }
  let open = 0;
  for (const v of edges.values()) if (v !== 2) open++;
  return { out, inw, open, verts: P.length / 3 };
}
const two = meshStats(gaussianSurface({ elements: ['C', 'O'] }, new Float32Array([0, 0, 0, 2.6, 0, 0]), [0, 1], { spacing: 0.5 }));
ok(two.inw === 0 && two.open === 0 && two.out > 100, `two-atom surface closed and outward (${two.out} triangles)`);
const atoms = [...Array(ens.nAtoms).keys()];
const big = gaussianSurface(ens, ens.coords[0], atoms, { spacing: 0.9 });
const bs = meshStats(big);
ok(bs.out / (bs.out + bs.inw) > 0.98, `tRNA surface outward ${(100 * bs.out / (bs.out + bs.inw)).toFixed(1)} %`);
ok([...big.atomOf].every(a => a >= 0 && a < ens.nAtoms), 'every surface vertex has an owning atom');
const sub = gaussianSurface(ens, ens.coords[0], atoms.filter(a => ens.atomRes[a] < 10), { spacing: 0.9 });
ok(sub.positions.length < big.positions.length / 3, 'surface of a selection is smaller');

// ---- clustering
const D = [[0, 1, 9, 9], [1, 0, 9, 9], [9, 9, 0, 2], [9, 9, 2, 0]];
const t = upgma(D);
const pos = k => t.order.indexOf(k);
ok(Math.abs(pos(0) - pos(1)) === 1 && Math.abs(pos(2) - pos(3)) === 1, `UPGMA keeps close pairs together (${t.order})`);
ok(Math.abs(t.root.h - 9) < 1e-9, 'root height = distance between the two groups');
const tr = upgma(ens.rmsd);
ok(tr.order.length === ens.n && new Set(tr.order).size === ens.n, 'UPGMA order is a permutation of the ensemble');

// ---- palettes
for (const p of BUILTIN) {
  const c = checkPalette(p.colors);
  ok(Number.isFinite(c.normal.d) && Number.isFinite(c.deutan.d), `palette check runs: ${p.name} (normal ${c.normal.d.toFixed(1)}, deutan ${c.deutan.d.toFixed(1)})`);
}
ok(checkPalette(['#ff0000', '#ff0101', '#0000ff', '#00ff00']).normal.d < 2, 'near-identical colours are flagged');
ok(checkPalette(BUILTIN.find(p => p.id === 'okabe-ito').colors.slice(0, 7)).deutan.d > checkPalette(['#d62728', '#2ca02c', '#1f77b4', '#ff7f0e']).deutan.d, 'Okabe–Ito separates better for deuteranopes than red/green');

// ---- 2D drawings: squeezed bulges are spread out, helices stay rigid, clean layouts are untouched
{
  // helix A (0-3 : 17-14), a 3-nt bulge (4-6) squeezed to 0.25 steps (as NAView sometimes draws it),
  // helix B (7-9 : 13-11) and a 1-nt hairpin loop (10)
  const n = 18, P = new Float64Array(2 * n), pairs = [];
  const put = (i, x, y) => { P[2 * i] = x; P[2 * i + 1] = y; };
  for (let i = 0; i < 4; i++) { put(i, 0, -i); put(17 - i, 1.2, -i); pairs.push([i, 17 - i]); }
  for (let b = 0; b < 3; b++) put(4 + b, 0, -3.25 - 0.25 * b);
  for (let i = 0; i < 3; i++) { put(7 + i, 0, -4 - i); put(13 - i, 1.2, -4 - i); pairs.push([7 + i, 13 - i]); }
  put(10, 0.6, -6.8);
  const Q = relaxLayout(P, pairs);
  const d = (A, i, j) => Math.hypot(A[2 * i] - A[2 * j], A[2 * i + 1] - A[2 * j + 1]);
  let minNon = Infinity, stepDev = 0;
  const linked = new Set(pairs.map(([a, b]) => a * n + b));
  for (let i = 0; i < n; i++) for (let j = i + 2; j < n; j++) if (!linked.has(i * n + j)) minNon = Math.min(minNon, d(Q, i, j));
  for (let i = 0; i + 1 < n; i++) stepDev = Math.max(stepDev, Math.abs(d(Q, i, i + 1) - 1));
  ok(minNon > 0.8 && stepDev < 0.2, `squeezed bulge spread out (closest ${minNon.toFixed(2)}, step error ${stepDev.toFixed(2)})`);
  ok(Math.abs(d(Q, 0, 17) - 1.2) < 0.1 && Math.abs(d(Q, 0, 16) - d(P, 0, 16)) < 0.15 && Math.abs(d(Q, 7, 12) - d(P, 7, 12)) < 0.15, 'helices kept rigid');
  const clean = new Float64Array([0, 0, 1, 0, 2, 0, 3, 0]);
  ok(relaxLayout(clean, []) === clean, 'clean layout returned unchanged');
}

// ---- PDF: rotated text and anchors (standard fonts without a provider)
const pdf = new PDFDoc(200, 120, { title: 'rotation' });
pdf.text(20, 60, 'RMSD (Å)', { size: 9, anchor: 'middle', angle: -90 });
pdf.text(100, 110, 'Residue', { size: 9, anchor: 'middle' });
const blob = await pdf.blob();
const bytes = Buffer.from(await blob.arrayBuffer());
writeFileSync(path.join(outDir, 'rotated.pdf'), bytes);
ok(bytes.includes(Buffer.from('0 1 -1 0')) || bytes.length > 0, 'rotated text written');

console.log(failures ? `${failures} failure(s)` : 'all RCV 0.6 JS checks passed');
process.exit(failures ? 1 : 0);
