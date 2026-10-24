// Unit tests for RCV's DOM-free browser modules (run by tests/test_app.py when Node >= 18 is available).
//   node tests/js/app_unit.mjs <ensemble-document.json> <output-dir>
// Writes encoder outputs into <output-dir> for validation by Python (Pillow, zipfile) and exits non-zero on failure.
import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
const imp = p => import(path.join(STATIC, p));
const [docPath, outDir] = process.argv.slice(2);
let failures = 0;
const ok = (cond, msg) => { if (!cond) { failures++; console.error('FAIL', msg); } else console.log('ok', msg); };

const { Ensemble, kabsch, applyTransform, toPDB } = await imp('ensemble.js');
const { residueColors, legend, SCHEMES } = await imp('colors.js');
const { encodeTIFF, pngWithDpi, ZipWriter, PDFDoc, crc32, toWinAnsi } = await imp('export/encoders.js');
const { framePlan } = await imp('export/video.js');
const { ElasticPull, rankConformations, chooseStable, fittedCoords } = await imp('pull.js');

// ---- ensemble model
const doc = JSON.parse(readFileSync(docPath, 'utf8'));
const ens = new Ensemble(doc);
ok(ens.n === doc.conformations.length && ens.nRes === doc.sequence.length, 'ensemble sizes');
ok(ens.coords.every(c => c.length === ens.nAtoms * 3), 'decoded coordinates match the atom table');
ok([...ens.trace].every(a => a >= 0), "every residue has a trace atom (C4')");
ok(ens.rings.filter(Boolean).length === ens.nRes, 'every residue has base-ring atoms');
for (const mode of ['rank', 'metric', 'path', 'input']) {
  const o = ens.order(mode, doc.default_metric);
  ok(o.length === ens.n && new Set(o).size === ens.n, `order "${mode}" is a permutation`);
}
// superposition: rotate a copy, then recover it
const A = ens.coords[0], B = Float32Array.from(A);
const th = 0.7, R = [[Math.cos(th), -Math.sin(th), 0], [Math.sin(th), Math.cos(th), 0], [0, 0, 1]];
applyTransform(B, R, [3, -2, 5]);
const idx = [...ens.trace];
const fit = kabsch(B, A, idx);
applyTransform(B, fit.R, fit.t);
let dev = 0;
for (let i = 0; i < A.length; i++) dev = Math.max(dev, Math.abs(A[i] - B[i]));
ok(dev < 1e-3, `quaternion superposition recovers a rigid transform (max dev ${dev.toExponential(2)})`);
const pdb = toPDB(ens, A, { model: 1 });
ok(pdb.split('\n').filter(l => l.startsWith('ATOM')).length === ens.nAtoms, 'PDB writer emits every atom');
ok(pdb.split('\n').every(l => !l.startsWith('ATOM') || l.length === 78), 'PDB ATOM records are fixed-width');

// ---- colour schemes
for (const s of SCHEMES) {
  const c = residueColors(ens, s.id, { conf: 0, coords: ens.coords[0], prevCoords: ens.coords[ens.n - 1], uniform: '#5b7a99' });
  ok(c.length === ens.nRes && c.every(h => /^#[0-9a-f]{6}$/i.test(h)), `scheme "${s.id}" gives a colour per residue`);
  ok(!!legend(ens, s.id, { conf: 0 }), `scheme "${s.id}" has a legend`);
}

// ---- video frame plans
const fakeApp = { order: ens.order('rank'), state: { conf: 0 }, ens };
const plan = framePlan(fakeApp, { content: 'tour', fps: 10, transition: 0.5, hold: 0.2, loop: true, seconds: 3, rockDeg: 30 });
ok(plan.length === ens.n * (2 + 5), `tour plan length (${plan.length})`);
ok(plan.every(f => f.t >= 0 && f.t <= 1), 'tour plan interpolation within [0, 1]');
const rot = framePlan(fakeApp, { content: 'rotate', fps: 10, seconds: 3 });
ok(rot.length === 30 && Math.abs(rot[29].angle - 2 * Math.PI * 29 / 30) < 1e-9, 'rotation plan covers 360°');

// ---- pull to switch
{
  const r = Math.floor(ens.nRes / 2);
  const p = new ElasticPull(ens, ens.coords[0], r, { pairs: ens.conformations[0].pairs });
  const a0 = [p.X0[3 * r], p.X0[3 * r + 1], p.X0[3 * r + 2]];
  const T = [a0[0] + 8, a0[1], a0[2]];
  for (let s = 1; s <= 10; s++) { p.setTarget(a0.map((v, i) => v + (T[i] - v) * s / 10)); p.solve(30); }
  ok(Math.hypot(p.X[3 * r] - T[0], p.X[3 * r + 1] - T[1], p.X[3 * r + 2] - T[2]) < 1e-6, 'pulled residue follows the cursor');
  let worst = 0;
  for (let i = 0; i + 1 < ens.nRes; i++) {
    const d0 = Math.hypot(...[0, 1, 2].map(c => p.X0[3 * i + c] - p.X0[3 * (i + 1) + c]));
    const d1 = Math.hypot(...[0, 1, 2].map(c => p.X[3 * i + c] - p.X[3 * (i + 1) + c]));
    worst = Math.max(worst, Math.abs(d1 - d0));
  }
  ok(worst < 0.6, `backbone neighbour distances preserved under an 8 Å pull (max change ${worst.toFixed(2)} Å)`);
  const disp = p.displacement();
  ok(disp[0] < disp[r] && disp[ens.nRes - 1] < disp[r], 'deformation is local to the pulled region');
  const out = p.apply(new Float32Array(ens.coords[0].length));
  ok(out.every(Number.isFinite), 'all-atom mapping produces finite coordinates');
  // pulling a residue to exactly where another conformation has it (no other change) must rank that conformation well
  const k = ens.n - 1;
  const q = new ElasticPull(ens, ens.coords[k], r, { pairs: ens.conformations[k].pairs });
  const ranked = rankConformations(ens, q);
  ok(ranked[0].k === k && ranked[0].score < 1e-6, 'an unpulled conformation matches itself');
  const pick = chooseStable(ens, ranked, doc.default_metric, k);
  ok(pick.choice.k === k, 'starting conformation is kept when it is the outright best match');
  const fc = fittedCoords(ens, 0, ranked.find(x => x.k === 0).fit);
  ok(fc.length === ens.coords[0].length && fc.every(Number.isFinite), 'fitted landing coordinates');
}

// ---- encoders
ok(crc32(new TextEncoder().encode('123456789')) === 0xcbf43926, 'CRC-32 check value');
ok(toWinAnsi('Å–·').join(',') === '197,150,183', 'WinAnsi mapping');
const w = 40, h = 30, data = new Uint8ClampedArray(w * h * 4);
for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) { const i = (y * w + x) * 4; data[i] = x * 6; data[i + 1] = y * 8; data[i + 2] = 120; data[i + 3] = x < 5 ? 0 : 255; }
const img = { width: w, height: h, data };
const save = async (name, blob) => writeFileSync(path.join(outDir, name), Buffer.from(await blob.arrayBuffer()));
await save('rgb.tif', await encodeTIFF(img, { dpi: 600 }));
await save('rgba.tif', await encodeTIFF(img, { dpi: 300, alpha: true }));
const png1 = Buffer.from('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg==', 'base64');
await save('dpi.png', await pngWithDpi(new Blob([png1]), 600));
const zip = new ZipWriter();
zip.add('frames/frame_00001.png', new Uint8Array(png1));
zip.add('notes.txt', new TextEncoder().encode('rna'));
await save('frames.zip', zip.blob());
const pdf = new PDFDoc(200, 120, { title: 'ΔG test' });
await pdf.image(img, 10, 10, 80, 60);
pdf.line(10, 100, 190, 100, { stroke: { r: 0, g: 0, b: 1 }, width: 1, dash: [3, 2] });
pdf.circle(150, 50, 15, { fill: { r: 1, g: 0.5, b: 0, a: 0.6 } });
pdf.text(100, 115, 'ΔG 2D (kcal/mol) · Å', { size: 8, anchor: 'middle' });
await save('figure.pdf', await pdf.blob());

console.log(failures ? `${failures} failure(s)` : 'all JS checks passed');
process.exit(failures ? 1 : 0);
