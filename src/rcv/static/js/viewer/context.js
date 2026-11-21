// What is bound to the RNA (document.context, see rcv/data.py): ligands as sticks, ions as spheres, and protein, DNA
// and other RNA chains as thin traces (one point per residue). The layer follows the shown conformation and morphs
// with it; it is drawn in the scene, so images, figures and videos include it.
import * as THREE from '../../vendor/three/three.module.min.js';
import { Cylinders, Spheres } from './geometry.js';
import { ELEMENT } from '../colors.js';

export const CONTEXT_COLORS = {
  carbon: '#b5519c',                                   // ligand carbons: a hue no RNA colour scheme uses
  other: '#8f98a3',
  protein: '#b08d6e', dna: '#6f93c4', rna: '#98a1ab',
  ions: { MG: '#3fa34d', NA: '#8e5bd0', K: '#7a4bb5', CA: '#5c8f2a', ZN: '#6d7f99', MN: '#9c6fb3', FE: '#c8642a', CO: '#d4739a', NI: '#4fae7b',
    CU: '#b8732e', CL: '#3fae49', SR: '#2f9d8a', BA: '#2e8b57', CD: '#c9a15a', CS: '#6a3fa0', RB: '#7a3fa0', TL: '#8a6d5a', PB: '#575961',
    PT: '#8f98a3', IR: '#3d5a9c', OS: '#3d6a9c', HG: '#9aa3ad', AU: '#d4a017', LI: '#b37fd6', BR: '#a1522c', I: '#7b3fa0', F: '#5fb34f' },
};
const KIND_LABEL = { ligand: 'ligand', ion: 'ion', protein: 'protein', dna: 'DNA', rna: 'other RNA chain' };
const R_STICK = 0.3, R_ION = 0.9, R_TRACE = 0.32;

export class ContextLayer {
  constructor(parent) {
    this.parent = parent;
    this.group = new THREE.Group();
    parent.add(this.group);
    this.ctx = null;
    this.shown = null;                 // Float32Array drawn now (atoms, then trace points)
    this.enabled = true;
  }

  get active() { return !!this.ctx && this.enabled; }
  get radius() { return this.active ? this.ctx.radius : 0; }

  build(ctx, material) {
    this.dispose();
    this.ctx = ctx;
    if (!ctx) return;
    this.material = material;
    const nIons = ctx.groups.filter(g => g.kind === 'ion').length;
    this.atoms = new Spheres(Math.max(1, ctx.nAtoms - nIons), material, 2);
    this.ions = new Spheres(Math.max(1, nIons), material, 2);
    this.bonds = new Cylinders(Math.max(1, ctx.bonds.length), material, 10);
    this.joints = new Spheres(Math.max(1, ctx.nPoints), material, 1);
    this.segs = new Cylinders(Math.max(1, ctx.segs.length), material, 6);
    this.ligandAtoms = [];
    this.ionAtoms = [];
    for (const g of ctx.groups) for (let a = g.first_atom; a < g.first_atom + g.n_atoms; a++) (g.kind === 'ion' ? this.ionAtoms : this.ligandAtoms).push(a);
    for (const m of [this.atoms, this.ions, this.bonds, this.joints, this.segs]) this.group.add(m.mesh);
    this.group.visible = this.enabled;
    this.shown = ctx.coords[0].slice();
    this._draw();
    this._colors();
  }

  setMaterial(material) {
    this.material = material;
    if (!this.ctx) return;
    for (const m of [this.atoms, this.ions, this.bonds, this.joints, this.segs]) m.mesh.material = material;
  }

  setEnabled(on) { this.enabled = !!on; this.group.visible = this.enabled; }

  // Shown coordinates: conformation a, or the blend a -> b at t (0..1).
  frame(a, b = a, t = 0) {
    const ctx = this.ctx;
    if (!ctx) return;
    const A = ctx.coords[a], B = ctx.coords[b];
    if (!A || !B) return;
    const d = this.shown;
    if (a === b || t <= 0) d.set(A);
    else for (let i = 0; i < d.length; i++) d[i] = A[i] + (B[i] - A[i]) * t;
    this._draw();
  }
  // From a snapshot of the shown coordinates (a morph that started mid-way) to conformation k.
  blendFrom(start, k, t) {
    const ctx = this.ctx;
    if (!ctx || !start || start.length !== this.shown.length || !ctx.coords[k]) return;
    const B = ctx.coords[k], d = this.shown;
    for (let i = 0; i < d.length; i++) d[i] = start[i] + (B[i] - start[i]) * t;
    this._draw();
  }
  snapshot() { return this.shown ? this.shown.slice() : null; }

  _draw() {
    const ctx = this.ctx, c = this.shown;
    if (!ctx) return;
    const at = (list, i, p) => { const a = 3 * list[i]; p[0] = c[a]; p[1] = c[a + 1]; p[2] = c[a + 2]; };
    this.atoms.set(this.ligandAtoms.length, (i, p) => at(this.ligandAtoms, i, p), () => R_STICK);
    this.ions.set(this.ionAtoms.length, (i, p) => at(this.ionAtoms, i, p), () => R_ION);
    const b = ctx.bonds;                                   // half-bond i: atom b[i] -> middle of the bond with b[i ^ 1]
    this.bonds.set(b.length, (i, seg) => {
      const p = 3 * b[i], q = 3 * b[i ^ 1];
      seg[0] = c[p]; seg[1] = c[p + 1]; seg[2] = c[p + 2];
      seg[3] = (c[p] + c[q]) / 2; seg[4] = (c[p + 1] + c[q + 1]) / 2; seg[5] = (c[p + 2] + c[q + 2]) / 2;
    }, R_STICK);
    const o = 3 * ctx.nAtoms;
    this.joints.set(ctx.nPoints, (i, p) => { p[0] = c[o + 3 * i]; p[1] = c[o + 3 * i + 1]; p[2] = c[o + 3 * i + 2]; }, () => R_TRACE);
    const s = ctx.segs;
    this.segs.set(s.length, (i, seg) => { const p = o + 3 * s[i]; for (let k = 0; k < 6; k++) seg[k] = c[p + k]; }, R_TRACE);
  }

  _colors() {
    const ctx = this.ctx, col = new THREE.Color(), lin = {};
    const of = hex => { if (!lin[hex]) { col.set(hex); lin[hex] = [col.r, col.g, col.b]; } return lin[hex]; };
    const atomHex = a => {
      const e = ctx.elements[a];
      if (ctx.groups[ctx.group[a]].kind === 'ion') return CONTEXT_COLORS.ions[e] || CONTEXT_COLORS.other;
      return e === 'C' ? CONTEXT_COLORS.carbon : (ELEMENT[e] || CONTEXT_COLORS.ions[e] || CONTEXT_COLORS.other);
    };
    const put = (hex, out) => { const v = of(hex); out[0] = v[0]; out[1] = v[1]; out[2] = v[2]; };
    this.atoms.setColors(this.ligandAtoms.length, (i, out) => put(atomHex(this.ligandAtoms[i]), out));
    this.ions.setColors(this.ionAtoms.length, (i, out) => put(atomHex(this.ionAtoms[i]), out));
    this.bonds.setColors(ctx.bonds.length, (i, out) => put(atomHex(ctx.bonds[i]), out));
    const traceHex = p => CONTEXT_COLORS[ctx.traces[ctx.pointTrace[p]].kind] || CONTEXT_COLORS.other;
    this.joints.setColors(ctx.nPoints, (i, out) => put(traceHex(i), out));
    this.segs.setColors(ctx.segs.length, (i, out) => put(traceHex(ctx.segs[i]), out));
  }

  // What is under the pointer: { kind, name, what, chain, color, score } (score: depth, as Viewer.pick) or null.
  pick(proj, x, y, unit) {
    const ctx = this.ctx;
    if (!this.active || !this.shown) return null;
    const c = this.shown;
    let best = null;
    const test = (i, radius, make) => {
      const [px, py, pz] = proj(c[3 * i], c[3 * i + 1], c[3 * i + 2]);
      const d = Math.hypot(px - x, py - y);
      if (d > Math.max(6, radius / unit)) return;
      const score = pz + d * 1e-4;
      if (!best || score < best.score) best = { score, make };
    };
    for (let a = 0; a < ctx.nAtoms; a++) {
      const g = ctx.groups[ctx.group[a]];
      test(a, g.kind === 'ion' ? R_ION : 0.6, () => ({ kind: g.kind, what: KIND_LABEL[g.kind], name: `${g.comp} ${g.num}`.trim(), chain: g.chain, atom: g.kind === 'ion' ? '' : ctx.names[a],
        color: g.kind === 'ion' ? (CONTEXT_COLORS.ions[ctx.elements[a]] || CONTEXT_COLORS.other) : CONTEXT_COLORS.carbon }));
    }
    for (let p = 0; p < ctx.nPoints; p++) {
      const t = ctx.traces[ctx.pointTrace[p]];
      test(ctx.nAtoms + p, 0.8, () => ({ kind: t.kind, what: KIND_LABEL[t.kind], name: `${t.comps[p - t.first]} ${t.nums[p - t.first]}`, chain: t.chain, atom: '', color: CONTEXT_COLORS[t.kind] }));
    }
    return best ? { ...best.make(), score: best.score } : null;
  }

  dispose() {
    for (const m of [this.atoms, this.ions, this.bonds, this.joints, this.segs]) {
      if (!m) continue;
      this.group.remove(m.mesh);
      m.mesh.geometry.dispose();
    }
    this.atoms = this.ions = this.bonds = this.joints = this.segs = null;
    this.ctx = null;
    this.shown = null;
  }
}
