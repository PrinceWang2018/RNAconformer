// Pull-to-switch: elastic deformation of the displayed conformation while one residue is dragged, and the
// choice of the nearest *stable* conformation (a computed, validated ensemble member) for the pulled shape.
// DOM-free so it can be unit-tested in Node.
import { fitRigid } from './ensemble.js';

const dist = (X, i, j) => Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]);

export class ElasticPull {
  // ens: Ensemble; coords: displayed all-atom coordinates (Float32Array); residue: dragged residue index;
  // pairs: base pairs [[i, j], ...] of the displayed conformation (stiff springs).
  constructor(ens, coords, residue, { pairs = [], cutoff = 13, tether = 0.012 } = {}) {
    this.ens = ens;
    this.residue = residue;
    this.coords0 = Float32Array.from(coords);
    const n = ens.nRes;
    this.n = n;
    this.atom = Int32Array.from(ens.trace);
    this.valid = [];
    this.X0 = new Float64Array(n * 3);
    for (let r = 0; r < n; r++) {
      const a = this.atom[r];
      if (a < 0) continue;
      this.valid.push(r);
      for (let d = 0; d < 3; d++) this.X0[3 * r + d] = coords[3 * a + d];
    }
    this.X = Float64Array.from(this.X0);
    const partner = new Map();
    for (const [i, j] of pairs) { partner.set(i, j); partner.set(j, i); }
    const si = [], sj = [], sl = [], ss = [];
    for (let a = 0; a < this.valid.length; a++) {
      for (let b = a + 1; b < this.valid.length; b++) {
        const i = this.valid[a], j = this.valid[b];
        const L = dist(this.X0, i, j);
        const gap = j - i;
        let s;
        if (gap === 1) s = 1.0;                                  // backbone
        else if (partner.get(i) === j) s = 0.8;                 // base pair
        else if (gap === 2) s = 0.6;
        else if (L < cutoff) s = 0.22 * Math.min(1, (6 / Math.max(L, 6)) ** 2);
        else continue;
        si.push(i); sj.push(j); sl.push(L); ss.push(s);
      }
    }
    this.si = Int32Array.from(si); this.sj = Int32Array.from(sj);
    this.sl = Float64Array.from(sl); this.ss = Float64Array.from(ss);
    this.tether = tether;
    this.target = [this.X0[3 * residue], this.X0[3 * residue + 1], this.X0[3 * residue + 2]];
    // local frames for mapping atoms: neighbours along the chain and the pairing partner
    this.frames = new Array(n);
    for (let r = 0; r < n; r++) {
      const f = [r - 1, r, r + 1, partner.get(r)].filter(x => x !== undefined && x >= 0 && x < n && this.atom[x] >= 0);
      this.frames[r] = f;
    }
  }

  setTarget(t) { this.target = [t[0], t[1], t[2]]; }

  solve(iters = 30) {
    const X = this.X, X0 = this.X0, d = this.residue, T = this.target, k = this.tether;
    const fix = () => { X[3 * d] = T[0]; X[3 * d + 1] = T[1]; X[3 * d + 2] = T[2]; };
    for (let it = 0; it < iters; it++) {
      for (const i of this.valid) {
        if (i === d) continue;
        for (let c = 0; c < 3; c++) X[3 * i + c] += k * (X0[3 * i + c] - X[3 * i + c]);
      }
      fix();
      for (let s = 0; s < this.si.length; s++) {
        const i = this.si[s], j = this.sj[s];
        const dx = X[3 * j] - X[3 * i], dy = X[3 * j + 1] - X[3 * i + 1], dz = X[3 * j + 2] - X[3 * i + 2];
        const len = Math.hypot(dx, dy, dz);
        if (len < 1e-9) continue;
        const C = len - this.sl[s];
        if (Math.abs(C) < 1e-6) continue;
        const wi = i === d ? 0 : 1, wj = j === d ? 0 : 1, ws = wi + wj;
        if (!ws) continue;
        const f = this.ss[s] * C / len / ws;
        X[3 * i] += dx * f * wi; X[3 * i + 1] += dy * f * wi; X[3 * i + 2] += dz * f * wi;
        X[3 * j] -= dx * f * wj; X[3 * j + 1] -= dy * f * wj; X[3 * j + 2] -= dz * f * wj;
      }
      fix();
    }
  }

  // Deformed trace as [[x, y, z], ...] (one per residue; undeformed position where no trace atom exists).
  trace() {
    const out = [];
    for (let r = 0; r < this.n; r++) out.push([this.X[3 * r], this.X[3 * r + 1], this.X[3 * r + 2]]);
    return out;
  }

  displacement() {
    const d = new Float64Array(this.n);
    for (const r of this.valid) d[r] = Math.hypot(this.X[3 * r] - this.X0[3 * r], this.X[3 * r + 1] - this.X0[3 * r + 1], this.X[3 * r + 2] - this.X0[3 * r + 2]);
    return d;
  }

  // Map the deformation onto every atom: each residue moves rigidly with its local frame.
  apply(out) {
    const ens = this.ens, c0 = this.coords0, X0 = this.X0, X = this.X;
    for (let r = 0; r < this.n; r++) {
      const res = ens.residues[r];
      const f = this.frames[r];
      if (!f.length) continue;
      const P = f.map(i => [X0[3 * i], X0[3 * i + 1], X0[3 * i + 2]]);
      const Q = f.map(i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]]);
      const w = f.map(i => (i === r ? 2 : 1));
      const { R, cp, cq } = fitRigid(P, Q, w);
      for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) {
        const x = c0[3 * a] - cp[0], y = c0[3 * a + 1] - cp[1], z = c0[3 * a + 2] - cp[2];
        out[3 * a] = R[0][0] * x + R[0][1] * y + R[0][2] * z + cq[0];
        out[3 * a + 1] = R[1][0] * x + R[1][1] * y + R[1][2] * z + cq[1];
        out[3 * a + 2] = R[2][0] * x + R[2][1] * y + R[2][2] * z + cq[2];
      }
    }
    return out;
  }
}

// Score every conformation against the pulled shape. Each conformation is fitted on the residues that the
// pull did not move, then scored on how well it reproduces the moved region and the dragged residue.
// alpha: weight of the directly manipulated region (dragged residue and its chain neighbours) against the
// weighted RMSD over the elastically moved region (which is only a guide, not physics).
// viewDir: unit view direction in model coordinates. A mouse drag moves the residue only within the screen
// plane, so for the dragged region the error along the view direction counts only `depthWeight`.
export function rankConformations(ens, pull, { alpha = 0.85, viewDir = null, depthWeight = 0.1 } = {}) {
  const X = pull.X, n = pull.n, d = pull.residue;
  const disp = pull.displacement();
  const maxD = Math.max(1e-6, ...disp);
  const valid = pull.valid;
  const fw = valid.map(r => 1 - 0.85 * disp[r] / maxD);
  const sw = valid.map(r => 0.2 + disp[r] / maxD);
  const target = valid.map(r => [X[3 * r], X[3 * r + 1], X[3 * r + 2]]);
  const out = [];
  for (let k = 0; k < ens.n; k++) {
    const C = ens.coords[k];
    const Y = valid.map(r => { const a = ens.trace[r]; return [C[3 * a], C[3 * a + 1], C[3 * a + 2]]; });
    const { R, cp, cq } = fitRigid(Y, target, fw);
    let s2 = 0, sws = 0, reach = 0, l2 = 0, ln = 0;
    Y.forEach((p, i) => {
      const x = p[0] - cp[0], y = p[1] - cp[1], z = p[2] - cp[2];
      const q = [R[0][0] * x + R[0][1] * y + R[0][2] * z + cq[0], R[1][0] * x + R[1][1] * y + R[1][2] * z + cq[1], R[2][0] * x + R[2][1] * y + R[2][2] * z + cq[2]];
      const e2 = (q[0] - target[i][0]) ** 2 + (q[1] - target[i][1]) ** 2 + (q[2] - target[i][2]) ** 2;
      s2 += sw[i] * e2; sws += sw[i];
      let ep2 = e2;
      if (viewDir) {
        const dz = (q[0] - target[i][0]) * viewDir[0] + (q[1] - target[i][1]) * viewDir[1] + (q[2] - target[i][2]) * viewDir[2];
        ep2 = e2 - (1 - depthWeight) * dz * dz;
      }
      if (valid[i] === d) reach = Math.sqrt(ep2);
      if (Math.abs(valid[i] - d) <= 1) { l2 += (valid[i] === d ? 2 : 1) * ep2; ln += valid[i] === d ? 2 : 1; }
    });
    const wrmsd = Math.sqrt(s2 / (sws || 1));
    const local = Math.sqrt(l2 / (ln || 1));
    const v = ens.conformations[k].validation;
    out.push({ k, score: alpha * local + (1 - alpha) * wrmsd, local, wrmsd, reach, valid: v ? v.passed !== false : true, fit: { R, cp, cq } });
  }
  out.sort((a, b) => a.score - b.score);
  return out;
}

// Among conformations that match about as well as the best one (within max(0.75 Å, 10 %)), prefer the most
// stable by the chosen metric — "relatively stable, not necessarily the lowest".
// `start`: the conformation the pull began from; it is kept only when it is itself the best match.
export function chooseStable(ens, ranked, metricKey, start = null) {
  const ok = ranked.filter(r => r.valid);
  const pool = ok.length ? ok : ranked;
  const best = pool[0];
  const tol = Math.max(0.75, 0.1 * best.score);
  const near = pool.filter(r => r.score <= best.score + tol && (r.k !== start || r === best));
  const m = ens.metricByKey[metricKey];
  let choice = best;
  if (m && (m.lower_is_better === true || m.lower_is_better === false)) {
    const sign = m.lower_is_better ? 1 : -1;
    const scored = near.filter(r => ens.value(r.k, metricKey) !== null);
    if (scored.length) choice = scored.reduce((a, b) => (sign * ens.value(b.k, metricKey) < sign * ens.value(a.k, metricKey) ? b : a));
  }
  return { choice, best, near, tolerance: tol };
}

// All-atom coordinates of conformation k placed in the frame used for matching (fitted on the residues the
// pull left in place), so the parts the user did not move stay where they are when the structure switches.
export function fittedCoords(ens, k, fit) {
  const C = ens.coords[k], out = new Float32Array(C.length), { R, cp, cq } = fit;
  for (let i = 0; i < C.length; i += 3) {
    const x = C[i] - cp[0], y = C[i + 1] - cp[1], z = C[i + 2] - cp[2];
    out[i] = R[0][0] * x + R[0][1] * y + R[0][2] * z + cq[0];
    out[i + 1] = R[1][0] * x + R[1][1] * y + R[1][2] * z + cq[1];
    out[i + 2] = R[2][0] * x + R[2][1] * y + R[2][2] * z + cq[2];
  }
  return out;
}
