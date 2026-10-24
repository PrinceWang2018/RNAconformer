// Client-side model of one ensemble document (schema rcv.document/1; legacy rnaconfswitch.document/1, rnafold3d.viewer/1).
import { b64ToFloat32 } from './util.js';

export const DOC_SCHEMAS = ['rcv.document/1', 'rnaconfswitch.document/1', 'rnafold3d.viewer/1'];
const PURINE_OUTLINE = ['C4', 'N9', 'C8', 'N7', 'C5', 'C6', 'N1', 'C2', 'N3'];
const SIX = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'];
const FIVE = ['C4', 'N9', 'C8', 'N7', 'C5'];      // same rotational sense as the purine outline

export class Ensemble {
  constructor(doc, opts = {}) {
    if (!doc || !DOC_SCHEMAS.includes(doc.schema)) throw new Error('not an RCV ensemble document');
    this.doc = doc;
    this.id = doc.id;
    this.serverId = opts.serverId || null;
    this.name = doc.name;
    this.sequence = doc.sequence;
    this.residues = doc.residues;
    this.nRes = doc.residues.length;
    this.atomNames = doc.atoms.names;
    this.elements = doc.atoms.elements;
    this.atomRes = Int32Array.from(doc.atoms.residue);
    this.nAtoms = this.atomNames.length;
    this.bonds = Int32Array.from(doc.bonds);
    this.conformations = doc.conformations;
    this.n = doc.conformations.length;
    this.coords = doc.conformations.map(c => b64ToFloat32(c.coords));
    this.rmsf = doc.ensemble.rmsf;
    this.rmsd = doc.ensemble.rmsd_c4;
    this.embedding = doc.ensemble.embedding;
    this.metrics = doc.metrics;
    this.metricByKey = Object.fromEntries(doc.metrics.map(m => [m.key, m]));
    this.terms = doc.terms || {};
    this._indexAtoms();
    this._bounds();
  }

  _indexAtoms() {
    const n = this.nRes;
    this.atomIndex = new Array(n);
    this.trace = new Int32Array(n).fill(-1);       // C4' (tube)
    this.c1 = new Int32Array(n).fill(-1);
    this.glyco = new Int32Array(n).fill(-1);       // N9 (purine) / N1 (pyrimidine)
    this.phos = new Int32Array(n).fill(-1);
    this.rings = new Array(n);                     // {outline:[idx], fans:[[idx...]]}
    for (let r = 0; r < n; r++) {
      const res = this.residues[r];
      const map = {};
      for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) map[this.atomNames[a]] = a;
      this.atomIndex[r] = map;
      const pick = (...names) => { for (const nm of names) if (nm in map) return map[nm]; return -1; };
      this.trace[r] = pick("C4'", "C3'", "C1'", 'P');
      this.c1[r] = pick("C1'");
      this.phos[r] = pick('P');
      const purine = 'N9' in map && 'C8' in map && 'N7' in map;
      this.glyco[r] = purine ? map.N9 : pick('N1', 'N9');
      const has = list => list.every(nm => nm in map);
      if (purine && has(PURINE_OUTLINE)) {
        this.rings[r] = { outline: PURINE_OUTLINE.map(nm => map[nm]), fans: [SIX.map(nm => map[nm]), FIVE.map(nm => map[nm])] };
      } else if (has(SIX)) {
        this.rings[r] = { outline: SIX.map(nm => map[nm]), fans: [SIX.map(nm => map[nm])] };
      } else this.rings[r] = null;
    }
  }

  _bounds() {
    let cx = 0, cy = 0, cz = 0, cnt = 0;
    for (const c of this.coords) for (let i = 0; i < c.length; i += 3) { cx += c[i]; cy += c[i + 1]; cz += c[i + 2]; cnt++; }
    cx /= cnt; cy /= cnt; cz /= cnt;
    let r2 = 0;
    for (const c of this.coords) for (let i = 0; i < c.length; i += 3) {
      const dx = c[i] - cx, dy = c[i + 1] - cy, dz = c[i + 2] - cz;
      r2 = Math.max(r2, dx * dx + dy * dy + dz * dz);
    }
    this.center = [cx, cy, cz];
    this.radius = Math.sqrt(r2) + 2.5;
    // framing radius: median over conformations of each conformation's extent around the common centre
    const ext = this.coords.map(c => {
      let m = 0;
      for (let i = 0; i < c.length; i += 3) { const dx = c[i] - cx, dy = c[i + 1] - cy, dz = c[i + 2] - cz; m = Math.max(m, dx * dx + dy * dy + dz * dz); }
      return Math.sqrt(m);
    }).sort((a, b) => a - b);
    this.fitRadius = Math.min(this.radius, ext[Math.floor(ext.length / 2)] + 2.0);
  }

  label(k) { return this.conformations[k]?.label ?? `conf ${k + 1}`; }
  // Short display name: imported models keep their own label (e.g. "Model 7"), predictions their conf id.
  shortId(k) {
    const c = this.conformations[k];
    if (!c) return `conf_${k + 1}`;
    if (this.doc.origin !== 'rnafold3d-engine' && c.label && c.label.length <= 28) return c.label;
    return c.id;
  }
  value(k, key) { const v = this.conformations[k]?.values?.[key]; return typeof v === 'number' ? v : null; }
  resLabel(r) { const res = this.residues[r]; return `${res.name}${res.num}`; }

  partner(k) {
    const p = new Int32Array(this.nRes).fill(-1);
    for (const [i, j] of this.conformations[k].pairs) { p[i] = j; p[j] = i; }
    return p;
  }
  pairInfo(k, r) {
    for (const [i, j, t, pk] of this.conformations[k].pairs) if (i === r || j === r) return { partner: i === r ? j : i, type: t, pk };
    return null;
  }

  residueCenter(coords, r) {
    const res = this.residues[r];
    let x = 0, y = 0, z = 0;
    for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) { x += coords[3 * a]; y += coords[3 * a + 1]; z += coords[3 * a + 2]; }
    const m = Math.max(1, res.n_atoms);
    return [x / m, y / m, z / m];
  }

  residueDisplacement(a, b) {
    const out = new Array(this.nRes).fill(0);
    if (!a || !b) return out;
    for (let r = 0; r < this.nRes; r++) {
      const res = this.residues[r];
      let s = 0;
      for (let i = res.first_atom; i < res.first_atom + res.n_atoms; i++) {
        const dx = a[3 * i] - b[3 * i], dy = a[3 * i + 1] - b[3 * i + 1], dz = a[3 * i + 2] - b[3 * i + 2];
        s += dx * dx + dy * dy + dz * dz;
      }
      out[r] = Math.sqrt(s / Math.max(1, res.n_atoms));
    }
    return out;
  }

  // Conformation display orders.
  order(mode, metricKey) {
    const idx = [...Array(this.n).keys()];
    if (mode === 'metric' && metricKey) {
      const m = this.metricByKey[metricKey];
      const sign = m && m.lower_is_better === false ? -1 : 1;
      return idx.sort((a, b) => sign * ((this.value(a, metricKey) ?? Infinity) - (this.value(b, metricKey) ?? Infinity)));
    }
    if (mode === 'rank') return idx.sort((a, b) => (this.conformations[a].rank ?? a) - (this.conformations[b].rank ?? b));
    if (mode === 'path') {
      // greedy nearest-neighbour walk through C4' RMSD space starting from the best-ranked conformation
      const start = this.order('rank')[0];
      const out = [start];
      const left = new Set(idx.filter(i => i !== start));
      while (left.size) {
        const last = out[out.length - 1];
        let best = -1, bd = Infinity;
        for (const j of left) { const d = this.rmsd[last][j]; if (d < bd) { bd = d; best = j; } }
        out.push(best);
        left.delete(best);
      }
      return out;
    }
    return idx;
  }

  // Superpose every conformation onto `ref` using the atoms of the given residues (quaternion method).
  superposeOn(residueList, ref = 0) {
    const atoms = [];
    for (const r of residueList) { const a = this.trace[r]; if (a >= 0) atoms.push(a); }
    if (atoms.length < 3) for (const r of residueList) { const res = this.residues[r]; for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) atoms.push(a); }
    if (atoms.length < 3) throw new Error('select at least 3 residues to superpose on');
    const target = this.coords[ref];
    for (let k = 0; k < this.n; k++) {
      if (k === ref) continue;
      const { R, t } = kabsch(this.coords[k], target, atoms);
      applyTransform(this.coords[k], R, t);
    }
    this._bounds();
  }

  rmsdBetween(a, b, atoms) {
    const A = this.coords[a], B = this.coords[b];
    let s = 0;
    for (const i of atoms) { const dx = A[3 * i] - B[3 * i], dy = A[3 * i + 1] - B[3 * i + 1], dz = A[3 * i + 2] - B[3 * i + 2]; s += dx * dx + dy * dy + dz * dz; }
    return Math.sqrt(s / Math.max(1, atoms.length));
  }
}

// Optimal rotation (Horn's quaternion method, Jacobi eigen-solver) mapping mobile[idx] onto target[idx].
export function kabsch(mobile, target, idx) {
  let mx = 0, my = 0, mz = 0, tx = 0, ty = 0, tz = 0;
  for (const i of idx) { mx += mobile[3 * i]; my += mobile[3 * i + 1]; mz += mobile[3 * i + 2]; tx += target[3 * i]; ty += target[3 * i + 1]; tz += target[3 * i + 2]; }
  const n = idx.length;
  mx /= n; my /= n; mz /= n; tx /= n; ty /= n; tz /= n;
  let Sxx = 0, Sxy = 0, Sxz = 0, Syx = 0, Syy = 0, Syz = 0, Szx = 0, Szy = 0, Szz = 0;
  for (const i of idx) {
    const ax = mobile[3 * i] - mx, ay = mobile[3 * i + 1] - my, az = mobile[3 * i + 2] - mz;
    const bx = target[3 * i] - tx, by = target[3 * i + 1] - ty, bz = target[3 * i + 2] - tz;
    Sxx += ax * bx; Sxy += ax * by; Sxz += ax * bz; Syx += ay * bx; Syy += ay * by; Syz += ay * bz; Szx += az * bx; Szy += az * by; Szz += az * bz;
  }
  const N = [
    [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
    [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
    [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
    [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
  ];
  const { values, vectors } = jacobi4(N);
  let best = 0;
  for (let i = 1; i < 4; i++) if (values[i] > values[best]) best = i;
  const [q0, q1, q2, q3] = [vectors[0][best], vectors[1][best], vectors[2][best], vectors[3][best]];
  const R = [
    [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
    [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
    [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
  ];
  const t = [tx - (R[0][0] * mx + R[0][1] * my + R[0][2] * mz), ty - (R[1][0] * mx + R[1][1] * my + R[1][2] * mz), tz - (R[2][0] * mx + R[2][1] * my + R[2][2] * mz)];
  return { R, t };
}

// Weighted rigid fit of point lists (arrays of [x, y, z]): Q ≈ R (P - cp) + cq. Weights default to 1.
export function fitRigid(P, Q, w = null) {
  const n = P.length;
  let W = 0;
  const cp = [0, 0, 0], cq = [0, 0, 0];
  for (let i = 0; i < n; i++) {
    const wi = w ? w[i] : 1;
    W += wi;
    for (let d = 0; d < 3; d++) { cp[d] += wi * P[i][d]; cq[d] += wi * Q[i][d]; }
  }
  for (let d = 0; d < 3; d++) { cp[d] /= W || 1; cq[d] /= W || 1; }
  if (n < 2) return { R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], cp, cq };
  let Sxx = 0, Sxy = 0, Sxz = 0, Syx = 0, Syy = 0, Syz = 0, Szx = 0, Szy = 0, Szz = 0;
  for (let i = 0; i < n; i++) {
    const wi = w ? w[i] : 1;
    const ax = P[i][0] - cp[0], ay = P[i][1] - cp[1], az = P[i][2] - cp[2];
    const bx = Q[i][0] - cq[0], by = Q[i][1] - cq[1], bz = Q[i][2] - cq[2];
    Sxx += wi * ax * bx; Sxy += wi * ax * by; Sxz += wi * ax * bz; Syx += wi * ay * bx; Syy += wi * ay * by;
    Syz += wi * ay * bz; Szx += wi * az * bx; Szy += wi * az * by; Szz += wi * az * bz;
  }
  return { R: quatRotation(Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz), cp, cq };
}

function quatRotation(Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz) {
  const N = [
    [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
    [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
    [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
    [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
  ];
  const { values, vectors } = jacobi4(N);
  let best = 0;
  for (let i = 1; i < 4; i++) if (values[i] > values[best]) best = i;
  const [q0, q1, q2, q3] = [vectors[0][best], vectors[1][best], vectors[2][best], vectors[3][best]];
  return [
    [q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
    [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
    [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3],
  ];
}

export function applyTransform(c, R, t) {
  for (let i = 0; i < c.length; i += 3) {
    const x = c[i], y = c[i + 1], z = c[i + 2];
    c[i] = R[0][0] * x + R[0][1] * y + R[0][2] * z + t[0];
    c[i + 1] = R[1][0] * x + R[1][1] * y + R[1][2] * z + t[1];
    c[i + 2] = R[2][0] * x + R[2][1] * y + R[2][2] * z + t[2];
  }
}

function jacobi4(A) {
  const a = A.map(r => r.slice());
  const v = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
  for (let sweep = 0; sweep < 50; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
    if (off < 1e-18) break;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) {
      if (Math.abs(a[p][q]) < 1e-15) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 4; k++) {
        const akp = a[k][p], akq = a[k][q];
        a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq;
      }
      for (let k = 0; k < 4; k++) {
        const apk = a[p][k], aqk = a[q][k];
        a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk;
      }
      for (let k = 0; k < 4; k++) {
        const vkp = v[k][p], vkq = v[k][q];
        v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq;
      }
    }
  }
  return { values: [a[0][0], a[1][1], a[2][2], a[3][3]], vectors: v };
}

// PDB writer for the displayed coordinates.
export function toPDB(ens, coords, { model = null, remarks = [] } = {}) {
  const lines = remarks.map(r => `REMARK 250 ${r}`.slice(0, 80));
  if (model !== null) lines.push(`MODEL     ${String(model).padStart(4)}`);
  let serial = 1;
  for (let r = 0; r < ens.nRes; r++) {
    const res = ens.residues[r];
    const resName = (res.comp || res.name).slice(0, 3).padStart(3);
    const num = parseInt(res.num, 10);
    const ins = String(res.num).replace(/^-?\d+/, '').slice(0, 1) || ' ';
    const chain = (res.chain || 'A').slice(0, 1);
    for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) {
      const name = ens.atomNames[a];
      const el = ens.elements[a];
      const nm = name.length >= 4 ? name.slice(0, 4) : (el.length === 1 ? ' ' + name : name).padEnd(4);
      const x = coords[3 * a].toFixed(3).padStart(8), y = coords[3 * a + 1].toFixed(3).padStart(8), z = coords[3 * a + 2].toFixed(3).padStart(8);
      lines.push(`ATOM  ${String(serial % 100000).padStart(5)} ${nm} ${resName} ${chain}${String(Number.isFinite(num) ? num : r + 1).padStart(4)}${ins}   ${x}${y}${z}  1.00  0.00          ${el.padStart(2)}`);
      serial++;
    }
  }
  lines.push('TER');
  if (model !== null) lines.push('ENDMDL');
  return lines.join('\n') + '\n';
}
