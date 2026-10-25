// Mechanism morph between two conformations: rigid modules (from the Pathway analysis) turn about their own centres
// along the shortest rotation and translate, instead of every atom moving on a straight line (which shrinks helices
// that rotate). In the staged form each module moves during its own time window, taken from the order of the
// base-pair events on the 2D refolding route, so the animation shows the mechanism step by step. Flexible residues
// (hinges, rearranging strands) blend the motions of the neighbouring modules along the chain plus their own
// residual displacement, so the chain stays connected; both end points are reproduced exactly.
// This is an illustration of the analysed mechanism, not a simulated trajectory.

const MIN_WINDOW = 0.35;

export function smooth(x) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  return x * x * (3 - 2 * x);
}

// q^s for a unit quaternion q = (w, x, y, z): the rotation by s times the angle, about the same axis (shortest arc).
export function quatPow(q, s) {
  let [w, x, y, z] = q;
  if (w < 0) { w = -w; x = -x; y = -y; z = -z; }
  const half = Math.acos(Math.min(1, Math.max(-1, w)));
  const sn = Math.sin(half);
  if (sn < 1e-9) return [1, 0, 0, 0];
  const h = half * s, k = Math.sin(h) / sn;
  return [Math.cos(h), x * k, y * k, z * k];
}

// rotate vector v by unit quaternion q (column-vector convention)
export function quatRotate(q, v) {
  const [w, x, y, z] = q;
  const tx = 2 * (y * v[2] - z * v[1]), ty = 2 * (z * v[0] - x * v[2]), tz = 2 * (x * v[1] - y * v[0]);
  return [v[0] + w * tx + (y * tz - z * ty), v[1] + w * ty + (z * tx - x * tz), v[2] + w * tz + (x * ty - y * tx)];
}

function windowOf(w) {
  if (!w) return null;
  let [a, b] = w;
  const c = (a + b) / 2, half = Math.max((b - a) / 2, MIN_WINDOW / 2);
  a = Math.max(0, c - half); b = Math.min(1, c + half);
  if (b - a < MIN_WINDOW) { if (a === 0) b = MIN_WINDOW; else a = 1 - MIN_WINDOW; }
  return [a, b];
}

// Build the morph. A, B: Float32Array coordinates (x, y, z per atom) of the two conformations; atomRes: residue per
// atom; motion: result.motion of the Pathway analysis (modules, blend); timing: result.timing (per residue
// first/last event time); staged: modules move in the order of the 2D events, else all together.
export function buildMechanism(A, B, atomRes, motion, timing, { staged = true, names = null } = {}) {
  const mods = motion?.modules || [];
  const nAtoms = A.length / 3;
  const nRes = motion?.blend?.length || 0;
  const win = mods.map(m => (staged ? windowOf(m.window) : null) || [staged ? 0.15 : 0, staged ? 0.85 : 1]);
  // residues: their own window (events) or the blend of their modules' windows
  const resWin = [];
  for (let r = 0; r < nRes; r++) {
    const own = staged && timing?.[r] && timing[r].first !== null && timing[r].first !== undefined ? windowOf([timing[r].first, timing[r].last]) : null;
    if (own && (motion.module_of?.[r] ?? -1) < 0) { resWin.push(own); continue; }
    const bl = motion.blend[r] || [];
    if (!bl.length) { resWin.push(staged ? [0.1, 0.9] : [0, 1]); continue; }
    let a = 0, b = 0;
    for (const [m, w] of bl) { a += w * win[m][0]; b += w * win[m][1]; }
    resWin.push([a, b]);
  }
  // per atom: rigid start positions relative to each blended module and the residual to B
  const blendOf = new Array(nAtoms);
  const residual = new Float32Array(A.length);
  for (let i = 0; i < nAtoms; i++) {
    const r = atomRes[i];
    const bl = (motion?.blend?.[r] || []).filter(([m]) => mods[m]);
    blendOf[i] = bl;
    const v = [A[3 * i], A[3 * i + 1], A[3 * i + 2]];
    let px = 0, py = 0, pz = 0, wsum = 0;
    for (const [m, w] of bl) {
      const M = mods[m];
      const d = [v[0] - M.centre_a[0], v[1] - M.centre_a[1], v[2] - M.centre_a[2]];
      const rot = quatRotate(M.quaternion, d);
      px += w * (rot[0] + M.centre_b[0]); py += w * (rot[1] + M.centre_b[1]); pz += w * (rot[2] + M.centre_b[2]);
      wsum += w;
    }
    if (wsum > 0) {
      residual[3 * i] = B[3 * i] - px / wsum; residual[3 * i + 1] = B[3 * i + 1] - py / wsum; residual[3 * i + 2] = B[3 * i + 2] - pz / wsum;
    } else {
      residual[3 * i] = B[3 * i] - v[0]; residual[3 * i + 1] = B[3 * i + 1] - v[1]; residual[3 * i + 2] = B[3 * i + 2] - v[2];
    }
  }
  const progress = (w, t) => smooth((t - w[0]) / Math.max(1e-6, w[1] - w[0]));
  // rigid-module estimate at time t (also the start of the refinement below)
  function rigidAt(t, out = new Float32Array(A.length)) {
    const qs = mods.map((M, m) => quatPow(M.quaternion, progress(win[m], t)));
    const ss = mods.map((M, m) => progress(win[m], t));
    const sr = resWin.map(w => progress(w, t));
    for (let i = 0; i < nAtoms; i++) {
      const bl = blendOf[i];
      const x = A[3 * i], y = A[3 * i + 1], z = A[3 * i + 2];
      let px = 0, py = 0, pz = 0, wsum = 0;
      for (const [m, w] of bl) {
        const M = mods[m];
        const rot = quatRotate(qs[m], [x - M.centre_a[0], y - M.centre_a[1], z - M.centre_a[2]]);
        const s = ss[m];
        px += w * (rot[0] + M.centre_a[0] + s * (M.centre_b[0] - M.centre_a[0]));
        py += w * (rot[1] + M.centre_a[1] + s * (M.centre_b[1] - M.centre_a[1]));
        pz += w * (rot[2] + M.centre_a[2] + s * (M.centre_b[2] - M.centre_a[2]));
        wsum += w;
      }
      if (wsum > 0) { px /= wsum; py /= wsum; pz /= wsum; } else { px = x; py = y; pz = z; }
      const s = sr[atomRes[i]] ?? smooth(t);
      out[3 * i] = px + s * residual[3 * i];
      out[3 * i + 1] = py + s * residual[3 * i + 1];
      out[3 * i + 2] = pz + s * residual[3 * i + 2];
    }
    return out;
  }
  if (!names) return { at: rigidAt, rigidAt, windows: win, residueWindows: resWin };

  // ---- refinement: distances between key atoms (C4', C1', P of every residue) interpolated with the residues'
  // progress, solved by sparse stress majorisation from the rigid estimate; atoms then follow their residue.
  const keyNames = new Set(["C4'", "C1'", 'P']);
  const resAtoms = [];
  for (let i = 0; i < nAtoms; i++) (resAtoms[atomRes[i]] ||= []).push(i);
  const key = [];                     // atom index of every key point
  const keyOfRes = [];
  for (let r = 0; r < resAtoms.length; r++) {
    const ks = (resAtoms[r] || []).filter(i => keyNames.has(names[i]));
    const pick = ks.length >= 3 ? ks : (resAtoms[r] || []).slice(0, 3);
    keyOfRes[r] = pick.map(i => { key.push(i); return key.length - 1; });
  }
  const P = key.length;
  const ka = new Float64Array(3 * P), kb = new Float64Array(3 * P);
  key.forEach((i, p) => { for (let c = 0; c < 3; c++) { ka[3 * p + c] = A[3 * i + c]; kb[3 * p + c] = B[3 * i + c]; } });
  const dist = (X, p, q) => { const dx = X[3 * p] - X[3 * q], dy = X[3 * p + 1] - X[3 * q + 1], dz = X[3 * p + 2] - X[3 * q + 2]; return Math.sqrt(dx * dx + dy * dy + dz * dz); };
  const keyRes = Int32Array.from(key.map(i => atomRes[i]));
  // sparse pair list (CSR): key-point pairs closer than 11 A in A or B, plus chain neighbours (weighted up)
  const deg = new Int32Array(P);
  const plist = [];
  for (let p = 0; p < P; p++) {
    for (let q = p + 1; q < P; q++) {
      const da = dist(ka, p, q), db = dist(kb, p, q);
      const chain = Math.abs(keyRes[p] - keyRes[q]) <= 1;
      if (Math.min(da, db) < 11 || chain) { plist.push(p, q, da, db, (chain ? 4 : 1) / Math.max(1, ((da + db) / 2) ** 2)); deg[p]++; deg[q]++; }
    }
  }
  const start = new Int32Array(P + 1);
  for (let p = 0; p < P; p++) start[p + 1] = start[p] + deg[p];
  const nq = new Int32Array(start[P]), nda = new Float64Array(start[P]), ndb = new Float64Array(start[P]), nw = new Float64Array(start[P]);
  const fill = Int32Array.from(start.subarray(0, P));
  for (let k = 0; k < plist.length; k += 5) {
    const p = plist[k], q = plist[k + 1];
    for (const [a, b] of [[p, q], [q, p]]) { const e = fill[a]++; nq[e] = b; nda[e] = plist[k + 2]; ndb[e] = plist[k + 3]; nw[e] = plist[k + 4]; }
  }
  const wsum = new Float64Array(P);
  for (let p = 0; p < P; p++) for (let e = start[p]; e < start[p + 1]; e++) wsum[p] += nw[e];

  function kabsch(Pts, Q, n) {           // rotation (row-major, column convention) and centroids for points 0..n-1
    const pc = [0, 0, 0], qc = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let c = 0; c < 3; c++) { pc[c] += Pts[3 * i + c] / n; qc[c] += Q[3 * i + c] / n; }
    const H = [0, 0, 0, 0, 0, 0, 0, 0, 0];
    for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) H[3 * a + b] += (Pts[3 * i + a] - pc[a]) * (Q[3 * i + b] - qc[b]);
    return { R: rotationFromCovariance(H), pc, qc };
  }
  function applyTo(T, src, i, dst, j) {      // dst[j..j+2] = T(src[i..i+2])
    const R = T.R, x = src[i] - T.pc[0], y = src[i + 1] - T.pc[1], z = src[i + 2] - T.pc[2];
    dst[j] = R[0] * x + R[1] * y + R[2] * z + T.qc[0];
    dst[j + 1] = R[3] * x + R[4] * y + R[5] * z + T.qc[1];
    dst[j + 2] = R[6] * x + R[7] * y + R[8] * z + T.qc[2];
  }
  // per residue: its key points in A and B (for placing its atoms on the solved key points)
  const resKeyA = [], resKeyB = [];
  for (let r = 0; r < resAtoms.length; r++) {
    const ks = keyOfRes[r] || [];
    resKeyA[r] = Float64Array.from(ks.flatMap(p => [ka[3 * p], ka[3 * p + 1], ka[3 * p + 2]]));
    resKeyB[r] = Float64Array.from(ks.flatMap(p => [kb[3 * p], kb[3 * p + 1], kb[3 * p + 2]]));
  }
  const cache = new Map();
  const STEPS = 36, SWEEPS = staged ? 32 : 20;
  const target = new Float64Array(start[P]);
  const tmpA = new Float64Array(3), tmpB = new Float64Array(3);
  function frame(k) {
    if (cache.has(k)) return cache.get(k);
    const t = k / STEPS;
    let out;
    if (k === 0) out = Float32Array.from(A);
    else if (k === STEPS) out = Float32Array.from(B);
    else {
      const guess = rigidAt(t);
      const X = new Float64Array(3 * P);
      for (let p = 0; p < P; p++) { const i = key[p]; X[3 * p] = guess[3 * i]; X[3 * p + 1] = guess[3 * i + 1]; X[3 * p + 2] = guess[3 * i + 2]; }
      const G = Float64Array.from(X);
      const sr = resWin.map(w => progress(w, t));
      for (let p = 0; p < P; p++) {
        const sp = sr[keyRes[p]];
        for (let e = start[p]; e < start[p + 1]; e++) { const s = (sp + sr[keyRes[nq[e]]]) / 2; target[e] = nda[e] + (ndb[e] - nda[e]) * s; }
      }
      for (let sweep = 0; sweep < SWEEPS; sweep++) {
        for (let p = 0; p < P; p++) {
          if (wsum[p] <= 0) continue;
          const x0 = X[3 * p], y0 = X[3 * p + 1], z0 = X[3 * p + 2];
          let wx = 0, wy = 0, wz = 0;
          for (let e = start[p]; e < start[p + 1]; e++) {
            const q = nq[e];
            const qx = X[3 * q], qy = X[3 * q + 1], qz = X[3 * q + 2];
            const dx = x0 - qx, dy = y0 - qy, dz = z0 - qz;
            const d = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6;
            const w = nw[e], f = target[e] / d;
            wx += w * (qx + f * dx); wy += w * (qy + f * dy); wz += w * (qz + f * dz);
          }
          X[3 * p] = wx / wsum[p]; X[3 * p + 1] = wy / wsum[p]; X[3 * p + 2] = wz / wsum[p];
        }
      }
      // keep the display frame of the rigid estimate
      const Tg = kabsch(X, G, P);
      const Xs = new Float64Array(3 * P);
      for (let p = 0; p < P; p++) applyTo(Tg, X, 3 * p, Xs, 3 * p);
      // atoms: blend of the residue's A and B atoms, each placed on the solved key points
      out = new Float32Array(A.length);
      for (let r = 0; r < resAtoms.length; r++) {
        const ks = keyOfRes[r];
        if (!resAtoms[r]) continue;
        const s = sr[r];
        if (!ks || ks.length < 3) {
          for (const i of resAtoms[r]) for (let c = 0; c < 3; c++) out[3 * i + c] = guess[3 * i + c];
          continue;
        }
        const sX = Float64Array.from(ks.flatMap(p => [Xs[3 * p], Xs[3 * p + 1], Xs[3 * p + 2]]));
        const TA = kabsch(resKeyA[r], sX, ks.length), TB = kabsch(resKeyB[r], sX, ks.length);
        for (const i of resAtoms[r]) {
          applyTo(TA, A, 3 * i, tmpA, 0);
          applyTo(TB, B, 3 * i, tmpB, 0);
          for (let c = 0; c < 3; c++) out[3 * i + c] = tmpA[c] + (tmpB[c] - tmpA[c]) * s;
        }
      }
    }
    cache.set(k, out);
    return out;
  }
  function at(t, out = new Float32Array(A.length)) {
    const x = Math.min(1, Math.max(0, t)) * STEPS, k = Math.floor(x), f = x - k;
    const F0 = frame(k), F1 = frame(Math.min(STEPS, k + 1));
    for (let i = 0; i < out.length; i++) out[i] = F0[i] + (F1[i] - F0[i]) * f;
    return out;
  }
  return { at, rigidAt, windows: win, residueWindows: resWin, steps: STEPS };
}

// Rotation matrix (row-major, acting on column vectors: q = R p) maximising tr(R^T H^T) for the covariance
// H = sum p q^T, via the quaternion method of Horn (1987) with power iteration on the 4x4 key matrix.
export function rotationFromCovariance(H) {
  const [Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz] = H;
  const N = [
    [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
    [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
    [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
    [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
  ];
  // largest eigenvector of the symmetric 4x4 key matrix by cyclic Jacobi rotations
  const a = N.map(r => r.slice());
  const V = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
  for (let sweep = 0; sweep < 12; sweep++) {
    let off = 0;
    for (let p = 0; p < 3; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
    if (off < 1e-22) break;
    for (let p = 0; p < 3; p++) {
      for (let q = p + 1; q < 4; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), sn = t * c;
        for (let k = 0; k < 4; k++) {
          const akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - sn * akq; a[k][q] = sn * akp + c * akq;
        }
        for (let k = 0; k < 4; k++) {
          const apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - sn * aqk; a[q][k] = sn * apk + c * aqk;
        }
        for (let k = 0; k < 4; k++) {
          const vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - sn * vkq; V[k][q] = sn * vkp + c * vkq;
        }
      }
    }
  }
  let best = 0;
  for (let k = 1; k < 4; k++) if (a[k][k] > a[best][best]) best = k;
  const v = [V[0][best], V[1][best], V[2][best], V[3][best]];
  const nv = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]) || 1;
  for (let k = 0; k < 4; k++) v[k] /= nv;
  const [w, x, y, z] = v;
  return [
    w * w + x * x - y * y - z * z, 2 * (x * y - w * z), 2 * (x * z + w * y),
    2 * (x * y + w * z), w * w - x * x + y * y - z * z, 2 * (y * z - w * x),
    2 * (x * z - w * y), 2 * (y * z + w * x), w * w - x * x - y * y + z * z,
  ];
}

// Bond-length distortion along a morph: the largest relative change of consecutive C4'-C4' distances (frames f).
export function chainStrain(frames, trace) {
  let worst = 0;
  const d = (X, a, b) => Math.hypot(X[3 * a] - X[3 * b], X[3 * a + 1] - X[3 * b + 1], X[3 * a + 2] - X[3 * b + 2]);
  const ref = [];
  for (let r = 0; r + 1 < trace.length; r++) ref.push(trace[r] >= 0 && trace[r + 1] >= 0 ? (d(frames[0], trace[r], trace[r + 1]) + d(frames[frames.length - 1], trace[r], trace[r + 1])) / 2 : null);
  for (const X of frames) {
    for (let r = 0; r + 1 < trace.length; r++) {
      if (ref[r] === null || ref[r] < 1e-6) continue;
      worst = Math.max(worst, Math.abs(d(X, trace[r], trace[r + 1]) - ref[r]) / ref[r]);
    }
  }
  return worst;
}
