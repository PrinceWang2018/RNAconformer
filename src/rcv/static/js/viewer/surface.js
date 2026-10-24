// Molecular surface: a smooth Gaussian ("blobby") surface around the atoms, meshed with naive surface nets.
//   field f(p) = Σ exp(-B (|p - a|² / R² - 1)), R = van der Waals radius, surface f = 1 (B = 2.2 gives a surface
//   just outside the van der Waals spheres that closes small crevices, like a solvent-excluded surface).
// Every vertex carries the atom that dominates the field there, so the surface takes the residue colours.
const VDW_R = { C: 1.7, N: 1.55, O: 1.52, P: 1.8, S: 1.8 };

export function gaussianSurface(ens, coords, atoms, { spacing = 0.9, blob = 2.2 } = {}) {
  if (!atoms.length) return null;
  const R = atoms.map(a => VDW_R[ens.elements[a]] || 1.7);
  const cut = R.map(r => r * Math.sqrt(1 + Math.log(200) / blob));       // contribution < 0.5 % beyond
  let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
  atoms.forEach((a, i) => {
    const x = coords[3 * a], y = coords[3 * a + 1], z = coords[3 * a + 2], c = cut[i];
    x0 = Math.min(x0, x - c); y0 = Math.min(y0, y - c); z0 = Math.min(z0, z - c);
    x1 = Math.max(x1, x + c); y1 = Math.max(y1, y + c); z1 = Math.max(z1, z + c);
  });
  const h = spacing;
  const nx = Math.ceil((x1 - x0) / h) + 2, ny = Math.ceil((y1 - y0) / h) + 2, nz = Math.ceil((z1 - z0) / h) + 2;
  if (nx * ny * nz > 1.2e7) return gaussianSurface(ens, coords, atoms, { spacing: spacing * 1.35, blob });
  const F = new Float32Array(nx * ny * nz), best = new Float32Array(nx * ny * nz), owner = new Int32Array(nx * ny * nz).fill(-1);
  const idx = (i, j, k) => i + nx * (j + ny * k);
  atoms.forEach((a, n) => {
    const ax = coords[3 * a], ay = coords[3 * a + 1], az = coords[3 * a + 2], r2 = R[n] * R[n], c = cut[n];
    const i0 = Math.max(0, Math.floor((ax - c - x0) / h)), i1 = Math.min(nx - 1, Math.ceil((ax + c - x0) / h));
    const j0 = Math.max(0, Math.floor((ay - c - y0) / h)), j1 = Math.min(ny - 1, Math.ceil((ay + c - y0) / h));
    const k0 = Math.max(0, Math.floor((az - c - z0) / h)), k1 = Math.min(nz - 1, Math.ceil((az + c - z0) / h));
    const c2 = c * c;
    for (let k = k0; k <= k1; k++) {
      const dz = z0 + k * h - az;
      for (let j = j0; j <= j1; j++) {
        const dy = y0 + j * h - ay, dyz = dy * dy + dz * dz;
        if (dyz > c2) continue;
        let p = idx(i0, j, k);
        for (let i = i0; i <= i1; i++, p++) {
          const dx = x0 + i * h - ax, d2 = dx * dx + dyz;
          if (d2 > c2) continue;
          const v = Math.exp(-blob * (d2 / r2 - 1));
          F[p] += v;
          if (v > best[p]) { best[p] = v; owner[p] = a; }
        }
      }
    }
  });
  // ---- surface nets: one vertex per cell crossing the surface, one quad per grid edge crossing it
  const cellV = new Int32Array((nx - 1) * (ny - 1) * (nz - 1)).fill(-1);
  const cidx = (i, j, k) => i + (nx - 1) * (j + (ny - 1) * k);
  const pos = [], own = [];
  const corner = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]];
  const edges = [[0, 1], [2, 3], [4, 5], [6, 7], [0, 2], [1, 3], [4, 6], [5, 7], [0, 4], [1, 5], [2, 6], [3, 7]];
  const val = new Float32Array(8);
  for (let k = 0; k < nz - 1; k++) for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
    let mask = 0;
    for (let c = 0; c < 8; c++) { const v = F[idx(i + corner[c][0], j + corner[c][1], k + corner[c][2])] - 1; val[c] = v; if (v > 0) mask |= 1 << c; }
    if (mask === 0 || mask === 255) continue;
    let sx = 0, sy = 0, sz = 0, cnt = 0;
    for (const [a, b] of edges) {
      if ((val[a] > 0) === (val[b] > 0)) continue;
      const t = val[a] / (val[a] - val[b]);
      sx += corner[a][0] + (corner[b][0] - corner[a][0]) * t;
      sy += corner[a][1] + (corner[b][1] - corner[a][1]) * t;
      sz += corner[a][2] + (corner[b][2] - corner[a][2]) * t;
      cnt++;
    }
    cellV[cidx(i, j, k)] = pos.length / 3;
    pos.push(x0 + (i + sx / cnt) * h, y0 + (j + sy / cnt) * h, z0 + (k + sz / cnt) * h);
    let bo = -1, bv = -Infinity;
    for (let c = 0; c < 8; c++) { const p = idx(i + corner[c][0], j + corner[c][1], k + corner[c][2]); if (F[p] > bv && owner[p] >= 0) { bv = F[p]; bo = owner[p]; } }
    own.push(bo);
  }
  const tri = [];
  const quad = (a, b, c, d, flip) => { if (a < 0 || b < 0 || c < 0 || d < 0) return; if (flip) tri.push(a, c, b, a, d, c); else tri.push(a, b, c, a, c, d); };
  for (let k = 1; k < nz - 1; k++) for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
    const inside = F[idx(i, j, k)] > 1;
    if (i < nx - 1 && (F[idx(i + 1, j, k)] > 1) !== inside)
      quad(cellV[cidx(i, j - 1, k - 1)], cellV[cidx(i, j, k - 1)], cellV[cidx(i, j, k)], cellV[cidx(i, j - 1, k)], !inside);
    if (j < ny - 1 && (F[idx(i, j + 1, k)] > 1) !== inside)
      quad(cellV[cidx(i - 1, j, k - 1)], cellV[cidx(i - 1, j, k)], cellV[cidx(i, j, k)], cellV[cidx(i, j, k - 1)], !inside);
    if (k < nz - 1 && (F[idx(i, j, k + 1)] > 1) !== inside)
      quad(cellV[cidx(i - 1, j - 1, k)], cellV[cidx(i, j - 1, k)], cellV[cidx(i, j, k)], cellV[cidx(i - 1, j, k)], !inside);
  }
  // ---- normals: gradient of the field (trilinear), pointing outwards
  const sample = (x, y, z) => {
    const fx = (x - x0) / h, fy = (y - y0) / h, fz = (z - z0) / h;
    const i = Math.min(nx - 2, Math.max(0, Math.floor(fx))), j = Math.min(ny - 2, Math.max(0, Math.floor(fy))), k = Math.min(nz - 2, Math.max(0, Math.floor(fz)));
    const u = fx - i, v = fy - j, w = fz - k;
    const c00 = F[idx(i, j, k)] * (1 - u) + F[idx(i + 1, j, k)] * u, c10 = F[idx(i, j + 1, k)] * (1 - u) + F[idx(i + 1, j + 1, k)] * u;
    const c01 = F[idx(i, j, k + 1)] * (1 - u) + F[idx(i + 1, j, k + 1)] * u, c11 = F[idx(i, j + 1, k + 1)] * (1 - u) + F[idx(i + 1, j + 1, k + 1)] * u;
    return (c00 * (1 - v) + c10 * v) * (1 - w) + (c01 * (1 - v) + c11 * v) * w;
  };
  const nv = pos.length / 3, nrm = new Float32Array(pos.length), e = h * 0.5;
  for (let q = 0; q < nv; q++) {
    const x = pos[3 * q], y = pos[3 * q + 1], z = pos[3 * q + 2];
    let gx = sample(x - e, y, z) - sample(x + e, y, z), gy = sample(x, y - e, z) - sample(x, y + e, z), gz = sample(x, y, z - e) - sample(x, y, z + e);
    const l = Math.hypot(gx, gy, gz) || 1;
    nrm[3 * q] = gx / l; nrm[3 * q + 1] = gy / l; nrm[3 * q + 2] = gz / l;
  }
  return { positions: Float32Array.from(pos), normals: nrm, atomOf: Int32Array.from(own), index: nv > 65535 ? Uint32Array.from(tri) : Uint16Array.from(tri), area: tri.length / 3 };
}
