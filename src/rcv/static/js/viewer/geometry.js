// Geometry builders for RNA representations. Every builder keeps its topology fixed and rewrites
// positions/normals in place, so conformations can be morphed at interactive frame rates.
import * as THREE from '../../vendor/three/three.module.min.js';

export const VDW = { C: 1.7, N: 1.55, O: 1.52, P: 1.8, S: 1.8 };

// ------------------------------------------------------------------ backbone tube (C4' spline)
export class Tube {
  constructor(ens, { radius = 1.0, radial = 14, sub = 8, maxGap = 9.5 } = {}) {
    this.ens = ens;
    this.radial = radial;
    this.sub = sub;
    this.maxGap = maxGap;
    this.radius = radius;
    this.radii = null;                       // optional per-residue radius (putty)
    this.capSteps = 4;
    this.geometry = new THREE.BufferGeometry();
    this.segments = null;
  }

  // Split residues with a trace atom into continuous segments (chain breaks > maxGap).
  _segments(coords) {
    const ens = this.ens;
    const segs = [];
    let cur = [];
    for (let r = 0; r < ens.nRes; r++) {
      const a = ens.trace[r];
      if (a < 0) { if (cur.length) segs.push(cur); cur = []; continue; }
      if (cur.length) {
        const b = ens.trace[cur[cur.length - 1]];
        const d = Math.hypot(coords[3 * a] - coords[3 * b], coords[3 * a + 1] - coords[3 * b + 1], coords[3 * a + 2] - coords[3 * b + 2]);
        if (d > this.maxGap) { segs.push(cur); cur = []; }
      }
      cur.push(r);
    }
    if (cur.length) segs.push(cur);
    return segs.filter(s => s.length >= 2);
  }

  _layout(segs) {
    // samples per segment: (n-1)*sub+1 along spline + capSteps at each end
    const R = this.radial;
    let nv = 0, ni = 0;
    this.layout = segs.map(s => {
      const ns = (s.length - 1) * this.sub + 1 + 2 * this.capSteps;
      const off = nv;
      nv += ns * R + 2;                       // + 2 pole vertices
      ni += (ns - 1) * R * 6 + 2 * R * 3;
      return { res: s, ns, off };
    });
    const pos = new Float32Array(nv * 3), nor = new Float32Array(nv * 3), col = new Float32Array(nv * 3);
    const idx = new Uint32Array(ni);
    let k = 0;
    for (const L of this.layout) {
      const { ns, off } = L;
      for (let i = 0; i < ns - 1; i++) for (let j = 0; j < R; j++) {
        const a = off + i * R + j, b = off + i * R + (j + 1) % R, c = a + R, d = b + R;
        idx[k++] = a; idx[k++] = b; idx[k++] = c;
        idx[k++] = b; idx[k++] = d; idx[k++] = c;
      }
      const p0 = off + ns * R, p1 = p0 + 1;
      for (let j = 0; j < R; j++) {
        idx[k++] = p0; idx[k++] = off + j; idx[k++] = off + (j + 1) % R;
        const last = off + (ns - 1) * R;
        idx[k++] = p1; idx[k++] = last + (j + 1) % R; idx[k++] = last + j;
      }
    }
    this.geometry.setIndex(new THREE.BufferAttribute(idx, 1));
    this.geometry.setAttribute('position', new THREE.BufferAttribute(pos, 3).setUsage(THREE.DynamicDrawUsage));
    this.geometry.setAttribute('normal', new THREE.BufferAttribute(nor, 3).setUsage(THREE.DynamicDrawUsage));
    this.geometry.setAttribute('color', new THREE.BufferAttribute(col, 3).setUsage(THREE.DynamicDrawUsage));
  }

  update(coords) {
    const segs = this._segments(coords);
    const key = segs.map(s => s[0] + ':' + s.length).join(',');
    if (key !== this.segKey) { this.segKey = key; this._layout(segs); }
    const pos = this.geometry.attributes.position.array, nor = this.geometry.attributes.normal.array;
    const R = this.radial, S = this.sub, C = this.capSteps;
    const ens = this.ens;
    for (const L of this.layout) {
      const n = L.res.length;
      const P = L.res.map(r => { const a = ens.trace[r]; return [coords[3 * a], coords[3 * a + 1], coords[3 * a + 2]]; });
      const pts = [], tans = [], rad = [];
      const radAt = t => {
        if (!this.radii) return this.radius;
        const i = Math.min(n - 2, Math.floor(t)), f = t - i;
        return this.radii[L.res[i]] * (1 - f) + this.radii[L.res[i + 1]] * f;
      };
      for (let i = 0; i < n - 1; i++) {
        const p0 = P[Math.max(0, i - 1)], p1 = P[i], p2 = P[i + 1], p3 = P[Math.min(n - 1, i + 2)];
        for (let s = 0; s < S || (i === n - 2 && s === S); s++) {
          const t = s / S;
          const [p, d] = catmull(p0, p1, p2, p3, t);
          pts.push(p); tans.push(d); rad.push(radAt(i + t));
        }
      }
      // rotation-minimising frames (double reflection)
      const ns = pts.length;
      const T = tans.map(normalize);
      const N = new Array(ns), B = new Array(ns);
      N[0] = normalize(perp(T[0]));
      for (let i = 0; i < ns - 1; i++) {
        const v1 = sub(pts[i + 1], pts[i]);
        const c1 = dot(v1, v1) || 1e-9;
        const rL = sub(N[i], scale(v1, 2 * dot(v1, N[i]) / c1));
        const tL = sub(T[i], scale(v1, 2 * dot(v1, T[i]) / c1));
        const v2 = sub(T[i + 1], tL);
        const c2 = dot(v2, v2) || 1e-9;
        N[i + 1] = normalize(sub(rL, scale(v2, 2 * dot(v2, rL) / c2)));
      }
      for (let i = 0; i < ns; i++) B[i] = cross(T[i], N[i]);
      // rounded caps: extra rings before the first and after the last sample
      const rings = [];
      for (let c = C; c >= 1; c--) {
        const phi = (c / C) * Math.PI / 2;
        rings.push({ p: add(pts[0], scale(T[0], -rad[0] * Math.sin(phi))), r: rad[0] * Math.cos(phi), T: T[0], N: N[0], B: B[0], nz: -Math.sin(phi) });
      }
      for (let i = 0; i < ns; i++) rings.push({ p: pts[i], r: rad[i], T: T[i], N: N[i], B: B[i], nz: 0 });
      for (let c = 1; c <= C; c++) {
        const phi = (c / C) * Math.PI / 2, e = ns - 1;
        rings.push({ p: add(pts[e], scale(T[e], rad[e] * Math.sin(phi))), r: rad[e] * Math.cos(phi), T: T[e], N: N[e], B: B[e], nz: Math.sin(phi) });
      }
      let v = L.off;
      for (const ring of rings) {
        const cz = Math.sqrt(Math.max(0, 1 - ring.nz * ring.nz));
        for (let j = 0; j < R; j++) {
          const a = (j / R) * Math.PI * 2, ca = Math.cos(a), sa = Math.sin(a);
          const dx = ring.N[0] * ca + ring.B[0] * sa, dy = ring.N[1] * ca + ring.B[1] * sa, dz = ring.N[2] * ca + ring.B[2] * sa;
          pos[3 * v] = ring.p[0] + ring.r * dx; pos[3 * v + 1] = ring.p[1] + ring.r * dy; pos[3 * v + 2] = ring.p[2] + ring.r * dz;
          nor[3 * v] = dx * cz + ring.T[0] * ring.nz; nor[3 * v + 1] = dy * cz + ring.T[1] * ring.nz; nor[3 * v + 2] = dz * cz + ring.T[2] * ring.nz;
          v++;
        }
      }
      const first = rings[0], last = rings[rings.length - 1];
      const r0 = rad[0], r1 = rad[ns - 1];
      const pole0 = add(pts[0], scale(T[0], -r0)), pole1 = add(pts[ns - 1], scale(T[ns - 1], r1));
      for (const [pp, tt] of [[pole0, scale(first.T, -1)], [pole1, last.T]]) {
        pos[3 * v] = pp[0]; pos[3 * v + 1] = pp[1]; pos[3 * v + 2] = pp[2];
        nor[3 * v] = tt[0]; nor[3 * v + 1] = tt[1]; nor[3 * v + 2] = tt[2];
        v++;
      }
      L.ringCount = rings.length;
    }
    this.geometry.attributes.position.needsUpdate = true;
    this.geometry.attributes.normal.needsUpdate = true;
    this.geometry.computeBoundingSphere();
    if (this._colors) this.setColors(this._colors);
  }

  // colors: Float32Array nRes*3 (linear RGB)
  setColors(colors) {
    this._colors = colors;
    if (!this.layout) return;
    const col = this.geometry.attributes.color.array;
    const R = this.radial, S = this.sub, C = this.capSteps;
    for (const L of this.layout) {
      const n = L.res.length;
      let v = L.off;
      const total = L.ringCount ?? ((n - 1) * S + 1 + 2 * C);
      for (let ri = 0; ri < total; ri++) {
        const sIdx = Math.min(Math.max(ri - C, 0), (n - 1) * S);
        const r = L.res[Math.min(n - 1, Math.round(sIdx / S))];
        for (let j = 0; j < R; j++, v++) { col[3 * v] = colors[3 * r]; col[3 * v + 1] = colors[3 * r + 1]; col[3 * v + 2] = colors[3 * r + 2]; }
      }
      const r0 = L.res[0], r1 = L.res[n - 1];
      for (const r of [r0, r1]) { col[3 * v] = colors[3 * r]; col[3 * v + 1] = colors[3 * r + 1]; col[3 * v + 2] = colors[3 * r + 2]; v++; }
    }
    this.geometry.attributes.color.needsUpdate = true;
  }
}

function catmull(p0, p1, p2, p3, t) {
  // centripetal Catmull-Rom (alpha = 0.5) via Barry-Goldman, with numerical derivative
  const f = tt => cr(p0, p1, p2, p3, tt);
  const p = f(t);
  const e = 1e-3;
  const a = f(Math.max(0, t - e)), b = f(Math.min(1, t + e));
  return [p, sub(b, a)];
}
function cr(p0, p1, p2, p3, t) {
  const d01 = Math.max(1e-4, Math.pow(dist(p0, p1), 0.5)), d12 = Math.max(1e-4, Math.pow(dist(p1, p2), 0.5)), d23 = Math.max(1e-4, Math.pow(dist(p2, p3), 0.5));
  const t0 = 0, t1 = d01, t2 = t1 + d12, t3 = t2 + d23;
  const tt = t1 + (t2 - t1) * t;
  const L = (a, b, ta, tb) => { const w = tb - ta < 1e-9 ? 0 : (tt - ta) / (tb - ta); return [a[0] + (b[0] - a[0]) * w, a[1] + (b[1] - a[1]) * w, a[2] + (b[2] - a[2]) * w]; };
  const A1 = L(p0, p1, t0, t1), A2 = L(p1, p2, t1, t2), A3 = L(p2, p3, t2, t3);
  const B1 = L(A1, A2, t0, t2), B2 = L(A2, A3, t1, t3);
  return L(B1, B2, t1, t2);
}
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const scale = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
function normalize(a) { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; }
function perp(t) { const a = Math.abs(t[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]; return cross(t, a); }

// ------------------------------------------------------------------ base slabs (filled rings)
export class Slabs {
  constructor(ens, { half = 0.32 } = {}) {
    this.ens = ens;
    this.half = half;
    this.geometry = new THREE.BufferGeometry();
    // vertex budget per residue: fans (top+bottom) + side quads
    this.items = [];
    let nv = 0;
    for (let r = 0; r < ens.nRes; r++) {
      const ring = ens.rings[r];
      if (!ring) continue;
      const fanTris = ring.fans.reduce((s, f) => s + f.length, 0);
      const count = fanTris * 3 * 2 + ring.outline.length * 6;
      this.items.push({ r, ring, off: nv, count });
      nv += count;
    }
    this.nv = nv;
    this.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(nv * 3), 3).setUsage(THREE.DynamicDrawUsage));
    this.geometry.setAttribute('normal', new THREE.BufferAttribute(new Float32Array(nv * 3), 3).setUsage(THREE.DynamicDrawUsage));
    this.geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(nv * 3), 3).setUsage(THREE.DynamicDrawUsage));
  }

  update(c) {
    const pos = this.geometry.attributes.position.array, nor = this.geometry.attributes.normal.array;
    const h = this.half;
    const P = i => [c[3 * i], c[3 * i + 1], c[3 * i + 2]];
    for (const it of this.items) {
      if (this.hidden && this.hidden.has(it.r)) {          // hidden residue: collapse its slab to a point
        const a = it.ring.outline[0];
        for (let v = it.off; v < it.off + it.count; v++) { pos[3 * v] = c[3 * a]; pos[3 * v + 1] = c[3 * a + 1]; pos[3 * v + 2] = c[3 * a + 2]; }
        continue;
      }
      const out = it.ring.outline.map(P);
      // Newell normal of the outline
      let n = [0, 0, 0];
      for (let i = 0; i < out.length; i++) {
        const a = out[i], b = out[(i + 1) % out.length];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]); n[1] += (a[2] - b[2]) * (a[0] + b[0]); n[2] += (a[0] - b[0]) * (a[1] + b[1]);
      }
      n = normalize(n);
      const up = scale(n, h), dn = scale(n, -h);
      let v = it.off;
      const put = (p, nn) => { pos[3 * v] = p[0]; pos[3 * v + 1] = p[1]; pos[3 * v + 2] = p[2]; nor[3 * v] = nn[0]; nor[3 * v + 1] = nn[1]; nor[3 * v + 2] = nn[2]; v++; };
      for (const fan of it.ring.fans) {
        const pts = fan.map(P);
        const cen = scale(pts.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / pts.length);
        for (let i = 0; i < pts.length; i++) {
          const a = pts[i], b = pts[(i + 1) % pts.length];
          put(add(cen, up), n); put(add(a, up), n); put(add(b, up), n);
          const nn = scale(n, -1);
          put(add(cen, dn), nn); put(add(b, dn), nn); put(add(a, dn), nn);
        }
      }
      const cen = scale(out.reduce((s, p) => add(s, p), [0, 0, 0]), 1 / out.length);
      for (let i = 0; i < out.length; i++) {
        const a = out[i], b = out[(i + 1) % out.length];
        let sn = normalize(cross(sub(b, a), n));
        const mid = scale(add(a, b), 0.5);
        if (dot(sn, sub(mid, cen)) < 0) sn = scale(sn, -1);
        const a1 = add(a, up), b1 = add(b, up), a0 = add(a, dn), b0 = add(b, dn);
        put(a0, sn); put(b0, sn); put(b1, sn);
        put(a0, sn); put(b1, sn); put(a1, sn);
      }
    }
    this.geometry.attributes.position.needsUpdate = true;
    this.geometry.attributes.normal.needsUpdate = true;
    this.geometry.computeBoundingSphere();
  }

  setColors(colors) {
    const col = this.geometry.attributes.color.array;
    for (const it of this.items) {
      const r = it.r;
      for (let v = it.off; v < it.off + it.count; v++) { col[3 * v] = colors[3 * r]; col[3 * v + 1] = colors[3 * r + 1]; col[3 * v + 2] = colors[3 * r + 2]; }
    }
    this.geometry.attributes.color.needsUpdate = true;
  }
}

// ------------------------------------------------------------------ instanced cylinders & spheres
const _m = new THREE.Matrix4();
const _q = new THREE.Quaternion();
const _v = new THREE.Vector3();
const _s = new THREE.Vector3();
const _up = new THREE.Vector3(0, 1, 0);

export class Cylinders {
  // pairs: Int32Array [a0,b0,a1,b1...] of point indices; each pair is drawn as one cylinder
  constructor(capacity, material, radial = 12) {
    const g = new THREE.CylinderGeometry(1, 1, 1, radial, 1, true);
    this.mesh = new THREE.InstancedMesh(g, material, Math.max(1, capacity));
    this.mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.mesh.frustumCulled = false;
    this.mesh.count = 0;
    this.capacity = capacity;
  }
  // segs: array-like of [x0,y0,z0,x1,y1,z1] produced by callback for i in [0,count)
  set(count, getSeg, radius) {
    const mesh = this.mesh;
    mesh.count = Math.min(count, this.capacity);
    const arr = mesh.instanceMatrix.array;
    const seg = new Float32Array(6);
    for (let i = 0; i < mesh.count; i++) {
      getSeg(i, seg);
      const dx = seg[3] - seg[0], dy = seg[4] - seg[1], dz = seg[5] - seg[2];
      const len = Math.hypot(dx, dy, dz) || 1e-6;
      _v.set(dx / len, dy / len, dz / len);
      _q.setFromUnitVectors(_up, _v);
      _s.set(radius, len, radius);
      _m.compose(_v.set((seg[0] + seg[3]) / 2, (seg[1] + seg[4]) / 2, (seg[2] + seg[5]) / 2), _q, _s);
      _m.toArray(arr, i * 16);
    }
    mesh.instanceMatrix.needsUpdate = true;
  }
  setColors(count, getColor) {
    const mesh = this.mesh;
    if (!mesh.instanceColor) mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(this.capacity * 3), 3);
    const a = mesh.instanceColor.array;
    const tmp = new Float32Array(3);
    for (let i = 0; i < Math.min(count, this.capacity); i++) { getColor(i, tmp); a[3 * i] = tmp[0]; a[3 * i + 1] = tmp[1]; a[3 * i + 2] = tmp[2]; }
    mesh.instanceColor.needsUpdate = true;
  }
}

export class Spheres {
  constructor(capacity, material, detail = 2) {
    const g = new THREE.IcosahedronGeometry(1, detail);
    this.mesh = new THREE.InstancedMesh(g, material, Math.max(1, capacity));
    this.mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.mesh.frustumCulled = false;
    this.mesh.count = 0;
    this.capacity = capacity;
  }
  set(count, getPos, getRadius) {
    const mesh = this.mesh;
    mesh.count = Math.min(count, this.capacity);
    const arr = mesh.instanceMatrix.array;
    const p = new Float32Array(3);
    for (let i = 0; i < mesh.count; i++) {
      getPos(i, p);
      const r = getRadius(i);
      const o = i * 16;
      arr[o] = r; arr[o + 1] = 0; arr[o + 2] = 0; arr[o + 3] = 0;
      arr[o + 4] = 0; arr[o + 5] = r; arr[o + 6] = 0; arr[o + 7] = 0;
      arr[o + 8] = 0; arr[o + 9] = 0; arr[o + 10] = r; arr[o + 11] = 0;
      arr[o + 12] = p[0]; arr[o + 13] = p[1]; arr[o + 14] = p[2]; arr[o + 15] = 1;
    }
    mesh.instanceMatrix.needsUpdate = true;
  }
  setColors(count, getColor) {
    const mesh = this.mesh;
    if (!mesh.instanceColor) mesh.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(this.capacity * 3), 3);
    const a = mesh.instanceColor.array;
    const tmp = new Float32Array(3);
    for (let i = 0; i < Math.min(count, this.capacity); i++) { getColor(i, tmp); a[3 * i] = tmp[0]; a[3 * i + 1] = tmp[1]; a[3 * i + 2] = tmp[2]; }
    mesh.instanceColor.needsUpdate = true;
  }
}
