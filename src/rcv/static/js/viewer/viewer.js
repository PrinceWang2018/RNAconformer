// 3D viewer: scene, representations, styles, picking, labels/overlays and offscreen export.
import * as THREE from '../../vendor/three/three.module.min.js';
import { TrackballControls } from './controls.js';
import { PostFX } from './post.js';
import { Tube, Slabs, Cylinders, Spheres, VDW } from './geometry.js';
import { ELEMENT } from '../colors.js';
import { svg } from '../util.js';

export const STYLES = {
  studio: { label: 'Studio', material: 'standard', roughness: 0.55, ao: 0.85, aoRadius: 6.0, aoIntensity: 1.6, outline: 0.85, outlineWidth: 1.4, fog: 0 },
  illustrative: { label: 'Illustrative', material: 'flat', ao: 0.9, aoRadius: 6.5, aoIntensity: 1.5, outline: 1.0, outlineWidth: 2.2, fog: 0 },
  glossy: { label: 'Glossy', material: 'physical', roughness: 0.28, ao: 0.6, aoRadius: 5, aoIntensity: 1.0, outline: 0.0, outlineWidth: 1.2, fog: 0 },
  depth: { label: 'Depth cue', material: 'standard', roughness: 0.6, ao: 0.8, aoRadius: 5.5, aoIntensity: 1.2, outline: 0.7, outlineWidth: 1.3, fog: 0.8 },
  plain: { label: 'Plain (fast)', material: 'lambert', ao: 0, aoRadius: 5, aoIntensity: 1, outline: 0, outlineWidth: 1, fog: 0 },
};

function makeMaterial(style, opts = {}) {
  const base = { vertexColors: !!opts.vertexColors, side: opts.side || THREE.FrontSide };
  if (opts.ghost) return new THREE.MeshLambertMaterial({ color: opts.color || 0x8a8f98, transparent: true, opacity: opts.opacity ?? 0.22, depthWrite: false });
  switch (style.material) {
    case 'flat': return new THREE.MeshBasicMaterial(base);
    case 'physical': return new THREE.MeshPhysicalMaterial({ ...base, roughness: style.roughness, metalness: 0, clearcoat: 0.65, clearcoatRoughness: 0.18 });
    case 'lambert': return new THREE.MeshLambertMaterial(base);
    default: return new THREE.MeshStandardMaterial({ ...base, roughness: style.roughness ?? 0.55, metalness: 0 });
  }
}

const REP = {
  'cartoon': { tube: true, slabs: true, glyco: true, atoms: false },
  'cartoon+sticks': { tube: true, slabs: false, glyco: false, atoms: 'thin' },
  'sticks': { tube: false, slabs: false, glyco: false, atoms: 'sticks' },
  'ballstick': { tube: false, slabs: false, glyco: false, atoms: 'ballstick' },
  'spacefill': { tube: false, slabs: false, glyco: false, atoms: 'spacefill' },
};

const _tmp = [0, 0, 0];

export class Viewer {
  constructor({ canvas, container, overlay, labels, onPick, onHover, onDblPick }) {
    this.canvas = canvas;
    this.container = container;
    this.overlayEl = overlay;
    this.labelsEl = labels;
    this.onPick = onPick || (() => {});
    this.onHover = onHover || (() => {});
    this.onDblPick = onDblPick || (() => {});
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: false, alpha: false, preserveDrawingBuffer: true, powerPreference: 'high-performance' });
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    this.renderer.toneMapping = THREE.NoToneMapping;
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.post = new PostFX(this.renderer);
    this.styles = STYLES;                    // style presets (exposed for preferences and plugins)
    this.scene = new THREE.Scene();
    this.pivot = new THREE.Group();          // rotation (controls) + pan
    this.inner = new THREE.Group();          // translation so rotation is about the ensemble centre
    this.pivot.add(this.inner);
    this.scene.add(this.pivot);
    this._lights();
    this.persp = new THREE.PerspectiveCamera(26, 1, 1, 1000);
    this.ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.1, 1000);
    this.projection = 'perspective';
    this.settings = { style: 'studio', background: '#ffffff', representation: 'cartoon', outlineColor: '#141414', quality: 'high', elementColors: true,
      alpha: { cartoon: 1, bases: 1, atoms: 1, surface: 1, overlay: 0.2 } };
    this.ens = null;
    this.coords = null;
    this.residueHex = null;
    this.selection = new Set();
    this.highlight = new Set();
    this.hidden = new Set();               // residues whose main representation is hidden (console: hide)
    this.extraSticks = new Set();          // residues always drawn as sticks (console: show sticks)
    this.extraSpheres = new Set();         // residues drawn as van der Waals spheres (console: show spheres)
    this.extraLabels = new Map();          // residue -> custom label text (null = residue name)
    this.labelMode = 'selection';
    this.overlays = new Map();
    this.renderHooks = new Set();          // fn(project, w, h) after every frame (annotation layer)
    this.animators = new Set();
    this.needsRender = true;
    this.controls = new TrackballControls(canvas, {
      onChange: () => { this.needsRender = true; },
      onStart: () => container.classList.add('rotating'),
      onEnd: () => container.classList.remove('rotating'),
    });
    this._bindPicking();
    this._resizeObs = new ResizeObserver(() => this.resize());
    this._resizeObs.observe(container);
    this.resize();
    this._last = performance.now();
    this.frameStats = { frames: 0, t0: performance.now(), fps: 0 };
    requestAnimationFrame(t => this._loop(t));
  }

  _lights() {
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0xb4afa5, 1.55));
    const key = new THREE.DirectionalLight(0xffffff, 2.3);
    key.position.set(-1.0, 1.35, 2.2);
    const fill = new THREE.DirectionalLight(0xffffff, 0.55);
    fill.position.set(1.6, -0.4, 1.0);
    const rim = new THREE.DirectionalLight(0xffffff, 0.35);
    rim.position.set(0.2, 0.8, -2.0);
    this.scene.add(key, fill, rim);
  }

  get camera() { return this.projection === 'orthographic' ? this.ortho : this.persp; }

  resize() {
    const w = Math.max(1, this.container.clientWidth), h = Math.max(1, this.container.clientHeight);
    this.renderer.setSize(w, h, false);
    const pr = this.renderer.getPixelRatio();
    this.post.setSize(w * pr, h * pr);
    this.cssSize = [w, h];
    this.needsRender = true;
  }

  // ------------------------------------------------------------------ ensemble & representations
  setEnsemble(ens, coords, { keepView = false } = {}) {
    this._disposeReps();
    this.ens = ens;
    this.coords = coords;
    this.selection.clear();
    this.highlight.clear();
    if (!ens) { this.needsRender = true; return; }
    if (!keepView) this.inner.position.set(-ens.center[0], -ens.center[1], -ens.center[2]);
    this._buildReps();
    if (!keepView) this.resetView(false);
    this.needsRender = true;
  }

  _materials() {
    const st = STYLES[this.settings.style] || STYLES.studio;
    this._mats?.forEach(m => m.dispose());
    // one material per kind of object, so each can have its own transparency
    this.matVC = makeMaterial(st, { vertexColors: true });          // cartoon tube
    this.matSlab = makeMaterial(st, { vertexColors: true });        // base slabs
    this.matGlyco = makeMaterial(st, {});                           // glycosidic sticks
    this.matInst = makeMaterial(st, {});                            // atoms and bonds
    this.matSurf = makeMaterial(st, { vertexColors: true, side: THREE.DoubleSide });
    this._mats = [this.matVC, this.matSlab, this.matGlyco, this.matInst, this.matSurf];
    this._applyAlpha();
  }

  _buildReps() {
    const ens = this.ens;
    this._materials();
    this.tube = new Tube(ens, { radius: 1.0 });
    this.tubeMesh = new THREE.Mesh(this.tube.geometry, this.matVC);
    this.tubeMesh.frustumCulled = false;
    this.slabs = new Slabs(ens);
    this.slabMesh = new THREE.Mesh(this.slabs.geometry, this.matSlab);
    this.slabMesh.frustumCulled = false;
    this.glyco = new Cylinders(ens.nRes, this.matGlyco, 10);
    this.atoms = new Spheres(ens.nAtoms, this.matInst, 2);
    this.bondCyl = new Cylinders(ens.bonds.length, this.matInst, 10);       // two half-bonds per bond
    this.selAtoms = new Spheres(ens.nAtoms, this.matInst, 2);
    this.selBonds = new Cylinders(ens.bonds.length, this.matInst, 10);
    this.vdwAtoms = new Spheres(ens.nAtoms, this.matInst, 3);
    this.slabs.hidden = this.hidden;
    for (const m of [this.tubeMesh, this.slabMesh, this.glyco.mesh, this.atoms.mesh, this.bondCyl.mesh, this.selAtoms.mesh, this.selBonds.mesh, this.vdwAtoms.mesh]) this.inner.add(m);
    this.ghostGroup = new THREE.Group();
    this.inner.add(this.ghostGroup);
    this.previewTube = null;
    this._applyRepVisibility();
    this.update();
  }

  _disposeReps() {
    this.setSurface(null);
    for (const m of [this.tubeMesh, this.slabMesh, this.glyco?.mesh, this.atoms?.mesh, this.bondCyl?.mesh, this.selAtoms?.mesh, this.selBonds?.mesh, this.vdwAtoms?.mesh]) {
      if (!m) continue;
      this.inner.remove(m);
      m.geometry.dispose();
    }
    this._clearGhosts();
    this._clearPreview();
    if (this.ghostGroup) this.inner.remove(this.ghostGroup);
    this.tube = this.slabs = this.glyco = this.atoms = this.bondCyl = this.selAtoms = this.selBonds = this.vdwAtoms = null;
    this.tubeMesh = this.slabMesh = null;
  }

  setStyle(style) {
    if (!STYLES[style]) return;
    this.settings.style = style;
    if (!this.ens) return;
    this._materials();
    this.tubeMesh.material = this.matVC;
    this.slabMesh.material = this.matSlab;
    this.glyco.mesh.material = this.matGlyco;
    for (const c of [this.atoms, this.bondCyl, this.selAtoms, this.selBonds, this.vdwAtoms]) c.mesh.material = this.matInst;
    if (this.surfaceMesh) this.surfaceMesh.material = this.matSurf;
    this.needsRender = true;
  }

  setRepresentation(rep) {
    if (!REP[rep]) return;
    this.settings.representation = rep;
    if (!this.ens) return;
    this._applyRepVisibility();
    this.update();
  }

  _applyRepVisibility() {
    const R = REP[this.settings.representation];
    this.tubeMesh.visible = R.tube;
    this.slabMesh.visible = R.slabs;
    this.glyco.mesh.visible = R.glyco;
    this.atoms.mesh.visible = !!R.atoms;
    this.bondCyl.mesh.visible = !!R.atoms && R.atoms !== 'spacefill';
    this.needsRender = true;
  }

  // Per-object opacity (1 = opaque): cartoon, bases, atoms, surface, overlay. Translucent objects do not write depth, so
  // ambient occlusion and outlines still come from the opaque objects behind them.
  setTransparency(key, opacity) {
    this.settings.alpha = { cartoon: 1, bases: 1, atoms: 1, surface: 1, overlay: 0.2, ...(this.settings.alpha || {}) };
    this.settings.alpha[key] = Math.max(0.05, Math.min(1, opacity));
    this._applyAlpha();
    this.needsRender = true;
  }
  _applyAlpha() {
    const A = { cartoon: 1, bases: 1, atoms: 1, surface: 1, overlay: 0.2, ...(this.settings.alpha || {}) };
    const set = (m, a) => {
      if (!m) return;
      const tr = a < 0.999;
      m.transparent = tr; m.opacity = a; m.depthWrite = !tr; m.needsUpdate = true;
    };
    set(this.matVC, A.cartoon); set(this.matSlab, A.bases); set(this.matGlyco, A.bases); set(this.matInst, A.atoms); set(this.matSurf, A.surface);
    if (this.ghostMat) { this.ghostMat.opacity = A.overlay; this.ghostMat.needsUpdate = true; }
    for (const [mesh, a] of [[this.tubeMesh, A.cartoon], [this.slabMesh, A.bases], [this.surfaceMesh, A.surface]]) if (mesh) mesh.renderOrder = a < 0.999 ? 5 : 0;
  }

  // Molecular surface from viewer/surface.js data (null removes it).
  setSurface(data) {
    if (this.surfaceMesh) { this.inner.remove(this.surfaceMesh); this.surfaceMesh.geometry.dispose(); this.surfaceMesh = null; }
    this._surfaceData = data;
    if (!data || !this.ens) { this.needsRender = true; return; }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(data.positions, 3));
    g.setAttribute('normal', new THREE.BufferAttribute(data.normals, 3));
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(data.positions.length), 3));
    g.setIndex(new THREE.BufferAttribute(data.index, 1));
    this.surfaceMesh = new THREE.Mesh(g, this.matSurf);
    this.surfaceMesh.frustumCulled = false;
    this.inner.add(this.surfaceMesh);
    this._applyAlpha();
    this._surfaceColors();
    this.needsRender = true;
  }
  setSurfaceVisible(on) { if (this.surfaceMesh) { this.surfaceMesh.visible = on; this.needsRender = true; } }
  _surfaceColors() {
    const m = this.surfaceMesh, d = this._surfaceData, lin = this.resLinear;
    if (!m || !d || !lin) return;
    const col = m.geometry.attributes.color.array, ens = this.ens;
    for (let v = 0; v < d.atomOf.length; v++) {
      const a = d.atomOf[v];
      if (a < 0) continue;
      if (this.surfaceAtomColors && this._atomColor) { this._atomColor(a, _tmp); col[3 * v] = _tmp[0]; col[3 * v + 1] = _tmp[1]; col[3 * v + 2] = _tmp[2]; continue; }
      const r = ens.atomRes[a];
      col[3 * v] = lin[3 * r]; col[3 * v + 1] = lin[3 * r + 1]; col[3 * v + 2] = lin[3 * r + 2];
    }
    m.geometry.attributes.color.needsUpdate = true;
  }

  setPutty(radii) {
    this._puttyRadii = radii;
    this._applyRadii();
  }

  // Per-residue tube radius: putty (flexibility) and hidden residues (radius 0).
  _applyRadii() {
    if (!this.tube) return;
    const putty = this._puttyRadii;
    if (!putty && !this.hidden.size) this.tube.radii = null;
    else this.tube.radii = Array.from({ length: this.ens.nRes }, (_, r) => (this.hidden.has(r) ? 0 : putty ? putty[r] : this.tube.radius));
    this.update();
  }

  // Console-level display sets (residue indices). Any argument left undefined keeps its current value.
  setDisplaySets({ hidden, sticks, spheres, labels, atomColors } = {}) {
    if (atomColors) { this.atomColorOverrides = new Map(atomColors); if (this.resLinear) { this._applyColors(); } }
    if (hidden) { this.hidden = new Set(hidden); if (this.slabs) this.slabs.hidden = this.hidden; this._atomSubset = null; }
    if (sticks) this.extraSticks = new Set(sticks);
    if (spheres) this.extraSpheres = new Set(spheres);
    if (labels) this.extraLabels = new Map(labels);
    if (hidden) this._applyRadii(); else this.update();
    this.needsRender = true;
  }

  _visibleAtomSubset() {
    if (!this.hidden.size) return null;
    if (this._atomSubset) return this._atomSubset;
    const ens = this.ens, atoms = [];
    for (let r = 0; r < ens.nRes; r++) {
      if (this.hidden.has(r)) continue;
      const res = ens.residues[r];
      for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) atoms.push(a);
    }
    const keep = new Uint8Array(ens.nAtoms);
    for (const a of atoms) keep[a] = 1;
    const bonds = [];
    for (let i = 0; i < ens.bonds.length; i += 2) if (keep[ens.bonds[i]] && keep[ens.bonds[i + 1]]) bonds.push(ens.bonds[i], ens.bonds[i + 1]);
    this._atomSubset = { atoms, bonds };
    return this._atomSubset;
  }

  // residueHex: array of '#rrggbb' per residue
  setResidueColors(residueHex) {
    this.residueHex = residueHex;
    if (!this.ens) return;
    const n = this.ens.nRes;
    const lin = new Float32Array(n * 3);
    const c = new THREE.Color();
    for (let i = 0; i < n; i++) { c.set(residueHex[i]); lin[3 * i] = c.r; lin[3 * i + 1] = c.g; lin[3 * i + 2] = c.b; }
    this.resLinear = lin;
    this._applyColors();
  }

  setResidueColorsLinear(lin) {
    this.resLinear = lin;
    this._applyColors();
  }

  _applyColors() {
    const lin = this.resLinear;
    if (!lin || !this.ens) return;
    const ens = this.ens;
    this.tube.setColors(lin);
    this.slabs.setColors(lin);
    this.glyco.setColors(ens.nRes, (i, out) => { out[0] = lin[3 * i]; out[1] = lin[3 * i + 1]; out[2] = lin[3 * i + 2]; });
    const elemLin = {};
    const c = new THREE.Color();
    for (const [k, hex] of Object.entries(ELEMENT)) { c.set(hex); elemLin[k] = [c.r, c.g, c.b]; }
    const atomHex = this.atomColorOverrides;
    const atomColor = (a, out) => {
      if (atomHex && atomHex.has(a)) { c.set(atomHex.get(a)); out[0] = c.r; out[1] = c.g; out[2] = c.b; return; }
      const e = ens.elements[a];
      if (this.settings.elementColors && e !== 'C' && elemLin[e]) { out[0] = elemLin[e][0]; out[1] = elemLin[e][1]; out[2] = elemLin[e][2]; return; }
      const r = ens.atomRes[a];
      out[0] = lin[3 * r]; out[1] = lin[3 * r + 1]; out[2] = lin[3 * r + 2];
    };
    this._atomColor = atomColor;
    this.atoms.setColors(ens.nAtoms, atomColor);
    const b = ens.bonds;
    this.bondCyl.setColors(b.length, (i, out) => atomColor(b[i], out));    // half-bond i belongs to atom b[i]
    this._updateSelectionReps();
    this._surfaceColors();
    this.needsRender = true;
  }

  setCoords(coords) {
    this.coords = coords;
    this.update();
  }

  update() {
    if (!this.ens || !this.coords) return;
    const ens = this.ens, c = this.coords, R = REP[this.settings.representation];
    if (R.tube) this.tube.update(c);
    if (R.slabs) this.slabs.update(c);
    if (R.glyco) {
      const rs = [];
      for (let r = 0; r < ens.nRes; r++) if (ens.trace[r] >= 0 && ens.glyco[r] >= 0 && ens.rings[r] && !this.hidden.has(r)) rs.push(r);
      this._glycoRes = rs;
      this.glyco.set(rs.length, (i, seg) => {
        const a = ens.trace[rs[i]], b = ens.glyco[rs[i]];
        seg[0] = c[3 * a]; seg[1] = c[3 * a + 1]; seg[2] = c[3 * a + 2]; seg[3] = c[3 * b]; seg[4] = c[3 * b + 1]; seg[5] = c[3 * b + 2];
      }, 0.36);
      if (this.resLinear) {
        const lin = this.resLinear;
        this.glyco.setColors(rs.length, (i, out) => { const r = rs[i]; out[0] = lin[3 * r]; out[1] = lin[3 * r + 1]; out[2] = lin[3 * r + 2]; });
      }
    }
    if (R.atoms) this._updateAtoms(this.atoms, this.bondCyl, this._visibleAtomSubset(), R.atoms);
    this._updateSelectionReps();
    this.needsRender = true;
  }

  _atomRadius(kind, a) {
    const e = this.ens.elements[a];
    if (kind === 'spacefill') return VDW[e] || 1.7;
    if (kind === 'ballstick') return 0.24 * (VDW[e] || 1.7) / 1.7 + 0.08;
    if (kind === 'thin') return 0.17;
    return 0.22;
  }

  _updateAtoms(spheres, cyl, subset, kind) {
    const ens = this.ens, c = this.coords;
    const atomList = subset ? subset.atoms : null;
    const count = atomList ? atomList.length : ens.nAtoms;
    spheres.set(count, (i, p) => { const a = atomList ? atomList[i] : i; p[0] = c[3 * a]; p[1] = c[3 * a + 1]; p[2] = c[3 * a + 2]; },
      i => this._atomRadius(kind, atomList ? atomList[i] : i));
    if (atomList && this._atomColor) spheres.setColors(count, (i, out) => this._atomColor(atomList[i], out));
    if (kind === 'spacefill') { cyl.mesh.count = 0; return; }
    const bonds = subset ? subset.bonds : ens.bonds;
    const nb = bonds.length;       // = 2 * number of bonds; half-bond i spans atom bonds[i] -> midpoint
    const r = kind === 'ballstick' ? 0.11 : this._atomRadius(kind, 0);
    cyl.set(nb, (i, seg) => {
      const a = bonds[i], b = bonds[i ^ 1];
      seg[0] = c[3 * a]; seg[1] = c[3 * a + 1]; seg[2] = c[3 * a + 2];
      seg[3] = (c[3 * a] + c[3 * b]) / 2; seg[4] = (c[3 * a + 1] + c[3 * b + 1]) / 2; seg[5] = (c[3 * a + 2] + c[3 * b + 2]) / 2;
    }, r);
    if (subset && this._atomColor) cyl.setColors(nb, (i, out) => this._atomColor(bonds[i], out));
  }

  setSelection(residues) {
    this.selection = new Set(residues);
    this._updateSelectionReps();
  }
  setHighlight(residues) {
    this.highlight = new Set(residues || []);
    this._updateSelectionReps();
  }

  _updateSelectionReps() {
    if (!this.ens || !this.selAtoms) return;
    const R = REP[this.settings.representation];
    // van der Waals spheres for "show spheres" residues (not needed when everything is spacefill already)
    if (this.vdwAtoms) {
      const sph = R.atoms === 'spacefill' ? [] : [...this.extraSpheres];
      const atoms = [];
      for (const r of sph) { const res = this.ens.residues[r]; for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) atoms.push(a); }
      this.vdwAtoms.mesh.visible = atoms.length > 0;
      if (atoms.length) {
        const c = this.coords;
        this.vdwAtoms.set(atoms.length, (i, p) => { const a = atoms[i]; p[0] = c[3 * a]; p[1] = c[3 * a + 1]; p[2] = c[3 * a + 2]; }, i => this._atomRadius('spacefill', atoms[i]) * 0.92);
        if (this._atomColor) this.vdwAtoms.setColors(atoms.length, (i, out) => this._atomColor(atoms[i], out));
      } else this.vdwAtoms.mesh.count = 0;
    }
    const show = new Set([...this.selection, ...this.highlight, ...this.extraSticks]);
    if (R.atoms || show.size === 0) {
      this.selAtoms.mesh.count = 0; this.selBonds.mesh.count = 0;
      this.selAtoms.mesh.visible = this.selBonds.mesh.visible = false;
      this.needsRender = true;
      return;
    }
    const ens = this.ens;
    const atoms = [];
    for (const r of show) { const res = ens.residues[r]; for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) atoms.push(a); }
    const inSet = new Set(atoms);
    const bonds = [];
    for (let i = 0; i < ens.bonds.length; i += 2) if (inSet.has(ens.bonds[i]) && inSet.has(ens.bonds[i + 1])) bonds.push(ens.bonds[i], ens.bonds[i + 1]);
    this.selAtoms.mesh.visible = this.selBonds.mesh.visible = true;
    this._updateAtoms(this.selAtoms, this.selBonds, { atoms, bonds }, 'sticks');
    this.needsRender = true;
  }

  // ------------------------------------------------------------------ ensemble overlay & preview
  setGhosts(list) {
    this._clearGhosts();
    if (!list || !this.ens) { this.needsRender = true; return; }
    this.ghostMat = makeMaterial(null, { ghost: true, color: 0x7f8794, opacity: this.settings.alpha?.overlay ?? 0.2 });
    for (const coords of list) {
      const t = new Tube(this.ens, { radius: 0.42, radial: 8, sub: 4 });
      t.update(coords);
      const m = new THREE.Mesh(t.geometry, this.ghostMat);
      m.renderOrder = 2;
      m.frustumCulled = false;
      this.ghostGroup.add(m);
    }
    this.needsRender = true;
  }
  _clearGhosts() {
    if (!this.ghostGroup) return;
    for (const m of [...this.ghostGroup.children]) { m.geometry.dispose(); this.ghostGroup.remove(m); }
    this.ghostMat?.dispose();
  }
  setPreview(coords) {
    this._clearPreview();
    if (!coords || !this.ens) { this.needsRender = true; return; }
    const t = new Tube(this.ens, { radius: 0.55, radial: 10, sub: 5 });
    t.update(coords);
    this.previewMat = makeMaterial(null, { ghost: true, color: 0xeb6834, opacity: 0.45 });
    this.previewTube = new THREE.Mesh(t.geometry, this.previewMat);
    this.previewTube.renderOrder = 3;
    this.previewTube.frustumCulled = false;
    this.inner.add(this.previewTube);
    this.needsRender = true;
  }
  _clearPreview() {
    if (!this.previewTube) return;
    this.inner.remove(this.previewTube);
    this.previewTube.geometry.dispose();
    this.previewMat?.dispose();
    this.previewTube = null;
  }

  // ------------------------------------------------------------------ camera
  _cameraFor(aspect, cam = null) {
    const ens = this.ens;
    const R = ens ? ens.radius : 30;
    const F = ens ? (this.fitAll ? ens.radius : ens.fitRadius) : 30;
    const zoom = this.controls.zoom;
    if (this.projection === 'orthographic') {
      const c = cam || this.ortho;
      const hh = (F * 1.04) / zoom;
      c.left = -hh * aspect; c.right = hh * aspect; c.top = hh; c.bottom = -hh;
      c.near = 0.1; c.far = R * 10;
      c.position.set(0, 0, R * 4);
      c.lookAt(0, 0, 0);
      c.updateProjectionMatrix();
      c.updateMatrixWorld(true);
      c.matrixWorldInverse.copy(c.matrixWorld).invert();
      this._depthRef = R * 4;
      this.controls.unitsPerPixel = (2 * hh) / Math.max(1, this.cssSize[1]);
      return c;
    }
    const c = cam || this.persp;
    const fov = c.fov * Math.PI / 180;
    const fit = Math.min(1, aspect);
    const dist = (F / Math.sin(fov / 2)) / fit * 1.04 / zoom;
    c.aspect = aspect;
    c.near = Math.max(0.5, dist - R * 1.3);
    c.far = dist + R * 1.3;
    c.position.set(0, 0, dist);
    c.lookAt(0, 0, 0);
    c.updateProjectionMatrix();
    c.updateMatrixWorld(true);
    c.matrixWorldInverse.copy(c.matrixWorld).invert();
    this._depthRef = dist;
    this.controls.unitsPerPixel = (2 * dist * Math.tan(fov / 2)) / Math.max(1, this.cssSize[1]);
    return c;
  }

  setProjection(p) { this.projection = p === 'orthographic' ? 'orthographic' : 'perspective'; this.needsRender = true; }

  resetView(render = true) {
    this.controls.reset();
    // a pleasant default orientation: principal axis horizontal
    if (this.ens && this.coords) {
      const q = principalOrientation(this.ens, this.coords);
      if (q) this.controls.quaternion.copy(q);
    }
    if (render) this.needsRender = true;
  }

  centerOn(point, animate = true) {
    // move pan so that the model point is at the screen centre
    const p = new THREE.Vector3(...point).sub(new THREE.Vector3(...this.ens.center)).applyQuaternion(this.controls.quaternion);
    const target = new THREE.Vector2(-p.x, -p.y);
    if (!animate) { this.controls.pan.copy(target); this.needsRender = true; return; }
    const from = this.controls.pan.clone();
    const t0 = performance.now();
    const fn = () => {
      const t = Math.min(1, (performance.now() - t0) / 450);
      const e = t < 0.5 ? 2 * t * t : 1 - Math.pow(-2 * t + 2, 2) / 2;
      this.controls.pan.lerpVectors(from, target, e);
      return t < 1;
    };
    this._addAnimator(fn);
  }

  _addAnimator(fn) {
    const wrap = dt => { const go = fn(dt); if (!go) this.animators.delete(wrap); return go; };
    this.animators.add(wrap);
    this.needsRender = true;
  }

  // ------------------------------------------------------------------ render loop
  _loop(t) {
    requestAnimationFrame(tt => this._loop(tt));
    const dt = Math.min(0.1, (t - this._last) / 1000);
    this._last = t;
    let active = this.controls.update(dt);
    for (const fn of [...this.animators]) { try { active = fn(dt) || active; } catch (e) { console.error(e); this.animators.delete(fn); } }
    if (active || this.needsRender) {
      this.needsRender = false;
      this.render();
    }
  }

  _syncPivot() {
    this.pivot.quaternion.copy(this.controls.quaternion);
    this.pivot.position.set(this.controls.pan.x, this.controls.pan.y, 0);
    this.scene.updateMatrixWorld(true);
  }

  _postSettings(extra = {}) {
    const st = this.styles[this.settings.style] || this.styles.studio;
    const R = this.ens ? this.ens.radius : 30;
    const ref = this._depthRef || R * 3;
    const q = this.settings.quality;
    return {
      background: this.settings.background, transparent: false,
      ao: q === 'low' ? 0 : st.ao, aoRadius: st.aoRadius, aoIntensity: st.aoIntensity,
      outline: st.outline, outlineWidth: st.outlineWidth, outlineColor: this.settings.outlineColor,
      fog: st.fog, fogNear: ref - R * 0.4, fogFar: ref + R * 1.1,
      fxaa: true, pixelScale: this.renderer.getPixelRatio(), ...extra,
    };
  }

  render() {
    const [w, h] = this.cssSize;
    const cam = this._cameraFor(w / h);
    this._syncPivot();
    if (!this.ens) {
      this.renderer.setRenderTarget(null);
      this.renderer.setClearColor(new THREE.Color(this.settings.background), 1);
      this.renderer.clear();
    } else {
      this.post.render(this.scene, cam, null, this._postSettings());
    }
    this._drawLabels(cam, w, h);
    this._drawOverlays(cam, w, h);
    if (this.renderHooks.size && this.ens) {
      const proj = this.projector(cam, w, h);
      for (const fn of this.renderHooks) { try { fn(proj, w, h); } catch (e) { console.error(e); } }
    }
    const fs = this.frameStats;
    fs.frames++;
    const now = performance.now();
    if (now - fs.t0 > 1000) { fs.fps = Math.round(fs.frames * 1000 / (now - fs.t0)); fs.frames = 0; fs.t0 = now; }
  }

  // ------------------------------------------------------------------ projection helpers
  projector(cam, w, h) {
    const m = new THREE.Matrix4().multiplyMatrices(cam.projectionMatrix, cam.matrixWorldInverse).multiply(this.inner.matrixWorld);
    const e = m.elements;
    return (x, y, z) => {
      const cx = e[0] * x + e[4] * y + e[8] * z + e[12];
      const cy = e[1] * x + e[5] * y + e[9] * z + e[13];
      const cz = e[2] * x + e[6] * y + e[10] * z + e[14];
      const cw = e[3] * x + e[7] * y + e[11] * z + e[15];
      return [(cx / cw * 0.5 + 0.5) * w, (1 - (cy / cw * 0.5 + 0.5)) * h, cz / cw];
    };
  }

  labelAnchor(r) {
    const ens = this.ens, c = this.coords;
    const a = ens.glyco[r] >= 0 ? ens.glyco[r] : ens.trace[r];
    if (a < 0) return ens.residueCenter(c, r);
    if (ens.rings[r]) {
      let x = 0, y = 0, z = 0;
      const o = ens.rings[r].outline;
      for (const i of o) { x += c[3 * i]; y += c[3 * i + 1]; z += c[3 * i + 2]; }
      return [x / o.length, y / o.length, z / o.length];
    }
    return [c[3 * a], c[3 * a + 1], c[3 * a + 2]];
  }

  labelResidues() {
    if (!this.ens) return [];
    const extra = [...this.extraLabels.keys()];
    let base;
    if (this.labelMode === 'none') base = [];
    else if (this.labelMode === 'all') base = [...Array(this.ens.nRes).keys()];
    else if (this.labelMode === 'every10') base = [...Array(this.ens.nRes).keys()].filter(r => r === 0 || (r + 1) % 10 === 0 || r === this.ens.nRes - 1);
    else base = [...this.selection];
    return [...new Set([...base, ...extra])].filter(r => !this.hidden.has(r) || this.extraLabels.has(r));
  }
  labelText(r) { return this.extraLabels.get(r) || this.ens.resLabel(r); }

  _drawLabels(cam, w, h) {
    const host = this.labelsEl;
    if (!host) return;
    const list = this.ens ? this.labelResidues() : [];
    while (host.children.length > list.length) host.lastChild.remove();
    while (host.children.length < list.length) host.append(document.createElement('div'));
    if (!list.length) return;
    const proj = this.projector(cam, w, h);
    list.forEach((r, k) => {
      const node = host.children[k];
      node.className = 'res-label';
      const p = this.labelAnchor(r);
      const [x, y, z] = proj(...p);
      node.textContent = this.labelText(r);
      node.style.left = `${x}px`;
      node.style.top = `${y}px`;
      node.style.display = z > 1 || z < -1 ? 'none' : '';
    });
  }

  // Overlays: functions returning primitives in model coordinates.
  //   {type:'line', a:[x,y,z], b:[x,y,z], color, dash:true, label:'12.3 Å', width}
  //   {type:'point', p:[x,y,z], color, r}
  addOverlay(id, fn) { this.overlays.set(id, fn); this.needsRender = true; }
  removeOverlay(id) { this.overlays.delete(id); this.needsRender = true; }
  overlayPrimitives() {
    const out = [];
    if (!this.ens || !this.coords) return out;
    for (const fn of this.overlays.values()) {
      try { out.push(...(fn({ coords: this.coords, ens: this.ens }) || [])); } catch (e) { console.error(e); }
    }
    return out;
  }
  _drawOverlays(cam, w, h) {
    const root = this.overlayEl;
    if (!root) return;
    root.replaceChildren();
    const prims = this.overlayPrimitives();
    if (!prims.length) return;
    const proj = this.projector(cam, w, h);
    for (const p of prims) {
      if (p.type === 'line') {
        const [x1, y1] = proj(...p.a), [x2, y2] = proj(...p.b);
        root.append(svg('line', { x1, y1, x2, y2, stroke: '#ffffff', 'stroke-width': (p.width || 2) + 3, 'stroke-linecap': 'round', opacity: 0.8 }));
        root.append(svg('line', { x1, y1, x2, y2, stroke: p.color || '#111', 'stroke-width': p.width || 2, 'stroke-dasharray': p.dash ? '6 4' : null, 'stroke-linecap': 'round' }));
        if (p.label) {
          const tx = (x1 + x2) / 2, ty = (y1 + y2) / 2 - 8;
          root.append(svg('text', { x: tx, y: ty, 'text-anchor': 'middle', 'font-size': 12, 'font-weight': 650, fill: '#111', stroke: '#fff', 'stroke-width': 3, 'paint-order': 'stroke', text: p.label }));
        }
      } else if (p.type === 'point') {
        const [x, y] = proj(...p.p);
        root.append(svg('circle', { cx: x, cy: y, r: p.r || 4, fill: p.color || '#111', stroke: '#fff', 'stroke-width': 1.5 }));
      }
    }
  }

  // Model coordinates under the mouse, on the plane through `refPoint` (model coordinates) facing the viewer.
  screenToModel(clientX, clientY, refPoint) {
    const [w, h] = this.cssSize;
    const cam = this._cameraFor(w / h);
    this._syncPivot();
    const rect = this.canvas.getBoundingClientRect();
    const nx = ((clientX - rect.left) / rect.width) * 2 - 1, ny = -(((clientY - rect.top) / rect.height) * 2 - 1);
    const ref = new THREE.Vector3(...refPoint).applyMatrix4(this.inner.matrixWorld).applyMatrix4(cam.matrixWorldInverse);
    const a = new THREE.Vector3(nx, ny, -1).applyMatrix4(cam.projectionMatrixInverse);
    const b = new THREE.Vector3(nx, ny, 1).applyMatrix4(cam.projectionMatrixInverse);
    const t = Math.abs(b.z - a.z) < 1e-9 ? 0 : (ref.z - a.z) / (b.z - a.z);
    const q = a.lerp(b, t).applyMatrix4(cam.matrixWorld);
    const inv = new THREE.Matrix4().copy(this.inner.matrixWorld).invert();
    q.applyMatrix4(inv);
    return [q.x, q.y, q.z];
  }

  // Viewing direction (unit vector toward the camera) in model coordinates.
  viewDirModel() {
    const v = new THREE.Vector3(0, 0, 1).applyQuaternion(this.controls.quaternion.clone().invert());
    return [v.x, v.y, v.z];
  }

  gpuName() {
    try {
      const gl = this.renderer.getContext();
      const ext = gl.getExtension('WEBGL_debug_renderer_info');
      return cleanGpuName(ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER));
    } catch { return ''; }
  }

  // ------------------------------------------------------------------ picking
  _bindPicking() {
    let down = null;
    this.canvas.addEventListener('pointerdown', e => { down = { x: e.clientX, y: e.clientY, t: performance.now() }; });
    this.canvas.addEventListener('pointerup', e => {
      if (!down) return;
      const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y);
      if (moved < 4 && e.button === 0) this.onPick(this.pick(e.clientX, e.clientY), e);
      down = null;
    });
    this.canvas.addEventListener('dblclick', e => this.onDblPick(this.pick(e.clientX, e.clientY), e));
    let pending = null;
    this.canvas.addEventListener('pointermove', e => {
      if (e.buttons) { this.onHover(null, e); return; }
      pending = e;
      if (this._hoverRaf) return;
      this._hoverRaf = requestAnimationFrame(() => { this._hoverRaf = 0; this.onHover(this.pick(pending.clientX, pending.clientY), pending); });
    });
    this.canvas.addEventListener('pointerleave', e => this.onHover(null, e));
  }

  pick(clientX, clientY) {
    if (!this.ens || !this.coords) return null;
    const rect = this.canvas.getBoundingClientRect();
    const x = clientX - rect.left, y = clientY - rect.top;
    const [w, h] = this.cssSize;
    const cam = this._cameraFor(w / h);
    this._syncPivot();
    const proj = this.projector(cam, w, h);
    const ens = this.ens, c = this.coords;
    const R = REP[this.settings.representation];
    const unit = this.controls.unitsPerPixel || 0.05;
    const radiusPx = Math.max(6, (R.atoms === 'spacefill' ? 1.6 : R.atoms ? 0.5 : 1.3) / unit);
    let best = null;
    for (let a = 0; a < ens.nAtoms; a++) {
      const [px, py, pz] = proj(c[3 * a], c[3 * a + 1], c[3 * a + 2]);
      const d = Math.hypot(px - x, py - y);
      if (d > radiusPx) continue;
      const score = pz + d * 1e-4;
      if (!best || score < best.score) best = { atom: a, residue: ens.atomRes[a], score, px, py };
    }
    return best;
  }

  // ------------------------------------------------------------------ export
  // Render the current view at a given pixel size into a 2D canvas (plus a projector for annotations).
  renderImage({ width, height, transparent = false, background = null, superSample = 0 } = {}) {
    if (!this.ens) throw new Error('nothing to render');
    const maxTex = this.renderer.capabilities.maxTextureSize;
    let ss = superSample || (width * height <= 9e6 ? 2 : 1);
    while (ss > 1 && (width * ss > maxTex || height * ss > maxTex || width * height * ss * ss > 6.4e7)) ss--;
    const W = width * ss, H = height * ss;
    if (W > maxTex || H > maxTex) throw new Error(`image too large for this GPU (max ${maxTex} px per side)`);
    const cam = this._cameraFor(width / height, this.projection === 'orthographic' ? this.ortho.clone() : this.persp.clone());
    // restore interactive camera parameters for the live view afterwards
    this._syncPivot();
    const post = new PostFX(this.renderer);
    post.setSize(W, H);
    const out = new THREE.WebGLRenderTarget(W, H, { type: THREE.UnsignedByteType, format: THREE.RGBAFormat, depthBuffer: false });
    const viewH = this.cssSize[1] * this.renderer.getPixelRatio();
    const settings = this._postSettings({ transparent, pixelScale: Math.max(1, H / Math.max(1, viewH)) * this.renderer.getPixelRatio(), fxaa: ss === 1 });
    if (background) settings.background = background;
    post.render(this.scene, cam, out, settings);
    const buf = new Uint8Array(W * H * 4);
    this.renderer.readRenderTargetPixels(out, 0, 0, W, H, buf);
    post.dispose();
    out.dispose();
    const big = document.createElement('canvas');
    big.width = W; big.height = H;
    const ctx = big.getContext('2d');
    const img = ctx.createImageData(W, H);
    const row = W * 4;
    for (let y = 0; y < H; y++) img.data.set(buf.subarray((H - 1 - y) * row, (H - y) * row), y * row);
    ctx.putImageData(img, 0, 0);
    let canvas = big;
    if (ss > 1) {
      canvas = document.createElement('canvas');
      canvas.width = width; canvas.height = height;
      const c2 = canvas.getContext('2d');
      c2.imageSmoothingEnabled = true;
      c2.imageSmoothingQuality = 'high';
      c2.drawImage(big, 0, 0, width, height);
    }
    this.needsRender = true;
    this._cameraFor(this.cssSize[0] / this.cssSize[1]);
    return { canvas, project: this.projector(cam, width, height), scale: height / this.cssSize[1] };
  }

  // Persistent offscreen pipeline for video capture.
  beginCapture(width, height) {
    const post = new PostFX(this.renderer);
    post.setSize(width, height);
    const out = new THREE.WebGLRenderTarget(width, height, { type: THREE.UnsignedByteType, format: THREE.RGBAFormat, depthBuffer: false });
    const buf = new Uint8Array(width * height * 4);
    const canvas = document.createElement('canvas');
    canvas.width = width; canvas.height = height;
    const ctx = canvas.getContext('2d');
    const img = ctx.createImageData(width, height);
    const viewH = this.cssSize[1] * this.renderer.getPixelRatio();
    const self = this;
    return {
      canvas,
      frame() {
        const cam = self._cameraFor(width / height, self.projection === 'orthographic' ? self.ortho.clone() : self.persp.clone());
        self._syncPivot();
        post.render(self.scene, cam, out, self._postSettings({ pixelScale: Math.max(1, height / Math.max(1, viewH)) * self.renderer.getPixelRatio(), fxaa: true }));
        self.renderer.readRenderTargetPixels(out, 0, 0, width, height, buf);
        const row = width * 4;
        for (let y = 0; y < height; y++) img.data.set(buf.subarray((height - 1 - y) * row, (height - y) * row), y * row);
        ctx.putImageData(img, 0, 0);
        self._cameraFor(self.cssSize[0] / self.cssSize[1]);
        return { canvas, project: self.projector(cam, width, height), scale: height / self.cssSize[1] };
      },
      end() { post.dispose(); out.dispose(); self.needsRender = true; },
    };
  }
}

// "ANGLE (Intel, Intel(R) Iris(R) Xe Graphics (0x00009A49) Direct3D11 ...)" -> "Intel Iris Xe Graphics"
export function cleanGpuName(raw) {
  let s = String(raw || '');
  const m = /^ANGLE \((.*)\)$/.exec(s);
  if (m) { const parts = m[1].split(', '); s = parts.length >= 2 ? parts[1] : parts[0]; }
  if (s.includes('Renderer: ')) s = s.split('Renderer: ')[1];
  s = s.replace(/\s*\(0x[0-9a-f]+\).*$/i, '').replace(/\s+(Direct3D|OpenGL|Vulkan)\S*.*$/i, '').replace(/\((R|TM)\)/gi, '');
  return s.replace(/\s+/g, ' ').trim();
}

// Orientation that lays the ensemble's principal axis horizontally (largest extent along screen X).
function principalOrientation(ens, coords) {
  const n = ens.nRes;
  const pts = [];
  for (let r = 0; r < n; r++) { const a = ens.trace[r]; if (a >= 0) pts.push([coords[3 * a], coords[3 * a + 1], coords[3 * a + 2]]); }
  if (pts.length < 3) return null;
  const m = [0, 0, 0];
  for (const p of pts) { m[0] += p[0]; m[1] += p[1]; m[2] += p[2]; }
  m[0] /= pts.length; m[1] /= pts.length; m[2] /= pts.length;
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) C[i][j] += (p[i] - m[i]) * (p[j] - m[j]);
  // power iteration for the two leading axes
  const lead = (M, orth) => {
    let v = [1, 0.5, 0.25];
    for (let it = 0; it < 60; it++) {
      let w = [0, 1, 2].map(i => M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2]);
      if (orth) { const d = w[0] * orth[0] + w[1] * orth[1] + w[2] * orth[2]; w = w.map((x, i) => x - d * orth[i]); }
      const l = Math.hypot(...w) || 1;
      v = w.map(x => x / l);
    }
    return v;
  };
  const e1 = lead(C, null);
  const e2 = lead(C, e1);
  const e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
  const mat = new THREE.Matrix4().set(e1[0], e1[1], e1[2], 0, e2[0], e2[1], e2[2], 0, e3[0], e3[1], e3[2], 0, 0, 0, 0, 1);
  return new THREE.Quaternion().setFromRotationMatrix(mat);
}
