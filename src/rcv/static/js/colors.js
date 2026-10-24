// Color schemes for the molecule. Palettes follow the validated categorical order of the app's
// data-viz palette; the nucleotide set (A green, C blue, G amber, U pink) passes colour-vision-deficiency
// separation for all pairs, and every scheme ships with a legend so identity never relies on colour alone.
import { hexToRgb, rgbToHex, clamp } from './util.js';

export const NUCLEOTIDE = { A: '#008300', C: '#2a78d6', G: '#eda100', U: '#e87ba4', N: '#8f8d86' };
export const NUCLEOTIDE_CLASSIC = { A: '#008300', C: '#2a78d6', G: '#eda100', U: '#e34948', N: '#8f8d86' };
export const SLOTS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'];
export const LOOP = '#c4c2ba';
export const LONE_PAIR = '#9b998f';
export const PAIRED = '#2a78d6';
export const PSEUDOKNOT = '#eb6834';
export const BLUE_RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6',
  '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'];
export const ORANGE_RAMP = ['#fde3d6', '#f9c4aa', '#f5a47f', '#f08456', '#eb6834', '#cf5422', '#a9421a', '#833313'];
export const VIRIDIS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'];
export const ELEMENT = { N: '#2f4fd1', O: '#e03131', P: '#f08c00', S: '#e8b800', F: '#5fb34f', CL: '#3fae49', BR: '#a1522c', SE: '#d08a2a', I: '#7b3fa0' };

export function mix(a, b, t) {
  const x = hexToRgb(a), y = hexToRgb(b);
  return rgbToHex([0, 1, 2].map(i => x[i] + (y[i] - x[i]) * t));
}

export function sampleRamp(ramp, t) {
  t = clamp(Number.isFinite(t) ? t : 0, 0, 1) * (ramp.length - 1);
  const i = Math.min(ramp.length - 2, Math.floor(t));
  return mix(ramp[i], ramp[i + 1], t - i);
}

export function helixColor(k) {
  if (k === -2) return LONE_PAIR;
  if (k < 0) return LOOP;
  const base = SLOTS[k % SLOTS.length];
  const tier = Math.floor(k / SLOTS.length) % 3;
  return tier === 0 ? base : tier === 1 ? mix(base, '#ffffff', 0.45) : mix(base, '#000000', 0.3);
}

// Structure elements of the shown conformation (categorical, fixed order).
export const ELEMENTS = [
  ['stem', 'Stem (paired)', '#2a78d6'], ['hairpin', 'Hairpin loop', '#eb6834'], ['internal', 'Internal loop', '#1baf7a'],
  ['bulge', 'Bulge', '#eda100'], ['multi', 'Multibranch loop', '#6250d6'], ['pk', 'Pseudoknot pair', '#e34948'], ['exterior', 'Exterior / ends', '#c4c2ba'],
];
const ELEMENT_COLOR = Object.fromEntries(ELEMENTS.map(([k, , c]) => [k, c]));
const ELEMENT_SLOT = { stem: 0, hairpin: 1, internal: 2, bulge: 3, multi: 6, pk: 7 };

// Active colour palette (Color ▸ Palettes): categorical slots and nucleotide colours used by every scheme.
export function applyPalette(pal) {
  if (!pal || !Array.isArray(pal.colors) || pal.colors.length < 4) return;
  SLOTS.splice(0, SLOTS.length, ...pal.colors);
  const nt = pal.nucleotides || { A: pal.colors[2] || pal.colors[0], C: pal.colors[0], G: pal.colors[3] || pal.colors[1], U: pal.colors[1] };
  for (const b of ['A', 'C', 'G', 'U']) if (nt[b]) NUCLEOTIDE[b] = nt[b];
  for (const row of ELEMENTS) {
    const slot = ELEMENT_SLOT[row[0]];
    if (slot !== undefined) row[2] = SLOTS[Math.min(slot, SLOTS.length - 1)] ?? row[2];
    ELEMENT_COLOR[row[0]] = row[2];
  }
}

// Loop type of every residue from the nested pairs of a conformation (pseudoknotted pairs marked separately).
export function structureElements(n, pairs) {
  const pt = new Int32Array(n).fill(-1);
  const out = new Array(n).fill('exterior');
  for (const [i, j, , pk] of pairs) {
    if (pk) { out[i] = out[j] = 'pk'; continue; }
    if (pt[i] < 0 && pt[j] < 0) { pt[i] = j; pt[j] = i; }
  }
  for (let i = 0; i < n; i++) if (pt[i] >= 0) out[i] = 'stem';
  // for each pair (p, q): classify the loop it closes
  for (let p = 0; p < n; p++) {
    const q = pt[p];
    if (q <= p) continue;
    const unpaired = [];
    let branches = 0, k = p + 1;
    const left = [], right = [];
    let seenBranch = false;
    while (k < q) {
      if (pt[k] > k) { branches++; seenBranch = true; k = pt[k] + 1; continue; }
      unpaired.push(k);
      (seenBranch ? right : left).push(k);
      k++;
    }
    const kind = branches === 0 ? 'hairpin' : branches === 1 ? (left.length && right.length ? 'internal' : 'bulge') : 'multi';
    for (const u of unpaired) if (out[u] !== 'pk') out[u] = kind;
  }
  return out;
}

export const SCHEMES = [
  { id: 'nucleotide', label: 'Nucleotide' },
  { id: 'helix', label: 'Helices' },
  { id: 'elements', label: 'Structure elements' },
  { id: 'pairing', label: 'Base pairing' },
  { id: 'pairfreq', label: 'Pairing frequency (ensemble)' },
  { id: 'switching', label: 'Switching residues' },
  { id: 'flexibility', label: 'Flexibility (RMSF)' },
  { id: 'displacement', label: 'Displacement from previous' },
  { id: 'position', label: 'Position 5′→3′' },
  { id: 'nucleotide-classic', label: 'Nucleotide (classic, U red)' },
  { id: 'uniform', label: 'Uniform' },
];

const custom = new Map();     // plugin-registered schemes
export function registerScheme(def) {
  if (!def || !def.id || typeof def.residueColors !== 'function') throw new Error('scheme needs id and residueColors()');
  custom.set(def.id, def);
  if (!SCHEMES.find(s => s.id === def.id)) SCHEMES.push({ id: def.id, label: def.label || def.id, plugin: true });
}

// ctx: { conf, prevCoords, coords, uniform }
export function residueColors(ens, scheme, ctx) {
  const n = ens.nRes;
  const out = new Array(n);
  const k = ctx.conf;
  if (custom.has(scheme)) {
    const cols = custom.get(scheme).residueColors(ens, ctx) || [];
    for (let i = 0; i < n; i++) out[i] = cols[i] || LOOP;
    return out;
  }
  switch (scheme) {
    case 'nucleotide':
    case 'nucleotide-classic': {
      const pal = scheme === 'nucleotide' ? NUCLEOTIDE : NUCLEOTIDE_CLASSIC;
      for (let i = 0; i < n; i++) out[i] = pal[ens.residues[i].name] || pal.N;
      return out;
    }
    case 'helix': {
      const h = ens.conformations[k].helix_of_residue || [];
      for (let i = 0; i < n; i++) out[i] = helixColor(h[i] ?? -1);
      return out;
    }
    case 'pairing': {
      out.fill(LOOP);
      for (const [i, j, , pk] of ens.conformations[k].pairs) { out[i] = out[j] = pk ? PSEUDOKNOT : PAIRED; }
      return out;
    }
    case 'elements': {
      const el = structureElements(n, ens.conformations[k].pairs);
      for (let i = 0; i < n; i++) out[i] = ELEMENT_COLOR[el[i]];
      return out;
    }
    case 'pairfreq': {
      const f = pairedFraction(ens);
      for (let i = 0; i < n; i++) out[i] = sampleRamp(BLUE_RAMP, f[i]);
      return out;
    }
    case 'switching': {
      const s = switchingClass(ens);
      for (let i = 0; i < n; i++) out[i] = s[i] === 'switch' ? '#eb6834' : s[i] === 'paired' ? PAIRED : LOOP;
      return out;
    }
    case 'flexibility': {
      const r = ens.rmsf;
      const max = Math.max(1e-6, ...r);
      for (let i = 0; i < n; i++) out[i] = sampleRamp(BLUE_RAMP.slice(1), r[i] / max);
      return out;
    }
    case 'displacement': {
      const d = ens.residueDisplacement(ctx.coords, ctx.prevCoords);
      const max = Math.max(2.0, ...d);
      for (let i = 0; i < n; i++) out[i] = sampleRamp(ORANGE_RAMP, d[i] / max);
      ctx.displacementMax = max;
      return out;
    }
    case 'position':
      for (let i = 0; i < n; i++) out[i] = sampleRamp(VIRIDIS, n > 1 ? i / (n - 1) : 0);
      return out;
    default:
      out.fill(ctx.uniform || '#5b7a99');
      return out;
  }
}

export function legend(ens, scheme, ctx) {
  if (!ens) return null;
  if (custom.has(scheme)) return custom.get(scheme).legend?.(ens, ctx) || { title: SCHEMES.find(s => s.id === scheme)?.label, items: [] };
  switch (scheme) {
    case 'nucleotide':
    case 'nucleotide-classic': {
      const pal = scheme === 'nucleotide' ? NUCLEOTIDE : NUCLEOTIDE_CLASSIC;
      const present = new Set(ens.residues.map(r => r.name));
      return { title: 'Nucleotide', items: ['A', 'C', 'G', 'U'].filter(b => present.has(b)).map(b => ({ color: pal[b], label: b })) };
    }
    case 'helix': {
      const h = ens.conformations[ctx.conf].helix_of_residue || [];
      const used = [...new Set(h.filter(v => v >= 0))].sort((a, b) => a - b);
      const items = used.map(f => {
        const fam = ens.doc.helix_families[f];
        return { color: helixColor(f), label: `${fam.name} ${fam.span[0] + 1}–${fam.span[1] + 1}` };
      });
      items.push({ color: LOOP, label: 'Unpaired' });
      if (h.includes(-2)) items.push({ color: LONE_PAIR, label: 'Isolated pair' });
      return { title: 'Helices (this conformation)', items };
    }
    case 'pairing':
      return { title: 'Base pairing', items: [{ color: PAIRED, label: 'Paired' }, { color: PSEUDOKNOT, label: 'Pseudoknot' }, { color: LOOP, label: 'Unpaired' }] };
    case 'elements': {
      const used = new Set(structureElements(ens.nRes, ens.conformations[ctx.conf].pairs));
      return { title: 'Structure elements (this conformation)', items: ELEMENTS.filter(([k]) => used.has(k)).map(([, label, color]) => ({ color, label })) };
    }
    case 'pairfreq':
      return { title: 'Fraction of conformations in which the residue is paired', ramp: BLUE_RAMP, min: '0%', max: '100%' };
    case 'switching': {
      const s = switchingClass(ens);
      return { title: 'Pairing across the ensemble', items: [{ color: '#eb6834', label: `Switching (${s.filter(x => x === 'switch').length})` }, { color: PAIRED, label: 'Always same partner' }, { color: LOOP, label: 'Always unpaired' }] };
    }
    case 'flexibility': {
      const max = Math.max(...ens.rmsf);
      return { title: 'Per-residue RMSF across the ensemble', ramp: BLUE_RAMP.slice(1), min: '0 Å', max: `${max.toFixed(1)} Å` };
    }
    case 'displacement':
      return { title: 'Displacement from previous conformation', ramp: ORANGE_RAMP, min: '0 Å', max: `${(ctx.displacementMax || 2).toFixed(1)} Å` };
    case 'position':
      return { title: 'Sequence position', ramp: VIRIDIS, min: '5′', max: '3′' };
    default:
      return { title: 'Uniform', items: [{ color: ctx.uniform || '#5b7a99', label: ens.name }] };
  }
}

export const BACKGROUNDS = [
  { id: 'white', label: 'White', color: '#ffffff' },
  { id: 'paper', label: 'Paper', color: '#f7f6f2' },
  { id: 'mist', label: 'Light gray', color: '#e9ebee' },
  { id: 'slate', label: 'Slate', color: '#2b2f36' },
  { id: 'dark', label: 'Dark', color: '#151514' },
  { id: 'black', label: 'Black', color: '#000000' },
];

// Fraction of conformations in which each residue is base-paired (cached per ensemble).
export function pairedFraction(ens) {
  if (ens._pairedFraction) return ens._pairedFraction;
  const f = new Array(ens.nRes).fill(0);
  for (const c of ens.conformations) { const seen = new Set(); for (const [i, j] of c.pairs) { seen.add(i); seen.add(j); } for (const r of seen) f[r]++; }
  return (ens._pairedFraction = f.map(v => v / Math.max(1, ens.n)));
}

// 'switch' when the pairing partner differs between conformations, else 'paired' / 'unpaired'.
export function switchingClass(ens) {
  if (ens._switching) return ens._switching;
  const parts = ens.conformations.map((_, k) => ens.partner(k));
  return (ens._switching = Array.from({ length: ens.nRes }, (_, r) => (parts.some(p => p[r] !== parts[0][r]) ? 'switch' : parts[0][r] >= 0 ? 'paired' : 'unpaired')));
}
