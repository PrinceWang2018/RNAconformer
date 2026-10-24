// Palette data and the colour-vision check (no DOM; also used by tests).
import { hexToRgb } from './util.js';

export const BUILTIN = [
  { id: 'rcv', name: 'RCV default', colors: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    nucleotides: { A: '#008300', C: '#2a78d6', G: '#eda100', U: '#e87ba4' } },
  { id: 'okabe-ito', name: 'Okabe–Ito (colour-blind safe)', colors: ['#0072B2', '#E69F00', '#009E73', '#F0E442', '#CC79A7', '#56B4E9', '#D55E00', '#000000'],
    nucleotides: { A: '#009E73', C: '#0072B2', G: '#E69F00', U: '#CC79A7' } },
  { id: 'tol-bright', name: 'Tol bright (colour-blind safe)', colors: ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB', '#000000'],
    nucleotides: { A: '#228833', C: '#4477AA', G: '#CCBB44', U: '#EE6677' } },
  { id: 'tol-muted', name: 'Tol muted (colour-blind safe)', colors: ['#332288', '#CC6677', '#117733', '#DDCC77', '#88CCEE', '#44AA99', '#AA4499', '#882255'],
    nucleotides: { A: '#117733', C: '#332288', G: '#DDCC77', U: '#CC6677' } },
  { id: 'pastel', name: 'Pastel', colors: ['#7aa6d6', '#f0a07a', '#86c9a5', '#f2cf6b', '#e8a3c1', '#a8d18d', '#a99ae0', '#e88f8e'],
    nucleotides: { A: '#86c9a5', C: '#7aa6d6', G: '#f2cf6b', U: '#e8a3c1' } },
  { id: 'gray', name: 'Grayscale (print)', colors: ['#2b2b2b', '#8c8c8c', '#555555', '#bdbdbd', '#404040', '#a6a6a6', '#6e6e6e', '#d4d4d4'],
    nucleotides: { A: '#2b2b2b', C: '#6e6e6e', G: '#a3a3a3', U: '#d0d0d0' } },
];

// ------------------------------------------------------------------ colour-vision check (OKLab distance × 100)
const lin = v => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
function oklab([r, g, b]) {
  const l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b, m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b, s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b;
  const l_ = Math.cbrt(l), m_ = Math.cbrt(m), s_ = Math.cbrt(s);
  return [0.2104542553 * l_ + 0.7936177850 * m_ - 0.0040720468 * s_, 1.9779984951 * l_ - 2.4285922050 * m_ + 0.4505937099 * s_, 0.0259040371 * l_ + 0.7827717662 * m_ - 0.8086757660 * s_];
}
// Machado et al. 2009, severity 1.0 (linear RGB)
const CVD = {
  protan: [[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]],
  deutan: [[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]],
  tritan: [[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]],
};
function sim(hex, kind) {
  const c = hexToRgb(hex).map(lin);
  if (!kind) return c;
  const M = CVD[kind];
  return M.map(row => Math.min(1, Math.max(0, row[0] * c[0] + row[1] * c[1] + row[2] * c[2])));
}
const dist = (a, b) => { const p = oklab(a), q = oklab(b); return Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) * 100; };

// Smallest pairwise distance of the first n colours, for normal vision and the three colour-vision deficiencies.
export function checkPalette(colors, n = 8) {
  const cols = colors.slice(0, n);
  const out = {};
  for (const kind of [null, 'protan', 'deutan', 'tritan']) {
    let worst = { d: Infinity };
    for (let i = 0; i < cols.length; i++) for (let j = i + 1; j < cols.length; j++) {
      const d = dist(sim(cols[i], kind), sim(cols[j], kind));
      if (d < worst.d) worst = { d, i, j };
    }
    out[kind || 'normal'] = worst;
  }
  return out;
}

