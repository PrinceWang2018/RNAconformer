// Small DOM, formatting and file helpers shared by the app.

export function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class') node.className = v;
    else if (k === 'text') node.textContent = v;
    else if (k === 'html') node.innerHTML = v;          // only used with trusted static markup
    else if (k.startsWith('on') && typeof v === 'function') node.addEventListener(k.slice(2), v);
    else if (k === 'style' && typeof v === 'object') { for (const [sk, sv] of Object.entries(v)) { if (sk.startsWith('--')) node.style.setProperty(sk, sv); else node.style[sk] = sv; } }
    else if (k === 'dataset') Object.assign(node.dataset, v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return node;
}

const SVGNS = 'http://www.w3.org/2000/svg';
export function svg(tag, attrs = {}, ...children) {
  const node = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'text') node.textContent = v;
    else if (k.startsWith('on') && typeof v === 'function') node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return node;
}

// Font of every exported figure, video and plot: Liberation Sans (SIL Open Font License, bundled in static/fonts;
// same metrics as the Helvetica/Arial family that journals expect), loaded as "RCV Figure Sans".
export const FIG_FONT = '"RCV Figure Sans", "Liberation Sans", Arimo, "Helvetica Neue", Helvetica, sans-serif';
export async function ensureFigureFont() {
  try { await Promise.all([document.fonts.load('12px "RCV Figure Sans"'), document.fonts.load('bold 12px "RCV Figure Sans"')]); } catch { /* fallback fonts */ }
}
export const figFont = (px, bold = false) => `${bold ? 'bold' : 'normal'} ${Math.max(1, px)}px ${FIG_FONT}`;

export const $ = (sel, root = document) => root.querySelector(sel);
export const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];
export const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
export const lerp = (a, b, t) => a + (b - a) * t;
export const easeInOut = t => (t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2);

export class Emitter {
  constructor() { this._h = new Map(); }
  on(name, fn) {
    if (!this._h.has(name)) this._h.set(name, new Set());
    this._h.get(name).add(fn);
    return () => this._h.get(name)?.delete(fn);
  }
  emit(name, detail) {
    for (const fn of [...(this._h.get(name) || [])]) {
      try { fn(detail); } catch (err) { console.error(`handler for ${name} failed`, err); }
    }
  }
}

export function fmt(v, digits = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return '–';
  if (typeof v !== 'number') return String(v);
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e5 || a < 1e-3)) return v.toExponential(2);
  if (Number.isInteger(v) && a < 1e5) return v.toLocaleString('en-US');
  const d = a >= 100 ? Math.min(digits, 1) : a >= 10 ? Math.min(digits, 2) : digits;
  return v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
}
// Metric values: counts stay integers, everything else gets fixed decimals (1 for |v| >= 100, else 2).
export function fmtVal(v, unit = 'x') {
  if (v === null || v === undefined || Number.isNaN(v)) return '–';
  if (unit === '' && Number.isInteger(v)) return v.toLocaleString('en-US');
  const d = Math.abs(v) >= 100 ? 1 : 2;
  return v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
}
export function fmtSigned(v, digits = 2) {
  if (v === null || v === undefined || Number.isNaN(v)) return '–';
  return (v > 0 ? '+' : v < 0 ? '−' : '±') + fmt(Math.abs(v), digits);
}

export function niceTicks(min, max, count = 4) {
  if (!Number.isFinite(min) || !Number.isFinite(max)) return [0, 1];
  if (min === max) { const p = Math.abs(min) * 0.1 || 1; min -= p; max += p; }
  const span = max - min;
  const step0 = span / Math.max(1, count);
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const err = step0 / mag;
  const step = (err >= 7.5 ? 10 : err >= 3.5 ? 5 : err >= 1.5 ? 2 : 1) * mag;
  const ticks = [];
  for (let v = Math.ceil(min / step - 1e-9) * step; v <= max + step * 1e-9; v += step) ticks.push(Math.abs(v) < step * 1e-9 ? 0 : v);
  return ticks;
}

export function b64ToFloat32(b64) {
  const bin = atob(b64);
  const buf = new ArrayBuffer(bin.length);
  const u8 = new Uint8Array(buf);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  return new Float32Array(buf);
}
export function float32ToB64(arr) {
  const u8 = new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength);
  let s = '';
  for (let i = 0; i < u8.length; i += 0x8000) s += String.fromCharCode.apply(null, u8.subarray(i, i + 0x8000));
  return btoa(s);
}

export function hexToRgb(hex) {
  const h = hex.replace('#', '');
  const n = parseInt(h.length === 3 ? h.split('').map(c => c + c).join('') : h, 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}
export function rgbToHex([r, g, b]) {
  return '#' + [r, g, b].map(v => Math.round(clamp(v, 0, 255)).toString(16).padStart(2, '0')).join('');
}
export function luminance(hex) {
  const [r, g, b] = hexToRgb(hex).map(v => { v /= 255; return v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

export function safeStorage() {
  try { const k = '__t'; localStorage.setItem(k, '1'); localStorage.removeItem(k); return localStorage; } catch { return null; }
}
export function loadPref(key, fallback) {
  try {
    let v = localStorage.getItem('rcv.' + key);
    if (v === null) v = localStorage.getItem('rnaconfswitch.' + key);      // preferences of version 0.4
    return v === null ? fallback : JSON.parse(v);
  } catch { return fallback; }
}
export function savePref(key, value) {
  try { localStorage.setItem('rcv.' + key, JSON.stringify(value)); } catch { /* storage unavailable */ }
}

// Save a Blob with the native save dialog when available (Chromium), otherwise as a download.
export async function saveBlob(blob, filename, types) {
  if (window.__native && window.__native.available()) return window.__native.save(blob, filename, types?.[0]?.description);
  if (window.showSaveFilePicker && !window.__noFilePicker) {
    try {
      const ext = filename.slice(filename.lastIndexOf('.'));
      const handle = await window.showSaveFilePicker({
        suggestedName: filename,
        types: types || [{ description: ext.slice(1).toUpperCase() + ' file', accept: { [blob.type || 'application/octet-stream']: [ext] } }],
      });
      const w = await handle.createWritable();
      await w.write(blob);
      await w.close();
      return handle.name;
    } catch (err) {
      if (err && err.name === 'AbortError') return null;
      // fall through to download (e.g. unsupported type in picker)
    }
  }
  const url = URL.createObjectURL(blob);
  const a = el('a', { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 30000);
  return filename;
}

export function readFileText(file) {
  return new Promise((resolve, reject) => {
    const r = new FileReader();
    r.onload = () => resolve(r.result);
    r.onerror = () => reject(r.error);
    r.readAsText(file);
  });
}

export function debounce(fn, ms) {
  let t = 0;
  return (...a) => { clearTimeout(t); t = setTimeout(() => fn(...a), ms); };
}

export function slug(s) {
  return String(s || 'rna').replace(/[^A-Za-z0-9_.-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 60) || 'rna';
}

export function sleep(ms) { return new Promise(r => setTimeout(r, ms)); }
export function nextFrame() { return new Promise(r => requestAnimationFrame(() => r())); }
