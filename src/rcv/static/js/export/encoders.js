// Dependency-free encoders for publication output: TIFF (Adobe Deflate + predictor, DPI tags),
// PNG with a pHYs (DPI) chunk, PDF (raster image + vector graphics/text) and ZIP (stored) archives.
// Works in browsers and in Node >= 18 (CompressionStream, Blob, Response are globals there too).

const CRC_TABLE = (() => {
  const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
  return t;
})();
export function crc32(bytes, crc = 0) {
  let c = (crc ^ 0xffffffff) >>> 0;
  for (let i = 0; i < bytes.length; i++) c = CRC_TABLE[(c ^ bytes[i]) & 255] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

export async function zlibDeflate(u8) {
  const stream = new Blob([u8]).stream().pipeThrough(new CompressionStream('deflate'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

// ------------------------------------------------------------------ TIFF
export async function encodeTIFF(img, { dpi = 300, alpha = false, software = 'RCV (RNA Conformation Viewer)' } = {}) {
  const { width: w, height: h, data } = img;
  const spp = alpha ? 4 : 3;
  const raw = new Uint8Array(w * h * spp);
  for (let y = 0; y < h; y++) {
    const p = [0, 0, 0, 0];
    for (let x = 0; x < w; x++) {
      const si = (y * w + x) * 4, di = (y * w + x) * spp;
      for (let c = 0; c < spp; c++) { const v = data[si + c]; raw[di + c] = (v - p[c]) & 255; p[c] = v; }
    }
  }
  const strip = await zlibDeflate(raw);
  const sw = new TextEncoder().encode(software + '\0');
  let off = 8;
  const stripOff = off; off += strip.length; off += off & 1;
  const bpsOff = off; off += spp * 2;
  const xresOff = off; off += 8;
  const yresOff = off; off += 8;
  const swOff = off; off += sw.length; off += off & 1;
  const ifdOff = off;
  const tags = [
    [256, 4, 1, w], [257, 4, 1, h], [258, 3, spp, bpsOff], [259, 3, 1, 8], [262, 3, 1, 2], [273, 4, 1, stripOff],
    [277, 3, 1, spp], [278, 4, 1, h], [279, 4, 1, strip.length], [282, 5, 1, xresOff], [283, 5, 1, yresOff],
    [284, 3, 1, 1], [296, 3, 1, 2], [305, 2, sw.length, swOff], [317, 3, 1, 2],
  ];
  if (alpha) tags.push([338, 3, 1, 2]);
  const total = ifdOff + 2 + tags.length * 12 + 4;
  const buf = new ArrayBuffer(total);
  const dv = new DataView(buf), u8 = new Uint8Array(buf);
  dv.setUint16(0, 0x4949); dv.setUint16(2, 42, true); dv.setUint32(4, ifdOff, true);
  u8.set(strip, stripOff);
  for (let i = 0; i < spp; i++) dv.setUint16(bpsOff + 2 * i, 8, true);
  dv.setUint32(xresOff, Math.round(dpi * 100), true); dv.setUint32(xresOff + 4, 100, true);
  dv.setUint32(yresOff, Math.round(dpi * 100), true); dv.setUint32(yresOff + 4, 100, true);
  u8.set(sw, swOff);
  dv.setUint16(ifdOff, tags.length, true);
  tags.forEach(([tag, type, count, value], i) => {
    const o = ifdOff + 2 + i * 12;
    dv.setUint16(o, tag, true); dv.setUint16(o + 2, type, true); dv.setUint32(o + 4, count, true);
    if (type === 3 && count === 1) dv.setUint16(o + 8, value, true); else dv.setUint32(o + 8, value, true);
  });
  dv.setUint32(ifdOff + 2 + tags.length * 12, 0, true);
  return new Blob([buf], { type: 'image/tiff' });
}

// ------------------------------------------------------------------ PNG pHYs (DPI)
export async function pngWithDpi(blob, dpi) {
  const src = new Uint8Array(await blob.arrayBuffer());
  const ppm = Math.round(dpi / 0.0254);
  const chunk = new Uint8Array(21);
  const dv = new DataView(chunk.buffer);
  dv.setUint32(0, 9);
  chunk.set([0x70, 0x48, 0x59, 0x73], 4);           // 'pHYs'
  dv.setUint32(8, ppm); dv.setUint32(12, ppm); chunk[16] = 1;
  dv.setUint32(17, crc32(chunk.subarray(4, 17)));
  const at = 8 + 25;                                 // after signature + IHDR
  const out = new Uint8Array(src.length + 21);
  out.set(src.subarray(0, at), 0); out.set(chunk, at); out.set(src.subarray(at), at + 21);
  return new Blob([out], { type: 'image/png' });
}

// ------------------------------------------------------------------ ZIP (stored)
export class ZipWriter {
  constructor() { this.parts = []; this.central = []; this.offset = 0; }
  add(name, bytes) {
    const nm = new TextEncoder().encode(name);
    const crc = crc32(bytes);
    const h = new DataView(new ArrayBuffer(30));
    h.setUint32(0, 0x04034b50, true); h.setUint16(4, 20, true); h.setUint16(6, 0x0800, true); h.setUint16(8, 0, true);
    h.setUint16(10, 0, true); h.setUint16(12, 0x21, true); h.setUint32(14, crc, true);
    h.setUint32(18, bytes.length, true); h.setUint32(22, bytes.length, true); h.setUint16(26, nm.length, true); h.setUint16(28, 0, true);
    this.parts.push(new Uint8Array(h.buffer), nm, bytes);
    const c = new DataView(new ArrayBuffer(46));
    c.setUint32(0, 0x02014b50, true); c.setUint16(4, 20, true); c.setUint16(6, 20, true); c.setUint16(8, 0x0800, true);
    c.setUint16(10, 0, true); c.setUint16(12, 0, true); c.setUint16(14, 0x21, true); c.setUint32(16, crc, true);
    c.setUint32(20, bytes.length, true); c.setUint32(24, bytes.length, true); c.setUint16(28, nm.length, true);
    c.setUint32(42, this.offset, true);
    this.central.push(new Uint8Array(c.buffer), nm);
    this.offset += 30 + nm.length + bytes.length;
  }
  blob() {
    const size = this.central.reduce((s, p) => s + p.length, 0);
    const e = new DataView(new ArrayBuffer(22));
    const n = this.central.length / 2;
    e.setUint32(0, 0x06054b50, true); e.setUint16(8, n, true); e.setUint16(10, n, true);
    e.setUint32(12, size, true); e.setUint32(16, this.offset, true);
    return new Blob([...this.parts, ...this.central, new Uint8Array(e.buffer)], { type: 'application/zip' });
  }
}

// ------------------------------------------------------------------ PDF
const WINANSI = { '€': 0x80, '‚': 0x82, '„': 0x84, '…': 0x85, '‘': 0x91, '’': 0x92, '“': 0x93, '”': 0x94, '•': 0x95, '–': 0x96, '—': 0x97, '™': 0x99 };
const ASCII_FALLBACK = { 'Δ': 'd', '′': "'", '″': '"', '−': '-', '→': '->', '←': '<-', '⁺': '+', '⁻': '-', 'ε': 'e', 'α': 'a', 'β': 'b', 'χ': 'chi', '≈': '~', '≤': '<=', '≥': '>=', '▼': 'v', '▲': '^', ' ': ' ', ' ': ' ' };
export function toWinAnsi(str) {
  const out = [];
  for (const ch of String(str)) {
    const code = ch.codePointAt(0);
    if (WINANSI[ch]) out.push(WINANSI[ch]);
    else if (code < 128 || (code >= 160 && code <= 255)) out.push(code);
    else if (ASCII_FALLBACK[ch]) for (const c of ASCII_FALLBACK[ch]) out.push(c.charCodeAt(0));
    else out.push(63);
  }
  return out;
}
function pdfString(str) {
  let s = '(';
  for (const b of toWinAnsi(str)) {
    if (b === 40 || b === 41 || b === 92) s += '\\' + String.fromCharCode(b);
    else if (b < 32 || b > 126) s += '\\' + b.toString(8).padStart(3, '0');
    else s += String.fromCharCode(b);
  }
  return s + ')';
}
const n3 = v => (Math.round(v * 1000) / 1000).toString();
// Symbol-font code points for Greek letters used in labels, and their advance widths (em)
const SYMBOL = { 'Δ': 'D', 'α': 'a', 'β': 'b', 'γ': 'g', 'δ': 'd', 'ε': 'e', 'χ': 'c', 'μ': 'm', 'σ': 's', 'τ': 't', 'κ': 'k', 'λ': 'l', 'θ': 'q', 'φ': 'f', 'ψ': 'y', 'ω': 'w' };
const SYMBOL_W = { D: 0.612, a: 0.631, b: 0.549, g: 0.411, d: 0.494, e: 0.439, c: 0.549, m: 0.576, s: 0.603, t: 0.439, k: 0.549, l: 0.549, q: 0.521, f: 0.521, y: 0.686, w: 0.686 };

let _measure = null;
function textWidth(str, size, bold) {
  if (typeof document === 'undefined') return str.length * size * 0.52;
  _measure = _measure || document.createElement('canvas').getContext('2d');
  _measure.font = `${bold ? 'bold ' : ''}${size}px "RCV Figure Sans", "Liberation Sans", Helvetica, sans-serif`;
  return _measure.measureText(str).width;
}

// Source of embedded fonts for PDF text: async (regularText, boldText) -> { regular, bold } (subsets of the bundled
// Liberation Sans from the local server) or null for the standard Helvetica fonts. Set once by the application.
let _fontProvider = null;
export function setPdfFontProvider(fn) { _fontProvider = fn; }
const WIN_SET = (() => { const s = new Set(); for (let c = 32; c < 127; c++) s.add(String.fromCharCode(c)); for (let c = 160; c < 256; c++) s.add(String.fromCharCode(c)); for (const ch of Object.keys(WINANSI)) s.add(ch); return s; })();
function b64bytes(b64) { const bin = atob(b64); const u = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return u; }
function pdfBytes(codes) {
  let s = '(';
  for (const b of codes) {
    if (b === 40 || b === 41 || b === 92) s += '\\' + String.fromCharCode(b);
    else if (b < 32 || b > 126) s += '\\' + b.toString(8).padStart(3, '0');
    else s += String.fromCharCode(b);
  }
  return s + ')';
}

export class PDFDoc {
  // Coordinates in the drawing API are top-left based, in points (1/72 inch).
  constructor(width, height, { title = 'RCV figure' } = {}) {
    this.W = width; this.H = height; this.title = title;
    this.ops = [];
    this.images = [];
    this.gstates = new Map();
  }
  _y(y) { return this.H - y; }
  _gs(fillA, strokeA) {
    const key = `${n3(fillA)}/${n3(strokeA)}`;
    if (!this.gstates.has(key)) this.gstates.set(key, { name: `GS${this.gstates.size + 1}`, ca: fillA, CA: strokeA });
    return this.gstates.get(key).name;
  }
  _style(s) {
    const o = [];
    if (s.fill) o.push(`${n3(s.fill.r)} ${n3(s.fill.g)} ${n3(s.fill.b)} rg`);
    if (s.stroke) o.push(`${n3(s.stroke.r)} ${n3(s.stroke.g)} ${n3(s.stroke.b)} RG`, `${n3(s.width ?? 1)} w`, '1 J 1 j');
    o.push(s.dash && s.dash.length ? `[${s.dash.map(n3).join(' ')}] 0 d` : '[] 0 d');
    const fa = s.fill ? (s.fill.a ?? 1) * (s.opacity ?? 1) : 1, sa = s.stroke ? (s.stroke.a ?? 1) * (s.strokeOpacity ?? s.opacity ?? 1) : 1;
    o.push(`/${this._gs(fa, sa)} gs`);
    return o.join(' ');
  }
  _paint(s) { return s.fill && s.stroke ? 'B' : s.fill ? 'f' : s.stroke ? 'S' : 'n'; }
  async image(img, x, y, w, h) {
    const { width, height, data } = img;
    const rgb = new Uint8Array(width * height * 3);
    const a = new Uint8Array(width * height);
    let hasAlpha = false;
    for (let i = 0, j = 0, k = 0; i < data.length; i += 4) {
      rgb[j++] = data[i]; rgb[j++] = data[i + 1]; rgb[j++] = data[i + 2];
      a[k++] = data[i + 3];
      if (data[i + 3] !== 255) hasAlpha = true;
    }
    const name = `Im${this.images.length + 1}`;
    this.images.push({ name, width, height, rgb: await zlibDeflate(rgb), alpha: hasAlpha ? await zlibDeflate(a) : null });
    this.ops.push(`q ${n3(w)} 0 0 ${n3(h)} ${n3(x)} ${n3(this._y(y + h))} cm /${name} Do Q`);
  }
  line(x1, y1, x2, y2, s) { this.ops.push(`q ${this._style({ ...s, fill: null })} ${n3(x1)} ${n3(this._y(y1))} m ${n3(x2)} ${n3(this._y(y2))} l S Q`); }
  rect(x, y, w, h, s) { this.ops.push(`q ${this._style(s)} ${n3(x)} ${n3(this._y(y + h))} ${n3(w)} ${n3(h)} re ${this._paint(s)} Q`); }
  circle(cx, cy, r, s) {
    const k = 0.5523 * r, Y = this._y(cy);
    this.ops.push(`q ${this._style(s)} ${n3(cx + r)} ${n3(Y)} m ${n3(cx + r)} ${n3(Y + k)} ${n3(cx + k)} ${n3(Y + r)} ${n3(cx)} ${n3(Y + r)} c ` +
      `${n3(cx - k)} ${n3(Y + r)} ${n3(cx - r)} ${n3(Y + k)} ${n3(cx - r)} ${n3(Y)} c ${n3(cx - r)} ${n3(Y - k)} ${n3(cx - k)} ${n3(Y - r)} ${n3(cx)} ${n3(Y - r)} c ` +
      `${n3(cx + k)} ${n3(Y - r)} ${n3(cx + r)} ${n3(Y - k)} ${n3(cx + r)} ${n3(Y)} c h ${this._paint(s)} Q`);
  }
  path(cmds, s) {
    const o = [];
    for (const c of cmds) {
      if (c[0] === 'M') o.push(`${n3(c[1])} ${n3(this._y(c[2]))} m`);
      else if (c[0] === 'L') o.push(`${n3(c[1])} ${n3(this._y(c[2]))} l`);
      else if (c[0] === 'C') o.push(`${n3(c[1])} ${n3(this._y(c[2]))} ${n3(c[3])} ${n3(this._y(c[4]))} ${n3(c[5])} ${n3(this._y(c[6]))} c`);
      else if (c[0] === 'Z') o.push('h');
    }
    this.ops.push(`q ${this._style(s)} ${o.join(' ')} ${this._paint(s)} Q`);
  }
  // Text is laid out when the file is written, once the fonts are known: embedded Liberation Sans subsets (F1 regular,
  // F2 bold, F3/F4 for characters outside WinAnsi such as Δ, ′, →) or the standard Helvetica / Symbol fonts.
  text(x, y, str, { size = 10, bold = false, color = { r: 0, g: 0, b: 0 }, anchor = 'start', opacity = 1, angle = 0 } = {}) {
    const s = String(str);
    const w = textWidth(s, size, bold);
    const dx = anchor === 'middle' ? -w / 2 : anchor === 'end' ? -w : 0;
    this.ops.push({ text: s, size, bold, x, y: this._y(y), dx, angle, color, gs: this._gs(opacity, 1) });
  }

  _textOp(o, fonts) {
    const runs = [];
    const push = (font, codes) => { if (runs.length && runs[runs.length - 1].font === font) runs[runs.length - 1].codes.push(...codes); else runs.push({ font, codes: [...codes] }); };
    for (const ch of o.text) {
      if (fonts) {
        const extra = (o.bold ? fonts.boldExtra : fonts.regExtra).get(ch);
        if (WIN_SET.has(ch)) push(o.bold ? 'F2' : 'F1', toWinAnsi(ch));
        else if (extra !== undefined) push(o.bold ? 'F4' : 'F3', [extra]);
        else push(o.bold ? 'F2' : 'F1', toWinAnsi(ch));
      } else {
        const sym = SYMBOL[ch];
        if (sym) push('F3', [sym.charCodeAt(0)]); else push(o.bold ? 'F2' : 'F1', toWinAnsi(ch));
      }
    }
    const body = runs.map(r => `/${r.font} ${n3(o.size)} Tf ${pdfBytes(r.codes)} Tj`).join(' ');
    // text matrix: rotation about the anchor point (SVG angles are clockwise, PDF counter-clockwise), then the anchor shift
    const a = -(o.angle || 0) * Math.PI / 180, c = Math.cos(a), s = Math.sin(a);
    const tx = o.x + c * (o.dx || 0), ty = o.y + s * (o.dx || 0);
    return `q /${o.gs} gs ${n3(o.color.r)} ${n3(o.color.g)} ${n3(o.color.b)} rg BT ${n3(c)} ${n3(s)} ${n3(-s)} ${n3(c)} ${n3(tx)} ${n3(ty)} Tm ${body} ET Q`;
  }

  async _fonts() {
    const provider = this.fontProvider === undefined ? _fontProvider : this.fontProvider;
    if (!provider) return null;
    let reg = '', bold = '';
    for (const o of this.ops) if (typeof o === 'object') { if (o.bold) bold += o.text; else reg += o.text; }
    if (!reg && !bold) return null;
    const uniq = s => [...new Set(s)].join('');
    let r = null;
    try { r = await provider(uniq(reg), uniq(bold)); } catch { r = null; }
    if (!r || !r.regular) return null;
    const extraMap = f => new Map((f?.extra || []).map((e, i) => [e.char, 33 + i]));
    return { reg: r.regular, bold: r.bold || r.regular, regExtra: extraMap(r.regular), boldExtra: extraMap(r.bold || r.regular) };
  }

  async blob() {
    const enc = new TextEncoder();
    const chunks = [];
    const offsets = [];
    let pos = 0;
    const push = b => { const u = typeof b === 'string' ? enc.encode(b) : b; chunks.push(u); pos += u.length; };
    const objs = [];             // [id, header-string, streamBytes|null]
    let id = 0;
    const next = () => ++id;
    const fonts = await this._fonts();
    const catalog = next(), pages = next(), page = next(), content = next(), f1 = next(), f2 = next(), f3 = next(), info = next();
    const f4 = fonts ? next() : null;
    const imgIds = this.images.map(im => ({ im, id: next(), smask: im.alpha ? next() : null }));
    const ops = this.ops.map(o => (typeof o === 'object' ? this._textOp(o, fonts) : o));
    const contentBytes = await zlibDeflate(enc.encode(ops.join('\n')));
    const xobj = imgIds.map(x => `/${x.im.name} ${x.id} 0 R`).join(' ');
    const gs = [...this.gstates.values()].map(g => `/${g.name} << /Type /ExtGState /ca ${n3(g.ca)} /CA ${n3(g.CA)} >>`).join(' ');
    objs.push([catalog, `<< /Type /Catalog /Pages ${pages} 0 R >>`]);
    objs.push([pages, `<< /Type /Pages /Kids [${page} 0 R] /Count 1 >>`]);
    objs.push([page, `<< /Type /Page /Parent ${pages} 0 R /MediaBox [0 0 ${n3(this.W)} ${n3(this.H)}] /Contents ${content} 0 R ` +
      `/Resources << /Font << /F1 ${f1} 0 R /F2 ${f2} 0 R /F3 ${f3} 0 R${f4 ? ` /F4 ${f4} 0 R` : ''} >> /XObject << ${xobj} >> /ExtGState << ${gs} >> >> >>`]);
    objs.push([content, `<< /Length ${contentBytes.length} /Filter /FlateDecode >>`, contentBytes]);
    if (fonts) {
      // embedded font subsets: one font program per style, shared by its WinAnsi font and its extra-character font
      const tag = () => Array.from({ length: 6 }, () => String.fromCharCode(65 + Math.floor(Math.random() * 26))).join('');
      const embed = async (f) => {
        const raw = b64bytes(f.ttf), packed = await zlibDeflate(raw);
        const file = next(), desc = next(), name = `${tag()}+${f.name}`;
        objs.push([file, `<< /Length ${packed.length} /Length1 ${raw.length} /Filter /FlateDecode >>`, packed]);
        objs.push([desc, `<< /Type /FontDescriptor /FontName /${name} /Flags 32 /FontBBox [${f.bbox.join(' ')}] /ItalicAngle ${f.italicAngle} ` +
          `/Ascent ${f.ascent} /Descent ${f.descent} /CapHeight ${f.capHeight} /StemV ${f.stemV} /FontFile2 ${file} 0 R >>`]);
        return { desc, name };
      };
      const R = await embed(fonts.reg), B = fonts.bold === fonts.reg ? R : await embed(fonts.bold);
      const win = (id, f, d) => objs.push([id, `<< /Type /Font /Subtype /TrueType /BaseFont /${d.name} /FirstChar 32 /LastChar 255 /Widths [${f.widths.join(' ')}] /Encoding /WinAnsiEncoding /FontDescriptor ${d.desc} 0 R >>`]);
      const extra = (id, f, d) => {
        const ex = f.extra || [];
        if (!ex.length) { objs.push([id, `<< /Type /Font /Subtype /TrueType /BaseFont /${d.name} /FirstChar 32 /LastChar 32 /Widths [278] /Encoding /WinAnsiEncoding /FontDescriptor ${d.desc} 0 R >>`]); return; }
        objs.push([id, `<< /Type /Font /Subtype /TrueType /BaseFont /${d.name} /FirstChar 33 /LastChar ${32 + ex.length} /Widths [${ex.map(e => e.width).join(' ')}] ` +
          `/Encoding << /Type /Encoding /BaseEncoding /WinAnsiEncoding /Differences [33 ${ex.map(e => '/' + e.glyph).join(' ')}] >> /FontDescriptor ${d.desc} 0 R >>`]);
      };
      win(f1, fonts.reg, R); win(f2, fonts.bold, B); extra(f3, fonts.reg, R); extra(f4, fonts.bold, B);
    } else {
      objs.push([f1, '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>']);
      objs.push([f2, '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>']);
      objs.push([f3, '<< /Type /Font /Subtype /Type1 /BaseFont /Symbol >>']);
    }
    const now = new Date();
    const d = `D:${now.getFullYear()}${String(now.getMonth() + 1).padStart(2, '0')}${String(now.getDate()).padStart(2, '0')}${String(now.getHours()).padStart(2, '0')}${String(now.getMinutes()).padStart(2, '0')}00`;
    objs.push([info, `<< /Title ${pdfString(this.title)} /Producer (RCV - RNA Conformation Viewer) /CreationDate (${d}) >>`]);
    for (const x of imgIds) {
      const sm = x.smask ? ` /SMask ${x.smask} 0 R` : '';
      objs.push([x.id, `<< /Type /XObject /Subtype /Image /Width ${x.im.width} /Height ${x.im.height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /FlateDecode /Length ${x.im.rgb.length}${sm} >>`, x.im.rgb]);
      if (x.smask) objs.push([x.smask, `<< /Type /XObject /Subtype /Image /Width ${x.im.width} /Height ${x.im.height} /ColorSpace /DeviceGray /BitsPerComponent 8 /Filter /FlateDecode /Length ${x.im.alpha.length} >>`, x.im.alpha]);
    }
    objs.sort((a, b) => a[0] - b[0]);
    push(new Uint8Array([...enc.encode('%PDF-1.4\n%'), 0xe2, 0xe3, 0xcf, 0xd3, 0x0a]));   // binary marker comment
    for (const [oid, head, stream] of objs) {
      offsets[oid] = pos;
      push(`${oid} 0 obj\n${head}\n`);
      if (stream) { push('stream\n'); push(stream); push('\nendstream\n'); }
      push('endobj\n');
    }
    const xref = pos;
    let x = `xref\n0 ${id + 1}\n0000000000 65535 f \n`;
    for (let i = 1; i <= id; i++) x += `${String(offsets[i]).padStart(10, '0')} 00000 n \n`;
    push(x);
    push(`trailer\n<< /Size ${id + 1} /Root ${catalog} 0 R /Info ${info} 0 R >>\nstartxref\n${xref}\n%%EOF\n`);
    return new Blob(chunks, { type: 'application/pdf' });
  }
}

// Walker backend that draws charts into a PDFDoc region (scale s, offset ox/oy, in points).
export function pdfBackend(doc, parseColor, ox, oy, s) {
  const st = (style, filled = true) => ({
    fill: filled ? parseColor(style.fill) : null, stroke: parseColor(style.stroke), width: style.strokeWidth * s,
    dash: style.dash ? style.dash.map(v => v * s) : null, opacity: style.opacity ?? 1, strokeOpacity: style.strokeOpacity ?? 1,
  });
  const X = v => ox + v * s, Y = v => oy + v * s;
  return {
    line: (x1, y1, x2, y2, style) => { const t = st(style, false); if (t.stroke) doc.line(X(x1), Y(y1), X(x2), Y(y2), t); },
    circle: (cx, cy, r, style) => { const t = st(style); if (t.fill || t.stroke) doc.circle(X(cx), Y(cy), r * s, t); },
    rect: (x, y, w, h, rx, style) => { const t = st(style); if (t.fill || t.stroke) doc.rect(X(x), Y(y), w * s, h * s, t); },
    path: (cmds, style) => {
      const t = st(style);
      if (!t.fill && !t.stroke) return;
      doc.path(cmds.map(c => (c[0] === 'Z' ? c : [c[0], ...c.slice(1).map((v, i) => (i % 2 === 0 ? X(v) : Y(v)))])), t);
    },
    text: (x, y, str, style) => {
      const c = parseColor(style.fill);
      if (!c) return;
      doc.text(X(x), Y(y), str, { size: style.fontSize * s, bold: Number(style.fontWeight) >= 600 || style.fontWeight === 'bold', color: c, anchor: style.anchor, opacity: c.a * (style.opacity ?? 1), angle: style.angle || 0 });
    },
  };
}
