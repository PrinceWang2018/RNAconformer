// Export of any chart drawn as SVG (right-panel charts, plot studio): standalone SVG with resolved styles and the
// figure font, vector PDF (embedded font subsets), PNG and TIFF at a chosen resolution, or a copy to the clipboard.
import { walkSVG, canvasBackend, parseColor } from '../panels/charts.js';
import { PDFDoc, pdfBackend, encodeTIFF, pngWithDpi } from './encoders.js';
import { saveBlob, slug, ensureFigureFont } from '../util.js';

const SVG_FONT = '"Liberation Sans", Arimo, "Helvetica Neue", Helvetica, sans-serif';
const PROPS = ['fill', 'stroke', 'stroke-width', 'stroke-dasharray', 'stroke-linecap', 'stroke-linejoin', 'opacity', 'fill-opacity', 'stroke-opacity', 'font-size', 'font-weight', 'text-anchor'];

function size(svgEl) {
  const vb = (svgEl.getAttribute('viewBox') || '').split(/[ ,]+/).map(Number);
  if (vb.length === 4 && vb[2] > 0) return { w: vb[2], h: vb[3], x: vb[0], y: vb[1] };
  const r = svgEl.getBoundingClientRect();
  return { w: r.width || 300, h: r.height || 150, x: 0, y: 0 };
}

// A self-contained SVG string: computed colours and sizes copied onto every element, no CSS variables or classes,
// invisible hit areas removed, the figure font named (editors substitute it if not installed).
export function standaloneSVG(svgEl, { title = '' } = {}) {
  const { w, h, x, y } = size(svgEl);
  const clone = svgEl.cloneNode(true);
  const src = [svgEl, ...svgEl.querySelectorAll('*')], dst = [clone, ...clone.querySelectorAll('*')];
  const drop = [];
  src.forEach((node, i) => {
    const out = dst[i];
    const cs = getComputedStyle(node);
    if (cs.display === 'none' || cs.visibility === 'hidden') { drop.push(out); return; }
    const tag = node.tagName.toLowerCase();
    const invisible = (cs.fill === 'none' || /rgba\([^)]*,\s*0\)/.test(cs.fill) || cs.fill === 'transparent') &&
      (cs.stroke === 'none' || /rgba\([^)]*,\s*0\)/.test(cs.stroke) || cs.stroke === 'transparent');
    if (['path', 'rect', 'line', 'circle', 'polyline'].includes(tag) && (invisible || node.classList.contains('hit'))) { drop.push(out); return; }
    for (const p of PROPS) {
      const v = cs.getPropertyValue(p);
      if (v && v !== 'normal' && !(p === 'opacity' && v === '1')) out.setAttribute(p, v);
    }
    if (tag === 'text') out.setAttribute('font-family', SVG_FONT);
    out.removeAttribute('class');
    out.removeAttribute('style');
    for (const a of [...out.attributes]) if (a.name.startsWith('data-') || a.name.startsWith('aria-') || a.name === 'role' || a.name === 'tabindex') out.removeAttribute(a.name);
  });
  for (const d of drop) d.remove();
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('viewBox', `${x} ${y} ${w} ${h}`);
  clone.setAttribute('width', `${w}pt`);
  clone.setAttribute('height', `${h}pt`);
  clone.setAttribute('font-family', SVG_FONT);
  if (title) { const t = document.createElementNS('http://www.w3.org/2000/svg', 'title'); t.textContent = title; clone.prepend(t); }
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(clone);
}

// Vector PDF, page = chart size in points (1 SVG unit = 1 pt) times `scale`.
export async function chartPDF(svgEl, { scale = 1, title = 'RCV chart', background = '#ffffff' } = {}) {
  await ensureFigureFont();
  const { w, h, x, y } = size(svgEl);
  const doc = new PDFDoc(w * scale, h * scale, { title });
  if (background) { const c = parseColor(hexToRgbStr(background)); if (c) doc.rect(0, 0, w * scale, h * scale, { fill: c }); }
  walkSVG(svgEl, pdfBackend(doc, parseColor, -x * scale, -y * scale, scale));
  return doc.blob();
}

// Raster: `dpi` with 1 SVG unit = 1 pt.
export async function chartCanvas(svgEl, { dpi = 600, background = '#ffffff' } = {}) {
  await ensureFigureFont();
  const { w, h, x, y } = size(svgEl);
  const s = dpi / 72;
  const c = document.createElement('canvas');
  c.width = Math.max(1, Math.round(w * s)); c.height = Math.max(1, Math.round(h * s));
  const ctx = c.getContext('2d');
  if (background) { ctx.fillStyle = background; ctx.fillRect(0, 0, c.width, c.height); }
  ctx.translate(-x * s, -y * s);
  walkSVG(svgEl, canvasBackend(ctx, s));
  return c;
}

export async function exportChart(svgEl, format, name, { dpi = 600, title = '', background = '#ffffff' } = {}) {
  let blob, ext = format;
  if (format === 'svg') blob = new Blob([standaloneSVG(svgEl, { title })], { type: 'image/svg+xml' });
  else if (format === 'pdf') blob = await chartPDF(svgEl, { title: title || name, background });
  else {
    const c = await chartCanvas(svgEl, { dpi, background: format === 'png-transparent' ? null : background });
    if (format === 'tiff') { blob = await encodeTIFF(c.getContext('2d').getImageData(0, 0, c.width, c.height), { dpi }); ext = 'tif'; }
    else { blob = await pngWithDpi(await new Promise(r => c.toBlob(r, 'image/png')), dpi); ext = 'png'; }
  }
  return saveBlob(blob, `${slug(name)}.${ext}`);
}

export async function copyChart(svgEl, dpi = 300) {
  const c = await chartCanvas(svgEl, { dpi });
  const blob = await new Promise(r => c.toBlob(r, 'image/png'));
  await navigator.clipboard.write([new ClipboardItem({ 'image/png': blob })]);
}

function hexToRgbStr(hex) {
  const n = parseInt(String(hex).replace('#', ''), 16);
  return `rgb(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255})`;
}

// Context-menu items for a chart.
export function chartMenuItems(svgEl, name, toast) {
  const run = fmt => async () => {
    try { const saved = await exportChart(svgEl, fmt, name, { title: name }); if (saved) toast(saved, 'ok', 'Chart exported'); }
    catch (err) { toast(err.message, 'error', 'Export failed'); }
  };
  return [
    { head: name },
    { label: 'Export as PDF (vector)…', action: run('pdf') },
    { label: 'Export as SVG (vector)…', action: run('svg') },
    { label: 'Export as PNG (600 dpi)…', action: run('png') },
    { label: 'Export as TIFF (600 dpi)…', action: run('tiff') },
    { sep: true },
    { label: 'Copy image', action: async () => { try { await copyChart(svgEl); toast('Chart copied as an image.', 'ok'); } catch (err) { toast('Copying images is not available here; export it instead.', 'error'); } } },
  ];
}
