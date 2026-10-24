// Multi-panel publication figure: several conformations in the same orientation, lettered A, B, C …, each with
// its name and energy, an optional 2D-structure inset, and one shared colour legend. PNG, TIFF or PDF.
import { encodeTIFF, pngWithDpi, PDFDoc } from './encoders.js';
import { drawLegend } from './video.js';
import { fmt, FIG_FONT, ensureFigureFont } from '../util.js';

export function defaultPanelConfs(app, max = 4) {
  const ens = app.ens;
  const seen = new Set(), out = [];
  for (const k of app.order) {                       // one representative per 2D structure first
    const key = ens.conformations[k].dbn;
    if (!seen.has(key)) { seen.add(key); out.push(k); }
    if (out.length >= max) break;
  }
  for (const k of app.order) { if (out.length >= Math.min(max, ens.n)) break; if (!out.includes(k)) out.push(k); }
  return out;
}

export async function composePanels(app, o) {
  await ensureFigureFont();
  const ens = app.ens, v = app.viewer;
  const n = o.confs.length;
  if (!n) throw new Error('choose at least one conformation');
  const cols = o.cols || (n <= 3 ? n : n === 4 ? 2 : 3);
  const rows = Math.ceil(n / cols);
  const pw = Math.round(o.width / cols), ph = Math.round(pw / (o.aspect || 1));
  const s = pw / 1000;                                   // type scale relative to a 1000 px panel
  const legendH = o.legend ? Math.round(150 * s * Math.min(1.6, Math.max(1, cols / 2))) : 0;
  const W = pw * cols, H = ph * rows + legendH;
  const canvas = document.createElement('canvas');
  canvas.width = W; canvas.height = H;
  const ctx = canvas.getContext('2d');
  const bg = o.transparent ? null : (o.background || v.settings.background);
  if (bg) { ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H); }
  const dark = bg ? app.isDarkColor(bg) : false;
  const ink = dark ? '#f0efec' : '#141414', sub = dark ? '#c3c2b7' : '#52514e';
  const saved = app.snapshotView();
  const key = app.state.metric, m = ens.metricByKey[key] || { label: key, unit: '' };
  try {
    for (let i = 0; i < n; i++) {
      const k = o.confs[i];
      const x = (i % cols) * pw, y = Math.floor(i / cols) * ph;
      app.display.set(ens.coords[k]);
      v.setCoords(app.display);
      v.setResidueColorsLinear(app.colorsFor(k, o.confs[i - 1] ?? null).lin);
      // with a 2D inset, the 3D model gets the left part of the panel and the 2D drawing its own column
      const w3 = o.inset2d ? Math.round(pw * 0.7) : pw;
      const r = v.renderImage({ width: w3, height: ph, transparent: !bg, background: bg });
      ctx.drawImage(r.canvas, x, y);
      if (o.annotations !== false && app.annotations) { ctx.save(); ctx.translate(x, y); app.annotations.drawTo(ctx, r.project, w3, ph, Math.max(0.25, r.scale), ens.coords[k]); ctx.restore(); }
      const font = (px, w = 600) => `${w >= 600 ? 'bold' : 'normal'} ${Math.round(px * s)}px ${FIG_FONT}`;
      ctx.textBaseline = 'alphabetic';
      if (o.letters) { ctx.fillStyle = ink; ctx.font = font(64, 750); ctx.fillText(String.fromCharCode(65 + i), x + 34 * s, y + 84 * s); }
      if (o.captions) {
        ctx.fillStyle = ink; ctx.font = font(34, 650);
        const tx = x + (o.letters ? 104 : 34) * s;
        ctx.fillText(ens.shortId(k), tx, y + 70 * s);
        const val = ens.value(k, key);
        ctx.fillStyle = sub; ctx.font = font(27, 500);
        const line2 = [val !== null ? `${m.label} ${fmt(val)}${m.unit ? ' ' + m.unit : ''}` : null, `${ens.conformations[k].pairs.length} bp`].filter(Boolean).join(' · ');
        ctx.fillText(line2, tx, y + 106 * s);
      }
      if (o.inset2d && app.ss2d) {
        const iw = pw - w3 - 20 * s, ih = Math.min(ph - 150 * s, iw * 2.2);
        app.ss2d.drawTo(ctx, [x + w3, y + (ph - ih) / 2 + 40 * s, iw, ih], { k, dark, frame: false, unit: 2.6 * s });
      }
      if (o.frame) { ctx.strokeStyle = dark ? 'rgba(255,255,255,0.18)' : 'rgba(0,0,0,0.14)'; ctx.lineWidth = Math.max(1, 2 * s); ctx.strokeRect(x + 1, y + 1, pw - 2, ph - 2); }
    }
    if (o.legend) {
      const lg = app.legendData();
      if (lg) drawLegend(ctx, lg, 34 * s, H - 30 * s, s * 1.3, ink, sub);
    }
  } finally {
    app.restoreView(saved);
    v.setResidueColorsLinear(app.curColors || app.colorsFor(app.state.conf, app.state.prev).lin);
    v.needsRender = true;
  }
  return canvas;
}

export async function exportPanels(app, o) {
  const canvas = await composePanels(app, o);
  const img = () => canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height);
  if (o.format === 'tiff') return { blob: await encodeTIFF(img(), { dpi: o.dpi, alpha: !!o.transparent }), ext: 'tif', canvas };
  if (o.format === 'pdf') {
    const pt = 72 / o.dpi;
    const doc = new PDFDoc(canvas.width * pt, canvas.height * pt, { title: `${app.ens.name} — ${o.confs.length} conformations` });
    await doc.image(img(), 0, 0, canvas.width * pt, canvas.height * pt);
    return { blob: await doc.blob(), ext: 'pdf', canvas };
  }
  const png = await new Promise(res => canvas.toBlob(res, 'image/png'));
  return { blob: await pngWithDpi(png, o.dpi), ext: 'png', canvas };
}
