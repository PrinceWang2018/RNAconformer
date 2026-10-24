// Publication figure export: high-resolution structure image (+ legend, labels, overlays) and an optional
// energy panel (profile, landscape, base-pair arcs) drawn as vectors in PDF or rasterised in PNG/TIFF.
import { FIG_FONT, ensureFigureFont } from '../util.js';
import { encodeTIFF, pngWithDpi, PDFDoc, pdfBackend } from './encoders.js';
import { walkSVG, canvasBackend, parseColor } from '../panels/charts.js';
import { drawLegend } from './video.js';

function drawLabels(ctx, app, r, background) {
  const v = app.viewer;
  const list = v.labelResidues();
  const s = Math.max(1, r.scale);
  ctx.save();
  ctx.font = `bold ${Math.round(12 * s)}px ${FIG_FONT}`;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  for (const res of list) {
    const [x, y, z] = r.project(...v.labelAnchor(res));
    if (z < -1 || z > 1) continue;
    const t = v.labelText ? v.labelText(res) : app.ens.resLabel(res);
    const w = ctx.measureText(t).width + 12 * s, h = 18 * s;
    const yy = y - h * 0.9;
    ctx.fillStyle = 'rgba(255,255,255,0.92)';
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(x - w / 2, yy - h / 2, w, h, 4 * s); else ctx.rect(x - w / 2, yy - h / 2, w, h);
    ctx.fill();
    ctx.fillStyle = '#111111';
    ctx.fillText(t, x, yy + 0.5);
  }
  ctx.restore();
}

function drawOverlayPrims(ctx, app, r) {
  const prims = app.viewer.overlayPrimitives();
  const s = Math.max(1, r.scale);
  ctx.save();
  for (const p of prims) {
    if (p.type === 'line') {
      const [x1, y1] = r.project(...p.a), [x2, y2] = r.project(...p.b);
      ctx.lineCap = 'round';
      ctx.strokeStyle = 'rgba(255,255,255,0.8)'; ctx.lineWidth = ((p.width || 2) + 3) * s; ctx.setLineDash([]);
      ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
      ctx.strokeStyle = p.color || '#111'; ctx.lineWidth = (p.width || 2) * s; ctx.setLineDash(p.dash ? [6 * s, 4 * s] : []);
      ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
      if (p.label) {
        ctx.setLineDash([]);
        ctx.font = `bold ${Math.round(12 * s)}px ${FIG_FONT}`;
        ctx.textAlign = 'center';
        ctx.lineWidth = 3 * s; ctx.strokeStyle = '#fff'; ctx.fillStyle = '#111';
        ctx.strokeText(p.label, (x1 + x2) / 2, (y1 + y2) / 2 - 8 * s);
        ctx.fillText(p.label, (x1 + x2) / 2, (y1 + y2) / 2 - 8 * s);
      }
    } else if (p.type === 'point') {
      const [x, y] = r.project(...p.p);
      ctx.beginPath(); ctx.arc(x, y, (p.r || 4) * s, 0, Math.PI * 2);
      ctx.fillStyle = p.color || '#111'; ctx.fill();
    }
  }
  ctx.restore();
}

// Panel blocks from the live energy panel (SVG) — measured in CSS px of the live charts.
function panelBlocks(app) {
  const ep = app.energyPanel;
  const ens = app.ens, m = ens.metricByKey[app.state.metric] || { label: app.state.metric };
  const blocks = [];
  const add = (title, node) => { if (node && node.getAttribute('viewBox')) { const vb = node.getAttribute('viewBox').split(' ').map(Number); blocks.push({ title, node, w: vb[2], h: vb[3] }); } };
  add(`${m.label}${m.unit ? ' (' + m.unit + ')' : ''} per conformation`, ep.profile);
  add(ep.landMode === 'map' ? 'Similarity map (MDS of C4′ RMSD)' : `${m.label} vs C4′ RMSD`, ep.land);
  add(`Base pairs: ${ens.shortId(app.state.conf)} (above) vs ${ens.shortId(app.state.prev ?? app.state.conf)} (below)`, ep.arcs.svg);
  const W = Math.max(...blocks.map(b => b.w));
  const gap = 26, titleH = 18;
  let H = 0;
  for (const b of blocks) { b.y = H + titleH; H += titleH + b.h + gap; }
  return { blocks, W, H: H - gap + 6 };
}

function withLightTheme(fn) {
  const root = document.documentElement;
  const prev = root.getAttribute('data-theme');
  root.setAttribute('data-theme', 'light');
  try { return fn(); } finally { if (prev) root.setAttribute('data-theme', prev); else root.removeAttribute('data-theme'); }
}

export async function composeFigure(app, o) {
  await ensureFigureFont();
  const r = app.viewer.renderImage({ width: o.width, height: o.height, transparent: o.transparent, background: o.background });
  const ctx = r.canvas.getContext('2d');
  if (o.labels) drawLabels(ctx, app, r, o.background);
  drawOverlayPrims(ctx, app, r);
  if (o.annotations !== false) app.annotations?.drawTo(ctx, r.project, o.width, o.height, Math.max(0.25, r.scale));
  if (o.legend) {
    const lg = app.legendData();
    const dark = !o.transparent && app.isDarkColor(o.background || app.viewer.settings.background);
    const s = o.height / 1080 * 1.1;
    if (lg) drawLegend(ctx, lg, 40 * s, o.height - 36 * s, s, dark ? '#f0efec' : '#141414', dark ? '#c3c2b7' : '#52514e');
  }
  return r.canvas;
}

export async function exportFigure(app, o) {
  const structure = await composeFigure(app, o);
  const panel = o.panel && app.ens ? panelBlocks(app) : null;
  const pxPerPt = o.dpi / 72;
  if (o.format === 'pdf') {
    const imgWpt = o.width / pxPerPt, imgHpt = o.height / pxPerPt;
    let pageW = imgWpt, s = 0;
    if (panel) { s = (imgHpt * 0.94) / panel.H; pageW += panel.W * s + 18; }
    const doc = new PDFDoc(pageW, imgHpt, { title: `${app.ens.name} — ${app.ens.shortId(app.state.conf)}` });
    const img = structure.getContext('2d').getImageData(0, 0, structure.width, structure.height);
    await doc.image(img, 0, 0, imgWpt, imgHpt);
    if (panel) {
      const ox = imgWpt + 12, oy = imgHpt * 0.03;
      withLightTheme(() => {
        for (const b of panel.blocks) {
          doc.text(ox, oy + (b.y - 6) * s, b.title, { size: 9.5 * s, bold: true, color: { r: 0.2, g: 0.2, b: 0.2 } });
          walkSVG(b.node, pdfBackend(doc, parseColor, ox, oy + b.y * s, s));
        }
      });
    }
    return { blob: await doc.blob(), ext: 'pdf' };
  }
  // raster formats
  let canvas = structure;
  if (panel) {
    const s = (o.height * 0.94) / panel.H;
    const out = document.createElement('canvas');
    out.width = Math.round(o.width + panel.W * s + 30 * (o.height / 1000));
    out.height = o.height;
    const c = out.getContext('2d');
    if (!o.transparent) { c.fillStyle = o.background || '#ffffff'; c.fillRect(0, 0, out.width, out.height); }
    c.drawImage(structure, 0, 0);
    const ox = o.width + 18 * (o.height / 1000), oy = o.height * 0.03;
    withLightTheme(() => {
      for (const b of panel.blocks) {
        c.save();
        c.fillStyle = '#333';
        c.font = `bold ${9.5 * s}px ${FIG_FONT}`;
        c.fillText(b.title, ox, oy + (b.y - 6) * s);
        c.translate(ox, oy + b.y * s);
        walkSVG(b.node, canvasBackend(c, s));
        c.restore();
      }
    });
    canvas = out;
  }
  if (o.format === 'tiff') {
    const img = canvas.getContext('2d').getImageData(0, 0, canvas.width, canvas.height);
    return { blob: await encodeTIFF(img, { dpi: o.dpi, alpha: !!o.transparent }), ext: 'tif' };
  }
  const png = await new Promise(res => canvas.toBlob(res, 'image/png'));
  return { blob: await pngWithDpi(png, o.dpi), ext: 'png' };
}
