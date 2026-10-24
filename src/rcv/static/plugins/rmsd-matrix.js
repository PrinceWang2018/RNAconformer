// Built-in plugin: pairwise C4′ RMSD heatmap of the ensemble. Hover a cell to compare two
// conformations; click it to morph from the row conformation to the column conformation.
const RAMP = ['#f3f8fe', '#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b'];

function rampColor(t) {
  t = Math.max(0, Math.min(1, t)) * (RAMP.length - 1);
  const i = Math.min(RAMP.length - 2, Math.floor(t)), f = t - i;
  const a = parseInt(RAMP[i].slice(1), 16), b = parseInt(RAMP[i + 1].slice(1), 16);
  const ch = s => Math.round(((a >> s) & 255) * (1 - f) + ((b >> s) & 255) * f);
  return `rgb(${ch(16)},${ch(8)},${ch(0)})`;
}

export default {
  id: 'rmsd-matrix',
  name: 'Similarity matrix',
  version: '1.0',
  description: 'Heatmap of pairwise C4′ RMSD between conformations; click a cell to morph between them.',
  activate(api) {
    const { el, svg, fmt } = api.util;
    let root = null;
    const render = () => {
      if (!root) return;
      root.replaceChildren();
      const ens = api.ensemble;
      root.append(el('div', { class: 'section-title', text: 'Similarity matrix (C4′ RMSD)' }));
      if (!ens) { root.append(el('div', { class: 'desc', text: 'No ensemble loaded.' })); return; }
      const order = api.order, n = order.length;
      const max = Math.max(1e-6, ...ens.rmsd.flat());
      const size = Math.max(200, Math.min(root.clientWidth - 8, 340));
      const pad = 30, cell = (size - pad) / n;
      const s = svg('svg', { viewBox: `0 0 ${size} ${size}`, width: '100%', class: 'chart', role: 'img', 'aria-label': 'Pairwise RMSD heatmap' });
      const cur = order.indexOf(api.conformation);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
        const v = ens.rmsd[order[i]][order[j]];
        const r = svg('rect', { x: pad + j * cell, y: pad + i * cell, width: cell + 0.3, height: cell + 0.3, fill: rampColor(v / max) });
        r.style.cursor = 'pointer';
        r.addEventListener('pointermove', e => api.ui.tooltip.show([api.ui.tipRow(`${fmt(v)} Å`, `${ens.shortId(order[i])} vs ${ens.shortId(order[j])}`),
          el('div', { class: 'tl', text: 'click: morph from row to column' })], e.clientX, e.clientY));
        r.addEventListener('pointerleave', () => api.ui.tooltip.hide());
        r.addEventListener('click', () => {
          api.goTo(order[i], { animate: false });
          requestAnimationFrame(() => api.goTo(order[j]));
        });
        s.append(r);
      }
      if (cur >= 0) {
        s.append(svg('rect', { x: pad, y: pad + cur * cell, width: n * cell, height: cell, fill: 'none', stroke: 'var(--text)', 'stroke-width': 1.5 }));
        s.append(svg('rect', { x: pad + cur * cell, y: pad, width: cell, height: n * cell, fill: 'none', stroke: 'var(--text)', 'stroke-width': 1.5 }));
      }
      const every = Math.max(1, Math.ceil(n / 10));
      for (let i = 0; i < n; i += every) {
        s.append(svg('text', { x: pad - 4, y: pad + (i + 0.5) * cell + 3.5, 'text-anchor': 'end', text: String(i + 1) }));
        s.append(svg('text', { x: pad + (i + 0.5) * cell, y: pad - 6, 'text-anchor': 'middle', text: String(i + 1) }));
      }
      root.append(s);
      const lg = el('div', { style: { display: 'flex', alignItems: 'center', gap: '8px', marginTop: '8px', fontSize: '11.5px', color: 'var(--text-2)' } },
        el('span', { text: '0 Å' }), el('span', { style: { flex: '1', height: '10px', borderRadius: '3px', background: `linear-gradient(90deg, ${RAMP.join(',')})` } }), el('span', { text: `${fmt(max)} Å` }));
      root.append(lg);
      if (cur >= 0) {
        const k = order[cur];
        const near = order.filter(j => j !== k && ens.rmsd[k][j] < 2.0).length;
        const nearest = order.filter(j => j !== k).sort((a, b) => ens.rmsd[k][a] - ens.rmsd[k][b])[0];
        root.append(el('div', { class: 'desc', text: `${ens.shortId(k)}: ${near} other conformation(s) within 2 Å` + (nearest !== undefined ? `; nearest is ${ens.shortId(nearest)} at ${fmt(ens.rmsd[k][nearest])} Å.` : '.') }));
      }
      root.append(el('div', { class: 'desc', text: 'Rows and columns follow the conformation order of the left panel. RMSD after optimal superposition of C4′ atoms.' }));
    };
    api.ui.addPanel({ id: 'panel', title: 'Similarity', render(body) { root = body; render(); } });
    for (const ev of ['ensemble', 'conformation', 'order', 'tab']) api.on(ev, render);
  },
};
