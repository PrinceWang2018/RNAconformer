// Built-in plugin: track residue–residue distances across conformations.
// Each tracked distance becomes a metric, so the energy panel can plot it (profile, funnel, map),
// and a dashed line with the live distance follows the structure while it morphs.
const ATOMS = [
  { id: "C4'", label: 'C4′' }, { id: 'P', label: 'P' }, { id: "C1'", label: 'C1′' },
  { id: 'glyco', label: 'N1/N9 (base)' }, { id: 'centre', label: 'Residue centre' },
];
const COLORS = ['#eb6834', '#6250d6', '#1baf7a', '#e34948', '#2a78d6'];

export default {
  id: 'distance-tracker',
  name: 'Distance tracker',
  version: '1.0',
  description: 'Measure distances between residues and follow them across conformations.',
  activate(api) {
    const { el, svg, fmt } = api.util;
    const tracks = [];
    let counter = 0;
    let root = null;

    const point = (ens, coords, r, atom) => {
      if (atom === 'centre') return ens.residueCenter(coords, r);
      const idx = atom === 'glyco' ? ens.glyco[r] : ens.atomIndex[r][atom];
      if (idx === undefined || idx < 0) return ens.residueCenter(coords, r);
      return [coords[3 * idx], coords[3 * idx + 1], coords[3 * idx + 2]];
    };
    const dist = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
    const label = (ens, t) => `${ens.resLabel(t.a)}–${ens.resLabel(t.b)}`;

    api.viewer.addOverlay('lines', ({ coords, ens }) => tracks.filter(t => t.ens === ens.id).map(t => {
      const a = point(ens, coords, t.a, t.atom), b = point(ens, coords, t.b, t.atom);
      return { type: 'line', a, b, color: t.color, dash: true, width: 2, label: `${dist(a, b).toFixed(1)} Å` };
    }));

    const add = () => {
      const ens = api.ensemble;
      const sel = api.selection;
      if (!ens || sel.length < 2) { api.ui.toast('Select two residues first (click one, Shift+click another).', 'info'); return; }
      const [a, b] = sel.slice(-2).sort((x, y) => x - y);
      const atom = root.querySelector('select[name=atom]').value;
      const t = { id: ++counter, ens: ens.id, a, b, atom, color: COLORS[(counter - 1) % COLORS.length] };
      t.key = `d${t.id}`;
      tracks.push(t);
      const atomLabel = ATOMS.find(x => x.id === atom).label;
      api.metrics.register({
        key: t.key, label: `Distance ${label(ens, t)} (${atomLabel})`, unit: 'Å', lower_is_better: null,
        description: `Distance between ${atomLabel} of ${label(ens, t)} in each conformation (Distance tracker plugin).`,
        compute: (e, k) => (e.id === t.ens ? dist(point(e, e.coords[k], t.a, t.atom), point(e, e.coords[k], t.b, t.atom)) : null),
      });
      api.viewer.requestRender();
      render();
    };

    const render = () => {
      if (!root) return;
      const ens = api.ensemble;
      const list = root.querySelector('.dt-list');
      list.replaceChildren();
      const mine = tracks.filter(t => ens && t.ens === ens.id);
      if (!mine.length) { list.append(el('div', { class: 'desc', text: 'No distances yet.' })); return; }
      for (const t of mine) {
        const vals = api.order.map(k => dist(point(ens, ens.coords[k], t.a, t.atom), point(ens, ens.coords[k], t.b, t.atom)));
        const cur = dist(point(ens, api.coords, t.a, t.atom), point(ens, api.coords, t.b, t.atom));
        const lo = Math.min(...vals), hi = Math.max(...vals);
        const W = 120, H = 28;
        const spark = svg('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H, class: 'chart' });
        const X = i => 3 + (vals.length > 1 ? i / (vals.length - 1) : 0.5) * (W - 6);
        const Y = v => H - 4 - ((v - lo) / ((hi - lo) || 1)) * (H - 8);
        spark.append(svg('polyline', { points: vals.map((v, i) => `${X(i)},${Y(v)}`).join(' '), fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 1.5 }));
        const ci = api.order.indexOf(api.conformation);
        if (ci >= 0) spark.append(svg('circle', { cx: X(ci), cy: Y(vals[ci]), r: 3.5, fill: t.color, stroke: 'var(--surface)', 'stroke-width': 1.5 }));
        const row = el('div', { style: { display: 'grid', gridTemplateColumns: '12px 1fr auto', gap: '4px 8px', alignItems: 'center', padding: '8px 0', borderBottom: '1px solid var(--grid)' } },
          el('span', { style: { width: '10px', height: '10px', borderRadius: '50%', background: t.color, display: 'inline-block' } }),
          el('div', {}, el('div', { style: { fontWeight: 600 }, text: `${label(ens, t)} · ${ATOMS.find(x => x.id === t.atom).label}` }),
            el('div', { class: 'muted', style: { fontSize: '11.5px' }, text: `now ${fmt(cur)} Å · range ${fmt(lo)}–${fmt(hi)} Å` })),
          el('div', { style: { display: 'flex', gap: '4px' } },
            el('button', { class: 'btn small', text: 'Plot', title: 'Show this distance in the energy panel charts', onclick: () => api.setMetric(`plugin.distance-tracker.${t.key}`) }),
            el('button', { class: 'btn small ghost', text: '✕', title: 'Remove', 'aria-label': 'Remove distance', onclick: () => { tracks.splice(tracks.indexOf(t), 1); api.viewer.requestRender(); render(); } })),
          el('span'), spark);
        list.append(row);
      }
    };

    api.ui.addPanel({
      id: 'panel', title: 'Distances',
      render(body) {
        root = body;
        const atomSel = el('select', { class: 'select small', name: 'atom', 'aria-label': 'Atom used for the distance' }, ...ATOMS.map(a => el('option', { value: a.id, text: a.label })));
        body.append(
          el('div', { class: 'section-title', text: 'Distance tracker' }),
          el('div', { class: 'desc', style: { marginBottom: '10px' }, text: 'Select two residues in the 3D view (click, then Shift+click), choose the atom, and add. The distance is followed through every conformation and can be plotted in the energy panel.' }),
          el('div', { style: { display: 'flex', gap: '6px', marginBottom: '8px' } }, atomSel, el('button', { class: 'btn small primary', text: 'Add from selection', onclick: add })),
          el('div', { class: 'dt-list' }));
        render();
      },
    });
    for (const ev of ['ensemble', 'conformation', 'order']) api.on(ev, render);
    let raf = 0;
    api.on('morph', () => { if (!raf) raf = requestAnimationFrame(() => { raf = 0; render(); }); });
  },
};
