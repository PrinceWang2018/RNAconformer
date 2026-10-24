// Two-sided arc diagram: base pairs of the current conformation (above the sequence axis) and of a
// comparison conformation (below). Pairs gained/lost between the two are emphasised.
import { svg, el } from '../util.js';
import { tooltip, tipRow } from '../ui.js';

export class ArcDiagram {
  constructor(app, host) {
    this.app = app;
    this.host = host;
    this.svg = svg('svg', { class: 'chart', role: 'img' });
    this.legend = el('div', { class: 'chips', style: { marginTop: '6px', fontSize: '11.5px', color: 'var(--text-2)', gap: '12px' } });
    host.append(this.svg, this.legend);
  }

  render(width) {
    const app = this.app, ens = app.ens;
    const s = this.svg;
    s.replaceChildren();
    this.legend.replaceChildren();
    if (!ens) return;
    const k = app.state.conf, c = app.state.prev ?? k;
    const n = ens.nRes;
    const W = Math.max(200, width), H = 150, base = 74, pad = 8;
    s.setAttribute('viewBox', `0 0 ${W} ${H}`);
    s.setAttribute('height', H);
    s.setAttribute('aria-label', `Arc diagram of base pairs in ${ens.shortId(k)} (above) and ${ens.shortId(c)} (below)`);
    const x = i => pad + (n > 1 ? i * (W - 2 * pad) / (n - 1) : (W - 2 * pad) / 2);
    const key = p => `${p[0]}-${p[1]}`;
    const cur = ens.conformations[k].pairs, cmp = ens.conformations[c].pairs;
    const curSet = new Set(cur.map(key)), cmpSet = new Set(cmp.map(key));
    const maxH = base - 10;
    const span = W - 2 * pad;
    const arc = (i, j, up) => {
      const x1 = x(i), x2 = x(j);
      const h = Math.min(maxH, Math.max(6, (x2 - x1) / span * maxH * 1.9 + 4));
      const cy = up ? base - h * 4 / 3 : base + h * 4 / 3;
      return `M ${x1} ${base} C ${x1} ${cy} ${x2} ${cy} ${x2} ${base}`;
    };
    const g = svg('g');
    // sequence axis
    g.append(svg('line', { x1: pad, x2: W - pad, y1: base, y2: base, stroke: 'var(--axis)', 'stroke-width': 1 }));
    const step = n <= 40 ? 1 : n <= 80 ? 5 : 10;
    for (let i = 0; i < n; i++) {
      if (n <= 40) g.append(svg('text', { x: x(i), y: base + 3.5, 'text-anchor': 'middle', 'font-size': 8.5, fill: 'var(--text-2)', stroke: 'var(--surface)', 'stroke-width': 2.5, 'paint-order': 'stroke', text: ens.residues[i].name }));
      else if ((i + 1) % step === 0) g.append(svg('line', { x1: x(i), x2: x(i), y1: base - 2, y2: base + 2, stroke: 'var(--axis)' }));
    }
    g.append(svg('text', { x: pad, y: 10, 'font-size': 10, fill: 'var(--text-3)', text: `${ens.shortId(k)} (current)` }));
    g.append(svg('text', { x: pad, y: H - 3, 'font-size': 10, fill: 'var(--text-3)', text: c === k ? 'same conformation' : `${ens.shortId(c)} (previous)` }));
    const hover = (p, label) => {
      let content = null;
      const build = () => {
        const t = p[2] ? `${p[2][0]}–${p[2][1]}` : 'pair';
        return [tipRow(`${ens.resLabel(p[0])}–${ens.resLabel(p[1])}`, ''),
          el('div', { class: 'tl', text: `${t} · ${label}${p[3] ? ' · pseudoknot' : ''}` })];
      };
      return {
        onpointerenter: e => { app.hoverResidues([p[0], p[1]]); content = build(); tooltip.show(content, e.clientX, e.clientY); },
        onpointermove: e => { if (content) tooltip.show(content, e.clientX, e.clientY); },
        onpointerleave: () => { content = null; app.hoverResidues([]); tooltip.hide(); },
      };
    };
    let gained = 0, lost = 0;
    for (const p of cur) {
      const isNew = c !== k && !cmpSet.has(key(p));
      if (isNew) gained++;
      const d = arc(p[0], p[1], true);
      g.append(svg('path', { d, fill: 'none', stroke: isNew ? 'var(--series-2)' : (p[3] ? 'var(--series-2)' : 'var(--series-1)'), 'stroke-width': isNew ? 2.2 : 1.4, opacity: isNew ? 1 : 0.85, 'stroke-dasharray': p[3] && !isNew ? '3 2' : null }));
      g.append(svg('path', { d, fill: 'none', stroke: 'transparent', 'stroke-width': 8, class: 'hit', ...hover(p, isNew ? 'gained' : 'shared') }));
    }
    if (c !== k) {
      for (const p of cmp) {
        const isLost = !curSet.has(key(p));
        if (isLost) lost++;
        const d = arc(p[0], p[1], false);
        g.append(svg('path', { d, fill: 'none', stroke: isLost ? 'var(--series-2)' : 'var(--axis)', 'stroke-width': isLost ? 2.2 : 1.2, 'stroke-dasharray': isLost ? '4 3' : null }));
        g.append(svg('path', { d, fill: 'none', stroke: 'transparent', 'stroke-width': 8, class: 'hit', ...hover(p, isLost ? 'lost in current' : 'shared') }));
      }
    }
    s.append(g);
    const key2 = (color, label, dash) => el('span', { class: 'lg', style: { display: 'inline-flex', alignItems: 'center', gap: '6px' } },
      el('span', { style: { display: 'inline-block', width: '16px', height: '0', borderTop: `2px ${dash ? 'dashed' : 'solid'} ${color}` } }), label);
    this.legend.append(key2('var(--series-1)', `${cur.length} pairs`), key2('var(--series-2)', `${gained} gained`), key2('var(--series-2)', `${lost} lost`, true));
  }
}
