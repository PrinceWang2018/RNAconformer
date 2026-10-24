// Left-hand ensemble/conformation lists and the bottom timeline scrubber.
import { el, svg, fmt, fmtVal, clamp } from '../util.js';
import { icon } from '../icons.js';
import { linear, extent } from './charts.js';
import { tooltip, tipRow, popupMenu } from '../ui.js';

export class EnsembleList {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    for (const ev of ['documents', 'ensemble']) app.on(ev, () => this.render());
  }
  render() {
    const app = this.app;
    this.root.replaceChildren();
    if (!app.ensembles.size) {
      this.root.append(el('div', { class: 'muted', style: { padding: '6px 10px 10px', fontSize: '12px' }, text: 'No ensembles yet. Open PDB/mmCIF files, an ensemble folder, or start a prediction.' }));
      return;
    }
    for (const ens of app.ensembles.values()) {
      const engine = ens.doc.origin === 'rnafold3d-engine';
      const item = el('div', { class: `ens-item ${app.ens === ens ? 'on' : ''}`, role: 'button', tabindex: 0, title: ens.doc.source?.path || ens.name },
        el('div', { class: 'meta' },
          el('div', { class: 'name', text: ens.name }),
          el('div', { class: 'sub' }, `${ens.n} conformation${ens.n === 1 ? '' : 's'} · ${ens.nRes} nt `, el('span', { class: `badge ${engine ? 'engine' : ''}`, text: engine ? 'predicted' : 'structures' }))),
        el('button', { class: 'btn ghost icon-only small x', title: 'Close ensemble', 'aria-label': 'Close ensemble', onclick: e => { e.stopPropagation(); app.closeEnsemble(ens.id); } }, icon('x')));
      item.addEventListener('click', () => app.selectEnsemble(ens.id));
      item.addEventListener('keydown', e => { if (e.key === 'Enter') app.selectEnsemble(ens.id); });
      this.root.append(item);
    }
  }
}

export class ConformationList {
  constructor(app, root) {
    this.app = app;
    this.root = root;
    for (const ev of ['ensemble', 'conformation', 'metric', 'order', 'values']) app.on(ev, () => this.render());
  }
  render() {
    const app = this.app, ens = app.ens;
    this.root.replaceChildren();
    if (!ens) return;
    const key = app.state.metric, m = ens.metricByKey[key] || {};
    const vals = ens.conformations.map((_, k) => ens.value(k, key));
    const [lo, hi] = extent(vals, 0);
    app.order.forEach((k, i) => {
      const c = ens.conformations[k];
      const v = vals[k];
      const valid = c.validation ? c.validation.passed : null;
      const t = v === null ? 0 : (v - lo) / ((hi - lo) || 1);
      const item = el('div', { class: `conf-item ${k === app.state.conf ? 'on' : ''} ${k === app.state.prev && k !== app.state.conf ? 'prev' : ''}`, role: 'option', 'aria-selected': k === app.state.conf ? 'true' : 'false', tabindex: -1 },
        el('div', { class: 'num', text: String(i + 1) }),
        el('div', { class: 'title', text: ens.shortId(k), title: c.label }),
        el('div', { class: 'val', text: fmtVal(v, m.unit), title: `${m.label || key}${m.unit ? ' (' + m.unit + ')' : ''}` }),
        el('div', { class: 'bar' }, el('i', { style: { width: `${Math.max(2, t * 100)}%` } })),
        el('div', { class: 'sub' },
          valid === null ? null : el('span', { class: valid ? 'state-ok' : 'state-bad', title: valid ? 'passes all-atom validation' : 'validation issues' }, icon(valid ? 'check' : 'warn')),
          c.explored ? el('span', { class: 'badge explored', title: 'Found by exploring a pull with the engine', text: `pull ${c.explored.pull}` }) : null,
          c.cluster && !c.explored ? el('span', { text: `cl ${c.cluster}` }) : null,
          el('span', { class: 'dbn', text: c.dbn, title: c.dbn })));
      item.addEventListener('click', () => app.goTo(k));
      item.addEventListener('contextmenu', e => {
        e.preventDefault();
        const cur = app.state.conf;
        popupMenu({ x: e.clientX, y: e.clientY }, [
          { label: `Go to ${ens.shortId(k)}`, action: () => app.goTo(k) },
          { label: k === cur ? 'Transition video…' : `Transition video ${ens.shortId(cur)} → ${ens.shortId(k)}…`, action: () => app.exportTransition(cur === k ? (app.state.prev ?? app.order[0]) : cur, k), disabled: ens.n < 2 },
          { label: `Compare with ${ens.shortId(cur)} (switch analysis)`, action: () => app.compareWith?.(k), disabled: k === cur },
          { sep: true },
          { label: 'Save this conformation as PDB…', action: () => app.exportPDB(false, k) },
        ]);
      });
      item.addEventListener('pointerenter', () => { if (k !== app.state.conf) app.previewConformation(k); });
      item.addEventListener('pointerleave', () => app.previewConformation(null));
      this.root.append(item);
    });
    this.root.querySelector('.conf-item.on')?.scrollIntoView({ block: 'nearest' });
  }
}

export class Timeline {
  constructor(app, track, posEl) {
    this.app = app;
    this.track = track;
    this.posEl = posEl;
    this.svg = svg('svg');
    track.append(this.svg);
    for (const ev of ['ensemble', 'conformation', 'metric', 'order', 'values']) app.on(ev, () => this.render());
    app.on('morph', m => this.updateHandle(m));
    new ResizeObserver(() => this.render()).observe(track);
    let dragging = false;
    const posAt = e => {
      const r = track.getBoundingClientRect();
      const n = this.app.order.length;
      return clamp(this.x.invert(e.clientX - r.left), 0, Math.max(0, n - 1));
    };
    track.addEventListener('pointerdown', e => {
      if (!this.app.ens) return;
      dragging = true;
      track.setPointerCapture(e.pointerId);
      this.app.scrub(posAt(e));
    });
    track.addEventListener('pointermove', e => {
      if (!this.app.ens || !this.x) return;
      if (dragging) { this.app.scrub(posAt(e)); return; }
      const p = Math.round(posAt(e)), k = this.app.order[p];
      const ens = this.app.ens, key = this.app.state.metric;
      tooltip.show([tipRow(fmtVal(ens.value(k, key), ens.metricByKey[key]?.unit), `${ens.shortId(k)} · position ${p + 1}`)], e.clientX, e.clientY);
    });
    track.addEventListener('pointerleave', () => tooltip.hide());
    const end = e => { if (!dragging) return; dragging = false; this.app.endScrub(posAt(e)); };
    track.addEventListener('pointerup', end);
    track.addEventListener('pointercancel', end);
    track.addEventListener('keydown', e => {
      if (e.key === 'ArrowLeft') { this.app.prev(); e.preventDefault(); }
      if (e.key === 'ArrowRight') { this.app.next(); e.preventDefault(); }
    });
  }

  render() {
    const app = this.app, ens = app.ens, s = this.svg;
    s.replaceChildren();
    const W = this.track.clientWidth || 400, H = this.track.clientHeight || 52;
    s.setAttribute('viewBox', `0 0 ${W} ${H}`);
    if (!ens) { this.posEl.textContent = '–'; this.x = null; return; }
    const n = app.order.length;
    this.x = linear(0, Math.max(1, n - 1), 10, W - 10);
    const key = app.state.metric;
    const vals = app.order.map(k => ens.value(k, key));
    const [lo, hi] = extent(vals, 0.1);
    const y = linear(lo, hi, H - 12, 6);
    const base = H - 12;
    s.append(svg('line', { x1: 10, x2: W - 10, y1: base + 0.5, y2: base + 0.5, stroke: 'var(--axis)' }));
    const pts = vals.map((v, i) => (v === null ? null : [this.x(i), y(v)])).filter(Boolean);
    if (pts.length > 1) {
      s.append(svg('path', { d: `M ${pts[0][0]} ${base} ` + pts.map(p => `L ${p[0]} ${p[1]}`).join(' ') + ` L ${pts[pts.length - 1][0]} ${base} Z`, fill: 'var(--series-1)', opacity: 0.1 }));
      s.append(svg('polyline', { points: pts.map(p => p.join(',')).join(' '), fill: 'none', stroke: 'var(--series-1)', 'stroke-width': 1.5, opacity: 0.8 }));
    }
    for (let i = 0; i < n; i++) s.append(svg('circle', { cx: this.x(i), cy: base, r: app.order[i] === app.state.conf ? 0 : 3, fill: 'var(--axis)' }));
    this.handle = svg('g', { 'pointer-events': 'none' },
      svg('line', { x1: 0, x2: 0, y1: 2, y2: base, stroke: 'var(--text)', 'stroke-width': 1.5 }),
      svg('circle', { cx: 0, cy: base, r: 6.5, fill: 'var(--accent)', stroke: 'var(--surface)', 'stroke-width': 2 }));
    s.append(this.handle);
    this.updateHandle(null);
  }

  updateHandle(m) {
    const app = this.app;
    if (!this.x || !app.ens) return;
    const order = app.order;
    let p;
    if (m && !m.done) {
      const a = order.indexOf(m.from), b = order.indexOf(m.to);
      p = a + (b - a) * m.t;
    } else p = order.indexOf(app.state.conf);
    const X = this.x(p);
    this.handle.setAttribute('transform', `translate(${X},0)`);
    const near = Math.round(p);
    this.posEl.textContent = m && !m.done && Math.abs(p - near) > 0.02 ? `${(p + 1).toFixed(1)} / ${order.length}` : `${near + 1} / ${order.length}`;
    this.track.setAttribute('aria-valuenow', String(near + 1));
    this.track.setAttribute('aria-valuemax', String(order.length));
    this.track.setAttribute('aria-valuemin', '1');
  }
}
