// Local compute jobs (prediction, pull exploration): polling, progress formatting and the viewer's jobs dock.
import { el } from './util.js';
import { icon } from './icons.js';

export function fmtDuration(s) {
  if (s === null || s === undefined || !Number.isFinite(s)) return '';
  s = Math.max(0, Math.round(s));
  if (s < 50) return `${Math.max(5, Math.ceil(s / 5) * 5)} s`;
  if (s < 3600) { const m = Math.round(s / 60); return `${m} min`; }
  const h = Math.floor(s / 3600), m = Math.round((s % 3600) / 60);
  return `${h} h ${m} min`;
}

// Poll a job until it finishes. Callbacks: onUpdate(info), onDone(info), onFail(info).
export function watchJob(api, id, { onUpdate, onDone, onFail, interval = 1200 } = {}) {
  let stopped = false, timer = 0;
  const tick = async () => {
    if (stopped) return;
    try {
      const info = await api.get(`/api/jobs/${id}`);
      onUpdate?.(info);
      if (info.status === 'done') { onDone?.(info); return; }
      if (info.status === 'failed' || info.status === 'cancelled') { onFail?.(info); return; }
    } catch (err) {
      onUpdate?.({ id, status: 'running', error: err.message });
    }
    timer = setTimeout(tick, interval);
  };
  tick();
  return () => { stopped = true; clearTimeout(timer); };
}

export function stageText(info) {
  const active = (info.stages || []).find(s => s.state === 'active') || (info.stages || []).filter(s => s.state === 'done').pop();
  return active ? `${active.label}${active.detail ? ' · ' + active.detail : ''}` : info.status;
}

// Compact list of running jobs in the 3D view.
export class JobsDock {
  constructor(app, host) {
    this.app = app;
    this.host = host;
    this.items = new Map();
  }

  track(info, { title, onDone } = {}) {
    const node = el('div', { class: 'job-item' });
    const bar = el('div', { class: 'progress' }, el('i'));
    const text = el('div', { class: 'job-text' });
    const cancel = el('button', { class: 'btn ghost icon-only small', title: 'Cancel', 'aria-label': 'Cancel job',
      onclick: () => this.app.api.post(`/api/jobs/${info.id}/cancel`, {}) }, icon('x'));
    node.append(el('div', { class: 'job-head' }, el('span', { class: 'spinner small' }), el('b', { text: title || info.title }), el('span', { class: 'grow' }), cancel), bar, text);
    this.host.append(node);
    this.host.hidden = false;
    const finish = () => { setTimeout(() => { node.remove(); this.items.delete(info.id); this.host.hidden = !this.items.size; }, 2500); };
    const stop = watchJob(this.app.api, info.id, {
      onUpdate: j => {
        bar.firstChild.style.width = `${Math.round((j.progress || 0) * 100)}%`;
        text.textContent = `${stageText(j)}${j.eta_s ? ' · about ' + fmtDuration(j.eta_s) + ' left' : ''}`;
      },
      onDone: async j => { node.classList.add('ok'); text.textContent = 'Finished'; finish(); await onDone?.(j); },
      onFail: j => { node.classList.add('bad'); text.textContent = j.status === 'cancelled' ? 'Cancelled' : (j.error || 'Failed'); finish(); },
    });
    this.items.set(info.id, { node, stop });
  }

  get running() { return this.items.size; }
}
