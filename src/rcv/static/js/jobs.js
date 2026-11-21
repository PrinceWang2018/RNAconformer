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
// A failed request changes nothing on screen (the last state stays). After `lostAfter` failures in a row the state is
// passed on with `lost: true`, so the page can say that the connection to the engine is gone, and polling goes on
// at a slower pace: if RCV answers again, the job is followed as before. A job the engine no longer knows (RCV was
// restarted) ends as failed.
export function watchJob(api, id, { onUpdate, onDone, onFail, interval = 1200, lostAfter = 5 } = {}) {
  let stopped = false, timer = 0, fails = 0, last = { id, status: 'running' };
  const tick = async () => {
    if (stopped) return;
    let wait = interval;
    try {
      const info = await api.get(`/api/jobs/${id}`);
      if (stopped) return;
      fails = 0;
      last = info;
      onUpdate?.(info);
      if (info.status === 'done') { onDone?.(info); return; }
      if (info.status === 'failed' || info.status === 'cancelled') { onFail?.(info); return; }
    } catch (err) {
      if (stopped) return;
      if (err.status === 404) {
        onFail?.({ ...last, status: 'failed', error: 'RCV was restarted, so this run is no longer followed. Results that were already saved are in the projects folder.' });
        return;
      }
      fails += 1;
      if (fails >= lostAfter) { onUpdate?.({ ...last, lost: true }); wait = Math.max(interval, 4000); }
    }
    timer = setTimeout(tick, wait);
  };
  tick();
  return () => { stopped = true; clearTimeout(timer); };
}

// Cancel a job. A run of more than a minute is only cancelled after a question (one click would discard it).
// Resolves true when the request to cancel was sent.
export async function cancelJob(app, info, exploration = false) {
  if ((info.elapsed_s || 0) > 60) {
    const { askConfirm } = await import('./ui.js');
    const ok = await askConfirm(exploration ? 'Cancel this exploration?' : 'Cancel this run?',
      `It has been running for ${fmtDuration(info.elapsed_s)}. What has been computed so far is discarded and cannot be continued later.`,
      exploration ? 'Cancel the exploration' : 'Cancel the run', 'Keep running');
    if (!ok) return false;
  }
  app.api.post(`/api/jobs/${info.id}/cancel`, {}).catch(() => { /* the watcher reports what happens next */ });
  return true;
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
    let latest = info;
    const cancel = el('button', { class: 'btn ghost icon-only small', title: 'Cancel', 'aria-label': 'Cancel job',
      onclick: async () => { if (await cancelJob(this.app, latest, true)) { cancel.disabled = true; text.textContent = 'Cancelling…'; } } }, icon('x'));
    node.append(el('div', { class: 'job-head' }, el('span', { class: 'spinner small' }), el('b', { text: title || info.title }), el('span', { class: 'grow' }), cancel), bar, text);
    this.host.append(node);
    this.host.hidden = false;
    const finish = () => { setTimeout(() => { node.remove(); this.items.delete(info.id); this.host.hidden = !this.items.size; }, 2500); };
    const stop = watchJob(this.app.api, info.id, {
      onUpdate: j => {
        latest = j;
        bar.firstChild.style.width = `${Math.round((j.progress || 0) * 100)}%`;
        text.textContent = j.lost ? 'Connection to the engine lost — trying to reconnect…' : j.cancelling ? 'Cancelling…'
          : `${stageText(j)}${j.eta_s ? ' · about ' + fmtDuration(j.eta_s) + ' left' : ''}`;
      },
      onDone: async j => { node.classList.add('ok'); text.textContent = 'Finished'; finish(); await onDone?.(j); },
      onFail: j => { node.classList.add('bad'); text.textContent = j.status === 'cancelled' ? 'Cancelled' : (j.error || 'Failed'); finish(); },
    });
    this.items.set(info.id, { node, stop });
  }

  get running() { return this.items.size; }
}
