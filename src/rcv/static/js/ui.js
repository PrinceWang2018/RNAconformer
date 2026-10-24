// UI primitives: tooltip, toasts, menus, dialogs, busy overlay.
import { el, $ } from './util.js';
import { icon } from './icons.js';

// ------------------------------------------------------------------ tooltip
const tip = () => $('#tooltip');
export const tooltip = {
  show(content, x, y) {
    const t = tip();
    t.replaceChildren(...(Array.isArray(content) ? content : [content]));
    t.hidden = false;
    const r = t.getBoundingClientRect();
    let left = x + 14, top = y + 14;
    if (left + r.width > window.innerWidth - 8) left = x - r.width - 14;
    if (top + r.height > window.innerHeight - 8) top = y - r.height - 14;
    t.style.left = `${Math.max(8, left)}px`;
    t.style.top = `${Math.max(8, top)}px`;
  },
  hide() { tip().hidden = true; },
};
export function tipRow(value, label, color) {
  return el('div', {}, color ? el('span', { class: 'key', style: { background: color } }) : null,
    el('span', { class: 'tv', text: value }), ' ', el('span', { class: 'tl', text: label || '' }));
}

// ------------------------------------------------------------------ toasts
export function toast(message, kind = 'info', title = '', ms = 4500) {
  const box = $('#toasts');
  const node = el('div', { class: `toast ${kind}`, role: kind === 'error' ? 'alert' : 'status' },
    el('span', { class: 'ti' }, icon(kind === 'error' ? 'warn' : kind === 'ok' ? 'check' : 'info')),
    el('div', { class: 'tb' }, title ? el('b', { text: title }) : null, el('span', { text: message })));
  box.append(node);
  setTimeout(() => node.remove(), kind === 'error' ? Math.max(ms, 8000) : ms);
  return node;
}

// ------------------------------------------------------------------ menus
let openMenu = null;
export function closeMenus() {
  if (openMenu) { openMenu.pop.remove(); openMenu.btn.classList?.remove('open'); openMenu = null; }
}
document.addEventListener('pointerdown', e => {
  if (openMenu && !openMenu.pop.contains(e.target) && !(openMenu.btn.contains && openMenu.btn.contains(e.target))) closeMenus();
});
document.addEventListener('keydown', e => { if (e.key === 'Escape') closeMenus(); });

// items: [{label, action, kbd, checked, disabled, sep, head}] or a function returning such a list
export function popupMenu(anchorBtn, items) {
  closeMenus();
  const list = typeof items === 'function' ? items() : items;
  const pop = el('div', { class: 'menu-pop', role: 'menu' });
  for (const it of list) {
    if (it.sep) { pop.append(el('div', { class: 'menu-sep' })); continue; }
    if (it.head) { pop.append(el('div', { class: 'menu-head', text: it.head })); continue; }
    const b = el('button', { class: 'menu-item', role: 'menuitem', disabled: it.disabled ? true : null },
      el('span', { class: 'check' }, it.checked ? icon('check') : null),
      el('span', { text: it.label }), it.kbd ? el('span', { class: 'kbd', text: it.kbd }) : null);
    b.addEventListener('click', () => { closeMenus(); it.action?.(); });
    pop.append(b);
  }
  document.body.append(pop);
  // anchor: a button (menu opens below it) or a point { x, y } (context menus)
  const r = anchorBtn.getBoundingClientRect ? anchorBtn.getBoundingClientRect() : { left: anchorBtn.x, bottom: anchorBtn.y - 4 };
  const pr = pop.getBoundingClientRect();
  pop.style.left = `${Math.max(4, Math.min(r.left, window.innerWidth - pr.width - 8))}px`;
  pop.style.top = `${Math.max(4, Math.min(r.bottom + 4, window.innerHeight - pr.height - 8))}px`;
  anchorBtn.classList?.add('open');
  openMenu = { pop, btn: anchorBtn };
  pop.querySelector('.menu-item:not(:disabled)')?.focus();
}

export function menuBar(root, menus) {
  root.replaceChildren();
  for (const m of menus) {
    const b = el('button', { class: 'menu-btn', text: m.label, 'aria-haspopup': 'menu' });
    b.addEventListener('click', () => { if (openMenu?.btn === b) closeMenus(); else popupMenu(b, m.items); });
    b.addEventListener('pointerenter', () => { if (openMenu && openMenu.btn !== b) popupMenu(b, m.items); });
    root.append(b);
  }
}

// ------------------------------------------------------------------ dialogs
export function dialog({ title, width = 560, body, actions = [], onClose }) {
  const dlg = el('dialog', { class: 'dlg', style: { '--dlg-w': `${width}px` } });
  const close = () => { dlg.close(); };
  const head = el('div', { class: 'dlg-head' }, el('h2', { text: title }),
    el('button', { class: 'btn ghost icon-only', type: 'button', 'aria-label': 'Close', onclick: close }, icon('x')));
  const bodyEl = el('div', { class: 'dlg-body' });
  const foot = el('div', { class: 'dlg-foot' });
  const btns = {};
  for (const a of actions) {
    if (a.spacer) { foot.append(el('span', { class: 'grow' })); continue; }
    if (a.node) { foot.append(a.node); continue; }
    const b = el('button', { class: `btn ${a.primary ? 'primary' : ''}`, type: 'button', text: a.label });
    b.addEventListener('click', async () => {
      if (a.close) { close(); return; }
      try { await a.action?.(b); } catch (err) { toast(err.message || String(err), 'error', 'Error'); }
    });
    btns[a.id || a.label] = b;
    foot.append(b);
  }
  dlg.append(el('div', { class: 'dlg-form' }, head, bodyEl, actions.length ? foot : null));
  if (typeof body === 'function') body(bodyEl); else if (body) bodyEl.append(body);
  dlg.addEventListener('close', () => { onClose?.(); dlg.remove(); });
  document.body.append(dlg);
  dlg.showModal();
  return { dlg, body: bodyEl, close, buttons: btns };
}

// Small text prompt (window.prompt is not available in every desktop web view).
export function askText(title, label, value = '', { pattern = null, hint = '' } = {}) {
  return new Promise(resolve => {
    const inp = el('input', { class: 'input', value, 'aria-label': label, spellcheck: 'false' });
    const msg = el('div', { class: 'help', text: hint });
    let result = null;
    const ok = () => {
      const v = inp.value.trim();
      if (!v) return;
      if (pattern && !pattern.test(v)) { msg.textContent = hint || 'invalid'; msg.classList.add('state-bad'); return; }
      result = v; d.close();
    };
    const d = dialog({ title, width: 420, body: b => b.append(el('div', { class: 'field' }, el('label', { text: label }), inp, msg)),
      actions: [{ label: 'Cancel', close: true }, { label: 'OK', primary: true, action: ok }], onClose: () => resolve(result) });
    inp.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); ok(); } });
    setTimeout(() => { inp.focus(); inp.select(); }, 0);
  });
}

export function field(label, control, help) {
  return el('div', { class: 'field' }, el('label', { text: label }), control, help ? el('div', { class: 'help', text: help }) : null);
}
export function select(options, value, attrs = {}) {
  const s = el('select', { class: 'select', ...attrs });
  for (const o of options) {
    if (o.group) {
      const g = el('optgroup', { label: o.group });
      for (const x of o.options) g.append(el('option', { value: x.value, text: x.label, selected: x.value === value ? true : null }));
      s.append(g);
    } else s.append(el('option', { value: o.value, text: o.label, selected: o.value === value ? true : null }));
  }
  return s;
}

// ------------------------------------------------------------------ busy overlay
let busyCount = 0;
export function busy(text) {
  busyCount++;
  $('#busy-text').textContent = text || 'Working…';
  $('#busy').hidden = false;
  return () => { busyCount = Math.max(0, busyCount - 1); if (!busyCount) $('#busy').hidden = true; };
}
