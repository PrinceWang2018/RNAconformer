// Desktop-application behaviour: nothing should feel like a web page. Text is not selectable outside text fields,
// the browser's right-click menu, image dragging, page zoom, browser shortcuts (print, find, reload, view source,
// developer tools, back/forward) and middle-click auto-scroll are replaced by the application's own behaviour.
// Plugin developers can re-enable the browser tools in Preferences (Developer mode).
import { closeMenus, popupMenu, toast } from './ui.js';
import { chartMenuItems } from './export/chartexport.js';

const EDITABLE = 'input:not([type=checkbox]):not([type=radio]):not([type=range]):not([type=color]):not([type=file]), textarea, [contenteditable="true"]';

export function installNativeFeel(app) {
  const dev = () => !!app.prefs.devMode;
  document.body.classList.add('native');

  // ---- right mouse button: our own menus only
  document.addEventListener('contextmenu', e => {
    if (dev() && e.shiftKey) return;                       // developer mode: Shift+right-click = browser menu
    if (e.defaultPrevented) return;                         // a component showed its own menu
    e.preventDefault();
    closeMenus();
    const t = e.target.nodeType === 3 ? e.target.parentElement : e.target;
    const field = t.closest?.(EDITABLE);
    if (field) { editMenu(field, e); return; }
    const chart = t.closest?.('svg.chart');
    if (chart) { popupMenu({ x: e.clientX, y: e.clientY }, chartMenuItems(chart, chartName(chart), toast)); return; }
    const table = t.closest?.('table');
    if (table) { popupMenu({ x: e.clientX, y: e.clientY }, tableMenu(table)); return; }
    const sel = t.closest?.('.selectable');
    if (sel) popupMenu({ x: e.clientX, y: e.clientY }, [{ label: 'Copy text', action: () => copyText(sel.innerText) }]);
  });

  // ---- left mouse button: no text selection or dragging outside text fields
  document.addEventListener('selectstart', e => {
    const t = e.target.nodeType === 3 ? e.target.parentElement : e.target;
    if (!t?.closest?.(`${EDITABLE}, .selectable`)) e.preventDefault();
  });
  document.addEventListener('dragstart', e => {
    const t = e.target.nodeType === 3 ? e.target.parentElement : e.target;
    if (!t?.closest?.('[draggable="true"]')) e.preventDefault();
  });
  // middle-click auto-scroll and the mouse's back / forward buttons
  document.addEventListener('mousedown', e => { if (e.button === 1) e.preventDefault(); });
  const blockNav = e => { if (e.button === 3 || e.button === 4) { e.preventDefault(); e.stopPropagation(); } };
  document.addEventListener('mouseup', blockNav, true);
  document.addEventListener('auxclick', e => { if (e.button !== 0) e.preventDefault(); });
  try {
    history.pushState({ rcv: 1 }, '', location.href);
    window.addEventListener('popstate', () => history.pushState({ rcv: 1 }, '', location.href));
  } catch { /* history not available */ }

  // ---- page zoom (Ctrl+wheel, pinch) is not a page feature here
  window.addEventListener('wheel', e => { if (e.ctrlKey) e.preventDefault(); }, { passive: false, capture: true });

  // ---- browser shortcuts
  window.addEventListener('keydown', e => {
    const k = e.key.toLowerCase(), ctrl = e.ctrlKey || e.metaKey;
    if (dev() && (e.key === 'F12' || (ctrl && e.shiftKey && ['i', 'j', 'c'].includes(k)))) return;
    const block = () => { e.preventDefault(); e.stopPropagation(); };
    if (e.key === 'F5' || e.key === 'F12' || e.key === 'F7' || e.key === 'F3' || e.key === 'BrowserBack' || e.key === 'BrowserForward') return block();
    if (e.altKey && (e.key === 'ArrowLeft' || e.key === 'ArrowRight' || e.key === 'Home')) return block();
    if (!ctrl) return;
    const inField = !!e.target.closest?.(EDITABLE);
    if (k === 'f' && !e.shiftKey) { block(); app.console?.setVisible(true); if (app.console) app.console.input.value = 'select '; return; }   // find = select
    if (k === 'p' && !e.shiftKey) { block(); if (app.ens) app.action('export-image'); return; }                                                  // print = export image
    // keys the application uses itself (it handles them outside text fields): only suppress the browser default
    if (['n', 'o', 's', 'e', 'h'].includes(k)) { e.preventDefault(); return; }
    if (['a', 'c', 'v', 'x', 'z', 'y'].includes(k) && inField) return;                                                       // normal text editing
    if (['r', 'u', 'g', 'd', 'j', 'l', 'k', 't', 'm', 'b', 'i'].includes(k) && !e.shiftKey) return block();
    if (e.shiftKey && ['i', 'j', 'c', 'delete', 'r', 'n', 'b', 'o', 't'].includes(k)) return block();
    if (k === '+' || k === '=' || k === '-' || k === '0') {           // page zoom -> view zoom
      block();
      const c = app.viewer?.controls;
      if (!c || app.page !== 'viewer') return;
      if (k === '0') app.viewer.resetView(); else c.zoom = Math.max(0.05, Math.min(40, c.zoom * (k === '-' ? 1 / 1.15 : 1.15)));
      app.viewer.needsRender = true;
    }
  }, true);

  // ---- text fields: no spell-check squiggles in sequences and commands
  document.addEventListener('focusin', e => { const t = e.target; if (t.matches?.(EDITABLE) && !t.hasAttribute('spellcheck')) t.setAttribute('spellcheck', 'false'); });
}

function chartName(svgEl) {
  const sec = svgEl.closest('.section');
  const title = sec?.querySelector('.section-title')?.firstChild?.textContent || svgEl.getAttribute('aria-label') || 'chart';
  return String(title).trim().slice(0, 60) || 'chart';
}

async function copyText(text) {
  try { await navigator.clipboard.writeText(text); toast('Copied.', 'ok', '', 1500); }
  catch { toast('Copying is not available here.', 'error'); }
}

function tableMenu(table) {
  const rows = [...table.querySelectorAll('tr')].map(tr => [...tr.children].map(td => td.innerText.replace(/\s+/g, ' ').trim()).join('\t'));
  return [{ label: 'Copy table (tab-separated, for Excel)', action: () => copyText(rows.join('\n')) },
    { label: 'Copy as CSV', action: () => copyText(rows.map(r => r.split('\t').map(c => (/[",]/.test(c) ? `"${c.replace(/"/g, '""')}"` : c)).join(',')).join('\n')) }];
}

// Cut / Copy / Paste / Select all for text fields (the application's own menu, not the browser's).
function editMenu(field, e) {
  field.focus();
  const ro = field.readOnly || field.disabled;
  const sel = typeof field.selectionStart === 'number' ? field.value.slice(field.selectionStart, field.selectionEnd) : String(window.getSelection());
  const exec = cmd => () => { field.focus(); document.execCommand(cmd); };
  popupMenu({ x: e.clientX, y: e.clientY }, [
    { label: 'Undo', kbd: 'Ctrl+Z', action: exec('undo'), disabled: ro },
    { label: 'Redo', kbd: 'Ctrl+Y', action: exec('redo'), disabled: ro },
    { sep: true },
    { label: 'Cut', kbd: 'Ctrl+X', action: exec('cut'), disabled: ro || !sel },
    { label: 'Copy', kbd: 'Ctrl+C', action: exec('copy'), disabled: !sel },
    { label: 'Paste', kbd: 'Ctrl+V', disabled: ro, action: async () => {
      field.focus();
      try { const text = await navigator.clipboard.readText(); document.execCommand('insertText', false, text); }
      catch { toast('Press Ctrl+V to paste here.', 'info', '', 2500); }
    } },
    { sep: true },
    { label: 'Select all', kbd: 'Ctrl+A', action: () => { field.focus(); field.select ? field.select() : document.execCommand('selectAll'); } },
  ]);
}
