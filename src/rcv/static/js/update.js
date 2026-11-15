// Check for updates: Help > Check for updates…, and once a day when an installed RCV starts (can be turned off).
// RCV never downloads or installs anything by itself: it says which version is available and, when the user
// chooses Download, opens the installer for this computer (Windows, macOS Apple Silicon / Intel, Linux .deb / .rpm)
// in the web browser. The local server (rcv/update.py) does the network request and picks the installer.
import { el } from './util.js';
import { dialog, toast } from './ui.js';
import { openExternal } from './dialogs.js';

const TITLE = 'Check for updates';

// Help > Check for updates…
export async function checkForUpdates(app) {
  toast('Checking for updates…', 'info', TITLE, 2500);
  let r;
  try { r = await app.api.get('/api/update'); }
  catch (err) { toast(err.message || String(err), 'error', TITLE, 8000); return; }
  if (r.newer) updateDialog(app, r);
  else toast(`You have the newest version of RCV (${r.current}).`, 'ok', TITLE);
}

// At start-up: quiet when there is nothing new, the check is off or was made today, or the computer is offline.
export async function autoCheckForUpdates(app) {
  if (!app.api.online || !app.info?.update?.auto) return;
  let r;
  try { r = await app.api.get('/api/update?auto=1'); } catch { return; }
  if (!r?.checked || !r.newer || r.skipped) return;
  // do not interrupt: wait until no other dialog is open (at most 10 minutes)
  for (let i = 0; i < 120 && document.querySelector('dialog[open]'); i++) await new Promise(res => setTimeout(res, 5000));
  if (!document.querySelector('dialog[open]')) updateDialog(app, r);
}

export function updateDialog(app, r) {
  const auto = el('input', { type: 'checkbox', checked: app.info?.update?.auto ? true : null });
  auto.addEventListener('change', async () => {
    try {
      const s = await app.api.post('/api/update/settings', { auto: auto.checked });
      if (app.info) app.info.update = s;
    } catch (err) { toast(err.message, 'error', TITLE); }
  });
  const where = r.installed
    ? (r.file ? el('p', {}, el('span', { text: `Installer for this computer (${r.platform_label}): ` }), el('code', { text: r.file, translate: 'no', style: { whiteSpace: 'nowrap' } }))
      : el('p', { text: 'This release has no installer for this computer; the download page lists all files.' }))
    : el('p', { text: 'You run RCV from source: update it with git pull, or install RCV with an installer from the download page.' });
  const target = r.installed && r.download ? r.download : r.page;
  const d = dialog({
    title: 'Update available', width: 540,
    body: b => b.append(
      el('p', {}, el('b', { text: `RCV ${r.latest} is available. You have version ${r.current}.` })),
      r.date ? el('div', { class: 'desc', text: `Released on ${r.date}.` }) : null,
      where,
      el('div', { class: 'desc', text: 'The download opens in your web browser; RCV does not install anything by itself. Close RCV before installing the new version over it: your settings, plugins and predictions are kept.' }),
      el('label', { class: 'check', style: { marginTop: '10px' } }, auto, 'Check for updates when RCV starts')),
    actions: [
      { label: 'Skip this version', action: async () => {
        try { const s = await app.api.post('/api/update/settings', { skip: r.latest }); if (app.info) app.info.update = s; } catch { /* still close */ }
        d.close();
      } },
      { spacer: true },
      { label: 'Release notes', action: () => openExternal(app, r.notes || r.page) },
      { label: 'Later', close: true },
      { label: r.installed && r.download ? 'Download' : 'Open download page', primary: true, action: () => { openExternal(app, target); d.close(); } },
    ],
  });
  return d;
}
