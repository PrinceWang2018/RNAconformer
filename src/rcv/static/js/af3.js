// AlphaFold 3 in RCV: settings of a local installation (used as the prediction engine on the start page) and
// opening AlphaFold Server results (.zip). AlphaFold 3 is not part of RCV; it needs its own installation.
import { el } from './util.js';
import { icon } from './icons.js';
import { dialog, field, toast, busy } from './ui.js';

export const AF3_SERVER = 'https://alphafoldserver.com/';

export async function loadAf3Settings(app) {
  if (!app.api.online) return null;
  try { return await app.api.get('/api/af3/settings'); } catch { return null; }
}

export function af3Configured(s) { return !!(s && (s.command?.trim() || s.script?.trim())); }

export function af3Summary(s) {
  if (!af3Configured(s)) return 'Not set up yet: choose your AlphaFold 3 installation in the settings.';
  if (s.command?.trim()) return `Custom command: ${s.command.trim().slice(0, 80)}${s.command.length > 80 ? '…' : ''}`;
  const parts = [s.script];
  parts.push(s.db_dir ? 'with MSA search' : 'without MSA (no databases)');
  parts.push(`${s.seeds || 1} seed${(s.seeds || 1) > 1 ? 's' : ''} × 5 samples`);
  return parts.join(' · ');
}

// Settings dialog: Python + run_alphafold.py + model parameters (+ databases), or a custom command (Docker, WSL).
export async function af3SettingsDialog(app, onSaved) {
  const s = (await loadAf3Settings(app)) || {};
  const inp = (value, label, ph = '') => el('input', { class: 'input mono', value: value || '', 'aria-label': label, placeholder: ph, spellcheck: 'false' });
  const py = inp(s.python, 'Python of the AlphaFold 3 environment', 'e.g. /home/me/miniconda3/envs/af3/bin/python');
  const script = inp(s.script, 'run_alphafold.py', 'e.g. /home/me/alphafold3/run_alphafold.py');
  const models = inp(s.model_dir, 'Model parameters folder', 'folder with af3.bin (from Google DeepMind)');
  const dbs = inp(s.db_dir, 'Databases folder (optional)', 'empty: run without MSA search');
  const seeds = el('input', { class: 'input', type: 'number', min: 1, max: 20, value: s.seeds || 1, 'aria-label': 'Seeds' });
  const extra = inp(s.extra, 'Extra flags', 'e.g. --flash_attention_implementation=xla');
  const cmd = el('textarea', { class: 'input mono', rows: 3, spellcheck: 'false', 'aria-label': 'Custom command',
    placeholder: 'docker run --rm --gpus all -v {input_dir}:/root/af_input -v {output_dir}:/root/af_output -v /data/af3_models:/root/models alphafold3 python run_alphafold.py --json_path=/root/af_input/{json_name} --model_dir=/root/models --output_dir=/root/af_output --norun_data_pipeline' });
  cmd.value = s.command || '';
  const report = el('ul', { class: 'af3-check' });
  const pick = (target, folder) => el('button', { class: 'btn', type: 'button', hidden: !app.native?.available(), onclick: async () => {
    const p = folder ? await app.native.chooseFolder() : (await app.native.chooseFiles('any'))[0];
    if (p) target.value = p;
  } }, 'Choose…');
  const row = (label, control, button, help) => el('div', { class: 'field' }, el('label', { text: label }),
    el('div', { class: 'row', style: { gap: '8px' } }, control, button || null), help ? el('div', { class: 'help', text: help }) : null);
  const values = () => ({ python: py.value.trim(), script: script.value.trim(), model_dir: models.value.trim(), db_dir: dbs.value.trim(),
    seeds: +seeds.value || 1, extra: extra.value.trim(), command: cmd.value.trim() });
  dialog({
    title: 'AlphaFold 3 on this computer', width: 640,
    body: b => b.append(
      el('p', { class: 'help', text: 'RCV can run your own AlphaFold 3 installation (Linux with an NVIDIA GPU; model parameters are granted by Google DeepMind). RCV writes the input, starts AlphaFold 3 and opens every sampled model as one ensemble coloured by pLDDT. Without a local installation, use the AlphaFold Server and open its result .zip in card 2.' }),
      row('Python of the AlphaFold 3 environment', py, pick(py, false)),
      row('run_alphafold.py', script, pick(script, false)),
      row('Model parameters folder', models, pick(models, true)),
      row('Databases folder (optional)', dbs, pick(dbs, true), 'Empty: AlphaFold 3 runs without MSA search (no 600 GB databases needed; usually less accurate for RNA).'),
      el('div', { class: 'grid2' }, field('Seeds (5 models each)', seeds), field('Extra flags', extra)),
      el('details', { class: 'wl-adv', open: s.command ? true : null }, el('summary', { text: 'Custom command (Docker, WSL, cluster)' }),
        el('div', { class: 'help', text: 'Replaces the fields above. Placeholders: {json}, {json_name}, {input_dir}, {output_dir}, {model_dir}, {db_dir}. The program must write its output folder into {output_dir}.' }),
        cmd),
      report),
    actions: [
      { label: 'Check', action: async () => {
        const done = busy('Checking the AlphaFold 3 installation…');
        try {
          const r = await app.api.post('/api/af3/check', values());
          report.replaceChildren(...r.items.map(i => el('li', { class: i.ok === null ? 'muted' : i.ok ? 'state-ok' : 'state-bad' },
            icon(i.ok === false ? 'warn' : i.ok ? 'check' : 'info'), el('span', { text: i.text }))));
        } finally { done(); }
      } },
      { spacer: true }, { label: 'Cancel', close: true },
      { label: 'Save', primary: true, action: async () => {
        const saved = await app.api.post('/api/af3/settings', values());
        toast('AlphaFold 3 settings saved.', 'ok', 'AlphaFold 3');
        onSaved?.(saved);
        document.querySelector('dialog.dlg[open]')?.close();
      } },
    ],
  });
}

// AlphaFold Server result (.zip) -> ensemble (bytes are sent to the local server, which reads the models).
export async function openAf3Zip(app, file) {
  const done = busy(`Reading the AlphaFold 3 models in ${file.name}…`);
  try {
    const r = await fetch('/api/af3/open', { method: 'POST', body: file,
      headers: { 'X-Viewer-Token': app.api.token || '', 'X-Name': encodeURIComponent(file.name), 'Content-Type': 'application/zip' } });
    const out = await r.json().catch(() => ({}));
    if (!r.ok) throw new Error(out.error || `could not open ${file.name}`);
    await app.loadServerDocument(out.id);
    toast(`${out.n_conformations} AlphaFold 3 models of ${out.name}, best ranking score first. Colour by Confidence to see the pLDDT.`, 'ok', 'AlphaFold 3', 7000);
  } catch (err) { toast(err.message, 'error', 'Could not open the AlphaFold 3 result'); } finally { done(); }
}
