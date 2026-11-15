// Recipes: named analyses — prediction settings and console steps — kept as small .rcvrecipe files that a lab can
// share. Tools > Recipes… lists the built-in and your own recipes; run one on the open ensemble, choose one for the
// next prediction (start page, Advanced settings: its steps run on the result), make one from the commands used
// on this ensemble, or import / export the file. The local server stores them (rcv/recipes.py).
import { el, slug, saveBlob, readFileText } from './util.js';
import { icon } from './icons.js';
import { dialog, toast, askConfirm } from './ui.js';
import { CommandError } from './console/commands.js';

const SAMPLING = [['', '(keep the form)'], ['fast', 'Fast'], ['standard', 'Standard'], ['thorough', 'Thorough']];

export class Recipes {
  constructor(app) { this.app = app; this.cache = null; }

  async list(refresh = false) {
    if (!this.app.api.online) return [];
    if (!this.cache || refresh) this.cache = await this.app.api.get('/api/recipes');
    return this.cache.recipes;
  }
  async find(key) {
    const k = String(key || '').trim().toLowerCase();
    return (await this.list()).find(r => r.id.toLowerCase() === k || r.name.toLowerCase() === k) || null;
  }

  // Runs the steps on the open ensemble through the console (plain-language lines work too); stops at an error.
  async run(recipe) {
    const app = this.app;
    if (!app.ens) { toast('Open an ensemble first: a recipe runs on the open ensemble.', 'error', 'Recipe'); return false; }
    if (!recipe.steps?.length) { toast(`"${recipe.name}" has no steps to run (it only holds prediction settings).`, 'info', 'Recipe'); return true; }
    app.console.print(`recipe "${recipe.name}" (${recipe.steps.length} steps)`, 'hint');
    let k = 0;
    for (const step of recipe.steps) {
      k++;
      if (!(await app.console.run(step, { record: false }))) {
        app.console.setVisible(true);
        toast(`Stopped at step ${k} of ${recipe.steps.length}: ${step} (see the console).`, 'error', recipe.name, 8000);
        return false;
      }
    }
    toast(`${recipe.steps.length} steps run on ${app.ens.name}. The commands are in the console and saved with the ensemble.`, 'ok', recipe.name, 5000);
    return true;
  }
}

function recipeFile(r) {
  const { rcv_recipe, name, description, author, prediction, steps } = r;
  return JSON.stringify({ rcv_recipe: rcv_recipe || 1, name, description, author, prediction, steps }, null, 2) + '\n';
}

// Tools > Recipes… (its own window: list on the left, the chosen recipe on the right)
export async function recipesDialog(app, selectId = null) {
  const R = app.recipes;
  let list = [];
  try { list = await R.list(true); } catch (err) { toast(err.message, 'error', 'Recipes'); return; }
  let cur = list.find(r => r.id === selectId) || list[0] || null;
  let dirty = false;

  const listBox = el('div', { class: 'rc-list', role: 'listbox', 'aria-label': 'Recipes' });
  const detail = el('div', { class: 'rc-detail' });
  const nameIn = el('input', { class: 'input', 'aria-label': 'Name' });
  const descIn = el('textarea', { class: 'input', rows: 3, 'aria-label': 'Description' });
  const authorIn = el('input', { class: 'input', 'aria-label': 'Author' });
  const modelsIn = el('input', { class: 'input', type: 'number', min: 1, max: 50, placeholder: '–', 'aria-label': 'Conformations to keep' });
  const sampIn = el('select', { class: 'select', 'aria-label': 'Sampling' }, ...SAMPLING.map(([v, t]) => el('option', { value: v, text: t })));
  const seedIn = el('input', { class: 'input', type: 'number', placeholder: '–', 'aria-label': 'Random seed' });
  const tplIn = el('select', { class: 'select', 'aria-label': 'Templates' }, ...[['', '(keep the form)'], ['auto', 'Use templates'], ['off', 'No templates']].map(([v, t]) => el('option', { value: v, text: t })));
  const stepsIn = el('textarea', { class: 'input mono rc-steps', rows: 12, spellcheck: 'false', 'aria-label': 'Steps' });
  const note = el('div', { class: 'desc' });
  for (const c of [nameIn, descIn, authorIn, modelsIn, sampIn, seedIn, tplIn, stepsIn]) c.addEventListener('input', () => { dirty = true; sync(); });

  const read = () => ({
    rcv_recipe: 1, name: nameIn.value.trim(), description: descIn.value.trim(), author: authorIn.value.trim(),
    prediction: Object.fromEntries(Object.entries({ models: modelsIn.value ? +modelsIn.value : null, sampling: sampIn.value || null,
      seed: seedIn.value ? +seedIn.value : null, templates: tplIn.value || null }).filter(([, v]) => v !== null)),
    steps: stepsIn.value.split(/\r?\n/).map(s => s.trim()).filter(Boolean),
  });
  const fill = r => {
    const p = r?.prediction || {};
    nameIn.value = r?.name || ''; descIn.value = r?.description || ''; authorIn.value = r?.author || '';
    modelsIn.value = p.models ?? ''; sampIn.value = p.sampling || ''; seedIn.value = p.seed ?? ''; tplIn.value = p.templates || '';
    stepsIn.value = (r?.steps || []).join('\n');
    const ro = !!r?.builtin;
    for (const c of [nameIn, descIn, authorIn, modelsIn, sampIn, seedIn, tplIn, stepsIn]) c.disabled = !r || ro;
    note.textContent = r?.error ? `This file could not be read: ${r.error}` : ro ? 'A built-in recipe. Duplicate it to change it.' : r?.path ? `Saved in ${r.path}` : 'Not saved yet.';
    dirty = false;
  };
  const drawList = () => {
    const group = (title, items) => items.length ? [el('div', { class: 'rc-group', text: title }), ...items.map(r => el('button', {
      class: `rc-item${cur && r.id === cur.id ? ' on' : ''}`, type: 'button', role: 'option', 'aria-selected': cur && r.id === cur.id ? 'true' : 'false',
      onclick: async () => { if (dirty && !(await askConfirm('Recipes', 'Discard the changes to this recipe?', 'Discard'))) return; cur = r; fill(r); drawList(); sync(); } },
      el('b', { text: r.name, translate: r.builtin ? null : 'no' }), el('span', { class: 'muted', text: r.error ? 'cannot be read' : `${r.steps?.length || 0} steps${Object.keys(r.prediction || {}).length ? ' · prediction settings' : ''}` })))] : [];
    listBox.replaceChildren(...group('Built-in', list.filter(r => r.builtin)), ...group('My recipes', list.filter(r => !r.builtin)));
  };
  const reload = async id => { list = await R.list(true); cur = list.find(r => r.id === id) || list[0] || null; fill(cur); drawList(); sync(); };
  const saveCur = async () => {
    const body = read();
    if (!body.name) { toast('Give the recipe a name.', 'error', 'Recipe'); return null; }
    const saved = await app.api.post('/api/recipes', { recipe: body, id: cur && !cur.builtin && cur.path ? cur.id : null });
    toast(`Saved in ${saved.path}`, 'ok', saved.name);
    await reload(saved.id);
    return saved;
  };
  const newRecipe = r => { cur = { id: '', builtin: false, ...r }; fill(cur); dirty = true; drawList(); sync(); nameIn.focus(); nameIn.select(); };

  const btn = (label, ic, fn, cls = 'btn') => el('button', { class: cls, type: 'button', onclick: async () => { try { await fn(); } catch (err) { toast(err.message || String(err), 'error', 'Recipe'); } } }, ic ? icon(ic) : null, label);
  const runBtn = btn('Run on this ensemble', 'play', async () => { const r = dirty || !cur.path && !cur.builtin ? { ...cur, ...read() } : cur; d.close(); await R.run(r); }, 'btn primary');
  const nextBtn = btn('Use for the next prediction', 'flask', async () => {
    const r = dirty ? { ...cur, ...read() } : cur;
    app.welcome.setRecipe(r);
    d.close();
    app.showPage('welcome');
    toast(`The next prediction uses "${r.name}": its settings are in the form, and its steps run on the result.`, 'info', 'Recipe', 6000);
  });
  const saveBtn = btn('Save', 'check', saveCur);
  const dupBtn = btn('Duplicate', 'copy', () => newRecipe({ ...read(), name: `${cur.name} (copy)` }));
  const expBtn = btn('Export file…', 'download', () => saveBlob(new Blob([recipeFile(dirty ? read() : cur)], { type: 'application/json' }), `${slug(cur.name)}.rcvrecipe`,
    [{ description: 'RCV recipe', accept: { 'application/json': ['.rcvrecipe'] } }]));
  const delBtn = btn('Delete', 'trash', async () => {
    if (!(await askConfirm('Delete recipe', `Delete the recipe "${cur.name}"? Its file is removed.`, 'Delete'))) return;
    await app.api.del(`/api/recipes/${encodeURIComponent(cur.id)}`);
    await reload(null);
  });
  const fileIn = el('input', { type: 'file', accept: '.rcvrecipe,.json', hidden: true });
  fileIn.addEventListener('change', async () => {
    const f = fileIn.files?.[0]; fileIn.value = '';
    if (!f) return;
    try { const r = await app.api.post('/api/recipes/import', { text: await readFileText(f) }); toast(`Imported into ${r.path}`, 'ok', r.name); await reload(r.id); }
    catch (err) { toast(err.message, 'error', f.name); }
  });
  const sync = () => {
    const has = !!cur, ens = !!app.ens;
    runBtn.disabled = !has || !ens || !!cur?.error;
    runBtn.title = ens ? '' : 'Open an ensemble first';
    nextBtn.disabled = !has || !!cur?.error;
    saveBtn.disabled = !has || cur.builtin || !dirty;
    dupBtn.disabled = !has || !!cur?.error;
    expBtn.disabled = !has || !!cur?.error;
    delBtn.disabled = !has || cur.builtin || !cur.path;
  };

  const f = (label, c) => el('div', { class: 'field' }, el('label', { text: label }), c);
  detail.append(
    f('Name', nameIn), f('Description', descIn),
    el('div', { class: 'section-title', text: 'Prediction settings (used when the recipe is chosen for a prediction)' }),
    el('div', { class: 'grid2' }, f('Conformations to keep', modelsIn), f('Sampling', sampIn)),
    el('div', { class: 'grid2' }, f('Random seed', seedIn), f('Templates', tplIn)),
    el('div', { class: 'field' }, el('label', { text: 'Steps: one console command per line, run in order on the ensemble (plain English works too)' }), stepsIn),
    f('Author', authorIn), note,
    el('div', { class: 'rc-actions' }, el('div', { class: 'rc-act-row' }, runBtn, nextBtn), el('div', { class: 'rc-act-row' }, saveBtn, dupBtn, expBtn, delBtn)));
  fill(cur); drawList(); sync();

  const d = dialog({
    title: 'Recipes', width: 980,
    body: b => b.append(
      el('div', { class: 'desc', text: 'A recipe is a named analysis: prediction settings and the steps that follow — conditions, colours, tabs, views, comparisons. Share it as a small .rcvrecipe file so everyone in the lab looks at an RNA the same way.' }),
      el('div', { class: 'rc-grid' }, el('div', { class: 'rc-left' }, listBox,
        el('div', { class: 'rc-left-actions' },
          btn('New recipe', 'plus', () => newRecipe({ name: 'New recipe', steps: [], prediction: {} })),
          btn('From this ensemble’s commands', 'terminal', () => {
            const log = app.ens?.commandLog || [];
            if (!log.length) { toast('No commands have been run on this ensemble yet (console, plain language or a recipe).', 'info', 'Recipe'); return; }
            newRecipe({ name: `${app.ens.name} analysis`, description: `Made from the commands used on ${app.ens.name}.`, steps: log.slice(-MAX_FROM_LOG), prediction: {} });
          }),
          btn('Import file…', 'upload', () => fileIn.click()),
          app.native?.available() ? btn('Open recipes folder', 'folder', () => app.native.reveal(`${R.cache?.folder}/x`)) : null)),
        detail), fileIn),
    actions: [{ label: 'Close', close: true }],
  });
  return d;
}
const MAX_FROM_LOG = 200;

// console: recipe list | recipe <name>
export function defineRecipeCommand(cmds, app) {
  cmds.def(['recipe'], { args: 'name | list', help: 'Run a recipe (Tools > Recipes) on the open ensemble, or list the recipes.', run: async (args, raw) => {
    const want = raw.trim();
    const all = await app.recipes.list(true);
    if (!want || want.toLowerCase() === 'list') return all.map(r => ({ text: `${r.name.padEnd(28)} ${r.steps?.length || 0} steps${r.builtin ? ' (built-in)' : ''}`, kind: 'out' }));
    const r = all.find(x => x.id.toLowerCase() === want.toLowerCase() || x.name.toLowerCase() === want.toLowerCase());
    if (!r) throw new CommandError(`no recipe "${want}" (recipe list)`);
    await app.recipes.run(r);
    return { text: `recipe "${r.name}" done`, kind: 'out' };
  } });
}
