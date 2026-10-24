// __NAME__ — an RCV interface plugin (template written by Plugins ▸ New plugin).
// Edit this file, then Plugins ▸ Reload window. Guide: Help ▸ Plugin developer guide (docs/PLUGINS.md).
// A plugin is an ES module whose default export has id, name, version, description and activate(api).

export default {
  id: '__ID__',                       // lower-case, unique; prefixes this plugin's metrics, schemes and panels
  name: '__NAME__',
  version: '1.0',
  description: 'What this plugin does, in one sentence.',

  activate(api) {
    if (api.version < 2) { api.ui.toast('This plugin needs RCV 0.6 or newer.', 'error'); return; }
    const { el, fmt } = api.util;

    // 1. a metric: one number per conformation, plotted in the energy panel and exported with the CSV table
    api.metrics.register({
      key: 'gc_pairs', label: 'G–C pairs', unit: '', lower_is_better: false,
      description: 'Number of G–C base pairs in each conformation (template plugin).',
      compute: (ens, k) => ens.conformations[k].pairs.filter(p => p[2] === 'CG').length,
    });

    // 2. a console command: type "gccount" in the console (`)
    api.commands.register({
      name: 'gccount', args: '[conformation number]', help: 'Count G–C pairs in a conformation',
      run: args => {
        const ens = api.ensemble;
        if (!ens) return 'no ensemble open';
        const k = args[0] ? api.order[(+args[0] || 1) - 1] : api.conformation;
        const n = ens.conformations[k].pairs.filter(p => p[2] === 'CG').length;
        return `${ens.shortId(k)}: ${n} G–C pairs`;
      },
    });

    // 3. a side-panel tab that follows the shown conformation
    let box = null;
    const render = () => {
      if (!box) return;
      const ens = api.ensemble;
      box.replaceChildren(ens
        ? el('div', { class: 'desc', text: `${ens.shortId(api.conformation)}: ${ens.conformations[api.conformation].pairs.length} base pairs, ${fmt(ens.value(api.conformation, api.metric))} (${api.metric})` })
        : el('div', { class: 'muted', text: 'No ensemble open.' }));
    };
    api.ui.addPanel({ id: 'main', title: '__NAME__', render: host => { box = el('div', { class: 'section' }); host.append(box); render(); } });
    api.on('conformation', render);
    api.on('ensemble', render);

    // 4. a menu item (Plugins menu)
    api.ui.addMenuItem({ label: '__NAME__: label the selection', action: () => {
      for (const r of api.selection) api.annotations.add({ kind: 'label', residue: r, text: api.ensemble.resLabel(r), color: api.palette.colors[1] });
    } });
  },
};
