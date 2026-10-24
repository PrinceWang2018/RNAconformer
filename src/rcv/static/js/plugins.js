// Plugin host. A plugin is an ES module whose default export is
//   { id, name, version, description, activate(api) }
// Built-in plugins live in static/plugins/, user plugins in the plugins folder of the user data folder (*.js).
// See docs/RCV.md (section Plugins) for the API reference.
import { el, svg, fmt, fmtSigned, loadPref, savePref } from './util.js';
import { icon } from './icons.js';
import { registerScheme } from './colors.js';
import { toast, tooltip, tipRow, dialog } from './ui.js';
import { kabsch } from './ensemble.js';
import { SLOTS, NUCLEOTIDE } from './colors.js';
import { exportChart, standaloneSVG } from './export/chartexport.js';

export const API_VERSION = 2;          // 2 (RCV 0.6): commands, annotations, palette, surface, export; version 1 calls unchanged

export class PluginHost {
  constructor(app) {
    this.app = app;
    this.plugins = [];
    this.menuItems = [];
    this.disabled = new Set(loadPref('disabledPlugins', []));
  }

  async loadAll(list) {
    for (const p of list) await this.load(p);
  }

  async load(meta) {
    const rec = { url: meta.url, file: meta.file, builtin: !!meta.builtin, enabled: !this.disabled.has(meta.url), status: 'pending' };
    this.plugins.push(rec);
    if (!rec.enabled) { rec.status = 'disabled'; rec.name = meta.file?.split(/[\\/]/).pop(); return; }
    try {
      const mod = await import(new URL(meta.url, location.href).href);
      const def = mod.default;
      if (!def || typeof def.activate !== 'function' || !def.id) throw new Error('default export must have id and activate(api)');
      Object.assign(rec, { id: def.id, name: def.name || def.id, version: def.version || '', description: def.description || '' });
      await def.activate(this.api(rec));
      rec.status = 'active';
    } catch (err) {
      rec.status = 'error';
      rec.error = err.message || String(err);
      console.error('plugin failed', meta.url, err);
      toast(`${meta.file || meta.url}: ${rec.error}`, 'error', 'Plugin failed to load');
    }
  }

  setEnabled(url, on) {
    if (on) this.disabled.delete(url); else this.disabled.add(url);
    savePref('disabledPlugins', [...this.disabled]);
  }

  api(rec) {
    const app = this.app;
    const host = this;
    return Object.freeze({
      version: API_VERSION,
      // read-only state & navigation
      get ensemble() { return app.ens; },
      get conformation() { return app.state.conf; },
      get previous() { return app.state.prev; },
      get order() { return [...app.order]; },
      get metric() { return app.state.metric; },
      get selection() { return [...app.state.selection]; },
      get coords() { return app.display; },
      on: (event, fn) => app.on(event, fn),
      goTo: (k, opts) => app.goTo(k, opts),
      select: (residues, mode = 'set') => app.select(residues, mode),
      setMetric: key => app.setMetric(key),
      // UI
      ui: {
        addPanel: ({ id, title, render }) => app.addPluginTab({ id: `${rec.id}:${id}`, title, render }),
        addMenuItem: ({ label, action, kbd }) => { host.menuItems.push({ label, action, kbd, plugin: rec.name }); app.rebuildMenus(); },
        toast: (msg, kind) => toast(msg, kind || 'info', rec.name),
        dialog, tooltip, tipRow, icon,
      },
      // 3D view
      viewer: {
        addOverlay: (id, fn) => app.viewer.addOverlay(`${rec.id}:${id}`, fn),
        removeOverlay: id => app.viewer.removeOverlay(`${rec.id}:${id}`),
        requestRender: () => { app.viewer.needsRender = true; },
        project: (x, y, z) => { const v = app.viewer; const [w, h] = v.cssSize; return v.projector(v.camera, w, h)(x, y, z); },
      },
      colors: { registerScheme: def => { registerScheme({ ...def, id: `${rec.id}:${def.id}` }); app.refreshSchemes(); } },
      metrics: { register: def => app.addClientMetric({ ...def, key: `plugin.${rec.id}.${def.key}`, plugin: rec.name }) },
      // console commands: run(args, raw) returns a string, an array of strings, or nothing (may be async)
      commands: {
        register: ({ name, args = '', help = '', run }) => {
          if (!/^[a-z][a-z0-9_-]*$/.test(name || '')) throw new Error('command names: lower-case letters, digits, _ and -');
          if (app.console.cmds.commands.has(name)) throw new Error(`command "${name}" exists already`);
          app.console.cmds.def([name], { args, help: `${help} [plugin ${rec.name}]`, run });
        },
        run: line => app.console.run(line, { echo: app.console.visible }),
      },
      annotations: {
        add: a => app.annotations.add(a),
        list: () => app.annotations.list().map(a => ({ ...a })),
        remove: id => app.annotations.remove(id),
      },
      palette: { get colors() { return [...SLOTS]; }, get nucleotides() { return { ...NUCLEOTIDE }; }, get name() { return app.palettes.active.name; } },
      display: {
        surface: (on, residues = null) => app.surface.set(on, { residues }),
        transparency: (object, opacity) => app.setTransparency(object, opacity),
        color: (residues, hex) => { const d = app.displayOf(); for (const r of residues) d.colors.set(r, hex); app.applyDisplay(); },
      },
      export: { chart: (svgEl, format, name, opts) => exportChart(svgEl, format, name, opts), svgText: (svgEl, opts) => standaloneSVG(svgEl, opts) },
      util: { el, svg, fmt, fmtSigned, kabsch },
    });
  }
}
