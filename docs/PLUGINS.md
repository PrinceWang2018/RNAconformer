# Writing plugins for RCV

RCV (RNA Conformation Viewer) can be extended without changing its source code. There are two kinds of plugins:

| Kind | File | Runs | Good for |
|---|---|---|---|
| **Interface plugin** | `name.js` (JavaScript module) | in the RCV window | panels, charts, colour schemes, console commands, menu items, 3D overlays, annotations, metrics computed from the shown data |
| **Analysis plugin** | `name.py` (Python) | in the local RCV engine, when an ensemble is opened | metrics that need numpy/scipy or heavier computation per conformation |

Plugin API version: **3** (RCV 0.6.1). Everything from versions 1 (RCV 0.4–0.5) and 2 (RCV 0.6.0) still works.
Version 3 adds per-ensemble metrics, `api.metrics.unregister` (section 2.7) and saved plugin state (section 2.9).

---

## 1. Quick start

1. **Plugins ▸ New plugin from template…** writes a commented template (JavaScript or Python) into your plugin
   folder and shows where it is. The JavaScript template already adds a metric, a console command, a side panel and
   a menu item.
2. **Plugins ▸ Open plugin folder** opens the folder in Explorer / Finder.
3. Edit the file in any text editor.
4. **Plugins ▸ Reload window** (JavaScript) or restart RCV (Python). **Plugins ▸ Manage plugins** lists every plugin
   with its status and any error, and lets you switch plugins on and off.

Plugin folders (created on first use):

| System | Folder |
|---|---|
| Windows | `%APPDATA%\RCV\plugins` |
| macOS | `~/Library/Application Support/RCV/plugins` |
| Linux | `~/.local/share/rcv/plugins` |

Plugins in the folders of earlier versions (`RNAconfSwitch`, `~/.rnafold3d/plugins`) are loaded too.

**Developer mode** (View ▸ Preferences): F12 / Ctrl+Shift+I open the browser developer tools of the RCV window
(JavaScript console, debugger, network log), and Shift+right-click shows the browser's own menu. Leave it off for
normal use.

---

## 2. Interface plugins (JavaScript)

A plugin is an ES module whose default export is an object with `id`, `name`, `version`, `description` and
`activate(api)`:

```js
export default {
  id: 'hello',                 // lower-case and unique
  name: 'Hello',
  version: '1.0',
  description: 'Says hello when the conformation changes.',
  activate(api) {
    api.on('conformation', ({ from, to }) => api.ui.toast(`switched to ${api.ensemble.shortId(to)}`));
  },
};
```

`activate` runs once when RCV starts (it may be `async`). Keep it fast; do heavy work lazily (when a panel opens or a
command runs). Check the API version if you rely on newer features: `if (api.version < 3) return;`.

### 2.1 Reading state

| Property | Meaning |
|---|---|
| `api.version` | plugin API version (3) |
| `api.ensemble` | the open ensemble (`null` when none) — see section 3 |
| `api.conformation` | index of the shown conformation |
| `api.previous` | index of the previous conformation (`null` at first) |
| `api.order` | conformation indices in the current order (rank, energy, similarity path or input) |
| `api.metric` | key of the energy metric shown in the energy panel |
| `api.selection` | selected residue indices (0-based, 5′→3′) |
| `api.coords` | the coordinates being displayed (`Float32Array` x,y,z per atom; changes during morphs) |

### 2.2 Events

`api.on(event, fn)`:

| Event | Payload | When |
|---|---|---|
| `ensemble` | the ensemble | an ensemble is opened, closed or selected |
| `conformation` | `{ from, to }` | a switch starts |
| `morph` | `{ from, to, t }` every frame; `{ done: true }` at the end | during morphs and scrubbing |
| `metric`, `order`, `values` | — | the energy metric, the order or metric values changed |
| `selection` | residue list | the selection changed |
| `colors` | — | residue colours changed |
| `hover` | residue list | the mouse is over residues (3D, 2D window, sequence bar) |
| `tool` | tool id | the mouse tool changed (`select`, `pull`, `paint`, `label`, `measure`) |
| `annotations`, `measurements`, `display`, `surface` | — | the user edited them |
| `tab` | tab id | a side-panel tab was shown |
| `pathway` | the Pathway result (barriers, routes, bottleneck, rigid modules, timing) | a Pathway analysis finished (RCV 0.7) |
| `conditions` | `{ temperature, sodium, ligand, useProbing }` | the conditions of the 2D analyses changed |
| `probing` | the ensemble | probing data were loaded, refitted or removed |

### 2.3 Acting

| Call | Effect |
|---|---|
| `api.goTo(k, { duration, animate })` | morph to conformation `k` |
| `api.select(residues, mode)` | select residues; `mode` = `set`, `add` or `toggle` |
| `api.setMetric(key)` | show a metric in the energy panel |
| `api.commands.run('color orange, resi 5-9')` | run any console command (section 2.8) |

### 2.4 User interface

| Call | Effect |
|---|---|
| `api.ui.addPanel({ id, title, render(host) })` | adds a tab to the right panel; `render` receives the empty tab element |
| `api.ui.addMenuItem({ label, action, kbd })` | adds an item to the Plugins menu |
| `api.ui.toast(message, kind)` | short notification; `kind` = `info`, `ok` or `error` |
| `api.ui.dialog({ title, width, body, actions })` | modal dialog (`body` is an element or `b => b.append(...)`; `actions` = `[{ label, primary, close, action }]`) |
| `api.ui.tooltip.show(content, x, y)` / `.hide()` | hover tooltip |
| `api.ui.tipRow(value, label)` | formatted tooltip line |
| `api.ui.icon(name)` | one of RCV's stroke icons as an SVG element |
| `api.util.el(tag, attrs, ...children)` / `api.util.svg(...)` | create HTML / SVG elements (`text:` sets text, `on…:` adds listeners) |
| `api.util.fmt(value, digits)` | number formatting used across RCV |

Panels that draw charts should use `<svg class="chart">`: they then get RCV's publication styling, and right-click
offers PDF/SVG/PNG/TIFF export automatically. To export from your own button use
`api.export.chart(svgElement, 'pdf' | 'svg' | 'png' | 'tiff', name, { dpi })`.

### 2.5 3D view

| Call | Effect |
|---|---|
| `api.viewer.addOverlay(id, ({ coords, ens }) => primitives)` | draw lines and points that follow the structure; redrawn every frame and included in image exports. Primitives: `{ type: 'line', a: [x,y,z], b: [x,y,z], color, width, dash, label }`, `{ type: 'point', p: [x,y,z], color, r }` |
| `api.viewer.removeOverlay(id)` | remove it |
| `api.viewer.project(x, y, z)` | model → screen pixels `[x, y, depth]` |
| `api.viewer.requestRender()` | redraw after changing something an overlay reads |
| `api.annotations.add({ kind, residue, text, color, size, arrow, box })` | add a label (`kind: 'label'`), a region label (`kind: 'region', residues: [...]`) or a note (`kind: 'note', x, y` as fractions of the view); returns the annotation |
| `api.annotations.list()` / `.remove(id)` | read / delete annotations |
| `api.display.color(residues, '#rrggbb')` | custom colour on top of the colour scheme |
| `api.display.surface(true, residues)` | molecular surface (all residues, or the given ones) |
| `api.display.transparency(object, opacity)` | `cartoon`, `bases`, `atoms`, `surface`, `overlay` or `all`; opacity 0.05–1 |

### 2.6 Colour schemes

```js
api.colors.registerScheme({
  id: 'purines', label: 'Purine / pyrimidine',
  residueColors: (ens, ctx) => ens.residues.map(r => ('AG'.includes(r.name) ? api.palette.colors[0] : api.palette.colors[1])),
  legend: () => ({ title: 'Base type', items: [{ color: api.palette.colors[0], label: 'Purine' }, { color: api.palette.colors[1], label: 'Pyrimidine' }] }),
});
```

`residueColors` returns one `#rrggbb` per residue; `ctx.conf` is the conformation being coloured (schemes may differ
per conformation). A legend is a list of `items` or a ramp `{ title, ramp: [colours], min, max }`. Use
`api.palette.colors` (the active palette, chosen in Color ▸ Palettes) so your scheme follows the user's palette.

### 2.7 Metrics

```js
api.metrics.register({
  key: 'loop_size', label: 'Apical loop size', unit: 'nt', lower_is_better: null,
  description: 'Unpaired nucleotides between the last pair and the next pair.',
  compute: (ens, k) => { /* return a number or null for conformation k */ },
});
```

Metrics appear in the energy panel (profile, funnel, map), the Details table, the plot studio and the CSV export.
The key is stored as `plugin.<plugin id>.<key>`. Registering the same key again replaces the metric.

Version 3:

| Call | Effect |
|---|---|
| `api.metrics.register({ ..., ensemble: ens })` | the metric belongs to that ensemble only (for metrics the user creates on one ensemble, such as a tracked distance); without `ensemble` it applies to every ensemble, also ones opened later |
| `api.metrics.unregister(key)` | removes the metric (and its values) from every ensemble and stops computing it for new ones |
| `api.metrics.unregister(key, ens)` | removes it from that ensemble only |

When the removed metric is shown in the energy panel, RCV switches to the ensemble's default metric. `ensemble` may
be the ensemble object or its `id`.

### 2.8 Console commands

```js
api.commands.register({
  name: 'loopsize', args: '[conformation]', help: 'Apical loop size of a conformation',
  run: (args, raw) => `loop size: ${/* ... */ 4} nt`,       // string, array of strings, or nothing; may be async
});
```

Commands are listed by `help`, take comma-separated arguments, can be used in scripts (`run file.rcvs`) and are
recorded in the ensemble's command log (saved in `.rcv` files) like built-in commands. `api.commands.run(line)` runs
any command, which gives plugins the whole console vocabulary (selection language included).

### 2.9 Saved state (version 3)

Things the user creates with your plugin (tracked distances, marked regions, settings per ensemble) can be saved
with the ensemble, like RCV's own measurements and annotations:

```js
const notes = new WeakMap();                      // ensemble -> your state
api.session.register({
  save: ens => notes.get(ens) ?? null,            // a JSON-compatible value, or null for "nothing to save"
  load: (ens, value) => { notes.set(ens, value ?? []); /* redraw, re-register metrics, ... */ },
});
// after every change the user makes:
api.session.changed();
```

* The state is saved in `.rcv` files (File ▸ Save) and restored when they are opened; a workspace keeps it for the
  shown ensemble, like measurements. It is kept when a pull exploration extends the ensemble with new conformations
  (same residues, new ensemble object).
* It is part of **undo / redo**: `api.session.changed()` records a step; undo and redo call `load` with an earlier
  value. `load(ens, null)` means "no state" — clear everything for that ensemble.
* `load` receives saved data that may come from an older version of your plugin or a damaged file: check it and
  keep only what is valid. Residue indices refer to the same ensemble, so they can be stored as they are.
* Keep the state per ensemble (a `WeakMap` keyed by the ensemble object works well) and keep it small: it is copied
  for every undo step.
* If `load` registers metrics, use the same keys that were saved, so a saved view that shows one of them in the
  energy panel opens with it again. The built-in Distance tracker (`static/plugins/distance-tracker.js`) is a
  complete example.

---

## 3. The ensemble object

`api.ensemble` (read-only for plugins):

| Member | Content |
|---|---|
| `name`, `sequence` | ensemble name, sequence (A/C/G/U) |
| `n`, `nRes`, `nAtoms` | conformations, residues, atoms |
| `residues[r]` | `{ name, comp, num, chain, first_atom, n_atoms }` (`num` = author numbering, may include an insertion code) |
| `atomNames[a]`, `elements[a]`, `atomRes[a]` | atom name, element, residue index |
| `atomIndex[r][name]` | atom index of a named atom in residue `r` |
| `trace[r]`, `c1[r]`, `glyco[r]`, `phos[r]` | atom index of C4′, C1′, N1/N9, P (−1 if absent) |
| `coords[k]` | `Float32Array` of conformation `k` (superposed; x,y,z per atom) |
| `conformations[k]` | `{ id, label, rank, pairs: [[i, j, type, pseudoknot]], stacks, dbn, input_dbn, values: { metric: number }, validation, cluster, explored }` |
| `metrics`, `metricByKey[key]` | metric descriptions `{ key, label, unit, group, lower_is_better, description }` |
| `rmsd[a][b]`, `rmsf[r]`, `embedding[k]` | pairwise C4′ RMSD (Å), per-residue RMSF (Å), 2D similarity map |
| `value(k, key)` | metric value or `null` |
| `shortId(k)`, `resLabel(r)` | display names, e.g. `conf_003`, `G34` |
| `partner(k)` | `Int32Array` pairing partner per residue (−1 unpaired) |
| `pairInfo(k, r)` | `{ partner, type, pk }` or `null` |
| `residueCenter(coords, r)`, `residueDisplacement(a, b)` | geometry helpers |

Residue and atom indices are 0-based; selections, colours and annotations use residue indices.

---

## 4. Analysis plugins (Python)

```python
PLUGIN = {
    'id': 'contacts', 'name': 'Tertiary contacts',
    'metrics': [{'key': 'contacts', 'label': 'Tertiary contacts', 'unit': '', 'lower_is_better': False,
                 'description': 'Residue pairs > 3 apart with heavy atoms within 4 Å.'}],
}

def compute(sequence, residues):
    # sequence: str 5'->3'; residues: list (one per nucleotide) of {atom_name: numpy array (3,)}, superposed
    return {'contacts': 12}
```

* `compute` runs for every conformation when an ensemble is opened (also for `.rcv` conversion); return numbers or
  `None`. Exceptions are reported as a warning and never stop loading.
* The values become metrics `plugin.<id>.<key>` everywhere (energy panel, tables, plots, CSV, `.rcv` files). RCV
  remembers which plugins were active, so a prediction folder is re-evaluated when the plugin set changes.
* numpy and scipy are available; installer builds contain only the packages RCV itself needs.
* Example: `examples/viewer_plugins/tertiary_contacts.py`.

---

## 5. Testing, sharing, good practice

* **Debugging:** enable developer mode, open the developer tools (F12) and use `console.log`; load errors are listed in
  Plugins ▸ Manage plugins. For Python plugins, errors appear as warnings when the ensemble opens.
* **Automated tests:** keep computations in plain functions (no DOM) and test them with Node.js, as RCV does in
  `tests/js/` (load an ensemble document exported with `rcv export-document`).
* **Sharing:** a plugin is a single file — share it, or publish it with a licence header. RCV is GPL-3.0; plugins are
  your own work and may use any licence.
* **Style:** use `api.palette.colors` and `<svg class="chart">` so figures match the rest of RCV; label units; never
  block the window for more than a moment (use `setTimeout` / `await` to split long loops).
* **Safety:** plugins run with your user rights, like any program. Install only plugins you trust.
