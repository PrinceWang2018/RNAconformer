// Command console (bottom drawer of the 3D view): type commands, see results, reuse history, paste scripts.
// Every command run on an ensemble is kept with it (saved in .rcv files), so a figure can be reproduced.
import { el, loadPref, savePref } from '../util.js';
import { icon } from '../icons.js';
import { CommandSet, COLOR_NAMES } from './commands.js';
import { SELECTION_KEYWORDS } from './selection.js';
import { SCHEMES } from '../colors.js';

export class ConsolePanel {
  constructor(app, host) {
    this.app = app;
    this.cmds = new CommandSet(app);
    this.history = loadPref('consoleHistory', []);
    this.hpos = this.history.length;
    this.out = el('div', { class: 'cs-out', role: 'log', 'aria-live': 'polite', translate: 'no' });   // commands and results stay in English
    this.input = el('input', { class: 'cs-in', spellcheck: 'false', autocomplete: 'off', 'aria-label': 'Command', placeholder: 'Type a command — e.g. select resi 20-25 · color orange, sele · distance 12/N1, 30/N3 · help' });
    this.root = el('div', { class: 'console', hidden: true },
      el('div', { class: 'cs-head', 'data-win-head': '' }, icon('terminal'), el('b', { text: 'Console' }), el('span', { class: 'muted', text: 'commands for selections, colours, labels, measurements, views and exports' }),
        el('span', { class: 'grow' }),
        el('button', { class: 'ss-btn', title: 'List commands', onclick: () => this.run('help') }, icon('info')),
        el('button', { class: 'ss-btn', title: 'Clear output', onclick: () => this.out.replaceChildren() }, icon('trash')),
        el('span', { class: 'win-btns', 'data-win': 'console' })),
      this.out,
      el('div', { class: 'cs-line' }, el('span', { class: 'cs-prompt', text: 'RCV›', translate: 'no' }), this.input));
    host.append(this.root);
    this.input.addEventListener('keydown', e => this._key(e));
    this.input.addEventListener('paste', e => {
      const text = e.clipboardData?.getData('text') || '';
      if (!text.includes('\n')) return;
      e.preventDefault();
      this.runScript(text);
    });
    app.registerSessionPart('commands', {
      save: ens => (ens?.commandLog || []).slice(-500),
      load: (ens, list) => { ens.commandLog = Array.isArray(list) ? list.slice() : []; },
    });
    this.print('RCV console — type help for the list of commands. Selections: resi 5-12, resn G and paired, helix H2, switching, rmsf > 3, within 5 of resi 20, sele.', 'hint');
  }

  setVisible(on) {
    if (!on && this.app.wm?.maxed === 'console') this.app.wm.restore(false);
    this.root.hidden = !on;
    savePref('consoleOpen', on);
    this.app.viewer.resize();
    this.app.rebuildMenus?.();
    this.app._syncToolbar?.();
    if (on) setTimeout(() => this.input.focus(), 0);
  }
  toggle() { this.setVisible(this.root.hidden); }
  get visible() { return !this.root.hidden; }

  print(text, kind = 'out') {
    this.out.append(el('div', { class: `cs-${kind}`, text }));
    while (this.out.children.length > 600) this.out.firstChild.remove();
    this.out.scrollTop = this.out.scrollHeight;
  }

  async run(line, { echo = true, record = true } = {}) {
    const src = String(line).trim();
    if (!src) return true;
    if (echo) this.print(`RCV› ${src}`, 'cmd');
    if (record) {
      if (this.history[this.history.length - 1] !== src) this.history.push(src);
      if (this.history.length > 300) this.history.splice(0, this.history.length - 300);
      savePref('consoleHistory', this.history);
      this.hpos = this.history.length;
    }
    try {
      this.app.history?.commit();
      const out = await this.cmds.execute(src);
      for (const o of out) this.print(o.text, o.kind);
      const ens = this.app.ens;
      if (ens && !/^(help|\?|info|ls|list|count|pairs|measurements|echo|print)\b/i.test(src)) (ens.commandLog ||= []).push(src);
      return true;
    } catch (err) {
      this.print(err.message || String(err), 'err');
      return false;
    }
  }

  async runScript(text) {
    const lines = String(text).split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#'));
    this.print(`running ${lines.length} command${lines.length === 1 ? '' : 's'}`, 'hint');
    for (const l of lines) if (!(await this.run(l))) { this.print('script stopped at the error above', 'err'); break; }
  }

  _key(e) {
    if (e.key === 'Enter') { const v = this.input.value; this.input.value = ''; this.run(v); e.preventDefault(); return; }
    if (e.key === 'ArrowUp') { if (this.hpos > 0) { this.hpos--; this.input.value = this.history[this.hpos] || ''; } e.preventDefault(); return; }
    if (e.key === 'ArrowDown') { if (this.hpos < this.history.length) { this.hpos++; this.input.value = this.history[this.hpos] || ''; } e.preventDefault(); return; }
    if (e.key === 'Escape') { this.setVisible(false); this.app.viewer.canvas.focus(); e.preventDefault(); return; }
    if (e.key === 'Tab') { e.preventDefault(); this._complete(); }
    e.stopPropagation();
  }

  // Tab completion: command names first, then selection keywords, names, conformations, colours.
  _complete() {
    const v = this.input.value, pos = this.input.selectionStart ?? v.length;
    const before = v.slice(0, pos);
    const m = /([A-Za-z0-9_'+-]*)$/.exec(before);
    const word = m ? m[1] : '';
    const isCmd = !/\s/.test(before.trim()) && !before.includes(',');
    const app = this.app, ens = app.ens;
    let pool;
    if (isCmd) pool = this.cmds.names();
    else {
      pool = [...SELECTION_KEYWORDS, ...Object.keys(COLOR_NAMES), ...SCHEMES.map(s => s.id), 'cartoon', 'sticks', 'spheres', 'ballstick', 'labels', 'everything', 'ensemble', 'putty', 'sequence', '2d'];
      if (ens) pool.push(...ens.conformations.map(c => c.id), ...[...app.displayOf().named.keys()], ...app.scenes.list().map(s => s.name), ...(ens.doc.helix_families || []).map(f => f.name));
    }
    const cand = [...new Set(pool)].filter(p => p.toLowerCase().startsWith(word.toLowerCase())).sort();
    if (!cand.length) return;
    let ins;
    if (cand.length === 1) ins = cand[0] + (isCmd ? ' ' : '');
    else {
      let pre = cand[0];
      for (const c of cand) while (!c.toLowerCase().startsWith(pre.toLowerCase())) pre = pre.slice(0, -1);
      ins = pre.length > word.length ? pre : word;
      if (pre.length <= word.length) this.print(cand.slice(0, 40).join('  ') + (cand.length > 40 ? ' …' : ''), 'hint');
    }
    this.input.value = before.slice(0, before.length - word.length) + ins + v.slice(pos);
    const np = before.length - word.length + ins.length;
    this.input.setSelectionRange(np, np);
  }
}
