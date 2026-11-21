// Plain-language front end of the console: turns a sentence into RCV commands, locally, with patterns (no
// language model, nothing leaves the computer). "colour by flexibility" -> color flexibility; "show the switch at
// 25 degrees with 150 mM salt" -> conditions 25, 150 mM + tab switch; "compare models 1 and 3" -> align 1, 3.
// The console runs a line as a command first and asks this module only when that fails, then shows and records
// the commands it ran (so the command log stays a reproducible script). English only, like the commands.
//
// translatePlain(text, { schemes, ligands }) -> { commands: [...] } or null when no pattern fits.

const COLORS = ['red', 'green', 'darkgreen', 'blue', 'navy', 'orange', 'yellow', 'gold', 'purple', 'violet', 'magenta', 'pink',
  'cyan', 'teal', 'brown', 'white', 'black', 'gray', 'grey', 'lightgray', 'lightgrey', 'salmon', 'wheat'];
const COLOR_RE = COLORS.join('|');
const NUM = '(-?\\d+(?:\\.\\d+)?)';
const NT = { g: 'G', guanine: 'G', guanines: 'G', gs: 'G', a: 'A', adenine: 'A', adenines: 'A', as: 'A', c: 'C', cytosine: 'C',
  cytosines: 'C', cytidine: 'C', cs: 'C', u: 'U', uracil: 'U', uracils: 'U', uridine: 'U', uridines: 'U', us: 'U' };
const TABS = { energy: 'energy', energies: 'energy', switch: 'switch', switching: 'switch', pathway: 'pathway', path: 'pathway',
  transition: 'pathway', details: 'details', detail: 'details' };
// colour schemes by the words people use for them
const SCHEME_WORDS = [
  [/modif|methyl|pseudouridin|inosin|edit|oxo/, 'modified'], [/flexib|rmsf|mobil|motion|dynamic/, 'flexibility'], [/switch/, 'switching'], [/frequen|ensemble pair/, 'pairfreq'],
  [/pair/, 'pairing'], [/helix|helices|stem/, 'helix'], [/element|secondary|loop/, 'elements'],
  [/displace|moved|movement|change/, 'displacement'], [/position|rainbow|5.?\s*(to|->)\s*3|order/, 'position'],
  [/classic/, 'nucleotide-classic'], [/nucleotide|base|sequence|residue type|type/, 'nucleotide'], [/uniform|single|one colou?r|plain/, 'uniform'],
];

// "residues 10 to 20", "nucleotide 5", "the guanines", "switching residues", "helix H2" -> selection expression
export function residuePhrase(s) {
  const t = ` ${String(s || '').toLowerCase().replace(/[.!?]+$/, '').trim()} `;
  let m;
  if ((m = /\b(?:residues?|nucleotides?|nts?|bases?|positions?|resi)\s+(-?\d+)\s*(?:-|–|to|through|thru|until)\s*(-?\d+)/.exec(t))) return `resi ${m[1]}-${m[2]}`;
  if ((m = /\b(?:residues?|nucleotides?|nts?|bases?|positions?|resi)\s+(-?\d+(?:\s*(?:,|and)\s*-?\d+)+)/.exec(t))) return `resi ${m[1].split(/\s*(?:,|and)\s*/).join('+')}`;
  if ((m = /\b(?:residue|nucleotide|nt|base|position|resi)\s+(-?\d+)/.exec(t))) return `resi ${m[1]}`;
  if ((m = /\b(?:helix|stem)\s+([a-z]?\d+)\b/.exec(t))) return `helix ${m[1].toUpperCase()}`;
  if (/\bswitch(ing)?\b/.test(t)) return 'switching';
  if (/\bpseudo-?knot/.test(t)) return 'pk';
  if (/\b(unpaired|single[- ]stranded|loops?)\b/.test(t)) return 'unpaired';
  if (/\b(paired|base[- ]paired|double[- ]stranded|helices|stems)\b/.test(t)) return 'paired';
  if (/\bflexible\b/.test(t)) return 'rmsf > 3';
  if (/\b(selection|selected|sele)\b/.test(t)) return 'sele';
  if (/\b(first) (residue|nucleotide)\b/.test(t)) return 'first';
  if (/\b(last) (residue|nucleotide)\b/.test(t)) return 'last';
  for (const [w, n] of Object.entries(NT)) if (new RegExp(`(^|\\s)(all\\s+|the\\s+|every\\s+)?${w}('?s)?(\\s|$)`).test(t) && (w.length > 2 || /\b(all|the|every)\s/.test(t) || t.trim() === w)) return `resn ${n}`;
  if (/\b(everything|all|whole|molecule|entire|it)\b/.test(t)) return 'all';
  return null;
}

function conditionsOf(c, ligands) {
  const parts = [];
  let m;
  if (/\bstandard conditions\b/.test(c)) return ['standard'];
  // the ligand first, and taken out of the text, so "theophylline 2 mM" is not read as salt
  for (const name of Object.keys(ligands || {})) {
    if (new RegExp(`\\b(?:no|without)\\s+(?:the\\s+)?(?:${name}|ligand)\\b`).test(c)) { parts.push('no ligand'); c = c.replace(new RegExp(name, 'g'), ''); continue; }
    if ((m = new RegExp(`\\b${name}\\b(?:\\s+(?:at\\s+)?${NUM}\\s*(µm|um|nm|mm)?)?`).exec(c))) {
      const unit = { um: 'µM', 'µm': 'µM', nm: 'nM', mm: 'mM' }[m[2]] || 'µM';
      parts.push(`${name} ${m[1] ? +m[1] : 10} ${unit}`);
      c = c.replace(m[0], ' ');
    }
  }
  if (!parts.length && /\b(no|without)\s+(the\s+)?ligand\b/.test(c)) parts.push('no ligand');
  if (/\bphysiological\b/.test(c)) parts.unshift('37', '150 mM');
  else if (/\broom temperature\b/.test(c)) parts.unshift('25');
  else if (/\bbody temperature\b/.test(c)) parts.unshift('37');
  if ((m = new RegExp(`${NUM}\\s*(?:°\\s*c?|degrees?(?:\\s*c(?:elsius)?)?|deg\\b|c\\b)`).exec(c))) parts.unshift(String(+m[1]));
  if ((m = new RegExp(`${NUM}\\s*(mm|m)\\b`).exec(c)) && (m[2] === 'mm' || /salt|na\b|sodium|nacl|kcl/.test(c))) parts.push(`${+m[1]} ${m[2] === 'mm' ? 'mM' : 'M'}`);
  return parts;
}

function schemeOf(phrase, schemes) {
  const p = phrase.toLowerCase();
  for (const s of schemes || []) if (p.trim() === s.id || p.includes(s.label.toLowerCase())) return s.id;
  for (const [re, id] of SCHEME_WORDS) if (re.test(p) && (!schemes || schemes.some(s => s.id === id))) return id;
  return null;
}

const nums = s => [...String(s).matchAll(/(?:^|[^\w/])(\d+)(?![\w/])/g)].map(m => m[1]);
const atomSpec = s => s.replace(/\/([a-z0-9']+)/g, (_, a) => '/' + a.toUpperCase());

// Each rule: [regex, (match, clause, opts) -> command(s) | null]. The first rule that gives commands wins.
const RULES = [
  [/^(help|commands|what can (i|you)|how do i|\?)/, () => 'help'],
  [/\b(?:fetch|download|get|load|open)\b.*?\b(?:pdb\s*(?:entry|id|structure)?\s*)?([0-9][a-z0-9]{3})\b/, m => `fetch ${m[1].toUpperCase()}`],
  // view motion before conditions ("turn 90 degrees" is not a temperature)
  [/\b(?:stop|end|halt)\s+(?:the\s+)?(?:rotat|spin|turn)/, () => 'spin off'],
  [/\b(?:rotate|turn)\b(?:\s+(?:the\s+)?(?:view|molecule|structure|it))?(?:\s+by)?\s+(-?\d+)/, (m, c) => {
    const ax = /\b(?:around|about|on|along)\s+(?:the\s+)?([xyz])\b/.exec(c) || /\b([xyz])[- ]axis\b/.exec(c);
    return `turn ${ax ? ax[1] : 'y'}, ${m[1]}`;
  }],
  [/\b(?:spin(?:ning)?|auto-?rotat\w*|keep rotating|start rotating|rotate continuously|turntable)\b/, () => 'spin on'],
  // movies and images before pathways and conformations
  [/\b(video|movie|animation|mp4|webm)\b/, (m, c) => {
    const n = nums(c);
    if (/\b(rotat|spin|turntable|360)/.test(c)) return 'video rotate';
    if (n.length >= 2) return `video transition ${n[0]}, ${n[1]}`;
    return 'video tour';
  }],
  [/\b(?:png|picture|image|snapshot|screenshot|figure|photo)\b/, (m, c) => /\b(save|export|take|make|create|get|render)\b/.test(c) ? 'png' : null],
  [/\b(?:save|export|write)\b/, (m, c) => {
    if (/\bworkspace|session\b/.test(c)) return 'save workspace';
    if (/\bpdb\b/.test(c)) return /\b(all|every|each)\b/.test(c) ? 'save pdb all' : 'save pdb';
    if (/\b(csv|table|energies|spreadsheet)\b/.test(c)) return 'save csv';
    if (/\b(rcv|file|ensemble|everything|it|results?)\b/.test(c)) return 'save rcv';
    return null;
  }],
  // comparisons between conformations
  [/\bmorph\b/, (m, c) => { const n = nums(c); return n.length >= 2 ? `morph ${n[0]}, ${n[1]}` : null; }],
  [/\b(?:pathway|refold|transition|how does|go from|turn into|path)\b/, (m, c) => {
    const n = nums(c);
    return n.length >= 2 ? `pathway ${n[0]}, ${n[1]}` : (/\b(tab|panel)\b/.test(c) ? 'tab pathway' : null);
  }],
  [/\b(?:compare|align|superpose|superimpose|rmsd|difference|differ)\b/, (m, c) => { const n = nums(c); return n.length >= 2 ? `align ${n[0]}, ${n[1]}` : null; }],
  [/\b(?:distance|how far|measure)\b/, (m, c) => {
    const at = [...c.matchAll(/\b(\d+(?:\/[a-z0-9']+)?)(?![\w])/g)].map(x => atomSpec(x[1]));
    if (at.length < 2) return null;
    return /\bangle\b/.test(c) && at.length >= 3 ? `angle ${at.slice(0, 3).join(', ')}` : `distance ${at[0]}, ${at[1]}`;
  }],
  // conditions of the 2D analyses (and "show the switch at 25 °C" -> also the Switch tab)
  [/(°|degrees?\b|\bdeg\b|\d\s*c\b|\bmm\b|\d\s*m\b|\bsalt\b|\bsodium\b|\bnacl\b|\bkcl\b|\bligand\b|\bphysiological\b|\broom temperature\b|\bbody temperature\b|\bstandard conditions\b|\btheophylline\b)/, (m, c, o) => {
    const p = conditionsOf(c, o.ligands);
    if (!p.length) return null;
    const out = [`conditions ${p.join(', ')}`];
    const tab = /\b(switch|switching)\b/.test(c) ? 'switch' : /\b(pathway|transition)\b/.test(c) ? 'pathway' : null;
    if (tab) out.push(`tab ${tab}`);
    return out;
  }],
  // playback and conformations
  [/\b(?:base pairs|which pairs|list (?:the )?pairs|what pairs|pairing)\b/, (m, c) => { const n = nums(c); return n.length ? `pairs ${n[0]}` : 'pairs'; }],
  [/^(?:stop|pause|halt)\b/, () => 'stop'],
  [/\b(?:play|animate|tour|cycle through|go through)\b/, () => 'play'],
  [/\bnext\b/, () => 'next'],
  [/\bprevious\b|\bprev\b|\bgo back\b/, () => 'prev'],
  [/\bfirst (?:conformation|conformer|model|structure|one)\b/, () => 'conf first'],
  [/\blast (?:conformation|conformer|model|structure|one)\b/, () => 'conf last'],
  [/\b(?:conformation|conformer|model|structure|state)\s*(?:number\s*)?#?(\d+)\b/, m => `conf ${m[1]}`],
  [/\b(?:sort|order|rank)\b.*\b(energy|rank|similar|similarity|path|input|original)\b/, m => `order ${{ similar: 'path', similarity: 'path', original: 'input' }[m[1]] || m[1]}`],
  // colours
  [/\b(?:reset|remove|clear|undo)\b.*\bcolou?rs?\b|\bdefault colou?rs?\b/, () => 'color default, all'],
  [/\b(?:background|bg)\b/, (m, c) => {
    if (/\bdark\b/.test(c)) return 'bg black';
    if (/\blight\b/.test(c)) return 'bg white';
    const k = new RegExp(`\\b(${COLOR_RE})\\b`).exec(c);
    return k ? `bg ${k[1]}` : null;
  }],
  [/\bcolou?r(?:ed|s)?\b.*?\bby\b\s+(?:the\s+)?(.+)$/, (m, c, o) => { const s = schemeOf(m[1], o.schemes); return s ? `color ${s}` : null; }],
  [new RegExp(`\\b(?:colou?r|make|paint|turn)\\b\\s+(.*?)\\s*(?:in\\s+|to\\s+)?\\b(${COLOR_RE})\\b`), m => {
    const sel = residuePhrase(m[1]) || (/^(it|this|them|everything|all)?$/.test(m[1].trim()) ? 'all' : null);
    return sel ? `color ${m[2]}, ${sel}` : null;
  }],
  // labels
  [/\b(?:remove|clear|hide|delete|no|turn off)\b.*\blabels?\b/, () => 'hide labels'],
  [/\blabels?\b.*\bevery\s+(?:10th|tenth|ten)\b/, () => 'set labels, every10'],
  [/\blabels?\b.*\b(?:all|every|each)\b\s+(?:residues?|nucleotides?|bases?)/, () => 'set labels, all'],
  [/\blabel\b\s+(.*)/, m => { const s = residuePhrase(m[1]); return s ? `label ${s}` : null; }],
  // view
  [/\breset\b.*\bview\b|\bfit\b.*\b(window|screen)\b|\bzoom out\b/, () => 'reset'],
  [/\b(?:zoom|focus|centre|center)\b(?:\s+in)?(?:\s+(?:on|to|at))?\s*(.*)$/, m => { const s = residuePhrase(m[1]); return s && s !== 'all' ? `zoom ${s}` : 'zoom'; }],
  // representations and parts of the window
  [/\b(?:open|show|go to|switch to|display)\b.*\b(energy|energies|switch|switching|pathway|details)\b\s*(?:tab|panel|analysis|plot|chart)?/, (m, c) => (/\b(tab|panel|analysis|plot|chart)\b/.test(c) || /^(open|go to|switch to)\b/.test(c) ? `tab ${TABS[m[1]]}` : null)],
  [/\b(show|display|turn on|add|hide|remove|turn off|close)\b/, (m, c) => {
    const on = /\b(show|display|turn on|add)\b/.test(m[1]);
    const v = on ? 'show' : 'hide';
    const sel = residuePhrase(c.replace(m[1], ''));
    const forSel = sel && !['all', 'paired', 'unpaired'].includes(sel) ? `, ${sel}` : '';
    if (/\b(2d|2-d|secondary structure|two-dimensional)\b/.test(c)) return `${v} 2d`;
    if (/\bsequence\b/.test(c)) return `${v} sequence`;
    if (/\b(ensemble|all conformations|every conformation|overlay|ghosts|superposition)\b/.test(c)) return `${v} ensemble`;
    if (/\b(putty|flexibility)\b/.test(c)) return `${v} putty`;
    if (/\bsurface\b/.test(c)) return `${v} surface${forSel}`;
    if (/\b(sticks|atoms|bases)\b/.test(c)) return `${v} sticks${forSel}`;
    if (/\b(spheres|space ?fill)\b/.test(c)) return `${v} spheres${forSel}`;
    if (/\blabels?\b/.test(c)) return `${v} labels${forSel}`;
    if (/\b(cartoon|backbone)\b/.test(c)) return `${v} cartoon${forSel}`;
    return null;
  }],
  [/\b(?:as|use|switch to|change to|draw|in)\b.*?\b(cartoon|sticks|spheres|space ?fill|ball[- ]and[- ]stick|ballstick|ball\+stick)\b/, m => {
    const r = m[1].replace(/\s/g, '');
    return `rep ${{ spacefill: 'spheres', 'ball-and-stick': 'ballstick', ballandstick: 'ballstick', 'ball+stick': 'ballstick' }[r] || r}`;
  }],
  [/\b(studio|illustrative|glossy|depth|plain)\b.*\b(style|look|rendering|mode)\b|\b(style|look|rendering)\b.*\b(studio|illustrative|glossy|depth|plain)\b/, m => `style ${m[1] || m[4]}`],
  // selections and questions
  [/\b(?:clear|remove|cancel)\b.*\bselection\b|\bdeselect\b|\bselect nothing\b|\bunselect\b/, () => 'deselect'],
  [/\b(?:select|pick|highlight|mark|find)\b\s+(.*)/, m => { const s = residuePhrase(m[1]); return s ? `select ${s}` : null; }],
  [/\bhow many\b\s+(.*)/, m => { const s = residuePhrase(m[1]); return s ? `count ${s}` : null; }],
  [/\b(?:info|summary|summarize|tell me about|what is (?:this|open|shown)|describe|overview)\b/, () => 'info'],
];

const FILLER = /^(?:please\s+|could you\s+|can you\s+|would you\s+|i want to\s+|i'd like to\s+|i would like to\s+|let me\s+|let's\s+|lets\s+|now\s+|also\s+|just\s+)+/;

function clean(s) {
  return String(s || '').toLowerCase().replace(/colour/g, 'color').replace(/[“”"]/g, '').replace(/\s+/g, ' ')
    .replace(/[.!?]+$/, '').trim().replace(FILLER, '').replace(/\s+please$/, '').trim();
}

function one(clause, opts) {
  const c = clean(clause);
  if (!c) return null;
  for (const [re, fn] of RULES) {
    const m = re.exec(c);
    if (!m) continue;
    const r = fn(m, c, opts);
    if (r) return [].concat(r);
  }
  return null;
}

export function translatePlain(text, opts = {}) {
  const src = String(text || '').trim();
  if (!src || src.length > 300) return null;
  // several requests in one line: "colour by flexibility and show the 2D structure", "...; then play"
  const parts = src.split(/\s*(?:;|\.\s+|,?\s+and then\s+|,?\s+then\s+|,\s*and\s+|\s+and\s+(?=(?:please\s+)?(?:show|hide|colou?r|select|label|zoom|play|go|open|save|turn|rotate|spin|make|compare|reset)\b))\s*/i).filter(Boolean);
  if (parts.length >= 2) {
    const out = [];
    for (const p of parts) { const r = one(p, opts); if (!r) { out.length = 0; break; } out.push(...r); }
    if (out.length) return { commands: out };
  }
  const whole = one(src, opts);
  return whole ? { commands: whole } : null;
}

// Examples shown by "help plain" and in the console's first line.
export const PLAIN_EXAMPLES = [
  'colour by flexibility', 'show the switch at 25 degrees with 150 mM salt', 'compare models 1 and 3',
  'select residues 10 to 20', 'make the guanines orange', 'show sticks for residues 22 to 25', 'label every 10th residue',
  'how does conformation 1 turn into 3', 'zoom in on helix H2', 'save a picture', 'make a video of the transition from 1 to 4',
];
