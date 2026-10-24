// RCV selection language (residue level), written for RNA conformational ensembles.
//
//   all | none | sele                    everything / nothing / the current interactive selection
//   resi 5  resi 5-12  resi 5-12+20+31-33   author residue numbers (as in the structure files)
//   pos 1-10                             sequence positions (1 = first residue shown)
//   resn A+G   A   C   G   U              nucleotide type
//   chain A                              chain identifier
//   paired | unpaired | pk               base-paired / unpaired / pseudoknotted in the shown conformation
//   helix H2 | helix 2                   residues of a helix family
//   switching                            pairing partner differs between conformations of the ensemble
//   changed                              pairing differs between the previous and the shown conformation
//   rmsf > 3   rmsf < 1.5                per-residue flexibility across the ensemble (Å)
//   moved > 2                            displacement between the previous and the shown conformation (Å)
//   within 6 of <sel>                    residues with an atom within 6 Å of the selection (shown coordinates)
//   name P                               residues that have this atom (atom-level commands use the atom)
//   first | last                         first / last residue
//   <name>                               a named selection
//   12-18                                shorthand for pos 12-18
// combined with  and / or / not  (also & | !) and parentheses.

const KEYWORDS = ['all', 'none', 'sele', 'resi', 'pos', 'resn', 'chain', 'paired', 'unpaired', 'pk', 'helix', 'switching', 'changed',
  'rmsf', 'moved', 'within', 'of', 'name', 'first', 'last', 'and', 'or', 'not', 'A', 'C', 'G', 'U'];
export const SELECTION_KEYWORDS = KEYWORDS;

function tokenize(src) {
  return [...String(src).matchAll(/\(|\)|&|\||!|>=|<=|>|<|[^\s()&|!<>]+/g)].map(m => m[0]);
}

// Parse "5-12+20+31-33" into ranges [[5,12],[20,20],[31,33]] (negative numbers allowed: -3-5).
function ranges(tok) {
  if (!tok) throw new Error('expected a number or range');
  return tok.split('+').map(part => {
    const m = /^(-?\d+)(?:-(-?\d+))?$/.exec(part);
    if (!m) throw new Error(`"${part}" is not a number or range`);
    const a = +m[1], b = m[2] !== undefined ? +m[2] : a;
    return [Math.min(a, b), Math.max(a, b)];
  });
}

export class SelectionError extends Error {}

// ctx: { ens, conf, prev, selection:Set, named:Map<string, Set>, coords:Float32Array }
export function evaluateSelection(src, ctx) {
  const text = String(src || '').trim();
  if (!text) throw new SelectionError('empty selection');
  const toks = tokenize(text);
  let i = 0;
  const ens = ctx.ens, n = ens.nRes;
  const peek = () => toks[i], next = () => toks[i++];
  const all = () => new Set([...Array(n).keys()]);
  const where = fn => { const s = new Set(); for (let r = 0; r < n; r++) if (fn(r)) s.add(r); return s; };
  const expectNum = () => { const t = next(); const v = parseFloat(t); if (!Number.isFinite(v)) throw new SelectionError(`expected a number, got "${t ?? 'end'}"`); return v; };
  const cmp = () => {
    const op = next();
    if (!['>', '<', '>=', '<='].includes(op)) throw new SelectionError(`expected > or < after the keyword, got "${op ?? 'end'}"`);
    const v = expectNum();
    return x => (op === '>' ? x > v : op === '<' ? x < v : op === '>=' ? x >= v : x <= v);
  };
  const partner = k => {
    const p = new Int32Array(n).fill(-1);
    for (const [a, b] of ens.conformations[k].pairs) { p[a] = b; p[b] = a; }
    return p;
  };

  const parseOr = () => {
    let s = parseAnd();
    while (peek() === 'or' || peek() === '|') { next(); const t = parseAnd(); s = new Set([...s, ...t]); }
    return s;
  };
  const parseAnd = () => {
    let s = parseNot();
    for (;;) {
      if (peek() === 'and' || peek() === '&') { next(); const t = parseNot(); s = new Set([...s].filter(r => t.has(r))); continue; }
      // implicit "and" is not allowed: it hides typos
      return s;
    }
  };
  const parseNot = () => {
    if (peek() === 'not' || peek() === '!') { next(); const s = parseNot(); return where(r => !s.has(r)); }
    return parseAtom();
  };
  const parseAtom = () => {
    const t = next();
    if (t === undefined) throw new SelectionError('selection ends too early');
    if (t === '(') { const s = parseOr(); if (next() !== ')') throw new SelectionError('missing ")"'); return s; }
    const low = t.toLowerCase();
    switch (low) {
      case 'all': case '*': return all();
      case 'none': return new Set();
      case 'sele': case 'selected': case 'selection': return new Set(ctx.selection || []);
      case 'resi': {
        const rs = ranges(next());
        return where(r => { const v = parseInt(ens.residues[r].num, 10); return rs.some(([a, b]) => v >= a && v <= b); });
      }
      case 'pos': { const rs = ranges(next()); return where(r => rs.some(([a, b]) => r + 1 >= a && r + 1 <= b)); }
      case 'resn': { const set = new Set(next().toUpperCase().split('+')); return where(r => set.has(ens.residues[r].name) || set.has((ens.residues[r].comp || '').toUpperCase())); }
      case 'chain': { const c = next(); return where(r => (ens.residues[r].chain || '') === c); }
      case 'paired': case 'stem': { const p = partner(ctx.conf); return where(r => p[r] >= 0); }
      case 'unpaired': case 'loop': { const p = partner(ctx.conf); return where(r => p[r] < 0); }
      case 'pk': case 'pseudoknot': { const s = new Set(); for (const [a, b, , pk] of ens.conformations[ctx.conf].pairs) if (pk) { s.add(a); s.add(b); } return s; }
      case 'helix': {
        const h = next();
        const fams = ens.doc.helix_families || [];
        const f = fams.find(x => x.name.toLowerCase() === String(h).toLowerCase()) || fams[parseInt(String(h).replace(/^h/i, ''), 10) - 1];
        if (!f) throw new SelectionError(`no helix "${h}" (helices: ${fams.map(x => x.name).join(', ') || 'none'})`);
        const hres = ens.conformations[ctx.conf].helix_of_residue || [];
        const inConf = where(r => hres[r] === f.id);
        if (inConf.size) return inConf;
        const s = new Set();            // helix absent from this conformation: use the family's span
        for (let r = f.span[0]; r <= f.span[1]; r++) s.add(r);
        return s;
      }
      case 'switching': {
        const parts = ens.conformations.map((_, k) => partner(k));
        return where(r => parts.some(p => p[r] !== parts[0][r]));
      }
      case 'changed': {
        if (ctx.prev === null || ctx.prev === undefined || ctx.prev === ctx.conf) return new Set();
        const a = partner(ctx.prev), b = partner(ctx.conf);
        return where(r => a[r] !== b[r]);
      }
      case 'rmsf': { const f = cmp(); return where(r => f(ens.rmsf[r])); }
      case 'moved': {
        const f = cmp();
        if (ctx.prev === null || ctx.prev === undefined) return new Set();
        const d = ens.residueDisplacement(ens.coords[ctx.conf], ens.coords[ctx.prev]);
        return where(r => f(d[r]));
      }
      case 'within': {
        const d = expectNum();
        if ((next() || '').toLowerCase() !== 'of') throw new SelectionError('expected "within <Å> of <selection>"');
        const ref = parseAtom();
        const c = ctx.coords || ens.coords[ctx.conf];
        const refAtoms = [];
        for (const r of ref) { const res = ens.residues[r]; for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) refAtoms.push(a); }
        const d2 = d * d;
        return where(r => {
          if (ref.has(r)) return true;
          const res = ens.residues[r];
          for (let a = res.first_atom; a < res.first_atom + res.n_atoms; a++) {
            for (const b of refAtoms) {
              const dx = c[3 * a] - c[3 * b], dy = c[3 * a + 1] - c[3 * b + 1], dz = c[3 * a + 2] - c[3 * b + 2];
              if (dx * dx + dy * dy + dz * dz <= d2) return true;
            }
          }
          return false;
        });
      }
      case 'name': { const names = new Set(next().split('+')); return where(r => [...names].some(nm => nm in ens.atomIndex[r])); }
      case 'first': return new Set([0]);
      case 'last': return new Set([n - 1]);
      default: break;
    }
    if (/^[ACGU]$/i.test(t)) { const b = t.toUpperCase(); return where(r => ens.residues[r].name === b); }
    if (ctx.named && ctx.named.has(t)) return new Set(ctx.named.get(t));
    if (/^-?\d+(-(-?\d+))?(\+-?\d+(-(-?\d+))?)*$/.test(t)) { const rs = ranges(t); return where(r => rs.some(([a, b]) => r + 1 >= a && r + 1 <= b)); }
    throw new SelectionError(`unknown word "${t}" (try: resi 5-10, resn G, paired, helix H1, rmsf > 2, within 5 of resi 12, or a selection name)`);
  };

  const s = parseOr();
  if (i < toks.length) throw new SelectionError(`unexpected "${toks[i]}" — combine conditions with and / or`);
  return s;
}

// Atom specifier for measurements: "23/N1" (author number / atom), "pos 5/P", or a selection (its C4' or
// first atom when it contains one residue). Returns an atom index.
export function resolveAtom(spec, ctx) {
  const ens = ctx.ens;
  const s = String(spec || '').trim();
  const m = /^(.*?)\s*\/\s*([A-Za-z0-9'"*]+)$/.exec(s);
  const atomName = m ? m[2].replace(/\*/g, "'").replace(/"/g, "''") : null;
  const resSel = m ? m[1] : s;
  let residues;
  if (/^-?\d+$/.test(resSel.trim())) {
    const v = parseInt(resSel, 10);
    residues = new Set(ens.residues.map((r, i) => (parseInt(r.num, 10) === v ? i : -1)).filter(i => i >= 0));
  } else residues = evaluateSelection(resSel, ctx);
  if (residues.size !== 1) throw new SelectionError(`"${spec}" must name exactly one residue (it matches ${residues.size})`);
  const r = [...residues][0];
  const idx = ens.atomIndex[r];
  if (atomName) {
    if (!(atomName in idx)) throw new SelectionError(`residue ${ens.resLabel(r)} has no atom ${atomName} (atoms: ${Object.keys(idx).join(' ')})`);
    return idx[atomName];
  }
  return ens.trace[r] >= 0 ? ens.trace[r] : ens.residues[r].first_atom;
}

// Compact text for a residue set: "G5-A12, U20" style using author numbers.
export function describeResidues(ens, set) {
  const list = [...set].sort((a, b) => a - b);
  if (!list.length) return 'nothing';
  const parts = [];
  for (let i = 0; i < list.length;) {
    let j = i;
    while (j + 1 < list.length && list[j + 1] === list[j] + 1) j++;
    parts.push(i === j ? ens.resLabel(list[i]) : `${ens.resLabel(list[i])}–${ens.resLabel(list[j])}`);
    i = j + 1;
  }
  return parts.length > 8 ? `${parts.slice(0, 8).join(', ')} … (${list.length} nt)` : `${parts.join(', ')} (${list.length} nt)`;
}
