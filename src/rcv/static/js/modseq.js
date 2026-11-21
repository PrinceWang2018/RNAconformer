// Modified nucleotides in a sequence line (mirrors rnafold3d/modifications.py): [code] in brackets for a modified
// nucleotide (MODOMICS short name, alias or PDB component: GGC[m6A]CU[Y]A, [A>I] for an A-to-I edited site), and * between
// two nucleotides for a phosphorothioate linkage (A*G). The catalogue comes from the server (/api/modifications).

let catalogue = null;            // { categories, entries }
let index = null;                // { exact: Map, lower: Map(lower -> Set(code)), byCode: Map }
let pending = null;

export function setCatalogue(cat) {
  catalogue = cat;
  const exact = new Map(), lower = new Map(), byCode = new Map();
  for (const e of cat?.entries || []) {
    byCode.set(e.code, e);
    for (const k of [e.code, ...(e.aliases || []), ...(e.ccd ? [e.ccd] : []), ...(e.ccd_also || [])]) {
      if (!exact.has(k)) exact.set(k, e.code);
      const lk = k.toLowerCase();
      if (!lower.has(lk)) lower.set(lk, new Set());
      lower.get(lk).add(e.code);
    }
  }
  index = { exact, lower, byCode };
  return index;
}

export function getCatalogue() { return catalogue; }
export function entry(code) { return index?.byCode.get(code) || null; }

// Fetch once per session (null when the engine is not reachable).
export function loadCatalogue(api) {
  if (catalogue) return Promise.resolve(catalogue);
  if (!pending) pending = api.get('/api/modifications').then(c => setCatalogue(c) && c).catch(() => { pending = null; return null; });
  return pending;
}

export function lookup(token) {
  if (!index) return null;
  const t = String(token).trim().replace(/ψ/g, 'Ψ');
  if (index.exact.has(t)) return index.byCode.get(index.exact.get(t));
  const hits = index.lower.get(t.toLowerCase());
  if (hits && hits.size === 1) return index.byCode.get([...hits][0]);
  const t2 = t.replace(/[ _]/g, '');
  return t2 !== t ? lookup(t2) : null;
}

export const hasNotation = text => /[[*]/.test(String(text || ''));

// One sequence line, upper case and T -> U outside the brackets. Digits outside the brackets are dropped: sequences
// copied from GenBank or a paper carry position numbers ("1 ggcagatctg 11 agcctgggag").
export function cleanLine(line) {
  return String(line || '').replace(/\s+/g, '').split(/(\[[^[\]]*\])/).map(p => (p.startsWith('[') ? p : p.replace(/\d+/g, '').toUpperCase().replace(/T/g, 'U'))).join('');
}

// -> { seq (parent bases; N for an unknown name or letter), mods [{pos, code}], problems [text], bad [letters],
//      unknown (catalogue not loaded) }
export function parse(text) {
  const s = String(text || '').replace(/\s+/g, '');
  const seq = [], mods = [], problems = [], bad = new Set();
  let ps = false, unknown = false;
  for (let k = 0; k < s.length;) {
    const ch = s[k];
    if (ch === '[') {
      const end = s.indexOf(']', k);
      if (end < 0) { problems.push(`unclosed "[" at character ${k + 1}`); break; }
      const tok = s.slice(k + 1, end);
      const e = lookup(tok);
      if (!index) { unknown = true; seq.push('N'); mods.push({ pos: seq.length, code: tok }); }
      else if (!e) { problems.push(`unknown modification [${tok}]`); seq.push('N'); }
      else if (e.kind === 'backbone') problems.push(`[${tok}] is a linkage: write * between two nucleotides`);
      else {
        seq.push(e.build);
        mods.push({ pos: seq.length, code: e.code });
        if (ps) mods.push({ pos: seq.length, code: 'ps' });
      }
      ps = false;
      k = end + 1;
      continue;
    }
    if (ch === '*') {
      if (!seq.length) problems.push('"*" (phosphorothioate) must stand between two nucleotides');
      ps = true; k++; continue;
    }
    if (ch >= '0' && ch <= '9') { k++; continue; }           // position numbers
    let up = ch.toUpperCase();
    if (up === 'T') up = 'U';
    if ('ACGU'.includes(up)) {
      seq.push(up);
      if (ps) mods.push({ pos: seq.length, code: 'ps' });
    } else { seq.push('N'); bad.add(ch); }
    ps = false; k++;
  }
  if (ps) problems.push('"*" at the end of the sequence');
  return { seq: seq.join(''), mods, problems, bad: [...bad], unknown };
}

export function format(seq, mods) {
  const by = new Map(), ps = new Set();
  for (const m of mods) (m.code === 'ps' ? ps.add(m.pos) : by.set(m.pos, m.code));
  let out = '';
  for (let i = 1; i <= seq.length; i++) {
    if (ps.has(i) && i > 1) out += '*';
    out += by.has(i) ? `[${by.get(i)}]` : seq[i - 1];
  }
  return out;
}

// "m6A ×1, Ψ ×2" (codes in order of first appearance; pseudouridine shown as Ψ).
export function summary(mods) {
  const count = new Map();
  for (const m of mods) count.set(m.code, (count.get(m.code) || 0) + 1);
  return [...count].map(([c, n]) => `${display(c)} ×${n}`).join(', ');
}

export const display = code => ({ Y: 'Ψ', m1Y: 'm1Ψ', m3Y: 'm3Ψ', Ym: 'Ψm', ps: 'PS' }[code] || code);

// Central A of the m6A consensus DRACH (D = A/G/U, R = A/G, H = A/C/U), 0-based positions.
export function drachSites(seq) {
  const out = [];
  for (let i = 2; i + 2 < seq.length; i++) {
    if (seq[i] === 'A' && 'AGU'.includes(seq[i - 2]) && 'AG'.includes(seq[i - 1]) && seq[i + 1] === 'C' && 'ACU'.includes(seq[i + 2])) out.push(i);
  }
  return out;
}
