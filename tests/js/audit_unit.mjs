// Start page and job watcher after the 0.7.6 audit (run by tests/test_v076_audit.py when Node is available):
// numbered sequences, the dot-bracket check, and the job watcher when requests fail.
import { fileURLToPath } from 'node:url';
import path from 'node:path';
import assert from 'node:assert/strict';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
globalThis.localStorage = { getItem: () => null, setItem: () => {} };
const node = () => ({ style: {}, dataset: {}, classList: { add() {}, remove() {}, toggle() {} }, append() {}, setAttribute() {}, addEventListener() {}, appendChild(x) { return x; } });
globalThis.document = { createElementNS: node, createElement: node, createTextNode: t => t, documentElement: { dataset: {} }, querySelector: () => null, addEventListener() {}, body: node() };
globalThis.window = globalThis.window || { addEventListener() {}, matchMedia: () => ({ matches: false, addEventListener() {} }) };
try { globalThis.navigator = globalThis.navigator || { languages: ['en'] }; } catch { /* read-only in newer Node */ }

const modseq = await import(path.join(STATIC, 'modseq.js'));
const { parseFastaRecords, parseSequenceInput, dotBracketProblem } = await import(path.join(STATIC, 'welcome.js'));
const { watchJob } = await import(path.join(STATIC, 'jobs.js'));

// ---- C1: position numbers are ignored, digits inside [brackets] stay
assert.equal(modseq.cleanLine('1 ggcagatctg 11 agcctgggag'), 'GGCAGAUCUGAGCCUGGGAG');
assert.equal(modseq.cleanLine('GGC[m2,2G]C 6 U'), 'GGC[m2,2G]CU');
let p = parseSequenceInput('        1 ggcagatctg agcctgggag ctctctgcc\n//\n');
assert.equal(p.seq, 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC');
assert.deepEqual(p.bad, []);
assert.deepEqual(modseq.parse('GG1CA').bad, []);
assert.deepEqual(parseSequenceInput('GGCAXGGC').bad, ['X']);
const recs = parseFastaRecords('>a\n1 ggcacuucgg 11 ugcc\n//\n>b\nGGCACUUCGGUGCC\n(((((....)))))\n');
assert.equal(recs.length, 2);
assert.equal(recs[0].sequence, 'GGCACUUCGGUGCC');
assert.equal(recs[0].problem, null);
assert.equal(recs[1].problem, null);

// ---- C2: the dot-bracket is checked before the run
assert.equal(dotBracketProblem('((((....))))', 12), null);
assert.equal(dotBracketProblem('((..[[..))..]]', 14), null);
assert.match(dotBracketProblem('((((....))))', 14), /12 characters/);
assert.match(dotBracketProblem('((....)))', 9), /"\)" at position 9 has no opening/);
assert.match(dotBracketProblem('((((....)))', 11), /"\(" at position 1 is never closed/);
assert.match(dotBracketProblem('((x(....))))', 12), /"x" at position 3/);
assert.match(dotBracketProblem('(((..)))', 8), /fewer than 3 nt/);
assert.match(parseFastaRecords('>a\nGGCACUUCGGUGCC\n>b\nGGCACUUCGGUGCC\n(((((....))).)\n')[1].problem, /never closed/);

// ---- B1: one failed request changes nothing; several say "lost"; the job is followed again afterwards
const sleep = ms => new Promise(r => setTimeout(r, ms));
{
  const script = ['ok', 'fail', 'ok', 'fail', 'fail', 'fail', 'ok', 'done'];
  let k = 0;
  const api = { get: async () => {
    const s = script[Math.min(k++, script.length - 1)];
    if (s === 'fail') throw new TypeError('Failed to fetch');
    return { id: 'j', status: s === 'done' ? 'done' : 'running', progress: k / 10, stages: [{ id: 's', state: 'active' }] };
  } };
  const updates = [];
  let done = null;
  watchJob(api, 'j', { interval: 5, lostAfter: 3, onUpdate: j => updates.push(j), onDone: j => { done = j; } });
  for (let i = 0; i < 400 && !done; i++) await sleep(10);
  assert.ok(done, 'the job was not followed to its end');
  assert.ok(updates.every(u => (u.stages || []).length === 1), 'a failed request emptied the progress card');
  assert.equal(updates.filter(u => u.lost).length, 1);                // after the third failure in a row only
  assert.equal(updates.at(-1).status, 'done');
}
{
  let n = 0;
  const api = { get: async () => { if (n++ === 0) return { id: 'j', status: 'running', title: 't' }; throw Object.assign(new Error('unknown job'), { status: 404 }); } };
  let failed = null;
  watchJob(api, 'j', { interval: 5, onFail: j => { failed = j; } });
  for (let i = 0; i < 100 && !failed; i++) await sleep(10);
  assert.equal(failed?.status, 'failed');
  assert.match(failed.error, /restarted/);
  assert.equal(failed.title, 't');
}

// ---- export of the secondary structure: dot-bracket, CT, BPSEQ (the same text as rcv/ss2dexport.py)
{
  const sec = await import(path.join(STATIC, 'export/secondary.js'));
  assert.deepEqual(sec.partners('((..[[..))..]]'), [9, 8, -1, -1, 13, 12, -1, -1, 1, 0, -1, -1, 5, 4]);
  const recs = [{ name: 'tar  demo', id: 'conf_001', seq: 'GGCACUUCGGUGCC', dbn: '(((((....)))))', energy: -7.1, also: ['conf_002'], numbers: null },
    { name: 'tar  demo', id: 'conf_003', seq: 'GGCACUUCGGUGCC', dbn: '.((((....)))).', energy: null, also: [], numbers: null }];
  assert.equal(sec.toDotBracket(recs), '>tar demo conf_001 dG=-7.10 kcal/mol also conf_002\nGGCACUUCGGUGCC\n(((((....)))))\n>tar demo conf_003\nGGCACUUCGGUGCC\n.((((....)))).\n');
  const ct = sec.toCT(recs).split('\n');
  assert.equal(ct[0], '   14  ENERGY = -7.10  tar demo conf_001');
  assert.equal(ct[1], '    1 G       0    2   14    1');
  assert.equal(ct[14], '   14 C      13    0    1   14');
  assert.equal(ct[15], '   14  tar demo conf_003');
  const bp = sec.toBpseq(recs).split('\n');
  assert.equal(bp[0], '# tar demo conf_001 dG=-7.10 kcal/mol');
  assert.equal(bp[1], '1 G 14');
  assert.equal(bp[6], '6 U 0');
  if (process.argv[2]) {                       // records written by the Python side: the text must be identical
    const { readFileSync } = await import('node:fs');
    const cases = JSON.parse(readFileSync(process.argv[2], 'utf8'));
    for (const c of cases.export) assert.equal(sec.FORMATS[c.format].write(c.records), c.text, `export differs from rcv ss (${c.format})`);
    // ---- C5: what is bound to the RNA is read with the document and moves with its conformation
    const { Ensemble } = await import(path.join(STATIC, 'ensemble.js'));
    const ens = new Ensemble(cases.complex), ctx = ens.context;
    assert.ok(ctx, 'the document of a complex has a context');
    assert.deepEqual([ctx.nAtoms, ctx.nPoints, ctx.bonds.length, ctx.segs.length, ctx.coords.length], [7, 8, 12, 7, 3]);
    assert.ok(ctx.radius > 5 && ctx.radius < 60);
    const dists = () => ens.coords.map((x, k) => { const a = ens.trace[4], c = ctx.coords[k]; return Math.hypot(c[18] - x[3 * a], c[19] - x[3 * a + 1], c[20] - x[3 * a + 2]); });
    const before = dists();
    assert.ok(Math.abs(before[0] - Math.hypot(5, 4)) < 0.01, 'the ion sits where the file puts it');
    ens.superposeOn([0, 1, 2, 3], 1);
    dists().forEach((v, k) => assert.ok(Math.abs(v - before[k]) < 1e-3, 'superposing moves what is bound with its conformation'));
    const plain = new Ensemble({ ...cases.complex, context: undefined });
    assert.equal(plain.context, null);
    assert.equal(new Ensemble({ ...cases.complex, context: { ...cases.complex.context, coords: cases.complex.context.coords.slice(1) } }).context, null);
  }
  const ens = { name: 'x', sequence: 'GGCACUUCGGUGCC', residues: [...'GGCACUUCGGUGCC'].map((b, i) => ({ name: b, num: String(i + 5) })),
    conformations: [{ dbn: '(((((....)))))' }, { dbn: '.((((....)))).' }, { dbn: '(((((....)))))' }],
    shortId: k => `c${k + 1}`, value: (k, key) => (key === 'dg2d' ? -k : null) };
  const d = sec.recordsOf(ens, 'distinct', [2, 0, 1], 0);
  assert.deepEqual(d.map(r => [r.id, r.also]), [['c3', ['c1']], ['c2', []]]);
  assert.equal(sec.recordsOf(ens, 'current', null, 1)[0].id, 'c2');
  assert.equal(sec.toCT(sec.recordsOf(ens, 'current', null, 0)).split('\n')[1], '    1 G       0    2   14    5');
}
console.log('all audit tests passed');
