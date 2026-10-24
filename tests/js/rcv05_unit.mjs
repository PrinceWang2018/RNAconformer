// Unit tests for RCV 0.5 browser modules (run by tests/test_rcv05.py when Node >= 18 is available).
//   node tests/js/rcv05_unit.mjs <ensemble-document.json>
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
const imp = p => import(path.join(STATIC, p));
let failures = 0;
const ok = (cond, msg) => { if (!cond) { failures++; console.error('FAIL', msg); } else console.log('ok', msg); };
const throws = (fn, re, msg) => { try { fn(); ok(false, msg + ' (no error)'); } catch (e) { ok(re.test(e.message), `${msg}: ${e.message}`); } };

const { Ensemble } = await imp('ensemble.js');
const S = await imp('console/selection.js');
const L = await imp('ss2dlayout.js');
const C = await imp('colors.js');
const V = await imp('export/video.js');
const K = await imp('console/commands.js');

const doc = JSON.parse(readFileSync(process.argv[2], 'utf8'));
const ens = new Ensemble(doc);
const ctx = { ens, conf: 0, prev: 1, selection: new Set([3, 4]), named: new Map([['acc', new Set([0, 1, 2])]]) };
const sel = q => S.evaluateSelection(q, ctx);

// ---- selection language
ok(sel('all').size === ens.nRes, 'all');
ok(sel('resi 1-7').size === 7 && sel('resi 1-7+66-72').size === 14, 'residue ranges');
ok(sel('pos 10-12 or G and pos 1-20').size === sel('pos 10-12').size + [...sel('G and pos 1-20')].filter(r => r < 9 || r > 11).length, '"and" binds tighter than "or"');
ok([...sel('resn A+G')].every(r => 'AG'.includes(ens.residues[r].name)), 'resn');
ok(sel('paired').size + sel('unpaired').size === ens.nRes, 'paired / unpaired partition');
ok(sel('not paired').size === sel('unpaired').size, 'not');
ok(sel('acc or sele').size === 5, 'named selection and current selection');
ok(sel('within 4 of resi 34').has(33), 'within');
ok(sel('(resi 1-5 or resi 70-76) and not unpaired').size <= 12, 'parentheses');
ok(sel('switching').size > 0 && sel('changed').size > 0, 'ensemble keywords');
ok(sel('rmsf > 0').size === ens.nRes, 'rmsf threshold');
ok(sel('helix H1').size >= 2, 'helix family');
throws(() => sel('resi 1-5 G'), /unexpected/, 'implicit and is rejected');
throws(() => sel('foo'), /unknown word/, 'unknown word');
throws(() => sel('(resi 1'), /missing/, 'unbalanced parenthesis');
throws(() => sel('rmsf 3'), /expected > or </, 'comparison needed');
ok(S.resolveAtom('34/N1', ctx) === ens.atomIndex[33].N1, 'atom specifier residue/atom');
ok(S.resolveAtom("pos 5/C4'", ctx) === ens.atomIndex[4]["C4'"], 'atom specifier with prime');
throws(() => S.resolveAtom('34/XX', ctx), /has no atom/, 'missing atom');
throws(() => S.resolveAtom('resi 1-3', ctx), /exactly one residue/, 'atom needs one residue');
ok(/G1–A5 \(5 nt\)/.test(S.describeResidues(ens, sel('resi 1-5'))), 'describe residues');

// ---- console argument parsing and colours
ok(JSON.stringify(K.splitArgs("d1, 12/C4', 30/N3")) === JSON.stringify(['d1', "12/C4'", '30/N3']), 'split arguments keeps primes');
ok(JSON.stringify(K.splitArgs('resi 34, "anticodon, loop"')) === JSON.stringify(['resi 34', 'anticodon, loop']), 'quoted argument');
ok(K.parseColor('orange') === '#eb6834' && K.parseColor('#abc') === '#aabbcc' && K.parseColor('nope') === null, 'colour names');

// ---- 2D layouts
for (const c of ens.conformations) {
  const P = L.layoutFor(ens, ens.conformations.indexOf(c));
  ok(P.length === 2 * ens.nRes && !P.some(Number.isNaN), `layout of ${c.id}`);
}
const dbn = '(((.(((....((((((.((......)).)))))).((((....))))(((((.......)))))))).)))....';
const st = [], pairs = [];
for (let i = 0; i < dbn.length; i++) { if (dbn[i] === '(') st.push(i); else if (dbn[i] === ')') pairs.push([st.pop(), i, 'CG', false]); }
const F = L.simpleLayout(dbn.length, L.pairTable(dbn.length, pairs));
let bad = 0, minD = Infinity;
for (let i = 0; i + 1 < dbn.length; i++) if (Math.abs(Math.hypot(F[2 * i + 2] - F[2 * i], F[2 * i + 3] - F[2 * i + 1]) - 1) > 1e-6) bad++;
for (let i = 0; i < dbn.length; i++) for (let j = i + 2; j < dbn.length; j++) minD = Math.min(minD, Math.hypot(F[2 * j] - F[2 * i], F[2 * j + 1] - F[2 * i + 1]));
ok(bad === 0 && minD > 0.99, 'fallback layout: unit backbone steps, no overlaps');
ok(L.nestedKey(dbn.length, pairs) === dbn, 'nested key');
const G = L.simpleLayout(dbn.length, L.pairTable(dbn.length, pairs.slice(0, 9)));
const f = L.fit2D(F, G);
const b0 = L.blendLayouts(F, G, 0, f), b1 = L.blendLayouts(F, G, 1, f);
ok(b0.every((v, i) => Math.abs(v - F[i]) < 1e-9) && b1.every((v, i) => Math.abs(v - G[i]) < 1e-9), 'layout blend hits both ends');

// ---- structure elements & ensemble schemes
const el = C.structureElements(dbn.length, pairs).map(x => ({ stem: 'S', hairpin: 'H', internal: 'I', bulge: 'B', multi: 'M', exterior: 'E', pk: 'P' })[x]).join('');
ok(el === 'SSSISSSMMMMSSSSSSISSHHHHHHSSISSSSSSMSSSSHHHHSSSSSSSSSHHHHHHHSSSSSSSSISSSEEEE', 'structure elements of a cloverleaf');
for (const id of ['elements', 'pairfreq', 'switching']) {
  const cols = C.residueColors(ens, id, { conf: 0 });
  ok(cols.length === ens.nRes && cols.every(c => /^#[0-9a-f]{6}$/i.test(c)), `scheme ${id}`);
  ok(!!C.legend(ens, id, { conf: 0 }).title, `legend ${id}`);
}

// ---- transition paths and frame plans
const p = V.transitionPath(ens, 0, ens.n - 1);
ok(p[0] === 0 && p[p.length - 1] === ens.n - 1 && new Set(p).size === p.length, 'transition path endpoints');
const app = { ens, order: [...Array(ens.n).keys()], state: { conf: 0, prev: null } };
const plan = V.framePlan(app, { content: 'transition', from: 0, to: ens.n - 1, route: 'path', fps: 10, transition: 2, hold: 0.5, boomerang: true, spinDeg: 20 });
ok(plan[0].from === 0 && plan[0].t === 0, 'transition starts on A');
ok(plan.some(fr => fr.phase === 'end' && fr.from === ens.n - 1), 'transition holds on B');
ok(plan[plan.length - 1].from === 0, 'boomerang returns to A');
ok(Math.abs(plan[0].angle + plan[plan.length - 1].angle) < 1e-9, 'slow turn is symmetric');

console.log(failures ? `${failures} failure(s)` : 'all RCV 0.5 JS checks passed');
process.exit(failures ? 1 : 0);
