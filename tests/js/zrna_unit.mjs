// Z-RNA detector of the viewer (src/rcv/static/js/zrna.js), run by tests/test_zrna.py when Node is available:
//   node tests/js/zrna_unit.mjs <cases.json>   cases: [{ name, seq, residues: [{atom: [x,y,z]}], pairs: [[i, j]] }]
// Prints one JSON line: { name: [{ pairs, residues, twists, synPurines, antiPyrimidines }] } for comparison with
// rnafold3d/zrna.py.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static/js');
globalThis.localStorage = { getItem: () => null, setItem: () => {} };
globalThis.document = { createElementNS: () => ({}), createElement: () => ({}), documentElement: { dataset: {} } };
const { detectZ } = await import(path.join(STATIC, 'zrna.js'));

const cases = JSON.parse(readFileSync(process.argv[2], 'utf8'));
const out = {};
for (const c of cases) {
  const atoms = (r, name) => c.residues[r]?.[name] ?? null;
  out[c.name] = detectZ(c.seq, atoms, c.pairs).map(z => ({
    pairs: z.pairs.map(([i, j]) => [i + 1, j + 1]), residues: z.residues.map(r => r + 1),
    twists: z.twists.map(t => Math.round(t * 10) / 10), synPurines: Math.round(z.synPurines * 100) / 100,
    antiPyrimidines: Math.round(z.antiPyrimidines * 100) / 100,
  }));
}
console.log(JSON.stringify(out));
