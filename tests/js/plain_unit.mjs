// Unit tests for the plain-language console (static/js/console/plain.js); run by tests/test_v075_workflow.py.
//   node tests/js/plain_unit.mjs
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const JS = path.resolve(here, '../../src/rcv/static/js');
const { translatePlain, residuePhrase, PLAIN_EXAMPLES } = await import(path.join(JS, 'console', 'plain.js'));
let failures = 0;
const SCHEMES = ['nucleotide', 'helix', 'elements', 'pairing', 'pairfreq', 'switching', 'flexibility', 'displacement', 'position',
  'nucleotide-classic', 'uniform'].map(id => ({ id, label: id }));
const opts = { schemes: SCHEMES, ligands: { theophylline: {} } };
const tr = s => translatePlain(s, opts)?.commands.join(' ; ') ?? null;
const eq = (input, want) => {
  const got = tr(input);
  if (got !== want) { failures++; console.error(`FAIL ${JSON.stringify(input)}\n   got  ${got}\n   want ${want}`); }
  else console.log('ok', input, '->', got);
};

eq('colour by flexibility', 'color flexibility');
eq('Color it by the switching residues', 'color switching');
eq('colour by helices', 'color helix');
eq('show the switch at 25 degrees with 150 mM salt', 'conditions 25, 150 mM ; tab switch');
eq('set the temperature to 30 °C', 'conditions 30');
eq('with theophylline at 2 mM', 'conditions theophylline 2 mM');
eq('add theophylline 10 µM at 25 C', 'conditions 25, theophylline 10 µM');
eq('without ligand', 'conditions no ligand');
eq('physiological conditions', 'conditions 37, 150 mM');
eq('standard conditions', 'conditions standard');
eq('compare models 1 and 3', 'align 1, 3');
eq('what is the RMSD between conformation 2 and 5?', 'align 2, 5');
eq('select residues 10 to 20', 'select resi 10-20');
eq('select nucleotides 5, 9 and 12', 'select resi 5+9+12');
eq('select the guanines', 'select resn G');
eq('highlight the switching residues', 'select switching');
eq('make the guanines orange', 'color orange, resn G');
eq('colour residues 10-20 red', 'color red, resi 10-20');
eq('make everything blue', 'color blue, all');
eq('reset the colours', 'color default, all');
eq('show sticks for residues 22 to 25', 'show sticks, resi 22-25');
eq('hide the 2D structure', 'hide 2d');
eq('show all conformations', 'show ensemble');
eq('show the flexibility putty', 'show putty');
eq('label every 10th residue', 'set labels, every10');
eq('remove the labels', 'hide labels');
eq('label residue 34', 'label resi 34');
eq('how does conformation 1 turn into 3', 'pathway 1, 3');
eq('show the pathway from 2 to 4', 'pathway 2, 4');
eq('zoom in on helix H2', 'zoom helix H2');
eq('zoom on residues 5 to 9', 'zoom resi 5-9');
eq('reset the view', 'reset');
eq('rotate 90 degrees around x', 'turn x, 90');
eq('turn the molecule by 45', 'turn y, 45');
eq('stop rotating', 'spin off');
eq('start spinning', 'spin on');
eq('go to model 5', 'conf 5');
eq('next conformation please', 'next');
eq('play the tour', 'play');
eq('stop', 'stop');
eq('morph from 1 to 4', 'morph 1, 4');
eq('distance between 12/n1 and 30/n3', 'distance 12/N1, 30/N3');
eq('how far is residue 5 from residue 20', 'distance 5, 20');
eq('save a picture', 'png');
eq('export all models as PDB', 'save pdb all');
eq('save the ensemble', 'save rcv');
eq('make a video of the transition from 1 to 4', 'video transition 1, 4');
eq('make a rotation movie', 'video rotate');
eq('fetch 1ANR', 'fetch 1ANR');
eq('download PDB entry 3d2g', 'fetch 3D2G');
eq('dark background', 'bg black');
eq('use spheres', 'rep spheres');
eq('show it as ball and stick', 'rep ballstick');
eq('open the switch tab', 'tab switch');
eq('how many guanines', 'count resn G');
eq('which base pairs does model 2 have', 'pairs 2');
eq('tell me about this ensemble', 'info');
eq('sort by energy', 'order energy');
eq('colour by flexibility and show the 2D structure', 'color flexibility ; show 2d');
eq('select residues 3 to 7; then zoom in on the selection', 'select resi 3-7 ; zoom sele');
eq('please clear the selection', 'deselect');
eq('the weather is nice today', null);
eq('banana', null);

if (residuePhrase('the loops') !== 'unpaired' || residuePhrase('pseudoknotted residues') !== 'pk') { failures++; console.error('FAIL residuePhrase'); }
for (const e of PLAIN_EXAMPLES) if (!tr(e)) { failures++; console.error('FAIL example not understood:', e); }

if (failures) { console.error(`${failures} failure(s)`); process.exit(1); }
console.log('all plain-language tests passed');
