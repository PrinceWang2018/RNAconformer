// The notes of a running prediction in every interface language (0.7.6, audit F5).
// node notes_unit.mjs <notes.json>: notes as rcv/jobs.py wrote them ({ text, keep, lead }; same: true = a line of
// another program, which must come back unchanged). Used by tests/test_v076_audit.py.
import { readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const STATIC = path.join(path.dirname(fileURLToPath(import.meta.url)), '..', '..', 'src', 'rcv', 'static');
globalThis.localStorage = { getItem: () => null, setItem: () => {}, removeItem: () => {} };
globalThis.location = { search: '' };
Object.defineProperty(globalThis, 'navigator', { value: { languages: ['en'] }, configurable: true });
globalThis.document = { documentElement: {}, body: null };
globalThis.MutationObserver = class { observe() {} };
globalThis.fetch = async url => ({ ok: true, json: async () => JSON.parse(readFileSync(path.join(STATIC, 'i18n', `${String(url).match(/i18n\/(.+)\.json/)[1]}.json`), 'utf8')) });
const i18n = await import(path.join(STATIC, 'js', 'i18n.js'));

const notes = JSON.parse(readFileSync(process.argv[2], 'utf8'));
let failures = 0;
const fail = msg => { failures++; console.error('FAIL', msg); };
for (const l of i18n.LANGUAGES) {
  if (l.code === 'en') continue;
  await i18n.setLanguage(l.code, { save: false });
  for (const n of notes) {
    const lead = n.lead && n.text.startsWith(n.lead) ? n.lead : '';
    const body = n.text.slice(lead.length);
    const out = i18n.tKeep(body, n.keep || []);
    if (/[-]/.test(out)) fail(`${l.code}: a marker is left in "${out}"`);
    for (const k of n.keep || []) if (!out.includes(k)) fail(`${l.code}: "${k}" was changed in "${out}"`);
    if (l.code === 'de' && /^ss\d+ \(suboptimal[,)]/.test(body)) continue;       // German: "suboptimal" as well
    if (n.same ? out !== body : out === body) fail(`${l.code}: ${n.same ? 'changed' : 'left in English'}: "${n.text.slice(0, 120)}"${n.same ? ` -> "${out.slice(0, 120)}"` : ''}`);
  }
}
await i18n.setLanguage('de', { save: false });
const de = i18n.tKeep('Prediction of 76 nt started on this computer with 4 CPU cores; results are saved in /home/some user/rcv runs/x.', ['/home/some user/rcv runs/x']);
if (de !== 'Vorhersage für 76 nt auf diesem Computer mit 4 CPU-Kern(en) gestartet; die Ergebnisse werden in /home/some user/rcv runs/x gespeichert.') fail('German: ' + de);
if (failures) { console.error(`${failures} failure(s)`); process.exit(1); }
console.log(`all ${notes.length} notes are translated in ${i18n.LANGUAGES.length - 1} languages`);
