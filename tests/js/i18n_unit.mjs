// Unit tests for the interface-language module (run by tests/test_rcv06_lang.py when Node >= 18 is available).
//   node tests/js/i18n_unit.mjs
// The module is loaded without a browser: only the translation functions are used, not the DOM walker.
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const STATIC = path.resolve(here, '../../src/rcv/static');
let failures = 0;
const ok = (cond, msg) => { if (!cond) { failures++; console.error('FAIL', msg); } else console.log('ok', msg); };

let current = null;
globalThis.localStorage = { getItem: () => null, setItem: () => {}, removeItem: () => {} };
globalThis.location = { search: '' };
Object.defineProperty(globalThis, 'navigator', { value: { languages: ['de-AT', 'en'] }, configurable: true });
globalThis.document = { documentElement: {}, body: null };
globalThis.MutationObserver = class { observe() {} };
globalThis.fetch = async url => {
  const code = String(url).match(/i18n\/(.+)\.json/)[1];
  current = JSON.parse(readFileSync(path.join(STATIC, 'i18n', `${code}.json`), 'utf8'));
  return { ok: true, json: async () => current };
};
const i18n = await import(path.join(STATIC, 'js', 'i18n.js'));

// system language matching
ok(i18n.matchLanguage(['de-AT']) === 'de', 'de-AT -> de');
ok(i18n.matchLanguage(['zh-HK']) === 'zh-TW' && i18n.matchLanguage(['zh-Hans-CN']) === 'zh-CN', 'Chinese variants');
ok(i18n.matchLanguage(['pt-PT']) === 'pt-BR' && i18n.matchLanguage(['nl-NL']) === 'en', 'Portuguese and unsupported');
ok(i18n.preferredLanguage() === 'de', 'system language is used when nothing was chosen');

// German
await i18n.setLanguage('de', { save: false });
const { t } = i18n;
ok(t('View') === 'Ansicht', 'plain menu label');
ok(t('  View ') === '  Ansicht ', 'surrounding spaces are kept');
ok(t('Undo (Ctrl+Z)') === 'Rückgängig (Strg+Z)', 'own translation of a label with its key');
ok(t('View (F9)') === 'Ansicht (F9)', 'keyboard hint after a known label');
ok(t('20 conformations') === '20 Konformationen', 'pattern with a number');
ok(t('Something new that has no translation') === 'Something new that has no translation', 'unknown text unchanged');
ok(t('GGCAGAUCUGAGCCUGGGAGCUCUCUGCC') === 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC', 'sequences unchanged');
ok(t('Deuteranopia: closest pair #2a78d6 / #6250d6 — ΔE 9.5') === 'Deuteranopie: ähnlichstes Paar #2a78d6 / #6250d6 — ΔE 9.5',
  'colour check, distinguishable');
ok(t('Normal vision: closest pair #eb6834 / #e34948 — ΔE 7.1 (hard to tell apart)')
  === 'Normales Sehen: ähnlichstes Paar #eb6834 / #e34948 — ΔE 7.1 (schwer zu unterscheiden)', 'colour check, too similar');
ok(t('Model 1: 0 other conformation(s) within 2 Å; nearest is Model 20 at 3.19 Å.').startsWith('Model 1: 0 weitere'),
  'filled-in names are kept');
ok(t('Conformations · Evaluator energy (kcal/mol)').startsWith('Konformationen · '), 'one value with a unit');
ok(i18n.tf('{0} conformations', 12) === '12 Konformationen', 'tf fills a pattern');

// Every language: menu names are translated; English comes back unchanged.
for (const l of i18n.LANGUAGES) {
  await i18n.setLanguage(l.code, { save: false });
  const v = t('View'), e = t('Edit');
  if (l.code === 'en') ok(v === 'View' && e === 'Edit', 'English is the source text');
  else ok(v !== 'View' && e !== 'Edit' && i18n.currentLanguage() === l.code, `${l.code}: menus translated (${e}/${v})`);
}

if (failures) { console.error(`${failures} failure(s)`); process.exit(1); }
console.log('all i18n tests passed');
