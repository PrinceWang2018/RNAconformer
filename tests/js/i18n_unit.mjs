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

// Lines of a running prediction are put together from parts (0.7.6, audit F2): a stage, its detail and the time left
await i18n.setLanguage('de', { save: false });
{
  const t2 = i18n.t;
  ok(t2('3D sampling · 12 of 32 trajectories') === '3D-Stichprobe · 12 von 32 Trajektorien', 'stage and detail (was half translated: the stage stayed English)');
  ok(t2('Relaxation and validation · about 40 s left') === 'Relaxation und Validierung · noch etwa 40 s', 'stage and time left');
  ok(t2('Pull-guided sampling · 3 of 12 trajectories · about 2 min left') === 'Zuggeführte Stichprobe · 3 von 12 Trajektorien · noch etwa 2 min', 'stage, detail and time left');
  ok(t2('40% · about 2 min left · 5 min elapsed') === '40 % · noch etwa 2 min · 5 min vergangen', 'progress line');
  ok(t2('2 hypotheses (seed structures only (ViennaRNA unavailable))') === '2 Hypothesen (nur Startstrukturen (ViennaRNA nicht verfügbar))', 'a translated value inside a pattern');
  ok(t2('101 × 32 (+8 Z-RNA) trajectories · about 5 h 44 min on 15 cores · at most 100 sequences per batch').endsWith('höchstens 100 Sequenzen pro Stapel'), 'estimate of a batch that is too large');
  ok(t2('variant "G1C": position 1 is G, not C') === 'Variante „G1C“: Position 1 ist G, nicht C', 'reason why a variant was not analysed');
  ok(t2('Compare two models') === 'Compare two models', 'a sentence is not taken for "{0} models"');
  ok(t2('Failed:') === 'Fehlgeschlagen:' && t2('some engine message without a translation') === 'some engine message without a translation', '"Failed:" is translated on its own');
}
await i18n.setLanguage('ja', { save: false });
ok(i18n.t('variant "G1C": position 1 is G, not C') === '変異体「G1C」: 位置 1 は C ではなく G です', 'placeholders may change their order');
// Durations (0.7.6, audit F6): "40 s" and "2 min" follow the language as "1 h 5 min" does; "{0} s" takes a number only
ok(i18n.t('40 s') === '40 秒' && i18n.t('2 min') === '2 分' && i18n.t('1.5 h') === '1.5 時間' && i18n.t('1 h 5 min') === '1 時間 5 分', 'durations in Japanese');
ok(i18n.t('40% · about 2 min left · 35 s elapsed') === '40% · 残り約 2 分 · 経過 35 秒', 'durations inside a progress line');
ok(i18n.t('Relaxation and validation · about 40 s left') === '緩和と検証 · 残り約 40 秒', 'a duration without three letters inside a line');
ok(i18n.t('ViennaRNA took 40 s') === 'ViennaRNA took 40 s' && i18n.t('tRNA s') === 'tRNA s', '"{0} s" does not catch a sentence that ends in " s"');
await i18n.setLanguage('ko', { save: false });
ok(i18n.t('40 s') === '40초' && i18n.t('2 min') === '2분', 'durations in Korean');
await i18n.setLanguage('de', { save: false });
ok(i18n.t('40 s') === '40 s' && i18n.t('2 min') === '2 min', 'German keeps the unit symbols');
// Data inside a text (audit F5): put aside while the text is translated, back unchanged
{
  const title = 'STRUCTURAL CHARACTERIZATION of the rev response element';
  const en = `Template 1I9F_A (${title}): sequence identity 97%, covers 22% of the sequence, E-value 0.028.`;
  ok(i18n.t(en) === en, 'a title in English prose stops the pattern when it is not named as data');
  ok(i18n.tKeep(en, [title]) === `Template 1I9F_A (${title}): Sequenzidentität 97%, deckt 22% der Sequenz ab, E-Wert 0.028.`, 'named as data, the title stays and the rest is translated');
  ok(i18n.tKeep('5 conformations', ['5 conformations']) === '5 conformations', 'a text that is all data is not translated');
  ok(i18n.tKeep('5 conformations', []) === '5 Konformationen' && i18n.tKeep('5 conformations', ['not in the text']) === '5 Konformationen', 'no data: as t()');
  ok(i18n.tKeep('Command: run $1 $& now', ['run $1 $& now']) === 'Befehl: run $1 $& now', 'data with "$" comes back as it was');
}
if (failures) { console.error(`${failures} failure(s)`); process.exit(1); }
console.log('all i18n tests passed');
