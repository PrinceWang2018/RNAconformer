// Interface languages.
//
// English is the source language. A translation file (static/i18n/<code>.json) maps English interface text to the
// chosen language, and this module translates the interface as it is drawn: text, tooltips (title), placeholders
// and accessible names, in the start page, menus, dialogs and panels. Text with numbers or names in it is matched
// with patterns such as "{0} conformations · {1} nt". Anything without a translation stays in English.
//
// Kept in English on purpose: figures and charts (SVG and canvas, so an exported figure looks the same in every
// language), the command console and its commands, structure and file data, and anything marked translate="no".
import { loadPref, savePref } from './util.js';

export const LANGUAGES = [
  { code: 'en', name: 'English', english: 'English' },
  { code: 'zh-CN', name: '简体中文', english: 'Chinese (Simplified)' },
  { code: 'zh-TW', name: '繁體中文', english: 'Chinese (Traditional)' },
  { code: 'ja', name: '日本語', english: 'Japanese' },
  { code: 'ko', name: '한국어', english: 'Korean' },
  { code: 'es', name: 'Español', english: 'Spanish' },
  { code: 'fr', name: 'Français', english: 'French' },
  { code: 'de', name: 'Deutsch', english: 'German' },
  { code: 'it', name: 'Italiano', english: 'Italian' },
  { code: 'pt-BR', name: 'Português (Brasil)', english: 'Portuguese (Brazil)' },
  { code: 'ru', name: 'Русский', english: 'Russian' },
];

const ATTRS = ['title', 'placeholder', 'aria-label', 'label'];
const SKIP_TAGS = new Set(['SCRIPT', 'STYLE', 'TEXTAREA', 'CODE', 'PRE', 'KBD', 'CANVAS', 'TITLE']);
const ATTRS_OK = new Set(['TEXTAREA']);   // the text is the user's, but the placeholder and tooltip are interface text
const SVG_NS = 'http://www.w3.org/2000/svg';

let lang = 'en';
let dict = null;                 // Map: English -> translation
let patterns = [];               // [{ re, out, anchor }]
let cache = new Map();           // English -> translation | null (no translation)
const textSrc = new WeakMap();   // text node -> English source of the translation shown
const textOut = new WeakMap();   // text node -> translation this module wrote
const attrSrc = new WeakMap();   // element -> { attr: English source }
const attrOut = new WeakMap();   // element -> { attr: translation written }
let observer = null;
const listeners = new Set();

// ------------------------------------------------------------------ language choice
export function supported(code) { return LANGUAGES.some(l => l.code === code); }

// Best match of the system languages, e.g. zh-HK -> zh-TW, pt-PT -> pt-BR, de-AT -> de.
export function matchLanguage(tags) {
  for (const raw of tags || []) {
    const t = String(raw || '').toLowerCase();
    if (!t) continue;
    if (t.startsWith('zh')) return /hant|tw|hk|mo/.test(t) ? 'zh-TW' : 'zh-CN';
    if (t.startsWith('pt')) return 'pt-BR';
    const base = t.split('-')[0];
    const hit = LANGUAGES.find(l => l.code.toLowerCase() === base);
    if (hit) return hit.code;
  }
  return 'en';
}

// ?lang=xx in the address > saved choice > system language.
export function preferredLanguage() {
  try {
    const q = new URLSearchParams(location.search).get('lang');
    if (q && supported(q)) return q;
  } catch { /* no URL */ }
  const saved = loadPref('lang', null);
  if (saved && supported(saved)) return saved;
  return matchLanguage(navigator.languages?.length ? navigator.languages : [navigator.language]);
}

export const currentLanguage = () => lang;
export const languageName = code => (LANGUAGES.find(l => l.code === code) || LANGUAGES[0]).name;
export function onLanguage(fn) { listeners.add(fn); return () => listeners.delete(fn); }

// ------------------------------------------------------------------ translation
const escRe = s => s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');

function compile(strings) {
  dict = new Map();
  patterns = [];
  cache = new Map();
  for (const [src, out] of Object.entries(strings || {})) {
    if (typeof out !== 'string' || !out || src.startsWith('_')) continue;
    if (/\{\d\}/.test(src)) {
      const parts = src.split(/(\{\d\})/);
      const order = [];
      const re = parts.map(p => { const m = /^\{(\d)\}$/.exec(p); if (m) { order.push(+m[1]); return '([\\s\\S]*?)'; } return escRe(p); }).join('');
      const lits = parts.filter(p => p && !/^\{\d\}$/.test(p));
      const anchor = [...lits].sort((a, b) => b.length - a.length)[0] || '';
      const weight = lits.reduce((n, p) => n + p.length, 0);        // more literal text = more specific
      patterns.push({ re: new RegExp(`^${re}$`), order, out, anchor, weight });
    } else dict.set(src, out);
  }
  patterns.sort((a, b) => b.weight - a.weight);
}

// Translation of one interface string (surrounding spaces are kept). English and unknown text come back unchanged.
export function t(s) {
  if (!dict || typeof s !== 'string') return s;
  const m = /^(\s*)([\s\S]*?)(\s*)$/.exec(s);
  const core = m[2];
  if (!core || !/[A-Za-z]{2}/.test(core)) return s;
  let out = cache.get(core);
  if (out === undefined) {
    out = lookup(core);
    if (out === null) {
      // a keyboard hint after the text: "Undo (Ctrl+Z)", "Pull a residue to switch conformation (G)"
      const k = KEY_HINT.exec(core);
      if (k) { const base = lookup(k[1]); if (base !== null) out = base + k[2]; }
    }
    if (cache.size > 6000) cache.clear();
    cache.set(core, out);
  }
  return out === null ? s : m[1] + out + m[3];
}

const KEY_HINT = /^([\s\S]*?\S)(\s*\((?:(?:Ctrl|Shift|Alt)\+[^()]+|F\d{1,2}|Esc|Space|[A-Z0-9`?←→])\))$/;
// English prose left in a filled-in value means the pattern caught a sentence it was not written for.
const PROSE = /[a-z]{3,}[\s,;:]+[a-z]{3,}/;

function lookup(core, depth = 0) {
  const hit = dict.get(core);
  if (hit !== undefined) return hit;
  for (const p of patterns) {
    if (p.anchor && !core.includes(p.anchor)) continue;
    const r = p.re.exec(core);
    if (!r) continue;
    const vals = {};
    let ok = true;
    p.order.forEach((k, i) => {
      const v = r[i + 1], vt = v.trim();
      const tv = dict.get(vt) ?? (depth < 2 && /[A-Za-z]{3}/.test(vt) ? (lookup(vt, depth + 1) ?? undefined) : undefined);
      if (tv === undefined && PROSE.test(v)) ok = false;
      vals[k] = tv === undefined ? v : v.replace(v.trim(), tv);
    });
    if (!ok) continue;
    return p.out.replace(/\{(\d)\}/g, (_, k) => vals[k] ?? '');
  }
  return null;
}

// Fill a pattern in code: tf('{0} conformations', 12)
export function tf(src, ...vals) {
  const s = src.replace(/\{(\d)\}/g, (_, k) => String(vals[+k] ?? ''));
  return t(s);
}

// ------------------------------------------------------------------ DOM
function skipped(node, forAttrs = false) {
  for (let e = node.nodeType === 1 ? node : node.parentNode; e && e.nodeType === 1; e = e.parentNode) {
    if (SKIP_TAGS.has(e.tagName) && !(forAttrs && e === node && ATTRS_OK.has(e.tagName))) return true;
    if (e.namespaceURI === SVG_NS) return true;
    if (e.getAttribute('translate') === 'no' || e.isContentEditable) return true;
  }
  return false;
}

function doText(node, force = false) {
  const cur = node.data;
  let src;
  if (textOut.get(node) === cur) { if (!force) return; src = textSrc.get(node) ?? cur; } else src = cur;
  const out = lang === 'en' ? src : t(src);
  if (out !== src) { textSrc.set(node, src); textOut.set(node, out); }
  else { textSrc.delete(node); textOut.delete(node); }
  if (out !== cur) node.data = out;
}

function doAttr(elm, name, force = false) {
  const cur = elm.getAttribute(name);
  if (cur === null) return;
  const outs = attrOut.get(elm) || {};
  const srcs = attrSrc.get(elm) || {};
  let src;
  if (outs[name] === cur) { if (!force) return; src = srcs[name] ?? cur; } else src = cur;
  const out = lang === 'en' ? src : t(src);
  if (out !== src) { srcs[name] = src; outs[name] = out; } else { delete srcs[name]; delete outs[name]; }
  attrSrc.set(elm, srcs); attrOut.set(elm, outs);
  if (out !== cur) elm.setAttribute(name, out);
}

function walk(root, force = false) {
  if (!root) return;
  if (root.nodeType === 3) { if (!skipped(root)) doText(root, force); return; }
  if (root.nodeType !== 1 && root.nodeType !== 9 && root.nodeType !== 11) return;
  if (root.nodeType === 1 && skipped(root, true)) return;
  const filter = {
    acceptNode(n) {
      if (n.nodeType === 3) return SKIP_TAGS.has(n.parentNode?.tagName) ? NodeFilter.FILTER_REJECT : NodeFilter.FILTER_ACCEPT;
      if (n.namespaceURI === SVG_NS || n.getAttribute('translate') === 'no') return NodeFilter.FILTER_REJECT;
      if (SKIP_TAGS.has(n.tagName) && !ATTRS_OK.has(n.tagName)) return NodeFilter.FILTER_REJECT;
      return NodeFilter.FILTER_ACCEPT;
    },
  };
  if (root.nodeType === 1) for (const a of ATTRS) if (root.hasAttribute(a)) doAttr(root, a, force);
  const w = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT | NodeFilter.SHOW_TEXT, filter);
  for (let n = w.nextNode(); n; n = w.nextNode()) {
    if (n.nodeType === 3) { if (n.data.trim()) doText(n, force); }
    else for (const a of ATTRS) if (n.hasAttribute(a)) doAttr(n, a, force);
  }
}

function startObserver() {
  if (observer) return;
  observer = new MutationObserver(muts => {
    if (lang === 'en' && !dict) return;
    for (const m of muts) {
      if (m.type === 'childList') { for (const n of m.addedNodes) walk(n); }
      else if (m.type === 'characterData') { if (!skipped(m.target)) doText(m.target); }
      else if (m.type === 'attributes') { if (!skipped(m.target, true)) doAttr(m.target, m.attributeName); }
    }
  });
  observer.observe(document.documentElement, { subtree: true, childList: true, characterData: true, attributes: true, attributeFilter: ATTRS });
}

async function loadStrings(code) {
  if (code === 'en') return null;
  const r = await fetch(`i18n/${code}.json`, { cache: 'no-cache' });
  if (!r.ok) throw new Error(`no translation file for ${code}`);
  const data = await r.json();
  return data.strings || data;
}

// Switch the interface language (also used at start-up). Saves the choice unless save is false.
export async function setLanguage(code, { save = true } = {}) {
  if (!supported(code)) code = 'en';
  let strings = null;
  try { strings = await loadStrings(code); } catch (err) { console.warn(err); code = 'en'; }
  lang = code;
  if (strings) compile(strings); else { dict = null; patterns = []; cache = new Map(); }
  document.documentElement.lang = code;
  if (save) savePref('lang', code);
  startObserver();
  walk(document.body, true);
  for (const fn of [...listeners]) { try { fn(code); } catch (err) { console.error(err); } }
  return code;
}

export function initLanguage() { return setLanguage(preferredLanguage(), { save: false }); }

// Translation coverage of the current page (for tests): English strings still visible in a translated interface.
export function untranslated(root = document.body) {
  const out = new Set();
  if (lang === 'en' || !root) return [];
  const w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
  for (let n = w.nextNode(); n; n = w.nextNode()) {
    const s = n.data.trim();
    if (s && /[A-Za-z]{3}/.test(s) && !skipped(n) && !textOut.has(n) && n.parentElement?.offsetParent !== null) out.add(s);
  }
  return [...out];
}
