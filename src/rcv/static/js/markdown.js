// Minimal Markdown renderer for RCV's own documentation shown inside the application (headings, paragraphs, lists,
// tables, fenced code, inline code, bold, italics, images, links). Builds DOM nodes; never uses innerHTML.
//   opts.image(src) -> URL or Promise<URL> for ![alt](src); without it images are left out.
//   opts.link(href) is called when a link is clicked; without it links are shown as text.
import { el } from './util.js';

export const slugify = s => String(s).toLowerCase().replace(/[`*]/g, '').replace(/[^\p{L}\p{N}\s-]/gu, '').trim().replace(/\s+/g, '-');

function inline(text, opts = {}) {
  const out = [];
  const re = /(`[^`]+`|\*\*[^*]+\*\*|\*[^*\s][^*]*\*|!\[[^\]]*\]\([^)]+\)|\[[^\]]+\]\([^)]+\))/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(text.slice(last, m.index));
    const s = m[0];
    if (s.startsWith('`')) out.push(el('code', { text: s.slice(1, -1) }));
    else if (s.startsWith('**')) out.push(el('b', {}, ...inline(s.slice(2, -2), opts)));
    else if (s.startsWith('![')) {
      const alt = s.slice(2, s.indexOf(']')), src = s.slice(s.indexOf('](') + 2, -1).trim().split(/\s+/)[0];
      if (opts.image) {
        const img = el('img', { alt, loading: 'lazy', draggable: 'false' });
        Promise.resolve(opts.image(src)).then(u => { if (u) img.src = u; else img.remove(); }).catch(() => img.remove());
        out.push(img);
      }
    } else if (s.startsWith('[')) {
      const label = s.slice(1, s.indexOf(']')), href = s.slice(s.indexOf('](') + 2, -1).trim();
      if (opts.link) out.push(el('a', { class: 'md-link', role: 'link', tabindex: 0, title: href, onclick: e => { e.preventDefault(); opts.link(href); } }, ...inline(label, opts)));
      else out.push(el('span', { class: 'md-link', text: label }));
    } else out.push(el('i', {}, ...inline(s.slice(1, -1), opts)));
    last = m.index + s.length;
  }
  if (last < text.length) out.push(text.slice(last));
  return out;
}

export function renderMarkdown(md, opts = {}) {
  const inl = t => inline(t, opts);
  const root = el('div', { class: 'md selectable' });
  const lines = String(md).replace(/\r/g, '').split('\n');
  let i = 0;
  const para = [];
  const flush = () => { if (para.length) { root.append(el('p', {}, ...inl(para.join(' ')))); para.length = 0; } };
  while (i < lines.length) {
    const l = lines[i];
    if (/^```/.test(l)) {
      flush();
      const code = [];
      i++;
      while (i < lines.length && !/^```/.test(lines[i])) code.push(lines[i++]);
      root.append(el('pre', {}, el('code', { text: code.join('\n') })));
      i++; continue;
    }
    const h = /^(#{1,4})\s+(.*)$/.exec(l);
    if (h) { flush(); root.append(el(`h${Math.min(5, h[1].length + 1)}`, { id: slugify(h[2]) }, ...inl(h[2]))); i++; continue; }
    if (/^---+\s*$/.test(l)) { flush(); root.append(el('hr')); i++; continue; }
    if (/^\|/.test(l)) {
      flush();
      const rows = [];
      while (i < lines.length && /^\|/.test(lines[i])) rows.push(lines[i++]);
      const cells = r => r.trim().replace(/^\||\|$/g, '').split(/(?<!\\)\|/).map(c => c.trim().replace(/\\\|/g, '|'));
      const body = rows.filter((r, k) => !(k === 1 && /^\|[\s:|-]+\|$/.test(r.trim())));
      const t = el('table', { class: 'table' }, el('thead', {}, el('tr', {}, ...cells(body[0]).map(c => el('th', {}, ...inl(c))))),
        el('tbody', {}, ...body.slice(1).map(r => el('tr', {}, ...cells(r).map(c => el('td', {}, ...inl(c)))))));
      root.append(t);
      continue;
    }
    if (/^\s*([*-]|\d+\.)\s+/.test(l)) {
      flush();
      const ordered = /^\s*\d+\./.test(l);
      const list = el(ordered ? 'ol' : 'ul');
      while (i < lines.length && /^\s*([*-]|\d+\.)\s+/.test(lines[i])) {
        let item = lines[i].replace(/^\s*([*-]|\d+\.)\s+/, '');
        i++;
        while (i < lines.length && /^\s{2,}\S/.test(lines[i]) && !/^\s*([*-]|\d+\.)\s+/.test(lines[i])) item += ' ' + lines[i++].trim();
        list.append(el('li', {}, ...inl(item)));
      }
      root.append(list);
      continue;
    }
    if (!l.trim()) { flush(); i++; continue; }
    para.push(l.trim());
    i++;
  }
  flush();
  return root;
}
