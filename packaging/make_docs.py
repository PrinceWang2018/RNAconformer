"""Render RCV's Markdown documentation as standalone HTML pages (installers and the website).

    python packaging/make_docs.py OUT_DIR [--site]

Writes OUT_DIR/index.html (user guide, docs/RCV.md), plugins.html (docs/PLUGINS.md), changelog.html
(CHANGELOG.md) and images/. --site adds the website header links. No third-party packages: the converter covers
the Markdown used in these files (headings, paragraphs, lists, tables, code blocks, inline code, bold, italics,
links, images, rules).
"""
from __future__ import annotations

import argparse
import html
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from rcv import WEBSITE, __version__  # noqa: E402

PAGES = [('index.html', ROOT / 'docs' / 'RCV.md', 'User guide'),
         ('plugins.html', ROOT / 'docs' / 'PLUGINS.md', 'Plugin developer guide'),
         ('changelog.html', ROOT / 'CHANGELOG.md', 'What is new')]
LINKS = {'RCV.md': 'index.html', 'PLUGINS.md': 'plugins.html', 'CHANGELOG.md': 'changelog.html'}
GITHUB = 'https://github.com/PrinceWang2018/RNAconformer/blob/master/'


def slugify(text):
    text = re.sub(r'[`*]', '', text.lower())
    text = re.sub(r'[^\w\s-]', '', text)
    return re.sub(r'\s+', '-', text.strip())


def link_target(href):
    if re.match(r'^(https?:|mailto:|#)', href):
        return href
    path, _, anchor = href.partition('#')
    name = path.split('/')[-1]
    if name in LINKS:
        return LINKS[name] + ('#' + anchor if anchor else '')
    if path.startswith('images/') or path.startswith('docs/images/'):
        return 'images/' + name
    return GITHUB + path.lstrip('./') + ('#' + anchor if anchor else '')     # other repository files


def inline(text):
    out, pos = [], 0
    pattern = re.compile(r'(`[^`]+`|\*\*[^*]+\*\*|\*[^*\s][^*]*\*|!\[[^\]]*\]\([^)]+\)|\[[^\]]+\]\([^)]+\))')
    for m in pattern.finditer(text):
        out.append(html.escape(text[pos:m.start()], quote=False))
        s = m.group(0)
        if s.startswith('`'):
            out.append(f'<code>{html.escape(s[1:-1], quote=False)}</code>')
        elif s.startswith('**'):
            out.append(f'<b>{inline(s[2:-2])}</b>')
        elif s.startswith('!['):
            alt, src = s[2:s.index(']')], s[s.index('](') + 2:-1].split()[0]
            out.append(f'<img src="{html.escape(link_target(src))}" alt="{html.escape(alt)}" loading="lazy">')
        elif s.startswith('['):
            label, href = s[1:s.index(']')], s[s.index('](') + 2:-1].strip()
            ext = ' target="_blank" rel="noopener"' if href.startswith('http') else ''
            out.append(f'<a href="{html.escape(link_target(href))}"{ext}>{inline(label)}</a>')
        else:
            out.append(f'<i>{inline(s[1:-1])}</i>')
        pos = m.end()
    out.append(html.escape(text[pos:], quote=False))
    return ''.join(out)


def markdown(md):
    lines = md.replace('\r', '').split('\n')
    out, para, toc, i = [], [], [], 0

    def flush():
        if para:
            out.append(f'<p>{inline(" ".join(para))}</p>')
            para.clear()
    while i < len(lines):
        line = lines[i]
        if line.startswith('```'):
            flush()
            code = []
            i += 1
            while i < len(lines) and not lines[i].startswith('```'):
                code.append(lines[i])
                i += 1
            out.append(f'<pre><code>{html.escape(chr(10).join(code), quote=False)}</code></pre>')
            i += 1
            continue
        h = re.match(r'^(#{1,4})\s+(.*)$', line)
        if h:
            flush()
            level, text = len(h.group(1)), h.group(2)
            sid = slugify(text)
            out.append(f'<h{level} id="{sid}">{inline(text)}</h{level}>')
            if level in (2, 3):
                toc.append((level, sid, re.sub(r'[`*]', '', text)))
            i += 1
            continue
        if re.match(r'^---+\s*$', line):
            flush()
            out.append('<hr>')
            i += 1
            continue
        if line.startswith('|'):
            flush()
            rows = []
            while i < len(lines) and lines[i].startswith('|'):
                rows.append(lines[i])
                i += 1
            cells = lambda r: [c.strip().replace('\\|', '|') for c in re.split(r'(?<!\\)\|', r.strip().strip('|'))]  # noqa: E731
            body = [r for k, r in enumerate(rows) if not (k == 1 and re.match(r'^\|[\s:|-]+\|$', r.strip()))]
            head = ''.join(f'<th>{inline(c)}</th>' for c in cells(body[0]))
            rest = ''.join('<tr>' + ''.join(f'<td>{inline(c)}</td>' for c in cells(r)) + '</tr>' for r in body[1:])
            out.append(f'<div class="tw"><table><thead><tr>{head}</tr></thead><tbody>{rest}</tbody></table></div>')
            continue
        if re.match(r'^\s*([*-]|\d+\.)\s+', line):
            flush()
            ordered = bool(re.match(r'^\s*\d+\.', line))
            items = []
            while i < len(lines) and re.match(r'^\s*([*-]|\d+\.)\s+', lines[i]):
                item = re.sub(r'^\s*([*-]|\d+\.)\s+', '', lines[i])
                i += 1
                while i < len(lines) and re.match(r'^\s{2,}\S', lines[i]) and not re.match(r'^\s*([*-]|\d+\.)\s+', lines[i]):
                    item += ' ' + lines[i].strip()
                    i += 1
                items.append(f'<li>{inline(item)}</li>')
            tag = 'ol' if ordered else 'ul'
            out.append(f'<{tag}>{"".join(items)}</{tag}>')
            continue
        if not line.strip():
            flush()
        else:
            para.append(line.strip())
        i += 1
    flush()
    return '\n'.join(out), toc


CSS = """
:root { --bg: #f7f6f2; --card: #fff; --text: #1d2430; --muted: #5d6675; --line: #e3e1da; --accent: #2a78d6; --code: #f1efe9; }
@media (prefers-color-scheme: dark) { :root { --bg: #16191e; --card: #1d2128; --text: #e7e9ee; --muted: #a2a9b5; --line: #2f343d; --accent: #6aa8f0; --code: #242932; } }
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--text); font: 15px/1.6 system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif; }
header { display: flex; align-items: center; gap: 14px; padding: 12px 22px; background: var(--card); border-bottom: 1px solid var(--line); position: sticky; top: 0; z-index: 2; }
header img { width: 26px; height: 26px; }
header b { font-size: 16px; letter-spacing: .02em; }
header nav { display: flex; gap: 16px; margin-left: auto; flex-wrap: wrap; }
header a { color: var(--muted); text-decoration: none; font-size: 14px; }
header a:hover, header a.on { color: var(--accent); }
.wrap { display: grid; grid-template-columns: 260px minmax(0, 1fr); gap: 32px; max-width: 1200px; margin: 0 auto; padding: 24px 22px 60px; }
aside { position: sticky; top: 72px; align-self: start; max-height: calc(100vh - 90px); overflow: auto; font-size: 13.5px; }
aside a { display: block; padding: 3px 8px; color: var(--muted); text-decoration: none; border-radius: 5px; }
aside a.l3 { padding-left: 20px; font-size: 13px; }
aside a:hover { background: var(--code); color: var(--text); }
main { background: var(--card); border: 1px solid var(--line); border-radius: 10px; padding: 8px 34px 34px; min-width: 0; }
h1 { font-size: 28px; margin: 22px 0 8px; } h2 { font-size: 22px; margin: 34px 0 10px; padding-top: 8px; border-top: 1px solid var(--line); }
h3 { font-size: 17px; margin: 24px 0 8px; } h4 { font-size: 15px; margin: 18px 0 6px; }
a { color: var(--accent); }
code { background: var(--code); padding: 1px 5px; border-radius: 4px; font: 13px/1.5 ui-monospace, Consolas, "Liberation Mono", monospace; }
pre { background: var(--code); padding: 12px 14px; border-radius: 8px; overflow: auto; } pre code { padding: 0; background: none; }
.tw { overflow-x: auto; } table { border-collapse: collapse; margin: 12px 0; font-size: 14px; }
th, td { border: 1px solid var(--line); padding: 6px 10px; text-align: left; vertical-align: top; } th { background: var(--code); }
img { max-width: 100%; height: auto; border: 1px solid var(--line); border-radius: 8px; margin: 8px 0; }
footer { text-align: center; color: var(--muted); font-size: 13px; padding: 0 0 30px; }
@media (max-width: 860px) { .wrap { grid-template-columns: 1fr; } aside { display: none; } main { padding: 4px 16px 24px; } }
"""


def page(title, body, toc, current, site):
    nav = ''.join(f'<a href="{f}"{" class=on" if f == current else ""}>{t}</a>' for f, _, t in PAGES)
    if site:
        nav = f'<a href="../">Home</a><a href="../#download">Download</a>' + nav
    else:
        nav += f'<a href="{WEBSITE}" target="_blank" rel="noopener">Website</a>'
    side = ''.join(f'<a class="l{lvl}" href="#{sid}">{html.escape(text)}</a>' for lvl, sid, text in toc)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)} – RCV {__version__}</title><link rel="icon" href="images/icon.svg"><style>{CSS}</style></head>
<body><header><img src="images/icon.svg" alt=""><b>RCV</b><span style="color:var(--muted)">RNA Conformation Viewer {__version__}</span><nav>{nav}</nav></header>
<div class="wrap"><aside>{side}</aside><main>{body}</main></div>
<footer>© 2026 Zixiang Wang · RCV is free software under the GPL-3.0 license · <a href="{WEBSITE}">{WEBSITE.split('//')[1].rstrip('/')}</a></footer>
</body></html>
"""


def build(out, site=False):
    out = Path(out)
    (out / 'images').mkdir(parents=True, exist_ok=True)
    for img in (ROOT / 'docs' / 'images').iterdir():
        if img.suffix.lower() in ('.png', '.svg', '.jpg', '.webp'):
            shutil.copy2(img, out / 'images' / img.name)
    shutil.copy2(ROOT / 'src' / 'rcv' / 'static' / 'img' / 'icon.svg', out / 'images' / 'icon.svg')
    for name, src, title in PAGES:
        body, toc = markdown(src.read_text(encoding='utf-8'))
        (out / name).write_text(page(title, body, toc, name, site), encoding='utf-8')
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('out')
    ap.add_argument('--site', action='store_true', help='links for the website (Home, Download)')
    a = ap.parse_args()
    print('documentation written to', build(a.out, a.site))
