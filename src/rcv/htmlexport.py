"""Interactive HTML export: one self-contained page (RCV Lite) with the ensemble, the 3D viewer and the 2D drawing.

The page needs no server and no installation: the viewer's JavaScript modules (and three.js, MIT licence) are
embedded base64-encoded and started from blob URLs, the ensemble document is embedded as JSON. Suitable for
supplementary material, reviewers and collaborators.
"""
from __future__ import annotations

import base64
import html
import json

from . import __version__
from . import paths as P

MODULES = ['js/lite.js', 'js/viewer/viewer.js', 'js/viewer/controls.js', 'js/viewer/post.js', 'js/viewer/geometry.js',
           'js/colors.js', 'js/util.js', 'js/ensemble.js', 'js/ss2dlayout.js', 'js/landscape.js', 'js/mechmorph.js',
           'vendor/three/three.module.min.js']

BOOT = r"""(function () {
  var mods = {};
  document.querySelectorAll('script[type="text/x-rcv-module"]').forEach(function (s) { mods[s.getAttribute('data-path')] = s.textContent; });
  function dec(b64) { var bin = atob(b64.replace(/\s+/g, '')); var u = new Uint8Array(bin.length); for (var i = 0; i < bin.length; i++) u[i] = bin.charCodeAt(i); return new TextDecoder().decode(u); }
  function resolve(base, spec) { var p = base.split('/'); p.pop(); spec.split('/').forEach(function (x) { if (x === '..') p.pop(); else if (x !== '.') p.push(x); }); return p.join('/'); }
  var urls = {};
  function url(path) {
    if (urls[path]) return urls[path];
    if (!(path in mods)) throw new Error('missing module ' + path);
    var src = dec(mods[path]).replace(/(\bfrom\s*|\bimport\s*)(['"])(\.{1,2}\/[^'"]+)\2/g, function (m, kw, q, spec) { return kw + q + url(resolve(path, spec)) + q; });
    urls[path] = URL.createObjectURL(new Blob([src], { type: 'text/javascript' }));
    return urls[path];
  }
  function fail(e) { var b = document.getElementById('err'); b.textContent = 'This page could not start the 3D viewer (' + (e && e.message || e) + '). Open it in a current Chrome, Edge, Firefox or Safari.'; b.hidden = false; }
  try { import(url('js/lite.js')).catch(fail); } catch (e) { fail(e); }
})();"""

STYLE = r"""
:root { color-scheme: light; --page: #f3f2ee; --surface: #fcfcfb; --surface-2: #f6f5f1; --border: rgba(11,11,11,.11); --text: #0b0b0b;
  --text-2: #52514e; --text-3: #898781; --accent: #2a78d6; --accent-weak: rgba(42,120,214,.12); --axis: #a8a69f; --border-strong: rgba(11,11,11,.22);
  --series-1: #2a78d6; --series-2: #eb6834; --font: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  --mono: ui-monospace, "Cascadia Mono", Consolas, monospace; }
@media (prefers-color-scheme: dark) { :root { color-scheme: dark; --page: #0d0d0d; --surface: #1a1a19; --surface-2: #151514; --border: rgba(255,255,255,.1);
  --text: #f0efec; --text-2: #c3c2b7; --text-3: #8f8d86; --accent: #3987e5; --accent-weak: rgba(57,135,229,.18); --axis: #6b6a65;
  --border-strong: rgba(255,255,255,.22); --series-1: #3987e5; --series-2: #f08456; } }
* { box-sizing: border-box; }
html, body { margin: 0; height: 100%; background: var(--page); color: var(--text); font: 13px/1.4 var(--font); }
body { display: grid; grid-template-rows: auto minmax(0, 1fr); overflow: hidden; }
header { display: flex; align-items: center; gap: 10px; padding: 8px 16px; background: var(--surface); border-bottom: 1px solid var(--border); min-width: 0; }
header svg { width: 26px; height: 26px; flex: none; }
header h1 { font-size: 15px; margin: 0; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
header .sub { color: var(--text-2); font-size: 12px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
header .grow { flex: 1; }
header .badge { font-size: 11px; color: var(--text-3); white-space: nowrap; }
main { display: grid; grid-template-columns: minmax(0, 1fr) 340px; min-height: 0; }
#vp { position: relative; min-width: 0; min-height: 0; overflow: hidden; }
#gl { width: 100%; height: 100%; display: block; touch-action: none; }
#ov { position: absolute; inset: 0; pointer-events: none; width: 100%; height: 100%; }
#lb { position: absolute; inset: 0; pointer-events: none; }
.res-label { position: absolute; transform: translate(-50%, -130%); font: 600 11px var(--font); background: rgba(255,255,255,.9); color: #111; padding: 1px 5px; border-radius: 4px; }
#legend { position: absolute; left: 14px; bottom: 12px; display: flex; flex-wrap: wrap; gap: 4px 10px; align-items: center; max-width: 70%; font-size: 12px; color: #141414; }
#legend .lg-title { width: 100%; font-size: 10.5px; text-transform: uppercase; letter-spacing: .06em; opacity: .7; }
#legend .lg { display: inline-flex; gap: 5px; align-items: center; } #legend .sw { width: 12px; height: 12px; border-radius: 3px; }
#legend .ramp { width: 120px; height: 10px; border-radius: 3px; }
#tip { position: fixed; z-index: 9; pointer-events: none; background: rgba(20,20,20,.9); color: #fff; padding: 4px 8px; border-radius: 6px; font-size: 12px; }
aside { background: var(--surface); border-left: 1px solid var(--border); overflow: auto; padding: 12px 14px; display: flex; flex-direction: column; gap: 12px; }
h2 { font-size: 10.5px; text-transform: uppercase; letter-spacing: .08em; color: var(--text-2); margin: 0 0 6px; }
.row { display: flex; gap: 6px; align-items: center; flex-wrap: wrap; }
button, select { font: inherit; color: inherit; background: var(--surface-2); border: 1px solid var(--border); border-radius: 6px; padding: 4px 10px; cursor: pointer; }
button:hover { border-color: var(--accent); }
label.check { display: inline-flex; gap: 5px; align-items: center; font-size: 12px; }
.seq { font: 11.5px/1.5 var(--mono); word-break: break-all; color: var(--text-2); }
#confs { display: flex; flex-direction: column; gap: 2px; max-height: 260px; overflow: auto; }
.conf { display: grid; grid-template-columns: 22px 1fr auto; gap: 0 8px; text-align: left; border-color: transparent; padding: 4px 8px; }
.conf.on { background: var(--accent-weak); border-color: var(--accent); }
.conf .n { color: var(--text-3); } .conf .v { font-variant-numeric: tabular-nums; }
.conf .dbn { grid-column: 2 / 4; font: 10.5px var(--mono); color: var(--text-3); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
#ss { width: 100%; height: 260px; display: block; background: var(--surface-2); border-radius: 8px; }
.cap { font-size: 11.5px; color: var(--text-3); margin-top: 4px; }
#profile { width: 100%; height: 90px; display: block; }
#methods { font-size: 12px; color: var(--text-2); }
#map { width: 100%; height: 220px; display: block; cursor: pointer; }
#pw-prof { width: 100%; height: 110px; display: block; }
#pw-phi { width: 100%; display: block; }
.pw-text { font-size: 12px; color: var(--text-2); display: flex; flex-direction: column; gap: 4px; }
.pw-text b { color: var(--text); }
#made { font-size: 11px; color: var(--text-3); margin-top: auto; }
#err { position: fixed; inset: 40% 10% auto; background: #b42d2d; color: #fff; padding: 14px 18px; border-radius: 8px; z-index: 10; }
[hidden] { display: none !important; }
@media (max-width: 820px) { main { grid-template-columns: 1fr; grid-template-rows: 55vh auto; overflow: auto; } body { overflow: auto; } aside { border-left: 0; border-top: 1px solid var(--border); } }
"""

BODY = """<header>{icon}<div style="min-width:0"><h1 id="title"></h1><div class="sub" id="subtitle"></div></div><span class="grow"></span>
<span class="badge">drag: rotate · right-drag: pan · wheel: zoom · ←/→: morph · space: play</span></header>
<main>
<div id="vp"><canvas id="gl" aria-label="3D structure"></canvas><svg id="ov"></svg><div id="lb"></div><div id="legend"></div></div>
<aside>
<section><h2>Conformations</h2>
<div class="row"><button id="prev" title="Previous (←)">◀</button><button id="play">Play</button><button id="next" title="Next (→)">▶</button><button id="reset">Reset view</button></div>
<div class="row" style="margin-top:6px"><select id="rep" aria-label="Representation"></select><select id="scheme" aria-label="Colours"></select>
<label class="check"><input type="checkbox" id="overlay"> overlay</label></div>
<div class="row" style="margin-top:6px"><select id="metric" aria-label="Energy" style="max-width:100%"></select></div>
<div id="confs" style="margin-top:6px"></div></section>
<section><svg id="profile"></svg><div class="cap" id="profile-cap"></div></section>
<section id="map-box"><h2>Ensemble map</h2><svg id="map" role="img" aria-label="Topology-aware ensemble map"></svg><div class="cap" id="map-cap"></div></section>
<section id="pw-box" hidden><h2>Pathway A → B</h2><div class="pw-text" id="pw-text"></div><svg id="pw-prof" role="img" aria-label="Energy along the refolding path"></svg><div class="row"><button id="pw-play">Play the mechanism</button></div><svg id="pw-phi" role="img" aria-label="phi-value scan"></svg><div class="cap" id="pw-cap"></div></section>
<section><h2>Secondary structure</h2><canvas id="ss"></canvas><div class="cap" id="ss-cap"></div></section>
<section><h2>Sequence</h2><div class="seq" id="sequence"></div></section>
<section id="methods-box"><h2>Methods</h2><div id="methods"></div></section>
<div id="made"></div>
</aside></main>
<div id="tip" hidden></div><div id="err" hidden></div>"""


def build_html(doc, view=None, title=None, methods=None, analysis=None):
    static = P.static_dir()
    if 'landscape' not in doc.get('ensemble', {}):          # documents from before RCV 0.7
        try:
            from .landscape import layout_doc
            doc = {**doc, 'ensemble': {**doc['ensemble'], 'landscape': layout_doc(doc)}}
        except Exception:
            pass
    data = json.dumps({'doc': doc, 'view': view or {}, 'title': title or doc.get('name'), 'methods': methods,
                       'analysis': analysis or {}, 'version': __version__}, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
    icon = (static / 'img' / 'icon.svg').read_text(encoding='utf-8')
    mods = []
    for rel in MODULES:
        b64 = base64.b64encode((static / rel).read_bytes()).decode('ascii')
        mods.append(f'<script type="text/x-rcv-module" data-path="{rel}">{b64}</script>')
    page_title = html.escape(f'{title or doc.get("name", "RNA ensemble")} · RCV')
    return ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
            f'<meta name="generator" content="RCV {__version__} (RNA Conformation Viewer)">\n'
            f'<title>{page_title}</title>\n<style>{STYLE}</style>\n</head>\n<body>\n'
            + BODY.format(icon=icon) + '\n'
            f'<script type="application/json" id="rcv-data">{data}</script>\n'
            + '\n'.join(mods) + f'\n<script>{BOOT}</script>\n'
            '<!-- three.js (embedded): MIT License, Copyright 2010-2023 three.js authors. RCV: GPL-3.0. -->\n'
            '</body>\n</html>\n')
