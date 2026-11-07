"""Four-panel comparison figure: RCV next to deep-learning RNA structure predictors.

  a  download and disk footprint (documented sizes; weights measured on the download servers),
  b  RCV's time per CASP16 target on a laptop CPU, without GPU (results/per_target.csv, results/timing.json),
  c  RNAs with several experimental conformations: how well the five models cover them (results/ensembles.json),
  d  what each tool offers for an RNA that moves (from each tool's documentation).

Writes results/rcv_comparison.svg (plain SVG, light and dark colours, hover titles). For print, add
`--pdf --png`: the SVG is rendered by a headless Edge/Chrome with the bundled Liberation Sans font.

    python benchmark/casp16/figure.py [--pdf] [--png] [--installer-mb 65,71]
"""
from __future__ import annotations

import argparse
import base64
import csv
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
FONTS = HERE.parents[1] / 'src' / 'rcv' / 'static' / 'fonts'

# ---------------------------------------------------------------------------------------------- data
# Sizes in MB. 'model': program/model weights (download); 'full': with the sequence databases each tool needs for
# its full-accuracy (MSA) mode, as stated in its documentation. Python/CUDA/PyTorch/JAX environments are not
# counted for the other tools, so their numbers are lower bounds. Checked 2026-11-07.
FOOTPRINT = [
    {'tool': 'RCV 0.7.3', 'kind': 'physics + templates · CPU', 'model': None, 'full': None,
     'note': 'complete installer: Python, engine, ViennaRNA, viewer and the 21.5 MB structure library'},
    {'tool': 'trRosettaRNA2', 'kind': 'deep learning', 'model': 108, 'full': 108 + 32_000, 'full_label': '32 GB',
     'note': 'weights params_trRNA2.tar.bz2 108 MB (measured); RNAcentral database ~32 GB uncompressed (README, optional)',
     'src': 'github.com/YangLab-SDU/trRosettaRNA2'},
    {'tool': 'RhoFold+', 'kind': 'deep learning', 'model': 508, 'full': 508 + 900_000, 'full_label': '~900 GB',
     'note': 'checkpoint 508 MB (measured); RNAcentral, Rfam and nt "totaling about 900GB" (README)',
     'src': 'github.com/ml4bio/RhoFold'},
    {'tool': 'RoseTTAFold2NA', 'kind': 'deep learning', 'model': 962, 'full_label': '≥ 482 GB',
     'full': 962 + 46_000 + 272_000 + 300 + 12_000 + 151_000,
     'note': 'weights RF2NA_apr23.tgz 962 MB (measured); UniRef30 46G, BFD 272G, Rfam 300M, RNAcentral 12G, nt 151G '
             '(README; PDB100 size not stated, so a lower bound)', 'src': 'github.com/uw-ipd/RoseTTAFold2NA'},
    {'tool': 'AlphaFold 3', 'kind': 'deep learning · GPU', 'model': 1021, 'full': 1021 + 630_000, 'full_label': '631 GB',
     'note': 'af3.bin.zst 1.02 GB (measured); genetic databases 630 GB unpacked, 252 GB download (installation.md)',
     'src': 'github.com/google-deepmind/alphafold3'},
]

# what each tool gives for an RNA with several states (documented outputs; '–' = not part of the documented output)
TOOLS = ['RCV 0.7.3', 'AlphaFold 3', 'trRosettaRNA2', 'RhoFold+', 'RoseTTAFold2NA']
HEAD = {'RCV 0.7.3': ('RCV', '0.7.3'), 'AlphaFold 3': ('AlphaFold', '3'), 'trRosettaRNA2': ('trRosetta-', 'RNA2'),
        'RhoFold+': ('RhoFold+', ''), 'RoseTTAFold2NA': ('RoseTTA-', 'Fold2NA')}
CAPABILITIES = [
    ('Conformations per run', ['ensemble¹', '5 per seed²', '1³', '1', '1']),
    ('Alternative 2D states, populations', ['✓', '–', '–', '–', '–']),
    ('Temperature, salt, ligand', ['✓', '–', '–', '–', '–']),
    ('Switching paths, barriers, rates', ['✓', '–', '–', '–', '–']),
    ('Fit to SHAPE/DMS probing data', ['✓', '–', '–', '–', '–']),
    ('Interactive ensemble viewer', ['✓', '–⁴', '–', '–', '–']),
    ('Runs without a GPU', ['✓', '–', 'not stated', '✓', 'not stated']),
    ('Sequence databases needed', ['none', '630 GB', '32 GB⁵', '~900 GB⁵', '≥ 481 GB']),
]
CAP_NOTES = ['¹ structurally distinct conformations from different 2D/template hypotheses (1–30, default 8)',
             '² diffusion samples of one prediction; more seeds give more samples',
             '³ the trRosettaRNA2 paper describes conformers; the default output is one model',
             '⁴ AlphaFold 3 results open in RCV (0.7.3), so RCV\'s analyses apply to them',
             '⁵ optional: single-sequence mode without databases is less accurate']

# ---------------------------------------------------------------------------------------------- style
W, H = 1240, 1010
CSS = """
.bg{fill:#fcfcfb}.t1{fill:#0b0b0b}.t2{fill:#52514e}.t3{fill:#898781}
.grid{stroke:#e1e0d9;stroke-width:1}.axis{stroke:#c3c2b7;stroke-width:1}.span{stroke:#c3c2b7;stroke-width:2}
.s1{fill:#2a78d6}.s2{fill:#eb6834}.s3{fill:#1baf7a}.oth{fill:#898781}.ring{stroke:#fcfcfb;stroke-width:2}
.hollow{fill:#fcfcfb;stroke:#2a78d6;stroke-width:2.2}.wash{fill:#2a78d6;fill-opacity:.08}.sep{stroke:#e1e0d9;stroke-width:1}
@media (prefers-color-scheme: dark){
.bg{fill:#1a1a19}.t1{fill:#f0efec}.t2{fill:#c3c2b7}.t3{fill:#898781}
.grid{stroke:#2c2c2a}.axis{stroke:#383835}.span{stroke:#4a4a45}
.s1{fill:#3987e5}.s2{fill:#d95926}.s3{fill:#199e70}.ring{stroke:#1a1a19}
.hollow{fill:#1a1a19;stroke:#3987e5}.wash{fill:#3987e5;fill-opacity:.12}.sep{stroke:#2c2c2a}}
"""


def esc(s):
    return str(s).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def T(x, y, s, cls='t1', size=12, anchor='start', weight=None, title=None):
    w = f' font-weight="{weight}"' if weight else ''
    t = f'<title>{esc(title)}</title>' if title else ''
    return f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}"{w}>{t}{esc(s)}</text>'


def dot(x, y, cls, r=5.5, title='', shape='circle'):
    t = f'<title>{esc(title)}</title>' if title else ''
    if shape == 'diamond':
        d = r * 1.25
        return f'<path class="{cls}" d="M{x:.1f} {y - d:.1f}L{x + d:.1f} {y:.1f}L{x:.1f} {y + d:.1f}L{x - d:.1f} {y:.1f}Z">{t}</path>'
    if shape == 'square':
        return f'<rect class="{cls}" x="{x - r:.1f}" y="{y - r:.1f}" width="{2 * r:.1f}" height="{2 * r:.1f}" rx="1.5">{t}</rect>'
    if shape == 'triangle':
        d = r * 1.25
        return f'<path class="{cls}" d="M{x:.1f} {y - d:.1f}L{x + d:.1f} {y + d * 0.8:.1f}L{x - d:.1f} {y + d * 0.8:.1f}Z">{t}</path>'
    return f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="{r}">{t}</circle>'


def fmt_mb(mb):
    if mb >= 1_000_000:
        return f'{mb / 1_000_000:.1f} TB'.replace('.0 TB', ' TB')
    if mb >= 1000:
        v = mb / 1000
        return f'{v:.0f} GB' if v >= 10 else f'{v:.1f} GB'
    return f'{mb:.0f} MB'


def panel_title(o, x, y, letter, title, sub):
    o.append(T(x, y, letter, 't1', 17, weight='bold'))
    o.append(T(x + 20, y, title, 't1', 14.5, weight='bold'))
    o.append(T(x + 20, y + 18, sub, 't2', 11.5))


# ---------------------------------------------------------------------------------------------- panels
def panel_size(o, x, y, w, installer):
    panel_title(o, x, y, 'a', 'Download and disk space', 'Log scale. Other tools: Python/GPU software not counted.')
    left, right = x + 150, x + w - 20
    lo, hi = 1, 6                                         # 10 MB .. 1 TB
    X = lambda mb: left + (right - left) * (math.log10(mb) - lo) / (hi - lo)
    top, rowh = y + 62, 50
    bottom = top + rowh * len(FOOTPRINT)
    for e, lab in zip(range(lo, hi + 1), ['10 MB', '100 MB', '1 GB', '10 GB', '100 GB', '1 TB']):
        o.append(f'<line class="grid" x1="{X(10 ** e):.1f}" x2="{X(10 ** e):.1f}" y1="{top - 6}" y2="{bottom}"/>')
        o.append(T(X(10 ** e), bottom + 16, lab, 't3', 11, 'middle'))
    # legend
    lx, ly = x + 20, bottom + 40
    for cls, shape, lab in [('s1 ring', 'circle', 'RCV: complete installer'), ('s3 ring', 'square', 'model weights'),
                            ('s2 ring', 'diamond', 'with sequence databases (full-accuracy mode)')]:
        o.append(dot(lx + 6, ly - 4, cls, 5, shape=shape))
        o.append(T(lx + 16, ly, lab, 't2', 11.5))
        lx += 26 + 6.2 * len(lab)
    for k, f in enumerate(FOOTPRINT):
        cy = top + rowh * k + rowh / 2
        o.append(T(left - 14, cy - 1, f['tool'], 't1', 12.5, 'end', weight='bold' if k == 0 else None))
        o.append(T(left - 14, cy + 13, f['kind'], 't3', 10.5, 'end'))
        if f['model'] is None:                            # RCV: one complete package
            lo_mb, hi_mb = installer
            mid = (lo_mb + hi_mb) / 2
            if hi_mb > lo_mb:
                o.append(f'<line class="span" x1="{X(lo_mb):.1f}" x2="{X(hi_mb):.1f}" y1="{cy:.1f}" y2="{cy:.1f}"/>')
            o.append(dot(X(mid), cy, 's1 ring', 6.5, f"RCV 0.7.3 installer: {lo_mb:.0f}–{hi_mb:.0f} MB, {f['note']}"))
            lab = f'{lo_mb:.0f}–{hi_mb:.0f} MB, everything included' if hi_mb > lo_mb else f'{mid:.0f} MB, everything included'
            o.append(T(X(hi_mb) + 12, cy + 4, lab, 't1', 11.5))
            continue
        o.append(f'<line class="span" x1="{X(f["model"]):.1f}" x2="{X(f["full"]):.1f}" y1="{cy:.1f}" y2="{cy:.1f}"/>')
        o.append(dot(X(f['model']), cy, 's3 ring', 5.5, f"{f['tool']}: {f['note']}", 'square'))
        o.append(dot(X(f['full']), cy, 's2 ring', 5.5, f"{f['tool']}: {f['note']}", 'diamond'))
        o.append(T(X(f['model']), cy - 11, fmt_mb(f['model']), 't3', 10.5, 'middle'))
        o.append(T(X(f['full']), cy - 11, f['full_label'], 't2', 10.5, 'middle'))


def panel_speed(o, x, y, w, rows, timing):
    panel_title(o, x, y, 'b', 'Time per RNA on a laptop CPU, no GPU',
                f"RCV 0.7.3, 5 models from 40 trajectories, CASP16 targets. {timing.get('cpu_short', 'Laptop CPU')}.")
    left, right, top, bottom = x + 56, x + w - 24, y + 62, y + 300
    pts = [(int(r['length']), float(r['rcv_wall_time_s']) / 60, r['target']) for r in rows]
    for t in timing.get('extra16', []):
        pts.append((t['length'], t['wall_time_s'] / 60, t['target']))
    four = [(t['length'], t['wall_time_s'] / 60, t['target']) for t in timing.get('runs', [])]
    top_v = max([p[1] for p in pts + four] + [1])
    step = 1 if top_v <= 7 else 2 if top_v <= 16 else 5
    ymax = int(math.ceil(top_v / step) * step)
    xmin, xmax = 50, 260
    X = lambda v: left + (right - left) * (v - xmin) / (xmax - xmin)
    Y = lambda v: bottom - (bottom - top) * v / ymax
    for v in range(0, ymax + 1, step):
        o.append(f'<line class="grid" x1="{left}" x2="{right}" y1="{Y(v):.1f}" y2="{Y(v):.1f}"/>')
        o.append(T(left - 8, Y(v) + 4, str(v), 't3', 11, 'end'))
    for v in range(50, 261, 50):
        o.append(T(X(v), bottom + 16, str(v), 't3', 11, 'middle'))
    o.append(f'<line class="axis" x1="{left}" x2="{right}" y1="{bottom}" y2="{bottom}"/>')
    o.append(T((left + right) / 2, bottom + 34, 'RNA length (nucleotides)', 't2', 11.5, 'middle'))
    o.append(f'<text class="t2" font-size="11.5" text-anchor="middle" transform="translate({x + 16},{(top + bottom) / 2}) '
             f'rotate(-90)">Minutes per RNA</text>')
    for L, m, tid in pts:
        o.append(dot(X(L), Y(m), 's1 ring', 5, f'{tid}, {L} nt: {m:.1f} min on 16 CPU cores'))
    for L, m, tid in four:
        o.append(dot(X(L), Y(m), 'hollow', 4.6, f'{tid}, {L} nt: {m:.1f} min on {timing.get("cores", 4)} CPU cores'))
    med = sorted(p[1] for p in pts)[len(pts) // 2]
    lx, ly = X(160), Y(ymax * 0.31)
    o.append(dot(lx + 5, ly, 's1 ring', 5))
    o.append(T(lx + 15, ly + 4, f'16 cores (median {med:.1f} min)', 't2', 11.5))
    if four:
        o.append(dot(lx + 5, ly + 20, 'hollow', 4.6))
        o.append(T(lx + 15, ly + 24, f'{timing.get("cores", 4)} cores, same settings', 't2', 11.5))
    notes = ['Hardware the others need (documentation):',
             'AlphaFold 3 — NVIDIA GPU, compute capability ≥ 8.0',
             '(A100/H100 tested), ≥ 64 GB RAM; genetic search on CPU.']
    if timing.get('peak_memory_gb'):
        notes.insert(0, f'RCV peak memory: {timing["peak_memory_gb"]:.1f} GB (247 nt, all {timing.get("peak_workers", 16)} '
                        'workers together).')
        notes.insert(1, '')
    for k, s in enumerate(notes):
        o.append(T(x + 20, bottom + 62 + 16 * k, s, 't2', 11.5))


def panel_states(o, x, y, w, ens):
    panel_title(o, x, y, 'c', 'RNAs solved in several conformations',
                'CASP16: best TM-score of the 5 models for each experimental conformation, averaged.')
    left, right = x + 170, x + w - 24
    rows = ens.get('multi', [])
    top, rowh = y + 92, 64
    bottom = top + rowh * len(rows)
    xmin, xmax = 0.2, 0.8
    X = lambda v: left + (right - left) * (v - xmin) / (xmax - xmin)
    for v in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
        o.append(f'<line class="grid" x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top - 6}" y2="{bottom}"/>')
        o.append(T(X(v), bottom + 16, f'{v:.1f}', 't3', 11, 'middle'))
    o.append(T((left + right) / 2, bottom + 34, 'Coverage of the experimental conformations (TM-score)', 't2', 11.5, 'middle'))
    lx, ly = x + 20, y + 62
    for cls, shape, lab in [('s1 ring', 'circle', 'RCV 0.7.3'), ('s2 ring', 'diamond', 'AlphaFold 3'),
                            ('oth ring', 'triangle', 'trRosettaRNA, Vfold, GuangzhouRNA (CASP16)')]:
        o.append(dot(lx + 6, ly - 4, cls, 5, shape=shape))
        o.append(T(lx + 16, ly, lab, 't2', 11.5))
        lx += 28 + 6.3 * len(lab)
    names = {'052': 'trRosettaRNA', '481': 'Vfold', '183': 'GuangzhouRNA'}
    for k, r in enumerate(rows):
        cy = top + rowh * k + rowh / 2
        o.append(T(left - 14, cy - 2, f"{r['target']} · {r['length']} nt", 't1', 12.5, 'end'))
        o.append(T(left - 14, cy + 12, r['family'], 't3', 10.5, 'end'))
        o.append(T(left - 14, cy + 25, f"{r['conformations']} conformations", 't3', 10.5, 'end'))
        m = r['methods']
        vals = [v['coverage'] for v in m.values()]
        o.append(f'<line class="span" x1="{X(min(vals)):.1f}" x2="{X(max(vals)):.1f}" y1="{cy:.1f}" y2="{cy:.1f}"/>')
        for g, nm in names.items():
            if g in m:
                o.append(dot(X(m[g]['coverage']), cy, 'oth ring', 4.5, f"{r['target']} {nm}: coverage {m[g]['coverage']:.2f}, "
                             f"spread of the 5 models {m[g]['spread']:.2f}", 'triangle'))
        if '304' in m:
            o.append(dot(X(m['304']['coverage']), cy, 's2 ring', 5.5, f"{r['target']} AlphaFold 3: coverage "
                         f"{m['304']['coverage']:.2f}, spread {m['304']['spread']:.2f}", 'diamond'))
        if 'rcv' in m:
            o.append(dot(X(m['rcv']['coverage']), cy, 's1 ring', 6, f"{r['target']} RCV 0.7.3: coverage "
                         f"{m['rcv']['coverage']:.2f}, spread {m['rcv']['spread']:.2f}"))
            o.append(T(X(m['rcv']['coverage']), cy - 11, f"{m['rcv']['coverage']:.2f}", 't2', 10.5, 'middle'))


def panel_capabilities(o, x, y, w):
    panel_title(o, x, y, 'd', 'For an RNA that moves', 'What each tool provides, from its documentation.')
    top, rowh = y + 88, 30
    c0 = x + 196
    colw = (w - 196 - 6) / len(TOOLS)
    o.append(f'<rect class="wash" x="{c0:.1f}" y="{top - 38}" width="{colw:.1f}" height="{rowh * len(CAPABILITIES) + 42}" rx="6"/>')
    for j, t in enumerate(TOOLS):
        l1, l2 = HEAD[t]
        o.append(T(c0 + colw * (j + 0.5), top - 20, l1, 't1', 11.5, 'middle', weight='bold' if j == 0 else None))
        if l2:
            o.append(T(c0 + colw * (j + 0.5), top - 6, l2, 't1', 11.5, 'middle', weight='bold' if j == 0 else None))
    for i, (cap, vals) in enumerate(CAPABILITIES):
        cy = top + rowh * i + 20
        o.append(f'<line class="sep" x1="{x + 20}" x2="{x + w - 10}" y1="{cy - 18:.1f}" y2="{cy - 18:.1f}"/>')
        o.append(T(x + 20, cy, cap, 't1', 11.5))
        for j, v in enumerate(vals):
            o.append(T(c0 + colw * (j + 0.5), cy, v, 't1' if v.startswith('✓') or j == 0 else 't2', 12 if v == '✓' else 11,
                       'middle', weight='bold' if v == '✓' else None))
    for k, s in enumerate(CAP_NOTES):
        o.append(T(x + 20, top + rowh * len(CAPABILITIES) + 20 + 15 * k, s, 't3', 10.5))


# ---------------------------------------------------------------------------------------------- figure
def build(installer):
    rows = list(csv.DictReader((RES / 'per_target.csv').open()))
    ens = json.loads((RES / 'ensembles.json').read_text()) if (RES / 'ensembles.json').exists() else {}
    timing = json.loads((RES / 'timing.json').read_text()) if (RES / 'timing.json').exists() else {}
    if timing.get('cpu'):
        timing['cpu_short'] = f"{timing['cpu'].replace('(R)', '').replace('(TM)', '')}, {timing.get('logical_cpus', 16)} cores"
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         'font-family="Liberation Sans, Arial, Helvetica, sans-serif" role="img" '
         'aria-label="RCV compared with deep-learning RNA structure predictors: size, speed and RNAs with several conformations">',
         f'<style>{CSS}</style>', f'<rect class="bg" width="{W}" height="{H}" rx="8"/>',
         T(24, 36, 'RCV next to deep-learning RNA 3D structure predictors', 't1', 19, weight='bold'),
         T(24, 58, 'Size, hardware and RNAs that move. RCV: physics + experimental templates, no deep learning, no GPU.',
           't2', 12.5)]
    panel_size(o, 24, 100, 590, installer)
    panel_speed(o, 636, 100, 580, rows, timing)
    panel_states(o, 24, 560, 590, ens)
    panel_capabilities(o, 636, 560, 580)
    o.append(T(24, H - 16, 'Sizes and requirements as documented by each tool (checked 7 Nov 2026); weight files measured on the '
               'download servers. Accuracy on the 15 CASP16 targets: benchmark/casp16/results.', 't3', 10.5))
    o.append('</svg>')
    return '\n'.join(o) + '\n'


def browser():
    for c in [os.environ.get('RCV_BROWSER', ''), '/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
              '/mnt/c/Program Files/Google/Chrome/Application/chrome.exe', shutil.which('chromium') or '',
              shutil.which('google-chrome') or '', shutil.which('microsoft-edge') or '']:
        if c and Path(c).exists():
            return c
    return None


def render(svg, pdf, png):
    """PDF (vector) and PNG (3x) via a headless Chromium browser, Liberation Sans embedded."""
    exe = browser()
    if not exe:
        raise SystemExit('no Edge/Chrome found for --pdf/--png (set RCV_BROWSER)')
    fonts = ''.join(f"@font-face{{font-family:'Liberation Sans';font-weight:{wt};src:url(data:font/ttf;base64,"
                    f"{base64.b64encode((FONTS / fn).read_bytes()).decode()})}}"
                    for fn, wt in [('LiberationSans-Regular.ttf', 400), ('LiberationSans-Bold.ttf', 700)])
    html = (f'<!DOCTYPE html><html><head><meta charset="utf-8"><style>{fonts}@page{{size:{W}px {H}px;margin:0}}'
            f'html,body{{margin:0;padding:0;background:#fcfcfb}}svg{{display:block}}</style></head><body>{svg}</body></html>')
    win = exe.startswith('/mnt/')
    tmpdir = Path(tempfile.mkdtemp(dir='/mnt/c/Windows/Temp' if win and Path('/mnt/c/Windows/Temp').exists() else None))
    page = tmpdir / 'figure.html'
    page.write_text(html, encoding='utf-8')
    path = lambda p: subprocess.run(['wslpath', '-w', str(p)], capture_output=True, text=True).stdout.strip() if win else str(p)
    url = 'file:///' + path(page).replace('\\', '/').lstrip('/')
    base = [exe, '--headless=new', '--disable-gpu', '--no-pdf-header-footer', '--hide-scrollbars']
    if pdf:
        subprocess.run(base + [f'--print-to-pdf={path(tmpdir / "figure.pdf")}', url], capture_output=True, timeout=120)
        shutil.copy(tmpdir / 'figure.pdf', pdf)
        print(pdf)
    if png:
        subprocess.run(base + [f'--window-size={W},{H}', '--force-device-scale-factor=3',
                               f'--screenshot={path(tmpdir / "figure.png")}', url], capture_output=True, timeout=120)
        shutil.copy(tmpdir / 'figure.png', png)
        print(png)
    shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--installer-mb', default='65,71', help='RCV installer size in MB, "low,high" over the platforms')
    ap.add_argument('--pdf', action='store_true')
    ap.add_argument('--png', action='store_true')
    a = ap.parse_args()
    lo, hi = (float(v) for v in (a.installer_mb.split(',') * 2)[:2])
    svg = build((lo, hi))
    out = RES / 'rcv_comparison.svg'
    out.write_text(svg, encoding='utf-8')
    print(out)
    if a.pdf or a.png:
        render(svg, RES / 'rcv_comparison.pdf' if a.pdf else None, RES / 'rcv_comparison.png' if a.png else None)


if __name__ == '__main__':
    main()
