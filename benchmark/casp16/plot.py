"""Dot plot of the CASP16 comparison (results/per_target.csv -> results/casp16_tm.svg).

One row per target; TM-score of model 1 for RCV (with and without templates), AlphaFold 3 and the best CASP16
group on that target.
Plain SVG (no plotting library) with light and dark colours, hover titles and marker shapes as a second cue.
"""
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

SERIES = [('rcv', 'RCV 0.7.3', 'circle'),
          ('rcvfree', 'RCV without templates', 'ring'),
          ('304', 'AlphaFold 3', 'diamond'),
          ('best', 'Best CASP16 group', 'triangle')]
CLS = {'rcv': 'rcv', 'rcvfree': 'rcvfree', '304': 'af3', 'best': 'best'}


def marker(kind, x, y, cls, title):
    t = f'<title>{title}</title>'
    if kind == 'circle':
        return f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="5.5">{t}</circle>'
    if kind == 'ring':
        return f'<circle class="{cls}" cx="{x:.1f}" cy="{y:.1f}" r="4.6">{t}</circle>'
    if kind == 'diamond':
        return f'<path class="{cls}" d="M{x:.1f} {y - 6.5:.1f} L{x + 6.5:.1f} {y:.1f} L{x:.1f} {y + 6.5:.1f} L{x - 6.5:.1f} {y:.1f} Z">{t}</path>'
    return f'<path class="{cls}" d="M{x:.1f} {y - 6.5:.1f} L{x + 6.5:.1f} {y + 5:.1f} L{x - 6.5:.1f} {y + 5:.1f} Z">{t}</path>'


def main():
    rows = list(csv.DictReader((RES / 'per_target.csv').open()))
    rows.sort(key=lambda r: -float(r['casp_best_top1_tm']))
    W, left, right, top, rowh = 760, 250, 30, 92, 26
    H = top + rowh * len(rows) + 56
    x0, x1 = left, W - right
    X = lambda v: x0 + (x1 - x0) * v
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
         'font-family="Liberation Sans, Arial, Helvetica, sans-serif" role="img" '
         'aria-label="TM-score of model 1 per CASP16 RNA target for RCV with and without templates, AlphaFold 3 and the best CASP16 group">',
         '<style>',
         '.bg{fill:#fcfcfb}.t1{fill:#0b0b0b}.t2{fill:#52514e}.t3{fill:#7a7974}.grid{stroke:#e6e5e0;stroke-width:1}',
         '.ref{stroke:#9a9993;stroke-width:1;stroke-dasharray:3 3}.span{stroke:#c9c8c2;stroke-width:2}',
         '.rcv{fill:#2a78d6}.af3{fill:#eb6834}.best{fill:#1baf7a}.m{stroke:#fcfcfb;stroke-width:2}',
         '.rcvfree{fill:#fcfcfb;stroke:#2a78d6;stroke-width:2.2}',
         '@media (prefers-color-scheme: dark){.bg{fill:#1a1a19}.t1{fill:#f0efec}.t2{fill:#c3c2b7}.t3{fill:#9a998f}',
         '.grid{stroke:#33332f}.ref{stroke:#6d6c66}.span{stroke:#4a4a45}.rcv{fill:#3987e5}.af3{fill:#d95926}',
         '.best{fill:#199e70}.m{stroke:#1a1a19}.rcvfree{fill:#1a1a19;stroke:#3987e5;stroke-width:2.2}}',
         '</style>',
         f'<rect class="bg" width="{W}" height="{H}" rx="8"/>',
         '<text class="t1" x="24" y="32" font-size="17" font-weight="bold">CASP16 RNA targets: TM-score of model 1</text>',
         '<text class="t2" x="24" y="52" font-size="12.5">Scored with US-align as in CASP16; higher is better. '
         'Without a template both RCV runs give the same model.</text>']
    lx = 24
    series = [x for x in SERIES if x[0] != 'rcvfree' or rows[0].get('rcvfree_top1_tm')]
    for key, label, kind in series:
        cls = CLS[key] + ('' if key == 'rcvfree' else ' m')
        o.append(marker(kind, lx + 6, 74, cls, label))
        o.append(f'<text class="t1" x="{lx + 18}" y="78" font-size="12.5">{label}</text>')
        lx += 30 + 7.1 * len(label)
    for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0):
        o.append(f'<line class="grid" x1="{X(v):.1f}" x2="{X(v):.1f}" y1="{top - 8}" y2="{top + rowh * len(rows)}"/>')
        o.append(f'<text class="t3" x="{X(v):.1f}" y="{top + rowh * len(rows) + 16}" font-size="11" text-anchor="middle">{v:.1f}</text>')
    o.append(f'<line class="ref" x1="{X(0.45):.1f}" x2="{X(0.45):.1f}" y1="{top - 8}" y2="{top + rowh * len(rows)}"/>')
    o.append(f'<text class="t3" x="{X(0.45) + 4:.1f}" y="{top + rowh * len(rows) + 32}" font-size="11">0.45 ≈ correct global fold</text>')
    o.append(f'<text class="t2" x="{(x0 + x1) / 2:.1f}" y="{top + rowh * len(rows) + 48}" font-size="12" text-anchor="middle">TM-score</text>')
    for k, r in enumerate(rows):
        y = top + rowh * k + rowh / 2
        name = f"{r['target']} · {r['family'].split(' (')[0]} · {r['length']} nt"
        o.append(f'<text class="t1" x="{left - 12}" y="{y + 4:.1f}" font-size="12" text-anchor="end">{name}</text>')
        vals = {'rcv': float(r['rcv_top1_tm']), '304': float(r['304_top1_tm']), 'best': float(r['casp_best_top1_tm'])}
        if r.get('rcvfree_top1_tm'):
            vals['rcvfree'] = float(r['rcvfree_top1_tm'])
        o.append(f'<line class="span" x1="{X(min(vals.values())):.1f}" x2="{X(max(vals.values())):.1f}" y1="{y:.1f}" y2="{y:.1f}"/>')
        tmpl = f" (template {r['rcv_template']})" if r.get('rcv_template') else ' (no template)'
        for key, label, kind in reversed(series):
            if key not in vals:
                continue
            cls = CLS[key] + ('' if key == 'rcvfree' else ' m')
            o.append(marker(kind, X(vals[key]), y, cls, f"{r['target']} {label}{tmpl if key == 'rcv' else ''}: TM-score {vals[key]:.2f}"))
    o.append('</svg>')
    (RES / 'casp16_tm.svg').write_text('\n'.join(o) + '\n', encoding='utf-8')
    print(RES / 'casp16_tm.svg')


if __name__ == '__main__':
    main()
