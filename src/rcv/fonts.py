"""Figure font for PDF export: Liberation Sans 2.1.5 (SIL Open Font License 1.1; metric-compatible with Arial and
Helvetica), bundled with RCV in static/fonts/ together with its licence. Exported PDFs embed subsets of it (fontTools),
as journals require embedded fonts. Without fontTools, PDFs fall back to the standard Helvetica font.
"""
from __future__ import annotations

import base64
import io
import threading
from pathlib import Path

_CACHE = {}
_LOCK = threading.Lock()


def font_path(bold=False):
    from . import paths as P
    p = P.static_dir() / 'fonts' / ('LiberationSans-Bold.ttf' if bold else 'LiberationSans-Regular.ttf')
    return p if p.is_file() else None


def _load(bold):
    with _LOCK:
        if bold not in _CACHE:
            path = font_path(bold)
            font = raw = None
            if path:
                try:
                    from fontTools.ttLib import TTFont
                    raw = path.read_bytes()
                    font = TTFont(io.BytesIO(raw))          # read-only: metrics and character map
                except Exception:
                    font = raw = None
            _CACHE[bold] = (path, font, raw)
        return _CACHE[bold]


def _winansi(code):
    try:
        ch = bytes([code]).decode('cp1252')
    except UnicodeDecodeError:
        return None
    return ch


def pdf_font(text, bold=False):
    """Subset of the figure font (regular or bold) for the given characters, with the metrics a PDF font dictionary
    needs. Returns None when fontTools is not available."""
    path, font, raw = _load(bold)
    if font is None:
        return None
    from fontTools import subset
    cmap = font.getBestCmap()
    upm = font['head'].unitsPerEm
    hmtx = font['hmtx']
    scale = 1000.0 / upm

    def width(u):
        g = cmap.get(u)
        return round(hmtx[g][0] * scale) if g else 0
    win = {c: _winansi(c) for c in range(32, 256)}
    win_chars = {ch for ch in win.values() if ch}
    used = set(text) | set(chr(c) for c in range(32, 127))
    extra = sorted(ch for ch in set(text) if ch not in win_chars and ord(ch) in cmap and ord(ch) > 32)[:90]
    unicodes = {ord(ch) for ch in used if ord(ch) in cmap} | {ord(ch) for ch in extra}
    opts = subset.Options()
    opts.layout_features = []
    opts.hinting = False
    opts.desubroutinize = True
    opts.drop_tables = list(opts.drop_tables) + ['FFTM']      # FontForge timestamp: not needed, avoids a warning
    opts.notdef_outline = True
    opts.name_IDs = ['*']
    opts.name_languages = ['*']
    opts.drop_tables += ['meta', 'DSIG', 'GSUB', 'GPOS', 'GDEF', 'kern', 'hdmx', 'VDMX', 'LTSH', 'JSTF', 'PCLT', 'gasp']
    from fontTools.ttLib import TTFont
    sub = TTFont(io.BytesIO(raw))
    s = subset.Subsetter(opts)
    s.populate(unicodes=unicodes)
    s.subset(sub)
    out = io.BytesIO()
    sub.save(out)
    head, os2, post, hhea = font['head'], font['OS/2'], font['post'], font['hhea']
    return {
        'name': 'LiberationSans-Bold' if bold else 'LiberationSans', 'ttf': base64.b64encode(out.getvalue()).decode('ascii'),
        'bbox': [round(v * scale) for v in (head.xMin, head.yMin, head.xMax, head.yMax)],
        'ascent': round(hhea.ascent * scale), 'descent': round(hhea.descent * scale),
        'capHeight': round(getattr(os2, 'sCapHeight', 0.716 * upm) * scale), 'italicAngle': float(post.italicAngle),
        'stemV': 120 if bold else 88, 'widths': [width(ord(win[c])) if win[c] else 0 for c in range(32, 256)],
        'extra': [{'char': ch, 'width': width(ord(ch)), 'glyph': f'uni{ord(ch):04X}'} for ch in extra],
    }
