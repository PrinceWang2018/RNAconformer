"""Draw the MSIX / Microsoft Store visual assets of RCV from packaging/icons/rcv.png (1024 px).

    python packaging/windows/msix/make_assets.py        # needs Pillow; writes packaging/windows/msix/Assets/

The PNGs are kept in the repository, so building the MSIX package (build_msix.ps1) needs no Pillow. Run this
again only when the icon changes. Sizes follow Microsoft's app icon guidance: tiles at scale 100/200/400 (the
package's resources.pri picks the right one), the app list / taskbar icon at target sizes 16-256 both plated and
unplated (`altform-unplated`: drawn without a coloured plate behind it), the Store logo and a splash screen.
"""
from __future__ import annotations

from pathlib import Path

from PIL import Image

HERE = Path(__file__).resolve().parent
ICON = HERE.parents[1] / 'icons' / 'rcv.png'
OUT = HERE / 'Assets'

# name -> base size (width, height) at scale 100; the icon fills `fill` of the shorter side
TILES = {
    'Square44x44Logo': ((44, 44), 1.0),
    'Square71x71Logo': ((71, 71), 0.62),
    'Square150x150Logo': ((150, 150), 0.62),
    'Square310x310Logo': ((310, 310), 0.62),
    'Wide310x150Logo': ((310, 150), 0.62),
    'StoreLogo': ((50, 50), 1.0),
    'SplashScreen': ((620, 300), 0.6),
    'FileIcon': ((44, 44), 1.0),
}
SCALES = (100, 200, 400)
TARGET_SIZES = (16, 20, 24, 30, 32, 36, 40, 48, 60, 64, 72, 80, 96, 256)


def render(icon, size, fill):
    w, h = size
    side = max(1, round(min(w, h) * fill))
    canvas = Image.new('RGBA', (w, h), (0, 0, 0, 0))
    art = icon.resize((side, side), Image.LANCZOS)
    canvas.alpha_composite(art, ((w - side) // 2, (h - side) // 2))
    return canvas


def main():
    icon = Image.open(ICON).convert('RGBA')
    icon = icon.crop(icon.getbbox())                   # the artwork without the transparent margin
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob('*.png'):
        old.unlink()
    n = 0
    for name, ((w, h), fill) in TILES.items():
        for scale in SCALES:
            size = (round(w * scale / 100), round(h * scale / 100))
            render(icon, size, fill).save(OUT / f'{name}.scale-{scale}.png', optimize=True)
            n += 1
    for name in ('Square44x44Logo', 'FileIcon'):
        for t in TARGET_SIZES:
            img = render(icon, (t, t), 1.0)
            img.save(OUT / f'{name}.targetsize-{t}.png', optimize=True)
            img.save(OUT / f'{name}.targetsize-{t}_altform-unplated.png', optimize=True)
            if name == 'Square44x44Logo':
                img.save(OUT / f'{name}.targetsize-{t}_altform-lightunplated.png', optimize=True)
            n += 3 if name == 'Square44x44Logo' else 2
    print(f'{n} images written to {OUT}')


if __name__ == '__main__':
    main()
