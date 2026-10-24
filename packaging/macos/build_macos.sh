#!/bin/bash
# Build RCV.app and a .dmg on macOS (run on the target architecture: Apple Silicon or Intel).
#   Requirements: Python 3.12 from python.org (universal2) or Homebrew; Xcode command-line tools.
#   Usage (from the repository root):  bash packaging/macos/build_macos.sh
#   (or let GitHub build it: .github/workflows/installers.yml builds Apple Silicon and Intel disk images)
# Optional signing/notarisation: export CODESIGN_ID="Developer ID Application: ..." and NOTARY_PROFILE=<keychain profile>.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PYTHON:-python3.12}
$PY -m venv .venv-build
. .venv-build/bin/activate
pip install --upgrade pip
pip install -r packaging/requirements-desktop.txt
pip install --no-deps -e .
python -c "import RNA, webview; from rcv import __version__; print('RCV', __version__, 'ViennaRNA', RNA.__version__)"
python -m PyInstaller packaging/rcv.spec --noconfirm --distpath dist --workpath build/pyinstaller
APP=dist/RCV.app
"$APP/Contents/MacOS/RCV" engine ss --sequence GGCACUUCGGUGCC
"$APP/Contents/MacOS/RCV" convert examples/viewer_demo/HIV1_TAR_NMR_1ANR.pdb -o build/selftest.rcv
"$APP/Contents/MacOS/RCV" info build/selftest.rcv
VERSION=$(python -c "from rcv import __version__; print(__version__)")
ARCH=$(uname -m)
if [ -n "${CODESIGN_ID:-}" ]; then
  codesign --deep --force --options runtime --timestamp --sign "$CODESIGN_ID" "$APP"
fi
# Disk image: RCV.app, a link to Applications (drag to install) and the user guide as web pages
python packaging/make_docs.py build/docs_html
STAGE=build/dmg
rm -rf "$STAGE" && mkdir -p "$STAGE"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
cp -R build/docs_html "$STAGE/RCV User Guide"
cp LICENSE "$STAGE/LICENSE.txt" && cp packaging/THIRD-PARTY-NOTICES.txt "$STAGE/"
mkdir -p dist/installer
DMG=dist/installer/RCV-$VERSION-macOS-$ARCH.dmg
rm -f "$DMG"
hdiutil create -volname "RCV $VERSION" -srcfolder "$STAGE" -ov -format UDZO "$DMG"
if [ -n "${CODESIGN_ID:-}" ] && [ -n "${NOTARY_PROFILE:-}" ]; then
  codesign --sign "$CODESIGN_ID" "$DMG"
  xcrun notarytool submit "$DMG" --keychain-profile "$NOTARY_PROFILE" --wait
  xcrun stapler staple "$DMG"
fi
echo "built $DMG"
