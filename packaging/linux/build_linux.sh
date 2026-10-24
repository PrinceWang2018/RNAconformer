#!/bin/bash
# Build the RCV Linux packages (x86_64): a PyInstaller one-folder build installed to /opt/rcv, as
#   dist/installer/RCV-<version>-Linux-amd64.deb    Debian, Ubuntu, Linux Mint, ...   (sudo apt install ./file.deb)
#   dist/installer/RCV-<version>-Linux-x86_64.rpm   Fedora, RHEL/Rocky/Alma, openSUSE (sudo dnf install ./file.rpm)
# Both add the `rcv` command, a menu entry (with "RCV with demo ensembles" and "RCV User Guide"), the icon, the
# .rcv file type and the user guide as web pages in /usr/share/doc/rcv/html.
#
#   Requirements: Python 3.12 with venv and a shared libpython (python.org, Docker "python" images or distro
#                 packages), dpkg-deb for the .deb and rpmbuild for the .rpm (on Debian/Ubuntu: apt install rpm).
#   Usage (from the repository root):  bash packaging/linux/build_linux.sh
#   (or let GitHub build it: .github/workflows/installers.yml builds on Debian 11 and test-installs the packages)
#
# The packages run on systems whose glibc is at least as new as the build machine's (Debian 11: glibc 2.31, so
# Ubuntu 20.04+, Debian 11+, Fedora 32+, RHEL/Rocky/Alma 9+). RCV opens in a Chrome, Chromium or Edge app window
# when one is installed, otherwise in the default browser; it quits when the window is closed.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PYTHON:-python3.12}
$PY -m venv .venv-build
. .venv-build/bin/activate
pip install --upgrade pip
pip install -r packaging/requirements-desktop.txt           # (pywebview is not used on Linux)
pip install --no-deps -e .
python -c "import RNA; from rcv import __version__; print('RCV', __version__, 'ViennaRNA', RNA.__version__)"
python -m PyInstaller packaging/rcv.spec --noconfirm --distpath dist --workpath build/pyinstaller

# Self-test of the build: engine, command-line tools and the interface server
APP=dist/RCV
mkdir -p build
"$APP/RCV" --version
"$APP/RCV" engine ss --sequence GGCACUUCGGUGCC
"$APP/RCV" convert examples/viewer_demo/HIV1_TAR_NMR_1ANR.pdb -o build/selftest.rcv
"$APP/RCV" info build/selftest.rcv
"$APP/RCV" --shell none --port 18765 --projects-dir build/selftest-projects > build/selftest-server.log 2>&1 &
PID=$!
for _ in $(seq 120); do grep -q 'running at' build/selftest-server.log && break; sleep 0.5; done
cat build/selftest-server.log
PORT=$(grep -o -m1 'localhost:[0-9]*' build/selftest-server.log | cut -d: -f2)
python -c "import sys, urllib.request; r = urllib.request.urlopen('http://localhost:$PORT/', timeout=30); print('server', r.status); sys.exit(r.status != 200)"
kill "$PID"; wait "$PID" || true

VERSION=$(python -c "from rcv import __version__; print(__version__)")
GLIBC=$(getconf GNU_LIBC_VERSION | awk '{print $2}')      # (no `| head`: with pipefail it can fail on SIGPIPE)
python packaging/make_docs.py build/docs_html

# File tree shared by both packages
STAGE=build/linux/root
rm -rf build/linux
mkdir -p "$STAGE/opt" "$STAGE/usr/bin" "$STAGE/usr/share/applications" "$STAGE/usr/share/mime/packages" \
         "$STAGE/usr/share/doc/rcv"
cp -a "$APP" "$STAGE/opt/rcv"
ln -s /opt/rcv/RCV "$STAGE/usr/bin/rcv"
install -m 644 packaging/linux/rcv.desktop "$STAGE/usr/share/applications/rcv.desktop"
install -m 644 packaging/linux/rcv-mime.xml "$STAGE/usr/share/mime/packages/rcv.xml"
for d in packaging/linux/icons/*/; do
  install -D -m 644 "$d/rcv.png" "$STAGE/usr/share/icons/hicolor/$(basename "$d")/apps/rcv.png"
done
cp -a build/docs_html "$STAGE/usr/share/doc/rcv/html"
install -m 644 LICENSE packaging/THIRD-PARTY-NOTICES.txt "$STAGE/usr/share/doc/rcv/"
cat > "$STAGE/usr/share/doc/rcv/copyright" <<EOF
RCV - RNA Conformation Viewer $VERSION
Copyright (C) 2026 Zixiang Wang <princewang2018@gmail.com>
Website: https://www.rnaconfviewer.com/   Source: https://github.com/PrinceWang2018/RNAconformer

RCV is free software under the GNU General Public License version 3 (see LICENSE in this folder).
/opt/rcv also contains Python, NumPy, SciPy, ViennaRNA and other software under their own licences; they are
listed, with the ViennaRNA licence and credit, in THIRD-PARTY-NOTICES.txt in this folder.
EOF
chmod -R u+rwX,go+rX,go-w "$STAGE"
find "$STAGE/usr/share" -type f -exec chmod 644 {} +
# Debugging information is not needed to run RCV (ViennaRNA's and Python's libraries shrink from ~85 to ~18 MB).
# Only the libraries: the RCV program itself carries the app's archive at its end, which strip would remove.
find "$STAGE/opt/rcv/_internal" -type f -name '*.so*' -exec sh -c 'for f; do strip --strip-debug "$f" 2>/dev/null || true; done' sh {} +
"$STAGE/opt/rcv/RCV" engine ss --sequence GGCACUUCGGUGCC     # still works after stripping
mkdir -p dist/installer

# .deb
DEBROOT=build/linux/deb
cp -a "$STAGE" "$DEBROOT"
mkdir -p "$DEBROOT/DEBIAN"
cat > "$DEBROOT/DEBIAN/control" <<EOF
Package: rcv
Version: $VERSION
Architecture: amd64
Maintainer: Zixiang Wang <princewang2018@gmail.com>
Installed-Size: $(du -sk "$STAGE" | cut -f1)
Depends: libc6 (>= $GLIBC), xdg-utils
Suggests: chromium | google-chrome-stable | microsoft-edge-stable
Section: science
Priority: optional
Homepage: https://www.rnaconfviewer.com/
Description: RNA Conformation Viewer: predict, view and compare RNA ensembles
 RCV predicts RNA conformational ensembles from a sequence (ViennaRNA secondary
 structures, 3D assembly and refinement) and shows them in 3D and 2D: morphs
 between conformations, plots, palettes, surfaces, annotations, figure and video
 export and plugins. Includes the demo ensembles and the user guide.
EOF
DEB=dist/installer/RCV-$VERSION-Linux-amd64.deb
dpkg-deb --root-owner-group -Zxz --build "$DEBROOT" "$DEB"
echo "built $DEB"

# .rpm
if command -v rpmbuild >/dev/null; then
  TOP=$(pwd)/build/linux/rpmbuild
  rpmbuild -bb packaging/linux/rcv-rpm.spec --define "_topdir $TOP" --define "rcv_version $VERSION" \
           --define "rcv_glibc $GLIBC" --define "rcv_stage $(pwd)/$STAGE"
  RPM=dist/installer/RCV-$VERSION-Linux-x86_64.rpm
  cp "$TOP/RPMS/x86_64/rcv-$VERSION-1.x86_64.rpm" "$RPM"
  echo "built $RPM"
else
  echo "rpmbuild not found: no .rpm built (on Debian/Ubuntu: sudo apt install rpm)"
fi
