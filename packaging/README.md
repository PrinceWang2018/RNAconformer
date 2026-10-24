# Building RCV (RNA Conformation Viewer) desktop installers

RCV ships as a self-contained desktop application: Python, the prediction engine (numpy, scipy,
ViennaRNA), the user interface and the demo data are bundled with PyInstaller, and the window is a native
web view (WebView2 on Windows, WKWebView on macOS) through pywebview; on Linux it is a Chrome, Chromium or Edge
app window (or the default browser). All computation runs on the user's computer; nothing is sent anywhere.

PyInstaller does not cross-compile: build the Windows installer on Windows, the macOS app on a Mac (Apple
Silicon and Intel need separate builds, or a CI runner of each kind) and the Linux packages on Linux.

| File | Purpose |
|---|---|
| `rcv_entry.py` | executable entry point (`freeze_support` first: the same executable runs the engine and its worker processes) |
| `rcv.spec` | PyInstaller recipe (one-folder build; `.app` bundle on macOS) |
| `requirements-desktop.txt` | build environment (Python 3.12, ViennaRNA 2.7.2 wheel, pywebview, PyInstaller, certifi for *Fetch from the PDB*, fontTools for embedding the figure font in PDFs) |
| `build_windows.ps1` | Windows: venv → PyInstaller → frozen-engine self-test → Inno Setup installer |
| `windows/RCV.iss` | Inno Setup script (Start-menu entry, optional desktop shortcut, optional `.rcv` file association, uninstaller) |
| `macos/build_macos.sh` | macOS: venv → PyInstaller `.app` → self-test → `.dmg` with an *Applications* link and the user guide (optional signing and notarisation) |
| `linux/build_linux.sh` | Linux: venv → PyInstaller → self-test (engine, tools, server) → `.deb` (dpkg-deb) and `.rpm` (rpmbuild) |
| `linux/rcv-rpm.spec`, `linux/rcv.desktop`, `linux/rcv-mime.xml`, `linux/icons/` | RPM recipe, menu entry (with *RCV with demo ensembles* and *RCV User Guide*), `.rcv` file type, icons 16–512 px |
| `make_docs.py` | user guide, plugin guide and changelog as web pages (installers: Start menu / disk image; website: `--site`) |
| `THIRD-PARTY-NOTICES.txt` | licences and credits of the bundled software (ViennaRNA, Python, NumPy, SciPy, …), installed next to RCV |
| `../.github/workflows/installers.yml` | builds all installers on GitHub (Windows x64, macOS Apple Silicon, macOS Intel, Linux `.deb` and `.rpm`) and test-installs the Linux packages |
| `icons/` | `rcv.ico`, `rcv.icns` and a 1024 px `rcv.png` rendered from `src/rcv/static/img/icon.svg` |

## Windows (10/11, x64)

```powershell
winget install Python.Python.3.12
winget install JRSoftware.InnoSetup          # optional, for the Setup.exe
powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1
```

Output: `dist\RCV\RCV.exe` (portable folder) and
`dist\installer\RCV-<version>-Windows-x64-Setup.exe`. The Microsoft Edge WebView2 runtime is part of
Windows 11 and current Windows 10; on older systems install the Evergreen runtime from Microsoft. Without
pywebview/WebView2 the app falls back to an Edge app window.

Unsigned installers trigger a SmartScreen warning; sign `RCV.exe` and the installer with a
code-signing certificate (`signtool sign /fd sha256 /tr <timestamp-url> /td sha256 ...`) before distribution.

## macOS (14+, Apple Silicon or Intel)

```bash
brew install python@3.12        # or the python.org installer
bash packaging/macos/build_macos.sh
```

Output: `dist/RCV.app` and `dist/installer/RCV-<version>-macOS-<arch>.dmg`. For distribution
outside your own Macs, set `CODESIGN_ID` (Developer ID Application certificate) and `NOTARY_PROFILE`
(`xcrun notarytool store-credentials`) before running the script; unsigned apps must be opened with
right-click ▸ Open the first time.

The app needs **macOS 14 (Sonoma) or later**: on a macOS 15 build machine pip installs the NumPy and SciPy wheels
built for macOS 14 (`macosx_14_0`), and `rcv.spec` sets `LSMinimumSystemVersion` to match, so older systems show
a clear message instead of failing at start. (Supporting macOS 11–13 would need a build machine with macOS 13 or
the older OpenBLAS wheels.)

## Linux (x86_64: .deb and .rpm)

```bash
sudo apt install python3.12-venv rpm      # Debian/Ubuntu; rpm only for the .rpm
bash packaging/linux/build_linux.sh
```

Output: `dist/RCV/RCV` (portable folder), `dist/installer/RCV-<version>-Linux-amd64.deb` and
`RCV-<version>-Linux-x86_64.rpm`. Both install RCV to `/opt/rcv`, the command `/usr/bin/rcv`, a menu entry
(with the actions *RCV with demo ensembles* and *RCV User Guide*), icons, the `.rcv` file type and the user guide
in `/usr/share/doc/rcv/html`. Install with `sudo apt install ./RCV-<version>-Linux-amd64.deb` or
`sudo dnf install ./RCV-<version>-Linux-x86_64.rpm`.

The bundle runs on systems whose glibc is at least as new as the build machine's, so build on an old
distribution: GitHub builds in a Debian 11 container (glibc 2.31), which covers Ubuntu 20.04+, Debian 11+,
Fedora 32+ and RHEL/Rocky/Alma 9+. pywebview is not bundled on Linux (its GTK/Qt web views are system
libraries): RCV opens in a Chrome, Chromium or Edge app window when one is installed, else in the default
browser, and quits when that window is closed.

## macOS without a Mac: GitHub Actions

`.github/workflows/installers.yml` builds the Windows installer, the macOS disk images for Apple Silicon
(`macos-15` runner) and Intel (`macos-15-intel`; GitHub keeps Intel Macs until autumn 2027) and the Linux
packages (Debian 11 container; then installed, run and removed on Ubuntu 20.04 and 24.04, Debian 12, Fedora,
AlmaLinux 9 and openSUSE Leap 15.6) on GitHub's machines. *Actions ▸ Build installers ▸ Run workflow* puts the installers on the run page as downloads; pushing
a version tag (`git tag v0.6.0 && git push origin v0.6.0`) also attaches them to the release of the website
repository when the secret `WEBSITE_RELEASE_TOKEN` (fine-grained token, *Contents: read and write* on
`PrinceWang2018/RNAconfViewer_website`) is set.

## Where the installers are published

The download links of www.rnaconfviewer.com have the form
`https://www.rnaconfviewer.com/download/v<version>/<file>`. Cloudflare Pages serves at most 25 MiB per file and
the installers are larger (Windows ≈ 41 MiB), so the files live in a **GitHub release of the website repository**
(tag `v<version>`) and the website's `_redirects` file forwards the `/download/...` addresses there. GitHub counts
the downloads, which feeds the website's *Statistics* tab.

## What was verified

* ViennaRNA 2.7.2 publishes wheels for Windows x64 (cp312), macOS arm64 and macOS x86_64 (10.14+), so the
  engine's only compiled third-party dependency is installable on both targets.
* **Windows 10 x64 (0.6.0, PyInstaller 6.22.3, Python 3.12.15, Inno Setup 6.7.3)**: the frozen app predicts
  2D structures with ViennaRNA and runs a full 3D prediction with worker processes (`RCV engine predict`),
  including a run name with non-ASCII characters on a Chinese (GBK) system; `RCV convert` / `RCV info` work;
  the bundled server passes the GUI scenarios for predictions and FASTA batches started in the app, exports,
  video, the user guide with its pictures, the window layout and all ten languages; the installer compiles
  (≈ 41 MiB).
* Not yet verified: a macOS build (use the GitHub workflow), the pywebview window and its native Open/Save
  dialogs, the `.rcv` file association after installation, code signing.
* macOS: double-clicking an `.rcv` file in Finder starts RCV, but pywebview does not pass the file to the app
  (Apple-event handling); open it with File ▸ Open until that is added. Windows passes it on the command line.

## Notes for maintainers

* **Visual C++ runtime (Windows)**: PyInstaller copies `MSVCP140.dll` & co. from the first folder on PATH that
  has them — possibly an old copy from another program (a Java JDK on the build machine did this), which makes
  ViennaRNA fail with "DLL initialization routine failed". `rcv.spec` replaces such copies with the newer ones
  from `System32` and prints `RCV: using ...` when it does.
* **UTF-8 mode**: frozen apps ignore `PYTHON*` environment variables, so `rcv.spec` sets the interpreter
  options `-X utf8` and `-u`: files, the engine's progress pipe and the logs are UTF-8 and unbuffered whatever
  the Windows code page.

* The bundle includes `rcv/templates/` (plugin templates), `docs/PLUGINS.md` (shown by Plugins ▸ Plugin developer
  guide), the figure font `static/fonts/LiberationSans-*.ttf` with its licence (SIL OFL 1.1 — keep
  `LICENSE-OFL.txt` next to the fonts) and the fontTools modules used for PDF subsetting. Without fontTools PDFs
  fall back to the standard PDF fonts.
* The Edge/Chrome app window (used when pywebview is missing) runs with its own profile
  (`%LOCALAPPDATA%\RCV\window-profile`, `--class=RCV`), so the taskbar shows the RCV icon, not the browser's,
  and settings stored by the window are separate from the user's browser. Deleting that folder resets them.

* Frozen builds locate resources through `rcv.paths.resource_dir()` (`sys._MEIPASS/rcv`);
  the engine's parameters are loaded with `importlib.resources` from `rnafold3d/data`.
* Predictions started from the app run `RCV engine predict ...` as a child process
  (`rcv.paths.engine_command`), so the window stays responsive and cancelling kills the run.
* User data: settings, logs and plugins in `%APPDATA%\RCV` / `~/Library/Application Support/RCV` /
  `~/.local/share/rcv`;
  predictions in `Documents/RCV/projects` (changeable in the app).
