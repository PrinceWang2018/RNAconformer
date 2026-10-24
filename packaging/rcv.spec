# -*- mode: python ; coding: utf-8 -*-
# PyInstaller recipe for RCV (one-folder build; a .app bundle on macOS).
#   pyinstaller packaging/rcv.spec --noconfirm
# Build on the target OS (Windows on Windows, macOS on macOS, Linux on Linux); PyInstaller does not cross-compile.
# Windows: packaging/build_windows.ps1, macOS: packaging/macos/build_macos.sh, Linux: packaging/linux/build_linux.sh.
import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules

ROOT = Path(SPECPATH).resolve().parent
SRC = ROOT / 'src'
ICONS = ROOT / 'packaging' / 'icons'
sys.path.insert(0, str(SRC))
from rcv import __version__  # noqa: E402

datas = [
    (str(SRC / 'rcv' / 'static'), 'rcv/static'),
    (str(SRC / 'rnafold3d' / 'data'), 'rnafold3d/data'),
    (str(ROOT / 'examples' / 'viewer_demo'), 'rcv/demo'),
    (str(SRC / 'rcv' / 'templates'), 'rcv/templates'),
    (str(ROOT / 'docs' / 'PLUGINS.md'), 'rcv/docs'),          # Help menu: plugin guide, user guide, what is new
    (str(ROOT / 'docs' / 'RCV.md'), 'rcv/docs'),
    (str(ROOT / 'docs' / 'images'), 'rcv/docs/images'),
    (str(ROOT / 'CHANGELOG.md'), 'rcv/docs'),
    (str(ROOT / 'LICENSE'), '.'),
    (str(ROOT / 'packaging' / 'THIRD-PARTY-NOTICES.txt'), '.'),
]
datas += collect_data_files('RNA')                    # ViennaRNA energy parameter files, if packaged as data
binaries = collect_dynamic_libs('RNA')
hiddenimports = (collect_submodules('rnafold3d') + collect_submodules('rcv') +
                 collect_submodules('fontTools.subset') + collect_submodules('fontTools.ttLib') +
                 ['RNA', 'threadpoolctl', 'scipy.optimize', 'scipy.linalg', 'scipy.spatial', 'webview', 'fontTools.ttLib'])

a = Analysis(
    [str(ROOT / 'packaging' / 'rcv_entry.py')],
    pathex=[str(SRC)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    excludes=['tkinter', 'matplotlib', 'PyQt5', 'PyQt6', 'PySide2', 'IPython', 'pytest'],
    noarchive=False,
)


def _newest_vc_runtime(binaries):
    """Windows: ship the newest copy of each Visual C++ runtime DLL.

    PyInstaller takes MSVCP140.dll & co. from the first folder on PATH that has them, which can be an old copy
    installed by another program (for example a Java JDK). Libraries built with a current compiler (ViennaRNA)
    then fail with "DLL initialization routine failed", so the copy in System32 is used when it is newer.
    """
    import os
    import pefile

    def version(path):
        try:
            pe = pefile.PE(str(path), fast_load=True)
            pe.parse_data_directories([pefile.DIRECTORY_ENTRY['IMAGE_DIRECTORY_ENTRY_RESOURCE']])
            ffi = pe.VS_FIXEDFILEINFO[0]
            return (ffi.FileVersionMS >> 16, ffi.FileVersionMS & 0xFFFF, ffi.FileVersionLS >> 16, ffi.FileVersionLS & 0xFFFF)
        except Exception:
            return (0,)
    names = {'msvcp140.dll', 'msvcp140_1.dll', 'msvcp140_2.dll', 'msvcp140_atomic_wait.dll', 'msvcp140_codecvt_ids.dll',
             'vcruntime140.dll', 'vcruntime140_1.dll', 'concrt140.dll', 'vcomp140.dll'}
    system32 = Path(os.environ.get('SystemRoot', r'C:\Windows')) / 'System32'
    out = []
    for dest, src, kind in binaries:
        cand = system32 / Path(dest).name
        if Path(dest).name.lower() in names and cand.is_file() and version(cand) > version(src):
            print(f'RCV: using {cand} ({".".join(map(str, version(cand)))}) instead of {src}')
            src = str(cand)
        out.append((dest, src, kind))
    return out


if sys.platform == 'win32':
    a.binaries = _newest_vc_runtime(a.binaries)

pyz = PYZ(a.pure)
# Interpreter options of the frozen app (it ignores PYTHON* environment variables): UTF-8 mode, so files and the
# engine's progress pipe are UTF-8 whatever the Windows code page (e.g. GBK), and unbuffered output, so progress
# lines reach the app as they are written.
options = [('X utf8', None, 'OPTION'), ('u', None, 'OPTION')]
exe = EXE(
    pyz, a.scripts, options,
    exclude_binaries=True,
    name='RCV',
    console=False,                                     # windowed app; the engine runs as a child process
    icon=str(ICONS / ('rcv.icns' if sys.platform == 'darwin' else 'rcv.ico')),
    upx=False,
)
coll = COLLECT(exe, a.binaries, a.datas, name='RCV', upx=False)

if sys.platform == 'darwin':
    app = BUNDLE(
        coll,
        name='RCV.app',
        icon=str(ICONS / 'rcv.icns'),
        bundle_identifier='org.rcv.app',
        version=__version__,
        info_plist={
            'CFBundleDisplayName': 'RCV',
            'CFBundleName': 'RCV',
            'CFBundleShortVersionString': __version__,
            'NSHighResolutionCapable': True,
            'LSMinimumSystemVersion': '14.0',          # the NumPy/SciPy wheels pip picks on macOS 15 need macOS 14
            'NSHumanReadableCopyright': 'Copyright © 2026 Zixiang Wang. GPL-3.0.',
            'CFBundleDocumentTypes': [
                {'CFBundleTypeName': 'RCV conformational ensemble', 'CFBundleTypeRole': 'Editor',
                 'LSItemContentTypes': ['org.rcv.ensemble'], 'LSHandlerRank': 'Owner'},
                {'CFBundleTypeName': 'RNA structure', 'CFBundleTypeRole': 'Viewer',
                 'LSItemContentTypes': ['public.data'], 'CFBundleTypeExtensions': ['pdb', 'cif', 'mmcif', 'ent']}],
            'UTExportedTypeDeclarations': [
                {'UTTypeIdentifier': 'org.rcv.ensemble', 'UTTypeDescription': 'RCV conformational ensemble',
                 'UTTypeConformsTo': ['public.data', 'public.zip-archive'],
                 'UTTypeTagSpecification': {'public.filename-extension': ['rcv'], 'public.mime-type': 'application/vnd.rcv+zip'}}],
        },
    )
