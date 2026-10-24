"""Desktop window for RCV (RNA Conformation Viewer).

Shells, in order of preference:
  webview  native window through pywebview (WebView2 on Windows, WKWebView on macOS, GTK/Qt on Linux),
           with native Open/Save dialogs. Used by the installer builds.
  window   chromeless app window of Microsoft Edge or Google Chrome (--app), e.g. in development or on WSL.
  browser  a tab in the default browser.
  none     no window; print the address (headless use, tests).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import webbrowser
from pathlib import Path

from . import APP_NAME


def _is_wsl():
    try:
        return 'microsoft' in Path('/proc/version').read_text().lower()
    except OSError:
        return False


def find_browser():
    """Path of a Chromium-based browser (Edge or Chrome) usable for app windows, or None."""
    candidates = []
    if _is_wsl():
        candidates += ['/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
                       '/mnt/c/Program Files/Microsoft/Edge/Application/msedge.exe',
                       '/mnt/c/Program Files/Google/Chrome/Application/chrome.exe',
                       '/mnt/c/Program Files (x86)/Google/Chrome/Application/chrome.exe']
    elif sys.platform == 'win32':
        for base in (os.environ.get('PROGRAMFILES(X86)', ''), os.environ.get('PROGRAMFILES', ''),
                     os.environ.get('LOCALAPPDATA', '')):
            candidates += [os.path.join(base, 'Microsoft', 'Edge', 'Application', 'msedge.exe'),
                           os.path.join(base, 'Google', 'Chrome', 'Application', 'chrome.exe')]
    elif sys.platform == 'darwin':
        candidates += ['/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
                       '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
                       '/Applications/Chromium.app/Contents/MacOS/Chromium']
    else:
        for exe in ('google-chrome', 'google-chrome-stable', 'chromium', 'chromium-browser', 'microsoft-edge'):
            w = shutil.which(exe)
            if w:
                candidates.append(w)
    for c in candidates:
        if c and os.path.exists(c):
            return c
    return None


def window_profile_dir(exe=None):
    """A browser profile used only by the RCV window: its own icon cache (no stale logo from earlier versions),
    no mixing with the user's browsing history, tabs or extensions. Returned as a path the browser understands."""
    from . import paths as P
    if _is_wsl() and exe and str(exe).startswith('/mnt/'):
        try:                                             # the Windows browser needs a Windows path
            out = subprocess.run(['cmd.exe', '/c', 'echo %LOCALAPPDATA%'], capture_output=True, text=True, timeout=10,
                                 cwd='/mnt/c' if os.path.isdir('/mnt/c') else None).stdout.strip()
            if out and '%' not in out:
                return out.rstrip('\\') + '\\RCV\\window-profile'
        except (OSError, subprocess.SubprocessError):
            pass
        return None
    if sys.platform == 'win32' and os.environ.get('LOCALAPPDATA'):
        return os.path.join(os.environ['LOCALAPPDATA'], 'RCV', 'window-profile')
    return str(P.user_data_dir() / 'window-profile')


def app_window_command(url, profile=True):
    exe = find_browser()
    if not exe:
        return None
    cmd = [exe, f'--app={url}', '--window-size=1600,1000', '--no-first-run', '--no-default-browser-check',
           '--disable-extensions', '--disable-sync', '--class=RCV']
    prof = window_profile_dir(exe) if profile else None
    if prof:
        cmd.insert(1, f'--user-data-dir={prof}')
    return cmd


def webview_available():
    """pywebview is installed and a display is available (always true on Windows and macOS)."""
    try:
        import webview  # noqa: F401
    except ImportError:
        return False
    if sys.platform.startswith('linux') and not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
        return False
    return True


def _dialog(kind):
    """File-dialog constant for pywebview 4-6 (FileDialog enum in newer releases)."""
    import webview
    fd = getattr(webview, 'FileDialog', None)
    if fd is not None:
        return getattr(fd, kind)
    return getattr(webview, f'{kind}_DIALOG')


class Bridge:
    """Functions callable from the page as window.pywebview.api.<name>(...) (desktop shell only)."""

    def __init__(self, state):
        self._state = state
        self._window = None

    def choose_save_path(self, suggested_name, description='File'):
        ext = Path(suggested_name).suffix
        types = (f'{description} (*{ext})', 'All files (*.*)') if ext else ('All files (*.*)',)
        res = self._window.create_file_dialog(_dialog('SAVE'), save_filename=suggested_name, file_types=types)
        path = res[0] if isinstance(res, (list, tuple)) and res else res
        if not path:
            return None
        path = str(path)
        self._state.approved_paths.add(path)
        return path

    def choose_open_files(self, kind='structures'):
        types = {'fasta': ('FASTA (*.fa;*.fasta;*.fna;*.txt)', 'All files (*.*)'),
                 'rcv': ('RCV ensembles (*.rcv)', 'All files (*.*)')}.get(kind, (
            'RCV files and structures (*.rcv;*.pdb;*.cif;*.mmcif;*.ent)', 'RCV ensembles (*.rcv)',
            'Structures (*.pdb;*.cif;*.mmcif;*.ent)', 'Workspaces (*.json)', 'All files (*.*)'))
        res = self._window.create_file_dialog(_dialog('OPEN'), allow_multiple=kind not in ('fasta', 'rcv'), file_types=types)
        return [str(p) for p in (res or [])]

    def read_text_file(self, path):
        """Text of a small file the user just chose in an Open dialog (FASTA input)."""
        p = Path(path)
        if p.stat().st_size > 20 * 1024 * 1024:
            raise ValueError('file too large')
        return p.read_text(encoding='utf-8', errors='replace')

    def choose_folder(self):
        res = self._window.create_file_dialog(_dialog('FOLDER'))
        return str(res[0]) if res else None

    def toggle_fullscreen(self):
        """Window ▸ Full screen (F11) for the native window."""
        self._window.toggle_fullscreen()
        return True

    def open_url(self, url):
        """Help ▸ RCV website and links in the user guide: opened in the system browser (web links only)."""
        url = str(url)
        if not re.match(r'^https?://', url):
            raise ValueError('only web links can be opened')
        webbrowser.open(url)
        return True

    def reveal(self, path):
        p = Path(path)
        target = p if p.is_dir() else p.parent
        if sys.platform == 'win32':
            os.startfile(str(target))                      # noqa: S606 (opens Explorer)
        elif sys.platform == 'darwin':
            subprocess.Popen(['open', str(target)])
        else:
            subprocess.Popen(['xdg-open', str(target)])
        return True


def run_webview(url, state, on_closed):
    import webview
    bridge = Bridge(state)
    window = webview.create_window(APP_NAME, url, js_api=bridge, width=1560, height=980, min_size=(1000, 660),
                                   background_color='#f3f2ee', text_select=True)
    bridge._window = window
    state.shell = 'webview'
    window.events.closed += on_closed
    storage = Path(os.environ['RCV_WEBVIEW_DATA']) if os.environ.get('RCV_WEBVIEW_DATA') else None
    webview.start(private_mode=False, storage_path=str(storage) if storage else None)   # blocks until closed


def open_external(url, mode):
    """Open the interface outside pywebview; returns the shell actually used."""
    if mode in ('auto', 'window'):
        cmd = app_window_command(url)
        if cmd:
            try:
                subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
                return 'window'
            except OSError:
                pass
    if _is_wsl() and shutil.which('cmd.exe'):
        subprocess.Popen(['cmd.exe', '/c', 'start', '', url.replace('&', '^&')], stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL, cwd='/mnt/c' if os.path.isdir('/mnt/c') else None)
        return 'browser'
    webbrowser.open(url)
    return 'browser'
