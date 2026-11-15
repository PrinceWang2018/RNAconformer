"""Check for a newer version of RCV, and pick the installer for this computer.

Nothing is downloaded or installed here: the app shows what is available and, when the user agrees, opens the
installer's address in the web browser. The check asks GitHub for the newest release of the website repository
(where the installers are published); when GitHub cannot be reached it reads the website's statistics file, which
lists the same releases and is refreshed daily. The request carries only "RCV/<version>" as its User-Agent.

Automatic checks (when RCV starts, at most once a day) are on in installed builds and off when running from source;
the user can turn them off. RCV_NO_UPDATE_CHECK=1 turns them off as well.
"""
from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

from . import WEBSITE, __version__
from . import paths as P

RELEASES_API = 'https://api.github.com/repos/PrinceWang2018/RNAconfViewer_website/releases/latest'
STATS_URL = WEBSITE + 'stats.json'
RELEASE_REPO = 'PrinceWang2018/RNAconfViewer_website'
AUTO_INTERVAL = 20 * 3600            # seconds between automatic checks

# platform key -> (label, end of the installer's file name)
PLATFORMS = {
    'windows-x64': ('Windows (64-bit)', '-Windows-x64-Setup.exe'),
    'macos-arm64': ('macOS (Apple Silicon)', '-macOS-arm64.dmg'),
    'macos-x86_64': ('macOS (Intel)', '-macOS-x86_64.dmg'),
    'linux-deb': ('Linux (.deb: Debian, Ubuntu, Linux Mint)', '-Linux-amd64.deb'),
    'linux-rpm': ('Linux (.rpm: Fedora, RHEL, openSUSE)', '-Linux-x86_64.rpm'),
}


def parse_version(v):
    """'v0.7.10' -> (0, 7, 10); text after the numbers (0.8.0rc1) is ignored."""
    m = re.match(r'v?(\d+(?:\.\d+)*)', str(v or '').strip())
    return tuple(int(x) for x in m.group(1).split('.')) if m else ()


def is_newer(latest, current):
    a, b = parse_version(latest), parse_version(current)
    n = max(len(a), len(b))
    return bool(a) and a + (0,) * (n - len(a)) > b + (0,) * (n - len(b))


def _os_release():
    for f in ('/etc/os-release', '/usr/lib/os-release'):
        try:
            text = Path(f).read_text(encoding='utf-8', errors='replace')
        except OSError:
            continue
        return {k: v.strip('"\'') for k, _, v in (ln.partition('=') for ln in text.splitlines()) if k}
    return {}


def _rosetta():
    """True for an Intel build running on an Apple Silicon Mac."""
    try:
        r = subprocess.run(['sysctl', '-n', 'sysctl.proc_translated'], capture_output=True, text=True, timeout=5)
        return r.stdout.strip() == '1'
    except (OSError, subprocess.SubprocessError):
        return False


def platform_key(system=None, machine=None, os_release=None, dpkg_installed=None):
    """Which installer fits this computer: a key of PLATFORMS, or None (no installer for this platform)."""
    system = system or sys.platform
    machine = (machine or platform.machine() or '').lower()
    if system.startswith('win'):
        return 'windows-x64' if machine in ('amd64', 'x86_64', 'arm64', 'aarch64', '') else None
    if system == 'darwin':
        return 'macos-arm64' if machine in ('arm64', 'aarch64') or (machine == 'x86_64' and _rosetta()) \
            else 'macos-x86_64'
    if system.startswith('linux'):
        if machine not in ('x86_64', 'amd64'):
            return None
        if dpkg_installed is None:
            dpkg_installed = Path('/var/lib/dpkg/info/rcv.list').is_file()
        if dpkg_installed:                              # RCV itself was installed from the .deb
            return 'linux-deb'
        rel = os_release if os_release is not None else _os_release()
        family = f"{rel.get('ID', '')} {rel.get('ID_LIKE', '')}".lower()
        if re.search(r'debian|ubuntu|mint|pop|elementary|kali|raspbian', family):
            return 'linux-deb'
        if re.search(r'fedora|rhel|centos|rocky|alma|suse|opensuse|mageia|amzn|ol\b', family):
            return 'linux-rpm'
        return 'linux-deb' if Path('/usr/bin/dpkg').exists() else 'linux-rpm' if Path('/usr/bin/rpm').exists() else None
    return None


def _get_json(url, timeout=15):
    import ssl
    import urllib.request
    contexts = [ssl.create_default_context()]
    try:
        import certifi
        contexts.append(ssl.create_default_context(cafile=certifi.where()))
    except ImportError:
        pass
    last = None
    for ctx in contexts:
        try:
            req = urllib.request.Request(url, headers={'User-Agent': f'RCV/{__version__}',
                                                       'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return json.loads(r.read(4 * 1024 * 1024).decode('utf-8'))
        except Exception as exc:                        # noqa: BLE001 - any failure: try the next source
            last = exc
    raise last


def latest_release(fetch=_get_json):
    """{'version', 'tag', 'date', 'assets': [file names], 'notes', 'source'} of the newest published release."""
    errors = []
    try:
        r = fetch(RELEASES_API)
        if r.get('tag_name'):
            return {'version': '.'.join(map(str, parse_version(r['tag_name']))), 'tag': r['tag_name'],
                    'date': (r.get('published_at') or '')[:10], 'assets': [a['name'] for a in r.get('assets', [])],
                    'notes': r.get('html_url') or WEBSITE + '#releases', 'source': 'github'}
        errors.append('GitHub: no release')
    except Exception as exc:                            # noqa: BLE001
        errors.append(f'GitHub: {getattr(exc, "reason", exc)}')
    try:
        stats = fetch(STATS_URL)
        rels = [r for r in stats.get('downloads', {}).get('releases', [])
                if r.get('repo') == RELEASE_REPO and parse_version(r.get('tag'))]
        if rels:
            r = max(rels, key=lambda r: parse_version(r['tag']))
            return {'version': '.'.join(map(str, parse_version(r['tag']))), 'tag': r['tag'],
                    'date': r.get('date', ''), 'assets': [a['name'] for a in r.get('assets', [])],
                    'notes': WEBSITE + '#releases', 'source': 'website'}
        errors.append('website: no release')
    except Exception as exc:                            # noqa: BLE001
        errors.append(f'website: {getattr(exc, "reason", exc)}')
    raise RuntimeError('could not check for updates (no internet connection?): ' + '; '.join(errors))


def check(current=__version__, fetch=_get_json, key='detect', installed=None):
    """Compare the newest release with this version and pick the installer for this computer."""
    rel = latest_release(fetch)
    key = platform_key() if key == 'detect' else key
    installed = P.frozen() if installed is None else installed
    label, suffix = PLATFORMS.get(key, ('', None))
    asset = next((a for a in rel['assets'] if suffix and a.endswith(suffix)), None)
    return {'current': current, 'latest': rel['version'], 'date': rel['date'], 'newer': is_newer(rel['version'], current),
            'platform': key, 'platform_label': label, 'installed': installed, 'file': asset,
            'download': f"{WEBSITE}download/{rel['tag']}/{asset}" if asset else None,
            'page': WEBSITE + '#download', 'notes': rel['notes'], 'source': rel['source'],
            'source_zip': 'https://github.com/PrinceWang2018/RNAconformer'}


# ------------------------------------------------------------------ settings (settings.json of the user data folder)
def settings():
    s = P.load_settings()
    auto = s.get('update_check')
    if auto is None:
        auto = P.frozen()
    if os.environ.get('RCV_NO_UPDATE_CHECK'):
        auto = False
    return {'auto': bool(auto), 'skip': s.get('update_skip'), 'last': s.get('update_last', 0)}


def save(auto=None, skip=None):
    values = {}
    if auto is not None:
        values['update_check'] = bool(auto)
    if skip is not None:
        values['update_skip'] = str(skip) or None
    if values:
        P.save_settings(values)
    return settings()


def auto_check(now=None, **kw):
    """The check made when RCV starts: None when it is turned off, was made recently, or failed."""
    s = settings()
    now = time.time() if now is None else now
    if not s['auto'] or now - float(s['last'] or 0) < AUTO_INTERVAL:
        return None
    P.save_settings({'update_last': now})
    try:
        r = check(**kw)
    except Exception:                                   # noqa: BLE001 - quiet when offline
        return None
    r['skipped'] = bool(s['skip']) and s['skip'] == r['latest']
    return r
