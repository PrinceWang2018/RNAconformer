"""File locations for RCV (RNA Conformation Viewer) on Windows, macOS and Linux, in development and frozen (installer) builds."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

APP = 'RCV'
LEGACY_APPS = ('RNAconfSwitch',)        # earlier name of the application (0.4)
PKG_DIR = Path(__file__).resolve().parent


def frozen():
    """True inside a PyInstaller build (Windows .exe / macOS .app)."""
    return bool(getattr(sys, 'frozen', False))


def resource_dir():
    if frozen():
        return Path(getattr(sys, '_MEIPASS', Path(sys.executable).parent)) / 'rcv'
    return PKG_DIR


def static_dir():
    return resource_dir() / 'static'


def demo_dir():
    """Demo ensembles: bundled with installers, examples/viewer_demo in a source checkout."""
    for p in (resource_dir() / 'demo', PKG_DIR.parents[1] / 'examples' / 'viewer_demo'):
        if p.is_dir():
            return p
    return None


def source_checkout():
    """Repository root when running from a source checkout, else None."""
    root = PKG_DIR.parents[1]
    return root if (root / 'pyproject.toml').exists() and not frozen() else None


def _data_dir(app):
    if sys.platform == 'win32':
        return Path(os.environ.get('APPDATA', Path.home() / 'AppData' / 'Roaming')) / app
    if sys.platform == 'darwin':
        return Path.home() / 'Library' / 'Application Support' / app
    return Path(os.environ.get('XDG_DATA_HOME', Path.home() / '.local' / 'share')) / app.lower()


def user_data_dir():
    env = os.environ.get('RCV_HOME') or os.environ.get('RNACONFSWITCH_HOME')
    return Path(env) if env else _data_dir(APP)


def legacy_data_dirs():
    """User data folders of earlier versions (settings and plugins are still read from there)."""
    if os.environ.get('RCV_HOME') or os.environ.get('RNACONFSWITCH_HOME'):
        return []
    return [d for d in (_data_dir(a) for a in LEGACY_APPS) if d.is_dir()]


def documents_dir(app=APP):
    docs = Path.home() / 'Documents'
    if sys.platform in ('win32', 'darwin') or docs.is_dir():
        return docs / app
    return Path.home() / app


def legacy_projects_dirs():
    """Prediction folders of earlier versions, still listed under recent predictions."""
    return [d / 'projects' for d in (documents_dir(a) for a in LEGACY_APPS) if (d / 'projects').is_dir()]


def plugin_dirs():
    dirs = [user_data_dir() / 'plugins'] + [d / 'plugins' for d in legacy_data_dirs() if (d / 'plugins').is_dir()]
    legacy = Path(os.environ.get('RNAFOLD3D_HOME', Path.home() / '.rnafold3d')) / 'plugins'
    if legacy.is_dir():
        dirs.append(legacy)
    return dirs


def settings_path():
    return user_data_dir() / 'settings.json'


def load_settings():
    for p in [settings_path()] + [d / 'settings.json' for d in legacy_data_dirs()]:
        try:
            return json.loads(p.read_text(encoding='utf-8'))
        except (OSError, ValueError):
            continue
    return {}


def save_settings(values):
    p = settings_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    current = load_settings()
    current.update(values)
    p.write_text(json.dumps(current, indent=2), encoding='utf-8')
    return current


def default_projects_dir():
    return Path(load_settings().get('projects_dir') or (documents_dir() / 'projects'))


def engine_command(args):
    """Command line that runs the rnafold3d engine CLI in a separate process (works in frozen builds)."""
    if frozen():
        return [sys.executable, 'engine', *args]
    return [sys.executable, '-m', 'rcv', 'engine', *args]


def engine_env():
    env = dict(os.environ)
    if not frozen():
        src = str(PKG_DIR.parent)
        env['PYTHONPATH'] = os.pathsep.join([src] + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
    env['PYTHONUNBUFFERED'] = '1'
    env['PYTHONIOENCODING'] = 'utf-8'
    for v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        env.setdefault(v, '1')
    return env
