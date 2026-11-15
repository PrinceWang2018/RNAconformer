"""Project library: every prediction run, .rcv file and structure the user has opened, plus the runs and .rcv files
found in the watched folders (the projects folder and any folder the user adds), with a thumbnail and the
parameters, searchable in File > Project library. Coming back after a month is one click.

The index is <user data>/library/index.json; thumbnails are small JPEGs next to it (thumbs/<key>.jpg), drawn by
the viewer when an item is shown. Nothing is copied or moved: the library points to the files where they are.
"""
from __future__ import annotations

import base64
import hashlib
import json
import os
import threading
import time
from pathlib import Path

from . import paths as P

MAX_ITEMS = 2000
MAX_THUMB = 400 * 1024
STRUCT_EXT = ('.pdb', '.ent', '.cif', '.mmcif')
_lock = threading.Lock()


def lib_dir():
    return P.user_data_dir() / 'library'


def key_of(path):
    return hashlib.sha1(str(Path(path).expanduser().resolve()).encode('utf-8', 'replace')).hexdigest()[:16]


def _load():
    try:
        d = json.loads((lib_dir() / 'index.json').read_text(encoding='utf-8'))
        return d if isinstance(d, dict) and isinstance(d.get('items'), dict) else {'items': {}, 'folders': []}
    except (OSError, ValueError):
        return {'items': {}, 'folders': []}


def _save(d):
    lib_dir().mkdir(parents=True, exist_ok=True)
    tmp = lib_dir() / 'index.json.part'
    tmp.write_text(json.dumps(d, indent=1, ensure_ascii=False), encoding='utf-8')
    tmp.replace(lib_dir() / 'index.json')


def run_parameters(manifest):
    """Short, readable parameters of a prediction run (from its manifest)."""
    if not isinstance(manifest, dict):
        return {}
    p = manifest.get('parameters') or {}
    inp = manifest.get('input') or {}
    tpl = inp.get('templates') if isinstance(inp.get('templates'), dict) else {}
    out = {'engine': (manifest.get('software') or {}).get('engine'), 'models': p.get('n_models'),
           'trajectories': p.get('trajectories_run') or p.get('trajectories'), 'seed': p.get('seed'),
           'temperature': p.get('temperature'), 'templates': len(tpl.get('hits') or []) if tpl else None,
           '2d': 'given' if inp.get('secondary_structure_mode') == 'user' else 'predicted' if inp.get('secondary_structure_mode') else None,
           'wall_time_s': (manifest.get('statistics') or {}).get('wall_time_s')}
    return {k: v for k, v in out.items() if v not in (None, '')}


def describe_doc(doc, manifest=None):
    """Library fields of an opened document: where it is, what it is, and the parameters."""
    src = doc.get('source') or {}
    pkg = src.get('path') if src.get('kind') == 'package' or (src.get('path') and Path(src['path']).is_dir()) else None
    if src.get('rcv') and pkg and Path(src['rcv']).parent == Path(pkg):
        path, kind = pkg, 'run'                          # the .rcv RCV keeps inside a prediction folder: the run
    elif src.get('rcv'):
        path, kind = src['rcv'], 'rcv'
    elif src.get('kind') == 'package' and src.get('path'):
        path, kind = src['path'], 'run'
    elif src.get('kind') == 'pdb' and src.get('path'):
        path, kind = src['path'], 'pdb'
    elif src.get('paths'):
        path, kind = src['paths'][0] if len(src['paths']) == 1 else str(Path(src['paths'][0]).parent), 'structure'
    elif src.get('path'):
        path, kind = src['path'], 'run' if Path(src['path']).is_dir() else 'structure'
    else:
        return None                                       # uploads and fetched-in-memory documents have no file
    params = run_parameters(manifest) if manifest else {}
    return {'path': str(Path(path).expanduser().resolve()), 'kind': kind, 'name': doc.get('name') or Path(path).stem,
            'length': len(doc.get('sequence') or '') or None, 'n_conformations': len(doc.get('conformations') or []) or None,
            'sequence': (doc.get('sequence') or '')[:80], 'origin': doc.get('origin'), 'params': params}


def record(doc, manifest=None):
    """Remember that a document was opened (called by the server for every document with a file)."""
    info = describe_doc(doc, manifest)
    if not info:
        return None
    k = key_of(info['path'])
    with _lock:
        d = _load()
        old = d['items'].get(k, {})
        d['items'][k] = {**old, **{x: v for x, v in info.items() if v not in (None, {}, '')}, 'key': k,
                         'opened': time.time(), 'opens': old.get('opens', 0) + 1,
                         'first_opened': old.get('first_opened') or time.time()}
        if len(d['items']) > MAX_ITEMS:                    # forget the least recently opened
            for kk, _ in sorted(d['items'].items(), key=lambda kv: kv[1].get('opened', 0))[:len(d['items']) - MAX_ITEMS]:
                d['items'].pop(kk, None)
        _save(d)
    return k


def save_thumbnail(path, data_url):
    k = key_of(path)
    b64 = str(data_url).split(',', 1)[-1]
    raw = base64.b64decode(b64.encode('ascii'), validate=False)
    if not raw or len(raw) > MAX_THUMB or raw[:2] != b'\xff\xd8':
        raise ValueError('a thumbnail is a JPEG of at most 400 kB')
    (lib_dir() / 'thumbs').mkdir(parents=True, exist_ok=True)
    (lib_dir() / 'thumbs' / f'{k}.jpg').write_bytes(raw)
    with _lock:
        d = _load()
        if k in d['items']:
            d['items'][k]['thumb'] = time.time()
            _save(d)
    return {'key': k}


def thumbnail(key):
    if not all(c in '0123456789abcdef' for c in str(key)):
        raise KeyError(key)
    f = lib_dir() / 'thumbs' / f'{key}.jpg'
    if not f.is_file():
        raise KeyError(key)
    return f.read_bytes()


def forget(path):
    k = key_of(path)
    with _lock:
        d = _load()
        d['items'].pop(k, None)
        _save(d)
    try:
        (lib_dir() / 'thumbs' / f'{k}.jpg').unlink()
    except OSError:
        pass
    return {'forgotten': k}


def folders():
    return [f for f in _load().get('folders', []) if isinstance(f, str)]


def set_folder(path, add=True):
    p = str(Path(path).expanduser().resolve())
    if add and not Path(p).is_dir():
        raise ValueError(f'not a folder: {p}')
    with _lock:
        d = _load()
        fs = [f for f in d.get('folders', []) if f != p]
        if add:
            fs.append(p)
        d['folders'] = fs
        _save(d)
    return fs


def _rcv_files(root, limit=300, max_depth=4):
    out = []
    root = Path(root)
    for dirpath, dirnames, filenames in os.walk(root):
        if len(Path(dirpath).relative_to(root).parts) >= max_depth:
            dirnames[:] = []
        dirnames[:] = [x for x in dirnames if not x.startswith('.') and x not in ('node_modules', '__pycache__', 'candidates', 'conformations', 'raw')]
        if (Path(dirpath) / 'metadata' / 'manifest.json').is_file():
            dirnames[:] = []                             # a prediction folder: listed as a run, its .rcv is a cache
            continue
        for f in filenames:
            if f.lower().endswith('.rcv'):
                out.append(Path(dirpath) / f)
                if len(out) >= limit:
                    return out
    return out


def items(watched_roots, find_packages):
    """Everything in the library: opened items plus what the watched folders hold, newest first."""
    d = _load()
    out = {k: dict(v) for k, v in d['items'].items()}
    roots = [str(r) for r in watched_roots] + folders()
    seen_roots = []
    for r in roots:
        if r in seen_roots or not Path(r).is_dir():
            continue
        seen_roots.append(r)
        for p in find_packages([r]):
            k = key_of(p['path'])
            base = out.get(k, {})
            out[k] = {'path': p['path'], 'kind': 'run', 'name': p.get('name'), 'length': p.get('length'),
                      'n_conformations': p.get('n_conformations'), 'sequence': p.get('sequence'),
                      'params': base.get('params') or {'engine': p.get('engine')}, 'created': p.get('created'),
                      **base, 'key': k, 'watched': True}
        for f in _rcv_files(r):
            k = key_of(f)
            base = out.get(k, {})
            out[k] = {'path': str(f.resolve()), 'kind': 'rcv', 'name': f.stem, **base, 'key': k, 'watched': True,
                      'modified': f.stat().st_mtime}
    thumbs = lib_dir() / 'thumbs'
    for k, v in out.items():
        v['exists'] = Path(v.get('path', '')).exists()
        v['thumb'] = (thumbs / f'{k}.jpg').is_file()
        v['when'] = max(v.get('opened') or 0, v.get('modified') or 0, _iso_ts(v.get('created')))
    return {'items': sorted(out.values(), key=lambda v: v['when'], reverse=True), 'folders': folders(),
            'default_folders': [str(r) for r in watched_roots]}


def _iso_ts(s):
    if not s:
        return 0
    try:
        import datetime
        return datetime.datetime.fromisoformat(str(s)).timestamp()
    except ValueError:
        return 0
