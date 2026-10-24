"""Local application core of RCV (Python standard library only).

The user interface is a web-technology front end shown in a desktop window (see desktop.py). This module
serves it and a small JSON API on the loopback interface (127.0.0.1) only. Every API call carries a
per-session token (header X-Viewer-Token) handed to the window in the URL fragment, and the Host header
must be a loopback name, so other web pages cannot use the API. All computation stays on this computer.
"""
from __future__ import annotations

import gzip
import importlib.util
import json
import mimetypes
import os
import platform
import re
import secrets
import sys
import threading
import time
import traceback
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse

from . import APP_NAME, WEBSITE, __version__
from . import paths as P
from . import rcvfile
from .data import STRUCT_EXT, build_document, collect_paths, collect_uploads, _package_confs
from .jobs import BatchJob, ExploreJob, PredictJob
from .methods import methods_text

MAX_BODY = 512 * 1024 * 1024
mimetypes.add_type('text/javascript', '.js')
mimetypes.add_type('text/javascript', '.mjs')
mimetypes.add_type('application/json', '.json')
mimetypes.add_type('application/manifest+json', '.webmanifest')
mimetypes.add_type('image/x-icon', '.ico')


class PythonPlugin:
    """Analysis plugin: a .py file defining PLUGIN = {'id','name','metrics':[...]} and compute(sequence, residues)."""

    def __init__(self, path):
        spec = importlib.util.spec_from_file_location(f'rcv_plugin_{path.stem}', path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        info = getattr(mod, 'PLUGIN')
        self.id = re.sub(r'[^a-z0-9_-]', '', str(info.get('id', path.stem)).lower()) or path.stem
        self.name = info.get('name', self.id)
        self.metrics = list(info.get('metrics', []))
        self.compute = getattr(mod, 'compute')
        self.path = str(path)


@dataclass
class Entry:
    """An open ensemble: its document plus the source conformations (needed to extend it later)."""
    doc: dict
    confs: list
    manifest: dict | None
    warnings: set
    source: dict
    name: str
    explorations: int = 0
    extra: list = field(default_factory=list)
    session: dict | None = None          # viewer state stored in an .rcv file
    rcv_path: str | None = None


def vienna_version():
    try:
        import RNA
        return getattr(RNA, '__version__', 'installed')
    except ImportError:
        return None


class AppState:
    def __init__(self, projects_dir, library_roots, plugin_dirs, verbose=False):
        self.token = secrets.token_urlsafe(24)
        self.lock = threading.Lock()
        self.entries = {}
        self.jobs = {}
        self.projects_dir = Path(projects_dir)
        self.library_roots = [Path(p) for p in library_roots]
        self.plugin_dirs = [Path(p) for p in plugin_dirs]
        self.verbose = verbose
        self.approved_paths = set()        # save targets chosen in native dialogs (desktop shell)
        self.shell = 'browser'
        self.last_seen = time.time()
        self.closing_at = None
        self.py_plugins, self.plugin_errors = self._load_py_plugins()

    # ---- plugins
    def _load_py_plugins(self):
        out, errors = [], []
        for d in self.plugin_dirs:
            if not d.is_dir():
                continue
            for f in sorted(d.glob('*.py')):
                try:
                    out.append(PythonPlugin(f))
                except Exception as exc:
                    errors.append(f'{f.name}: {type(exc).__name__}: {exc}')
        return out, errors

    def js_plugins(self):
        items = [{'url': f'plugins/{f.name}', 'builtin': True, 'file': f.name}
                 for f in sorted((P.static_dir() / 'plugins').glob('*.js'))]
        for k, d in enumerate(self.plugin_dirs):
            if d.is_dir():
                items += [{'url': f'user-plugins/{k}/{f.name}', 'builtin': False, 'file': str(f)} for f in sorted(d.glob('*.js'))]
        return items

    # ---- documents
    def invalidate_library(self):
        self._library_cache = None

    def _store(self, entry):
        self.invalidate_library()
        with self.lock:
            self.entries[entry.doc['id']] = entry
        return summary(entry.doc)

    def open_paths(self, paths, chain=None, name=None):
        paths = [Path(p).expanduser() for p in paths]
        rcvs = [p for p in paths if p.suffix.lower() == '.rcv']
        if rcvs:
            if len(paths) > 1:
                raise ValueError('open .rcv files one at a time (each holds a complete ensemble)')
            return self.open_rcv(rcvs[0].read_bytes(), rcvs[0], name)
        if len(paths) == 1 and not chain and not name:
            cached = self._cached_rcv(paths[0])
            if cached:
                return cached
        nm, confs, manifest, warnings, source = collect_paths(paths, chain)
        doc = build_document(name or nm, confs, manifest, warnings, source, self.py_plugins)
        return self._store(Entry(doc, confs, manifest, set(warnings), source, name or nm))

    def _cached_rcv(self, pkg):
        """A prediction folder saved by RCV also holds <name>.rcv: open that (instant) instead of re-evaluating."""
        man = pkg / 'metadata' / 'manifest.json'
        if not man.is_file():
            return None
        for f in sorted(pkg.glob('*.rcv')):
            try:
                if f.stat().st_mtime < man.stat().st_mtime:
                    continue
                doc, session, info = rcvfile.read(f)
            except (OSError, ValueError, KeyError):
                continue
            if (doc.get('source') or {}).get('path') != str(pkg.resolve()) or \
                    sorted(doc.get('plugins') or []) != sorted(p.id for p in self.py_plugins):
                continue
            return self._entry_from_rcv(doc, session, info, f, manifest=json.loads(man.read_text()))
        return None

    def open_rcv(self, data, path=None, name=None):
        doc, session, info = rcvfile.unpack(data)
        if name:
            doc['name'] = name
        return self._entry_from_rcv(doc, session, info, path)

    def _entry_from_rcv(self, doc, session, info, path=None, manifest=None):
        confs = rcvfile.conformations_from_doc(doc)
        source = dict(doc.get('source') or {})
        if path:
            source['rcv'] = str(Path(path).resolve())
        doc['source'] = source
        manifest = manifest or info.get('provenance_manifest')
        entry = Entry(doc, confs, manifest, set(doc.get('warnings') or []), source, doc.get('name') or 'ensemble',
                      session=session, rcv_path=str(path) if path else None)
        out = self._store(entry)
        if session:
            out['session'] = session
        return out

    def pack_rcv(self, doc, session=None, thumbnail=None, entry=None):
        return rcvfile.pack(doc, session=session, thumbnail=thumbnail, methods=methods_text(doc, vienna_version()),
                            provenance_manifest=entry.manifest if entry else None, app_version=__version__)

    def open_uploads(self, files, chain=None, name=None):
        nm, confs, manifest, warnings, source = collect_uploads(files, chain)
        doc = build_document(name or nm, confs, manifest, warnings, source, self.py_plugins)
        return self._store(Entry(doc, confs, manifest, set(warnings), source, name or nm))

    def extend_document(self, entry, package_dir):
        """Append the conformers of an exploration run to an open ensemble (new document id)."""
        manifest = json.loads((Path(package_dir) / 'metadata' / 'manifest.json').read_text())
        warnings = set()
        new = _package_confs(Path(package_dir), manifest, None, warnings)
        entry.explorations += 1
        n = entry.explorations
        base_rank = max((c.meta.get('rank', 0) for c in entry.confs), default=len(entry.confs))
        added = []
        for c in new:
            meta = dict(c.meta)
            cid = f'pull{n}_{meta.get("id", "conf")}'
            meta.update(id=cid, rank=base_rank + meta.get('rank', 1),
                        explored={'pull': n, 'package': str(package_dir), **(meta.get('pull') or {})})
            c.meta = meta
            c.label = f'Pull {n} · {meta["rank"] - base_rank}'
            added.append(cid)
        confs = entry.confs + new
        doc = build_document(entry.name, confs, entry.manifest, entry.warnings | warnings, entry.source, self.py_plugins)
        old_id = entry.doc['id']
        entry.doc, entry.confs = doc, confs
        with self.lock:
            self.entries.pop(old_id, None)
            self.entries[doc['id']] = entry
        return {**summary(doc), 'replaces': old_id}, added


def summary(doc):
    return {'id': doc['id'], 'name': doc['name'], 'n_conformations': len(doc['conformations']),
            'length': len(doc['sequence']), 'origin': doc['origin'], 'source': doc.get('source')}


def entry_summary(entry):
    out = summary(entry.doc)
    if entry.session:
        out['session'] = entry.session
    return out


# ----------------------------------------------------------------------------- library & browsing
def find_packages(roots, limit=400, max_depth=6):
    seen, out = set(), []
    for root in roots:
        root = Path(root).expanduser()
        if not root.is_dir():
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            depth = len(Path(dirpath).relative_to(root).parts)
            dirnames[:] = [d for d in dirnames if not d.startswith('.') and d not in
                           ('node_modules', '__pycache__', 'candidates', 'conformations', 'raw', 'explorations')]
            if depth >= max_depth:
                dirnames[:] = []
            if Path(dirpath).name == 'metadata' and 'manifest.json' in filenames:
                pkg = Path(dirpath).parent.resolve()
                if pkg in seen:
                    continue
                seen.add(pkg)
                try:
                    m = json.loads((Path(dirpath) / 'manifest.json').read_text())
                except (OSError, ValueError):
                    continue
                if not str(m.get('schema', '')).startswith('rnafold3d.ensemble'):
                    continue
                out.append({'path': str(pkg), 'name': m.get('input', {}).get('name') or pkg.name,
                            'length': m.get('input', {}).get('length'), 'n_conformations': len(m.get('conformers', [])),
                            'engine': m.get('software', {}).get('engine'), 'created': m.get('created'),
                            'sequence': (m.get('input', {}).get('sequence') or '')[:60], 'root': str(root)})
                dirnames[:] = []
                if len(out) >= limit:
                    break
    out.sort(key=lambda x: x.get('created') or '', reverse=True)
    return out


def browse(path):
    p = Path(path).expanduser() if path else P.documents_dir()
    if not p.exists():
        p = Path.home()
    p = p.resolve()
    if p.is_file():
        p = p.parent
    entries = []
    try:
        for e in sorted(p.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            if e.name.startswith('.'):
                continue
            try:
                is_dir = e.is_dir()
            except OSError:
                continue
            if is_dir:
                entries.append({'name': e.name, 'path': str(e), 'type': 'dir',
                                'package': (e / 'metadata' / 'manifest.json').exists()})
            elif e.suffix.lower() in STRUCT_EXT + ('.rcv',):
                entries.append({'name': e.name, 'path': str(e), 'type': 'file', 'size': e.stat().st_size,
                                'rcv': e.suffix.lower() == '.rcv'})
    except PermissionError:
        raise ValueError(f'permission denied: {p}')
    places = [{'name': 'Projects folder', 'path': str(P.default_projects_dir())}, {'name': 'Home folder', 'path': str(Path.home())}]
    if P.demo_dir():
        places.append({'name': 'Demo data', 'path': str(P.demo_dir())})
    if P.source_checkout():
        places.append({'name': 'Source folder', 'path': str(P.source_checkout())})
    win_users = Path('/mnt/c/Users')
    if win_users.is_dir():
        for u in sorted(win_users.iterdir()):
            if (u / 'Downloads').is_dir() and u.name not in ('Public', 'Default', 'All Users', 'Default User'):
                places.append({'name': f'Windows: {u.name}', 'path': str(u)})
    return {'path': str(p), 'parent': str(p.parent) if p.parent != p else None, 'entries': entries[:2000],
            'places': places, 'is_package': (p / 'metadata' / 'manifest.json').exists()}


def fetch_pdb(pdb_id, folder):
    """Download an entry from the RCSB PDB as mmCIF (only when the user asks for it with the fetch command)."""
    import ssl
    import urllib.error
    import urllib.request
    pid = str(pdb_id).strip().upper()
    if not re.fullmatch(r'[0-9][A-Z0-9]{3}', pid):
        raise ValueError('a PDB id has 4 characters, e.g. 1ANR')
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / f'{pid}.cif'
    if target.is_file() and target.stat().st_size > 1000:
        return target
    url = f'https://files.rcsb.org/download/{pid}.cif'
    contexts = [ssl.create_default_context()]
    try:
        import certifi
        contexts.append(ssl.create_default_context(cafile=certifi.where()))
    except ImportError:
        pass
    last = None
    for ctx in contexts:
        try:
            req = urllib.request.Request(url, headers={'User-Agent': f'RCV/{__version__} (RNA Conformation Viewer)'})
            with urllib.request.urlopen(req, timeout=30, context=ctx) as r:
                data = r.read(200 * 1024 * 1024)
            break
        except urllib.error.HTTPError as exc:
            raise ValueError(f'{pid} could not be downloaded from the RCSB PDB ({exc.code})')
        except (urllib.error.URLError, OSError) as exc:
            last = exc
    else:
        raise ValueError(f'no connection to the RCSB PDB ({getattr(last, "reason", last)}); fetch needs internet access')
    tmp = target.with_suffix('.part')
    tmp.write_bytes(data)
    tmp.replace(target)
    return target


# Documents shown by the Help menu: repository paths in a source checkout, rcv/docs/<path below docs> in installers.
DOC_FILES = {'plugins': 'docs/PLUGINS.md', 'guide': 'docs/RCV.md', 'changelog': 'CHANGELOG.md'}


def doc_file(rel):
    src = P.source_checkout()
    sub = rel[5:] if rel.startswith('docs/') else rel
    for f in ((src / rel) if src else None, P.resource_dir() / 'docs' / sub):
        if f and f.is_file():
            return f
    return None


def system_info(state):
    return {'app': APP_NAME, 'version': __version__, 'cpus': os.cpu_count() or 1,
            'platform': f'{platform.system()} {platform.release()}', 'machine': platform.machine(),
            'python': platform.python_version(), 'vienna': vienna_version(), 'frozen': P.frozen(),
            'shell': state.shell, 'projects_dir': str(state.projects_dir)}


# ----------------------------------------------------------------------------- HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = f'{APP_NAME}/{__version__}'
    protocol_version = 'HTTP/1.1'

    @property
    def state(self) -> AppState:
        return self.server.state

    def log_message(self, fmt, *args):
        if self.state.verbose:
            sys.stderr.write('[app] ' + (fmt % args) + '\n')

    def _host_ok(self):
        host = (self.headers.get('Host') or '').rsplit(':', 1)[0].strip('[]').lower()
        return host in ('localhost', '127.0.0.1', '::1')

    def _send(self, code, body=b'', ctype='application/json', extra=None):
        if isinstance(body, str):
            body = body.encode('utf-8')
        gz = len(body) > 4096 and 'gzip' in (self.headers.get('Accept-Encoding') or '') and \
            ctype.startswith(('application/json', 'text/'))
        if gz:
            body = gzip.compress(body, 5)
        self.send_response(code)
        self.send_header('Content-Type', ctype + ('; charset=utf-8' if ctype.startswith(('text/', 'application/json')) else ''))
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Referrer-Policy', 'no-referrer')
        if gz:
            self.send_header('Content-Encoding', 'gzip')
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != 'HEAD':
            self.wfile.write(body)

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj, separators=(',', ':'), allow_nan=False))

    def _error(self, code, msg):
        self._json({'error': msg}, code)

    def _consume_body(self):
        """Read the request body once, before anything else: a body left unread on a keep-alive connection would be
        parsed as the start of the next request (e.g. '{}GET ...')."""
        try:
            n = int(self.headers.get('Content-Length') or 0)
        except ValueError:
            n = 0
        if n > MAX_BODY:
            self.close_connection = True
            self._body_cache = None
            return False
        self._body_cache = self.rfile.read(n) if n > 0 else b''
        return True

    def _raw_body(self):
        raw = getattr(self, '_body_cache', None)
        if raw is None:
            raise ValueError('request too large')
        if (self.headers.get('Content-Encoding') or '') == 'gzip':
            raw = gzip.decompress(raw)
        return raw

    def _rcv_bytes(self):
        """An .rcv file sent in the body, or named by X-Path (a file the user chose in a native Open dialog)."""
        path = unquote(self.headers.get('X-Path', ''))
        if path:
            f = Path(path).expanduser()
            if f.suffix.lower() != '.rcv' or not f.is_file():
                raise ValueError('not an .rcv file')
            return f.read_bytes()
        return self._raw_body()

    def _body(self):
        raw = self._raw_body()
        return json.loads(raw or b'{}')

    def _static(self, rel):
        if rel in ('', '/'):
            rel = 'index.html'
        rel = unquote(rel).lstrip('/')
        dev = os.environ.get('RCV_DEVDIR') or os.environ.get('RNACONFSWITCH_DEVDIR')
        if rel.startswith('_dev/') and dev:                     # developer test harness (tools/viewer_harness)
            base = Path(dev).resolve()
            target = (base / rel[5:]).resolve()
            if base not in target.parents or not target.is_file():
                return self._error(404, 'not found')
            return self._send(200, target.read_bytes(), mimetypes.guess_type(str(target))[0] or 'application/octet-stream')
        static = P.static_dir().resolve()
        if rel.startswith('user-plugins/'):
            parts = rel.split('/')
            if len(parts) != 3 or not parts[1].isdigit() or int(parts[1]) >= len(self.state.plugin_dirs):
                return self._error(404, 'not found')
            base = self.state.plugin_dirs[int(parts[1])]
            target = (base / parts[2]).resolve()
            if target.parent != base.resolve() or target.suffix != '.js':
                return self._error(404, 'not found')
        else:
            target = (static / rel).resolve()
            if static not in target.parents:
                return self._error(404, 'not found')
        if not target.is_file():
            return self._error(404, 'not found')
        ctype = mimetypes.guess_type(str(target))[0] or 'application/octet-stream'
        extra = {}
        if target.name == 'index.html':
            extra['Content-Security-Policy'] = ("default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
                                                "img-src 'self' data: blob:; media-src 'self' blob:; connect-src 'self'; "
                                                "worker-src 'self' blob:; object-src 'none'; base-uri 'none'")
        self._send(200, target.read_bytes(), ctype, extra)

    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        self._dispatch('GET')

    def do_POST(self):
        self._dispatch('POST')

    def do_DELETE(self):
        self._dispatch('DELETE')

    def _dispatch(self, method):
        if not self._consume_body():
            return self._error(413, 'request too large')
        if not self._host_ok():
            return self._error(403, 'forbidden host')
        url = urlparse(self.path)
        path = url.path
        if not path.startswith('/api/'):
            if method not in ('GET', 'HEAD'):
                return self._error(405, 'method not allowed')
            return self._static(path)
        query = parse_qs(url.query)
        token = self.headers.get('X-Viewer-Token', '')
        if not token and path == '/api/bye':           # sendBeacon cannot set headers
            token = (query.get('token') or [''])[0]
        if not secrets.compare_digest(token, self.state.token):
            return self._error(401, 'missing or invalid session token')
        self.state.last_seen = time.time()
        if path != '/api/bye':
            self.state.closing_at = None
        try:
            self._api(method, path[5:], query)
        except (ValueError, FileNotFoundError, KeyError, PermissionError, RuntimeError) as exc:
            msg = str(exc) if not isinstance(exc, KeyError) else f'not found: {exc}'
            self._error(400, msg)
        except Exception as exc:
            if self.state.verbose:
                traceback.print_exc()
            self._error(500, f'{type(exc).__name__}: {exc}')

    def _api(self, method, route, query):
        st = self.state
        parts = [p for p in route.split('/') if p]
        if method == 'GET' and parts == ['info']:
            return self._json({
                'app': APP_NAME, 'version': __version__,
                'capabilities': {'vienna': vienna_version() is not None, 'predict': True, 'browse': True,
                                 'native_dialogs': st.shell == 'webview'},
                'system': system_info(st),
                'documents': [entry_summary(e) for e in st.entries.values()],
                'jobs': [j.info() for j in st.jobs.values()],
                'plugins': st.js_plugins(),
                'python_plugins': [{'id': p.id, 'name': p.name, 'path': p.path} for p in st.py_plugins],
                'plugin_errors': st.plugin_errors,
                'plugin_dirs': [str(p) for p in st.plugin_dirs],
                'demo': bool(P.demo_dir()), 'website': WEBSITE})
        if parts and parts[0] == 'documents':
            if method == 'GET' and len(parts) == 1:
                return self._json([summary(e.doc) for e in st.entries.values()])
            if len(parts) == 3 and parts[2] == 'methods' and method == 'GET':
                entry = st.entries.get(parts[1])
                if entry is None:
                    return self._error(404, 'unknown document')
                from .methods import citation_text
                return self._json({'methods': methods_text(entry.doc, vienna_version()), 'citation': citation_text()})
            if len(parts) == 2:
                entry = st.entries.get(parts[1])
                if entry is None:
                    return self._error(404, 'unknown document')
                if method == 'GET':
                    return self._json(entry.doc)
                if method == 'DELETE':
                    with st.lock:
                        st.entries.pop(parts[1], None)
                    return self._json({'ok': True})
        if method == 'POST' and parts == ['open']:
            body = self._body()
            if not body.get('paths'):
                raise ValueError('no paths given')
            return self._json(st.open_paths(body['paths'], body.get('chain') or None, body.get('name') or None))
        if method == 'POST' and parts == ['open-demo']:
            demo = P.demo_dir()
            if not demo:
                raise ValueError('no demo data in this installation')
            out = []
            for name in ('HIV1_TAR_NMR_1ANR.pdb', 'HIV1_TAR_prediction', '1EHZ_tRNA_prediction'):
                if (demo / name).exists():
                    out.append(st.open_paths([demo / name]))
            return self._json(out)
        if method == 'POST' and parts == ['rcv', 'pack']:
            # current ensemble (+ viewer state, thumbnail) -> .rcv bytes; the window saves them where the user chooses
            body = self._body()
            entry = st.entries.get(body.get('document') or '')
            doc = entry.doc if entry else body.get('doc')
            if not doc:
                raise ValueError('unknown ensemble')
            thumb = body.get('thumbnail')
            if thumb and thumb.startswith('data:'):
                thumb = thumb.split(',', 1)[1]
            import base64 as _b64
            data = st.pack_rcv(doc, session=body.get('session'), thumbnail=_b64.b64decode(thumb) if thumb else None,
                               entry=entry)
            stats = {'bytes': len(data), 'pdb_bytes': rcvfile.pdb_size_estimate(len(doc['conformations']),
                                                                                len(doc['atoms']['names']))}
            return self._send(200, data, rcvfile.MIME, {'X-RCV-Stats': json.dumps(stats),
                                                        'Access-Control-Expose-Headers': 'X-RCV-Stats'})
        if method == 'POST' and parts == ['rcv', 'open']:
            data = self._raw_body()
            name = unquote(self.headers.get('X-Name', '')) or None
            out = st.open_rcv(data)
            if name:
                out['file'] = name
            return self._json(out)
        if method == 'POST' and parts == ['rcv', 'info']:
            # Tools > Inspect: describe an .rcv file without opening it
            doc, session, info = rcvfile.unpack(self._rcv_bytes())
            key = doc.get('default_metric')
            metric = next((m for m in doc.get('metrics', []) if m['key'] == key), {})
            return self._json({
                'name': doc.get('name'), 'sequence': doc.get('sequence'), 'n_conformations': len(doc['conformations']),
                'n_atoms': len(doc['atoms']['names']), 'origin': doc.get('origin'), **{k: info[k] for k in
                ('format', 'saved', 'app_version', 'bytes', 'pdb_bytes_estimate', 'thumbnail', 'methods', 'contents')},
                'has_session': bool(session), 'metric': {'key': key, 'label': metric.get('label', key), 'unit': metric.get('unit', '')},
                'conformations': [{'id': c['id'], 'label': c.get('label'), 'dbn': c.get('dbn'), 'value': c['values'].get(key)}
                                  for c in doc['conformations']],
                'provenance': doc.get('provenance'), 'source': doc.get('source')})
        if method == 'POST' and parts == ['rcv', 'extract']:
            # Tools > Extract: .rcv -> one multi-model PDB, or a ZIP with one PDB per conformation
            import io as _io
            import zipfile as _zip
            doc, _, _ = rcvfile.unpack(self._rcv_bytes())
            if (self.headers.get('X-Mode') or 'zip') == 'multi':
                return self._send(200, rcvfile.pdb_text(doc), 'chemical/x-pdb')
            buf = _io.BytesIO()
            with _zip.ZipFile(buf, 'w', _zip.ZIP_DEFLATED) as z:
                for k, c in enumerate(doc['conformations']):
                    z.writestr(f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', c['id'])}.pdb", rcvfile.pdb_text(doc, k))
            return self._send(200, buf.getvalue(), 'application/zip')
        if method == 'POST' and parts == ['pdf-font']:
            # subsets of the bundled figure font (Liberation Sans, OFL) for vector PDF figures; 'available': false = Helvetica
            from .fonts import pdf_font
            body = self._body()
            reg = pdf_font(str(body.get('regular', ''))[:20000], False)
            bold = pdf_font(str(body.get('bold', ''))[:20000], True) if body.get('bold') else None
            return self._json({'available': reg is not None, 'regular': reg, 'bold': bold})
        if method == 'POST' and parts == ['methods']:
            from .methods import citation_text
            doc = self._body().get('doc') or {}
            return self._json({'methods': methods_text(doc, vienna_version()), 'citation': citation_text()})
        if method == 'POST' and parts == ['upload']:
            body = self._body()
            if not body.get('files'):
                raise ValueError('no files uploaded')
            return self._json(st.open_uploads(body['files'], body.get('chain') or None, body.get('name') or None))
        if method == 'GET' and parts == ['library']:
            now = time.time()
            cache = getattr(st, '_library_cache', None)
            if not cache or now - cache[0] > 20 or (query.get('refresh') or [''])[0]:
                cache = (now, find_packages(st.library_roots))
                st._library_cache = cache
            return self._json({'roots': [str(r) for r in st.library_roots], 'packages': cache[1]})
        if method == 'GET' and parts == ['browse']:
            return self._json(browse((query.get('path') or [''])[0]))
        if method == 'POST' and parts == ['settings']:
            body = self._body()
            if body.get('projects_dir'):
                d = Path(body['projects_dir']).expanduser()
                d.mkdir(parents=True, exist_ok=True)
                st.projects_dir = d
                if d not in st.library_roots:
                    st.library_roots.insert(0, d)
                P.save_settings({'projects_dir': str(d)})
            return self._json(system_info(st))
        if method == 'POST' and parts in (['predict'], ['explore'], ['batch']):
            job = {'predict': PredictJob, 'explore': ExploreJob, 'batch': BatchJob}[parts[0]](st, self._body())
            with st.lock:
                st.jobs[job.id] = job
            job.start()
            return self._json(job.info())
        if parts and parts[0] == 'jobs':
            if method == 'GET' and len(parts) == 1:
                return self._json([j.info() for j in st.jobs.values()])
            job = st.jobs.get(parts[1]) if len(parts) >= 2 else None
            if job is None:
                return self._error(404, 'unknown job')
            if method == 'GET' and len(parts) == 2:
                return self._json(job.info())
            if method == 'POST' and len(parts) == 3 and parts[2] == 'cancel':
                job.cancel()
                return self._json(job.info())
        if method == 'POST' and parts == ['write-file']:
            # desktop shell: write bytes to a path the user picked in a native Save dialog
            target = unquote(self.headers.get('X-Path', ''))
            if not target or target not in st.approved_paths:
                raise PermissionError('this path was not chosen in a save dialog')
            data = self._raw_body()
            Path(target).write_bytes(data)
            st.approved_paths.discard(target)
            return self._json({'path': target, 'bytes': len(data)})
        if method == 'POST' and parts == ['read-json']:
            # sessions chosen in a native Open dialog (desktop shell); JSON files only
            target = Path(self._body().get('path', ''))
            if target.suffix.lower() != '.json' or not target.is_file() or target.stat().st_size > 300 * 1024 * 1024:
                raise ValueError('not a readable .json file')
            return self._json({'path': str(target), 'text': target.read_text(encoding='utf-8')})
        if method == 'POST' and parts == ['export-html']:
            from .htmlexport import build_html
            body = self._body()
            entry = st.entries.get(body.get('document') or '')
            doc = entry.doc if entry else body.get('doc')
            if not doc:
                raise ValueError('unknown ensemble')
            page = build_html(doc, body.get('view'), body.get('title'), methods_text(doc, vienna_version()))
            return self._send(200, page, 'text/html')
        if method == 'POST' and parts == ['switch']:
            # 2D-level switch analysis (ViennaRNA): state populations and the single-mutation scan between two states
            if vienna_version() is None:
                raise ValueError('the switch analysis needs ViennaRNA (pip install ViennaRNA)')
            from .switch import analyse, ensemble_states
            body = self._body()
            entry = st.entries.get(body.get('document') or '')
            doc = entry.doc if entry else body.get('doc')
            if not doc:
                raise ValueError('unknown ensemble')
            states = ensemble_states(doc)
            for extra in body.get('extra_states') or []:
                if isinstance(extra, str) and len(extra) == len(doc['sequence']) and extra not in [s['dbn'] for s in states]:
                    states.append({'dbn': extra, 'conformations': [], 'source': 'added'})
            a, b = body.get('a'), body.get('b')
            res = analyse(doc['sequence'], states, float(body.get('temperature', 37.0)),
                          int(a) if a is not None else None, int(b) if b is not None else None)
            return self._json(res)
        if method == 'POST' and parts == ['plugins', 'new']:
            # Plugins > New plugin: write a commented template into the user plugin folder
            body = self._body()
            kind = 'py' if body.get('kind') == 'py' else 'js'
            name = re.sub(r'[^A-Za-z0-9 _-]+', '', str(body.get('name') or 'My plugin')).strip()[:40] or 'My plugin'
            pid = re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-') or 'my-plugin'
            folder = st.plugin_dirs[0]
            folder.mkdir(parents=True, exist_ok=True)
            target = folder / f'{pid}.{kind}'
            n = 2
            while target.exists():
                target = folder / f'{pid}-{n}.{kind}'
                n += 1
            tpl = (P.resource_dir() / 'templates' / f'plugin_template.{kind}').read_text(encoding='utf-8')
            target.write_text(tpl.replace('__NAME__', name).replace('__ID__', target.stem), encoding='utf-8')
            return self._json({'path': str(target), 'folder': str(folder)})
        if method == 'GET' and len(parts) == 2 and parts[0] == 'docs' and parts[1] in DOC_FILES:
            # Help > User guide, What is new, Plugin developer guide (Markdown shown inside the app)
            f = doc_file(DOC_FILES[parts[1]])
            if f:
                return self._send(200, f.read_text(encoding='utf-8'), 'text/markdown')
            return self._error(404, 'this document is not installed')
        if method == 'GET' and len(parts) == 3 and parts[:2] == ['docs', 'images']:
            name = parts[2]
            f = doc_file('docs/images/' + name) if re.fullmatch(r'[A-Za-z0-9_.-]+\.(png|svg|jpg|webp)', name) else None
            if not f:
                return self._error(404, 'no such image')
            ctype = {'.png': 'image/png', '.svg': 'image/svg+xml', '.jpg': 'image/jpeg', '.webp': 'image/webp'}[f.suffix.lower()]
            return self._send(200, f.read_bytes(), ctype)
        if method == 'POST' and parts == ['read-text']:
            # command scripts for the console (run <file>); plain text files only
            target = Path(self._body().get('path', '')).expanduser()
            if target.suffix.lower() not in ('.txt', '.rcvs', '.cmd') or not target.is_file() or target.stat().st_size > 2 * 1024 * 1024:
                raise ValueError('scripts must be .rcvs, .cmd or .txt files under 2 MB')
            return self._json({'path': str(target), 'text': target.read_text(encoding='utf-8', errors='replace')})
        if method == 'POST' and parts == ['fetch']:
            body = self._body()
            path = fetch_pdb(body.get('id', ''), st.projects_dir / 'fetched')
            return self._json(st.open_paths([path], body.get('chain') or None))
        if method == 'POST' and parts == ['client-log']:
            body = self._body()
            msg = str(body.get('message', ''))[:2000].replace('\n', ' | ')
            sys.stderr.write(f"[client {str(body.get('level', 'info'))[:8]}] {msg}\n")
            sys.stderr.flush()
            return self._json({'ok': True})
        if method == 'POST' and parts == ['bye']:
            # the window was closed (or reloaded): quit a few seconds later unless it reconnects
            st.closing_at = time.time() + 4.0
            return self._json({'ok': True})
        return self._error(404, f'unknown endpoint {method} /api/{route}')


class AppServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, addr, state):
        super().__init__(addr, Handler)
        self.state = state


def make_server(port=8765, verbose=False, projects_dir=None):
    projects = Path(projects_dir) if projects_dir else P.default_projects_dir()
    roots = [projects]
    if not projects_dir:                                 # predictions made with earlier versions (RNAconfSwitch)
        roots += [d for d in P.legacy_projects_dirs() if d != projects]
    if P.demo_dir():
        roots.append(P.demo_dir())
    if P.source_checkout():
        roots.append(P.source_checkout() / 'runs')      # not the whole checkout: benchmark runs are slow to scan
    state = AppState(projects, roots, P.plugin_dirs(), verbose)
    server = None
    for p in ([port] if port == 0 else [port] + list(range(port + 1, port + 20)) + [0]):
        try:
            server = AppServer(('127.0.0.1', p), state)
            break
        except OSError:
            continue
    if server is None:
        raise RuntimeError('could not open a local port')
    return server, state
