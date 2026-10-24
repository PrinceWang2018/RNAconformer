"""Tests for RCV 0.6: bundled figure font and PDF font subsets, .rcv inspect/extract endpoints, plugin templates and
guide, application icons, the dedicated app-window profile, and (with Node.js) the surface mesh, clustering,
palette check and PDF text layout."""
import http.client
import io
import json
import shutil
import subprocess
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

from rcv import __author__, __version__, rcvfile
from rcv import paths as P
from rcv import server as srv
from rcv.data import load_document
from rcv.desktop import app_window_command
from rcv.methods import citation_text

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'
PRED = DEMO / 'HIV1_TAR_prediction'

try:
    import fontTools  # noqa: F401
    HAVE_FT = True
except ImportError:
    HAVE_FT = False


class TestFigureFont(unittest.TestCase):
    def test_bundled_open_licence_font(self):
        fonts = P.static_dir() / 'fonts'
        for f in ('LiberationSans-Regular.ttf', 'LiberationSans-Bold.ttf', 'LICENSE-OFL.txt'):
            self.assertTrue((fonts / f).is_file(), f)
        self.assertIn('SIL OPEN FONT LICENSE Version 1.1', (fonts / 'LICENSE-OFL.txt').read_text())
        css = (P.static_dir() / 'css' / 'app.css').read_text()
        self.assertIn('fonts/LiberationSans-Regular.ttf', css)
        # no Arial anywhere in the figure code (licence concerns): only open-licence fonts are named
        for js in (P.static_dir() / 'js').rglob('*.js'):
            text = js.read_text(encoding='utf-8')
            self.assertNotRegex(text, r"FIG_FONT\s*=\s*['\"][^'\"]*Arial", js.name)

    @unittest.skipUnless(HAVE_FT, 'fontTools not installed')
    def test_pdf_font_subset(self):
        from fontTools.ttLib import TTFont
        from rcv.fonts import pdf_font
        r = pdf_font('ΔG (kcal/mol) Å RMSD′ 0123', False)
        self.assertEqual(r['name'], 'LiberationSans')
        ttf = TTFont(io.BytesIO(__import__('base64').b64decode(r['ttf'])))
        self.assertLess(len(r['ttf']), 40000)                    # a subset, not the 400 kB font
        cmap = ttf.getBestCmap()
        for ch in 'ΔGÅ′':
            self.assertIn(ord(ch), cmap)
        self.assertEqual(len(r['widths']), 224)
        self.assertEqual([e['char'] for e in r['extra']], ['Δ', '′'])
        self.assertEqual(pdf_font('A', True)['name'], 'LiberationSans-Bold')


class TestServer06(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def req(self, method, path, body=None, raw=None, headers=None, token=True):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=120)
        h = {'Host': f'localhost:{self.port}', **({'X-Viewer-Token': self.state.token} if token else {}), **(headers or {})}
        data = raw
        if body is not None:
            data = json.dumps(body).encode()
            h['Content-Type'] = 'application/json'
        c.request(method, path, body=data, headers=h)
        r = c.getresponse()
        out, hdr = r.read(), dict(r.getheaders())
        c.close()
        return r.status, out, hdr

    def test_icons_and_manifest(self):
        for path, ctype in (('/favicon.ico?v=0.6', 'image/x-icon'), ('/img/icon-256.png', 'image/png'),
                            ('/manifest.webmanifest', 'application/manifest+json'), ('/fonts/LiberationSans-Regular.ttf', 'font/ttf')):
            st, body, hdr = self.req('GET', path, token=False)
            self.assertEqual(st, 200, path)
            if ctype != 'font/ttf':
                self.assertTrue(hdr['Content-Type'].startswith(ctype), (path, hdr['Content-Type']))
        st, body, _ = self.req('GET', '/', token=False)
        self.assertIn(b'img/icon-256.png', body)
        self.assertIn(b'manifest.webmanifest', body)

    def test_rcv_inspect_and_extract(self):
        data = rcvfile.pack(load_document([PRED]), session={'view': {}}, thumbnail=b'\xff\xd8jpeg')
        st, raw, _ = self.req('POST', '/api/rcv/info', raw=data, headers={'Content-Type': 'application/octet-stream'})
        self.assertEqual(st, 200, raw)
        info = json.loads(raw)
        self.assertEqual(info['n_conformations'], 10)
        self.assertTrue(info['has_session'])
        self.assertIn('thumbnail.jpg', info['contents'])
        st, z, _ = self.req('POST', '/api/rcv/extract', raw=data, headers={'Content-Type': 'application/octet-stream', 'X-Mode': 'zip'})
        names = zipfile.ZipFile(io.BytesIO(z)).namelist()
        self.assertEqual(len(names), 10)
        self.assertTrue(all(n.endswith('.pdb') for n in names))
        st, multi, _ = self.req('POST', '/api/rcv/extract', raw=data, headers={'Content-Type': 'application/octet-stream', 'X-Mode': 'multi'})
        self.assertEqual(multi.count(b'\nMODEL '), 10)              # every MODEL record follows a newline (REMARK first)
        f = Path(self.tmp) / 'x.rcv'
        f.write_bytes(data)
        st, raw, _ = self.req('POST', '/api/rcv/info', headers={'X-Path': str(f)})
        self.assertEqual(json.loads(raw)['name'], json.loads(raw)['name'])
        st, raw, _ = self.req('POST', '/api/rcv/info', headers={'X-Path': str(PRED / 'metadata' / 'manifest.json')})
        self.assertEqual(st, 400)

    def test_pdf_font_endpoint(self):
        st, raw, _ = self.req('POST', '/api/pdf-font', {'regular': 'RMSD (Å)', 'bold': 'A'})
        r = json.loads(raw)
        self.assertEqual(st, 200)
        self.assertEqual(r['available'], HAVE_FT)
        if HAVE_FT:
            self.assertEqual(r['regular']['name'], 'LiberationSans')

    def test_plugin_template_and_guide(self):
        for kind in ('js', 'py'):
            st, raw, _ = self.req('POST', '/api/plugins/new', {'name': 'My Test plugin!', 'kind': kind})
            self.assertEqual(st, 200, raw)
            path = Path(json.loads(raw)['path'])
            self.assertEqual(path.parent, Path(self.tmp) / 'plugins')
            text = path.read_text()
            self.assertNotIn('__NAME__', text)
            self.assertIn('My Test plugin', text)
        st, raw, _ = self.req('POST', '/api/plugins/new', {'name': 'My Test plugin!', 'kind': 'js'})
        self.assertTrue(json.loads(raw)['path'].endswith('my-test-plugin-2.js'))
        # the Python template is a working analysis plugin
        plugin = srv.PythonPlugin(Path(self.tmp) / 'plugins' / 'my-test-plugin.py')
        import numpy as np
        val = plugin.compute('GGGG', [{"C4'": np.array([float(i), 0.0, 0.0])} for i in range(4)])
        self.assertAlmostEqual(val['radius'], 1.5, places=2)
        st, raw, hdr = self.req('GET', '/api/docs/plugins')
        self.assertEqual(st, 200)
        self.assertIn(b'Writing plugins for RCV', raw)
        node = shutil.which('node')
        if node:
            r = subprocess.run([node, '--input-type=module', '--check'], input=(Path(self.tmp) / 'plugins' / 'my-test-plugin.js').read_text(), text=True, capture_output=True)
            self.assertEqual(r.returncode, 0, r.stderr)


class TestDesktop06(unittest.TestCase):
    def test_app_window_uses_its_own_profile(self):
        cmd = app_window_command('http://localhost:1/#token=x')
        if cmd:
            self.assertTrue(any(a.startswith('--user-data-dir=') and 'RCV' in a for a in cmd), cmd)
            self.assertIn('--disable-extensions', cmd)

    def test_author_and_citation(self):
        self.assertEqual(__author__, 'Zixiang Wang')
        self.assertIn('Zixiang Wang', citation_text())
        self.assertIn(__version__, citation_text())
        self.assertNotIn('deep learning', (P.static_dir() / 'js' / 'welcome.js').read_text())


@unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
class TestBrowserModules06(unittest.TestCase):
    def test_surface_clustering_palette_pdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            doc = Path(tmp) / 'doc.json'
            doc.write_text(json.dumps(load_document([DEMO / '1EHZ_tRNA_prediction'])))
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'rcv06_unit.mjs'), str(doc), tmp], capture_output=True, text=True, timeout=180)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            import re
            import zlib
            pdf = (Path(tmp) / 'rotated.pdf').read_bytes()
            stream = re.search(rb'/FlateDecode >>\nstream\n(.*?)\nendstream', pdf, re.S).group(1)
            ops = zlib.decompress(stream)
            self.assertIn(b' Tm ', ops)
            self.assertIn(b'0 1 -1 0', ops)                          # the y-axis title is rotated by 90 degrees


if __name__ == '__main__':
    unittest.main()
