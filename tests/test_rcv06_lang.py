"""RCV 0.6: interface languages, window layout and the in-app help documents."""
import http.client
import json
import re
import shutil
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path

from rcv import server as srv

ROOT = Path(__file__).resolve().parents[1]
STATIC = ROOT / 'src' / 'rcv' / 'static'
I18N = STATIC / 'i18n'
PLACEHOLDER = re.compile(r'\{\d\}')


def load(code):
    return json.loads((I18N / f'{code}.json').read_text(encoding='utf-8'))


class TestTranslationFiles(unittest.TestCase):
    def test_language_list_matches_files(self):
        js = (STATIC / 'js' / 'i18n.js').read_text(encoding='utf-8')
        codes = re.findall(r"\{ code: '([A-Za-z-]+)'", js)
        self.assertEqual(codes[0], 'en')                       # English is the source language, no file
        self.assertEqual(sorted(codes[1:]), sorted(p.stem for p in I18N.glob('*.json')))
        self.assertGreaterEqual(len(codes), 11)

    def test_same_strings_and_placeholders_everywhere(self):
        ref = None
        for f in sorted(I18N.glob('*.json')):
            d = json.loads(f.read_text(encoding='utf-8'))
            self.assertEqual(d['_language'], f.stem)
            s = d['strings']
            if ref is None:
                ref = set(s)
            self.assertEqual(sorted(set(s) ^ ref), [], f'{f.name}: strings differ from the other languages')
            for src, out in s.items():
                self.assertTrue(out.strip(), f'{f.name}: empty translation of {src!r}')
                # a translation may leave out a placeholder (the plural "s" in "ensemble{1}"), never invent one
                self.assertLessEqual(set(PLACEHOLDER.findall(out)), set(PLACEHOLDER.findall(src)),
                                     f'{f.name}: placeholders of {src!r} -> {out!r}')
        self.assertGreater(len(ref), 1000)

    def test_patterns_can_be_split(self):
        # Two placeholders side by side ("ΔE {3}{4}") cannot be split reliably, so the text would stay English.
        # The methods paragraph is meant to stay in English (it is pasted into manuscripts).
        allowed = {'Figures show {0}{1} rendered with RCV ({2} style); residues are coloured by {3}{4}.'}
        bad = [k for k in load('de')['strings'] if re.search(r'\{\d\}\{\d\}', k) and k not in allowed]
        self.assertEqual(bad, [])

    def test_file_dialog_places_are_translated(self):
        # the names the server gives the shortcut folders of the Open dialog
        src = (ROOT / 'src' / 'rcv' / 'server.py').read_text(encoding='utf-8')
        names = set(re.findall(r"\{'name': '([^']+)', 'path'", src)) | {'Windows: {0}'}
        self.assertIn('Projects folder', names)
        for code in ('zh-CN', 'ru'):
            s = load(code)['strings']
            self.assertEqual(sorted(n for n in names if n not in s), [], code)

    @unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
    def test_translation_matcher(self):
        r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'i18n_unit.mjs')], capture_output=True, text=True,
                           timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout[-3000:] + r.stderr[-3000:])
        self.assertIn('all i18n tests passed', r.stdout)


class TestWindowLayout(unittest.TestCase):
    def test_every_part_is_in_the_window_menu(self):
        js = (STATIC / 'js' / 'windows.js').read_text(encoding='utf-8')
        ids = re.findall(r"\{ id: '(\w+)', label: '([^']+)'", js)
        self.assertEqual({i for i, _ in ids} >= {'left', 'right', 'ss2d', 'console', 'seqbar', 'toolbar', 'timeline',
                                                  'statusbar'}, True)
        s = load('zh-CN')['strings']
        self.assertEqual([lab for _, lab in ids if lab not in s], [])          # labels have translations


class TestHelpDocuments(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def get(self, path):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=60)
        c.request('GET', path, headers={'Host': f'localhost:{self.port}', 'X-Viewer-Token': self.state.token})
        r = c.getresponse()
        out = r.status, r.read(), r.getheader('Content-Type')
        c.close()
        return out

    def test_guide_changelog_and_plugin_guide(self):
        for name, first in (('guide', '# '), ('changelog', '# Changelog'), ('plugins', '# ')):
            st, raw, ctype = self.get(f'/api/docs/{name}')
            self.assertEqual(st, 200, name)
            self.assertTrue(ctype.startswith('text/markdown'))
            self.assertTrue(raw.decode('utf-8').startswith(first), name)
        st, _, _ = self.get('/api/docs/secrets')
        self.assertEqual(st, 404)

    def test_guide_images(self):
        md = (ROOT / 'docs' / 'RCV.md').read_text(encoding='utf-8')
        names = sorted(set(re.findall(r'\]\(images/([^)\s]+)\)', md)))
        self.assertTrue(names)
        for n in names:
            st, raw, ctype = self.get(f'/api/docs/images/{n}')
            self.assertEqual(st, 200, n)
            self.assertTrue(ctype.startswith('image/'), n)
        for bad in ('..%2Fserver.py', 'x.py', '.hidden.png'):
            st, _, _ = self.get(f'/api/docs/images/{bad}')
            self.assertEqual(st, 404, bad)


if __name__ == '__main__':
    unittest.main()
