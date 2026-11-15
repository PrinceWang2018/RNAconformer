"""RCV 0.7.5 workflow: project library, recipes and the plain-language console."""
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from rcv import library as L
from rcv import recipes as R
from rcv import server as srv

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'


class _Home(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.env = mock.patch.dict(os.environ, {'RCV_HOME': str(Path(self.tmp.name) / 'home')})
        self.env.start()

    def tearDown(self):
        self.env.stop()
        self.tmp.cleanup()


class TestRecipes(_Home):
    def test_builtin_recipes_are_valid(self):
        names = [r['name'] for r in R.list_recipes()]
        self.assertIn('Riboswitch report', names)
        for r in R.BUILTIN:
            R.validate(r)
            self.assertTrue(r['id'].startswith('builtin-'))

    def test_validate(self):
        r = R.validate({'name': ' Lab view ', 'steps': 'color helix\n\n  tab switch ', 'prediction': {'models': '99', 'sampling': 'fast', 'templates': 'no'}})
        self.assertEqual(r['name'], 'Lab view')
        self.assertEqual(r['steps'], ['color helix', 'tab switch'])
        self.assertEqual(r['prediction'], {'models': 50, 'sampling': 'fast', 'templates': 'off'})
        for bad in ({}, {'name': 'x'}, {'name': 'x', 'steps': ['a'], 'prediction': {'sampling': 'slow'}}, [], {'name': 'x', 'steps': ['a'], 'rcv_recipe': 9}):
            with self.assertRaises((ValueError, TypeError)):
                R.validate(bad)

    def test_save_import_delete(self):
        a = R.save({'name': 'My TAR view', 'steps': ['color switching']})
        self.assertTrue(Path(a['path']).is_file() and a['path'].endswith('.rcvrecipe'))
        b = R.save({'name': 'My TAR view', 'steps': ['color helix']})          # same name -> a second file
        self.assertNotEqual(a['id'], b['id'])
        a2 = R.save({'name': 'My TAR view', 'steps': ['color pairing']}, rid=a['id'])
        self.assertEqual(a2['id'], a['id'])
        self.assertEqual(R.get(a['id'])['steps'], ['color pairing'])
        c = R.save({'name': 'Riboswitch report', 'steps': ['tab switch']}, rid='builtin-riboswitch-report')
        self.assertFalse(c['id'].startswith('builtin-'))                       # a changed built-in becomes the user's own
        shared = json.dumps({'rcv_recipe': 1, 'name': 'Shared', 'description': 'from a colleague', 'steps': ['show putty']})
        d = R.import_text(shared)
        self.assertEqual(R.get('shared')['description'], 'from a colleague')
        with self.assertRaises(ValueError):
            R.import_text('not json')
        R.delete(d['id'])
        with self.assertRaises(KeyError):
            R.get(d['id'])
        with self.assertRaises(ValueError):
            R.delete('builtin-riboswitch-report')
        with self.assertRaises(ValueError):
            R.delete('../evil')
        # a broken file is listed with its error, not hidden and not fatal
        (R.recipes_dir() / 'broken.rcvrecipe').write_text('{', encoding='utf-8')
        self.assertTrue(any(r['id'] == 'broken' and r.get('error') for r in R.list_recipes()))


class TestLibrary(_Home):
    def setUp(self):
        super().setUp()
        self.work = Path(self.tmp.name) / 'work'
        shutil.copytree(DEMO / 'HIV1_TAR_prediction', self.work / 'runs' / 'HIV1_TAR_prediction')
        shutil.copy(DEMO / 'HIV1_TAR_NMR_1ANR.pdb', self.work / 'tar.pdb')

    def test_opened_documents_are_recorded(self):
        st = srv.AppState(self.work / 'projects', [], [])
        run = st.open_paths([self.work / 'runs' / 'HIV1_TAR_prediction'])
        st.open_paths([self.work / 'tar.pdb'])
        st.open_paths([self.work / 'tar.pdb'])
        idx = L.items([], srv.find_packages)
        by_kind = {v['kind']: v for v in idx['items']}
        self.assertEqual(set(by_kind), {'run', 'structure'})
        r = by_kind['run']
        self.assertEqual((r['name'], r['length'], r['n_conformations']), (run['name'], 29, run['n_conformations']))
        self.assertTrue(r['params'].get('engine', '').startswith('rnafold3d'))
        self.assertEqual(r['params']['models'], 10)
        self.assertEqual(by_kind['structure']['opens'], 2)
        self.assertTrue(all(v['exists'] for v in idx['items']))
        # thumbnail (JPEG only) for a document, then forget
        L.save_thumbnail(r['path'], 'data:image/jpeg;base64,/9j/4AAQSkZJRgABAQ==')
        self.assertTrue(next(v for v in L.items([], srv.find_packages)['items'] if v['kind'] == 'run')['thumb'])
        self.assertTrue(L.thumbnail(L.key_of(r['path'])).startswith(b'\xff\xd8'))
        with self.assertRaises(ValueError):
            L.save_thumbnail(r['path'], 'data:image/png;base64,iVBORw0KGgo=')
        with self.assertRaises(KeyError):
            L.thumbnail('../../etc')
        L.forget(r['path'])
        self.assertEqual([v['kind'] for v in L.items([], srv.find_packages)['items']], ['structure'])

    def test_watched_folders_and_missing_files(self):
        idx = L.items([self.work / 'runs'], srv.find_packages)
        self.assertEqual([(v['kind'], v.get('watched')) for v in idx['items']], [('run', True)])
        extra = self.work / 'shared'
        (extra / 'sub').mkdir(parents=True)
        (extra / 'sub' / 'old.rcv').write_bytes(b'PK')
        L.set_folder(extra)
        self.assertEqual(L.folders(), [str(extra.resolve())])
        kinds = sorted(v['kind'] for v in L.items([self.work / 'runs'], srv.find_packages)['items'])
        self.assertEqual(kinds, ['rcv', 'run'])
        with self.assertRaises(ValueError):
            L.set_folder(self.work / 'nope')
        L.set_folder(extra, add=False)
        self.assertEqual(L.folders(), [])
        # a recorded file that was deleted stays listed, marked missing
        st = srv.AppState(self.work / 'projects', [], [])
        st.open_paths([self.work / 'tar.pdb'])
        (self.work / 'tar.pdb').unlink()
        gone = [v for v in L.items([], srv.find_packages)['items'] if v['kind'] == 'structure']
        self.assertEqual([v['exists'] for v in gone], [False])

    def test_run_and_its_cached_rcv_are_one_item(self):
        doc = {'name': 'x', 'sequence': 'GGAA', 'conformations': [{}],
               'source': {'kind': 'package', 'path': str(self.work / 'runs' / 'HIV1_TAR_prediction'),
                          'rcv': str(self.work / 'runs' / 'HIV1_TAR_prediction' / 'x.rcv')}}
        self.assertEqual(L.describe_doc(doc)['kind'], 'run')
        self.assertIsNone(L.describe_doc({'name': 'up', 'source': {'kind': 'upload'}}))


class TestPlainLanguage(unittest.TestCase):
    @unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
    def test_translations(self):
        r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'plain_unit.mjs')], capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout[-3000:] + r.stderr[-3000:])
        self.assertIn('all plain-language tests passed', r.stdout)

    def test_console_falls_back_to_plain_language(self):
        js = (ROOT / 'src/rcv/static/js/console/panel.js').read_text(encoding='utf-8')
        self.assertIn('translatePlain(src', js)
        cmds = (ROOT / 'src/rcv/static/js/console/commands.js').read_text(encoding='utf-8')
        self.assertIn('takes no arguments', cmds)        # "reset the colours" must not reset the view


if __name__ == '__main__':
    unittest.main()
