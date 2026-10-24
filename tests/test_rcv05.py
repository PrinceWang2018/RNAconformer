"""Tests for RCV 0.5: the .rcv file format, 2D layouts, batch input, switch analysis, methods text, interactive
HTML export, the new API endpoints and command-line tools, and (with Node.js) the console selection language,
2D layout maths and transition paths."""
import base64
import http.client
import io
import json
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

import numpy as np

from rcv import rcvfile
from rcv import server as srv
from rcv.data import build_document, load_document
from rcv.jobs import BatchJob, parse_fasta
from rcv.methods import methods_text

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'
NMR = DEMO / 'HIV1_TAR_NMR_1ANR.pdb'
PRED = DEMO / 'HIV1_TAR_prediction'

try:
    import RNA  # noqa: F401
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False


class TestRcvFormat(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pred = load_document([PRED])
        cls.nmr = load_document([NMR])

    def test_round_trip_is_exact_to_pdb_precision(self):
        for doc in (self.pred, self.nmr):
            data = rcvfile.pack(doc, session={'view': {'conf': 2}}, thumbnail=b'\x89PNG\r\n', methods='m')
            back, session, info = rcvfile.unpack(data)
            err = np.abs(rcvfile.doc_coords(doc) - rcvfile.doc_coords(back)).max()
            self.assertLessEqual(err, 0.0005 + 1e-6)
            self.assertEqual(back['bonds'], doc['bonds'])
            self.assertEqual({k: v for k, v in back.items() if k != 'conformations'},
                             {k: v for k, v in doc.items() if k != 'conformations'})
            for a, b in zip(doc['conformations'], back['conformations']):
                self.assertEqual({k: v for k, v in a.items() if k != 'coords'}, {k: v for k, v in b.items() if k != 'coords'})
            self.assertEqual(session, {'view': {'conf': 2}})
            self.assertEqual(info['methods'], 'm')
            self.assertTrue(info['thumbnail'])

    def test_smaller_than_pdb_and_self_describing(self):
        data = rcvfile.pack(self.nmr)
        pdb = rcvfile.pdb_size_estimate(len(self.nmr['conformations']), len(self.nmr['atoms']['names']))
        self.assertLess(len(data) * 8, pdb)                       # at least 8x smaller than the PDB models
        z = zipfile.ZipFile(io.BytesIO(data))
        self.assertEqual(z.infolist()[0].filename, 'mimetype')
        self.assertEqual(z.infolist()[0].compress_type, zipfile.ZIP_STORED)
        self.assertEqual(z.read('mimetype'), rcvfile.MIME.encode())
        self.assertIn(b'rcv-coords/1', z.read('README.txt'))
        self.assertTrue(rcvfile.is_rcv(data))
        man = json.loads(z.read('rcv.json'))
        self.assertEqual(man['format'], 'rcv/1')
        self.assertEqual(man['n_conformations'], 20)

    def test_damaged_and_foreign_files_are_rejected(self):
        data = bytearray(rcvfile.pack(self.pred))
        z = zipfile.ZipFile(io.BytesIO(bytes(data)))
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w') as out:
            for info in z.infolist():
                raw = z.read(info.filename)
                if info.filename == 'coords.bin':
                    raw = raw[:-1] + bytes([raw[-1] ^ 0xFF])
                out.writestr(info, raw)
        with self.assertRaisesRegex(ValueError, 'damaged'):
            rcvfile.unpack(buf.getvalue())
        with self.assertRaisesRegex(ValueError, 'not an .rcv'):
            rcvfile.unpack(b'not a zip')
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w') as out:
            out.writestr('rcv.json', json.dumps({'format': 'rcv/9'}))
            out.writestr('coords.bin', b'')
        with self.assertRaisesRegex(ValueError, 'rcv/9'):
            rcvfile.unpack(buf.getvalue())

    def test_rebuild_from_rcv_keeps_energies(self):
        back, _, info = rcvfile.unpack(rcvfile.pack(self.pred, provenance_manifest={'x': 1}))
        confs = rcvfile.conformations_from_doc(back)
        doc = build_document(back['name'], confs, None, set(), back['source'])
        for a, b in zip(self.pred['conformations'], doc['conformations']):
            for k, v in a['values'].items():
                if isinstance(v, float):
                    self.assertAlmostEqual(v, b['values'][k], delta=0.02, msg=k)
        self.assertEqual(info['provenance_manifest'], {'x': 1})

    def test_pdb_text(self):
        text = rcvfile.pdb_text(self.pred, 0)
        self.assertEqual(sum(1 for l in text.splitlines() if l.startswith('ATOM')), len(self.pred['atoms']['names']))
        atoms = [l for l in text.splitlines() if l.startswith('ATOM')]
        self.assertTrue(all(len(l) in (78, 80) and l[76:78].strip() in ('C', 'N', 'O', 'P') for l in atoms))
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in atoms])
        self.assertLess(np.abs(xyz - rcvfile.doc_coords(self.pred)[0]).max(), 0.0006)
        multi = rcvfile.pdb_text(self.pred)
        self.assertEqual(multi.count('MODEL '), len(self.pred['conformations']))


class TestDocumentAdditions(unittest.TestCase):
    def test_layouts_and_checksums(self):
        doc = load_document([NMR])
        self.assertIn(Path(NMR).name, doc['source']['checksums'])
        if HAVE_VIENNA:
            keys = {c['layout_key'] for c in doc['conformations']}
            self.assertTrue(keys <= set(doc['layouts2d']))
            for k, xy in doc['layouts2d'].items():
                self.assertEqual(len(xy), 2 * len(doc['sequence']))
                pts = np.array(xy).reshape(-1, 2)
                self.assertGreater(pts[0, 1], pts[:, 1].mean())          # 5' end drawn at the bottom

    def test_methods_text(self):
        doc = load_document([PRED])
        text = methods_text(doc, '2.7.2')
        for word in ('rnafold3d engine', 'ViennaRNA 2.7.2', 'trajectories', 'C4′', 'not experimental free energies'):
            self.assertIn(word, text)
        text2 = methods_text(load_document([NMR]))
        self.assertIn('HIV1_TAR_NMR_1ANR.pdb', text2)


class TestBatchInput(unittest.TestCase):
    def test_fasta_parsing(self):
        recs = parse_fasta('>a x\nGGGAAACCC\n>b\nggg aaa ccc\n(((...)))\n>a\nACGTACGTAC\n')
        self.assertEqual([r['name'] for r in recs], ['a', 'b', 'a_2'])
        self.assertEqual(recs[1]['structure'], '(((...)))')
        self.assertEqual(recs[2]['sequence'], 'ACGUACGUAC')
        self.assertEqual(parse_fasta('GGGGAAACCCC')[0]['name'], 'seq1')

    def test_batch_validation(self):
        state = srv.AppState(Path(tempfile.mkdtemp()), [], [])
        with self.assertRaisesRegex(ValueError, 'unsupported letters'):
            BatchJob(state, {'fasta': '>a\nGGGGAAACCCC\n>b\nGGGXAAACCCC\n'})
        with self.assertRaisesRegex(ValueError, 'no sequences'):
            BatchJob(state, {'fasta': '>a\n'})
        job = BatchJob(state, {'fasta': '>a\nGGGGAAACCCC\n>b\nGGCACUUCGGUGCC\n', 'models': 2, 'sampling': 'fast', 'cores': 1})
        self.assertEqual([s['label'] for s in job.stages], ['a · 11 nt', 'b · 14 nt'])
        self.assertGreater(job.eta_s, 0)


@unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
class TestSwitchAnalysis(unittest.TestCase):
    def test_states_populations_and_scan(self):
        from rcv.switch import analyse, ensemble_states, compatible
        doc = load_document([PRED])
        states = ensemble_states(doc)
        self.assertEqual(len(states), 2)
        self.assertEqual(sum(len(s['conformations']) for s in states), len(doc['conformations']))
        r = analyse(doc['sequence'], states, 37.0, 0, 1)
        self.assertAlmostEqual(sum(s['relative'] for s in r['states']), 1.0, places=6)
        self.assertTrue(all(0 < s['p_ensemble'] <= 1 for s in r['states']))
        sc = r['scan']
        self.assertEqual(len(sc['matrix']), len(doc['sequence']))
        self.assertTrue(all(len(row) == 3 for row in sc['matrix']))
        # a mutation that breaks a pair only present in A must shift the balance towards B
        A, B = sc['dbn_a'], sc['dbn_b']
        only_a = [i for i in range(len(A)) if A[i] != '.' and B[i] == '.']
        if only_a:
            i = only_a[0]
            worst = max(sc['matrix'][i].values(), key=lambda v: v['ddg'])
            self.assertGreater(worst['ddg'], 0)
        self.assertTrue(all(x['ddg'] < 0 for x in sc['towards_a']))
        self.assertTrue(all(x['ddg'] > 0 for x in sc['towards_b']))
        self.assertEqual(compatible('GAAAC', '(...)'), '(...)')
        self.assertEqual(compatible('GAAAA', '(...)'), '.....')


class TestHtmlExport(unittest.TestCase):
    def test_page_is_self_contained(self):
        from rcv.htmlexport import MODULES, build_html
        doc = load_document([PRED])
        page = build_html(doc, {'conf': 1}, 'TAR </script> test', 'methods')
        self.assertTrue(page.startswith('<!doctype html>'))
        for m in MODULES:
            self.assertIn(f'data-path="{m}"', page)
        data = page.split('<script type="application/json" id="rcv-data">', 1)[1].split('</script>', 1)[0]
        parsed = json.loads(data)
        self.assertEqual(parsed['title'], 'TAR </script> test')
        self.assertEqual(len(parsed['doc']['conformations']), 10)
        self.assertNotIn('src="http', page)                         # nothing loaded from the network
        self.assertIn('MIT License', page)


class TestServerAdditions(unittest.TestCase):
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

    def req(self, method, path, body=None, raw=None, headers=None):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=120)
        h = {'Host': f'localhost:{self.port}', 'X-Viewer-Token': self.state.token, **(headers or {})}
        data = raw
        if body is not None:
            data = json.dumps(body).encode()
            h['Content-Type'] = 'application/json'
        c.request(method, path, body=data, headers=h)
        r = c.getresponse()
        out = r.read()
        hdr = dict(r.getheaders())
        c.close()
        return r.status, out, hdr

    def test_rcv_pack_open_and_fast_reopen(self):
        st, raw, _ = self.req('POST', '/api/open', {'paths': [str(PRED)]})
        doc_id = json.loads(raw)['id']
        session = {'format': 'rcv.view/1', 'view': {'conf': 3, 'colorScheme': 'helix'}}
        st, data, hdr = self.req('POST', '/api/rcv/pack', {'document': doc_id, 'session': session,
                                                           'thumbnail': 'data:image/png;base64,' + base64.b64encode(b'png').decode()})
        self.assertEqual(st, 200, data[:200])
        stats = json.loads(hdr['X-RCV-Stats'])
        self.assertEqual(stats['bytes'], len(data))
        self.assertGreater(stats['pdb_bytes'], 8 * len(data))
        st, raw, _ = self.req('POST', '/api/rcv/open', raw=data, headers={'Content-Type': 'application/octet-stream', 'X-Name': 'tar.rcv'})
        self.assertEqual(st, 200, raw)
        s = json.loads(raw)
        self.assertEqual(s['session'], session)
        self.assertEqual(s['n_conformations'], 10)
        # opened from disk too, and extensible like any ensemble
        f = Path(self.tmp) / 'tar.rcv'
        f.write_bytes(data)
        st, raw, _ = self.req('POST', '/api/open', {'paths': [str(f)]})
        self.assertEqual(st, 200, raw)
        entry = self.state.entries[json.loads(raw)['id']]
        self.assertEqual(len(entry.confs), 10)
        self.assertEqual(entry.rcv_path, str(f))
        st, raw, _ = self.req('POST', '/api/open', {'paths': [str(f), str(NMR)]})
        self.assertEqual(st, 400)
        st, raw, _ = self.req('GET', f'/api/documents/{doc_id}/methods')
        self.assertIn('rnafold3d engine', json.loads(raw)['methods'])

    def test_keep_alive_after_post_without_body_use(self):
        # regression: POST /api/open-demo sends '{}' that the handler does not need; the next request on the same
        # connection must still parse ('501 Unsupported method ({}GET)' before the fix)
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=120)
        h = {'Host': f'localhost:{self.port}', 'X-Viewer-Token': self.state.token, 'Content-Type': 'application/json'}
        for _ in range(2):
            c.request('POST', '/api/bye', body=b'{}', headers=h)
            r = c.getresponse(); r.read()
            self.assertEqual(r.status, 200)
            c.request('GET', '/api/info', headers=h)
            r = c.getresponse(); r.read()
            self.assertEqual(r.status, 200)
        c.request('POST', '/api/info', body=b'{"x":1}', headers={**h, 'Host': 'evil.example.com'})
        r = c.getresponse(); r.read()
        self.assertEqual(r.status, 403)
        c.request('GET', '/api/info', headers=h)
        r = c.getresponse(); r.read()
        self.assertEqual(r.status, 200)
        c.close()
        self.state.closing_at = None

    def test_html_switch_and_guards(self):
        st, raw, _ = self.req('POST', '/api/open', {'paths': [str(PRED)]})
        doc_id = json.loads(raw)['id']
        st, page, hdr = self.req('POST', '/api/export-html', {'document': doc_id, 'view': {'conf': 0}})
        self.assertEqual(st, 200)
        self.assertTrue(hdr['Content-Type'].startswith('text/html'))
        self.assertIn(b'x-rcv-module', page)
        if HAVE_VIENNA:
            st, raw, _ = self.req('POST', '/api/switch', {'document': doc_id, 'a': 0, 'b': 1})
            self.assertEqual(st, 200, raw)
            self.assertIn('scan', json.loads(raw))
        st, raw, _ = self.req('POST', '/api/read-text', {'path': str(NMR)})
        self.assertEqual(st, 400)
        script = Path(self.tmp) / 'fig.rcvs'
        script.write_text('select resi 1-5\n')
        st, raw, _ = self.req('POST', '/api/read-text', {'path': str(script)})
        self.assertEqual(json.loads(raw)['text'], 'select resi 1-5\n')
        st, raw, _ = self.req('POST', '/api/fetch', {'id': 'x'})               # validated before any download
        self.assertEqual(st, 400)
        st, raw, _ = self.req('POST', '/api/batch', {'fasta': '>a\nACGU\n'})
        self.assertEqual(st, 400)


class TestCommandLine(unittest.TestCase):
    def test_convert_info_extract_html(self):
        env = {'PYTHONPATH': str(ROOT / 'src')}
        import os
        env = {**os.environ, **env}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'nmr.rcv'
            r = subprocess.run([sys.executable, '-m', 'rcv', 'convert', str(NMR), '-o', str(out)], capture_output=True, text=True, env=env, timeout=300)
            self.assertEqual(r.returncode, 0, r.stderr)
            info = json.loads(r.stdout.strip().splitlines()[-1])
            self.assertGreater(info['ratio'], 8)
            r = subprocess.run([sys.executable, '-m', 'rcv', 'info', str(out)], capture_output=True, text=True, env=env, timeout=120)
            self.assertIn('20 conformations', r.stdout)
            r = subprocess.run([sys.executable, '-m', 'rcv', 'extract', str(out), '-o', str(Path(tmp) / 'pdb')], capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(len(list((Path(tmp) / 'pdb').glob('*.pdb'))), 20, r.stderr)
            r = subprocess.run([sys.executable, '-m', 'rcv', 'html', str(out), '-o', str(Path(tmp) / 'p.html')], capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertGreater((Path(tmp) / 'p.html').stat().st_size, 500_000)


@unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
class TestBrowserModules05(unittest.TestCase):
    def test_selection_layout_transition_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            doc = Path(tmp) / 'doc.json'
            doc.write_text(json.dumps(load_document([DEMO / '1EHZ_tRNA_prediction'])))
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'rcv05_unit.mjs'), str(doc)], capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)


if __name__ == '__main__':
    unittest.main()
