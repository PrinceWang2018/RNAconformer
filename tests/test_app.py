"""Tests for the RCV application: multi-model reading, evaluator, ensemble documents, local API,
compute jobs, pull restraints, and (when Node.js is installed) the browser-side encoders, data model and pull logic."""
import base64
import http.client
import json
import re
import shutil
import subprocess
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

import numpy as np

from rnafold3d.build import atom_order
from rnafold3d.evaluate import evaluate
from rnafold3d.structure_io import read_models, write_pdb
from rcv import server as srv
from rcv.data import load_document, load_uploads
from rcv.desktop import app_window_command
from rcv.jobs import Job, PredictJob, pull_restraints

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'
NMR = DEMO / 'HIV1_TAR_NMR_1ANR.pdb'
PRED = DEMO / 'HIV1_TAR_prediction'


def _model_text(chain, shift=0.0):
    res = [{k: v + shift for k, v in r.atoms.items() if k in atom_order(r.base)} for r in chain.residues]
    with tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False) as t:
        name = t.name
    write_pdb(name, res, chain.sequence)
    text = Path(name).read_text()
    Path(name).unlink()
    return text


class TestMultiModel(unittest.TestCase):
    def test_nmr_file_gives_one_entry_per_model(self):
        models = read_models(NMR)
        self.assertEqual(len(models), 20)
        self.assertTrue(all(len(m[1]['A'].residues) == 29 for m in models))
        a, b = models[0][1]['A'].residues[5].atoms["C4'"], models[7][1]['A'].residues[5].atoms["C4'"]
        self.assertGreater(float(np.linalg.norm(a - b)), 0.01)     # models really differ

    def test_in_memory_text_and_single_model(self):
        models = read_models(NMR)
        text = _model_text(models[0][1]['A'])
        one = read_models(text=text, filename='x.pdb')
        self.assertEqual(len(one), 1)
        self.assertEqual(one[0][1]['A'].sequence, models[0][1]['A'].sequence)


class TestEvaluator(unittest.TestCase):
    def test_terms_add_up_and_pairs_found(self):
        ch = read_models(NMR)[0][1]['A']
        ev = evaluate(ch.sequence, [r.atoms for r in ch.residues])
        v = ev['values']
        self.assertGreaterEqual(v['n_pairs'], 8)                      # TAR stem + upper stem
        total = sum(v[k] or 0.0 for k in ('dg2d', 'elec', 'stack3d', 'steric'))
        self.assertAlmostEqual(v['eval_total'], round(total, 2), places=2)
        self.assertEqual(len(ev['dbn']), len(ch.sequence))
        self.assertGreater(v['elec'], 0.0)
        self.assertTrue(ev['validation']['passed'])


class TestViewerDocument(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nmr = load_document([NMR])
        cls.pred = load_document([PRED])

    def test_nmr_document(self):
        d = self.nmr
        self.assertEqual(d['schema'], 'rcv.document/1')
        self.assertEqual(len(d['conformations']), 20)
        self.assertEqual(d['default_metric'], 'eval_total')
        self.assertEqual(d['conformations'][0]['label'], 'Model 1')
        n_atoms = len(d['atoms']['names'])
        for c in d['conformations']:
            self.assertEqual(len(base64.b64decode(c['coords'])), n_atoms * 3 * 4)
        m = np.array(d['ensemble']['rmsd_c4'])
        self.assertTrue(np.allclose(m, m.T) and np.allclose(np.diag(m), 0.0))
        self.assertGreater(m.max(), 1.0)
        self.assertEqual(len(d['ensemble']['embedding']), 20)
        self.assertEqual(len(d['ensemble']['rmsf']), 29)
        self.assertGreaterEqual(len(d['bonds']) // 2, n_atoms - 1)
        self.assertGreaterEqual(len(d['helix_families']), 2)

    def test_ensemble_is_centred(self):
        d = self.nmr
        n = len(d['atoms']['names'])
        xyz = [np.frombuffer(base64.b64decode(c['coords']), '<f4').reshape(n, 3) for c in d['conformations']]
        c4 = [i for i, a in enumerate(d['atoms']['names']) if a == "C4'"]
        centre = np.mean([x[c4].mean(axis=0) for x in xyz], axis=0)
        self.assertLess(float(np.linalg.norm(centre)), 1.0)       # ensemble centred on the origin

    def test_engine_package_document(self):
        d = self.pred
        self.assertEqual(d['origin'], 'rnafold3d-engine')
        self.assertEqual(d['default_metric'], 'score')
        score = next(m for m in d['metrics'] if m['key'] == 'score')
        for c in d['conformations']:
            parts = sum(c['values'].get(t, 0.0) or 0.0 for t in score['terms'])
            self.assertAlmostEqual(parts, c['values']['score'], delta=0.01)
            self.assertTrue(c['validation']['passed'])

    def test_uploads_and_sequence_mismatch(self):
        models = read_models(NMR)
        files = [{'name': f'm{k}.pdb', 'text': _model_text(models[k][1]['A'])} for k in (0, 3)]
        d = load_uploads(files)
        self.assertEqual(len(d['conformations']), 2)
        bad = _model_text(models[1][1]['A']).replace(' G A   1', ' C A   1')
        with self.assertRaises(ValueError):
            load_uploads([files[0], {'name': 'bad.pdb', 'text': bad}])


class TestServer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.state = state
        cls.server = srv.AppServer(('127.0.0.1', 0), state)
        cls.port = cls.server.server_address[1]
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def req(self, method, path, body=None, token=True, host=None):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=60)
        headers = {'Host': host or f'localhost:{self.port}'}
        if token:
            headers['X-Viewer-Token'] = self.state.token
        data = None
        if body is not None:
            data = json.dumps(body).encode()
            headers['Content-Type'] = 'application/json'
        c.request(method, path, body=data, headers=headers)
        r = c.getresponse()
        raw = r.read()
        c.close()
        return r.status, raw

    def test_static_and_security(self):
        st, raw = self.req('GET', '/', token=False)
        self.assertEqual(st, 200)
        self.assertIn(b'RCV', raw)
        self.assertEqual(self.req('GET', '/api/info', token=False)[0], 401)
        self.assertEqual(self.req('GET', '/api/info', host='evil.example.com')[0], 403)
        self.assertEqual(self.req('GET', '/../server.py', token=False)[0], 404)
        self.assertEqual(self.req('GET', '/%2e%2e/server.py', token=False)[0], 404)
        self.assertEqual(self.req('GET', '/vendor/three/three.module.min.js', token=False)[0], 200)

    def test_open_documents_and_endpoints(self):
        st, raw = self.req('GET', '/api/info')
        info = json.loads(raw)
        self.assertEqual(st, 200)
        self.assertGreaterEqual(len(info['plugins']), 2)
        st, raw = self.req('POST', '/api/open', {'paths': [str(NMR)]})
        self.assertEqual(st, 200, raw)
        doc_id = json.loads(raw)['id']
        st, raw = self.req('GET', f'/api/documents/{doc_id}')
        self.assertEqual(json.loads(raw)['schema'], 'rcv.document/1')
        st, raw = self.req('GET', '/api/library')
        self.assertTrue(any(p['name'].startswith('HIV1_TAR') for p in json.loads(raw)['packages']))
        st, raw = self.req('GET', '/api/browse?path=' + str(DEMO))
        names = [e['name'] for e in json.loads(raw)['entries']]
        self.assertIn('HIV1_TAR_NMR_1ANR.pdb', names)
        self.assertEqual(self.req('POST', '/api/client-log', {'message': 'test'})[0], 200)
        st, raw = self.req('POST', '/api/open', {'paths': ['/nonexistent/x.pdb']})
        self.assertEqual(st, 400)
        st, raw = self.req('POST', '/api/predict', {'sequence': 'ACGUX'})
        self.assertEqual(st, 400)
        self.assertEqual(self.req('DELETE', f'/api/documents/{doc_id}')[0], 200)

    def test_app_window_command(self):
        cmd = app_window_command('http://localhost:1/#token=x')
        if cmd:
            self.assertTrue(any(a.startswith('--app=') for a in cmd))

    def test_desktop_endpoints_are_guarded(self):
        st, raw = self.req('GET', '/api/info')
        info = json.loads(raw)
        self.assertEqual(info['app'], 'RNA Conformation Viewer')
        self.assertGreaterEqual(info['system']['cpus'], 1)
        # files can only be written to paths approved by a native Save dialog
        st, _ = self.req('POST', '/api/write-file', {'x': 1})
        self.assertEqual(st, 400)
        st, _ = self.req('POST', '/api/read-json', {'path': str(NMR)})
        self.assertEqual(st, 400)
        st, _ = self.req('POST', '/api/explore', {'document': 'nope', 'conformation': 0, 'residue': 0, 'deformed': []})
        self.assertEqual(st, 400)
        st, raw = self.req('POST', '/api/predict', {'sequence': 'ACGU'})
        self.assertEqual(st, 400)
        self.assertIn('too short', json.loads(raw)['error'])

    def test_extend_document_with_exploration_package(self):
        s = self.state.open_paths([PRED])
        entry = self.state.entries[s['id']]
        n0 = len(entry.doc['conformations'])
        summary, added = self.state.extend_document(entry, PRED)        # any engine package works as "new" conformers
        self.assertEqual(len(added), n0)
        doc = self.state.entries[summary['id']].doc
        self.assertEqual(len(doc['conformations']), 2 * n0)
        self.assertTrue(all(c['explored'] for c in doc['conformations'][n0:]))
        self.assertEqual(summary['replaces'], s['id'])


class TestCli(unittest.TestCase):
    def test_export_document(self):
        from rcv.app import main
        from rnafold3d.cli import main as engine_main
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'doc.json'
            self.assertEqual(main(['export-document', str(PRED), '--out', str(out)]), 0)
            self.assertEqual(json.loads(out.read_text())['schema'], 'rcv.document/1')
            out2 = Path(tmp) / 'doc2.json'
            self.assertEqual(engine_main(['export-view', str(PRED), '--out', str(out2)]), 0)   # engine CLI forwards


class TestJobs(unittest.TestCase):
    def test_pull_restraints_tie_pulled_residue_to_unmoved_anchors(self):
        n = 30
        x0 = np.array([[3.5 * i, 0.0, 0.0] for i in range(n)])
        x1 = x0.copy()
        x1[14:17, 1] += np.array([6.0, 10.0, 6.0])               # residue 15 (0-based) pulled sideways
        rest, info = pull_restraints(x0, x1, 15)
        self.assertEqual(info['pulled'][0], 16)
        anchors = {r['j'] for r in rest}
        self.assertTrue(all(abs(a - 15) >= 4 for a in anchors))
        self.assertGreaterEqual(len(anchors), 3)
        for r in rest:
            self.assertAlmostEqual(r['target'], float(np.linalg.norm(x1[r['i']] - x1[r['j']])), places=1)

    def test_engine_log_parsing(self):
        state = type('S', (), {'projects_dir': tempfile.gettempdir()})()
        job = PredictJob(state, {'sequence': 'GGCACUUCGGUGCC', 'models': 2, 'sampling': 'fast', 'cores': 1})
        for line in ['[00:00:01] x: 14 nt; 1 secondary-structure hypotheses (ViennaRNA)',
                     '[00:00:01]   ss1 (((((....))))) -7.3 kcal/mol [vienna_mfe]',
                     '[00:00:01] running 6 trajectories on 1 worker(s)',
                     '[00:00:20]   progress 3/6 trajectories finished',
                     "[00:00:40]   traj   4 [ss1] rejected score3d=None ['severe_clash'] 12.0s",
                     '[00:00:41] 5/6 valid models; exported 2 conformers to /x in 40s',
                     '[00:00:41] WARNING: only 2 distinct valid conformers']:
            job.on_line(line)
        stages = {s['id']: s['state'] for s in job.stages}
        self.assertEqual(stages['ss'], 'done')
        self.assertEqual(stages['select'], 'done')
        self.assertEqual(job.result['valid'], 5)
        texts = ' '.join(n['text'] for n in job.notes)
        self.assertIn('1 secondary-structure hypothesis from ViennaRNA', texts)
        self.assertIn('rejected by validation: severe clash', texts)
        self.assertTrue(any(n['level'] == 'warning' for n in job.notes))
        self.assertGreater(job.progress, 0.9)


@unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
class TestBrowserModules(unittest.TestCase):
    def test_js_units_and_encoders(self):
        from PIL import Image
        with tempfile.TemporaryDirectory() as tmp:
            doc = Path(tmp) / 'doc.json'
            doc.write_text(json.dumps(load_document([PRED])))
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'app_unit.mjs'), str(doc), tmp],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            rgb = Image.open(Path(tmp) / 'rgb.tif')
            rgb.load()
            self.assertEqual((rgb.mode, rgb.size), ('RGB', (40, 30)))
            self.assertEqual(rgb.getpixel((10, 3)), (60, 24, 120))
            self.assertAlmostEqual(rgb.info['dpi'][0], 600, places=0)
            rgba = Image.open(Path(tmp) / 'rgba.tif')
            rgba.load()
            self.assertEqual(rgba.mode, 'RGBA')
            self.assertEqual(rgba.getpixel((1, 1))[3], 0)
            with Image.open(Path(tmp) / 'dpi.png') as png:
                self.assertAlmostEqual(png.info['dpi'][0], 600, places=0)
            z = zipfile.ZipFile(Path(tmp) / 'frames.zip')
            self.assertIsNone(z.testzip())
            self.assertEqual(z.read('notes.txt'), b'rna')
            pdf = (Path(tmp) / 'figure.pdf').read_bytes()
            xref = int(re.search(rb'startxref\n(\d+)', pdf).group(1))
            rows = pdf[xref:].split(b'\n')
            n = int(rows[1].split()[1])
            for i in range(1, n):
                off = int(rows[2 + i][:10])
                self.assertTrue(pdf[off:].startswith(f'{i} 0 obj'.encode()), f'xref entry {i}')


if __name__ == '__main__':
    unittest.main()
