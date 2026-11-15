"""0.7.5: the compound screening window — rcv/screening.py, the /api/screen routes and the two-state screen
behind them (rnafold3d/screening/screen.py). Step 4 and part of step 5 of docs/NEXT_STEPS.md §1.7.

The run here is deliberately tiny (two small compounds, the shortest search there is) so that the test is
about the plumbing — library in, site chosen, job run, table out, pose back to the viewer — and not about
the quality of the scoring function, which is what `benchmark/screening/` measures.
"""
import http.client
import json
import shutil
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rcv import screening, server as srv            # noqa: E402
from rcv.data import load_document                  # noqa: E402
from rcv.rcvfile import doc_coords                  # noqa: E402

DEMO = ROOT / 'examples' / 'viewer_demo'
EXAMPLE = ROOT / 'src' / 'rcv' / 'static' / 'examples' / 'demo-compounds.sdf'


class TestLibrary(unittest.TestCase):
    def test_the_example_file_reads_as_four_compounds(self):
        s = screening.parse_library(EXAMPLE.read_text(encoding='utf-8'), name='demo')
        self.assertEqual(s['count'], 4)
        self.assertEqual(s['rejected_total'], 0)
        self.assertEqual(s['name'], 'demo')
        self.assertEqual([e['name'] for e in s['examples']][:2], ['adenine', 'theophylline'])
        self.assertTrue(10 <= s['atoms'][0] <= s['atoms'][1] <= s['atoms'][2] <= 20)
        self.assertTrue(100 < s['weight'][0] <= s['weight'][1] < 300)
        self.assertEqual(s['charge'][0], 0.0)                     # purines and xanthines: neutral at pH 7
        self.assertEqual(len(screening.library_molecules(s['handle'])), 4)

    def test_a_broken_record_is_counted_with_its_reason(self):
        text = EXAMPLE.read_text(encoding='utf-8') + 'not a molfile at all\n$$$$\n'
        s = screening.parse_library(text, name='mixed')
        self.assertEqual(s['count'], 4)
        self.assertEqual(s['rejected_total'], 1)
        self.assertEqual(sum(s['rejected'].values()), 1)

    def test_an_unknown_handle_says_so(self):
        with self.assertRaises(KeyError):
            screening.library_molecules('nosuchhandle')


class TestSite(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = load_document([DEMO / '1EHZ_tRNA_prediction'])
        cls.coords = np.asarray(doc_coords(cls.doc))

    def test_a_pocket_and_a_selection_both_give_a_box(self):
        site, info = screening.site_of(self.doc, self.coords, 0, {'mode': 'pocket', 'pocket': 1})
        self.assertEqual(info['mode'], 'pocket')
        self.assertEqual(info['pocket'], 1)
        self.assertGreater(len(info['residues']), 3)
        self.assertGreater(site.size, 6.0)
        chosen, info2 = screening.site_of(self.doc, self.coords, 0,
                                          {'mode': 'residues', 'residues': info['residues'][:5]})
        self.assertEqual(info2['mode'], 'residues')
        self.assertGreater(info2['n_atoms'], 10)
        self.assertGreater(chosen.size, 6.0)
        from rnafold3d.screening.screen import MAX_SIZE
        self.assertLessEqual(chosen.size, MAX_SIZE)
        with self.assertRaises(ValueError):
            screening.site_of(self.doc, self.coords, 0, {'mode': 'residues', 'residues': []})

    def test_the_other_state_is_measured_at_the_same_nucleotides(self):
        site, info = screening.site_of(self.doc, self.coords, 0, {'mode': 'pocket', 'pocket': 1})
        other, info_b = screening.matching_site(self.doc, self.coords, 1, site, info['residues'])
        self.assertEqual(info_b['mode'], 'residues')
        self.assertEqual(info_b['residues'], info['residues'])
        # the same nucleotides in another conformation: a box of a similar size, in a place of its own
        self.assertLess(abs(other.size - site.size), 0.5 * site.size + 8.0)
        # with no nucleotides to go by at all, the same box is used rather than a different pocket
        from rnafold3d.screening.screen import Site
        plain = Site(centre=site.centre, size=site.size)
        same, info_c = screening.matching_site(self.doc, self.coords, 1, plain, [])
        self.assertEqual(info_c['mode'], 'box')
        self.assertEqual(same.centre, site.centre)


class TestScreenRoutes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        st, raw = cls._req('POST', '/api/open', {'paths': [str(DEMO / 'HIV1_TAR_NMR_1ANR.pdb')]})
        assert st == 200, raw[:200]
        cls.doc = json.loads(raw)['id']

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def _req(cls, method, path, body=None, raw_body=None, headers=None):
        c = http.client.HTTPConnection('127.0.0.1', cls.port, timeout=600)
        h = {'Host': f'localhost:{cls.port}', 'X-Viewer-Token': cls.state.token,
             'Content-Type': 'application/json', **(headers or {})}
        payload = raw_body if raw_body is not None else (json.dumps(body).encode() if body is not None else None)
        c.request(method, path, body=payload, headers=h)
        r = c.getresponse()
        out = r.read()
        c.close()
        return r.status, out

    def post(self, path, body=None, raw_body=None, headers=None):
        st, raw = self._req('POST', path, body, raw_body, headers)
        self.assertEqual(st, 200, raw[:400])
        return json.loads(raw)

    def test_a_whole_run_through_the_routes(self):
        lib = self.post('/api/screen/library', raw_body=EXAMPLE.read_bytes(),
                        headers={'Content-Type': 'application/octet-stream', 'X-Name': 'demo-compounds.sdf'})
        self.assertEqual(lib['count'], 4)
        site = self.post('/api/screen/site', {'document': self.doc, 'conformation': 0,
                                              'site': {'mode': 'pocket', 'pocket': 1}})
        self.assertIn('size', site)
        self.assertEqual(len(site['centre']), 3)
        job = self.post('/api/screen', {'document': self.doc, 'a': 0, 'b': 1, 'library': lib['handle'],
                                        'site': {'mode': 'pocket', 'pocket': 1}, 'exhaustiveness': 'quick',
                                        'n_poses': 1, 'keep_poses': 2, 'jobs': 2, 'limit': 2})
        self.assertEqual(job['kind'], 'screen')
        for _ in range(600):
            st, raw = self._req('GET', f'/api/screen/{job["id"]}')
            self.assertEqual(st, 200, raw[:300])
            info = json.loads(raw)
            if info['status'] in ('done', 'failed', 'cancelled'):
                break
            time.sleep(0.5)
        self.assertEqual(info['status'], 'done', info.get('error'))
        self.assertEqual(info['result']['total'], 2)
        self.assertTrue(info['result']['two_state'])
        self.assertEqual([s['state'] for s in info['stages']], ['done'] * 5)

        table = json.loads(self._req('GET', f'/api/screen/{job["id"]}/results')[1])
        self.assertEqual(table['total'], 2)
        rows = table['rows']
        self.assertEqual(len(rows), 2)
        for r in rows:
            self.assertIsNotNone(r['score_a'])
            self.assertIsNotNone(r['score_b'])
            # each number is rounded on its own, so they agree to within that rounding
            self.assertAlmostEqual(r['selectivity'], r['score_b'] - r['score_a'], delta=0.011)
        self.assertGreaterEqual(rows[0]['selectivity'], rows[1]['selectivity'])     # ranked, best first
        self.assertEqual(sorted(table['sites']), ['a', 'b'])
        self.assertEqual(table['parameters']['exhaustiveness'], 'quick')

        st, csv = self._req('GET', f'/api/screen/{job["id"]}/csv')
        self.assertEqual(st, 200)
        lines = csv.decode().strip().split('\n')
        self.assertEqual(len(lines), 3)
        self.assertIn('selectivity', lines[0])
        st, sdf = self._req('GET', f'/api/screen/{job["id"]}/sdf?state=a')
        self.assertEqual(st, 200)
        self.assertEqual(sdf.decode().count('$$$$'), 2)
        self.assertIn('RCV_selectivity', sdf.decode())

        pose = self.post(f'/api/screen/{job["id"]}/pose', {'index': rows[0]['index'], 'state': 'a'})
        self.assertEqual(len(pose['xyz']), rows[0]['atoms'])
        self.assertEqual(len(pose['elements']), rows[0]['atoms'])
        self.assertTrue(pose['residues'])                     # the nucleotides the pose touches
        self.assertEqual(pose['conformation'], 0)
        self.assertTrue(all(0 <= i < len(pose['xyz']) and 0 <= j < len(pose['xyz']) for i, j in pose['bonds']))

    def test_bad_requests_are_refused(self):
        self.assertEqual(self._req('GET', '/api/screen/nosuchjob')[0], 404)
        lib = self.post('/api/screen/library', raw_body=EXAMPLE.read_bytes(),
                        headers={'Content-Type': 'application/octet-stream'})
        self.assertEqual(self._req('POST', '/api/screen', {'document': self.doc, 'library': lib['handle'],
                                                           'a': 99})[0], 400)
        self.assertEqual(self._req('POST', '/api/screen', {'document': 'missing'})[0], 400)
        self.assertEqual(self._req('POST', '/api/screen', {'document': self.doc, 'library': 'nope'})[0], 400)
        self.assertEqual(self._req('POST', '/api/screen/library', raw_body=b'')[0], 400)
        self.assertEqual(self._req('POST', '/api/screen/site',
                                   {'document': self.doc, 'conformation': 0,
                                    'site': {'mode': 'residues', 'residues': []}})[0], 400)


if __name__ == '__main__':
    unittest.main()
