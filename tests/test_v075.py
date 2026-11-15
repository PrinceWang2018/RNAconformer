"""0.7.5: pocket finding (rnafold3d/screening/pockets.py), the pocket list in the Details tab and
/api/insight/pockets. The first step of the conformation-selective compound screen (docs/NEXT_STEPS.md §1)."""
import http.client
import json
import shutil
import sys
import tempfile
import threading
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d.screening import pockets as P        # noqa: E402
from rcv import insight, server as srv              # noqa: E402
from rcv.data import load_document                  # noqa: E402
from rcv.rcvfile import doc_coords                  # noqa: E402

DEMO = ROOT / 'examples' / 'viewer_demo'
SAMPLES = ROOT / 'benchmark' / 'data' / 'samples'


def shell(n=900, radius=10.0):
    """`n` carbon atoms spread over a sphere: a closed cavity of known volume inside."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return radius * np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], 1)


def _doc(name='1EHZ_tRNA_prediction'):
    if name not in _doc.cache:
        _doc.cache[name] = load_document([DEMO / name])
    return _doc.cache[name]


_doc.cache = {}


class TestPocketGeometry(unittest.TestCase):
    def test_closed_cavity_has_the_volume_of_the_cavity(self):
        xyz = shell(900, 10.0)
        found = P.find(xyz, elements=['C'] * len(xyz), min_volume=50)
        self.assertEqual(len(found), 1)
        p = found[0]
        want = 4 / 3 * np.pi * (10.0 - P.VDW['C'] - P.PROBE) ** 3      # the space a probe centre can reach
        self.assertAlmostEqual(p['volume'], want, delta=0.1 * want)
        self.assertGreater(p['burial'], 0.9)                           # enclosed in every direction
        self.assertGreater(p['enclosure'], 0.5)    # not 1: the far wall of a 20 Å sphere is past the 12 Å scan
        self.assertEqual(p['id'], 1)
        self.assertEqual(p['score'], 1.0 if p['burial'] == 1.0 else p['burial'])   # no pairing given: geometry alone
        self.assertFalse(p['wide'])
        self.assertEqual(p['points'].shape, (p['n_points'], 3))
        self.assertLess(np.linalg.norm(p['center']), 1.0)              # the cavity is centred on the origin

    def test_an_open_surface_has_no_pocket(self):
        g = np.arange(-12.0, 12.1, 2.0)
        xyz = np.array([[x, y, 0.0] for x in g for y in g])            # a flat sheet of atoms
        self.assertEqual(P.find(xyz, elements=['C'] * len(xyz)), [])

    def test_the_result_does_not_depend_on_the_frame_or_the_call(self):
        xyz = shell(600, 11.0)
        rng = np.random.default_rng(7)
        q = rng.normal(size=(3, 3))
        R = np.linalg.qr(q)[0] * np.sign(np.linalg.det(np.linalg.qr(q)[0]))
        a = P.find(xyz, elements=['C'] * len(xyz), min_volume=50)
        b = P.find(xyz, elements=['C'] * len(xyz), min_volume=50)
        c = P.find(xyz @ R + 31.0, elements=['C'] * len(xyz), min_volume=50)
        self.assertEqual(a[0]['n_points'], b[0]['n_points'])
        self.assertEqual(a[0]['residues'], b[0]['residues'])
        self.assertAlmostEqual(a[0]['volume'], c[0]['volume'], delta=0.03 * a[0]['volume'])

    def test_burial_threshold_decides_what_counts_as_a_site(self):
        doc = _doc()
        xyz = np.asarray(doc_coords(doc))[0]
        kw = dict(elements=doc['atoms']['elements'], residue_of=doc['atoms']['residue'])
        shallow = P.find(xyz, **kw)                     # the default: 5 of 7, chosen on known ligand sites
        deep = P.find(xyz, burial_min=6, **kw)
        closed = P.find(xyz, burial_min=7, **kw)
        self.assertEqual(P.BURIAL_MIN, 5)
        self.assertGreaterEqual(sum(p['volume'] for p in shallow), sum(p['volume'] for p in deep))
        self.assertGreaterEqual(sum(p['volume'] for p in deep), sum(p['volume'] for p in closed))
        for p in shallow:
            self.assertGreaterEqual(p['volume'], P.MIN_VOLUME)
            self.assertTrue(set(p['residues']) <= set(range(len(doc['sequence']))))
            self.assertTrue(p['residues'])

    def test_lining_residues_are_the_nucleotides_at_the_pocket(self):
        doc = _doc()
        xyz = np.asarray(doc_coords(doc))[0]
        res_of = np.asarray(doc['atoms']['residue'])
        found = P.find(xyz, elements=doc['atoms']['elements'], residue_of=res_of)
        self.assertTrue(found)
        p = found[0]
        for r in p['residues']:
            d = np.linalg.norm(p['points'][:, None, :] - xyz[res_of == r][None], axis=2)
            self.assertLessEqual(d.min(), P.CONTACT + 1e-9)
        self.assertEqual(sorted(int(r) for r in p['contacts']), p['residues'])
        self.assertTrue(all(c >= 1 for c in p['contacts'].values()))
        self.assertLessEqual(sum(p['contacts'].values()), int(np.isin(res_of, p['residues']).sum()))

    def test_score_ranks_enclosed_loop_sites_above_helix_grooves(self):
        groove = {'burial': 0.9, 'loop_fraction': 0.0, 'volume': 400.0}
        loop = {'burial': 0.9, 'loop_fraction': 1.0, 'volume': 400.0}
        small = {'burial': 0.9, 'loop_fraction': 1.0, 'volume': 150.0}
        blind = {'burial': 0.9, 'loop_fraction': None, 'volume': 400.0}
        self.assertAlmostEqual(P._score(groove), 0.9 * 0.4, places=3)
        self.assertAlmostEqual(P._score(loop), 0.9, places=3)
        self.assertAlmostEqual(P._score(small), 0.9 * 150 / 300, places=3)
        self.assertAlmostEqual(P._score(blind), 0.9, places=3)
        self.assertLess(P._score(groove), P._score(loop))

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            P.find(np.zeros((3, 3)))
        with self.assertRaises(ValueError):
            P.find(np.zeros(9))
        with self.assertRaises(ValueError):                   # volume_near needs the grid points
            P.volume_near([{'volume': 1.0, 'residues': [0], 'points': None}], [0], np.zeros((2, 3)), [0, 0])
        self.assertEqual(P.volume_near([], [0], np.zeros((2, 3)), [0, 0]), (0.0, []))
        self.assertEqual(P.radii(['c', 'XX', 'P']).tolist(), [P.VDW['C'], P.DEFAULT_VDW, P.VDW['P']])


class TestPocketComparison(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = _doc()
        cls.xyz = np.asarray(doc_coords(cls.doc))
        cls.res_of = cls.doc['atoms']['residue']
        kw = dict(elements=cls.doc['atoms']['elements'], residue_of=cls.res_of)
        cls.a = P.find(cls.xyz[0], flexible=insight._flexible(cls.doc, 0), **kw)
        cls.b = P.find(cls.xyz[1], flexible=insight._flexible(cls.doc, 1), **kw)

    def test_volume_near_measures_the_same_site_in_both_conformations(self):
        site = self.a[0]['residues']
        v_own, parts = P.volume_near(self.a, site, self.xyz[0], self.res_of)
        self.assertGreater(v_own, 0.5 * self.a[0]['volume'])            # most of the pocket is at its own lining
        self.assertLessEqual(v_own, sum(p['volume'] for p in self.a) + 1e-6)
        self.assertEqual(parts[0][0], 0)
        self.assertEqual(P.volume_near(self.a, [], self.xyz[0], self.res_of), (0.0, []))

    def test_a_site_open_in_one_conformation_only(self):
        rows = P.compare(self.a, self.b, self.xyz[0], self.xyz[1], self.res_of)
        self.assertTrue(rows)
        for r in rows:
            self.assertIn(r['side'], ('A', 'B'))
            self.assertIn(r['open_in'], ('A', 'B', 'both'))
            self.assertGreater(max(r['volume_a'], r['volume_b']), 0)
            self.assertAlmostEqual(r['difference'], (r['volume_a'] - r['volume_b']) / max(r['volume_a'], r['volume_b']), places=3)
            self.assertEqual(r['open_in'], 'both' if abs(r['difference']) <= 0.25 else ('A' if r['difference'] > 0 else 'B'))
            self.assertTrue(r['residues'])
        self.assertEqual(sorted(abs(r['difference']) for r in rows)[::-1],
                         [abs(r['difference']) for r in rows])          # largest change first

    def test_the_comparison_is_symmetric(self):
        ab = P.compare(self.a, self.b, self.xyz[0], self.xyz[1], self.res_of)
        ba = P.compare(self.b, self.a, self.xyz[1], self.xyz[0], self.res_of)
        back = {tuple(r['residues']): r for r in ba}
        seen = 0
        for r in ab:
            m = back.get(tuple(r['residues']))
            if m is None:
                continue
            seen += 1
            self.assertAlmostEqual(r['difference'], -m['difference'], places=3)
            self.assertAlmostEqual(r['volume_a'], m['volume_b'], places=1)
            self.assertEqual(r['open_in'], {'A': 'B', 'B': 'A', 'both': 'both'}[m['open_in']])
        self.assertGreaterEqual(seen, len(self.a))                      # every site of A is measured both ways


class TestPocketsOfADocument(unittest.TestCase):
    def test_pockets_of_a_conformation(self):
        doc = _doc()
        r = insight.pockets(doc, 0, doc_coords(doc))
        json.dumps(r)                                                   # the answer goes to the browser as JSON
        self.assertEqual((r['conformation'], r['other']), (0, None))
        self.assertEqual(r['params']['burial_min'], P.BURIAL_MIN)
        self.assertEqual(r['n_pockets'], len(r['pockets']))
        n = len(doc['sequence'])
        for p in r['pockets']:
            self.assertGreaterEqual(p['volume'], P.MIN_VOLUME)
            self.assertTrue(0.0 <= p['burial'] <= 1.0 and 0.0 <= p['loop_fraction'] <= 1.0 and 0.0 <= p['score'] <= 1.0)
            self.assertNotIn('points', p)                               # the grid points stay on the server
            self.assertEqual(p['residues'], sorted(set(p['residues'])))
            self.assertEqual([x for s in p['segments'] for x in range(s['start'], s['end'] + 1)], p['residues'])
            self.assertTrue(all(0 <= x < n for x in p['residues']))
        self.assertEqual([p['score'] for p in r['pockets']], sorted((p['score'] for p in r['pockets']), reverse=True))

    def test_loop_fraction_follows_the_pairing_of_that_conformation(self):
        doc = _doc()
        f = insight._flexible(doc, 0)
        self.assertEqual(len(f), len(doc['sequence']))
        self.assertTrue(set(np.unique(f)) <= {0.0, 1.0})
        dbn = insight.conf_dbn(doc, 0)
        self.assertEqual([i for i, v in enumerate(f) if v == 0.0], [i for i, c in enumerate(dbn) if c != '.'])
        r = insight.pockets(doc, 0, doc_coords(doc))
        for p in r['pockets']:
            open_res = [x for x in p['residues'] if dbn[x] == '.']
            self.assertEqual(p['loop_fraction'] > 0, bool(open_res))

    def test_comparison_of_two_conformations(self):
        doc = _doc()
        r = insight.pockets(doc, 0, doc_coords(doc), other=2)
        json.dumps(r)
        self.assertEqual((r['conformation'], r['other']), (0, 2))
        self.assertIn('n_pockets_other', r)
        for p in r['pockets']:
            self.assertIn(p['open_in'], ('A', 'B', 'both'))
            self.assertGreaterEqual(p['volume_site'], 0.0)
            self.assertAlmostEqual(p['difference'], (p['volume_site'] - p['volume_other'])
                                   / max(p['volume_site'], p['volume_other']), places=2)
        for c in r['changes']:
            self.assertTrue(c['pocket'] is not None or c['other_pocket'] is not None)
            self.assertTrue(c['segments'])
        self.assertEqual(insight.pockets(doc, 1, doc_coords(doc), other=1)['other'], None)   # itself is no comparison
        with self.assertRaises(ValueError):
            insight.pockets(doc, 99, doc_coords(doc))
        with self.assertRaises(ValueError):
            insight.pockets(doc, 0, doc_coords(doc), other=99)

    @unittest.skipUnless((SAMPLES / '2GDI_A' / 'reference.pdb').exists(), 'sample structures not in the tree')
    def test_a_ligand_bound_crystal_structure_has_a_ligand_sized_site(self):
        # 2GDI is the TPP riboswitch crystallised with its ligand; the normalised reference holds the RNA only,
        # so the site the ligand occupied is empty space in the file. What is checked here is that the best
        # pocket of such a structure is the size of a real binding site and lined by a junction rather than by
        # a helix groove — not which nucleotides bind thiamine pyrophosphate, which this file cannot show.
        doc = load_document([SAMPLES / '2GDI_A' / 'reference.pdb'])
        r = insight.pockets(doc, 0, doc_coords(doc))
        self.assertTrue(r['pockets'])
        best = r['pockets'][0]
        self.assertGreater(best['volume'], 200.0)
        self.assertLess(best['volume'], P.WIDE_VOLUME)
        self.assertGreater(best['score'], 0.4)
        self.assertGreater(best['loop_fraction'], 0.3)                  # a junction site, not a helix groove
        self.assertGreater(len(best['residues']), 10)
        self.assertGreater(len(r['pockets'][0]['segments']), 1)         # it is closed by more than one strand


class TestPocketsRoute(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        st, raw = cls._req('POST', '/api/open', {'paths': [str(DEMO / '1EHZ_tRNA_prediction')]})
        assert st == 200, raw[:200]
        cls.doc = json.loads(raw)['id']

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    @classmethod
    def _req(cls, method, path, body=None):
        c = http.client.HTTPConnection('127.0.0.1', cls.port, timeout=300)
        h = {'Host': f'localhost:{cls.port}', 'X-Viewer-Token': cls.state.token, 'Content-Type': 'application/json'}
        c.request(method, path, body=json.dumps(body).encode() if body is not None else None, headers=h)
        r = c.getresponse()
        out = r.read()
        c.close()
        return r.status, out

    def post(self, path, body):
        st, raw = self._req('POST', path, body)
        self.assertEqual(st, 200, raw[:300])
        return json.loads(raw)

    def test_route(self):
        r = self.post('/api/insight/pockets', {'document': self.doc, 'conformation': 1})
        self.assertEqual((r['conformation'], r['other']), (1, None))
        self.assertIn('pockets', r)
        again = self.post('/api/insight/pockets', {'document': self.doc, 'conformation': 1})
        self.assertEqual(again, r)                                      # cached, and the analysis is deterministic
        c = self.post('/api/insight/pockets', {'document': self.doc, 'conformation': 1, 'other': 0})
        self.assertEqual(c['other'], 0)
        self.assertIn('changes', c)
        s = self.post('/api/insight/pockets', {'document': self.doc, 'conformation': 1, 'burial_min': 6})
        self.assertEqual(s['params']['burial_min'], 6)
        self.assertLessEqual(sum(p['volume'] for p in s['pockets']),
                             sum(p['volume'] for p in r['pockets']) + 1e-6)
        self.assertEqual(self._req('POST', '/api/insight/pockets', {'document': self.doc, 'conformation': 99})[0], 400)
        self.assertEqual(self._req('POST', '/api/insight/pockets', {'document': 'missing'})[0], 400)


if __name__ == '__main__':
    unittest.main()
