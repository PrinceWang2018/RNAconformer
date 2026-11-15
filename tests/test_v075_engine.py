"""0.7.5 engine and data: motif fragments, pseudoknot hypotheses, domains of long RNAs, the base-base potential, the
faster nested-subset routine, and SAXS / FRET / cryo-EM data with maximum-entropy reweighting (rcv/expdata.py)."""
import base64
import http.client
import json
import shutil
import struct
import sys
import tempfile
import threading
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import basepot as BP                  # noqa: E402
from rnafold3d import domains as DOM                 # noqa: E402
from rnafold3d import fragments as FR                # noqa: E402
from rnafold3d.motifs import nested_greedy           # noqa: E402
from rnafold3d.templates import load_library         # noqa: E402
from rcv import expdata as XD, server as srv         # noqa: E402
from rcv.rcvfile import doc_coords                   # noqa: E402

DEMO = ROOT / 'examples' / 'viewer_demo'
try:
    import RNA  # noqa: F401
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False


def _brute_nested(pairs):
    """Reference: remove the pair with most crossings, recomputing all counts each time."""
    p = sorted(pairs)
    while True:
        cnt = [sum((a < c < b < d) or (c < a < d < b) for c, d in p) for a, b in p]
        if not p or max(cnt) == 0:
            return p
        p.pop(int(np.argmax(cnt)))


class TestNestedAndFragments(unittest.TestCase):
    def test_nested_greedy_matches_reference(self):
        rng = np.random.default_rng(3)
        for _ in range(30):
            n = 60
            res = rng.permutation(n)[:30]
            pairs = sorted({(int(min(a, b)), int(max(a, b))) for a, b in zip(res[::2], res[1::2]) if abs(a - b) > 3})
            self.assertEqual(nested_greedy(pairs), _brute_nested(list(pairs)))

    def test_hairpin_loops(self):
        #           ((((....))))..((((...))))
        pairs = [(0, 11), (1, 10), (2, 9), (3, 8), (14, 24), (15, 23), (16, 22), (17, 21)]
        self.assertEqual(FR.hairpin_loops(pairs), [(3, 8), (17, 21)])

    @unittest.skipUnless(load_library() is not None, 'template library not installed')
    def test_gnra_hairpin_gets_a_gnra_fragment(self):
        lib = load_library()
        seq = 'GGCGGAAACGCC'                         # GC-closed GAAA tetraloop
        pairs = [(0, 11), (1, 10), (2, 9), (3, 8)]
        allowed = lib.allowed('2010-01-01')
        cands = FR.hairpin_candidates(seq, 3, 8, lib, allowed)
        self.assertTrue(cands)
        self.assertEqual(cands[0][0], 1.0)            # an identical loop and closing pair exist
        dr, tors, loops, used = FR.hairpin_fragments(seq, pairs, lib, allowed, np.random.default_rng(1))
        self.assertEqual(loops, [(3, 8)])
        self.assertTrue(all(r in tors for r in range(4, 8)))
        self.assertEqual(used[0]['source_loop'][0] + used[0]['source_loop'][2:], 'GAA')   # G N R A
        self.assertTrue(dr['i'] and len(dr['i']) == len(dr['d']) == len(dr['w']))
        m = FR.merge_restraints(dr, {'i': [0], 'j': [5], 'd': [9.0], 'tol': [1.0]})
        self.assertEqual(len(m['i']), len(dr['i']) + 1)
        self.assertEqual(m['w'][-1], 1.0)


class TestDomains(unittest.TestCase):
    def test_split_at_exterior_loop_and_multiloop(self):
        def hp(i, j, k=5):
            return [(i + t, j - t) for t in range(k)]
        # two 60-nt domains side by side, then one 200-nt helix closing a multiloop of two 70-nt branches
        pairs = hp(0, 59) + hp(62, 121) + hp(130, 329, 4) + hp(140, 209) + hp(215, 284)
        doms = DOM.split(340, pairs, max_len=160, min_len=25)
        # the two 60-nt domains fit together (122 <= 160); the 200-nt helix is split at its multiloop, whose two
        # branches again fit together into one domain
        self.assertEqual(doms, [(0, 122), (140, 285)])
        self.assertEqual(DOM.split(340, pairs, max_len=100, min_len=25), [(0, 60), (62, 122), (140, 210), (215, 285)])

    def test_torsions_of_an_ideal_stem(self):
        from rnafold3d.build import ResidueBuilder, load_geometry, build_stem
        rb = ResidueBuilder(load_geometry())
        seq = 'GGGGAAAACCCC'
        atoms, _ = build_stem(rb, seq, [(0, 11), (1, 10), (2, 9), (3, 8)], refine=False)
        res = [atoms[r] for r in (0, 1, 2, 3)]
        tor = DOM.torsions_of(res)
        af = rb.g.aform
        for name, col in (('alpha', 0), ('gamma', 2), ('epsilon', 4)):
            self.assertAlmostEqual(tor[1][col], af[name], delta=2.0)


class TestBasePotential(unittest.TestCase):
    @unittest.skipUnless(BP.table() is not None, 'potential table not built')
    def test_table_and_score(self):
        tab = BP.table()
        self.assertEqual(tab.shape, (BP.N_TYPES, 2, BP.N_D, BP.H_BINS))
        self.assertTrue(np.isfinite(tab).all() and np.abs(tab).max() <= 3.0)
        from rnafold3d.domains import read_model
        res = read_model(DEMO / 'HIV1_TAR_NMR_1ANR.pdb') if (DEMO / 'HIV1_TAR_NMR_1ANR.pdb').exists() else None
        if res:
            seq = ''.join('G' if 'O6' in r else 'A' if 'N6' in r else 'U' if 'O4' in r else 'C' for r in res)
            e = BP.score(seq, res, [])
            self.assertTrue(np.isfinite(e))
            self.assertNotEqual(e, 0.0)


@unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
class TestPseudoknots(unittest.TestCase):
    def test_h_type_pseudoknot_of_1rnk(self):
        from rnafold3d import pseudoknots as PK
        seq = 'GGCGCAGUGGGCUAGCGCCACUCAAAAGGCCCAU'      # the 1RNK H-type pseudoknot
        h = PK.hypothesis(seq)
        self.assertIsNotNone(h)
        self.assertIn('[', h['dbn'])
        self.assertGreaterEqual(len(h['restraints']), 4)
        self.assertLess(h['energy_kcal_mol'], h['pseudoknot']['mfe_kcal_mol'])
        for i, j in h['restraints']:
            self.assertTrue(any(a < i < b < j or i < a < j < b for a, b in h['pairs']))
        dr = PK.helix_restraints(seq, h['restraints'])
        self.assertEqual(len(dr['i']), len(set(r for p in h['restraints'] for r in p)) * (len(set(r for p in h['restraints'] for r in p)) - 1) // 2)

    def test_simple_hairpin_has_none(self):
        from rnafold3d import pseudoknots as PK
        self.assertIsNone(PK.hypothesis('GGGGAAACCCC'))

    def test_prediction_with_a_pseudoknot_hypothesis(self):
        from rnafold3d.predict import predict, PredictConfig
        out = Path(tempfile.mkdtemp())
        try:
            man = predict('GGCGCAGUGGGCUAGCGCCACUCAAAAGGCCCAU', out, config=PredictConfig(n_models=2, trajectories=3, workers=1,
                                                                                  templates='off', motifs='off', pseudoknots='auto'), log=None)
            hyps = man['input']['secondary_structures']
            self.assertEqual(hyps[0]['source'], 'vienna_pseudoknot')
            self.assertTrue(man['input']['secondary_structure_method']['pseudoknots']['found'])
            self.assertGreaterEqual(man['statistics']['valid'], 1)
        finally:
            shutil.rmtree(out, ignore_errors=True)


def _mrc(P, res=6.0, vox=1.5):
    lo = P.min(0) - 12
    n = np.ceil((P.max(0) + 12 - lo) / vox).astype(int)
    g = np.stack(np.meshgrid(*[np.arange(k) for k in n], indexing='ij'), -1).reshape(-1, 3) * vox + lo
    rho = XD._simulate(P, g, res).reshape(n)
    data = np.transpose(rho, (2, 1, 0)).astype('<f4')
    h = bytearray(1024)
    struct.pack_into('<4i', h, 0, n[0], n[1], n[2], 2)
    struct.pack_into('<3i', h, 28, *n)
    struct.pack_into('<3f', h, 40, *(n * vox))
    struct.pack_into('<3i', h, 64, 1, 2, 3)
    struct.pack_into('<3f', h, 196, *lo)
    h[208:212] = b'MAP '
    return bytes(h) + data.tobytes()


class TestExperimentalData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        st = srv.AppState(Path(cls.tmp) / 'runs', [], [])
        d = st.open_paths([DEMO / 'HIV1_TAR_prediction'])
        cls.doc = st.entries[d['id']].doc
        cls.C = doc_coords(cls.doc)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_saxs_profile_and_reweighting(self):
        el, X = XD._heavy(self.doc, self.C)
        q = np.linspace(0.01, 0.4, 60)
        I = XD.saxs_profile(el, X[0], q)
        slope = np.polyfit(q[:8] ** 2, np.log(I[:8]), 1)[0]
        self.assertAlmostEqual(float(np.sqrt(-3 * slope)), XD._rg(self.doc, self.C, 0), delta=1.0)   # Guinier Rg
        I3 = XD.saxs_profile(el, X[3], q)
        y = I3 * (1 + 0.02 * np.random.default_rng(0).normal(size=len(q)))
        text = '# q I sigma\n' + '\n'.join(f'{a} {b} {c}' for a, b, c in zip(q * 10, y, 0.02 * I3))   # in 1/nm
        data = XD.parse_saxs(text)
        self.assertIn('converted', data['unit'])
        r = XD.fit(self.doc, self.C, [{'kind': 'saxs', **data}])
        agr = r['datasets'][0]['agreement']
        self.assertEqual(int(np.argmin(agr)), 3)
        self.assertEqual(int(np.argmax(r['reweighting']['weights'])), 3)
        self.assertLess(r['reweighting']['chi2'][0], r['reweighting']['chi2_prior'][0])

    def test_fret_parse_and_fit(self):
        n = len(self.doc['sequence'])
        d = XD.parse_fret('donor acceptor value\n1 29 E=0.6 R0=54\n5 20 30 3  # comment\n', n)
        self.assertEqual(d['n_points'], 2)
        self.assertAlmostEqual(d['pairs'][0]['distance'], 54 * (1 / 0.6 - 1) ** (1 / 6), places=1)
        with self.assertRaises(ValueError):
            XD.parse_fret('1 1 30', n)
        r = XD.fit(self.doc, self.C, [{'kind': 'fret', **d}])
        self.assertEqual(len(r['datasets'][0]['agreement']), len(self.doc['conformations']))
        self.assertEqual(len(r['reweighting']['weights']), len(self.doc['conformations']))

    def test_map_parse_place_and_reweight(self):
        from scipy.spatial.transform import Rotation
        el, X = XD._heavy(self.doc, self.C)
        R = Rotation.from_euler('xyz', [40, -70, 15], degrees=True).as_matrix()
        raw = _mrc(X[2] @ R + np.array([30.0, -5.0, 12.0]))
        m = XD.parse_map(base64.b64encode(raw).decode(), 'synthetic', resolution=6.0)
        self.assertGreater(m['n_points'], 500)
        r = XD.fit(self.doc, self.C, [{'kind': 'map', **m}])
        cc = r['datasets'][0]['agreement']
        self.assertEqual(int(np.argmax(cc)), 2)
        self.assertGreater(cc[2], 0.95)
        self.assertGreater(r['reweighting']['weights'][2], 0.5)
        with self.assertRaises(ValueError):
            XD.parse_map(b'\0' * 100)

    def test_joint_maxent(self):
        rng = np.random.default_rng(5)
        O1, O2 = rng.normal(size=(6, 20)), rng.normal(size=(6, 4))
        w_true = np.array([0.7, 0.3, 0, 0, 0, 0])
        sets = [{'O': O1, 'y': w_true @ O1, 'sigma': np.full(20, 0.05), 'scale': 'none'},
                {'O': O2, 'y': 2.0 * (w_true @ O2) + 1.0, 'sigma': np.full(4, 0.05), 'scale': 'linear'}]
        r = XD.maxent(sets, 6, theta=0.5)
        self.assertGreater(r['weights'][0] + r['weights'][1], 0.9)
        self.assertTrue(all(a < b for a, b in zip(r['chi2'], r['chi2_prior'])))


class TestExpdataServer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        st, raw = cls._req('POST', '/api/open', {'paths': [str(DEMO / 'HIV1_TAR_prediction')]})
        cls.id = json.loads(raw)['id']

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

    def test_parse_and_fit_routes(self):
        st, raw = self._req('POST', '/api/expdata/parse', {'kind': 'fret', 'text': '1 29 40\n3 25 35 4', 'n': 29, 'name': 'f.txt'})
        self.assertEqual(st, 200, raw[:300])
        d = json.loads(raw)
        self.assertEqual((d['kind'], d['name'], d['n_points']), ('fret', 'f.txt', 2))
        st, raw = self._req('POST', '/api/expdata/fit', {'document': self.id, 'datasets': [d]})
        self.assertEqual(st, 200, raw[:300])
        r = json.loads(raw)
        self.assertIn('reweighting', r)
        st, raw = self._req('POST', '/api/expdata/parse', {'kind': 'xray', 'text': ''})
        self.assertEqual(st, 400)


if __name__ == '__main__':
    unittest.main()
