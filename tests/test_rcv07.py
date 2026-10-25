"""Tests for RCV 0.7: solution conditions (salt, ligand), transition analysis (paths, network flux, bottleneck,
phi-value scan, rigid modules), the topology-aware ensemble map, probing data (parsing, agreement, reweighting),
the confidence model's features and labels, the new server endpoints, the HTML page, and (with Node.js) the
mechanism morph and the map drawing."""
import http.client
import json
import math
import shutil
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path

import numpy as np

from rcv import server as srv
from rcv import landscape, probing
from rcv.data import load_document
from rcv.rcvfile import doc_coords

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'
TAR = 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC'
THEO = 'GGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACC'

try:
    import RNA  # noqa: F401
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False


def _docs():
    if not hasattr(_docs, 'cache'):
        _docs.cache = {'tar': load_document([DEMO / 'HIV1_TAR_prediction']), 'trna': load_document([DEMO / '1EHZ_tRNA_prediction']),
                       'nmr': load_document([DEMO / 'HIV1_TAR_NMR_1ANR.pdb'])}
    return _docs.cache


@unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
class TestConditions(unittest.TestCase):
    def test_validation_and_salt(self):
        from rcv.conditions import Conditions
        with self.assertRaises(ValueError):
            Conditions.from_json({'temperature': 150})
        with self.assertRaises(ValueError):
            Conditions.from_json({'sodium': 10})
        with self.assertRaises(ValueError):
            Conditions.from_json({'ligand': {'motif_seq': 'GAU', 'motif_struct': '(.'}})
        s = '((((((...((((......))))))))))'
        e1 = Conditions.from_json({}).energy(Conditions.from_json({}).fold_compound(TAR), TAR, s)
        low = Conditions.from_json({'sodium': 0.1})
        e2 = low.energy(low.fold_compound(TAR), TAR, s)
        self.assertGreater(e2, e1 + 0.5)                       # less salt destabilises the helix
        self.assertIn('0.1 M Na+', low.describe())

    def test_ligand_motif_bonus_matches_fold(self):
        from rcv.conditions import Conditions
        c = Conditions.from_json({'ligand': {'preset': 'theophylline', 'conc_uM': 100}})
        self.assertAlmostEqual(c.binding_energy(), -c.kT * math.log1p(100 / 0.32), places=6)
        fc = c.fold_compound(THEO)
        s, e = fc.mfe()
        self.assertEqual(len(c.motif_sites(THEO)), 1)
        self.assertAlmostEqual(c.energy(fc, THEO, s), e, places=2)       # eval + motif bonus == constrained MFE
        self.assertEqual(c.motif_bonus('.' * len(THEO), c.motif_sites(THEO)), 0.0)
        self.assertEqual(Conditions.from_json({'ligand': {'preset': 'theophylline', 'conc_uM': 0}}).binding_energy(), 0.0)

    def test_dose_response(self):
        from rcv.conditions import Conditions
        from rcv.switch import dose_response, compatible
        bound = compatible(THEO, '((((...((((((((.....)))))...)))...))))')
        other = compatible(THEO, '((((((((....))))))))..................')
        d = dose_response(THEO, [{'dbn': bound}, {'dbn': other}], Conditions.from_json({'ligand': {'preset': 'theophylline'}}), [0, 1, 100])
        share = [r['relative'][0] for r in d['rows']]
        self.assertTrue(share[0] <= share[1] <= share[2])
        self.assertEqual(d['rows'][1]['bound'], [True, False])


@unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
class TestTransition(unittest.TestCase):
    def test_two_and_three_state_networks(self):
        from rcv.transition import transition_path_theory
        kT = 0.6
        r = transition_path_theory([0.0, 1.0], {(0, 1): 5.0}, 0, 1, kT, k0=1.0)
        self.assertAlmostEqual(r['rate'], math.exp(-5.0 / kT), delta=1e-9 * math.exp(-5 / kT))
        # A - I - B in series: the committor of I follows from its two rates
        r = transition_path_theory([0.0, -1.0, 0.5], {(0, 2): 4.0, (2, 1): 3.0}, 0, 1, kT, k0=1.0)
        kIA, kIB = math.exp(-(4.0 - 0.5) / kT), math.exp(-(3.0 - 0.5) / kT)
        self.assertAlmostEqual(r['committor'][2], kIB / (kIA + kIB), places=3)
        self.assertEqual(r['pathways'][0]['nodes'], [0, 2, 1])
        self.assertAlmostEqual(r['pathways'][0]['flux_share'], 1.0, places=3)
        # mean first-passage time with a deep minimum beside the direct route (solved by hand)
        r = transition_path_theory([-17.0, -18.4, -22.4], {(0, 1): -17.0, (1, 2): -18.4, (0, 2): -17.0}, 0, 1, kT, k0=1.0)
        kIB, kIA = math.exp(-4.0 / kT), math.exp(-5.4 / kT)
        tau = (1 + 1 / (kIA + kIB)) / (2 - kIA / (kIA + kIB))
        self.assertAlmostEqual(r['mfpt'] / tau, 1.0, places=6)

    def test_tar_pair(self):
        from rcv.conditions import Conditions
        from rcv import transition as T
        doc = _docs()['tar']
        r = T.analyse(doc, 0, 4, Conditions.from_json({}), coords=doc_coords(doc))
        self.assertEqual(r['A']['dbn'], '((((((...((((......))))))))))')
        self.assertAlmostEqual(r['direct']['barrier'], 4.1, places=1)
        self.assertLessEqual(r['best']['barrier'], r['direct']['barrier'] + 1e-9)
        self.assertEqual(r['bottleneck']['unit'], [[12, 19]])          # C13-G20 opens before the new pair forms
        self.assertGreater(r['rate'], 0)
        self.assertTrue(r['motion']['modules'])
        json.dumps(r)                                                     # plain JSON types only
        again = T.analyse(doc, 0, 4, Conditions.from_json({}), coords=doc_coords(doc))
        self.assertEqual(r['effective_barrier'], again['effective_barrier'])

    def test_network_and_trap(self):
        from rcv.conditions import Conditions
        from rcv import transition as T
        doc = _docs()['trna']
        r = T.analyse(doc, 0, 4, Conditions.from_json({}))
        self.assertGreater(len(r['nodes']), 2)
        q = r['tpt']['committor']
        self.assertTrue(all(0 <= x <= 1 for x in q))
        self.assertEqual((q[0], q[1]), (0.0, 1.0))
        shares = [p['flux_share'] for p in r['tpt']['pathways']]
        self.assertLessEqual(sum(shares), 1.0 + 1e-6)
        self.assertLessEqual(r['effective_barrier'], r['direct']['barrier'] + 1e-6)   # more routes, not slower
        self.assertEqual(r['nodes'][0]['conformations'][0], 0)                # basins map back to the conformations
        self.assertEqual(r['nodes'][1]['conformations'][0], 4)

    def test_basins(self):
        from rcv.conditions import Conditions
        from rcv import transition as T
        # NMR models of one fold: frayed variants that relax into the same 2D basin, so no barrier separates them
        nmr = _docs()['nmr']
        r = T.analyse(nmr, 0, 11, Conditions.from_json({}))
        self.assertTrue(r['same_basin'])
        self.assertEqual(r['A']['basin'], r['B']['basin'])
        self.assertGreater(r['B']['above_basin'], 0)
        self.assertIsNone(r['rate'])
        self.assertIsNone(r['tpt'])
        phases = [p['phase'] for p in r['best']['profile']]
        self.assertEqual(phases[0], 'relax_a')
        self.assertEqual(phases[-1], 'relax_b')
        self.assertEqual(r['best']['profile'][0]['s'], r['A']['dbn'])
        self.assertEqual(r['best']['profile'][-1]['s'], r['B']['dbn'])
        json.dumps(r)
        # a minimum is its own basin
        tar = _docs()['tar']
        r = T.analyse(tar, 0, 4, Conditions.from_json({}))
        self.assertEqual((r['A']['above_basin'], r['B']['above_basin']), (0.0, 0.0))
        self.assertFalse(r['same_basin'])

    def test_trap(self):
        from rcv.transition import find_trap, transition_path_theory
        kT = 0.6
        # A (0) and B (1) with a deep minimum (2) hanging off A: it has committor ~0 and lies 3 kcal/mol below A
        E = [-10.0, -11.0, -13.0]
        edges = {(0, 1): -4.0, (0, 2): -6.0}
        tpt = transition_path_theory(E, edges, 0, 1, kT, k0=1.0)
        trap = find_trap(E, tpt['committor'], tpt['occupancy'])
        self.assertEqual((trap['node'], trap['kind'], trap['depth']), (2, 'a-side', 3.0))
        self.assertGreater(trap['time_share'], 0.9)                                     # it holds the RNA
        self.assertEqual(trap, {**find_trap(E, tpt['committor']), 'time_share': trap['time_share'], 'nodes': [2]})
        self.assertAlmostEqual(sum(tpt['occupancy']), 1.0, places=3)
        # occupancy of the hand-solved network: the deep minimum beside the route holds the RNA (stochastic
        # simulation of the same network: 0.9986)
        r = transition_path_theory([-17.0, -18.4, -22.4], {(0, 1): -17.0, (1, 2): -18.4, (0, 2): -17.0}, 0, 1, kT, k0=1.0)
        self.assertAlmostEqual(r['occupancy'][2], 0.9986, places=3)
        # two near-degenerate wells split the time; together they are the trap
        E4 = [-10.0, -11.0, -13.0, -13.0]
        t4 = transition_path_theory(E4, {(0, 1): -4.0, (0, 2): -6.0, (2, 3): -11.0}, 0, 1, kT, k0=1.0)
        tr = find_trap(E4, t4['committor'], t4['occupancy'])
        self.assertEqual(sorted(tr['nodes']), [2, 3])
        self.assertGreater(tr['time_share'], 0.9)
        # the trap slows the switch: the mean first-passage time exceeds that of the bare A-B crossing
        bare = transition_path_theory(E[:2], {(0, 1): -4.0}, 0, 1, kT, k0=1.0)
        self.assertGreater(tpt['mfpt'], 1.5 * bare['mfpt'])
        # a shallow side minimum is no trap
        self.assertIsNone(find_trap([-10.0, -11.0, -10.2], [0.0, 1.0, 0.1]))
        # a deep minimum between the basins is
        self.assertEqual(find_trap([-10.0, -11.0, -12.0], [0.0, 1.0, 0.5])['kind'], 'deep')

    def test_phi_scan(self):
        from rcv.conditions import Conditions
        from rcv.transition import phi_scan
        sc = phi_scan(TAR, ['((((((...((((......))))))))))', '((((((...(((.((...)))))))))))'], Conditions.from_json({}))
        self.assertEqual(len(sc['matrix']), len(TAR))
        self.assertTrue(all(len(row) == 3 for row in sc['matrix']))
        for row in sc['matrix']:
            for v in row.values():
                if v['phi'] is not None:
                    self.assertGreaterEqual(abs(v['ddg_eq']), 0.6)
                    self.assertAlmostEqual(v['phi'], v['ddg_ts'] / v['ddg_eq'], delta=0.02)
        self.assertTrue(sc['faster'] and sc['slower'])

    def test_rigid_module_recovers_a_rotation(self):
        from rcv.transition import rigid_modules
        doc = _docs()['tar']
        X = doc_coords(doc)[0].astype(float)
        res = doc['residues']
        moving = [k for r in range(17, 29) for k in range(res[r]['first_atom'], res[r]['first_atom'] + res[r]['n_atoms'])]
        c = X[moving].mean(axis=0)
        ang = math.radians(40)
        R = np.array([[math.cos(ang), -math.sin(ang), 0], [math.sin(ang), math.cos(ang), 0], [0, 0, 1]])
        Y = X.copy()
        Y[moving] = (X[moving] - c) @ R.T + c + np.array([2.0, 0, 0])
        pairs = [p[:2] for p in doc['conformations'][0]['pairs'] if not (p[0] >= 17 and p[1] >= 17) and not (p[0] < 17 and p[1] >= 17)]
        mods = rigid_modules(doc, X, Y, pairs, pairs)
        rot = [m for m in mods['modules'] if set(range(20, 27)) <= set(m['residues'])]
        self.assertTrue(rot, mods['modules'])
        self.assertAlmostEqual(rot[0]['rotation'], 40.0, delta=2.0)
        self.assertAlmostEqual(abs(rot[0]['screw']['axis'][2]), 1.0, places=2)


class TestLandscape(unittest.TestCase):
    def test_shapes(self):
        self.assertEqual(landscape.shape_of('((((((...((((......))))))))))'), '[]')
        self.assertEqual(landscape.shape_of('((((...))))..((((...))))'), '[][]')
        self.assertEqual(landscape.shape_of('(((..((...))..((...))..)))'), '[[][]]')
        self.assertEqual(landscape.shape_of('.....'), '_')

    def test_layout_guarantees(self):
        rng = np.random.default_rng(1)
        centres = np.array([[0, 0, 0], [12, 0, 0], [6, 9, 0]], float)
        dbns = ['((((....))))........', '((((....))))((....))', '((..))((..))((..))..']          # shapes [], [][], [][][]
        P, S = [], []
        for k in range(45):
            P.append(centres[k % 3] + rng.normal(0, 3, 3))
            S.append(dbns[k % 3])
        D = np.linalg.norm(np.array(P)[:, None] - np.array(P)[None], axis=2)
        L = landscape.layout(D, S)
        X = np.array(L['positions'])
        self.assertEqual(len(L['classes']), 3)
        for c in L['classes']:
            d = np.linalg.norm(X[c['members']] - np.array(c['center']), axis=1)
            self.assertTrue((d <= c['radius'] + 1e-6).all())
        cs = L['classes']
        for i in range(len(cs)):
            for j in range(i + 1, len(cs)):
                gap = np.linalg.norm(np.array(cs[i]['center']) - np.array(cs[j]['center'])) - cs[i]['radius'] - cs[j]['radius']
                self.assertGreaterEqual(gap, -1e-6)
        q, m = L['quality']['map'], L['quality']['mds']
        self.assertGreaterEqual(q['topology_purity'], m['topology_purity'])       # neighbours share their topology
        self.assertGreaterEqual(q['trustworthiness'], m['trustworthiness'] - 0.02)  # without losing neighbourhoods
        self.assertEqual(L, landscape.layout(D, S))                     # deterministic

    def test_document_has_map(self):
        doc = _docs()['trna']
        L = doc['ensemble']['landscape']
        self.assertEqual(len(L['positions']), len(doc['conformations']))
        self.assertIn('trustworthiness', L['quality']['map'])


class TestProbing(unittest.TestCase):
    def test_parse_formats(self):
        self.assertEqual(probing.parse('1\t0.5\n2\t-999\n3\t1.2', 3)[0], [0.5, None, 1.2])
        self.assertEqual(probing.parse('1 0.5 0.01 G\n3 0.7 0.0 C', 3)[0], [0.5, None, 0.7])
        self.assertEqual(probing.parse('pos,nt,reactivity\n1,G,0.4\n2,C,nan', 2)[0], [0.4, None])
        self.assertEqual(probing.parse('0.1\n0.2', 2)[1], 'one value per line')
        with self.assertRaises(ValueError):
            probing.parse('a,b\n1,2', 2)

    def test_normalise_and_agreement(self):
        doc = _docs()['tar']
        U = probing.unpaired_matrix(doc)
        vals = [0.1 + 1.5 * u for u in U[0]]
        nv, f = probing.normalise(vals, doc['sequence'])
        self.assertGreater(f, 0)
        self.assertAlmostEqual(probing.agreement(doc, nv)[0], 1.0)
        dms, _ = probing.normalise(vals, doc['sequence'], 'DMS')
        self.assertTrue(all(v is None for v, b in zip(dms, doc['sequence']) if b in 'GU'))

    def test_reweighting_moves_towards_the_data(self):
        doc = _docs()['tar']
        U = probing.unpaired_matrix(doc)
        from rcv.switch import ensemble_states
        st = ensemble_states(doc)
        w = np.zeros(len(U))
        w[st[1]['conformations']] = 1 / len(st[1]['conformations'])
        y = 0.2 + 1.2 * (w @ U)
        r = probing.reweight(doc, list(y), theta=0.5)
        share = sum(r['weights'][k] for k in st[1]['conformations'])
        self.assertGreater(share, 0.6)                                    # prior share 0.2
        self.assertLess(r['chi2'], r['chi2_prior'])
        self.assertAlmostEqual(sum(r['weights']), 1.0, places=5)


class TestConfidenceFeatures(unittest.TestCase):
    def test_features_and_labels(self):
        from rnafold3d import qa
        from rnafold3d.metrics import read_numbered_pdb
        model = read_numbered_pdb(DEMO / 'HIV1_TAR_prediction' / 'conformations' / 'conf_001.pdb')
        keys = sorted(model)
        seq = ''.join(model[k].base for k in keys)
        F = qa.residue_features(seq, [model[k] for k in keys])
        self.assertEqual(F.shape, (len(keys), len(qa.feature_names())))
        self.assertFalse(np.isnan(F).any())
        self.assertGreater(F[:, 6].sum(), 10)                             # paired nucleotides found
        self.assertEqual(set(qa.lddt_per_residue(model, model).values()), {1.0})
        if qa.available():
            p = qa.predict(seq, [model[k] for k in keys])
            self.assertTrue(0 < p['global'] < 1)
            self.assertEqual(len(p['per_residue']), len(keys))


class TestServer07(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        st, raw = cls._req('POST', '/api/open', {'paths': [str(DEMO / '1EHZ_tRNA_prediction')]})
        cls.doc_id = json.loads(raw)['id']

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

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_switch_transition_scan(self):
        r = self.post('/api/switch', {'document': self.doc_id, 'conditions': {'temperature': 25, 'sodium': 0.15}, 'a': 0, 'b': 1})
        self.assertIn('25 °C, 0.15 M Na+', r['conditions_text'])
        t = self.post('/api/transition', {'document': self.doc_id, 'a': 1, 'b': 0, 'conditions': {}})
        self.assertIn('bottleneck', t)
        self.assertIn('motion', t)
        route = [t['nodes'][k]['dbn'] for k in t['best']['route']]
        sc = self.post('/api/transition/scan', {'document': self.doc_id, 'route': route})
        self.assertEqual(len(sc['matrix']), 76)
        st, raw = self._req('POST', '/api/transition', {'document': self.doc_id, 'a': 0, 'b': 99})
        self.assertEqual(st, 400)

    def test_landscape_probing_confidence_html(self):
        L = self.post('/api/landscape', {'document': self.doc_id})
        self.assertEqual(len(L['positions']), 5)
        p = self.post('/api/probing/parse', {'text': '\n'.join(f'{i + 1} {0.1 * (i % 7)}' for i in range(76)), 'sequence': 'A' * 76})
        self.assertEqual(p['n_values'], 76)
        f = self.post('/api/probing/fit', {'document': self.doc_id, 'values': p['normalised']})
        self.assertEqual(len(f['agreement']), 5)
        self.assertEqual(len(f['reweighting']['weights']), 5)
        c = self.post('/api/confidence', {'document': self.doc_id})
        if c['available']:
            self.assertEqual(len(c['conformations']), 5)
        st, raw = self._req('POST', '/api/export-html', {'document': self.doc_id, 'analysis': {'pathway': None}})
        self.assertEqual(st, 200)
        self.assertIn(b'js/landscape.js', raw)
        self.assertIn(b'js/mechmorph.js', raw)
        self.assertIn(b'"landscape"', raw)


@unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
class TestBrowserModules07(unittest.TestCase):
    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_mechanism_morph_and_map(self):
        from rcv.conditions import Conditions
        from rcv import transition as T
        doc = _docs()['trna']
        r = T.analyse(doc, 1, 0, Conditions.from_json({}), coords=doc_coords(doc))
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / 'case.json').write_text(json.dumps({'doc': doc, 'result': r}))
            out = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'rcv07_unit.mjs'), str(Path(tmp) / 'case.json')], capture_output=True, text=True, timeout=300)
            self.assertEqual(out.returncode, 0, out.stdout + out.stderr)


if __name__ == '__main__':
    unittest.main()
