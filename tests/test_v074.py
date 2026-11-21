"""0.7.4: structure-aware template search, probing/alignment-guided 2D, per-template model order, and the advanced
analyses inside the Switch, Pathway and Details tabs (rcv/insight.py, /api/insight/*)."""
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

from rnafold3d import templates as T                 # noqa: E402
from rnafold3d import ss as SS                        # noqa: E402
from rnafold3d.predict import select_diverse          # noqa: E402
from rcv import insight, server as srv                # noqa: E402
from rcv.conditions import Conditions                 # noqa: E402
from rcv.data import load_document                    # noqa: E402
from rcv.rcvfile import doc_coords                    # noqa: E402
from rcv.switch import ensemble_states, pairs_of      # noqa: E402

DEMO = ROOT / 'examples' / 'viewer_demo'
SAMPLES = ROOT / 'benchmark' / 'data' / 'samples'
try:
    import RNA  # noqa: F401
    HAVE_VIENNA = True
except ImportError:
    HAVE_VIENNA = False
TRNA = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'      # yeast tRNA-Phe (1EHZ)
TAR = 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC'


def _docs():
    if not hasattr(_docs, 'cache'):
        _docs.cache = {k: load_document([DEMO / v]) for k, v in
                       (('tar', 'HIV1_TAR_prediction'), ('trna', '1EHZ_tRNA_prediction'), ('nmr', 'HIV1_TAR_NMR_1ANR.pdb'))}
    return _docs.cache


class TestStructureAlignment(unittest.TestCase):
    def test_score_matrix_core_is_the_sequence_alignment(self):
        q, t = 'GGCAGAUCUGAGCC', 'AAGGCAGAUGUGAGCCUU'
        qc = T.CODE[np.frombuffer(q.encode(), np.uint8)]
        tc = T.CODE[np.frombuffer(t.encode(), np.uint8)]
        self.assertEqual(T.align(q, t), T._sw_matrix(T._score_rows(qc)[:, tc]))

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_profile_and_states(self):
        p = T.pair_profile(TAR)
        self.assertEqual(p.shape, (len(TAR), 3))
        self.assertTrue(np.allclose(p.sum(1), 1.0, atol=1e-6))
        self.assertGreater(p[1, 0], 0.8)            # 5' strand of the stem pairs downstream
        self.assertGreater(p[-2, 1], 0.8)           # 3' strand pairs upstream

    @unittest.skipUnless(HAVE_VIENNA and T.LIBRARY_FILE.exists(), 'needs ViennaRNA and the structure library')
    def test_remote_trna_homologs_are_found_and_extend_the_alignment(self):
        lib = T.load_library()
        kw = dict(exclude=('1EHZ',), max_hits=8)                   # the structure of this very tRNA is left out
        seq_hits = T.search(TRNA, lib, **kw)
        all_hits = T.search(TRNA, lib, structure=True, **kw)
        self.assertGreaterEqual(len(all_hits), len(seq_hits))
        by = {h['found_by'] for h in all_hits}
        self.assertTrue(by & {'sequence + structure', 'sequence, aligned with structure'}, by)
        for h in all_hits:
            if h['found_by'] != 'sequence':
                self.assertGreaterEqual(h['coverage'], T.STRUCT_MIN_COVERAGE)
                self.assertGreaterEqual(h['pairs_compatible'], T.STRUCT_MIN_COMPAT * h['pairs_mapped'])
            q = [a for a, _ in h['aligned']]
            t = [b for _, b in h['aligned']]
            self.assertEqual(q, sorted(set(q)))
            self.assertEqual(t, sorted(set(t)))
        self.assertGreaterEqual(max(h['coverage'] for h in all_hits), max(h['coverage'] for h in seq_hits) if seq_hits else 0)

    @unittest.skipUnless(HAVE_VIENNA and T.LIBRARY_FILE.exists(), 'needs ViennaRNA and the structure library')
    def test_unrelated_sequence_gets_no_structure_hit(self):
        rng = np.random.default_rng(11)
        lib = T.load_library()
        for _ in range(3):
            q = ''.join(rng.choice(list('ACGU'), 80))
            self.assertEqual([h for h in T.search(q, lib, structure=True) if h['found_by'] != 'sequence'], [])


class TestGuided2D(unittest.TestCase):
    def test_read_reactivities(self):
        d = Path(tempfile.mkdtemp())
        try:
            (d / 'a.txt').write_text('# SHAPE\n' + '\n'.join(['0.1', 'NA', '-999', '1.5'] * 3) + '\n')
            v = SS.read_reactivities(d / 'a.txt', 12)
            self.assertEqual(v[:4], [0.1, None, None, 1.5])
            (d / 'b.txt').write_text('\n'.join(f'{i}\t{0.2 * i}' for i in range(1, 11)) + '\n')
            v = SS.read_reactivities(d / 'b.txt', 12)
            self.assertEqual(v[0], 0.2)
            self.assertIsNone(v[11])
            with self.assertRaises(ValueError):
                SS.read_reactivities(d / 'a.txt', 30)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_read_alignment(self):
        d = Path(tempfile.mkdtemp())
        try:
            (d / 'x.fa').write_text('>h1\nGGCAGAU-CUGAGCCUGGGAGCUCUCUGCC\n>query\nGGCAGAUCTGAGCCUGGGAGCUC\nUCUGCC\n'.replace('CTG', 'CUG', 1))
            (d / 'x.fa').write_text('>h1\nGGCAGAUC-GAGCCUGGGAGCUCUCUGCC\n>query\n' + TAR + '\n')
            aln = SS.read_alignment(d / 'x.fa', TAR)
            self.assertEqual(aln[0], TAR)
            (d / 'x.sto').write_text('# STOCKHOLM 1.0\nh1  ' + TAR.replace('G', 'A', 1) + '\nq   ' + TAR.replace('U', 'T') + '\n//\n')
            self.assertEqual(SS.read_alignment(d / 'x.sto', TAR)[0], TAR)
            with self.assertRaises(ValueError):
                SS.read_alignment(d / 'x.fa', 'ACGU' * 5)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_probing_changes_the_hypotheses_and_keeps_true_energies(self):
        seq = 'GGGAAACUUCGGUUUCCCAAAAGGGAAACCGAAGUUUCCC'       # two alternative folds
        base, _ = SS.predict_structures(seq, 1)
        other = '.' * len(seq)
        vals = [0.0 if ch != '.' else 2.0 for ch in base[0]['dbn']]
        flipped = [2.0 - v for v in vals]                       # data that contradict the MFE
        h, meta = SS.predict_structures(seq, 1, reactivities=flipped)
        self.assertIn('probing', meta)
        self.assertNotEqual(h[0]['dbn'], base[0]['dbn'])
        self.assertIn('rank_energy_kcal_mol', h[0])
        fc = RNA.fold_compound(seq)
        self.assertAlmostEqual(h[0]['energy_kcal_mol'], round(fc.eval_structure(h[0]['dbn']), 2), places=2)
        same, _ = SS.predict_structures(seq, 1, reactivities=vals)
        self.assertEqual(same[0]['dbn'], base[0]['dbn'])
        self.assertNotEqual(other, same[0]['dbn'])

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_alignment_consensus_is_the_first_hypothesis(self):
        h1 = list(TAR)
        h1[1], h1[27] = 'A', 'U'
        h1[3], h1[25] = 'G', 'C'
        hyps, meta = SS.predict_structures(TAR, 3, alignment=[TAR, ''.join(h1)])
        self.assertIn('alignment_consensus', hyps[0]['source'])
        self.assertEqual(meta['alignment']['sequences'], 2)
        self.assertEqual([h['id'] for h in hyps], [f'ss{k}' for k in range(1, len(hyps) + 1)])
        for i, j in hyps[0]['pairs']:
            self.assertIn(TAR[i] + TAR[j], ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'))


class TestModelOrder(unittest.TestCase):
    def test_best_model_of_each_template_comes_first(self):
        rng = np.random.default_rng(0)
        cands = []
        for t_rank in range(3):
            centre = rng.normal(0, 30, (20, 3))
            for pos in range(3):
                cands.append({'_c4': centre + rng.normal(0, 3.0, (20, 3)), 'total_score': 0.0,
                              'rank_key': (0, pos, t_rank, float(pos)), 'template': t_rank})
        reps = select_diverse(cands, 3, 2.0)
        self.assertEqual([cands[i]['template'] for i in reps], [0, 1, 2])


class TestInsight(unittest.TestCase):
    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_rescue_designs_restore_the_balance(self):
        doc = _docs()['trna']
        st = ensemble_states(doc)
        self.assertGreater(len(st), 1)
        x = insight.experiments(doc, st, 0, 1, Conditions(), doc_coords(doc))
        self.assertTrue(x['rescue'])
        seq = doc['sequence']
        own = {'A': set(pairs_of(st[0]['dbn'])) - set(pairs_of(st[1]['dbn'])), 'B': set(pairs_of(st[1]['dbn'])) - set(pairs_of(st[0]['dbn']))}
        for d in x['rescue']:
            self.assertIn(tuple(d['pair']), own[d['state']])
            self.assertEqual({d['first']['pos'], d['second']['pos']}, set(d['pair']))
            self.assertGreaterEqual(abs(d['ddg_first']), 1.0)
            self.assertLessEqual(abs(d['ddg_double']), 0.5 * abs(d['ddg_first']) + 1e-9)
            self.assertEqual((d['ddg_first'] > 0), d['state'] == 'A')       # breaking a helix of A favours B
            self.assertLess(abs(d['frac_a_double'] - d['frac_a_wt']), abs(d['frac_a_first'] - d['frac_a_wt']))
            m = d['sequence_double']
            self.assertEqual(sum(a != b for a, b in zip(seq, m)), 2)
            i, j = d['pair']
            self.assertIn(m[i] + m[j], ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'))
        seg = x['probing']['segments']
        self.assertTrue(seg and all(s['reactive_in'] in 'AB' for s in seg))
        for s in seg:                                   # reactive = unpaired in that state, paired in the other
            dbn = st[0 if s['reactive_in'] == 'A' else 1]['dbn']
            self.assertTrue(all(dbn[k] == '.' for k in range(s['start'], s['end'] + 1)))
        for lab in x['labels']:
            self.assertGreaterEqual(lab['j'] - lab['i'], 8)
            self.assertGreater(abs(lab['fret_a'] - lab['fret_b']), 0.15)
        with self.assertRaises(ValueError):
            insight.experiments(doc, st, 0, 0, Conditions())

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_variants(self):
        seq = TAR
        st = [{'dbn': '((((((...((((......))))))))))', 'conformations': []}]
        r = insight.variants(seq, st, Conditions(), 'G1C, 30a' if False else 'C13G, G1C+C29G')
        self.assertEqual([x['variant'] for x in r['rows']], ['C13G', 'G1C+C29G'])
        self.assertEqual(r['rows'][1]['positions'], [0, 28])
        self.assertEqual(len(r['rows'][0]['profile_diff']), len(seq))
        self.assertIn(r['rows'][0]['effect'], ('small', 'moderate', 'large'))
        # a compensated pair changes the fold less than a broken stem pair
        broken = insight.variants(seq, st, Conditions(), 'G2C')['rows'][0]
        fixed = insight.variants(seq, st, Conditions(), 'G2C+C28G')['rows'][0]
        self.assertLess(fixed['pairs_changed'], max(1, broken['pairs_changed']) + 1)
        self.assertLessEqual(fixed['change'], broken['change'] + 1e-6)
        for bad in ('X5A', 'A1C', 'G99A', 'G1G', ''):
            with self.assertRaises(ValueError):
                insight.parse_variants(bad, seq)
        self.assertEqual(insight.parse_variants('g1t', seq)[0][0], 'G1U')
        # 0.7.6: one mistyped entry no longer stops the others; 34G>U notation; a variant entered twice counts once
        r = insight.variants(seq, st, Conditions(), 'C13G, X5A, G99A, 13C>G, n.2G>C')
        self.assertEqual([x['variant'] for x in r['rows']], ['C13G', 'G2C'])
        self.assertEqual([k['text'] for k in r['skipped']], ['X5A', 'G99A'])
        self.assertIn('outside the sequence (1-29)', r['skipped'][1]['reason'])
        with self.assertRaises(ValueError):
            insight.variants(seq, st, Conditions(), 'X5A')
        with self.assertRaises(ValueError):
            insight.parse_variants('G1C+G1A', seq)
        # positions are the residue numbers RCV shows: a file numbered 17-45 is addressed with those numbers
        nums = [str(i) for i in range(17, 17 + len(seq))]
        v = insight.parse_variants('G17C, 29G', seq, numbers=nums)
        self.assertEqual([(x[0], x[2]) for x in v], [('G17C', [0]), ('C29G', [12])])
        with self.assertRaises(ValueError):
            insight.parse_variants('G1C', seq, numbers=nums)
        # numbers that are not unique (two chains joined) or carry insertion codes: positions count along the sequence
        self.assertEqual(insight.parse_variants('G1C', seq, numbers=['5'] * len(seq))[0][2], [0])
        self.assertEqual(insight.parse_variants('G1C', seq, numbers=[str(i) for i in range(1, len(seq) + 1)])[0][2], [0])

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_cotranscriptional(self):
        doc = _docs()['trna']
        st = ensemble_states(doc)
        x = insight.cotranscriptional(doc['sequence'], Conditions(), st, 0, 1)
        n = len(doc['sequence'])
        self.assertEqual(x['rows'][-1]['length'], n)
        self.assertEqual(x['rows'][-1]['dbn'], x['final'])
        self.assertEqual(x['rows'][-1]['transient'], 0)
        for r in x['rows']:
            self.assertEqual(len(r['code']), r['length'])
            self.assertEqual(len(r['dbn']), r['length'])
            self.assertEqual(r['code'].count('t'), 2 * r['transient'])
            self.assertTrue(0.0 <= r['a'] <= 1.0 and 0.0 <= r['b'] <= 1.0)
        self.assertTrue(all(e['broken'] >= 3 for e in x['events']))
        self.assertIn('peak_a', x)

    def test_motifs_of_trna_phe(self):
        ref = SAMPLES / '1EHZ_A' / 'reference.pdb'
        if not ref.exists():
            self.skipTest('benchmark sample 1EHZ not present')
        doc = load_document([ref])
        m = insight.motifs(doc, 0, doc_coords(doc))
        labels = {x['label'] for x in m['motifs']}
        self.assertEqual(m['counts'].get('hairpin loop'), 3)
        self.assertIn('4-way junction', labels)
        junction = next(x for x in m['motifs'] if x['type'] == 'junction')
        self.assertGreaterEqual(junction['coaxial'], 1)                 # acceptor stem on the T stem
        self.assertIn('U8·A14', labels)                                 # reverse Hoogsteen pair of the tRNA core
        self.assertIn('U54·A58', labels)                                # T-loop
        for x in m['motifs']:
            self.assertTrue(x['residues'] == sorted(x['residues']) and x['detail'])

    def test_motifs_names(self):
        self.assertEqual(insight._hairpin_name('GAAA'), 'GNRA tetraloop')
        self.assertEqual(insight._hairpin_name('UUCG'), 'UNCG tetraloop')
        self.assertEqual(insight._hairpin_name('AUC'), 'triloop')
        doc = _docs()['tar']
        m = insight.motifs(doc, 0, doc_coords(doc))
        self.assertTrue(any(x['type'] in ('bulge', 'internal loop') for x in m['motifs']))
        self.assertTrue(any(x['type'] == 'hairpin loop' for x in m['motifs']))

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_target_sites(self):
        doc = _docs()['trna']
        x = insight.target_sites(doc, Conditions(), 16)
        n = len(doc['sequence'])
        self.assertEqual(len(x['profile']), n - 16 + 1)
        self.assertTrue(all(v >= 0 for v in x['profile']))
        self.assertEqual(len(x['unpaired']), n)
        best = x['best']
        self.assertEqual(best[0]['dg_open'], min(x['profile']))
        comp = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
        for b in best:
            self.assertEqual(''.join(comp[ch] for ch in reversed(b['antisense_dna'])), doc['sequence'][b['start']:b['end'] + 1])
        self.assertTrue(all(abs(p['start'] - q['start']) >= 8 for k, p in enumerate(best) for q in best[:k]))
        with self.assertRaises(ValueError):
            insight.target_sites(doc, Conditions(), 4)

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_target_sites_one_pass_equals_constrained_folding(self):
        """0.7.6: the opening energy of every window comes from one pass over the ensemble; it must be the
        difference of the ensemble free energies with and without the window forced single-stranded."""
        import RNA
        doc = _docs()['trna']
        seq, width = doc['sequence'], 18
        for cond in (Conditions(), Conditions.from_json({'temperature': 25, 'sodium': 0.15}, len(seq))):
            x = insight.target_sites(doc, cond, width)
            md = cond.model()
            fc = RNA.fold_compound(seq, md)
            _s, mfe = fc.mfe()
            fc.exp_params_rescale(mfe)
            _, g0 = fc.pf()
            for s0 in (0, 7, 30, len(seq) - width):
                f2 = RNA.fold_compound(seq, md)
                for t in range(s0, s0 + width):
                    f2.hc_add_up(t + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
                f2.exp_params_rescale(mfe)
                _, g1 = f2.pf()
                self.assertAlmostEqual(x['profile'][x['starts'].index(s0)], max(0.0, g1 - g0), delta=0.011)
            self.assertEqual(x['starts'], list(range(len(seq) - width + 1)))

    def test_compare(self):
        d = _docs()
        X = doc_coords(d['tar'])
        same = insight.compare(d['tar'], X, d['tar'], X, 0, 0)
        self.assertLess(same['rmsd'], 1e-6)
        self.assertAlmostEqual(same['tm'], 1.0, places=6)
        self.assertAlmostEqual(same['local_agreement_mean'], 1.0, places=6)
        self.assertEqual(same['moved'], [])
        x = insight.compare(d['tar'], X, d['nmr'], doc_coords(d['nmr']))
        self.assertEqual(x['matched'], len(d['tar']['sequence']))
        self.assertEqual(len(x['displacement']), len(d['tar']['sequence']))
        self.assertGreater(x['rmsd'], 0.5)
        self.assertLessEqual(x['cross_rmsd_min'], x['rmsd'] + 1e-6)
        self.assertEqual(len(x['nearest_b']), min(40, len(d['nmr']['conformations'])))
        with self.assertRaises(ValueError):
            insight.compare(d['tar'], X, d['trna'], doc_coords(d['trna']))

    def test_compare_finds_the_moved_part_and_maps_a_deletion(self):
        d = _docs()['trna']
        X = doc_coords(d).copy()
        n = len(d['sequence'])
        first = np.array([r['first_atom'] for r in d['residues']] + [X.shape[1]])
        Y = X.copy()
        Y[0, first[60]:first[70]] += np.array([0.0, 25.0, 0.0])          # move residues 61-70 of conformation 1
        x = insight.compare(d, X, d, Y, 0, 0)
        self.assertTrue(x['moved'])
        self.assertTrue(55 <= x['moved'][0]['start'] <= 62 and 67 <= x['moved'][0]['end'] <= 71, x['moved'])
        self.assertLess(max(v for v in x['displacement'][:50]), 3.0)      # the rest stays superposed
        ia, ib = insight._sequence_map(d['sequence'], d['sequence'][:30] + d['sequence'][31:])
        self.assertEqual(len(ia), n - 1)
        self.assertEqual((ia[40], ib[40]), (41, 40))


class TestServer074(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()
        cls.ids = {}
        for key, name in (('trna', '1EHZ_tRNA_prediction'), ('tar', 'HIV1_TAR_prediction'), ('nmr', 'HIV1_TAR_NMR_1ANR.pdb')):
            st, raw = cls._req('POST', '/api/open', {'paths': [str(DEMO / name)]})
            cls.ids[key] = json.loads(raw)['id']

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

    def test_motifs_and_compare_routes(self):
        m = self.post('/api/insight/motifs', {'document': self.ids['trna'], 'conformation': 1})
        self.assertEqual(m['conformation'], 1)
        self.assertIn('motifs', m)
        c = self.post('/api/insight/compare', {'document': self.ids['tar'], 'other': self.ids['nmr'], 'ka': 0, 'kb': 2})
        self.assertEqual((c['ka'], c['kb']), (0, 2))
        self.assertEqual(self._req('POST', '/api/insight/compare', {'document': self.ids['tar']})[0], 400)
        self.assertEqual(self._req('POST', '/api/insight/compare', {'document': self.ids['tar'], 'other': self.ids['trna']})[0], 400)
        self.assertEqual(self._req('POST', '/api/insight/motifs', {'document': self.ids['tar'], 'conformation': 99})[0], 400)
        self.assertEqual(self._req('POST', '/api/insight/nothing', {'document': self.ids['tar']})[0], 404)
        self.assertEqual(self._req('POST', '/api/insight/motifs', {'document': 'missing'})[0], 400)

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_switch_and_pathway_routes(self):
        doc = self.ids['trna']
        e = self.post('/api/insight/experiments', {'document': doc, 'a': 0, 'b': 1, 'conditions': {'temperature': 30}})
        self.assertIn('30 °C', e['conditions_text'])
        self.assertTrue(e['rescue'])
        self.assertEqual(self._req('POST', '/api/insight/experiments', {'document': doc, 'a': 0, 'b': 0})[0], 400)
        v = self.post('/api/insight/variants', {'document': doc, 'text': 'G1C, 30A', 'a': 0, 'b': 1})
        self.assertEqual(len(v['rows']), 2)
        self.assertTrue(v['two_state'])
        st, raw = self._req('POST', '/api/insight/variants', {'document': doc, 'text': 'Z9Q'})
        self.assertEqual(st, 400)
        self.assertIn('cannot read variant', json.loads(raw)['error'])
        g = self.post('/api/insight/cotranscription', {'document': doc, 'conf_a': 0, 'conf_b': 1})
        self.assertEqual(g['length'], 76)
        states = ensemble_states(_docs()['trna'])
        if not any(0 in s['conformations'] and 1 in s['conformations'] for s in states):
            self.assertIn('peak_a', g)
            self.assertIn('peak_b', g)
        g0 = self.post('/api/insight/cotranscription', {'document': doc})
        self.assertNotIn('peak_a', g0)
        t = self.post('/api/insight/targets', {'document': doc, 'width': 20})
        self.assertEqual(t['width'], 20)
        self.assertEqual(len(t['best'][0]['antisense_dna']), 20)


if __name__ == '__main__':
    unittest.main()
