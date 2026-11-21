"""0.7.6: modified nucleotides (rnafold3d/modifications.py) — sequence notation, ViennaRNA soft/hard constraints for
the 2D, 3D decoration of the parent bases, PDB/mmCIF output, the RCV document, jobs and the server API."""
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

from rnafold3d import modifications as M          # noqa: E402

try:
    import RNA  # noqa: F401
    HAVE_VIENNA = True
except ImportError:                                  # pragma: no cover
    HAVE_VIENNA = False


class TestNotation(unittest.TestCase):
    def test_catalogue_is_consistent(self):
        self.assertGreaterEqual(len(M.ENTRIES), 90)
        codes = [e['code'] for e in M.ENTRIES]
        self.assertEqual(len(codes), len(set(codes)))
        for e in M.ENTRIES:
            self.assertIn(e['parent'], 'ACGUN', e['code'])
            self.assertIn(e['category'], M.CATEGORIES, e['code'])
            self.assertTrue(e['name'] and e['note'], e['code'])
        cat = M.catalogue()
        self.assertEqual({r['code'] for r in cat['entries']}, set(codes))
        json.dumps(cat)

    def test_lookup_by_alias_ccd_and_name(self):
        self.assertEqual(M.lookup('m6A')['code'], 'm6A')
        self.assertEqual(M.lookup('6MZ')['code'], 'm6A')
        self.assertEqual(M.lookup('PSU')['code'], 'Y')
        self.assertEqual(M.lookup('Ψ')['code'], 'Y')
        self.assertEqual(M.lookup('A>I')['code'], 'I')
        self.assertEqual(M.lookup('8oxoG')['code'], '8oxoG')
        self.assertIsNone(M.lookup('not-a-modification'))

    def test_parse_format_round_trip(self):
        text = 'GGC[m6A]CU[Y]A*G[I]Cu'
        seq, mods, problems = M.parse(text)
        self.assertEqual(problems, [])
        self.assertEqual(seq, 'GGCACUUAGGCU')            # inosine is built on G
        self.assertEqual(M.normalize(mods, seq), mods)
        self.assertIn({'pos': 4, 'code': 'm6A'}, mods)
        self.assertIn({'pos': 9, 'code': 'ps'}, mods)
        self.assertEqual(M.format_sequence(seq, mods), 'GGC[m6A]CU[Y]A*G[I]CU')

    def test_problems(self):
        for text, word in [('GG[xyz]C', 'unknown'), ('GG[m6A', 'unclosed'), ('*GGC', 'between'), ('GGC*', 'end'),
                           ('GGXC', 'letter'), ('G[ps]C', 'linkage')]:
            self.assertTrue(any(word in p for p in M.parse(text)[2]), text)

    def test_normalize_checks_the_base(self):
        with self.assertRaisesRegex(ValueError, 'needs A'):
            M.normalize([{'pos': 1, 'code': 'm6A'}], 'GGC')
        with self.assertRaisesRegex(ValueError, 'outside'):
            M.normalize([{'pos': 9, 'code': 'm6A'}], 'AAA')
        with self.assertRaisesRegex(ValueError, 'two modifications'):
            M.normalize([{'pos': 1, 'code': 'm6A'}, {'pos': 1, 'code': 'm1A'}], 'AAA')
        # a sugar modification and a phosphorothioate may share a position
        self.assertEqual(len(M.normalize([{'pos': 2, 'code': 'Am'}, {'pos': 2, 'code': 'ps'}], 'AAA')), 2)


@unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
class TestEnergies(unittest.TestCase):
    HAIRPIN = 'GGGAGCUUCGGCUCCC'

    def _mfe(self, seq, mods):
        fc, info = M.fold_compound(seq, mods)
        s, e = fc.mfe()
        return s, e, info

    def test_blocked_base_stays_unpaired(self):
        s0, _, _ = self._mfe(self.HAIRPIN, [])
        self.assertEqual(s0[1], '(')
        s, _, info = self._mfe(self.HAIRPIN, [{'pos': 2, 'code': 'm1G'}])
        self.assertEqual(s[1], '.')
        self.assertEqual(info['blocked'], [2])

    def test_measured_parameters_change_the_energy(self):
        seq = 'GGGAUUCUUCGGAAUCCC'
        _, e0, _ = self._mfe(seq, [])
        _, e1, info = self._mfe(seq, [{'pos': 5, 'code': 'Y'}, {'pos': 6, 'code': 'Y'}])
        self.assertEqual(info['measured'], [5, 6])
        self.assertLess(e1, e0)                          # pseudouridine stabilises the stem

    def test_8oxoG_energies_follow_the_melting_data(self):
        """Choi et al., Nucleic Acids Res. 2017, 45:2099, Table 1: the 12-mer duplex AAGAGGGAUGAC with 8-oxoG at G6,
        G6+G7 and G5+G6+G7, and with A opposite the lesion. Free-energy differences to the unmodified duplex from
        the melting temperatures (two-state model, nearest-neighbour enthalpy): +3.0, +6.1, +10.4 and +3.8 kcal/mol.
        The pair with C follows these numbers. The pair with A does not: by the author's decision it is the
        strongest pair of 8-oxoG (energies of a U-A pair); the wobble pair with U is the weakest."""
        full = '(' * 12 + ')' * 12

        def energy(seq, sites):
            fc, _ = M.fold_compound(seq, [{'pos': k, 'code': '8oxoG'} for k in sites])
            return fc.eval_structure(full)
        c, a = 'AAGAGGGAUGAC&GUCAUCCCUCUU', 'AAGAGGGAUGAC&GUCAUCACUCUU'
        ref = energy(c, [])
        for seq, sites, measured, tol in ((c, [6], 3.0, 0.5), (c, [6, 7], 6.1, 0.5), (c, [5, 6, 7], 10.4, 1.5)):
            self.assertAlmostEqual(energy(seq, sites) - ref, measured, delta=tol, msg=f'{sites} in {seq}')
        u = 'AAGAGGGAUGAC&GUCAUCUCUCUU'
        with_c, with_a, with_u = energy(c, [6]), energy(a, [6]), energy(u, [6])
        self.assertLess(with_a, with_c - 0.4)                # A is the dominant partner ...
        self.assertLess(with_c, with_u - 1.0)                # ... and U the weakest
        self.assertAlmostEqual(with_u - energy(u, []), 3.0, delta=0.05)      # the wobble is weakened like the pair with C
        e = M.BY_CODE['8oxoG']['energy'][1]
        for up in M.CANON:                                   # in every neighbour context
            for down in M.CANON:
                two = {p: M._stack(up, ref) + M._stack(ref, down) + 2 * d for p, (ref, d) in e.items()}
                self.assertLess(two['A'], two['C'], (up, down))
                self.assertLess(two['C'], two['U'], (up, down))
        # the stack shared by two lesions carries both corrections, whichever strands the lesions are in
        for seq, sites in (('GGAGCGUCC&GGACGCUCC', [4, 14]), ('GGACGAUCC&GGAUCGUCC', [5, 15])):
            fc, _ = M.fold_compound(seq, [{'pos': k, 'code': '8oxoG'} for k in sites])
            fc0, _ = M.fold_compound(seq, [])
            s9 = '(' * 9 + ')' * 9
            self.assertAlmostEqual(fc.eval_structure(s9) - fc0.eval_structure(s9), 6.0, delta=0.05)

    def test_8oxoG_pairs_with_A(self):
        seq = 'GGGAGCUUCGGAUCCC'                         # G5 faces A12: a mismatch for G
        mods = [{'pos': 5, 'code': '8oxoG'}]
        info = M.fold_compound(seq, mods)[1]
        self.assertEqual(info['estimated'], [5])
        _, e_g, _ = self._mfe(seq, [])
        _, e_ox, _ = self._mfe(seq, mods)
        self.assertNotEqual(e_g, e_ox)

    def test_inosine_pairs_with_C(self):
        info = M.fold_compound('GGGIGCUUCGGCCCCC'.replace('I', 'G'), [{'pos': 4, 'code': 'I'}])[1]
        self.assertEqual(info['measured'], [4])

    def test_split_pairs(self):
        seq = 'GGGAGCUUCGGAUCCC'
        pairs = [(0, 15), (1, 14), (2, 13), (4, 11)]
        rigid, special, dropped = M.split_pairs(seq, pairs, [{'pos': 5, 'code': '8oxoG'}, {'pos': 2, 'code': 'm1G'}])
        self.assertIn((4, 11), special)
        self.assertIn((1, 14), dropped)
        self.assertIn((0, 15), rigid)


class TestPrediction(unittest.TestCase):
    """One small prediction with modifications through the whole engine, then the RCV document."""

    @classmethod
    def setUpClass(cls):
        if not HAVE_VIENNA:
            raise unittest.SkipTest('ViennaRNA not installed')
        from rnafold3d.predict import predict, PredictConfig
        cls.tmp = Path(tempfile.mkdtemp())
        cls.seq_text = 'GG[m6A]CUUC[Y]GAAG[m1G]AG*U[Gm]CC'
        cls.man = predict(cls.seq_text, cls.tmp / 'run', config=PredictConfig(
            n_models=1, trajectories=2, workers=1, templates='off', motifs='off', pseudoknots='off'), log=None)
        cls.out = cls.tmp / 'run'

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_manifest(self):
        mods = self.man['input']['modifications']
        codes = [m['code'] for m in mods['sites']]
        self.assertEqual(sorted(codes), sorted(['m6A', 'Y', 'm1G', 'ps', 'Gm']))
        self.assertTrue((self.out / 'input' / 'sequence_modified.txt').exists())

    def _pdb(self):
        files = sorted(self.out.rglob('*.pdb'))
        self.assertTrue(files)
        return files[0].read_text()

    def test_pdb_records_and_geometry(self):
        text = self._pdb()
        self.assertIn('MODRES', text)
        het = [l for l in text.splitlines() if l.startswith('HETATM')]
        comps = {l[17:20].strip() for l in het}
        for c in ('6MZ', 'PSU', '1MG', 'OMG'):
            self.assertIn(c, comps)
        xyz = {}
        for l in het:
            xyz.setdefault((l[17:20].strip(), int(l[22:26])), {})[l[12:16].strip()] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        m6a = next(v for (c, _), v in xyz.items() if c == '6MZ')
        self.assertAlmostEqual(float(np.linalg.norm(m6a['C9'] - m6a['N6'])), 1.45, delta=0.08)   # C9: the methyl of 6MZ
        omg = next(v for (c, _), v in xyz.items() if c == 'OMG')
        self.assertAlmostEqual(float(np.linalg.norm(omg['CM2'] - omg["O2'"])), 1.42, delta=0.08)
        psu = next(v for (c, _), v in xyz.items() if c == 'PSU')
        self.assertAlmostEqual(float(np.linalg.norm(psu['C5'] - psu["C1'"])), 1.50, delta=0.12)  # C-glycosidic bond
        sp = [l for l in text.splitlines() if l[12:16].strip() == 'S1P']
        self.assertEqual(len(sp), 1)

    def test_mmcif(self):
        cifs = sorted(self.out.rglob('*.cif'))
        if not cifs:
            self.skipTest('no mmCIF written')
        text = cifs[0].read_text()
        self.assertIn('mon_nstd_parent_comp_id', text)
        self.assertIn('(6MZ)', text)

    def test_document(self):
        from rcv.data import load_document
        doc = load_document([self.out])
        mods = {m['index']: m for m in doc['modifications']}
        self.assertEqual(mods[2]['code'], 'm6A')
        self.assertEqual(mods[2]['category'], 'methyl')
        self.assertEqual(doc['residues'][2].get('mod'), 'm6A')
        self.assertTrue(any(m['code'] == 'Y' for m in mods.values()))
        self.assertEqual(doc['sequence'][:5], 'GGACU')


class TestModificationsInTheSearch(unittest.TestCase):
    """Pucker and syn/anti preferences, modification atoms as riders of the sampled model, local relaxation."""
    SEQ = 'GGGAGCUUCGGCUCCCAAUGAGGACU'
    PAIRS = [(0, 15), (1, 14), (2, 13), (3, 12)]

    def _model(self, want):
        from rnafold3d.fold import FoldModel
        mods = M.normalize([{'pos': p, 'code': c} for c, p in want.items()], self.SEQ)
        return FoldModel(self.SEQ, self.PAIRS, modifications=mods), mods

    def test_preferences(self):
        pr = M.sampling_preferences([{'pos': 1, 'code': 'D'}, {'pos': 2, 'code': '+A'}, {'pos': 3, 'code': '8oxoG'},
                                     {'pos': 4, 'code': 'Um'}, {'pos': 5, 'code': 'm6A'}])
        self.assertEqual(pr['pucker'][0][0], 'C2endo')
        self.assertEqual(pr['pucker'][1], ('C3endo', 1.0))
        self.assertEqual(pr['chi'][2][0], 'syn')
        self.assertEqual(pr['chi'][1][0], 'anti')
        self.assertNotIn(4, pr['pucker'])

    def test_riders_follow_the_parent_atoms(self):
        m, mods = self._model({'m2,2G': 5, 'Um': 7, 'm6A': 17, '8oxoG': 20, 'ps': 22, 'dA': 24, 'm1Y': 26})
        names = {(r['res'] + 1, r['name']) for r in m.riders}
        for key in [(5, 'CN2'), (7, 'CM2'), (17, 'CM6'), (20, 'O8'), (22, 'S1P'), (26, 'CM1')]:
            self.assertIn(key, names)
        self.assertTrue(m.ghost[m.idx[(23, "O2'")]])          # dA: no O2'
        paired = {i for p in self.PAIRS for i in p}
        for seed in range(3):
            m.sample_torsions(np.random.default_rng(seed))
            m.set_loops_active(True)
            m.build()
            m.update_riders()
            res = m.residues_atoms()
            for rd in m.riders:
                r = rd['res']
                at = {k: v.copy() for k, v in res[r].items()}
                if rd['group'] == 'phosphate':
                    ref = at['P'] + 1.96 * (at['OP1'] - at['P']) / np.linalg.norm(at['OP1'] - at['P'])
                else:
                    ops = [op for x in mods if x['pos'] - 1 == r for op in M.BY_CODE[x['code']]['atoms']
                           if not (op[0] == 'swap' and op[3] == 'P')]
                    M._place(at, ops, None, r in paired)
                    ref = at[rd['name']]
                self.assertLess(float(np.linalg.norm(ref - m.xyz[rd['k']])), 0.05, (seed, r + 1, rd['name']))
        # the riders take part in the contact scoring
        n, ov = m.rider_contacts()
        self.assertGreaterEqual(n, 0)

    def test_draws_follow_the_preferences(self):
        m, _ = self._model({'D': 8, '+A': 17, '8oxoG': 20})
        rng = np.random.default_rng(0)
        c2 = sum(m.pucker_of_row(7, m.draw_row(rng, 7)) == 'C2endo' for _ in range(400)) / 400
        self.assertGreater(c2, 0.55)
        self.assertTrue(all(m.pucker_of_row(16, m.draw_row(rng, 16)) == 'C3endo' for _ in range(100)))
        syn = sum(-20 <= m.draw_row(rng, 19)[6] <= 120 for _ in range(400)) / 400
        self.assertGreater(syn, 0.45)

    def test_unmodified_model_is_unchanged(self):
        from rnafold3d.fold import FoldModel
        m = FoldModel(self.SEQ, self.PAIRS)
        self.assertFalse(m.riders)
        lib = m.geom.loop_library
        a, b = np.random.default_rng(5), np.random.default_rng(5)
        for _ in range(20):
            self.assertTrue(np.array_equal(m.draw_row(a, 7), lib[b.integers(len(lib))]))

    def test_chi_prior_prefers_syn(self):
        from rnafold3d.build import load_geometry
        from rnafold3d.loops import TorsionPrior
        P = TorsionPrior(load_geometry())
        t = P.chi_table('syn', 1.5)
        self.assertLess(t[62 + 180], t[-160 + 180] - 1.0)
        a = P.chi_table('anti', 1.0)
        self.assertGreater(a[62 + 180], P('chi', 62.0) - 1e-9)

    @unittest.skipUnless(HAVE_VIENNA, 'ViennaRNA not installed')
    def test_trajectory_with_a_syn_purine_in_a_loop(self):
        from dataclasses import asdict
        from rnafold3d.predict import _run_trajectory, PredictConfig
        s, mods, _ = M.parse('GGCGCAAG[8oxoG]AAGCGCC')
        pairs = [(0, 15), (1, 14), (2, 13), (3, 12)]
        cfg = asdict(PredictConfig(templates='off', motifs='off'))
        cfg['modifications'] = M.normalize(mods, s)
        out = _run_trajectory({'sequence': s, 'hypothesis': {'id': 'ss1', 'pairs': pairs}, 'index': 0, 'seed': 3,
                               'config': cfg})
        self.assertEqual(out['status'], 'valid', out.get('reasons'))
        conf = out['modification_conformations'][0]
        self.assertEqual(conf['code'], '8oxoG')
        self.assertEqual(conf['glycosidic'], 'syn')
        self.assertIn('contacts_after', out['modification_relax'])


class TestModificationInteractions(unittest.TestCase):
    """Hydrogen bonds, the pseudouridine water bridge and stacking made by the modifications (modinteract.py)."""
    SEQ = 'GGGAGCUUCGGCUCCCAAUGAGGACUGCAU'
    PAIRS = [(0, 15), (1, 14), (2, 13), (3, 12)]

    def _mi(self, want):
        from rnafold3d.fold import FoldModel
        from rnafold3d.modinteract import ModInteractions
        mods = M.normalize([{'pos': p, 'code': c} for c, p in want.items()], self.SEQ)
        m = FoldModel(self.SEQ, self.PAIRS, modifications=mods)
        return m, ModInteractions(m)

    def _roles(self, m, mi, r, name):
        k = m.idx[(r - 1, name)]
        if k not in set(mi.p_idx.tolist()):
            return None
        j = mi.p_idx.tolist().index(k)
        return bool(mi.p_acc[j]), bool(mi.p_don[j])

    def test_roles_follow_the_chemistry(self):
        m, mi = self._mi({'Y': 7, 'm1Y': 8, 't6A': 17, '8oxoG': 20, 'm1G': 22, 'm2,2G': 5, 'D': 26})
        q = {(t[5] + 1, t[6]): (t[1], t[2]) for t in mi.q}
        self.assertEqual(q[(7, 'N1')], (True, False))          # Ψ: new N1-H donor
        self.assertNotIn((8, 'N1'), q)                         # m1Ψ: the methyl takes it
        self.assertEqual(q[(20, 'N7')], (True, False))         # 8-oxo-G: N7-H
        self.assertEqual(q[(17, 'N11')], (True, False))        # t6A ureido N-H
        self.assertEqual(q[(17, 'O13')], (True, True))         # threonine hydroxyl
        self.assertIsNone(self._roles(m, mi, 22, 'N1'))        # m1G: no N1-H
        self.assertIsNone(self._roles(m, mi, 5, 'N2'))         # m2,2G: no amino H left
        self.assertEqual([int(m.atom_res[k]) + 1 for k, _, _ in mi.wb], [7])
        stack = dict((r + 1, s) for r, s in mi.stack)
        self.assertLess(stack[26], 0)                          # dihydrouridine stacks less
        self.assertGreater(stack[7], 0)

    def test_helical_residues_have_no_stacking_term(self):
        from rnafold3d.fold import FoldModel
        from rnafold3d.modinteract import ModInteractions
        mods = M.normalize([{'pos': 13, 'code': 'Y'}], self.SEQ)          # U13 pairs with A4
        mi = ModInteractions(FoldModel(self.SEQ, self.PAIRS, modifications=mods))
        self.assertEqual(mi.stack, [])                          # its stacking is in the 2D parameters

    def test_energy_and_report(self):
        from rnafold3d.search import Scorer
        m, _ = self._mi({'Y': 7, 't6A': 17, 'D': 26})
        m.sample_torsions(np.random.default_rng(0))
        m.set_loops_active(True)
        m.build()
        sc = Scorer(m)
        self.assertIsNotNone(sc.modint)
        d = sc(detail=True)
        self.assertIn('modification', d)
        self.assertLessEqual(d['modification'], 0.0)
        rep = sc.modint.report(sc.base_centroids())
        self.assertEqual({r['residue'] for r in rep} >= {7, 17, 26}, True)
        json.dumps(rep)
        from rnafold3d.fold import FoldModel
        self.assertIsNone(Scorer(FoldModel(self.SEQ, self.PAIRS)).modint)

    def test_stacking_geometry(self):
        """Two stacked helix bases count as one stacked face."""
        from rnafold3d.fold import FoldModel
        from rnafold3d.modinteract import ModInteractions
        from rnafold3d.search import Scorer
        mods = M.normalize([{'pos': 17, 'code': 'm6A'}], self.SEQ)
        m = FoldModel(self.SEQ, self.PAIRS, modifications=mods)
        m.sample_torsions(np.random.default_rng(0))
        m.set_loops_active(True)
        m.build()
        mi = ModInteractions(m)
        mi.stack = [(1, 0.5)]                                   # probe a helical residue: stacked on both faces
        (_, _, faces, partners), = mi._stacking(Scorer(m).base_centroids())
        self.assertGreater(faces, 1.5)
        self.assertEqual(sorted(j for j, _ in partners)[:2], [0, 2])


class TestTemplateFit(unittest.TestCase):
    """Engine 2.5: a fit onto the template atoms that fails validation is repeated with the residues at fault
    released; models whose fit stays rejected rank after those with an accepted fit."""

    def test_release_grows_around_the_issues(self):
        from rnafold3d.predict import _release_for_refit, TEMPLATE_FIT_MAX_RELEASED
        issues = [{'residue': 11, 'type': 'chain_break'}, {'residue': 40, 'type': 'bond_angle'}]      # 1-based
        r0 = _release_for_refit(issues, 0, set(), 100)
        self.assertEqual(r0, {9, 10, 11, 38, 39, 40})
        r1 = _release_for_refit(issues, 1, r0, 100)
        self.assertEqual(r1, set(range(8, 13)) | set(range(37, 42)))
        self.assertIsNone(_release_for_refit(issues, 0, r1, 100))                     # nothing new to release
        self.assertEqual(_release_for_refit([{'residue': 1}], 0, set(), 100), {0, 1})  # chain ends are respected
        self.assertEqual(_release_for_refit([{'residue': 100}], 0, set(), 100), {98, 99})
        # too much of the chain: given up (the model then keeps its first fit, marked as not accepted)
        many = [{'residue': k} for k in range(1, 30, 3)]
        self.assertIsNone(_release_for_refit(many, 0, set(), 60))
        self.assertLessEqual(len(_release_for_refit(many, 0, set(), 100)), TEMPLATE_FIT_MAX_RELEASED * 100)

    def test_models_with_a_rejected_fit_come_last(self):
        from rnafold3d.predict import _order_template_models
        mk = lambda score, accepted: {'score_3d': score, 'template_fit': {'accepted': accepted}}
        by_t = {0: [mk(-50.0, False), mk(-10.0, True), mk(-30.0, True)], 1: [mk(-5.0, True), mk(-80.0, False)]}
        _order_template_models(by_t, True)
        order = sorted((r for rs in by_t.values() for r in rs), key=lambda r: r['rank_key'])
        self.assertEqual([(r['rank_key'][3], r['score_3d']) for r in order],
                         [(0, -30.0), (1, -5.0), (0, -10.0), (1, -80.0), (0, -50.0)])
        self.assertTrue(all(r['template_fit']['accepted'] for r in order[:3]))
        # a single fit attempt (template_fit_attempts=1): the order of 0.7.5, by score within each template
        by_t = {0: [mk(-50.0, False), mk(-10.0, True)], 1: [mk(-5.0, True), mk(-80.0, False)]}
        _order_template_models(by_t, False)
        order = sorted((r for rs in by_t.values() for r in rs), key=lambda r: r['rank_key'])
        self.assertEqual([(r['rank_key'][3], r['score_3d']) for r in order], [(0, -50.0), (1, -80.0), (0, -10.0), (1, -5.0)])
        # a model without a fit (no template atoms) is not treated as rejected
        by_t = {0: [{'score_3d': 1.0}, mk(-9.0, False)]}
        _order_template_models(by_t, True)
        self.assertEqual(by_t[0][0]['rank_key'][:3], (0, 0, 0))


class TestJobsAndAlphaFold(unittest.TestCase):
    def test_sequence_with_modifications(self):
        from rcv.jobs import sequence_with_modifications
        seq, mods = sequence_with_modifications('gg[m6A]cu')
        self.assertEqual(seq, 'GGACU')
        self.assertEqual(mods, [{'pos': 3, 'code': 'm6A'}])
        with self.assertRaisesRegex(ValueError, 'unknown modification'):
            sequence_with_modifications('GG[nope]C')
        with self.assertRaisesRegex(ValueError, 'unsupported letters'):
            sequence_with_modifications('GGXC')
        seq, mods = sequence_with_modifications('GGACU', [{'pos': 4, 'code': 'm5C'}])
        self.assertEqual(mods, [{'pos': 4, 'code': 'm5C'}])

    def test_fasta_keeps_tokens(self):
        from rcv.jobs import parse_fasta
        recs = parse_fasta('>a\nGG[m6A]\nCU\n>b\nACGU\n')
        self.assertEqual(recs[0]['sequence'], 'GG[m6A]CU')

    def test_alphafold_input(self):
        from rcv.alphafold import input_json
        js = input_json('x', 'GGACU', modifications=[{'pos': 3, 'code': 'm6A'}])
        rna = js['sequences'][0]['rna']
        self.assertEqual(rna['modifications'], [{'modificationType': '6MZ', 'basePosition': 3}])


class TestServerApi(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from rcv import server as srv
        cls.tmp = tempfile.mkdtemp()
        state = srv.AppState(Path(cls.tmp) / 'runs', [], [])
        cls.state = state
        cls.server = srv.AppServer(('127.0.0.1', 0), state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def req(self, method, path, body=None):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=60)
        headers = {'Host': f'localhost:{self.port}', 'X-Viewer-Token': self.state.token}
        data = None
        if body is not None:
            data = json.dumps(body).encode()
            headers['Content-Type'] = 'application/json'
        c.request(method, path, body=data, headers=headers)
        r = c.getresponse()
        raw = r.read()
        c.close()
        return r.status, json.loads(raw or b'null')

    def test_catalogue_and_parse(self):
        st, cat = self.req('GET', '/api/modifications')
        self.assertEqual(st, 200)
        self.assertIn('methyl', cat['categories'])
        self.assertTrue(any(e['code'] == '8oxoG' for e in cat['entries']))
        st, r = self.req('POST', '/api/modifications/parse', {'sequence': 'GG[m6A]CU[nope]'})
        self.assertEqual(st, 200)
        self.assertEqual(r['modifications'][0], {'pos': 3, 'code': 'm6A'})
        self.assertTrue(r['problems'])


class TestFrontendFiles(unittest.TestCase):
    def test_editor_strings_are_translated(self):
        i18n = ROOT / 'src' / 'rcv' / 'static' / 'i18n'
        for f in i18n.glob('*.json'):
            d = json.loads(f.read_text())['strings']
            for k in ('Modified nucleotides…', 'Use in the sequence', 'm6A at DRACH motifs', 'blocks pairing'):
                self.assertIn(k, d, f.name)


if __name__ == '__main__':
    unittest.main()
