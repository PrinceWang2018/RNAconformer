"""0.7.6: Z-RNA - identification of left-handed helices (rnafold3d/zrna.py, rcv/static/js/zrna.js) and their
prediction by the engine (Z-form stems, A->Z estimate, Z hypothesis)."""
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rnafold3d import zrna as Z
from rnafold3d.build import ResidueBuilder, load_geometry, build_stem
from rnafold3d.structure_io import read_models

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / 'tests' / 'data' / 'zrna_2gxb_duplex.pdb'
DUPLEX = [(i, 11 - i) for i in range(6)]


def load_fixture():
    _, chains, _ = read_models(FIXTURE)[0]
    ch = list(chains.values())[0]
    return ''.join(r.base for r in ch.residues), [{k: np.asarray(v, float) for k, v in r.atoms.items()} for r in ch.residues]


def built(seq, pairs, form):
    out, rep = build_stem(ResidueBuilder(load_geometry()), seq, pairs, form=form)
    return [out.get(r, {}) for r in range(len(seq))], rep


class TestDetection(unittest.TestCase):
    def test_crystal_z_rna_is_found(self):
        seq, res = load_fixture()
        regs = Z.detect(seq, res, DUPLEX)
        self.assertEqual(len(regs), 1)
        z = regs[0]
        self.assertEqual(len(z['pairs']), 6)
        self.assertEqual(z['syn_purines'], 1.0)
        self.assertLess(z['mean_twist'], -15)
        # the dinucleotide repeat: one small and one large left-handed step
        self.assertTrue(all(t < -40 for t in z['twists'][1::2]))
        json.dumps(regs)

    def test_a_form_is_not_z(self):
        seq = 'CGCGCGCGCGCG'
        res, _ = built(seq, DUPLEX, 'A')
        self.assertEqual(Z.detect(seq, res, DUPLEX), [])
        twists = [Z.step_twist(res[i], res[j], res[k], res[l]) for (i, j), (k, l) in zip(DUPLEX, DUPLEX[1:])]
        self.assertTrue(all(30 < t < 37 for t in twists))

    def test_needs_syn_purines_and_three_pairs(self):
        seq, res = load_fixture()
        self.assertEqual(Z.detect(seq, res, DUPLEX[:2]), [])          # two pairs are never enough
        chi = [-160.0] * 12                                           # all anti: a left-handed twist alone is not Z
        self.assertEqual(Z.detect(seq, res, DUPLEX, chi=chi), [])

    def test_broken_pairs_are_ignored(self):
        seq, res = load_fixture()
        moved = [dict(a) for a in res]
        for k in moved[11]:                                           # pull residue 12 away from its partner
            moved[11][k] = moved[11][k] + np.array([12.0, 0, 0])
        regs = Z.detect(seq, moved, DUPLEX)
        self.assertTrue(all((1, 12) not in [tuple(p) for p in z['pairs']] for z in regs))


class TestZStem(unittest.TestCase):
    def test_built_z_stem_matches_the_crystal(self):
        seq, ref = load_fixture()
        res, rep = built(seq, DUPLEX, 'Z')
        self.assertEqual(rep['form'], 'Z')
        self.assertLess(rep['pair_fit_rmsd'], 0.3)
        self.assertGreaterEqual(rep['min_contact'] or 9, 2.5)
        regs = Z.detect(seq, res, DUPLEX)
        self.assertEqual(len(regs[0]['pairs']), 6)
        from rnafold3d.geomutil import kabsch
        common = [(r, n) for r in range(12) for n in res[r] if n in ref[r] and n not in ('OP1', 'OP2')]
        A = np.array([res[r][n] for r, n in common])
        B = np.array([ref[r][n] for r, n in common])
        R, t = kabsch(A, B)
        rmsd = float(np.sqrt(np.mean(np.sum((A @ R.T + t - B) ** 2, axis=1))))
        self.assertLess(rmsd, 2.0)                                    # 1.4 A; A-form is 8 A away

    def test_fold_model_builds_z_helix(self):
        from rnafold3d.fold import FoldModel
        seq = 'CGCGCGCGAAAACGCGCGCG'
        pairs = [(i, 19 - i) for i in range(8)]
        m = FoldModel(seq, pairs, z_pairs=pairs)
        self.assertEqual(m.z_residues, {r for p in pairs for r in p})
        self.assertEqual(m.puck[1], 'C3endo')                       # syn G
        self.assertEqual(m.puck[0], 'C2endo')                       # anti C
        self.assertEqual(m.stem_reports[0]['form'], 'Z')
        self.assertEqual(FoldModel(seq, pairs).stem_reports[0]['form'], 'A')


class TestPropensity(unittest.TestCase):
    HAIRPIN = 'CGCGCGAAAACGCGCG'
    PAIRS = [(i, 15 - i) for i in range(6)]

    def test_salt_modifications_and_binder(self):
        dg = lambda **kw: Z.candidates(self.HAIRPIN, self.PAIRS, **kw)[0]['dG_kcal_mol']
        self.assertGreater(dg(salt=0.15), 4.0)                        # A-form at physiological salt
        self.assertLess(dg(salt=6.0), 0.5)                            # left-handed near 6 M (Hall 1984)
        m8g = [{'pos': 4, 'code': 'm8G'}, {'pos': 14, 'code': 'm8G'}]
        self.assertLess(dg(salt=0.15, mods=m8g), dg(salt=0.15) - 5)
        self.assertLess(dg(salt=0.15, binder=True), Z.Z_PROPOSE_MAX)

    def test_only_alternating_watson_crick_segments(self):
        self.assertEqual(Z.candidates('GGGGGGAAAACCCCCC', self.PAIRS, salt=6.0), [])     # no alternation
        self.assertEqual(Z.candidates('UGUGUGAAAAGCGCGC', self.PAIRS, salt=6.0), [])     # G-U wobbles
        ca = Z.candidates('CACACAAAAAUGUGUG', self.PAIRS, salt=6.0)[0]                    # CA/UG repeat
        self.assertEqual(len(ca['pairs']), 6)
        self.assertGreater(ca['dG_kcal_mol'], Z.candidates(self.HAIRPIN, self.PAIRS, salt=6.0)[0]['dG_kcal_mol'])

    def test_junctions_cost_and_open_a_pair(self):
        seq = 'GGAGCGCGCGCGAAAGCGCGCGCGCUCC'
        pairs = [(i, 27 - i) for i in range(12)]
        c = Z.candidates(seq, pairs, salt=6.0)[0]
        self.assertEqual(c['pairs'][0], (3, 24))                     # GCGC... from pair 4 (A3 breaks the alternation)
        self.assertEqual(c['n_junctions'], 1)                        # the other end is the hairpin loop
        h = Z.hypothesis(seq, {'id': 'ss1', 'pairs': pairs, 'dbn': '', 'energy_kcal_mol': -20.0, 'source': 'mfe'}, c)
        self.assertEqual(h['zrna']['opened_pairs'], [(3, 26)])         # the pair at the A-Z junction opens
        self.assertIn((0, 27), h['pairs'])
        self.assertTrue(set(h['z_pairs']) <= set(map(tuple, h['pairs'])))
        self.assertEqual(h['dbn'].count('('), len(h['pairs']))
        self.assertAlmostEqual(h['energy_kcal_mol'], -20.0 + c['dG_kcal_mol'], places=2)

    def test_m8g_is_a_syn_modification(self):
        from rnafold3d import modifications as M
        mods = M.normalize([{'pos': 4, 'code': 'm8G'}, {'pos': 6, 'code': '8BrG'}], self.HAIRPIN)
        chi = M.sampling_preferences(mods)['chi']
        self.assertEqual(chi[3][0], 'syn')
        self.assertEqual(chi[5][0], 'syn')


class TestPrediction(unittest.TestCase):
    def test_z_hypothesis_end_to_end(self):
        from rnafold3d.predict import predict, PredictConfig
        tmp = tempfile.mkdtemp()
        try:
            cfg = PredictConfig(n_models=4, trajectories=4, templates='off', motifs='off', zrna='auto',
                                zrna_binder=True, anneal_steps=300)
            predict('CGCGCGCGAAAACGCGCGCG', tmp, name='z', config=cfg, log=lambda *a: None)
            man = json.loads((Path(tmp) / 'metadata' / 'manifest.json').read_text())
            z = man['input']['zrna']
            self.assertEqual(z['mode'], 'auto')
            zid = z['hypothesis']
            self.assertEqual(man['parameters']['trajectory_allocation'][zid], 2)
            zc = [c for c in man['conformers'] if c['secondary_structure_id'] == zid]
            self.assertTrue(zc)
            self.assertTrue(zc[0]['zrna_regions'])
            self.assertTrue(all(c['zrna_regions'] == [] for c in man['conformers'] if c['secondary_structure_id'] != zid))
            # with a Z-binding protein the estimate favours Z (dG < 0): the left-handed model is the first one,
            # although the 3D score, calibrated on A-form helices, scores its stem far worse
            self.assertLess(z['candidates'][0]['dG_kcal_mol'], 0)
            first = man['conformers'][0]
            self.assertEqual(first['secondary_structure_id'], zid)
            self.assertIn('helix_form_score_excluded', first)
            self.assertGreater(first['z_fraction_estimate'], 0.5)
            a_form = [c for c in man['conformers'] if c['secondary_structure_id'] != zid]
            self.assertTrue(a_form and all('helix_form_score_excluded' not in c for c in a_form))
        finally:
            shutil.rmtree(tmp)

    def test_default_is_off(self):
        from rnafold3d.predict import PredictConfig
        self.assertEqual(PredictConfig().zrna, 'off')


class TestDuplexFiles(unittest.TestCase):
    """Z-RNA in the PDB is mostly two-stranded: the viewer shows strands together with the chain option A+B."""

    def two_chain_file(self, tmp):
        out = []
        for line in FIXTURE.read_text().splitlines():
            if line.startswith('ATOM') and int(line[22:26]) > 6:
                line = line[:21] + 'B' + line[22:]
            out.append(line)
        f = Path(tmp) / 'duplex.pdb'
        f.write_text('\n'.join(out) + '\n')
        return f

    def test_join_chains_and_pairs_across_the_break(self):
        from rcv.data import load_document
        tmp = tempfile.mkdtemp()
        try:
            f = self.two_chain_file(tmp)
            one = load_document([f])
            self.assertEqual(one['sequence'], 'CGCGCG')                  # default: the longest chain only
            both = load_document([f], chain_id='A+B')
            self.assertEqual(both['sequence'], 'CGCGCGCGCGCG')
            pairs = [tuple(p[:2]) for p in both['conformations'][0]['pairs']]
            self.assertIn((5, 6), pairs)                                  # across the strand break, |i-j| = 1
            self.assertEqual(len(pairs), 6)
            dg = both['conformations'][0]['values']['dg2d']
            self.assertLess(dg, -5)                                       # co-folding energy, not "infeasible"
        finally:
            shutil.rmtree(tmp)

    def test_hairpins_still_need_a_loop(self):
        from rnafold3d.annotate import chain_breaks
        seq, res = load_fixture()
        from types import SimpleNamespace as NS
        self.assertEqual(chain_breaks([NS(index=k + 1, base=b, atoms=a) for k, (b, a) in enumerate(zip(seq, res))]), {6})


@unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
class TestViewerDetectorParity(unittest.TestCase):
    """zrna.js (live frames) and zrna.py (manifest) must give the same regions."""

    def test_same_regions(self):
        seq, res = load_fixture()
        rng = np.random.default_rng(1)
        cases = []
        for k, sd in enumerate((0.0, 0.3, 0.3, 0.5, 0.5, 0.8)):
            noisy = [{a: (v + rng.normal(0, sd, 3)).tolist() for a, v in r.items()} for r in res]
            cases.append({'name': f'crystal{k}', 'seq': seq, 'residues': noisy, 'pairs': DUPLEX})
        a_res, _ = built('CGCGCGCGCGCG', DUPLEX, 'A')
        cases.append({'name': 'aform', 'seq': 'CGCGCGCGCGCG', 'residues': [{a: v.tolist() for a, v in r.items()} for r in a_res], 'pairs': DUPLEX})
        tmp = tempfile.mkdtemp()
        try:
            f = Path(tmp) / 'cases.json'
            f.write_text(json.dumps(cases))
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'zrna_unit.mjs'), str(f)], capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            js = json.loads(r.stdout)
        finally:
            shutil.rmtree(tmp)
        found = 0
        for c in cases:
            py = Z.detect(c['seq'], [{a: np.array(v) for a, v in x.items()} for x in c['residues']], c['pairs'])
            self.assertEqual([[list(p) for p in z['pairs']] for z in py], [z['pairs'] for z in js[c['name']]], c['name'])
            for a, b in zip(py, js[c['name']]):
                self.assertTrue(np.allclose(a['twists'], b['twists'], atol=0.11))
                self.assertEqual(a['syn_purines'], b['synPurines'])
            found += bool(py)
        self.assertGreaterEqual(found, 3)
        self.assertEqual(js['aform'], [])


if __name__ == '__main__':
    unittest.main()
