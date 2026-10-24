"""Tests for the v2 engine (builder, fold tree, loop closure, validation, metrics, pipeline)."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rnafold3d.geomutil import place, dihedral, angle, kabsch, rotation_about_axis, frame, from_frame, to_frame
from rnafold3d.build import load_geometry, ResidueBuilder, build_stem, build_strand, atom_order
from rnafold3d.fold import FoldModel, check_pairs
from rnafold3d.loops import TorsionPrior, close_linker, ccd_close
from rnafold3d.checks import validate_model
from rnafold3d.annotate import dotbracket_to_pairs, nested_subset, pairs_to_dotbracket, canonical_pairs
from rnafold3d.structure_io import parse_cif, write_pdb, write_mmcif, read_mmcif, read_pdb
from rnafold3d.metrics import kabsch_rmsd, ermsd, c4_rmsd
from rnafold3d.structure_io import Residue

TRNA = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
TRNA_DB = '(((((((..((((........)))).(((((.......))))).....(((((.......))))))))))))....'


def pairs0(db):
    return [(i - 1, j - 1) for i, j in dotbracket_to_pairs(db)]


class GeometryPrimitives(unittest.TestCase):
    def test_place_reproduces_internal_coordinates(self):
        a, b, c = np.array([0, 0, 0.]), np.array([1.5, 0, 0]), np.array([2.1, 1.4, 0.3])
        for t in (-179.0, -60.0, 0.0, 75.0, 180.0):
            d = place(a, b, c, 1.53, 111.0, t)
            self.assertAlmostEqual(np.linalg.norm(d - c), 1.53, places=9)
            self.assertAlmostEqual(angle(b, c, d), 111.0, places=6)
            self.assertAlmostEqual(((dihedral(a, b, c, d) - t + 180) % 360) - 180, 0.0, places=6)

    def test_kabsch_and_frames(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(7, 3))
        y = x @ rotation_about_axis([1, -2, 0.5], 0.9).T + [3, 1, -2]
        r, t = kabsch(x, y)
        self.assertLess(np.abs(x @ r.T + t - y).max(), 1e-10)
        loc = to_frame(x, x[0], x[1], x[2])
        back = from_frame(loc, y[0], y[1], y[2])
        self.assertLess(np.abs(back - y).max(), 1e-10)
        self.assertAlmostEqual(abs(np.linalg.det(frame(x[0], x[1], x[2]))), 1.0, places=10)


class Builder(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rb = ResidueBuilder(load_geometry())

    def test_strand_links_are_exact(self):
        res = build_strand(self.rb, list('GACU' * 3))
        g = self.rb.g
        for a, b in zip(res, res[1:]):
            self.assertAlmostEqual(np.linalg.norm(a["O3'"] - b['P']), g.bond("O3'-P"), places=6)
            self.assertAlmostEqual(angle(a["C3'"], a["O3'"], b['P']), g.angle("C3'-O3'-P"), places=4)

    def test_ideal_duplex_is_valid_and_self_consistent(self):
        seq = 'GGCAUCGG' + 'AAAA' + 'CCGAUGCC'
        pairs = [(k, 19 - k) for k in range(8)]
        model, rep = build_stem(self.rb, seq, pairs)
        self.assertLess(rep['pair_fit_rmsd'], 0.3)
        # all pairs recovered by the independent annotator
        residues = [Residue(i, seq[i], seq[i], model[i]) for i in sorted(model)]
        found = {(i, j) for i, j, _ in canonical_pairs(residues)}
        self.assertEqual(found, set(pairs))

    def test_geometry_provenance_excludes_benchmark(self):
        raw = load_geometry().raw
        excluded = set(raw['provenance']['filters']['excluded_benchmark_ids'])
        used = {e['pdb_id'] for e in raw['provenance']['entries']}
        self.assertFalse(used & excluded)
        self.assertGreater(len(used), 50)


class FoldTree(unittest.TestCase):
    def test_trna_topology_and_kinematics(self):
        m = FoldModel(TRNA, pairs0(TRNA_DB))
        kinds = m.summary()['linkers']
        self.assertEqual(kinds, {'tree': 3, 'hairpin': 3, 'closure': 1, 'tail3': 1})
        rng = np.random.default_rng(1)
        m.sample_torsions(rng)
        m.set_loops_active(False)
        m.build()
        for k in m.active_dofs:
            before = m.measured_torsion(k)
            m.rotate(k, 31.0)
            after = m.measured_torsion(k)
            self.assertAlmostEqual(((after - before - 31.0 + 180) % 360) - 180, 0.0, places=6)
        for r in range(m.n - 1):
            a, b = m.idx[(r, "O3'")], m.idx[(r + 1, 'P')]
            if m.present[a] and m.present[b]:
                self.assertAlmostEqual(np.linalg.norm(m.xyz[a] - m.xyz[b]), m.rb.b_O3P, places=6)

    def test_rejects_pseudoknots_and_short_hairpins(self):
        with self.assertRaises(ValueError):
            check_pairs(20, [(0, 10), (5, 15)])
        with self.assertRaises(ValueError):
            check_pairs(10, [(2, 4)])

    def test_hairpin_closure_by_ccd(self):
        seq, db = 'GGGCAUUCGAUGCCC'[:14] + 'C', '((((.......))))'
        m = FoldModel(seq, pairs0(db))
        rng = np.random.default_rng(3)
        m.set_loops_active(False)
        m.sample_torsions(rng)
        m.build()
        prior = TorsionPrior(m.geom)
        rep = close_linker(m, m.loop_linkers[0], rng, prior)
        self.assertTrue(rep['closed'])
        self.assertLess(rep['gap_rmsd'], 0.1)
        self.assertEqual(int(m.present.sum()), m.n_atoms - 3)    # all atoms except 5'-terminal phosphate


class Validation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rb = ResidueBuilder(load_geometry())
        cls.seq = 'GGCAUCGGCC'
        cls.res = build_strand(cls.rb, list(cls.seq))

    def test_ideal_strand_passes(self):
        rep = validate_model(self.seq, self.res)
        self.assertTrue(rep['passed'], rep['issues'])

    def test_detects_chain_break_clash_and_mirror_image(self):
        broken = [dict(r) for r in self.res]
        for k in list(broken[5]):
            broken[5][k] = broken[5][k] + np.array([0.0, 0.0, 0.6])
        rep = validate_model(self.seq, broken)
        self.assertIn('chain_break', rep['issue_counts'])
        mirror = [{k: v * np.array([1, 1, -1]) for k, v in r.items()} for r in self.res]
        self.assertIn('chirality', validate_model(self.seq, mirror)['issue_counts'])
        clash = [dict(r) for r in self.res]
        clash[8] = {k: v - clash[8]["C1'"] + clash[2]["C1'"] + 0.5 for k, v in clash[8].items()}
        self.assertGreater(validate_model(self.seq, clash)['severe_clashes'], 0)


class Annotation(unittest.TestCase):
    def test_nested_subset_and_dotbracket(self):
        pairs = [(1, 20), (2, 19), (5, 25), (6, 24), (10, 15)]
        kept, removed = nested_subset(pairs)
        self.assertEqual(len(kept), 3)
        db = pairs_to_dotbracket(30, pairs)
        self.assertEqual(sorted(dotbracket_to_pairs(db)), sorted(pairs))


class Metrics(unittest.TestCase):
    def test_rigid_invariance(self):
        rng = np.random.default_rng(2)
        rb = ResidueBuilder(load_geometry())
        res = build_strand(rb, list('GACUGACUGA'))
        ref = {i + 1: Residue(i + 1, b, b, r) for i, (b, r) in enumerate(zip('GACUGACUGA', res))}
        R = rotation_about_axis([0.3, 1, -1], 2.0)
        moved = {k: Residue(k, v.base, v.comp_id, {a: x @ R.T + 4.0 for a, x in v.atoms.items()}) for k, v in ref.items()}
        self.assertLess(c4_rmsd(ref, moved)[0], 1e-8)
        self.assertLess(ermsd(ref, moved)[0], 1e-8)


class IO(unittest.TestCase):
    def test_cif_parser_quotes_and_loops(self):
        text = """data_X
_struct.title 'A "quoted" title'
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
ATOM "C4'" G
ATOM P U
#
"""
        d = parse_cif(text)
        self.assertEqual(d['_atom_site']['label_atom_id'], ["C4'", 'P'])
        self.assertEqual(d['_struct']['title'], ['A "quoted" title'])

    def test_pdb_and_mmcif_round_trip(self):
        rb = ResidueBuilder(load_geometry())
        seq = 'GCAU'
        res = build_strand(rb, list(seq))
        with tempfile.TemporaryDirectory() as d:
            write_pdb(Path(d) / 'x.pdb', res, seq)
            write_mmcif(Path(d) / 'x.cif', res, seq)
            p = next(iter(read_pdb(Path(d) / 'x.pdb')[0].values()))
            c = next(iter(read_mmcif(Path(d) / 'x.cif')[0].values()))
            self.assertEqual(p.sequence, seq)
            self.assertEqual(c.sequence, seq)
            self.assertLess(np.abs(c.residues[2].atoms["C4'"] - res[2]["C4'"]).max(), 1e-3)


class Pipeline(unittest.TestCase):
    def test_end_to_end_is_valid_and_deterministic(self):
        from rnafold3d.predict import predict, PredictConfig
        seq, db = 'GGACUUCGGUCC', '((((....))))'
        cfg = PredictConfig(n_models=2, trajectories=3, seed=7, anneal_steps=200, workers=1)
        with tempfile.TemporaryDirectory() as d:
            m1 = predict(seq, Path(d) / 'a', dbn=db, config=cfg, log=None)
            m2 = predict(seq, Path(d) / 'b', dbn=db, config=cfg, log=None)
            self.assertGreaterEqual(len(m1['conformers']), 1)
            for c in m1['conformers']:
                self.assertTrue(c['validation']['passed'])
                self.assertTrue((Path(d) / 'a' / c['files']['pdb']).exists())
                self.assertTrue((Path(d) / 'a' / c['files']['mmcif']).exists())
            for sub in ('metadata/manifest.json', 'metadata/scores.csv', 'metadata/clusters.json',
                        'metadata/validation.json', 'config.json', 'logs/run.log', 'input/sequence.fasta'):
                self.assertTrue((Path(d) / 'a' / sub).exists(), sub)
            a = (Path(d) / 'a' / 'conformations' / 'conf_001.pdb').read_text()
            b = (Path(d) / 'b' / 'conformations' / 'conf_001.pdb').read_text()
            self.assertEqual(a, b)
            self.assertEqual(json.loads((Path(d) / 'a' / 'metadata' / 'manifest.json').read_text())['schema'],
                             'rnafold3d.ensemble/1')


if __name__ == '__main__':
    unittest.main()


class EnsemblePackage(unittest.TestCase):
    def test_import_external_conformers(self):
        from rnafold3d.ensemble_io import import_ensemble, load_package
        rb = ResidueBuilder(load_geometry())
        seq = 'GGCAUCG'
        a = build_strand(rb, list(seq))
        R = rotation_about_axis([1, 0, 0], 0.3)
        b = [{k: v @ R.T for k, v in r.items()} for r in a]
        with tempfile.TemporaryDirectory() as d:
            for name, res in (('a.pdb', a), ('b.pdb', b)):
                write_pdb(Path(d) / name, res, seq)
            man = import_ensemble([Path(d) / 'a.pdb', Path(d) / 'b.pdb'], Path(d) / 'pkg', name='t')
            self.assertEqual(len(man['conformers']), 2)
            self.assertTrue(all(c['validation']['passed'] for c in man['conformers']))
            self.assertLess(man['pairwise_c4_rmsd']['matrix'][0][1], 1e-3)   # rigid copy
            self.assertEqual(load_package(Path(d) / 'pkg')['schema'], 'rnafold3d.ensemble/1')


class UserInput(unittest.TestCase):
    def test_noncanonical_and_pseudoknot_pairs_are_reported(self):
        from rnafold3d.predict import predict, PredictConfig
        seq = 'GGGCAAAAGCCCAAGAAAUUCGGG'
        dbn = '(((..[[..)))....]]......'      # contains G-A/A-A style non-canonical pairs
        cfg = PredictConfig(n_models=1, trajectories=1, seed=1, anneal_steps=100, workers=1)
        with tempfile.TemporaryDirectory() as d:
            m = predict(seq, Path(d) / 'r', dbn=dbn, config=cfg, log=None)
            self.assertIn('non-canonical', ' '.join(m['warnings']))
