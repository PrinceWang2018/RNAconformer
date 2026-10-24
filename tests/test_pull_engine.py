"""Engine support for RCV pull exploration: distance restraints used as a temporary bias, then removed."""
import tempfile
import unittest

import numpy as np

from rnafold3d.fold import FoldModel
from rnafold3d.predict import PredictConfig, predict
from rnafold3d.search import Scorer

SEQ = 'GGCACUUCGGUGCC'
PAIRS = [(0, 13), (1, 12), (2, 11), (3, 10), (4, 9)]


class TestPullRestraints(unittest.TestCase):
    def _model(self):
        model = FoldModel(SEQ, PAIRS, explicit_loops=True)
        model.set_loops_active(True)
        model.sample_torsions(np.random.default_rng(3))
        model.build()
        return model

    def test_flat_bottom_energy_and_bias_switch(self):
        model = self._model()
        I = model.idx
        d = float(np.linalg.norm(model.xyz[I[(6, "C4'")]] - model.xyz[I[(0, "C4'")]]))
        at_rest = Scorer(model, pull_restraints=[{'i': 6, 'j': 0, 'target': d, 'tol': 1.0}])
        self.assertAlmostEqual(at_rest.pull_energy(), 0.0, places=6)
        self.assertAlmostEqual(at_rest.pull_deviation(), 0.0, places=6)
        pulled = Scorer(model, pull_restraints=[{'i': 6, 'j': 0, 'target': d + 5.0, 'tol': 1.0}])
        self.assertAlmostEqual(pulled.pull_energy(), 16.0, places=4)          # (5 - 1)^2
        biased = pulled(detail=True)
        self.assertAlmostEqual(biased['pull_bias'], pulled.w.pull * 16.0, places=3)
        pulled.pull_scale = 0.0                                                # bias removed for relaxation
        unbiased = pulled(detail=True)
        self.assertEqual(unbiased['pull_bias'], 0.0)
        self.assertAlmostEqual(biased['total'] - unbiased['total'], biased['pull_bias'], places=3)

    def test_pull_guided_prediction_reports_unbiased_scores(self):
        cfg = PredictConfig(n_models=2, trajectories=2, seed=5, workers=2, relax_steps=400,
                            seed_structures=('(((((....)))))',),
                            pull_restraints=({'i': 6, 'j': 0, 'a': "C4'", 'b': "C4'", 'target': 24.0, 'tol': 1.5},))
        with tempfile.TemporaryDirectory() as tmp:
            m = predict(SEQ, tmp, name='pull', config=cfg, log=None)
            self.assertTrue(m['conformers'])
            hyp = m['input']['secondary_structures'][0]
            self.assertEqual(hyp['source'], 'seed')
            for c in m['conformers']:
                self.assertTrue(c['validation']['passed'])
                self.assertIn('rms_relaxed_A', c['pull'])
                self.assertEqual(c['score_terms']['pull_bias'], 0.0)


if __name__ == '__main__':
    unittest.main()
