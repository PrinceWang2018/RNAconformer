import unittest
import numpy as np
from rnafold3d.core import parse
from rnafold3d.assembly import assemble
from rnafold3d.sampling import evaluate
from rnafold3d.ensemble_sampling import PoolConfig, loop_runs, rotate_loop, run_pool
from rnafold3d.ensemble import EnsembleConfig

class TestM24(unittest.TestCase):
    def setUp(self):
        self.model=parse('GGGAAACCC','(((...)))')
        self.assembly=assemble(self.model,1,42,10)[0]
    def test_loop_runs(self):
        self.assertEqual(loop_runs(self.model),((3,4,5),))
    def test_loop_rotation_preserves_virtual_bonds_and_stem(self):
        state=evaluate(self.model,self.assembly,seed=42)
        self.assertIsNotNone(state)
        trial=rotate_loop(self.model,state,np.random.default_rng(4),0.8)
        self.assertIsNotNone(trial)
        old={s.residue:np.array(s.C4) for s in state.sites}
        new={s.residue:np.array(s.C4) for s in trial.sites}
        for i in (0,1,2,6,7,8):np.testing.assert_allclose(old[i],new[i])
        for i in range(2,6):self.assertAlmostEqual(np.linalg.norm(old[i]-old[i+1]),np.linalg.norm(new[i]-new[i+1]),places=7)
    def test_config_validation(self):
        with self.assertRaises(ValueError):PoolConfig(snapshot_interval=0)
        with self.assertRaises(ValueError):PoolConfig(loop_move_probability=1.1)
    def test_pool_anneal_snapshots(self):
        result=run_pool(self.model,[self.assembly],PoolConfig(steps=8,snapshot_interval=2,seed=42),ensemble_config=EnsembleConfig(max_clashes=50))
        self.assertEqual(result['m2_4']['snapshot_count'],6)
        self.assertGreaterEqual(result['selected_count'],1)
        self.assertGreater(result['m2_4']['runs'][0]['attempted']['loop'],0)
    def test_pool_remc(self):
        result=run_pool(self.model,[self.assembly],PoolConfig(steps=4,snapshot_interval=2,mode='remc',exchange_interval=2,seed=42),ensemble_config=EnsembleConfig(max_clashes=50))
        self.assertEqual(result['m2_4']['snapshot_count'],13)
        self.assertEqual(result['m2_4']['runs'][0]['swap_attempts'],3)

if __name__=='__main__':unittest.main()
