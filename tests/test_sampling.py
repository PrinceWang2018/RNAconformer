import unittest
import numpy as np
from rnafold3d.core import parse
from rnafold3d.assembly import assemble
from rnafold3d.sampling import SamplingConfig, evaluate, perturb, run_annealing, run_replica_exchange

class TestSampling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=parse('GGGAAACCC','(((...)))')
        cls.assembly=assemble(cls.model,n_candidates=1,seed=42,proposals=10)[0]
    def test_complete_energy(self):
        state=evaluate(self.model,self.assembly,seed=3)
        self.assertIsNotNone(state)
        self.assertEqual(len(state.sites),9)
    def test_rigid_stem_geometry(self):
        moved=perturb(self.assembly,np.random.default_rng(2),SamplingConfig())
        for a,b in zip(self.assembly.stem_sites,moved.stem_sites):
            x=np.asarray([s.C4 for s in a]);y=np.asarray([s.C4 for s in b])
            self.assertTrue(np.allclose(np.linalg.norm(x[:,None]-x[None,:],axis=-1),np.linalg.norm(y[:,None]-y[None,:],axis=-1)))
    def test_anneal_reproducible(self):
        cfg=SamplingConfig(steps=8,seed=7)
        a=run_annealing(self.model,self.assembly,cfg)
        b=run_annealing(self.model,self.assembly,cfg)
        self.assertEqual(a['trace'],b['trace'])
        self.assertLessEqual(a['best'].energy,a['trace'][0]+1e-9)
    def test_remc_accounting(self):
        cfg=SamplingConfig(steps=8,seed=8,exchange_interval=2)
        r=run_replica_exchange(self.model,self.assembly,cfg)
        self.assertEqual(len(r['traces']),4)
        self.assertEqual(len(r['traces'][0]),9)
        self.assertLessEqual(r['swaps'],r['swap_attempts'])
    def test_invalid_settings(self):
        with self.assertRaises(ValueError):SamplingConfig(temperature=-1)
if __name__=='__main__':unittest.main()
