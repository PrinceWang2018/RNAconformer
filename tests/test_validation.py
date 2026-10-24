import unittest
import numpy as np
from rnafold3d.core import parse
from rnafold3d.assembly import assemble
from rnafold3d.validation import aligned_rmsd, validate_candidate, filter_candidates, ValidationConfig

class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.model=parse('A'*9,'(((...)))')
        self.candidates=assemble(self.model,3,seed=11,proposals=10)
    def test_rmsd_rigid_invariance(self):
        a={0:np.array([0.,0.,0.]),1:np.array([1.,0.,0.]),2:np.array([0.,1.,0.])}
        b={k:np.array([-v[1],v[0],v[2]])+np.array([10.,2.,1.]) for k,v in a.items()}
        self.assertLess(aligned_rmsd(a,b),1e-8)
    def test_rmsd_mismatch(self):
        with self.assertRaises(ValueError):aligned_rmsd({0:np.zeros(3)},{1:np.zeros(3)})
    def test_validation_report(self):
        report,coords,closures=validate_candidate(self.model,self.candidates[0])
        self.assertEqual(report.candidate_index,0)
        self.assertEqual(report.n_closures,len(closures))
        self.assertTrue(set(coords).issubset(set(range(9))))
    def test_filter_deterministic(self):
        a=filter_candidates(self.model,self.candidates)
        b=filter_candidates(self.model,self.candidates)
        self.assertEqual(a[1],b[1])
    def test_strict_missing_tails(self):
        m=parse('A'*11,'.(((...))).')
        a=assemble(m,1,seed=1,proposals=1)[0]
        report,_,_=validate_candidate(m,a)
        self.assertIn('missing_residues',report.reasons)
    def test_duplicate_suppression(self):
        cfg=ValidationConfig(rmsd_threshold=1e6,require_complete=False,clash_distance=0.01,closure_tolerance=100.)
        selected,reports,_=filter_candidates(self.model,self.candidates,cfg)
        self.assertLessEqual(len(selected),1)
    def test_invalid_settings(self):
        with self.assertRaises(ValueError):filter_candidates(self.model,self.candidates,max_candidates=0)
        with self.assertRaises(ValueError):validate_candidate(self.model,self.candidates[0],ValidationConfig(clash_distance=-1))

if __name__=='__main__':unittest.main()
