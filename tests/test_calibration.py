import unittest
import numpy as np
from rnafold3d.geometry import NucleotideSites
from rnafold3d.calibration import regularize_loop_frames,validate_guide_geometry

class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.sites=tuple(NucleotideSites(i,(float(i)*6,2,0),(float(i)*6,0,0),(float(i)*6,0,3)) for i in range(5))
    def test_preserves_stem_and_c4(self):
        out=regularize_loop_frames(self.sites,{0,4})
        self.assertEqual(out[0],self.sites[0]);self.assertEqual(out[4],self.sites[4])
        self.assertEqual([s.C4 for s in out],[s.C4 for s in self.sites])
    def test_normalized_loop_offsets(self):
        out=regularize_loop_frames(self.sites,{0,4})
        for s in out[1:4]:
            self.assertAlmostEqual(np.linalg.norm(np.subtract(s.P,s.C4)),3)
            self.assertAlmostEqual(np.linalg.norm(np.subtract(s.N,s.C4)),3.6)
    def test_valid(self):
        r=validate_guide_geometry(regularize_loop_frames(self.sites,{0,4}),{0,4})
        self.assertTrue(r.passed_guide_checks)
        self.assertEqual(r.n_interfaces,2)
    def test_broken_interface(self):
        broken=list(self.sites)
        broken[1]=NucleotideSites(1,(50,2,0),(50,0,0),(50,0,3))
        r=validate_guide_geometry(broken,{0,4})
        self.assertGreater(r.n_bad_interfaces,0)
    def test_duplicate(self):
        with self.assertRaises(ValueError):validate_guide_geometry(self.sites+(self.sites[0],))
    def test_nonfinite(self):
        bad=list(self.sites);bad[2]=NucleotideSites(2,(float('nan'),0,0),bad[2].C4,bad[2].N)
        with self.assertRaises(ValueError):validate_guide_geometry(bad)
if __name__=='__main__':unittest.main()
