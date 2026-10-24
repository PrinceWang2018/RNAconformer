import unittest
import numpy as np
from rnafold3d import parse,build_topology,aform_three_site,rigid_transform_sites,stem_endpoint_frames

class M1Tests(unittest.TestCase):
    def model(self,db):return parse('A'*len(db),db)
    def test_hairpin(self):
        t=build_topology(self.model('(((...)))'))
        self.assertEqual(t.junctions[1].kind,'hairpin')
        self.assertEqual(t.junctions[1].segments,((3,4,5),))
    def test_internal(self):
        t=build_topology(self.model('((..((...))..))'))
        self.assertEqual(t.junctions[1].kind,'internal_loop')
        self.assertEqual(tuple(map(len,t.junctions[1].segments)),(2,2))
    def test_bulge(self):
        t=build_topology(self.model('((.((...)) ))'.replace(' ','')))
        self.assertEqual(t.junctions[1].kind,'bulge')
    def test_multibranch(self):
        t=build_topology(self.model('((..((...))..((...))..))'))
        self.assertEqual(t.junctions[1].kind,'multibranch_junction')
        self.assertEqual(len(t.junctions[1].child_stems),2)
        self.assertEqual(len(t.junctions[1].segments),3)
    def test_external(self):
        t=build_topology(self.model('..((...))...((...))..'))
        self.assertEqual(len(t.junctions[0].child_stems),2)
        self.assertEqual(tuple(map(len,t.junctions[0].segments)),(2,3,2))
    def test_pseudoknot_flag(self):
        t=build_topology(self.model('(([[))]]'))
        self.assertFalse(t.supported_planar)
        with self.assertRaises(ValueError):build_topology(self.model('(([[))]]'),strict=True)
    def test_three_site(self):
        m=self.model('(((...)))');s=aform_three_site(m.stems[0]);self.assertEqual(len(s),6)
        self.assertEqual(set(x.residue for x in s),{0,1,2,6,7,8})
        for x in s:
            self.assertTrue(np.isfinite(np.array([x.P,x.C4,x.N])).all())
        shifted=rigid_transform_sites(s,np.eye(3),[1,2,3]);np.testing.assert_allclose(np.array(shifted[0].C4)-np.array(s[0].C4),[1,2,3])
        self.assertEqual(len(stem_endpoint_frames(s,m.stems[0])),2)
if __name__=='__main__':unittest.main()
