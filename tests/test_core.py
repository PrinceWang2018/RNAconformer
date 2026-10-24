import unittest
import numpy as np
from rnafold3d import parse,aform_stem,transform,endpoint_reachable

class CoreTests(unittest.TestCase):
    def test_hairpin(self):
        m=parse('GGGAAACCC','(((...)))')
        self.assertEqual(len(m.stems),1)
        self.assertEqual(len(m.pairs),3)
        self.assertEqual(len(m.elements),1)
        self.assertEqual(m.elements[0].kind,'hairpin_or_internal')
    def test_multistem(self):
        m=parse('GGGAAACCCGGGAAACCC','(((...)))(((...)))')
        self.assertEqual(len(m.stems),2)
    def test_pseudoknot_brackets(self):
        m=parse('GGGAAACCC','(((.[[)))') if False else parse('GGGGCCCC','(([[))]]')
        self.assertEqual(len(m.pairs),4)
    def test_invalid(self):
        with self.assertRaises(ValueError): parse('GGG','(((')
        with self.assertRaises(ValueError): parse('GGG','...(')
    def test_geometry(self):
        m=parse('GGGAAACCC','(((...)))'); ids,xyz=aform_stem(m.stems[0])
        self.assertEqual(xyz.shape,(6,3))
        self.assertEqual(set(ids),{0,1,2,6,7,8})
        np.testing.assert_allclose(transform(xyz,np.eye(3),np.zeros(3)),xyz)
    def test_reach(self):
        self.assertTrue(endpoint_reachable([0,0,0],[10,0,0],2))
        self.assertFalse(endpoint_reachable([0,0,0],[20,0,0],2))
if __name__=='__main__':unittest.main()
