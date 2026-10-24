import unittest
from rnafold3d.atomistic import Atom
from rnafold3d.stereochemistry import validate_stereochemistry,dihedral

class TestStereochemistry(unittest.TestCase):
    def test_missing_atoms_fails(self):
        report=validate_stereochemistry([{}], 'A')
        self.assertFalse(report['passed'])
        self.assertTrue(any(x['kind']=='missing_atom' for x in report['issues']))
    def test_broken_link_fails(self):
        r={"O3'":Atom("O3'",'O',(0,0,0)), 'P':Atom('P','P',(0,0,0))}
        report=validate_stereochemistry([r,r], 'AA')
        self.assertTrue(any(x['kind']=='phosphodiester_distance' for x in report['issues']))
    def test_nonfinite_fails(self):
        report=validate_stereochemistry([{'P':Atom('P','P',(float('nan'),0,0))}], 'A')
        self.assertTrue(any(x['kind']=='nonfinite' for x in report['issues']))
    def test_torsion(self):
        def a(x):return Atom('X','C',x)
        v=dihedral(a((0,1,0)),a((0,0,0)),a((1,0,0)),a((1,0,1)))
        self.assertAlmostEqual(abs(v),90,places=5)
    def test_sequence_length(self):
        with self.assertRaises(ValueError):validate_stereochemistry([{}],'AA')
if __name__=='__main__':unittest.main()
