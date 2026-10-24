import unittest
import numpy as np
from rnafold3d.atomistic import Atom, validate_atoms
from rnafold3d.torsion_refinement import refine_torsions, TorsionConfig, _rotatable_edges

class TorsionTests(unittest.TestCase):
    def residue(self,shift):
        positions={'P':(0.,0.,0.),'OP1':(-.5,1.,0.),'OP2':(-.5,-1.,0.),"O5'":(1.5,0.,0.),"C5'":(2.5,1.,0.),"C4'":(3.5,1.,1.),"O4'":(3.,2.,1.),"C1'":(3.,2.,2.),"C2'":(4.,2.,2.),"C3'":(4.5,1.,2.),"O3'":(5.5,1.,2.),"O2'":(4.,3.,2.),'N9':(2.,2.,3.),'C8':(1.,2.,3.),'N7':(1.,3.,3.),'C5':(2.,4.,3.),'C4':(2.5,3.,3.),'N3':(3.5,3.,3.),'C2':(4.,4.,3.),'N1':(3.5,5.,3.),'C6':(2.5,5.,3.),'N6':(2.,6.,3.)}
        return {k:Atom(k,'P' if k=='P' else 'N' if k.startswith('N') else 'O' if k.startswith('O') else 'C',tuple(np.array(v)+shift)) for k,v in positions.items()}
    def test_has_backbone_rotors(self):
        a=self.residue(np.zeros(3));b=self.residue(np.array([5.,0.,4.]))
        names=[x[0] for x in _rotatable_edges([a,b],'AA')]
        self.assertIn('gamma',names);self.assertIn('epsilon',names);self.assertIn('zeta',names)
    def test_bond_lengths_preserved_for_rotor(self):
        a=self.residue(np.zeros(3));b=self.residue(np.array([5.,0.,4.]))
        out,report=refine_torsions([a,b],'AA',TorsionConfig(max_nfev=12))
        self.assertGreater(report['torsions'],0)
        for i in range(2):
            for x,y in [('P',"O5'"),("O5'","C5'"),("C5'","C4'"),("C3'","O3'")]:
                d0=np.linalg.norm(np.array([a,b][i][x].xyz)-[a,b][i][y].xyz)
                d1=np.linalg.norm(np.array(out[i][x].xyz)-out[i][y].xyz)
                self.assertAlmostEqual(d0,d1,places=6)
        self.assertLessEqual(report['objective_after'],report['objective_before']+1e-7)
    def test_rejects_invalid(self):
        with self.assertRaises(ValueError):refine_torsions([self.residue(np.zeros(3))],'GG')
        with self.assertRaises(ValueError):refine_torsions([{'P':Atom('P','P',(0,0,0))}],'A')
    def test_single_residue(self):
        a=self.residue(np.zeros(3))
        out,report=refine_torsions([a],'A',TorsionConfig(max_nfev=8))
        self.assertEqual(len(out),1)
        self.assertTrue(np.isfinite(report['objective_after']))
if __name__=='__main__':unittest.main()
