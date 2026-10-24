import unittest
import numpy as np
from rnafold3d.atomistic import Atom,validate_atoms
from rnafold3d.refinement import refine_rigid_chain,RefinementConfig

class RefinementTests(unittest.TestCase):
    def residue(self,shift):
        coords={'P':(0.,0.,0.),"O5'":(1.,0.,0.),"C4'":(0.,2.,0.),"C3'":(0.,2.,1.),"O3'":(0.,2.,2.),'N9':(1.,3.,0.)}
        return {k:Atom(k,'P' if k=='P' else ('N' if k=='N9' else 'O' if k.startswith('O') else 'C'),tuple(np.array(v)+shift)) for k,v in coords.items()}
    def test_reduces_broken_link(self):
        a=self.residue(np.array([0.,0.,0.]));b=self.residue(np.array([0.,2.,7.]))
        original=validate_atoms([a,b])['link_distances'][0]
        refined,report=refine_rigid_chain([a,b],RefinementConfig(max_nfev=160))
        new=validate_atoms(refined)['link_distances'][0]
        self.assertLess(abs(new-1.61),abs(original-1.61))
        self.assertEqual(len(refined),2)
    def test_preserves_intra_residue_geometry(self):
        a=self.residue(np.zeros(3));b=self.residue(np.array([0.,2.,7.]))
        refined,_=refine_rigid_chain([a,b])
        for before,after in zip([a,b],refined):
            keys=list(before)
            for i in range(len(keys)):
                for j in range(i+1,len(keys)):
                    u,v=keys[i],keys[j]
                    self.assertAlmostEqual(np.linalg.norm(np.array(before[u].xyz)-before[v].xyz),np.linalg.norm(np.array(after[u].xyz)-after[v].xyz),places=6)
    def test_reject_missing_atoms(self):
        with self.assertRaises(ValueError):refine_rigid_chain([{'P':Atom('P','P',(0,0,0))}])
    def test_guard_length(self):
        a=self.residue(np.zeros(3))
        with self.assertRaises(ValueError):refine_rigid_chain([a]*41)
if __name__=='__main__':unittest.main()
