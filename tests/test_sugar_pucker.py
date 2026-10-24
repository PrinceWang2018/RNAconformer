import unittest
import numpy as np
from rnafold3d.atomistic import Atom
from rnafold3d.sugar_pucker import dihedral,pseudorotation,refine_sugar_pucker,refine_chain_puckers,PuckerConfig


def residue():
    # Non-planar ring with chemically named substituent placeholders.
    names=["O4'","C1'","C2'","C3'","C4'"]
    pts=[(1.2,0.,0.15),(.37,1.14,-.1),(-.97,.71,.13),(-.97,-.71,-.12),(.37,-1.14,.08)]
    d={k:Atom(k,k[0],p) for k,p in zip(names,pts)}
    for k,p in {"O2'":(-1.5,1.2,.3),"O3'":(-1.6,-1.1,-.3),"C5'":(.6,-2.0,.4)}.items():d[k]=Atom(k,k[0],p)
    return d

class SugarPuckerTests(unittest.TestCase):
    def test_dihedral(self):
        self.assertTrue(np.isfinite(dihedral([1,0,0],[0,0,0],[0,1,0],[0,1,1])))
    def test_pseudorotation(self):
        p=pseudorotation(residue());self.assertTrue(0<=p['phase_deg']<360)
    def test_refine(self):
        r=residue();out,report=refine_sugar_pucker(r,PuckerConfig(max_nfev=30))
        self.assertLessEqual(report['objective_after'],report['objective_before']+1e-7)
        self.assertEqual(out["O2'"],r["O2'"])
    def test_chain(self):
        out,reports=refine_chain_puckers([residue()],PuckerConfig(max_nfev=15))
        self.assertEqual(len(out),1);self.assertEqual(len(reports),1)
    def test_incomplete(self):
        r=residue();del r["C3'"]
        with self.assertRaises(ValueError):refine_sugar_pucker(r)

if __name__=='__main__':unittest.main()
