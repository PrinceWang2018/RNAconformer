import unittest
from dataclasses import replace
import numpy as np
from rnafold3d.core import parse
from rnafold3d.geometry import NucleotideSites, aform_three_site
from rnafold3d.energy import EnergyConfig, score_sites

class EnergyTests(unittest.TestCase):
    def setUp(self):
        self.model=parse('GGGAAACCC','(((...)))')
        self.sites=list(aform_three_site(self.model.stems[0]))

    def test_decomposition_finite(self):
        e=score_sites(self.model,self.sites,require_complete=False)
        self.assertAlmostEqual(e.total,sum((e.bonded,e.pairing,e.orientation,e.stacking,e.steric,e.electrostatic)))
        self.assertEqual(e.n_pairs,3)
        self.assertEqual(e.n_stacks,4)

    def test_pair_displacement_penalty(self):
        a=score_sites(self.model,self.sites,require_complete=False)
        moved=[]
        for s in self.sites:
            if s.residue==8:
                moved.append(replace(s,N=tuple(np.array(s.N)+[100,0,0])))
            else: moved.append(s)
        b=score_sites(self.model,moved,require_complete=False)
        self.assertGreater(b.pairing,a.pairing)

    def test_stacking_direction(self):
        a=score_sites(self.model,self.sites,require_complete=False)
        self.assertLessEqual(a.stacking,0)

    def test_steric_clash(self):
        model=parse('AAAAAAAA','........')
        def make(i,x):return NucleotideSites(i,(x,0,1),(x,0,0),(x+1,0,0))
        far=[make(0,0),make(5,10)]
        close=[make(0,0),make(5,0.5)]
        self.assertGreater(score_sites(model,close,require_complete=False).steric,
                           score_sites(model,far,require_complete=False).steric)

    def test_screened_electrostatics(self):
        model=parse('AAAAAAAA','........')
        def make(i,x):return NucleotideSites(i,(x,0,1),(x,0,0),(x+1,0,0))
        near=score_sites(model,[make(0,0),make(5,5)],require_complete=False)
        far=score_sites(model,[make(0,0),make(5,25)],require_complete=False)
        self.assertGreater(near.electrostatic,far.electrostatic)

    def test_invalid_and_missing(self):
        with self.assertRaises(ValueError):score_sites(self.model,self.sites)
        with self.assertRaises(ValueError):EnergyConfig(screening_length=0)
        with self.assertRaises(ValueError):score_sites(self.model,self.sites+[self.sites[0]],require_complete=False)

    def test_rigid_invariance(self):
        e1=score_sites(self.model,self.sites,require_complete=False)
        moved=[replace(s,P=tuple(np.array(s.P)+[7,2,9]),C4=tuple(np.array(s.C4)+[7,2,9]),N=tuple(np.array(s.N)+[7,2,9])) for s in self.sites]
        e2=score_sites(self.model,moved,require_complete=False)
        self.assertAlmostEqual(e1.total,e2.total,places=9)

if __name__=='__main__':unittest.main()
