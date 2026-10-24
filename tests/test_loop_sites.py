import unittest
import numpy as np
from rnafold3d.core import parse
from rnafold3d.assembly import assemble
from rnafold3d.loop_sites import build_full_chain_guides
from rnafold3d.geometry import NucleotideSites
from rnafold3d.energy import score_sites

class LoopSitesTests(unittest.TestCase):
    def setUp(self):
        self.model=parse('A'*24,'((..((...))..((...))..))')
        self.candidate=assemble(self.model,1,seed=42,proposals=20)[0]
    def test_stem_unchanged(self):
        out=build_full_chain_guides(self.model,self.candidate)
        by_id={s.residue:s for s in out.sites}
        for stem in self.candidate.stem_sites:
            for s in stem:self.assertEqual(by_id[s.residue],s)
    def test_no_duplicate_and_finite(self):
        out=build_full_chain_guides(self.model,self.candidate)
        ids=[s.residue for s in out.sites]
        self.assertEqual(len(ids),len(set(ids)))
        for s in out.sites:self.assertTrue(np.isfinite(np.array([s.P,s.C4,s.N])).all())
    def test_complete_flag_consistency(self):
        out=build_full_chain_guides(self.model,self.candidate)
        self.assertEqual(out.complete, not out.missing_residues and all(c.converged for c in out.closures) and out.max_junction_link_error<=0.05)
    def test_complete_scores_when_available(self):
        out=build_full_chain_guides(self.model,self.candidate)
        if out.complete:self.assertEqual(score_sites(self.model,out.sites).n_residues,len(self.model.sequence))
    def test_terminal_tails_reported(self):
        model=parse('A'*11,'..(((...)))')
        candidate=assemble(model,1,seed=2,proposals=0)[0]
        out=build_full_chain_guides(model,candidate)
        self.assertIn(0,out.missing_residues)
        self.assertIn(1,out.missing_residues)
        self.assertFalse(out.complete)
    def test_reproducible(self):
        a=build_full_chain_guides(self.model,self.candidate,seed=7)
        b=build_full_chain_guides(self.model,self.candidate,seed=7)
        self.assertEqual(a,b)

if __name__=='__main__':unittest.main()
