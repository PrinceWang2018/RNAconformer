import unittest
import numpy as np
from rnafold3d.core import parse
from rnafold3d.assembly import assemble
from rnafold3d.closure import close_chain,close_assembly

class AssemblyClosureTests(unittest.TestCase):
    def test_deterministic_assembly(self):
        db='((..((...))..((...))..))'
        m=parse('A'*len(db),db)
        a=assemble(m,n_candidates=3,seed=42,proposals=20)
        b=assemble(m,n_candidates=3,seed=42,proposals=20)
        self.assertEqual(len(a),3)
        self.assertEqual(a[0].translations,b[0].translations)
        self.assertTrue(np.isfinite([x.placement_score for x in a]).all())
        self.assertEqual(len(a[0].stem_sites),3)
    def test_pseudoknot_rejected(self):
        m=parse('A'*8,'(([[))]]')
        with self.assertRaises(ValueError):assemble(m)
    def test_loop_closed(self):
        r=close_chain([0,0,0],[8,0,0],(1,2,3),step=6.)
        self.assertTrue(r.converged,r)
        self.assertEqual(len(r.xyz),3)
    def test_unreachable(self):
        r=close_chain([0,0,0],[30,0,0],(1,2),step=6.)
        self.assertEqual(r.reason,'unreachable')
        self.assertFalse(r.converged)
    def test_zero_residue(self):
        self.assertTrue(close_chain([0,0,0],[6,0,0],(),step=6.).converged)
        self.assertFalse(close_chain([0,0,0],[10,0,0],(),step=6.).converged)
    def test_hairpin_closure_reports(self):
        m=parse('G'*9,'(((...)))')
        a=assemble(m,n_candidates=1,seed=1,proposals=2)[0]
        c=close_assembly(m,a)
        self.assertEqual(len(c),1)
        self.assertEqual(c[0].residues,(3,4,5))
    def test_invalid(self):
        with self.assertRaises(ValueError):close_chain([0,0],[1,0,0],(1,))
        with self.assertRaises(ValueError):assemble(parse('A'*9,'(((...)))'),n_candidates=0)

if __name__=='__main__':unittest.main()
