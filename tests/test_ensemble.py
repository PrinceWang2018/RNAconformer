import unittest
from types import SimpleNamespace
from rnafold3d.core import parse
from rnafold3d.geometry import NucleotideSites
from rnafold3d.ensemble import select_ensemble, EnsembleConfig

class EnsembleTests(unittest.TestCase):
    def setUp(self):
        self.model=parse('AAAA','....')
    def state(self, energy, shift=0., incomplete=False):
        n=3 if incomplete else 4
        sites=tuple(NucleotideSites(i,(i*6.+shift,1.,0.),(i*6.+shift,0.,0.),(i*6.+shift,0.,1.)) for i in range(n))
        return SimpleNamespace(energy=energy,breakdown={'total':energy},sites=sites)
    def test_duplicate_and_best(self):
        out=select_ensemble(self.model,[{'source':'a','state':self.state(4)}, {'source':'b','state':self.state(2,10)}])
        self.assertEqual(out['selected_count'],1)
        self.assertEqual(out['selected_ids'],[1])
    def test_incomplete_rejected(self):
        out=select_ensemble(self.model,[{'source':'bad','state':self.state(1,incomplete=True)}])
        self.assertEqual(out['selected_count'],0)
        self.assertIn('incomplete_chain',out['records'][0]['reasons'])
    def test_invalid_config(self):
        with self.assertRaises(ValueError):EnsembleConfig(rmsd_threshold=0)
    def test_empty(self):
        self.assertEqual(select_ensemble(self.model,[])['selected_count'],0)
    def test_different_shapes(self):
        a=self.state(3)
        sites=list(a.sites)
        sites[-1]=NucleotideSites(3,(18,7,0),(18,6,0),(18,6,1))
        b=SimpleNamespace(energy=4,breakdown={},sites=tuple(sites))
        out=select_ensemble(self.model,[{'source':'a','state':a},{'source':'b','state':b}],EnsembleConfig(rmsd_threshold=0.5,max_link_error=10))
        self.assertEqual(out['selected_count'],2)
if __name__=='__main__':unittest.main()
