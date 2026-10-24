import tempfile,unittest
from pathlib import Path
from rnafold3d.atomistic import Atom,fit_residue,reconstruct,validate_atoms,write_pdb,export_ensemble
from rnafold3d.geometry import NucleotideSites

class AtomisticTests(unittest.TestCase):
    def setUp(self):
        self.atoms={k:Atom(k,'P' if k=='P' else 'C',(float(i%3),float(i//3),float((i*i)%5))) for i,k in enumerate(['P',"C4'",'N9',"O3'","O5'","C3'","C5'","C1'",'C2','C4','C5','C6','N1','N3','N6','N7','C8','O2'])}
        self.site=NucleotideSites(0,self.atoms['P'].xyz,self.atoms["C4'"].xyz,self.atoms['N9'].xyz)
    def test_fit(self):
        _,err=fit_residue(self.atoms,self.site,'A');self.assertLess(err,1e-8)
    def test_reconstruct(self):
        result,errors=reconstruct('A',[self.site],{'A':self.atoms});self.assertEqual(len(result),1);self.assertLess(errors[0],1e-8)
    def test_missing_template(self):
        with self.assertRaises(ValueError):reconstruct('G',[self.site],{})
    def test_incomplete_chain(self):
        with self.assertRaises(ValueError):reconstruct('AA',[self.site],{'A':self.atoms})
    def test_bad_anchor(self):
        bad=NucleotideSites(0,(0,0,0),(0,0,0),(0,0,0))
        with self.assertRaises(ValueError):reconstruct('A',[bad],{'A':self.atoms})
    def test_validation_single(self):
        rep=validate_atoms([self.atoms]);self.assertTrue(rep['passed_basic_checks'])
    def test_pdb(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'a.pdb';write_pdb(path,'A',[self.atoms],validate_atoms([self.atoms]));self.assertIn('ATOM',path.read_text())
    def test_export(self):
        with tempfile.TemporaryDirectory() as d:
            t=Path(d)/'a.pdb';write_pdb(t,'A',[self.atoms],validate_atoms([self.atoms]))
            result=export_ensemble('A',[[self.site]],{'A':t},Path(d)/'out')
            self.assertEqual(len(result['conformers']),1)
            self.assertTrue((Path(d)/'out'/'manifest.json').exists())
if __name__=='__main__':unittest.main()
