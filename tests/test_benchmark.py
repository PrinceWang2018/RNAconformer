import unittest
import tempfile
from pathlib import Path
import numpy as np
from rnafold3d.benchmark import aligned_rmsd, read_pdb, benchmark


def pdb(path, coords, residues=('A','C','G','U')):
    with open(path,'w') as f:
        for i,(xyz,base) in enumerate(zip(coords,residues),1):
            x,y,z=xyz
            f.write(f"ATOM   {i:4d}  C4' {base:>3s} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n")
        f.write('END\n')

class BenchmarkTests(unittest.TestCase):
    def test_rigid_invariance(self):
        a=np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3.]],float)
        b=a@np.array([[0,-1,0],[1,0,0],[0,0,1]])+np.array([10,4,-3])
        self.assertLess(aligned_rmsd(a,b),1e-8)
    def test_nonzero_error(self):
        a=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1.]],float)
        b=a.copy();b[0]+=np.array([0,0,2])
        self.assertGreater(aligned_rmsd(a,b),0.1)
    def test_reference_and_missing_residue(self):
        with tempfile.TemporaryDirectory() as d:
            ref=Path(d)/'ref.pdb'; good=Path(d)/'good.pdb'; bad=Path(d)/'bad.pdb'
            coords=np.array([[0,0,0],[2,0,0],[2,2,0],[2,2,2.]])
            pdb(ref,coords);pdb(good,coords+3);pdb(bad,coords[:2])
            result=benchmark(ref,[good,bad])
            self.assertEqual(result['n_evaluable'],1)
            self.assertTrue(result['success_at_k'])
            self.assertEqual(result['per_conformer'][1]['status'],'insufficient_coverage')
    def test_missing_reference(self):
        with self.assertRaises(OSError): read_pdb('/not/a/real/file.pdb')
