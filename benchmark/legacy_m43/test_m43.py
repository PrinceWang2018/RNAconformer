import tempfile,unittest,json
from pathlib import Path
from run_m43 import prepare_reference

CIF='''data_test
loop_
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.label_comp_id
_atom_site.label_atom_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.label_alt_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
R 1 G "C4'" 0 0 0 1 . X 11 ?
R 2 C "C4'" 1 0 0 1 . X 12 ?
R 3 A "C4'" 2 0 0 1 . X 13 ?
'''
class TestReference(unittest.TestCase):
 def test_reference_mapping(self):
  with tempfile.TemporaryDirectory() as td:
   root=Path(td);(root/'ref.cif').write_text(CIF)
   rep=prepare_reference(root/'ref.cif','R',root/'out')
   self.assertEqual(rep['observed_canonical_residues'],3)
   self.assertEqual(rep['c4_coverage'],1)
   self.assertEqual(json.loads((root/'out'/'residue_mapping.json').read_text())['1']['auth_seq_id'],'11')
   self.assertIn("C4'",(root/'out'/'reference.pdb').read_text())
if __name__=='__main__':unittest.main()
