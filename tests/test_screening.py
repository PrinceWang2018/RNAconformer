"""0.7.5: the compound-screening core — compounds (`chem`), the RNA as a receptor (`receptor`), the scoring
function (`score`), its grid (`grid`) and the search (`dock`).

No network and no downloads: the compounds here are built in the test, the receptor is a sample structure
from the tree. The benchmark that redocks real complexes is `benchmark/screening/redock.py`, which needs
the structures from RCSB and is not a unit test.
"""
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d.screening import chem, dock, grid, score          # noqa: E402
from rnafold3d.screening.receptor import Receptor, read_complex  # noqa: E402

SAMPLES = ROOT / 'benchmark' / 'data' / 'samples'

# a glycine zwitterion drawn by hand (heavy atoms, with the hydrogens left out): N-CH2-COOH
GLYCINE = """glycine
  test

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0
    1.4700    0.0000    0.0000 C   0  0
    2.0000    1.4000    0.0000 C   0  0
    1.3000    2.4000    0.0000 O   0  0
    3.3000    1.5000    0.0000 O   0  0
  1  2  1  0
  2  3  1  0
  3  4  2  0
  3  5  1  0
M  END
$$$$
"""
# benzamide: a planar aromatic ring with a rotatable bond to an amide
BENZAMIDE = """benzamide
  test

  9  9  0  0  0  0  0  0  0  0999 V2000
    0.0000    1.3950    0.0000 C   0  0
    1.2080    0.6975    0.0000 C   0  0
    1.2080   -0.6975    0.0000 C   0  0
    0.0000   -1.3950    0.0000 C   0  0
   -1.2080   -0.6975    0.0000 C   0  0
   -1.2080    0.6975    0.0000 C   0  0
    0.0000    2.8950    0.0000 C   0  0
    1.2000    3.5000    0.0000 O   0  0
   -1.2000    3.5000    0.0000 N   0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  1  7  1  0
  7  8  2  0
  7  9  1  0
M  END
$$$$
"""


class TestCompounds(unittest.TestCase):
    def test_a_molfile_is_read_with_its_rings_and_flags(self):
        m = chem.read_sdf(BENZAMIDE)[0]
        self.assertEqual(len(m), 9)
        self.assertEqual(m.name, 'benzamide')
        self.assertEqual(m.formula.startswith('C7'), True)
        self.assertEqual(len(m.rings), 1)
        self.assertEqual(int(m.aromatic.sum()), 6)             # the ring, not the amide
        self.assertTrue(m.aromatic[:6].all() and not m.aromatic[6:].any())
        self.assertTrue(m.hydrophobic[1:6].all())              # ring carbons away from the amide
        self.assertFalse(m.hydrophobic[6])                     # the carbonyl carbon is not apolar
        self.assertTrue(m.acceptor[7])                         # the carbonyl oxygen
        self.assertEqual(m.rotatable, [(0, 6)])                # only the ring-to-amide bond turns
        self.assertAlmostEqual(float(m.charge.sum()), 0.0, places=6)

    def test_an_acid_and_a_base_are_charged_at_ph_7(self):
        m = chem.read_sdf(GLYCINE)[0]
        self.assertAlmostEqual(float(m.charge.sum()), 0.0, places=6)      # +1 amine, −1 carboxylate
        self.assertGreater(m.charge[0], 0.9)                              # the nitrogen
        self.assertLess(m.charge[3] + m.charge[4], -0.9)                  # the two oxygens share −1
        self.assertAlmostEqual(m.charge[3], m.charge[4], places=6)
        neutral = chem.prepare(m.elements, m.xyz, m.bonds, ph=None)
        self.assertAlmostEqual(float(neutral.charge.sum()), 0.0, places=6)
        self.assertLess(abs(neutral.charge[0]), 0.5)                      # no perception: no formal charge

    def test_only_an_aromatic_nitrogen_with_a_ring_double_bond_is_a_cation(self):
        """A pyridinium is +1; an N-methylated ring amide (the theophylline and flavin case) is neutral.

        Both are six-membered, planar, aromatic by the rule in `prepare`, and both have a nitrogen with
        three heavy bonds. What separates them is whether that nitrogen holds a double bond inside the
        ring, and getting it wrong puts a whole unit of charge on a neutral compound in front of a
        polyanion.
        """
        import numpy as np
        ring = np.array([[np.cos(a), np.sin(a), 0.0] for a in np.arange(6) * np.pi / 3]) * 1.39

        def hexagon(elements, ring_orders, extra=(), hydrogens=None):
            xyz = np.vstack([ring] + [np.asarray(x, float).reshape(1, 3) for _e, x, _b in extra])
            els = list(elements) + [e for e, _x, _b in extra]
            bonds = [(i, (i + 1) % 6, o) for i, o in enumerate(ring_orders)]
            bonds += [(len(ring) + k, b[0], b[1]) for k, (_e, _x, b) in enumerate(extra)]
            return chem.prepare(els, xyz, bonds, hydrogens=hydrogens)

        methyl = ('C', ring[0] * 2.05, (0, 1))                 # on the nitrogen, which is ring atom 0
        pyridinium = hexagon(['N'] + ['C'] * 5, [2, 1, 2, 1, 2, 1], [methyl],
                             hydrogens=[0, 1, 1, 1, 1, 1, 3])
        self.assertTrue(pyridinium.aromatic[:6].all())
        self.assertAlmostEqual(float(pyridinium.charge.sum()), 1.0, places=6)
        self.assertGreater(pyridinium.charge[0], 0.9)

        # the same ring, but the nitrogen's two ring bonds are single and the double bond is the exocyclic
        # carbonyl of the neighbouring carbon: a 1-methyl-2-pyridone, as in a theophylline
        pyridone = hexagon(['N'] + ['C'] * 5, [1, 1, 2, 1, 2, 1],
                           [methyl, ('O', ring[1] * 2.1, (1, 2))], hydrogens=[0, 0, 1, 1, 1, 1, 3, 0])
        self.assertTrue(pyridone.aromatic[:6].all())
        self.assertAlmostEqual(float(pyridone.charge.sum()), 0.0, places=6)
        self.assertLess(abs(pyridone.charge[0]), 0.5)

        # an aniline amine keeps its pair in the ring and is not protonated at pH 7 either
        aniline = hexagon(['C'] * 6, [2, 1, 2, 1, 2, 1], [('N', ring[0] * 2.05, (0, 1))],
                          hydrogens=[0, 1, 1, 1, 1, 1, 2])
        self.assertAlmostEqual(float(aniline.charge.sum()), 0.0, places=6)

    def test_bonds_from_distance_find_the_same_molecule(self):
        m = chem.read_sdf(BENZAMIDE)[0]
        guessed = chem.prepare(m.elements, m.xyz)                          # as a PDB ligand arrives
        self.assertEqual(sorted((i, j) for i, j, _o in guessed.bonds),
                         sorted(tuple(sorted((i, j))) for i, j, _o in m.bonds))
        self.assertEqual(int(guessed.aromatic.sum()), 6)
        self.assertEqual(guessed.rotatable, [(0, 6)])
        # with no hydrogens in the file every polar atom is treated as donor and acceptor
        self.assertTrue(guessed.donor[8] and guessed.acceptor[8])

    def test_torsions_turn_a_branch_and_keep_the_molecule_rigid(self):
        m = chem.read_sdf(BENZAMIDE)[0]
        root, branches = chem.torsion_tree(m)
        self.assertEqual(len(branches), 1)
        self.assertEqual(len(root), 6)                                     # the ring stays still
        self.assertEqual(sorted(branches[0][2]), [6, 7, 8])
        turned = chem.apply_torsions(m, [np.pi / 3], (root, branches))
        self.assertTrue(np.allclose(turned[:6], m.xyz[:6]))                # the root did not move
        self.assertFalse(np.allclose(turned[7], m.xyz[7]))
        for i, j, _o in m.bonds:                                           # bond lengths unchanged
            self.assertAlmostEqual(float(np.linalg.norm(turned[i] - turned[j])),
                                   float(np.linalg.norm(m.xyz[i] - m.xyz[j])), places=6)
        self.assertTrue(np.allclose(chem.apply_torsions(m, [2 * np.pi], (root, branches)), m.xyz, atol=1e-6))

    def test_sdf_round_trip(self):
        m = chem.read_sdf(GLYCINE)[0]
        text = chem.write_sdf(m, properties=[{'score': -1.5}])
        back = chem.read_sdf(text)[0]
        self.assertEqual(back.elements, m.elements)
        self.assertTrue(np.allclose(back.xyz, m.xyz, atol=1e-4))
        self.assertEqual(back.properties['score'], '-1.5')
        self.assertAlmostEqual(float(back.charge.sum()), 0.0, places=6)
        many = chem.write_sdf(m, xyz=[m.xyz, m.xyz + 1.0], properties=[{'pose': 1}, {'pose': 2}])
        self.assertEqual(len(chem.read_sdf(many)), 2)

    def test_a_broken_record_does_not_stop_the_file(self):
        out = chem.read_sdf(GLYCINE + 'not a molfile\n$$$$\n' + BENZAMIDE)
        self.assertEqual(len(out), 3)
        self.assertEqual([len(m) for m in out], [5, 0, 9])
        self.assertIn('error', out[1].properties)


class TestReceptor(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = Receptor.from_pdb(SAMPLES / '1Y26_A' / 'reference.pdb')

    def test_every_nucleotide_carries_the_charge_of_its_phosphate(self):
        rec = self.rec
        charges = [round(float(rec.charge[rec.residue == r].sum()), 6) for r in sorted(set(rec.residue))]
        self.assertEqual(sorted(set(charges)), [-1.0, 0.0])          # 0 only where there is no phosphate
        self.assertEqual(charges.count(0.0), 1)
        self.assertAlmostEqual(float(rec.charge.sum()), -(len(charges) - 1), places=6)

    def test_types_come_from_the_atom_names(self):
        rec = self.rec
        name = {(int(r), n): k for k, (r, n) in enumerate(zip(rec.residue, rec.name))}
        for r in sorted(set(rec.residue))[:40]:
            base = rec.resname[name[(r, "C1'")]] if (r, "C1'") in name else None
            if base == 'G':
                self.assertTrue(rec.donor[name[(r, 'N1')]] and rec.acceptor[name[(r, 'O6')]])
                self.assertTrue(rec.acceptor[name[(r, 'N7')]] and not rec.donor[name[(r, 'N7')]])
            if base == 'C' and (r, 'C5') in name:
                self.assertTrue(rec.hydrophobic[name[(r, 'C5')]])     # the one apolar carbon of a pyrimidine
            if (r, "O2'") in name:
                self.assertTrue(rec.donor[name[(r, "O2'")]] and rec.acceptor[name[(r, "O2'")]])
        self.assertFalse(rec.hydrophobic[[name[(r, 'P')] for r in sorted(set(rec.residue))[1:5]]].any())

    def test_every_base_has_a_plane(self):
        rec = self.rec
        self.assertEqual(len(rec.ring_centre), len(set(rec.residue)))
        self.assertTrue(np.allclose(np.linalg.norm(rec.ring_normal, axis=1), 1.0))
        d = np.linalg.norm(rec.ring_centre[:, None] - rec.ring_centre[None], axis=2)
        self.assertLess(np.median(np.sort(d, axis=1)[:, 1]), 6.0)     # stacked neighbours are close

    def test_a_complex_is_split_into_receptor_and_compounds(self):
        rna = [ln for ln in (SAMPLES / '1Y26_A' / 'reference.pdb').read_text().splitlines()
               if ln.startswith('ATOM')]
        ligand = ['HETATM%5d  C%-2d LIG X 500    %8.3f%8.3f%8.3f  1.00  0.00           C' % (9000 + k, k, x, 0.0, 0.0)
                  for k, x in enumerate((60.0, 61.5, 63.0, 64.5, 66.0, 67.5))]
        water = ['HETATM 9100  O   HOH X 600    %8.3f%8.3f%8.3f  1.00  0.00           O' % (70.0, 0.0, 0.0)]
        rec, ligs = read_complex('\n'.join(rna + ligand + water) + '\n')
        self.assertEqual(len(ligs), 1)
        self.assertEqual(ligs[0].properties['residue'], 'LIG')
        self.assertEqual(len(ligs[0]), 6)
        self.assertEqual(len(rec), len(rna))                           # the water and the ligand are not in it
        far = ['HETATM 9201  C1  LIG X 700    %8.3f%8.3f%8.3f  1.00  0.00           C' % (900.0, 0.0, 0.0)]
        rec2, ligs2 = read_complex('\n'.join(rna + far) + '\n')
        self.assertEqual(ligs2, [])                                    # one atom is not a compound

    def test_linked_residues_are_one_compound(self):
        # an aminoglycoside is deposited as its sugars, one residue each: bonded to each other, not to the
        # chain, so they are one compound and not three
        rna = [ln for ln in (SAMPLES / '1Y26_A' / 'reference.pdb').read_text().splitlines()
               if ln.startswith('ATOM')]
        het = []
        for n, (res, num, x0) in enumerate((('SUG', 500, 60.0), ('SUG', 501, 67.5), ('OTH', 502, 200.0))):
            for k, dx in enumerate((0.0, 1.5, 3.0, 4.5, 6.0, 7.2)):
                het.append('HETATM%5d  C%-2d %s X %3d    %8.3f%8.3f%8.3f  1.00  0.00           C'
                           % (9000 + 10 * n + k, k, res, num, x0 + dx, 0.0, 0.0))
        rec, ligs = read_complex('\n'.join(rna + het) + '\n')
        self.assertEqual(len(ligs), 2)
        self.assertEqual(len(ligs[0]), 12)                             # the two sugars, joined
        self.assertEqual(ligs[0].properties['residues'], ['SUG', 'SUG'])
        self.assertEqual(ligs[1].properties['residues'], ['OTH'])


class TestScoring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = Receptor.from_pdb(SAMPLES / '1Y26_A' / 'reference.pdb')
        cls.mol = chem.read_sdf(BENZAMIDE)[0]

    def test_a_compound_inside_the_rna_is_punished_and_one_far_away_is_nothing(self):
        rec, m = self.rec, self.mol
        far = m.copy_with(m.xyz + rec.xyz.mean(0) + 60.0)
        self.assertAlmostEqual(score.score(far.xyz, far, rec), 0.0, places=6)
        inside = m.copy_with(m.xyz - m.xyz.mean(0) + rec.xyz[0])
        self.assertGreater(score.score(inside.xyz, inside, rec), 5.0)
        terms = score.pair_terms(inside.xyz, inside, rec)
        self.assertGreater(terms['repulsion'], 1.0)

    def test_stacking_is_face_to_face_and_not_edge_on(self):
        rec, m = self.rec, self.mol
        c, n = rec.ring_centre[10], rec.ring_normal[10]
        ring = m.xyz[:6] - m.xyz[:6].mean(0)
        normal = np.linalg.svd(ring - ring.mean(0))[2][2]
        R = _rotation_between(normal, n)
        stacked = ring @ R.T + c + 3.6 * n                       # the ring on the face of the base
        face = score.stacking(_pad(stacked, m), m, rec)
        edge = score.stacking(_pad(ring @ _rotation_between(normal, np.cross(n, [1.0, 0, 0])).T + c + 3.6 * n, m),
                              m, rec)
        self.assertGreater(face, 3.0)
        self.assertLess(edge, 0.25 * face)
        away = score.stacking(_pad(ring @ R.T + c + 9.0 * n, m), m, rec)
        self.assertLess(away, 0.05 * face)

    def test_the_terms_move_in_the_right_direction(self):
        # one 2′-hydroxyl as the whole receptor, and one acceptor walked towards it: the hydrogen-bond term
        # appears at contact and is gone a hydrogen bond and a half away
        lone = Receptor([[0.0, 0.0, 0.0]], ['O'], [0], ["O2'"], ['G'])
        self.assertTrue(lone.donor[0] and lone.acceptor[0])
        probe = chem.prepare(['O', 'C'], [[0.0, 0, 0], [1.4, 0, 0]])
        out = [score.pair_terms(probe.xyz - probe.xyz[0] + np.array([d, 0, 0]), probe, lone)['hbond']
               for d in (2.6, 3.0, 3.3, 5.0)]
        self.assertAlmostEqual(out[0], 2.0, places=6)         # at 2.6 Å, both directions count in full
        self.assertGreater(out[0], out[1])
        self.assertGreater(out[1], out[2])
        self.assertEqual(out[3], 0.0)
        far = score.pair_terms(probe.xyz + 50.0, probe, lone)
        self.assertEqual(sum(abs(v) for v in far.values()), 0.0)

    def test_charge_sign_decides_the_electrostatics(self):
        rec = self.rec
        here = rec.xyz[rec.charge.argmin()] + np.array([4.0, 0.0, 0.0])      # next to a phosphate oxygen
        plus = chem.prepare(['N', 'C'], [[0, 0, 0], [1.4, 0, 0]], formal=[1.0, 0.0])
        minus = chem.prepare(['O', 'C'], [[0, 0, 0], [1.4, 0, 0]], formal=[-1.0, 0.0])
        e_plus = score.pair_terms(plus.xyz - plus.xyz[0] + here, plus, rec)['electrostatic']
        e_minus = score.pair_terms(minus.xyz - minus.xyz[0] + here, minus, rec)['electrostatic']
        self.assertLess(e_plus, 0)
        self.assertGreater(e_minus, 0)


class TestGridAndSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = Receptor.from_pdb(SAMPLES / '1Y26_A' / 'reference.pdb')
        from rnafold3d.screening import pockets as P
        found = P.find(cls.rec.xyz, elements=cls.rec.element, residue_of=cls.rec.residue, top=1)
        cls.pocket = found[0]
        cls.grid = grid.Grid.for_pocket(cls.rec, pocket=cls.pocket, spacing=0.4)
        cls.mol = chem.read_sdf(BENZAMIDE)[0]
        cls.mol = cls.mol.copy_with(cls.mol.xyz - cls.mol.xyz.mean(0) + np.asarray(cls.pocket['center']))

    def test_the_grid_agrees_with_the_scoring_function(self):
        g, m, rec = self.grid, self.mol, self.rec
        prep = g.prepare_ligand(m)
        direct, bd = score.score(m.xyz, m, rec, breakdown=True)
        gridded, bg = g.energy(m.xyz, m, prep, breakdown=True)
        self.assertEqual(bg['outside'], 0.0)
        for term in ('gauss1', 'gauss2', 'hbond', 'stacking', 'electrostatic', 'desolvation'):
            # the electrostatic map is a 1/r² field, the one term whose curvature the grid feels: a fixed
            # tolerance on top of the relative one, since the value itself is often near zero
            self.assertAlmostEqual(bd['raw'][term], bg['raw'][term],
                                   delta=0.1 + 0.06 * abs(bd['raw'][term]), msg=term)
        self.assertAlmostEqual(direct, gridded, delta=0.5 + 0.15 * abs(direct))
        # the fast path folds the weights into the float32 maps, the term path weights the sums afterwards:
        # the same number up to the rounding of a float32
        self.assertAlmostEqual(g.energy(m.xyz, m, prep), gridded, delta=1e-4 * max(1.0, abs(gridded)))

    def test_many_poses_at_once_give_the_same_answer(self):
        g, m = self.grid, self.mol
        prep = g.prepare_ligand(m)
        poses = np.array([m.xyz + np.array([dx, 0.0, 0.0]) for dx in (0.0, 0.7, 1.4)])
        many = g.energies(poses, m, prep)
        self.assertEqual(many.shape, (3,))
        for k, p in enumerate(poses):
            self.assertAlmostEqual(float(many[k]), g.energy(p, m, prep), places=6)

    def test_a_pose_outside_the_box_is_rejected(self):
        g, m = self.grid, self.mol
        prep = g.prepare_ligand(m)
        self.assertGreater(g.energy(m.xyz + 100.0, m, prep, prep['weights']), 50.0)

    def test_the_gradient_follows_the_energy(self):
        d = dock.Docker(self.mol, self.grid)
        rng = np.random.default_rng(4)
        p = d.start(rng)
        e, gr = d.energy_and_gradient(p)
        self.assertAlmostEqual(e, d.energy(p), places=6)
        self.assertEqual(gr.shape, p.shape)
        step = np.zeros_like(p)
        step[0] = 0.05
        self.assertAlmostEqual((d.energy(p + step) - e) / 0.05, gr[0], delta=0.5 + 0.2 * abs(gr[0]))

    def test_docking_is_deterministic_and_stays_in_the_pocket(self):
        g, m, rec = self.grid, self.mol, self.rec
        a = dock.dock(m, g, receptor=rec, exhaustiveness=(2, 4), seed=11)
        b = dock.dock(m, g, receptor=rec, exhaustiveness=(2, 4), seed=11)
        c = dock.dock(m, g, receptor=rec, exhaustiveness=(2, 4), seed=12)
        self.assertTrue(a and len(a) <= 9)
        self.assertTrue(np.allclose(a[0].xyz, b[0].xyz))
        self.assertAlmostEqual(a[0].score, b[0].score, places=9)
        self.assertEqual([round(p.score, 6) for p in a], sorted(round(p.score, 6) for p in a))
        self.assertFalse(np.allclose(a[0].xyz, c[0].xyz) and a[0].score == c[0].score)
        for p in a:
            self.assertTrue(np.all(p.xyz.mean(0) > g.lo) and np.all(p.xyz.mean(0) < g.hi))
            self.assertLess(p.score, 10.0)
            for i, j, _o in m.bonds:                       # the compound keeps its own shape
                self.assertAlmostEqual(float(np.linalg.norm(p.xyz[i] - p.xyz[j])),
                                       float(np.linalg.norm(m.xyz[i] - m.xyz[j])), places=5)

    def test_rmsd_allows_the_symmetry_of_the_compound(self):
        m = chem.read_sdf(BENZAMIDE)[0]
        flipped = m.xyz.copy()
        flipped[[1, 2, 3, 4, 5]] = m.xyz[[5, 4, 3, 2, 1]]           # the ring turned over on itself
        self.assertGreater(dock.rmsd(flipped, m.xyz), 1.0)
        self.assertAlmostEqual(dock.rmsd(flipped, m.xyz, m.elements), 0.0, places=6)
        self.assertAlmostEqual(dock.rmsd(m.xyz + 1.0, m.xyz), float(np.sqrt(3)), places=6)
        with self.assertRaises(ValueError):
            dock.rmsd(m.xyz, m.xyz[:3])


def _pad(ring_xyz, m):
    """The molecule's coordinates with its ring put where `ring_xyz` says (for the stacking tests)."""
    out = np.array(m.xyz, float)
    out[:6] = ring_xyz
    out[6:] = ring_xyz.mean(0)
    return out


def _rotation_between(a, b):
    """A rotation that takes the unit vector `a` to the unit vector `b`."""
    a = np.asarray(a, float) / np.linalg.norm(a)
    b = np.asarray(b, float) / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(a @ b)
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)



class TestTwoStateScreen(unittest.TestCase):
    """The screen itself: a box per state, one score per state, and the order the window shows."""

    @classmethod
    def setUpClass(cls):
        from rnafold3d.screening import screen as SC
        cls.SC = SC
        cls.rec = Receptor.from_pdb(SAMPLES / '1Y26_A' / 'reference.pdb')
        from rnafold3d.screening import pockets as P
        cls.pockets = P.find(cls.rec.xyz, elements=cls.rec.element, residue_of=cls.rec.residue, top=2)
        cls.mol = chem.read_sdf(BENZAMIDE)[0]

    def test_a_pocket_gives_a_box_and_a_wide_one_is_cropped(self):
        SC = self.SC
        site = SC.Site.from_pocket(self.pockets[0])
        self.assertGreater(site.size, 6.0)
        self.assertEqual(len(site.centre), 3)
        self.assertTrue(set(site.residues) <= set(self.pockets[0]['residues']))
        self.assertEqual(site.cropped, site.size >= SC.MAX_SIZE)
        huge = SC.Site(centre=(0.0, 0.0, 0.0), size=200.0)
        self.assertTrue(huge.cropped)
        self.assertEqual(huge.size, SC.MAX_SIZE)       # a pocket spanning the molecule is not screened whole
        self.assertIn('cropped', site.as_dict())

    def test_every_state_is_scored_and_selectivity_is_b_minus_a(self):
        SC = self.SC
        a = SC.Site.from_atoms(self.rec.xyz[np.isin(self.rec.residue, self.pockets[0]['residues'][:6])])
        b = SC.Site(centre=tuple(np.asarray(a.centre) + np.array([6.0, 0.0, 0.0])), size=a.size)
        sc = SC.Screen(self.rec, a, self.rec, b, exhaustiveness=(1, 2), n_poses=1)
        row = sc.one(self.mol, 0, keep_poses=True)
        self.assertEqual(sorted(sc.states), ['a', 'b'])
        self.assertIsNotNone(row.score_a)
        self.assertIsNotNone(row.score_b)
        self.assertAlmostEqual(row.selectivity, row.score_b - row.score_a, places=9)
        self.assertEqual(row.worse, max(row.score_a, row.score_b))
        self.assertEqual(len(row.poses['a']), len(self.mol))
        d = row.as_dict(poses=True)
        self.assertEqual(d['atoms'], len(self.mol))
        self.assertIn('poses', d)
        self.assertNotIn('poses', row.as_dict())        # a table of 100,000 rows carries no coordinates

    def test_one_state_only_leaves_selectivity_empty(self):
        SC = self.SC
        a = SC.Site.from_pocket(self.pockets[0])
        row = SC.Screen(self.rec, a, exhaustiveness=(1, 2), n_poses=1).one(self.mol)
        self.assertIsNone(row.selectivity)
        self.assertIsNone(row.score_b)
        self.assertIsNotNone(row.worse)

    def test_the_order_is_selectivity_then_the_worse_score(self):
        SC = self.SC
        rows = [{'name': 'weak but selective', 'selectivity': 2.0, 'worse': -3.0},
                {'name': 'strong and selective', 'selectivity': 2.0, 'worse': -8.0},
                {'name': 'strong, not selective', 'selectivity': 0.5, 'worse': -9.0},
                {'name': 'not placed in B', 'selectivity': None, 'worse': -11.0}]
        order = [r['name'] for r in sorted(rows, key=SC.sort_key)]
        self.assertEqual(order, ['strong and selective', 'weak but selective', 'strong, not selective',
                                 'not placed in B'])

    def test_a_library_runs_in_one_process_and_in_several(self):
        SC = self.SC
        site = SC.Site.from_pocket(self.pockets[0])
        sc = SC.Screen(self.rec, site, exhaustiveness=(1, 2), n_poses=1, seed=3)
        mols = [self.mol, chem.read_sdf(GLYCINE)[0]]
        seen = []
        one = SC.run(sc, mols, jobs=1, progress=lambda d, t, r: seen.append((d, t)))
        self.assertEqual(seen, [(1, 2), (2, 2)])
        many = SC.run(sc, mols, jobs=2)
        self.assertEqual([r['name'] for r in one], [r['name'] for r in many])
        for x, y in zip(one, many):                     # the same seed, the same answer, however it is split
            self.assertAlmostEqual(x['score_a'], y['score_a'], places=6)
        s = SC.summarise(one)
        self.assertEqual(s['compounds'], 2)
        self.assertEqual(s['failed'], 0)
        self.assertEqual(len(s['top']), 2)

    def test_a_compound_too_big_for_a_screen_is_reported_not_dropped(self):
        SC = self.SC
        big = chem.prepare(['C'] * (SC.MAX_ATOMS + 1),
                           np.arange(3 * (SC.MAX_ATOMS + 1)).reshape(-1, 3) * 1.6, bonds=[], name='chain')
        row = SC.Screen(self.rec, SC.Site.from_pocket(self.pockets[0])).one(big)
        self.assertIn('too big', row.error)
        self.assertIsNone(row.score_a)
        scored = {'selectivity': -5.0, 'worse': -1.0}        # even a bad score sorts above no score at all
        self.assertEqual(sorted([row.as_dict(), scored], key=SC.sort_key)[0], scored)


if __name__ == '__main__':
    unittest.main()
