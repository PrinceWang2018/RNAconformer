"""0.7.2 regressions: clash-free stems next to G-U wobbles; one-line errors from the rcv command-line tools."""
import contextlib
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from rnafold3d.build import ResidueBuilder, load_geometry, build_stem, PAIR_HBONDS   # noqa: E402


def _is_base(a):
    return "'" not in a and not a.startswith(('P', 'OP'))


def worst_contact(model):
    """Shortest distance between atoms of different residues that the validator would call a clash candidate."""
    keys = [(r, a) for r in sorted(model) for a in model[r]]
    pts = np.array([model[r][a] for r, a in keys])
    worst = np.inf
    for i, j in cKDTree(pts).query_pairs(2.6):
        (ra, na), (rb, nb) = keys[i], keys[j]
        if ra == rb or (abs(ra - rb) == 1 and not (_is_base(na) and _is_base(nb))):
            continue
        worst = min(worst, float(np.linalg.norm(pts[i] - pts[j])))
    return worst


class TestWobbleStems(unittest.TestCase):
    """Before 0.7.2 about 5% of stems containing a G-U wobble were built with two atoms closer than 2.2 A
    (a squeezed neighbouring pair, overlapping wobble bases or stacked neighbours). Stems are rigid during the
    search, so every trajectory on such a 2D hypothesis was rejected (e.g. the MFE fold of CASP16 R1205)."""

    CASES = [('UUGG', 'UCAG'),            # R1205 stem 24-27/40-43: G26-C41 O6-N4 was 2.05 A, U40/C41 2.12 A
             ('CUUUUG', 'UGGAAG'),        # terminal wobble: bases overlapped (1.47 A)
             ('GCUUG', 'UGGGC'),
             ('GUUG', 'UGGC')]

    @classmethod
    def setUpClass(cls):
        cls.rb = ResidueBuilder(load_geometry())

    def test_no_clash_and_pairs_kept(self):
        for s1, s2 in self.CASES:
            seq = s1 + s2
            n = len(s1)
            pairs = [(i, 2 * n - 1 - i) for i in range(n)]
            model, rep = build_stem(self.rb, seq, pairs)
            with self.subTest(stem=f'{s1}/{s2}'):
                self.assertGreaterEqual(worst_contact(model), 2.2)
                for i, j in pairs:
                    a, b = seq[i], seq[j]
                    hb = PAIR_HBONDS.get((a, b)) or [(y, x) for x, y in PAIR_HBONDS[(b, a)]]
                    for x, y in hb:
                        self.assertGreater(float(np.linalg.norm(model[i][x] - model[j][y])), 2.4)
                self.assertLess(rep['pair_fit_rmsd'], 0.6)


class TestToolErrors(unittest.TestCase):
    def run_tool(self, *argv):
        from rcv.app import main
        err = io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
            code = main(list(argv))
        return code, err.getvalue()

    def test_missing_file_is_one_line(self):
        code, err = self.run_tool('info', 'no_such_file.rcv')
        self.assertEqual(code, 2)
        self.assertEqual(err.count('\n'), 1)
        self.assertIn('rcv info: error:', err)
        self.assertNotIn('Traceback', err)

    def test_folder_given_to_info_suggests_convert(self):
        with tempfile.TemporaryDirectory() as d:
            code, err = self.run_tool('info', d)
        self.assertEqual(code, 2)
        self.assertIn('rcv convert', err)

    def test_not_an_rcv_file(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'x.pdb'
            p.write_text('ATOM\n')
            code, err = self.run_tool('extract', str(p), '-o', str(Path(d) / 'out'))
        self.assertEqual(code, 2)
        self.assertIn('not an .rcv file', err)

    def test_debug_variable_shows_traceback(self):
        os.environ['RCV_DEBUG'] = '1'
        try:
            with self.assertRaises(OSError):
                self.run_tool('info', 'no_such_file.rcv')
        finally:
            del os.environ['RCV_DEBUG']


if __name__ == '__main__':
    unittest.main()
