"""0.7.3: template search (alignment), template hypotheses, junction motifs, AlphaFold 3 import and local runs."""
import io
import json
import os
import sys
import tempfile
import time
import unittest
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rnafold3d import templates as T          # noqa: E402
from rnafold3d import motifs as M             # noqa: E402
from rnafold3d.structure_io import read_models  # noqa: E402

DEMO_PDB = ROOT / 'examples' / 'viewer_demo' / 'HIV1_TAR_NMR_1ANR.pdb'


def brute_sw(q, t):
    n, m = len(q), len(t)
    H = np.zeros((n + 1, m + 1), int)
    s = lambda a, b: (-1 if 'N' in (a, b) else (T.MATCH if a == b else T.MISMATCH))
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            v = max(0, H[i - 1][j - 1] + s(q[i - 1], t[j - 1]))
            for k in range(j):
                v = max(v, H[i][k] - T.GAP_OPEN - T.GAP_EXT * (j - k))
            for k in range(i):
                v = max(v, H[k][j] - T.GAP_OPEN - T.GAP_EXT * (i - k))
            H[i][j] = v
    return int(H.max())


class FakeLib:
    def __init__(self, seqs):
        self.n_chains = len(seqs)
        self.off = np.cumsum([0] + [len(x) for x in seqs])
        self.codes = T.CODE[np.frombuffer(''.join(seqs).encode(), np.uint8)]


class TestAlignment(unittest.TestCase):
    def test_scan_and_traceback_match_brute_force(self):
        rng = np.random.default_rng(3)
        for trial in range(12):
            q = ''.join(rng.choice(list('ACGU'), rng.integers(6, 20)))
            ts = [''.join(rng.choice(list('ACGU'), rng.integers(4, 26))) for _ in range(3)]
            ts[1] = ts[1][:4] + q[1:-1].replace('A', 'G', 1) + ts[1][4:]
            per = T.sw_scan(q, FakeLib(ts))
            for t, score in zip(ts, per):
                self.assertEqual(int(score), brute_sw(q, t))
                sc, path = T.align(q, t)
                self.assertEqual(sc, int(score))
                # the path re-scores to the reported score
                tot, prev = 0, None
                for a, b in path:
                    tot += T.MATCH if q[a] == t[b] else T.MISMATCH
                    if prev is not None:
                        g = (a - prev[0] - 1) + (b - prev[1] - 1)
                        if g:
                            tot -= T.GAP_OPEN + T.GAP_EXT * g
                    prev = (a, b)
                self.assertEqual(tot, sc)


class TestJunctions(unittest.TestCase):
    def test_ports_and_strands(self):
        # 3-way junction: closing pair (0, 39); children (3, 14) and (17, 36); strands 2, 2, 2
        pairs = [(0, 39), (1, 38), (2, 37), (3, 14), (4, 13), (5, 12), (17, 36), (18, 35), (19, 34)]
        js = M.junctions(pairs)
        self.assertEqual(len(js), 1)
        J = js[0]
        self.assertEqual(J['ports'], [(37, 2), (3, 14), (17, 36)])
        self.assertEqual(J['strands'], [0, 2, 0])
        self.assertEqual(J['depth'], [3, 3, 3])
        self.assertEqual(M.strand_residues(J, 1), [15, 16])

    def test_nested_greedy_removes_crossing(self):
        kept = M.nested_greedy([(0, 10), (1, 9), (5, 20), (12, 18)])
        self.assertNotIn((5, 20), kept)
        self.assertEqual(len(kept), 3)


@unittest.skipUnless(T.load_library() is not None, 'template library not installed')
class TestTemplateLibrary(unittest.TestCase):
    def test_trna_finds_trna_templates_before_cutoff(self):
        seq = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'   # yeast tRNA-Phe
        lib = T.load_library()
        hits = T.search(seq, lib, cutoff='2000-01-01', max_hits=3)
        self.assertTrue(hits)
        self.assertTrue(all(h['release'] < '2000-01-01' for h in hits))
        hy = T.hypothesis(seq, hits[0], lib)
        self.assertGreaterEqual(len(hy['pairs']), 15)          # cloverleaf
        self.assertGreater(len(hy['distance_restraints']['i']), 200)
        # the uncovered part is folded by ViennaRNA with the template pairs enforced (numpy ints broke this)
        self.assertIsNotNone(hy['energy_kcal_mol'])
        ta = hy['template_atoms']
        self.assertEqual(ta['names'][1], "C4'")
        self.assertGreater(len(ta['atoms']), 0.8 * len(hits[0]['aligned']))

    def test_partial_template_folds_the_rest(self):
        # tRNA-Phe with an unrelated 30-nt hairpin inserted into the anticodon loop: half of the insert is
        # not covered by any tRNA template and must still be folded (not left as one long loop)
        trna = 'GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA'
        seq = trna[:36] + 'GGCACUUCGGUGCCAAGGACUUCGGUCCUU' + trna[36:]
        lib = T.load_library()
        hits = T.search(seq, lib, cutoff='2000-01-01', max_hits=1)
        self.assertTrue(hits)
        hy = T.hypothesis(seq, hits[0], lib)
        self.assertIsNotNone(hy['energy_kcal_mol'])
        inserted = {i for p in hy['pairs'] for i in p if 36 <= i < 66}
        self.assertGreaterEqual(len(inserted), 8)

    def test_cutoff_excludes_everything_early(self):
        lib = T.load_library()
        self.assertEqual(T.search('GGGAAACCCAGGGAAACCCUUU', lib, cutoff='1950-01-01'), [])


# --------------------------------------------------------------------------------------- Cartesian refinement
def _heavy(ch):
    from rnafold3d.build import atom_order
    return [{a: np.array(v, float) for a, v in r.atoms.items() if a in atom_order(r.base)} for r in ch.residues]


class TestCartesian(unittest.TestCase):
    def setUp(self):
        models = read_models(DEMO_PDB)
        self.a = next(iter(models[0][1].values()))
        self.b = next(iter(models[5][1].values()))
        self.seq = ''.join(r.base for r in self.a.residues)

    def test_pull_onto_template_atoms_keeps_covalent_geometry(self):
        from rnafold3d.cartesian import refine_cartesian, template_targets
        start, other = _heavy(self.a), _heavy(self.b)
        names = ['P', "C4'", "C1'", 'NG', 'C2', 'C4', 'C6']
        rot = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])        # any frame: targets are superposed
        tmpl = {}
        for r, at in enumerate(other):
            g = 'N9' if self.seq[r] in 'AG' else 'N1'
            tmpl[r] = [list(at[g if nm == 'NG' else nm] @ rot.T + 7.0) if (g if nm == 'NG' else nm) in at else None
                       for nm in names]
        targets = template_targets(self.seq, start, tmpl, names)
        self.assertGreater(len(targets), 5 * len(self.seq))
        out, st = refine_cartesian(self.seq, start, targets)
        self.assertLess(st['rms_to_targets_A'], 0.5)
        self.assertLess(st['rms_to_targets_A'], 0.3 * st['rms_to_targets_start_A'] + 0.1)
        self.assertLess(st['max_covalent_deviation_A'], 0.15)
        for r in range(len(self.seq)):      # bonds kept
            self.assertAlmostEqual(np.linalg.norm(out[r]["C1'"] - out[r]["C2'"]),
                                   np.linalg.norm(start[r]["C1'"] - start[r]["C2'"]), delta=0.05)

    def test_regularize_repairs_a_strained_linkage(self):
        from rnafold3d.build import load_geometry
        from rnafold3d.cartesian import regularize
        res = _heavy(self.a)
        for nm in ('P', 'OP1', 'OP2'):                 # stretch the O3'(9)-P(10) bond by 0.4 A
            res[10][nm] = res[10][nm] + 0.4 * (res[10]['P'] - res[9]["O3'"]) / np.linalg.norm(res[10]['P'] - res[9]["O3'"])
        self.assertGreater(np.linalg.norm(res[10]['P'] - res[9]["O3'"]), 1.9)
        out, st = regularize(self.seq, res, load_geometry())
        self.assertAlmostEqual(float(np.linalg.norm(out[10]['P'] - out[9]["O3'"])), 1.602, delta=0.05)
        self.assertLess(st['atoms_moved_rms_A'], 0.3)


# --------------------------------------------------------------------------------------- AlphaFold 3
def _af3_cif(residues, seq, plddt):
    lines = ['data_fold_test', '#', '_entity_poly.entity_id 1', '_entity_poly.type polyribonucleotide',
             f'_entity_poly.pdbx_seq_one_letter_code_can {seq}', '#', 'loop_'] + [f'_atom_site.{k}' for k in (
        'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_seq_id', 'label_asym_id',
        'label_entity_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv',
        'auth_seq_id', 'auth_asym_id', 'pdbx_PDB_model_num')]
    k = 1
    for i, r in enumerate(residues, 1):
        for nm, xyz in r.atoms.items():
            q = f'"{nm}"' if "'" in nm else nm
            lines.append(f'ATOM {k} {nm[0]} {q} . {r.base} {i} A 1 ? {xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f} 1.00 '
                         f'{plddt[i - 1]:.2f} {i} A 1')
            k += 1
    return '\n'.join(lines + ['#']) + '\n'


def make_server_zip(path, n=5):
    models = read_models(DEMO_PDB)
    ch = next(iter(models[0][1].values()))
    seq = ''.join(r.base for r in ch.residues)
    with zipfile.ZipFile(path, 'w') as z:
        for k in range(n):
            ch = next(iter(models[k % len(models)][1].values()))
            plddt = [50 + k * 5 + (i % 7) for i in range(len(seq))]
            z.writestr(f'fold_tar_test/fold_tar_test_model_{k}.cif', _af3_cif(ch.residues, seq, plddt))
            z.writestr(f'fold_tar_test/fold_tar_test_summary_confidences_{k}.json',
                       json.dumps({'ptm': 0.5 + 0.01 * k, 'iptm': None, 'ranking_score': [0.3, 0.9, 0.5, 0.7, 0.1][k],
                                   'fraction_disordered': 0.0, 'has_clash': 0.0}))
        z.writestr('fold_tar_test/fold_tar_test_job_request.json', json.dumps([{'name': 'TAR test'}]))
        z.writestr('fold_tar_test/terms_of_use.md', 'terms')
    return seq


class TestAlphaFold3(unittest.TestCase):
    def test_server_zip(self):
        from rcv import alphafold
        from rcv.data import build_document
        from rcv.server import confidence_of
        with tempfile.TemporaryDirectory() as d:
            zp = Path(d) / 'fold_tar_test.zip'
            seq = make_server_zip(zp)
            name, confs, warnings, info = alphafold.collect_zip(zp)
            self.assertEqual(name, 'TAR test')
            self.assertEqual(len(confs), 5)
            self.assertEqual([c.meta['af3']['ranking_score'] for c in confs], [0.9, 0.7, 0.5, 0.3, 0.1])
            self.assertEqual(len(confs[0].meta['plddt']), len(seq))
            doc = build_document(name, confs, None, warnings, info)
            self.assertEqual(doc['origin'], 'alphafold3')
            self.assertEqual(doc['default_metric'], 'af3_ranking')
            self.assertIn('plddt', doc['conformations'][0]['values'])
            conf = confidence_of(doc)
            self.assertEqual(conf['info']['source'], 'alphafold3')
            self.assertAlmostEqual(conf['conformations'][0]['per_residue'][0], 0.55, places=2)
            # the same through the generic path reader (Open files / rcv convert)
            from rcv.data import collect_paths
            nm2, confs2, _, _, src = collect_paths([zp])
            self.assertEqual(src['kind'], 'alphafold3')
            self.assertEqual(len(confs2), 5)

    def test_local_output_folder(self):
        from rcv import alphafold
        with tempfile.TemporaryDirectory() as d:
            zp = Path(d) / 'x.zip'
            make_server_zip(zp, n=2)
            files = alphafold._files_from_zip(zp.read_bytes())
            job = Path(d) / 'out' / 'tar_test'
            for k, (s, smp) in enumerate([(1, 0), (1, 1)]):
                sub = job / f'seed-{s}_sample-{smp}'
                sub.mkdir(parents=True)
                (sub / f'tar_test_seed-{s}_sample-{smp}_model.cif').write_bytes(files[f'fold_tar_test/fold_tar_test_model_{k}.cif'])
                (sub / f'tar_test_seed-{s}_sample-{smp}_summary_confidences.json').write_bytes(
                    files[f'fold_tar_test/fold_tar_test_summary_confidences_{k}.json'])
            (job / 'tar_test_model.cif').write_bytes(files['fold_tar_test/fold_tar_test_model_1.cif'])
            (job / 'tar_test_ranking_scores.csv').write_text('seed,sample,ranking_score\n1,0,0.3\n1,1,0.9\n')
            self.assertTrue(alphafold.is_af3_dir(job))
            name, confs, _, _ = alphafold.collect_dir(job)
            self.assertEqual(len(confs), 2)                      # the top-level copy of the best model is skipped
            self.assertEqual(confs[0].label, 'AF3 seed 1 · sample 1')

    def test_command_and_input(self):
        from rcv import alphafold
        j = alphafold.input_json('x', 'GGGAAACCC', seeds=2)
        self.assertEqual(j['modelSeeds'], [1, 2])
        self.assertEqual(j['sequences'][0]['rna']['unpairedMsa'], '')
        cmd = alphafold.command({'python': 'py', 'script': 'run_alphafold.py', 'model_dir': '/m'}, '/in/x.json', '/out')
        self.assertIn('--norun_data_pipeline', cmd)
        self.assertIn('--json_path=/in/x.json', cmd)
        cmd = alphafold.command({'command': 'docker run -v {input_dir}:/i img --json_path=/i/{json_name} --output_dir={output_dir}'},
                                '/in/x.json', '/out')
        self.assertEqual(cmd[:3], ['docker', 'run', '-v'])
        self.assertIn('--json_path=/i/x.json', cmd)

    def test_local_job_with_stand_in_program(self):
        """A stand-in for run_alphafold.py writes AlphaFold 3 output; the job opens it as an ensemble."""
        from rcv import server as srv
        from rcv.jobs import AlphaFold3Job
        with tempfile.TemporaryDirectory() as d:
            zp = Path(d) / 'x.zip'
            make_server_zip(zp, n=3)
            fake = Path(d) / 'run_alphafold.py'
            fake.write_text(f'''
import json, sys, zipfile, pathlib
args = dict(a[2:].split('=', 1) for a in sys.argv[1:] if '=' in a)
job = json.load(open(args['json_path']))
assert '--norun_data_pipeline' in sys.argv and job['sequences'][0]['rna']['unpairedMsa'] == ''
print('Running model inference with seed 1', flush=True)
z = zipfile.ZipFile({str(zp)!r})
out = pathlib.Path(args['output_dir']) / job['name'].lower()
for k in range(3):
    sub = out / f'seed-1_sample-{{k}}'
    sub.mkdir(parents=True)
    (sub / 'model.cif').write_bytes(z.read(f'fold_tar_test/fold_tar_test_model_{{k}}.cif'))
    (sub / 'summary_confidences.json').write_bytes(z.read(f'fold_tar_test/fold_tar_test_summary_confidences_{{k}}.json'))
print('Fold job done', flush=True)
''')
            state = srv.AppState(Path(d) / 'runs', [], [])
            seq = ''.join(r.base for r in next(iter(read_models(DEMO_PDB)[0][1].values())).residues)
            job = AlphaFold3Job(state, {'name': 'tar', 'sequence': seq, 'engine': 'alphafold3',
                                        'af3': {'python': sys.executable, 'script': str(fake), 'model_dir': d, 'db_dir': ''}})
            job.start()
            for _ in range(600):
                if job.status in ('done', 'failed'):
                    break
                time.sleep(0.1)
            self.assertEqual(job.status, 'done', job.error)
            self.assertEqual(job.document['n_conformations'], 3)
            self.assertTrue(Path(job.result['rcv']).exists())


if __name__ == '__main__':
    unittest.main()
