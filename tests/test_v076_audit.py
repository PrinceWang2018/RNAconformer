"""RCV 0.7.6, fixes from the audit of 2026-11-20 (RCV_0.7.6_audit.md): cancelling stops the whole computation,
structures that did not open, readable error messages, the local API, numbered sequences and the dot-bracket check."""
import ast
import http.client
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path

import numpy as np

from rnafold3d import structure_io as sio
from rnafold3d.annotate import check_dotbracket
from rnafold3d.build import atom_order
from rnafold3d.predict import normalize_sequence
from rnafold3d.structure_io import read_models, write_mmcif, write_pdb
from rcv import jobs
from rcv import server as srv
from rcv.data import collect_uploads, load_uploads

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / 'examples' / 'viewer_demo'
NMR = DEMO / 'HIV1_TAR_NMR_1ANR.pdb'
SRC = str(ROOT / 'src')


def demo_residues(n=8):
    ch = next(iter(read_models(NMR)[0][1].values()))
    res = [{k: v for k, v in r.atoms.items() if k in atom_order(r.base)} for r in ch.residues[:n]]
    return res, ch.sequence[:n]


def cif_text(residues, seq, comps=None):
    with tempfile.TemporaryDirectory() as d:
        write_mmcif(Path(d) / 'x.cif', residues, seq, comps=comps)
        return (Path(d) / 'x.cif').read_text()


def pdb_text(residues, seq, comps=None):
    with tempfile.TemporaryDirectory() as d:
        write_pdb(Path(d) / 'x.pdb', residues, seq, comps=comps)
        return (Path(d) / 'x.pdb').read_text()


def atom_rows(text):
    return [l for l in text.splitlines() if l.startswith(('ATOM ', 'HETATM '))]


# ------------------------------------------------------------------------------------------- A1
@unittest.skipUnless(os.name == 'posix', 'process inspection uses ps')
class TestCancelStopsTheWorkers(unittest.TestCase):
    SEQ = 'GGGAGCUCAGCCUUCACUGCAAUGGGCCGUAGACCCAUAGGCGAAACCGAAGCGCUAGUGCAAUAGGAGGGAUAGGCGAUUCCCGCCUAGGGAUAUGCUCACGAAGCUAGGCUCCC'

    def start(self, out):
        env = dict(os.environ, PYTHONPATH=SRC + os.pathsep + os.environ.get('PYTHONPATH', ''))
        p = subprocess.Popen([sys.executable, '-m', 'rnafold3d', 'predict', '--sequence', self.SEQ, '--models', '2',
                              '--trajectories', '4', '--workers', '2', '--templates', 'off', '--out', out],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, bufsize=1)
        for line in p.stdout:
            if 'running' in line and 'trajectories' in line:
                break
        deadline = time.time() + 30
        kids = []
        while time.time() < deadline and len(kids) < 2:
            time.sleep(0.5)
            kids = jobs._descendants(p.pid)
        self.assertGreaterEqual(len(kids), 2, 'the engine did not start its workers')
        return p, kids

    def gone(self, pids, within=8.0):
        deadline = time.time() + within
        while time.time() < deadline:
            table = jobs._process_table()
            left = [k for k in pids if k in table]
            if not left:
                return True
            time.sleep(0.3)
        for k in left:
            try:
                os.kill(k, signal.SIGKILL)
            except OSError:
                pass
        return False

    def test_cancel_ends_engine_and_workers(self):
        with tempfile.TemporaryDirectory() as d:
            p, kids = self.start(d)
            t0 = time.time()
            jobs.stop_process_tree(p)
            self.assertIsNotNone(p.poll())
            self.assertLess(time.time() - t0, 12)
            self.assertTrue(self.gone(kids), 'worker processes were still running after Cancel')
            p.stdout.close()

    def test_workers_stop_when_the_engine_is_killed(self):
        with tempfile.TemporaryDirectory() as d:
            p, kids = self.start(d)
            p.kill()
            p.wait(10)
            self.assertTrue(self.gone(kids), 'workers outlived the killed engine process')
            p.stdout.close()

    def test_launcher_script_and_its_child_are_stopped(self):
        with tempfile.TemporaryDirectory() as d:
            script = Path(d) / 'launch.sh'
            script.write_text('#!/bin/sh\nsleep 300 &\nwait\n')
            script.chmod(0o755)
            p = subprocess.Popen([str(script)])
            time.sleep(0.6)
            kids = jobs._descendants(p.pid)
            self.assertEqual(len(kids), 1)
            jobs.stop_process_tree(p, grace=3, all_at_once=True)
            self.assertTrue(self.gone(kids, 4))


# ------------------------------------------------------------------------------------------- A2
class TestExploreKeepsModifications(unittest.TestCase):
    def test_job_options_are_added_to_the_modifications(self):
        from rnafold3d import cli
        seen = {}

        def fake_predict(seq, out, name='query', dbn=None, config=None, log=None):
            seen.update(seq=seq, config=config)
            return {'conformers': [1], 'statistics': {}}
        with tempfile.TemporaryDirectory() as d:
            job = Path(d) / 'job.json'
            job.write_text(json.dumps({'seed_structures': ['((((....))))'], 'pull_restraints': [], 'relax_steps': 10,
                                       'max_structures': 2}))
            old, cli.predict = cli.predict, fake_predict
            try:
                with open(os.devnull, 'w') as null:
                    stdout, sys.stdout = sys.stdout, null
                    try:
                        code = cli.main(['predict', '--sequence', 'GGC[m6A]CUUCG[Y]GCC', '--out', d, '--job-json', str(job)])
                    finally:
                        sys.stdout = stdout
            finally:
                cli.predict = old
        self.assertEqual(code, 0)
        cfg = seen['config']
        self.assertEqual([(m['pos'], m['code']) for m in cfg.modifications], [(4, 'm6A'), (10, 'Y')])
        self.assertEqual(cfg.seed_structures, ('((((....))))',))
        self.assertEqual(cfg.max_structures, 2)


# -------------------------------------------------------------------------------------- A3, A4
class TestStructuresThatDidNotOpen(unittest.TestCase):
    def two_models_by_alternate(self, comps=None):
        """A file like 1ZEV: the alternate conformers A and C are given as models 1 and 2."""
        res, seq = demo_residues(6)
        text = cif_text(res, seq, comps)
        head, rows = text[:text.index('ATOM ' if 'ATOM ' in text else 'HETATM ')], atom_rows(text)

        def model(alt, num):
            out = []
            for r in rows:
                f = r.split()
                f[4], f[-1] = alt, str(num)
                out.append(' '.join(f))
            return out
        return head + '\n'.join(model('A', 1) + model('C', 2)) + '\n#\n', seq

    def test_models_labelled_by_alternate_conformer(self):
        text, seq = self.two_models_by_alternate()
        models = read_models(text=text, filename='alt.cif')
        self.assertEqual(len(models), 2)
        for _mid, chains, _meta in models:
            ch = next(iter(chains.values()))
            self.assertEqual(''.join(r.base for r in ch.residues), seq)
            self.assertGreater(len(ch.residues[2].atoms), 15)
        doc = load_uploads([{'name': 'alt.cif', 'text': text}])
        self.assertEqual(len(doc['conformations']), 2)

    def test_first_alternate_of_each_atom(self):
        res, seq = demo_residues(4)
        rows = atom_rows(pdb_text(res, seq))
        out = []
        for r in rows:                                   # conformer B before A for one residue: B is the first
            if r[22:26].strip() == '2':
                out.append(r[:16] + 'B' + r[17:30] + f'{float(r[30:38]) + 9:8.3f}' + r[38:])
                out.append(r[:16] + 'A' + r[17:])
            else:
                out.append(r)
        ch = next(iter(read_models(text='\n'.join(out) + '\nEND\n', filename='b.pdb')[0][1].values()))
        self.assertEqual(len(ch.residues), 4)
        self.assertEqual(len(ch.residues[1].atoms), len(res[1]))
        self.assertAlmostEqual(ch.residues[1].atoms["C4'"][0], res[1]["C4'"][0] + 9, places=2)

    def test_empty_models_are_left_out(self):
        res, seq = demo_residues(4)
        body = '\n'.join(atom_rows(pdb_text(res, seq)))
        water = 'HETATM    1  O   HOH A 101       1.000   2.000   3.000  1.00  0.00           O'
        text = f'MODEL        1\n{body}\nENDMDL\nMODEL        2\n{water}\nENDMDL\nMODEL        3\n{body}\nENDMDL\nEND\n'
        self.assertEqual([m[0] for m in read_models(text=text, filename='m.pdb')], ['1', '3'])

    def test_mirror_image_rna_and_parents_from_the_table(self):
        res, seq = demo_residues(6)                      # 2G32: L-nucleotides, no parent declared in the file
        comps = ['0' + b for b in seq]
        text = cif_text(res, seq, comps)
        text = re.sub(r"(0[ACGU]) 'RNA linking' [ACGU]", r"\1 'L-RNA linking' ?", text)
        ch = next(iter(read_models(text=text, filename='l.cif')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)
        self.assertEqual([r.comp_id for r in ch.residues], comps)
        ch = next(iter(read_models(text=pdb_text(res, seq, comps), filename='l.pdb')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)

    def test_dna_parent_inside_an_rna_chain(self):
        res, seq = demo_residues(6)                      # 3P4A: every residue declares a DNA parent
        names = {'A': 'XA1', 'C': 'XC1', 'G': 'XG1', 'U': 'XU1'}
        comps = [names[b] for b in seq]
        text = cif_text(res, seq, comps)
        text = re.sub(r"(X([ACGU])1) 'RNA linking' [ACGU]", lambda m: f"{m.group(1)} 'DNA linking' D{m.group(2)}", text)
        self.assertIn("'DNA linking' DG", text)
        ch = next(iter(read_models(text=text, filename='f.cif')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)
        # the same residues in a DNA chain stay out
        dna = text.replace('polyribonucleotide', 'polydeoxyribonucleotide')
        models = read_models(text=dna, filename='d.cif')
        self.assertEqual(models[0][1], {})
        self.assertEqual(sum(len(v) for v in models[0][2]['others']['dna'].values()), len(seq))

    def test_parent_read_from_the_base_atoms(self):
        res, seq = demo_residues(6)
        names = {'A': 'QA9', 'C': 'QC9', 'G': 'QG9', 'U': 'QU9'}
        comps = [names[b] for b in seq]
        text = re.sub(r"(Q[ACGU]9) 'RNA linking' [ACGU]", r"\1 'RNA linking' ?", cif_text(res, seq, comps))
        ch = next(iter(read_models(text=text, filename='q.cif')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)
        self.assertEqual(sio.base_from_atoms({"C1'", 'N9', 'C8', 'C4', 'N6'}), 'A')
        self.assertEqual(sio.base_from_atoms({"C1'", 'N9', 'C8', 'C4', 'O6', 'N2'}), 'G')
        self.assertEqual(sio.base_from_atoms({"C1'", 'N1', 'C6', 'C5', 'C4', 'N4'}), 'C')
        self.assertEqual(sio.base_from_atoms({"C1'", 'N1', 'C6', 'C5', 'C4', 'O4'}), 'U')
        self.assertIsNone(sio.base_from_atoms({'N9', 'C8', 'C4', 'N6'}))          # a free base: no sugar
        self.assertIsNone(sio.base_from_atoms({"C1'", "C2'", "O4'"}))             # no base

    def test_modres_and_seqres_of_pdb_files(self):
        res, seq = demo_residues(6)
        names = {'A': 'QA9', 'C': 'QC9', 'G': 'QG9', 'U': 'QU9'}
        comps = [names[b] for b in seq]
        text = pdb_text(res, seq, comps)                 # write_pdb adds MODRES records with the RNA parent
        ch = next(iter(read_models(text=text, filename='m.pdb')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)
        body = '\n'.join(l for l in text.splitlines() if not l.startswith('MODRES'))
        self.assertEqual(read_models(text=body, filename='n.pdb')[0][1], {})       # unknown components, no header
        seqres = 'SEQRES   1 A    6  ' + ' '.join(comps) + '\n'
        ch = next(iter(read_models(text=seqres + body, filename='s.pdb')[0][1].values()))
        self.assertEqual(''.join(r.base for r in ch.residues), seq)                # polymer residues: base from atoms
        # a DNA parent counts as RNA only with a 2' substituent
        dna = '\n'.join(f'MODRES XXXX {c} A {i:4d}   D{b}  X' for i, (c, b) in enumerate(zip(comps, seq), 1)) + '\n'
        ch = next(iter(read_models(text=dna + body, filename='d.pdb')[0][1].values()))
        self.assertEqual(len(ch.residues), len(seq))
        no2 = '\n'.join(l for l in body.splitlines() if l[12:16].strip() != "O2'")
        self.assertEqual(read_models(text=dna + no2, filename='e.pdb')[0][1], {})

    def test_free_base_is_not_a_residue(self):
        res, seq = demo_residues(5)
        rows = atom_rows(pdb_text(res, seq))
        base = [a for a in res[1] if "'" not in a and not a.startswith(('P', 'OP'))]
        lig = [f'HETATM{900 + k:5d} {a:<4s} ADE A 101    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {a[0]:>2s}'
               for k, (a, (x, y, z)) in enumerate((a, res[1][a] + 30.0) for a in base)]
        models = read_models(text='\n'.join(rows + lig) + '\nEND\n', filename='lig.pdb')
        ch = next(iter(models[0][1].values()))
        self.assertEqual(len(ch.residues), 5)
        self.assertEqual([k[3] for k in models[0][2]['others']['groups']], ['ADE'])


# ------------------------------------------------------------------------------------------- A5
class TestEngineFailureMessages(unittest.TestCase):
    def job(self):
        return jobs.PredictJob(type('S', (), {'projects_dir': tempfile.gettempdir()})(),
                               {'sequence': 'GGCACUUCGGUGCC', 'models': 2, 'cores': 1})

    def test_no_valid_model_is_said_in_words(self):
        job = self.job()
        for line in ['[00:00:41] 0/6 valid models; exported 0 conformers to /x in 40s', '{', ' "out": "/x",',
                     ' "conformers": 0,', ' "valid": 0,', ' "rejected": 6,', ' "errors": 0,', '}']:
            job.lines.append(line)
            job.on_line(line)
        msg = job.engine_failure(2)
        self.assertIn('passed the all-atom validation', msg)
        self.assertIn('more sampling', msg)
        self.assertNotIn('"rejected"', msg)

    def test_the_engines_own_error_line_is_shown(self):
        job = self.job()
        job.lines += ['[00:00:01] x: 14 nt', 'rnafold3d predict: error: the dot-bracket is not balanced: ")" at position 9 '
                      'has no opening "("']
        self.assertTrue(job.engine_failure(2).startswith('the dot-bracket is not balanced'))

    def test_json_fragments_are_not_shown(self):
        job = self.job()
        job.lines += ['Traceback (most recent call last):', 'MemoryError', '{', ' "valid": 3,', '}']
        msg = job.engine_failure(1)
        self.assertIn('MemoryError', msg)
        self.assertNotIn('"valid"', msg)


# -------------------------------------------------------------------------------------- A6, A7
class TestJobOptions(unittest.TestCase):
    state = type('S', (), {'projects_dir': tempfile.gettempdir()})()

    def test_batch_passes_zrna_on(self):
        job = jobs.BatchJob(self.state, {'fasta': '>a\nGGCACUUCGGUGCC\n>b\nCGCGCGCGCGCGAAAACGCGCGCGCGCG\n', 'models': 2,
                                         'zrna': 'auto', 'zrna_salt': 4.0, 'zrna_binder': True, 'cores': 1,
                                         'probing_text': '1 0.2\n'})
        self.assertEqual(job.params['zrna'], 'auto')
        self.assertEqual(job.unused, ['probing data'])
        one = jobs.PredictJob(self.state, dict(job.params, name='b', sequence=job.records[1]['sequence']))
        self.assertEqual((one.zrna, one.zrna_salt, one.zrna_binder), ('auto', 4.0, True))
        plain = jobs.BatchJob(self.state, {'fasta': '>a\nGGCACUUCGGUGCC\n>b\nGGCACUUCGGUGCC\n', 'cores': 1})
        self.assertNotIn('zrna', plain.params)
        self.assertGreater(sum(job.estimates), 0)

    def test_seed_zero_is_a_seed(self):
        self.assertEqual(jobs.PredictJob(self.state, {'sequence': 'GGCACUUCGGUGCC', 'seed': 0}).seed, 0)
        self.assertEqual(jobs.PredictJob(self.state, {'sequence': 'GGCACUUCGGUGCC'}).seed, 42)
        with self.assertRaisesRegex(ValueError, 'random seed'):
            jobs.PredictJob(self.state, {'sequence': 'GGCACUUCGGUGCC', 'seed': -5})
        js = (ROOT / 'src' / 'rcv' / 'static' / 'js' / 'welcome.js').read_text(encoding='utf-8')
        self.assertNotIn('seedIn.value || 42', js)

    def test_numbers_are_reported_in_words(self):
        for params, word in (({'models': 'abc'}, 'conformations to keep must be a whole number'),
                             ({'cores': 'many'}, 'CPU cores must be a whole number'),
                             ({'seed': 'x'}, 'random seed must be a whole number'),
                             ({'zrna': 'auto', 'zrna_salt': 'salty'}, 'must be a number')):
            with self.assertRaisesRegex(ValueError, word):
                jobs.PredictJob(self.state, dict({'sequence': 'GGCACUUCGGUGCC'}, **params))
        job = jobs.PredictJob(self.state, {'sequence': 'GGCACUUCGGUGCC'})           # nothing given: the defaults
        self.assertEqual((job.models, job.seed), (8, 42))
        self.assertGreaterEqual(job.cores, 1)


# ------------------------------------------------------------------------------------------- A8
class TestExtractStaysInItsFolder(unittest.TestCase):
    def test_names_from_the_file_cannot_leave_the_folder(self):
        from rcv import rcvfile
        from rcv.app import _extract
        from rcv.data import load_document
        doc = load_document([NMR])
        doc['conformations'] = doc['conformations'][:3]
        doc['conformations'][0]['id'] = '../../escaped'
        doc['conformations'][1]['id'] = '/abs/olute'
        doc['conformations'][2]['id'] = '..'
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / 'in' / 'x.rcv'
            f.parent.mkdir()
            rcvfile.write(f, doc)
            out = Path(d) / 'deep' / 'er' / 'out'
            with open(os.devnull, 'w') as null:
                stdout, sys.stdout = sys.stdout, null
                try:
                    self.assertEqual(_extract([str(f), '-o', str(out)]), 0)
                finally:
                    sys.stdout = stdout
            written = sorted(p.relative_to(d).as_posix() for p in Path(d).rglob('*.pdb'))
            self.assertEqual(len(written), 3)
            self.assertTrue(all(w.startswith('deep/er/out/') and w.count('/') == 3 for w in written), written)


# ------------------------------------------------------------------------------------------- B4
class TestFileMessages(unittest.TestCase):
    def message(self, name, text):
        with self.assertRaises(ValueError) as ctx:
            collect_uploads([{'name': name, 'text': text}])
        return str(ctx.exception)

    def test_truncated_pdb_line(self):
        res, seq = demo_residues(4)
        rows = atom_rows(pdb_text(res, seq))
        rows[10] = rows[10][:34]
        msg = self.message('cut.pdb', '\n'.join(rows) + '\n')
        self.assertIn('cut.pdb', msg)
        self.assertIn('line 11', msg)
        self.assertIn('not a complete ATOM record', msg)
        self.assertNotIn('could not convert', msg)

    def test_nan_coordinates(self):
        res, seq = demo_residues(4)
        rows = atom_rows(pdb_text(res, seq))
        rows[5] = rows[5][:30] + '     nan' + rows[5][38:]
        msg = self.message('nan.pdb', '\n'.join(rows) + '\n')
        self.assertIn('line 6', msg)
        self.assertIn('not a number', msg)
        text = cif_text(res, seq)
        first = atom_rows(text)[3]
        f = first.split()
        f[10] = 'nan'
        self.assertIn('atom 4', self.message('nan.cif', text.replace(first, ' '.join(f))))

    def test_incomplete_mmcif(self):
        res, seq = demo_residues(4)
        text = cif_text(res, seq)
        cut = text[:text.index('_atom_site.Cartn_x')] + text[text.index('_atom_site.Cartn_y'):]
        msg = self.message('part.cif', cut)
        self.assertIn('_atom_site.Cartn_x', msg)
        self.assertIn('incomplete', msg)
        self.assertNotIn("'Cartn_x'", msg)
        self.assertIn('no atom table', self.message('none.cif', 'data_x\n_entry.id x\n'))

    def test_file_without_rna_says_what_it_holds(self):
        ala = ['ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N',
               'ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00           C',
               'ATOM      3  C   ALA A   1      10.520   6.009  -4.116  1.00  0.00           C',
               'ATOM      4  O   ALA A   1       9.407   6.415  -4.405  1.00  0.00           O',
               'ATOM      5  CA  GLY A   2      10.700   5.400  -2.900  1.00  0.00           C',
               "ATOM      6  C4'  DA B   1       1.000   2.000   3.000  1.00  0.00           C",
               'HETATM    7 MG    MG A 101       4.000   2.000   3.000  1.00  0.00          MG',
               'HETATM    8  O   HOH A 201       7.000   2.000   3.000  1.00  0.00           O']
        msg = self.message('prot.pdb', '\n'.join(ala) + '\nEND\n')
        self.assertIn('no RNA found', msg)
        self.assertIn('protein (2 residues)', msg)
        self.assertIn('DNA (1 nucleotide)', msg)
        self.assertIn('MG', msg)
        self.assertEqual(sio.describe_others(sio.new_others()), '')

    def test_engine_command_line_numbers(self):
        from rnafold3d import cli
        for argv, word in ((['--models', '0'], '--models must be 1 or more'), (['--seed', '-1'], '--seed must be 0 or more'),
                           (['--modifications', 'abc:m6A'], 'POSITION:CODE')):
            with tempfile.TemporaryDirectory() as d:
                with tempfile.TemporaryFile('w+') as err:
                    stderr, sys.stderr = sys.stderr, err
                    try:
                        code = cli.main(['predict', '--sequence', 'GGCACUUCGGUGCC', '--out', d] + argv)
                    finally:
                        sys.stderr = stderr
                    err.seek(0)
                    self.assertEqual(code, 2)
                    self.assertIn(word, err.read())


# -------------------------------------------------------------------------------------- B5, B6
class TestServerHardening(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.state = srv.AppState(Path(cls.tmp) / 'runs', [DEMO], [Path(cls.tmp) / 'plugins'])
        cls.server = srv.AppServer(('127.0.0.1', 0), cls.state)
        cls.port = cls.server.server_address[1]
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        for j in cls.state.jobs.values():
            j.cancel()
        cls.server.shutdown()
        cls.server.server_close()
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def req(self, method, path, body=None, token=True, raw=None):
        c = http.client.HTTPConnection('127.0.0.1', self.port, timeout=60)
        headers = {'Host': f'localhost:{self.port}'}
        if token:
            headers['X-Viewer-Token'] = self.state.token if token is True else token
        data = raw
        if body is not None:
            data = json.dumps(body).encode()
            headers['Content-Type'] = 'application/json'
        c.request(method, path, body=data, headers=headers)
        r = c.getresponse()
        out = r.read()
        c.close()
        return r.status, json.loads(out or b'null')

    def test_body_must_be_a_json_object(self):
        st, r = self.req('POST', '/api/predict', [1, 2, 3])
        self.assertEqual(st, 400)
        self.assertIn('JSON object', r['error'])
        st, r = self.req('POST', '/api/predict', raw=b'{not json')
        self.assertEqual(st, 400)
        self.assertIn('not valid JSON', r['error'])

    def test_wrong_values_answer_400_with_words(self):
        st, r = self.req('POST', '/api/predict', {'sequence': 'GGCACUUCGGUGCC', 'models': 'abc'})
        self.assertEqual(st, 400)
        self.assertIn('whole number', r['error'])
        st, r = self.req('POST', '/api/predict', {'sequence': 'GGCACUUCGGUGCC', 'structure': '(((((....)))).'})
        self.assertEqual(st, 400)
        self.assertIn('never closed', r['error'])
        st, r = self.req('POST', '/api/batch', {'fasta': '>a\nGGCACUUCGGUGCC\n>b\nGGCACUUCGGUGCC\n', 'models': None,
                                                'cores': 'x'})
        self.assertEqual(st, 400)
        self.assertIn('CPU cores', r['error'])

    def test_token_is_checked_before_the_body_is_read(self):
        import socket
        s = socket.create_connection(('127.0.0.1', self.port), timeout=10)
        s.sendall((f'POST /api/predict HTTP/1.1\r\nHost: localhost:{self.port}\r\nContent-Type: application/json\r\n'
                   f'Content-Length: 400000000\r\n\r\n').encode())               # the body is never sent
        s.settimeout(5)
        answer = s.recv(4096).decode('latin-1')
        s.close()
        self.assertTrue(answer.startswith('HTTP/1.1 401'), answer[:80])
        self.assertIn('Connection: close', answer)
        self.assertEqual(self.req('GET', '/api/info', token='wrong-tökén')[0], 401)
        self.assertEqual(self.req('GET', '/api/ping')[1], {'ok': True})

    def test_window_gone_only_without_jobs(self):
        from rcv.app import window_gone, IDLE_QUIT_S
        st = type('S', (), {})()
        st.jobs, st.last_seen = {}, 1000.0
        self.assertFalse(window_gone(st, now=1000.0 + IDLE_QUIT_S - 1))
        self.assertTrue(window_gone(st, now=1000.0 + IDLE_QUIT_S + 1))
        st.jobs = {'a': type('J', (), {'status': 'running'})()}
        self.assertFalse(window_gone(st, now=1000.0 + 10 * IDLE_QUIT_S))
        st.jobs = {'a': type('J', (), {'status': 'done'})()}
        self.assertTrue(window_gone(st, now=1000.0 + 10 * IDLE_QUIT_S))


# -------------------------------------------------------------------------------------- C1, C2
class TestSequenceAndStructureInput(unittest.TestCase):
    def test_numbered_sequences(self):
        genbank = '        1 ggcagatctg agcctgggag ctctctgcc\n//\n'
        self.assertEqual(jobs.sequence_with_modifications(genbank.replace('//', ''))[0], 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC')
        self.assertEqual(normalize_sequence('1 ggcagatctg 11 agcctgggag'), 'GGCAGAUCUGAGCCUGGGAG')
        recs = jobs.parse_fasta('>tar\n' + genbank)
        self.assertEqual(recs[0]['sequence'], 'GGCAGAUCUGAGCCUGGGAGCUCUCUGCC')
        seq, mods = jobs.sequence_with_modifications('1 GGC[m6A]C 6 U[m2,2G]GCC')        # digits inside brackets stay
        self.assertEqual(seq, 'GGCACUGGCC')
        self.assertEqual([(m['pos'], m['code']) for m in mods], [(4, 'm6A'), (7, 'm2,2G')])
        with self.assertRaisesRegex(ValueError, 'unsupported letters'):
            jobs.sequence_with_modifications('GGCAXGGC')

    def test_dot_bracket_is_checked(self):
        self.assertIsNone(check_dotbracket('((((....))))', 12))
        self.assertIsNone(check_dotbracket('((..[[..))..]]'))
        self.assertIn('12 characters but the sequence has 14', check_dotbracket('((((....))))', 14))
        self.assertIn('")" at position 9 has no opening', check_dotbracket('(((....)))))'[1:]))
        self.assertIn('"(" at position 1 is never closed', check_dotbracket('((((....)))'))
        self.assertIn('character "x" at position 3', check_dotbracket('((x(....))))'))
        self.assertIn('hairpin loop of 2 nt', check_dotbracket('(((..)))'))
        state = type('S', (), {'projects_dir': tempfile.gettempdir()})()
        with self.assertRaisesRegex(ValueError, 'hairpin loop of 2 nt'):
            jobs.PredictJob(state, {'sequence': 'GGGGAACCCCAA', 'structure': '((((..))))..'})
        with self.assertRaisesRegex(ValueError, 'b: the dot-bracket is not balanced'):
            jobs.BatchJob(state, {'records': [{'name': 'a', 'sequence': 'GGCACUUCGGUGCC', 'structure': '(((((....)))))'},
                                              {'name': 'b', 'sequence': 'GGCACUUCGGUGCC', 'structure': '(((((....))).)'[::-1]}]})

    @unittest.skipUnless(shutil.which('node'), 'Node.js not installed')
    def test_the_window_checks_the_same(self):
        from rcv import ss2dexport
        from rcv.data import load_document
        doc = load_document([DEMO / 'HIV1_TAR_prediction'])
        cases = [{'format': f, 'records': ss2dexport.records(doc, w), 'text': ss2dexport.text(doc, f, w)}
                 for f in ss2dexport.FORMATS for w in ('distinct', 'all')]
        with tempfile.TemporaryDirectory() as d:
            bound = load_uploads([{'name': 'complex.pdb', 'text': complex_pdb(shift=5.0)}])
            (Path(d) / 'cases.json').write_text(json.dumps({'export': cases, 'complex': bound}))
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'audit_unit.mjs'), str(Path(d) / 'cases.json')],
                               capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout[-3000:] + r.stderr[-3000:])
        self.assertIn('all audit tests passed', r.stdout)


# ------------------------------------------------------------------- export of the 2D structure
class TestSecondaryStructureExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from rcv.data import load_document
        cls.doc = load_document([DEMO / 'HIV1_TAR_prediction'])

    def test_dot_bracket(self):
        from rcv import ss2dexport
        text = ss2dexport.text(self.doc, 'dbn', 'all')
        lines = text.splitlines()
        n = len(self.doc['conformations'])
        self.assertEqual(len(lines), 3 * n)
        for k in range(n):
            head, seq, dbn = lines[3 * k:3 * k + 3]
            self.assertTrue(head.startswith('>'))
            self.assertIn(self.doc['conformations'][k]['id'], head)
            self.assertEqual(seq, self.doc['sequence'])
            self.assertEqual(dbn, self.doc['conformations'][k]['dbn'])
            self.assertIsNone(check_dotbracket(dbn, len(seq)))
        distinct = ss2dexport.records(self.doc, 'distinct')
        self.assertEqual(len(distinct), len({c['dbn'] for c in self.doc['conformations']}))
        self.assertEqual(sum(1 + len(r['also']) for r in distinct), n)
        one = ss2dexport.text(self.doc, 'dbn', 'current', 2)
        self.assertEqual(one.count('>'), 1)
        self.assertIn(self.doc['conformations'][2]['dbn'], one)
        with self.assertRaisesRegex(ValueError, 'does not exist'):
            ss2dexport.text(self.doc, 'dbn', 'current', 999)
        with self.assertRaisesRegex(ValueError, 'unknown format'):
            ss2dexport.text(self.doc, 'xyz')

    def test_ct_and_bpseq_hold_the_same_pairs(self):
        from rcv import ss2dexport
        seq, n = self.doc['sequence'], len(self.doc['sequence'])
        for k, conf in enumerate(self.doc['conformations'][:4]):
            pairs = sorted((p[0] + 1, p[1] + 1) for p in conf['pairs'])
            ct = ss2dexport.text(self.doc, 'ct', 'current', k).splitlines()
            self.assertEqual(int(ct[0].split()[0]), n)
            rows = [l.split() for l in ct[1:]]
            self.assertEqual(len(rows), n)
            self.assertEqual(''.join(r[1] for r in rows), seq)
            self.assertEqual([int(r[0]) for r in rows], list(range(1, n + 1)))
            self.assertEqual([int(r[2]) for r in rows], list(range(0, n)))
            self.assertEqual([int(r[3]) for r in rows], list(range(2, n + 1)) + [0])
            self.assertEqual(sorted((int(r[0]), int(r[4])) for r in rows if int(r[4]) > int(r[0])), pairs)
            bp = [l.split() for l in ss2dexport.text(self.doc, 'bpseq', 'current', k).splitlines() if not l.startswith('#')]
            self.assertEqual(''.join(r[1] for r in bp), seq)
            self.assertEqual(sorted((int(r[0]), int(r[2])) for r in bp if int(r[2]) > int(r[0])), pairs)

    def test_ct_keeps_the_numbering_of_the_file(self):
        import copy
        from rcv import ss2dexport
        doc = copy.deepcopy(self.doc)
        for i, r in enumerate(doc['residues']):
            r['num'] = str(i + 17)
        self.assertEqual(ss2dexport.text(doc, 'ct', 'current', 0).splitlines()[1].split()[-1], '17')
        doc['residues'][3]['num'] = '19A'                                          # insertion code: plain 1..n instead
        self.assertEqual(ss2dexport.text(doc, 'ct', 'current', 0).splitlines()[1].split()[-1], '1')

    def test_pseudoknot_brackets(self):
        from rcv import ss2dexport
        self.assertEqual(ss2dexport.partners('((..[[..))..]]'), [9, 8, -1, -1, 13, 12, -1, -1, 1, 0, -1, -1, 5, 4])
        self.assertEqual(ss2dexport.partners('..)('), [-1, -1, -1, -1])

    def test_command_line(self):
        from rcv.app import main
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'tar.ct'
            with open(os.devnull, 'w') as null:
                stdout, sys.stdout = sys.stdout, null
                try:
                    self.assertEqual(main(['ss', str(NMR), '-o', str(out), '--conformation', '3']), 0)
                finally:
                    sys.stdout = stdout
            lines = out.read_text().splitlines()
            self.assertEqual(len(lines), 30)                                       # CT from the file name
            self.assertIn('Model 3', lines[0])


# ------------------------------------------------------------------- C5: what is bound to the RNA
def complex_pdb(n_models=3, n_res=10, ligand=True, ion=True, protein=True, shift=0.0):
    """Models of the NMR demo with a six-atom ligand, a magnesium ion, water and a short protein chain placed at the
    same spot relative to nucleotide 5 in every model; each model is also moved as a whole (shift * model number)."""
    models = read_models(NMR)[:n_models]
    out = []
    for mi, (_mid, chains, _meta) in enumerate(models):
        ch = next(iter(chains.values()))
        res = [{k: v for k, v in r.atoms.items() if k in atom_order(r.base)} for r in ch.residues[:n_res]]
        seq = ch.sequence[:n_res]
        move = np.array([shift * mi, -2 * shift * mi, shift * mi])
        rows = [l for l in atom_rows(pdb_text([{a: x + move for a, x in r.items()} for r in res], seq))]
        o = res[4]["C4'"] + move
        k = 9000

        def het(rec, name, comp, chain, num, xyz, elem):
            nonlocal k
            k += 1
            nm = name if len(name) == 4 else ' ' + name.ljust(3)
            return f'{rec:<6s}{k:5d} {nm} {comp:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {elem:>2s}'
        if ligand:
            for i, (nm, el) in enumerate((('C1', 'C'), ('C2', 'C'), ('N3', 'N'), ('C4', 'C'), ('O5', 'O'), ('C6', 'C'))):
                rows.append(het('HETATM', nm, 'LIG', 'A', 101, o + [6 + 1.4 * np.cos(i * np.pi / 3), 1.4 * np.sin(i * np.pi / 3), 3], el))
            rows.append(het('HETATM', 'H1', 'LIG', 'A', 101, o + [6, 0, 4], 'H'))
        if ion:
            rows.append(het('HETATM', 'MG', 'MG', 'A', 102, o + [-5, 4, 0], 'MG'))
        rows.append(het('HETATM', 'O', 'HOH', 'A', 201, o + [-5, 6, 0], 'O'))
        if protein:
            for i in range(8):
                rows.append(het('ATOM', 'CA', 'ALA', 'B', i + 1, o + [12 + 1.5 * np.cos(i * 1.75), 1.5 * np.sin(i * 1.75), 1.5 * i - 6], 'C'))
        out += [f'MODEL     {mi + 1:4d}'] + rows + ['ENDMDL']
    return '\n'.join(out) + '\nEND\n'


def context_xyz(doc, k):
    import base64
    return np.frombuffer(base64.b64decode(doc['context']['coords'][k]), '<f4').reshape(-1, 3)


def rna_xyz(doc, k):
    import base64
    return np.frombuffer(base64.b64decode(doc['conformations'][k]['coords']), '<f4').reshape(-1, 3)


class TestBatchCommand(unittest.TestCase):
    """`rcv batch` prints every note of the batch (F1: it went silent after 60) and Ctrl+C stops the computation."""

    def _run(self, job_class, argv_extra=()):
        import contextlib
        import io
        from unittest import mock
        from rcv import app, jobs
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / 'two.fa'
            fasta.write_text('>a\nGGGGAAAACCCC\n', encoding='utf-8')
            made = []

            def make(state, params):
                made.append(job_class(state, 'test batch'))
                made[-1].result = {'out': tmp}
                made[-1].items = [{'name': 'a', 'status': 'done', 'n_conformations': 1, 'rcv': 'a.rcv'}]
                return made[-1]
            out, err = io.StringIO(), io.StringIO()
            with mock.patch.object(jobs, 'BatchJob', make), contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = app._batch([str(fasta), '-o', tmp, *argv_extra])
        return code, out.getvalue(), made[0]

    def test_every_note_is_printed(self):
        from rcv.jobs import Job

        class Talkative(Job):
            def run(self):
                for k in range(100):
                    self.note(f'first burst {k}')
                time.sleep(1.4)                          # the command looks once a second
                for k in range(100):
                    self.note(f'second burst {k}')

        code, out, _job = self._run(Talkative)
        lines = [l for l in out.splitlines() if l.startswith('[')]
        self.assertEqual(code, 0)
        self.assertEqual(len(lines), 200)
        self.assertEqual([l.split('] ', 1)[1] for l in lines],
                         [f'first burst {k}' for k in range(100)] + [f'second burst {k}' for k in range(100)])
        self.assertEqual(json.loads(out.splitlines()[-1])['status'], 'done')

    def test_ctrl_c_cancels_the_job(self):
        from unittest import mock
        from rcv import app
        from rcv.jobs import Job
        cancelled = threading.Event()

        class Long(Job):
            def run(self):
                self.note('working')
                cancelled.wait(20)

            def cancel(self):
                self.cancelled = True
                cancelled.set()

        real_sleep, calls = time.sleep, []

        def sleep(seconds):
            calls.append(seconds)
            if len(calls) == 2:
                raise KeyboardInterrupt
            real_sleep(min(seconds, 0.2))
        with mock.patch.object(app.time, 'sleep', sleep):
            code, out, job = self._run(Long)
        self.assertEqual(code, 130)
        self.assertTrue(cancelled.is_set())
        self.assertEqual(job.status, 'cancelled')
        self.assertIn('working', out)
        self.assertEqual(json.loads(out.splitlines()[-1])['status'], 'cancelled')


class TestNotesInOtherLanguages(unittest.TestCase):
    """Audit F5: the notes of a running prediction are translated; paths, titles and names in them stay."""
    I18N = ROOT / 'src' / 'rcv' / 'static' / 'i18n'

    def test_a_note_names_its_data(self):
        job = jobs.Job(None, 'x')
        job.note('Saved my rna.rcv (1.20 MB): all conformations.', keep=['my rna.rcv', '', 'not in the text', 'my rna.rcv'])
        job.note('5 conformations in 12 s.', lead='trna_1: ')
        job.note('plain')
        job.note(f'results are saved in {Path("/a b/c")}.', 'warning', keep=[Path('/a b/c')])
        a, b, c, d = job.notes
        self.assertEqual(a['keep'], ['my rna.rcv'])
        self.assertEqual((b['text'], b['lead']), ('trna_1: 5 conformations in 12 s.', 'trna_1: '))   # `rcv batch` prints text
        self.assertEqual(set(c), {'t', 'level', 'text'})
        self.assertEqual((d['keep'], d['level']), ([str(Path('/a b/c'))], 'warning'))
        json.dumps(job.notes)

    def _texts(self, path, wanted):
        """The text pieces (12 characters or more) written in the first argument of the calls that `wanted` accepts."""
        def pieces(node):
            if isinstance(node, ast.Constant):
                return [node.value.strip()] if isinstance(node.value, str) and len(node.value.strip()) >= 12 else []
            if isinstance(node, ast.Subscript):                    # item["n_conformations"]: a key, not text
                return pieces(node.value)
            return [x for child in ast.iter_child_nodes(node) for x in pieces(child)]
        tree = ast.parse(path.read_text(encoding='utf-8'))
        return [x for node in ast.walk(tree) if isinstance(node, ast.Call) and node.args and wanted(node.func) for x in pieces(node.args[0])]

    def test_every_note_text_has_a_translation(self):
        keys = list(json.loads((self.I18N / 'de.json').read_text(encoding='utf-8'))['strings'])
        own = self._texts(ROOT / 'src' / 'rcv' / 'jobs.py', lambda f: isinstance(f, ast.Attribute) and f.attr == 'note')
        engine = self._texts(ROOT / 'src' / 'rnafold3d' / 'predict.py', lambda f: isinstance(f, ast.Attribute) and f.attr == 'append'
                             and isinstance(f.value, ast.Name) and f.value.id == 'warnings')
        self.assertGreater(len(own), 40)
        self.assertGreaterEqual(len(engine), 8)
        other_program = {'AlphaFold 3:'}                           # in front of a line that AlphaFold 3 printed
        missing = [c for c in own + engine if c not in other_program and not any(c in k for k in keys)]
        self.assertEqual(missing, [], 'note texts without an entry in the language files')

    def test_every_2d_source_of_the_engine_has_words(self):
        # audit F8: "ss1 (vienna_mfe_probing, ...)" was shown; every label the engine gives a 2D hypothesis is worded
        src = ROOT / 'src' / 'rnafold3d'
        labels = set()
        for name in ('ss.py', 'pseudoknots.py', 'predict.py'):
            text = (src / name).read_text(encoding='utf-8')
            labels |= set(re.findall(r"'source': '([A-Za-z_ ]+)'", text))
            labels |= {m for pair in re.findall(r"'source': '([a-z_]+)' if \w+ else '([a-z_]+)'", text) for m in pair}
        self.assertGreaterEqual(labels, {'vienna_mfe', 'vienna_mfe_probing', 'vienna_suboptimal', 'vienna_suboptimal_probing',
                                         'alignment_consensus', 'vienna_pseudoknot', 'seed', 'user', 'domain assembly'})
        keys = json.loads((self.I18N / 'de.json').read_text(encoding='utf-8'))['strings']
        for label in sorted(labels):
            words = jobs.hypothesis_source(label)
            self.assertNotIn('_', words, label)
            self.assertIn(words, keys, f'no translation entry for the source "{words}" ({label})')
        self.assertEqual(jobs.hypothesis_source('template 8UYJ_A+1I9F_A'), 'from experimental structure 8UYJ_A+1I9F_A')
        self.assertEqual(jobs.hypothesis_source('Z-RNA form of ss2 (template 1EHZ_A)'), 'Z-RNA form of ss2 (from experimental structure 1EHZ_A)')
        self.assertEqual(jobs.hypothesis_source('a_new_label'), 'a_new_label')
        self.assertIn('Z-RNA form of {0} ({1})', keys)

    def test_real_notes_in_every_language(self):
        job = jobs.PredictJob.__new__(jobs.PredictJob)
        jobs.Job.__init__(job, None, 'x')
        job.n_nt = 76
        title = "BtCoV-HKU5 5' proximal stem-loop 5, conformation 4"
        for line in ['[00:00:01] template search: 3 significant hit(s) in 0.2 s',
                     '[00:00:01] template search: 1 significant hit(s) in 0.4 s (released before 2024-05-01)',
                     '[00:00:01] template search: 0 significant hit(s) in 0.2 s',
                     '[00:00:01] template search: 0 significant hit(s) in 0.2 s (released before 2024-05-01)',
                     f"[00:00:01]   template 8UYJ_A E=6.2e-55 identity 100% coverage 100% pairs 46/46 [sequence] '{title}'",
                     f"[00:00:01]   template 8UYE_A E=6.2e-55 identity 100% coverage 100% pairs 46/46 [sequence, aligned with structure] '{title}'",
                     f"[00:00:01]   template 8UYM_A E=7.3e-32 identity 87% coverage 40% pairs 40/44 '{title}'",
                     '[00:00:01] x: 76 nt; 6 secondary-structure hypotheses (ViennaRNA)',
                     '[00:00:01] x: 76 nt; 1 secondary-structure hypotheses (user-supplied)',
                     '[00:00:01] x: 76 nt; 2 secondary-structure hypotheses (seed structures only (ViennaRNA unavailable))',
                     '[00:00:01]   ss1 ((((....)))) -7.3 kcal/mol [template 8UYJ_A]', '[00:00:01]   ss2 ((((....)))) -7.3 kcal/mol [vienna_mfe]',
                     '[00:00:01]   ss3 ((((....)))) -6.1 kcal/mol [vienna_suboptimal]', '[00:00:01]   ss4 ((((....)))) None kcal/mol [user]',
                     '[00:00:01]   ss5 ((((....)))) -6.1 kcal/mol [seed]', '[00:00:01]   ss6 ((((....)))) -6.1 kcal/mol [domain assembly]',
                     '[00:00:01]   ss7 ((((....)))) -7.9 kcal/mol [vienna_mfe_probing]', '[00:00:01]   ss8 ((((....)))) -7.1 kcal/mol [vienna_suboptimal_probing]',
                     '[00:00:01]   ss9 ((((....)))) -7.0 kcal/mol [alignment_consensus]', '[00:00:01]   ss10 (([[..)).]] -5.2 kcal/mol [vienna_pseudoknot]',
                     '[00:00:01]   ss11 ((((....)))) -2.2 kcal/mol [Z-RNA form of ss2 (vienna_mfe)]',
                     '[00:00:02] running 32 trajectories on 4 workers', '[00:00:02] running 1 trajectories on 1 worker',
                     "[00:01:02]   traj   3 [ss1] rejected score3d=None ['loop_closure_failed'] 93.1s",
                     '[00:01:03] 30/32 valid models; exported 5 conformers', '[00:01:03] 1/1 valid models; exported 1 conformers',
                     'WARNING: sequence longer than the 200-nt design range; results are untested at this size',
                     'WARNING: only 3 distinct valid conformers (requested 5); 12 of 32 trajectories produced valid models',
                     'WARNING: 4 crossing (pseudoknot) pairs are not part of the nested fold tree; they were dropped',
                     "WARNING: 1 non-canonical pairs in the user 2D were dropped (only A-U, G-C, G-U helices are modelled): [(3, 20, 'AG')]"]:
            job.on_line(line)
        self.assertEqual(len(job.notes), 30)
        self.assertEqual([n['keep'] for n in job.notes if 'keep' in n], [[title]] * 3 + [['loop closure failed']])
        texts = [n['text'] for n in job.notes]
        self.assertIn(f'Template 8UYJ_A ({title}): sequence identity 100%, covers 100% of the sequence, E-value 6.2e-55.', texts)   # F7
        self.assertIn('ss7 (minimum free energy with probing data, ΔG -7.9 kcal/mol): ((((....))))', texts)                         # F8
        self.assertIn('ss11 (Z-RNA form of ss2 (minimum free energy), ΔG -2.2 kcal/mol): ((((....))))', texts)
        self.assertFalse([t for t in texts if re.search(r'[a-z]_[a-z]', t.replace('8UYJ_A', '')) and 'loop' not in t], 'an internal label is shown')
        out = Path('/home/some user/RCV predictions/my rna_1')
        job.note(f'Prediction of 76 nt started on this computer with {jobs.plural(4, "CPU core")}; results are saved in {out}.', keep=[out])
        job.note(f'Batch of {jobs.plural(1, "sequence")} started; each prediction uses {jobs.plural(1, "CPU core")}. '
                 f'Results are saved in {out}.', keep=[out])
        job.note(f'{jobs.plural(5, "conformation")} in 118 s.', lead='trna_1: ')
        job.note(f'failed: {jobs.PredictJob.no_model_message}', 'error', lead='third ')
        job.note('Cancelled by the user.', 'warning')
        other = [{'text': 'AlphaFold 3: Running model inference with seed 1 took 42.1 seconds', 'same': True,
                  'keep': ['Running model inference with seed 1 took 42.1 seconds']},
                 {'text': 'fourth failed: the engine ended with exit code 1: something broke badly', 'lead': 'fourth ', 'same': True}]
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / 'notes.json').write_text(json.dumps(job.notes + other), encoding='utf-8')
            r = subprocess.run(['node', str(ROOT / 'tests' / 'js' / 'notes_unit.mjs'), str(Path(d) / 'notes.json')],
                               capture_output=True, text=True, encoding='utf-8')
        self.assertEqual(r.returncode, 0, r.stdout[-1500:] + r.stderr[-3000:])


class TestHtmlPageModules(unittest.TestCase):
    """The interactive HTML page carries its own copy of the viewer modules: every module one of them imports
    must be in the list (0.7.6: viewer.js began to import viewer/context.js and the page no longer started)."""

    def test_every_imported_module_is_bundled(self):
        import posixpath
        import re
        from rcv import htmlexport
        static = Path(htmlexport.__file__).parent / 'static'
        bundled = set(htmlexport.MODULES)
        missing = []
        for m in htmlexport.MODULES:
            if m.startswith('vendor/'):
                continue
            src = (static / m).read_text(encoding='utf-8')
            for spec in re.findall(r"""(?:from|import)\s*\(?\s*['"](\.{1,2}/[^'"]+)['"]""", src):
                target = posixpath.normpath(posixpath.join(posixpath.dirname(m), spec))
                if target not in bundled:
                    missing.append(f'{m} -> {target}')
        self.assertEqual(missing, [])


class TestBoundMolecules(unittest.TestCase):
    def test_ligand_ion_and_protein_are_kept_as_context(self):
        doc = load_uploads([{'name': 'complex.pdb', 'text': complex_pdb(shift=7.0)}])
        ctx = doc['context']
        self.assertEqual(len(doc['conformations']), 3)
        self.assertEqual(len(doc['sequence']), 10)                                  # the RNA itself is unchanged
        self.assertEqual([(g['comp'], g['kind'], g['n_atoms']) for g in ctx['groups']], [('LIG', 'ligand', 6), ('MG', 'ion', 1)])
        self.assertEqual(ctx['atoms']['names'], ['C1', 'C2', 'N3', 'C4', 'O5', 'C6', 'MG'])        # no hydrogen, no water
        self.assertEqual(ctx['atoms']['elements'][-1], 'MG')
        self.assertEqual(len(ctx['bonds']), 12)                                     # the six bonds of the ring
        self.assertTrue(all(a < 6 and b < 6 for a, b in zip(ctx['bonds'][::2], ctx['bonds'][1::2])))
        self.assertEqual([(t['kind'], t['chain'], t['n']) for t in ctx['traces']], [('protein', 'B', 8)])
        self.assertEqual(ctx['traces'][0]['comps'], ['ALA'] * 8)
        self.assertEqual(ctx['n_points'], 8)
        self.assertEqual(ctx['summary'], '1 ligand (LIG), 1 ion (MG), protein (1 chain B; 8 residues)')
        self.assertEqual(len(ctx['coords']), 3)
        # the models were moved apart in the file: after superposition everything bound keeps its place at the RNA
        c4 = doc['atoms']['names'].index("C4'", doc['residues'][4]['first_atom'])
        d = [np.linalg.norm(context_xyz(doc, k) - rna_xyz(doc, k)[c4], axis=1) for k in range(3)]
        self.assertEqual(d[0].shape, (15,))
        for k in (1, 2):
            np.testing.assert_allclose(d[k], d[0], atol=2e-3)
        self.assertAlmostEqual(float(d[0][6]), float(np.linalg.norm([-5, 4, 0])), places=2)
        self.assertFalse(any('Bound molecules' in w for w in doc['warnings']))

    def test_rna_only_files_have_no_context(self):
        doc = load_uploads([{'name': 'rna.pdb', 'text': complex_pdb(ligand=False, ion=False, protein=False)}])
        self.assertNotIn('context', doc)
        from rcv.data import load_document
        self.assertNotIn('context', load_document([DEMO / 'HIV1_TAR_prediction']))

    def test_only_what_every_conformation_has(self):
        both = [{'name': 'a.pdb', 'text': complex_pdb(1)}, {'name': 'b.pdb', 'text': complex_pdb(1, ion=False)}]
        doc = load_uploads(both)
        self.assertEqual([g['comp'] for g in doc['context']['groups']], ['LIG'])
        self.assertIn('Only the bound molecules that every conformation has are shown.', doc['warnings'])
        doc = load_uploads([both[0], {'name': 'c.pdb', 'text': complex_pdb(1, ligand=False, ion=False, protein=False)}])
        self.assertNotIn('context', doc)
        self.assertTrue(any('not every conformation has them' in w for w in doc['warnings']))

    def test_the_nearest_groups_when_there_are_too_many(self):
        from rcv import data
        old = data.CONTEXT_MAX_ATOMS, data.CONTEXT_MAX_POINTS
        try:
            data.CONTEXT_MAX_ATOMS = 6
            doc = load_uploads([{'name': 'complex.pdb', 'text': complex_pdb(2)}])
            self.assertEqual([g['comp'] for g in doc['context']['groups']], ['LIG'])            # nearer than the ion
            self.assertTrue(any('1 of 2 ligands and ions are shown' in w for w in doc['warnings']))
            data.CONTEXT_MAX_ATOMS, data.CONTEXT_MAX_POINTS = old[0], 4
            doc = load_uploads([{'name': 'complex.pdb', 'text': complex_pdb(2)}])
            self.assertEqual(doc['context']['traces'], [])
            self.assertEqual(len(doc['context']['groups']), 2)
        finally:
            data.CONTEXT_MAX_ATOMS, data.CONTEXT_MAX_POINTS = old

    def test_other_rna_chains_are_traces(self):
        text = complex_pdb(1, ligand=False, ion=False, protein=False)
        rows = [l for l in text.splitlines() if l.startswith('ATOM')]
        second = [l[:21] + 'C' + l[22:30] + f'{float(l[30:38]) + 40:8.3f}' + l[38:] for l in rows if int(l[22:26]) <= 6]
        doc = load_uploads([{'name': 'two.pdb', 'text': '\n'.join(rows + second) + '\nEND\n'}])
        self.assertEqual(len(doc['sequence']), 10)
        self.assertEqual([(t['kind'], t['chain'], t['n']) for t in doc['context']['traces']], [('rna', 'C', 6)])
        self.assertIn('other RNA (1 chain C; 6 nucleotides)', doc['context']['summary'])
        joined = load_uploads([{'name': 'two.pdb', 'text': '\n'.join(rows + second) + '\nEND\n'}], chain_id='A+C')
        self.assertEqual(len(joined['sequence']), 16)                               # both strands shown: nothing is left over
        self.assertNotIn('context', joined)

    def test_mmcif_ligand(self):
        res, seq = demo_residues(6)
        text = cif_text(res, seq).rstrip().rstrip('#').rstrip()
        o = res[2]["C4'"]
        extra = [f'HETATM {900 + i} {el} {nm} . SAM B 2 . ? {o[0] + 5 + i:.3f} {o[1]:.3f} {o[2] + 2:.3f} 1.00 0.00 301 SAM A {nm} 1'
                 for i, (nm, el) in enumerate((('N1', 'N'), ('C2', 'C'), ('SD', 'S')))]
        extra.append(f'HETATM 950 K K . K C 3 . ? {o[0]:.3f} {o[1] + 6:.3f} {o[2]:.3f} 1.00 0.00 302 K A K 1')
        extra.append(f'HETATM 951 O O . HOH D 4 . ? {o[0]:.3f} {o[1] + 8:.3f} {o[2]:.3f} 1.00 0.00 401 HOH A O 1')
        doc = load_uploads([{'name': 'lig.cif', 'text': text + '\n' + '\n'.join(extra) + '\n#\n'}])
        self.assertEqual(doc['sequence'], seq)
        self.assertEqual([(g['comp'], g['num'], g['chain'], g['kind']) for g in doc['context']['groups']],
                         [('SAM', '301', 'A', 'ligand'), ('K', '302', 'A', 'ion')])
        self.assertEqual(doc['context']['atoms']['elements'], ['N', 'C', 'S', 'K'])

    def test_rcv_files_keep_the_bound_molecules(self):
        from rcv import rcvfile
        from rcv.data import build_document
        doc = load_uploads([{'name': 'complex.pdb', 'text': complex_pdb()}])
        back, _session, _info = rcvfile.unpack(rcvfile.pack(doc))
        self.assertEqual(back['context'], doc['context'])
        confs = rcvfile.conformations_from_doc(back)
        self.assertEqual([g['comp'] for g in confs[1].context['groups']], ['LIG', 'MG'])
        again = build_document('again', confs)                                      # e.g. after a pull exploration
        self.assertEqual(again['context']['summary'], doc['context']['summary'])
        np.testing.assert_allclose(context_xyz(again, 2), context_xyz(doc, 2), atol=5e-3)
        confs[2].context = None                                                     # a new conformation without them
        self.assertNotIn('context', build_document('mixed', confs))

    def test_the_page_draws_and_follows_them(self):
        js = (ROOT / 'src' / 'rcv' / 'static' / 'js').resolve()
        self.assertIn('parseContext(doc.context', (js / 'ensemble.js').read_text(encoding='utf-8'))
        self.assertIn("label: 'Bound molecules (ligands, ions, other chains)'", (js / 'app.js').read_text(encoding='utf-8'))
        self.assertIn('viewer.setContextFrame(state.conf)', (js / 'lite.js').read_text(encoding='utf-8'))


if __name__ == '__main__':
    unittest.main()
