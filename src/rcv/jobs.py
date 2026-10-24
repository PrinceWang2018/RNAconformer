"""Local compute jobs: structure prediction from sequence and pull-guided exploration.

Jobs run the rnafold3d engine in a separate process on this computer's CPU cores (the same executable
in frozen builds, ``python -m rcv engine`` otherwise), report stages, progress, an estimated
time remaining and human-readable notes, and hand the finished ensemble to the viewer.
"""
from __future__ import annotations

import json
import os
import re
import secrets
import subprocess
import threading
import time
from pathlib import Path

import numpy as np

from . import __version__, rcvfile
from .methods import methods_text
from .paths import engine_command, engine_env

SAMPLING = {'fast': (2, 6), 'standard': (4, 8), 'thorough': (8, 16)}     # trajectories per model, minimum


def estimate_trajectory_seconds(n):
    """Rough single-trajectory wall time on one core (seconds), refined from live progress."""
    return 45.0 + 1.0 * n + 0.004 * n * n


def plural(n, word, many=None):
    n = int(n)
    return f'{n} {word if n == 1 else (many or word + "s")}'


def safe_name(text, default='rna'):
    name = re.sub(r'[^A-Za-z0-9_.-]+', '_', text or '').strip('_.')
    return (name or default)[:60]


class Job:
    kind = 'job'
    STAGES = ()

    def __init__(self, state, title):
        self.state = state
        self.id = secrets.token_hex(6)
        self.title = title
        self.status = 'queued'
        self.stages = [{'id': sid, 'label': label, 'state': 'pending', 'detail': ''} for sid, label in self.STAGES]
        self.progress = 0.0
        self.eta_s = None
        self.notes = []
        self.lines = []
        self.document = None
        self.result = {}
        self.error = None
        self.proc = None
        self.cancelled = False
        self.started = time.time()
        self.finished = None
        self._lock = threading.Lock()

    # ---------------------------------------------------------------- reporting
    def note(self, text, level='info'):
        with self._lock:
            self.notes.append({'t': round(time.time() - self.started, 1), 'level': level, 'text': str(text)[:500]})

    def stage(self, sid, state, detail=None):
        for s in self.stages:
            if s['id'] == sid:
                s['state'] = state
                if detail is not None:
                    s['detail'] = detail
            elif state == 'active' and s['state'] == 'active':
                s['state'] = 'done'

    def info(self):
        with self._lock:
            return {'id': self.id, 'kind': self.kind, 'title': self.title, 'status': self.status,
                    'stages': [dict(s) for s in self.stages], 'progress': round(self.progress, 3),
                    'eta_s': None if self.eta_s is None else round(self.eta_s), 'notes': list(self.notes[-60:]),
                    'elapsed_s': round((self.finished or time.time()) - self.started, 1), 'document': self.document,
                    'result': self.result, 'error': self.error, 'log_tail': self.lines[-15:]}

    # ---------------------------------------------------------------- lifecycle
    def start(self):
        threading.Thread(target=self._main, daemon=True, name=f'job-{self.id}').start()

    def cancel(self):
        self.cancelled = True
        p = self.proc
        if p and p.poll() is None:
            p.terminate()
            try:
                p.wait(5)
            except subprocess.TimeoutExpired:
                p.kill()

    def _main(self):
        self.status = 'running'
        try:
            self.run()
            if self.cancelled:
                raise RuntimeError('cancelled')
            self.status = 'done'
            self.progress = 1.0
            self.eta_s = 0
        except Exception as exc:                      # report every failure to the user
            if self.cancelled:
                self.status = 'cancelled'
                self.note('Cancelled by the user.', 'warning')
            else:
                self.status = 'failed'
                self.error = str(exc)
                self.note(str(exc), 'error')
            for s in self.stages:
                if s['state'] == 'active':
                    s['state'] = 'error' if self.status == 'failed' else 'pending'
        finally:
            self.finished = time.time()

    def run(self):
        raise NotImplementedError

    # ---------------------------------------------------------------- engine process
    def run_engine(self, args):
        flags = getattr(subprocess, 'CREATE_NO_WINDOW', 0) if os.name == 'nt' else 0   # no console flash on Windows
        self.proc = subprocess.Popen(engine_command(args), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     text=True, encoding='utf-8', errors='replace', env=engine_env(), bufsize=1,
                                     creationflags=flags)
        for raw in self.proc.stdout:
            line = raw.rstrip()
            if not line:
                continue
            with self._lock:
                self.lines.append(line)
                if len(self.lines) > 400:
                    del self.lines[:100]
            try:
                self.on_line(line)
            except Exception:                         # never let progress parsing kill a job
                pass
        code = self.proc.wait()
        if self.cancelled:
            raise RuntimeError('cancelled')
        if code != 0:
            tail = [l for l in self.lines[-6:] if l.strip()]
            raise RuntimeError('the prediction engine stopped with an error: ' + (' | '.join(tail) or f'exit code {code}'))

    # Engine log lines -> stages, progress, notes (format of rnafold3d.predict).
    def on_line(self, line):
        msg = line.split('] ', 1)[-1].strip() if line.startswith('[') else line.strip()
        m = re.match(r'(.+): (\d+) nt; (\d+) secondary-structure hypotheses \((.+)\)', msg)
        if m:
            self.stage('ss', 'done', f'{plural(m.group(3), "hypothesis", "hypotheses")} ({m.group(4)})')
            self.note(f'{plural(m.group(3), "secondary-structure hypothesis", "secondary-structure hypotheses")} from {m.group(4)}.')
            self.progress = max(self.progress, 0.05)
            return
        m = re.match(r'(ss\d+) (\S+) (\S+) kcal/mol \[(\w+)\]', msg)
        if m:
            src = {'vienna_mfe': 'minimum free energy', 'vienna_suboptimal': 'suboptimal', 'user': 'your input',
                   'seed': 'current conformation'}.get(m.group(4), m.group(4))
            e = '' if m.group(3) == 'None' else f', ΔG {m.group(3)} kcal/mol'
            self.note(f'{m.group(1)} ({src}{e}): {m.group(2)}')
            return
        m = re.match(r'running (\d+) trajectories on (\d+) worker', msg)
        if m:
            self.total, self.workers = int(m.group(1)), int(m.group(2))
            self.sample_t0 = time.time()
            self.stage('sample', 'active', f'0 of {self.total} trajectories')
            self.note(f'Sampling {plural(self.total, "independent trajectory", "independent trajectories")} on '
                      f'{plural(self.workers, "CPU core")}.')
            per = estimate_trajectory_seconds(getattr(self, 'n_nt', 50))
            self.eta_s = per * -(-self.total // max(1, self.workers)) + 10
            return
        m = re.match(r'progress (\d+)/(\d+) trajectories finished', msg)
        if m:
            k, n = int(m.group(1)), int(m.group(2))
            el = time.time() - getattr(self, 'sample_t0', self.started)
            self.progress = 0.08 + 0.82 * k / n
            self.eta_s = el / k * (n - k) + 8 if k else None
            self.stage('sample', 'active', f'{k} of {n} trajectories')
            return
        m = re.match(r'traj\s+(\d+) \[(ss\d+)\] (\w+)\s+score3d=(\S+)\s*(.*)', msg)
        if m:
            self.stage('select', 'active', 'validating and clustering')
            if m.group(3) != 'valid':
                why = re.sub(r"[\[\]']", '', m.group(5).rsplit(' ', 1)[0]).replace('_', ' ') or m.group(3)
                self.note(f'Trajectory {m.group(1)} ({m.group(2)}) rejected by validation: {why}.', 'warning')
            return
        m = re.match(r'(\d+)/(\d+) valid models; exported (\d+) conformers', msg)
        if m:
            self.result.update(valid=int(m.group(1)), trajectories=int(m.group(2)), conformers=int(m.group(3)))
            self.stage('select', 'done', f'{m.group(1)} of {m.group(2)} valid · {plural(m.group(3), "distinct conformation")}')
            self.note(f'{m.group(1)} of {m.group(2)} trajectories passed all-atom validation; '
                      f'{plural(m.group(3), "structurally distinct conformation")} kept.')
            self.progress = max(self.progress, 0.93)
            return
        if msg.startswith('WARNING:'):
            self.note(msg[8:].strip(), 'warning')
        elif 'error:' in msg:
            self.note(msg, 'error')


class PredictJob(Job):
    kind = 'predict'
    STAGES = (('ss', 'Secondary structures'), ('sample', '3D sampling'), ('select', 'Validation and selection'),
              ('load', 'Opening the viewer'))

    def __init__(self, state, params):
        seq = re.sub(r'\s+', '', params.get('sequence', '')).upper().replace('T', 'U')
        if not seq or not re.fullmatch(r'[ACGU]+', seq):
            bad = sorted(set(seq) - set('ACGU'))
            raise ValueError(f'the sequence may contain only A, C, G, U (T is read as U); found {"".join(bad) or "nothing"}')
        if len(seq) < 8:
            raise ValueError('the sequence is too short (at least 8 nt)')
        if len(seq) > 400:
            raise ValueError('sequences longer than 400 nt are outside the engine design range')
        structure = (params.get('structure') or '').strip()
        if structure and len(structure) != len(seq):
            raise ValueError(f'the dot-bracket has {len(structure)} characters but the sequence has {len(seq)} nt')
        name = safe_name(params.get('name'), 'rna')
        super().__init__(state, f'{name} · {len(seq)} nt')
        self.sequence, self.name, self.structure, self.n_nt = seq, name, structure, len(seq)
        self.models = max(1, min(50, int(params.get('models', 8))))
        per, low = SAMPLING.get(params.get('sampling', 'standard'), SAMPLING['standard'])
        self.trajectories = int(params.get('trajectories') or max(low, per * self.models))
        self.seed = int(params.get('seed', 42))
        cpu = os.cpu_count() or 1
        self.cores = max(1, min(cpu, int(params.get('cores') or max(1, cpu - 1))))
        self.out = Path(state.projects_dir) / f'{name}_{time.strftime("%Y%m%d-%H%M%S")}'
        per_traj = estimate_trajectory_seconds(len(seq))
        self.eta_s = per_traj * -(-self.trajectories // self.cores) + 15
        self.result = {'out': str(self.out), 'trajectories': self.trajectories, 'cores': self.cores}

    def run(self):
        self.out.parent.mkdir(parents=True, exist_ok=True)
        self.stage('ss', 'active', 'ViennaRNA' if not self.structure else 'from your dot-bracket')
        self.note(f'Prediction of {self.n_nt} nt started on this computer with {plural(self.cores, "CPU core")}; '
                  f'results are saved in {self.out}.')
        if self.n_nt > 200:
            self.note('Sequences above 200 nt are outside the tested range; expect long run times and lower accuracy.',
                      'warning')
        args = ['predict', '--sequence', self.sequence, '--name', self.name, '--models', str(self.models),
                '--trajectories', str(self.trajectories), '--seed', str(self.seed), '--workers', str(self.cores),
                '--out', str(self.out)]
        if self.structure:
            args += ['--structure', self.structure]
        self.run_engine(args)
        if not (self.out / 'metadata' / 'manifest.json').exists():
            raise RuntimeError('the engine finished without producing an ensemble')
        self.stage('load', 'active', 'computing energies and ensemble statistics')
        summary = self.state.open_paths([self.out])
        self.document = summary
        entry = self.state.entries.get(summary['id'])
        if entry is not None:
            try:
                from .server import vienna_version
                target = self.out / f'{self.name}.rcv'
                size = rcvfile.write(target, entry.doc, methods=methods_text(entry.doc, vienna_version()),
                                     provenance_manifest=entry.manifest, app_version=__version__)
                entry.rcv_path = str(target)
                self.result.update(rcv=str(target), rcv_bytes=size)
                self.note(f'Saved {target.name} ({size / 1e6:.2f} MB): all conformations, energies and 2D structures '
                          'in one file.')
            except Exception as exc:                     # the ensemble is open anyway; saving can be repeated
                self.note(f'Could not write the .rcv file: {exc}', 'warning')
        self.stage('load', 'done', plural(summary['n_conformations'], 'conformation'))
        if self.result.get('valid') == 0:
            raise RuntimeError('no trajectory produced a valid model; try more sampling or a known secondary structure')


# ------------------------------------------------------------------------------------------- pull
def pull_restraints(original, deformed, residue, max_pulled=3, max_anchors=5, tol=1.5):
    """Distance restraints (C4'-C4') that encode an interactively deformed shape.

    original, deformed: (n, 3) C4' coordinates before/after the drag; residue: the dragged residue (0-based).
    The dragged residue and up to two other strongly moved residues are tied to well-spread anchor residues
    that barely moved, with target distances taken from the deformed shape.
    """
    X0, X1 = np.asarray(original, float), np.asarray(deformed, float)
    n = len(X0)
    disp = np.linalg.norm(X1 - X0, axis=1)
    pulled = [int(residue)]
    for i in np.argsort(-disp):
        if len(pulled) >= max_pulled or disp[i] < 3.0:
            break
        if all(abs(int(i) - p) >= 3 for p in pulled):
            pulled.append(int(i))
    limit = max(1.0, 0.2 * float(disp.max()))
    cand = [i for i in range(n) if disp[i] <= limit and all(abs(i - p) >= 4 for p in pulled)]
    if len(cand) < 3:
        cand = [int(i) for i in np.argsort(disp)[: max(3, n // 3)] if all(abs(int(i) - p) >= 2 for p in pulled)]
    anchors = []
    if cand:
        start = max(cand, key=lambda i: float(np.linalg.norm(X1[i] - X1[residue])))
        anchors.append(start)
        while len(anchors) < min(max_anchors, len(cand)):
            far = max((i for i in cand if i not in anchors),
                      key=lambda i: min(float(np.linalg.norm(X1[i] - X1[a])) for a in anchors))
            anchors.append(far)
    out = []
    for p in pulled:
        for a in anchors:
            out.append({'i': p, 'j': a, 'a': "C4'", 'b': "C4'", 'target': round(float(np.linalg.norm(X1[p] - X1[a])), 2),
                        'tol': tol})
    return out, {'pulled': [p + 1 for p in pulled], 'anchors': [a + 1 for a in sorted(anchors)],
                 'max_displacement_A': round(float(disp.max()), 2)}


class ExploreJob(Job):
    """Sample new stable conformations that realise an interactive pull, then add them to the ensemble."""
    kind = 'explore'
    STAGES = (('ss', 'Secondary structures'), ('sample', 'Pull-guided sampling'), ('select', 'Relaxation and validation'),
              ('load', 'Adding to the ensemble'))

    def __init__(self, state, params):
        entry = state.entries.get(params.get('document'))
        if entry is None:
            raise ValueError('unknown ensemble')
        doc = entry.doc
        k = int(params['conformation'])
        residue = int(params['residue'])
        deformed = np.asarray(params['deformed'], float)
        n = len(doc['sequence'])
        if deformed.shape != (n, 3):
            raise ValueError('deformed trace must give one C4\' position per residue')
        names = doc['atoms']['names']
        c4 = []
        for r in doc['residues']:
            sl = names[r['first_atom']:r['first_atom'] + r['n_atoms']]
            if "C4'" not in sl:
                raise ValueError("every residue needs a C4' atom for pull exploration")
            c4.append(r['first_atom'] + sl.index("C4'"))
        import base64
        xyz = np.frombuffer(base64.b64decode(doc['conformations'][k]['coords']), '<f4').reshape(-1, 3)
        original = xyz[c4]
        restraints, info = pull_restraints(original, deformed, residue)
        super().__init__(state, f'Explore pull on {doc["residues"][residue]["name"]}{doc["residues"][residue]["num"]}')
        self.entry, self.k, self.residue, self.restraints, self.pull_info = entry, k, residue, restraints, info
        self.n_nt = n
        self.sequence = doc['sequence']
        self.seed_dbn = doc['conformations'][k]['dbn']
        cpu = os.cpu_count() or 1
        self.cores = max(1, min(cpu, int(params.get('cores') or max(1, cpu - 1))))
        self.trajectories = max(2, min(64, int(params.get('trajectories', 8))))
        self.models = max(1, min(self.trajectories, int(params.get('models', 4))))
        self.out = Path(state.projects_dir) / 'explorations' / f'{safe_name(entry.name)}_pull_{time.strftime("%Y%m%d-%H%M%S")}'
        self.eta_s = estimate_trajectory_seconds(n) * 1.15 * -(-self.trajectories // self.cores) + 20
        self.result = {'out': str(self.out), 'pull': info, 'restraints': len(restraints)}

    def run(self):
        self.out.mkdir(parents=True, exist_ok=True)
        self.stage('ss', 'active', 'current 2D + ViennaRNA alternatives')
        self.note(f'Pulling residue(s) {", ".join(map(str, self.pull_info["pulled"]))} relative to anchors '
                  f'{", ".join(map(str, self.pull_info["anchors"]))} (largest displacement '
                  f'{self.pull_info["max_displacement_A"]} Å): {len(self.restraints)} temporary distance restraints.')
        self.note('Each trajectory is sampled with the restraints as a bias, then relaxed without them, so the result '
                  'is a stable state of the unbiased model near the pulled shape.')
        job = {'seed_structures': [self.seed_dbn], 'pull_restraints': self.restraints, 'relax_steps': 1500,
               'max_structures': 2}
        (self.out / 'pull_job.json').write_text(json.dumps(job, indent=1))
        args = ['predict', '--sequence', self.sequence, '--name', 'pull', '--models', str(self.models),
                '--trajectories', str(self.trajectories), '--seed', str(int(time.time()) % 100000),
                '--workers', str(self.cores), '--out', str(self.out / 'run'), '--job-json', str(self.out / 'pull_job.json')]
        self.run_engine(args)
        self.stage('load', 'active', 'merging new conformations')
        summary, added = self.state.extend_document(self.entry, self.out / 'run')
        self.document = summary
        self.result['added'] = added
        self.stage('load', 'done', f'{len(added)} new conformation(s)')
        self.note(f'{len(added)} new conformation(s) added to the ensemble: {", ".join(added)}.')


# ------------------------------------------------------------------------------------------ batch
def parse_fasta(text):
    """FASTA (or one bare sequence) -> [{'name', 'sequence', 'structure'}]; a dot-bracket line after a sequence is
    taken as its known secondary structure."""
    records, cur = [], None
    for raw in str(text or '').splitlines():
        line = raw.strip()
        if not line or line.startswith((';', '#')):
            continue
        if line.startswith('>'):
            cur = {'name': line[1:].strip().split()[0] if line[1:].strip() else '', 'sequence': '', 'structure': ''}
            records.append(cur)
            continue
        if cur is None:
            cur = {'name': '', 'sequence': '', 'structure': ''}
            records.append(cur)
        token = line.split()[0]
        if re.fullmatch(r'[().\[\]{}<>]+', token) and cur['sequence']:
            cur['structure'] += token
        else:
            cur['sequence'] += re.sub(r'\s+', '', line).upper().replace('T', 'U')
    seen = {}
    for k, r in enumerate(records, 1):
        base = safe_name(r['name'], f'seq{k}')
        n = seen.get(base, 0)
        seen[base] = n + 1
        r['name'] = base if n == 0 else f'{base}_{n + 1}'
    return [r for r in records if r['sequence']]


class BatchJob(Job):
    """Predict several sequences one after another (each uses all chosen CPU cores) and save one .rcv per sequence."""
    kind = 'batch'
    MAX_SEQUENCES = 100

    def __init__(self, state, params):
        records = params.get('records') or parse_fasta(params.get('fasta', ''))
        # names become folder names: sanitise and de-duplicate whatever the client sent
        seen, clean = {}, []
        for k, r in enumerate(records, 1):
            base = safe_name(str(r.get('name') or ''), f'seq{k}')
            n = seen.get(base, 0)
            seen[base] = n + 1
            clean.append({'name': base if n == 0 else f'{base}_{n + 1}',
                          'sequence': re.sub(r'\s+', '', str(r.get('sequence') or '')).upper().replace('T', 'U'),
                          'structure': str(r.get('structure') or '').strip()})
        records = [r for r in clean if r['sequence']]
        if not records:
            raise ValueError('no sequences found in the FASTA input')
        if len(records) > self.MAX_SEQUENCES:
            raise ValueError(f'{len(records)} sequences; a batch can hold at most {self.MAX_SEQUENCES}')
        problems = []
        for r in records:
            seq = r['sequence']
            bad = sorted(set(seq) - set('ACGU'))
            if bad:
                problems.append(f'{r["name"]}: unsupported letters {"".join(bad)}')
            elif not 8 <= len(seq) <= 400:
                problems.append(f'{r["name"]}: {len(seq)} nt (8–400 nt supported)')
            elif r.get('structure') and len(r['structure']) != len(seq):
                problems.append(f'{r["name"]}: dot-bracket length {len(r["structure"])} differs from {len(seq)} nt')
        if problems:
            raise ValueError('; '.join(problems[:6]) + (f' (+{len(problems) - 6} more)' if len(problems) > 6 else ''))
        self.STAGES = tuple((f'item{k}', f'{r["name"]} · {len(r["sequence"])} nt') for k, r in enumerate(records))
        super().__init__(state, f'Batch of {plural(len(records), "sequence")}')
        self.records = records
        self.params = {k: params.get(k) for k in ('models', 'sampling', 'cores', 'seed')}
        cpu = os.cpu_count() or 1
        self.cores = max(1, min(cpu, int(params.get('cores') or max(1, cpu - 1))))
        self.out = Path(state.projects_dir) / f'batch_{time.strftime("%Y%m%d-%H%M%S")}'
        models = max(1, min(50, int(params.get('models', 8))))
        per, low = SAMPLING.get(params.get('sampling', 'standard'), SAMPLING['standard'])
        traj = max(low, per * models)
        self.estimates = [estimate_trajectory_seconds(len(r['sequence'])) * -(-traj // self.cores) + 25 for r in records]
        self.eta_s = sum(self.estimates)
        self.items = [{'name': r['name'], 'length': len(r['sequence']), 'sequence': r['sequence'], 'status': 'queued'}
                      for r in records]
        self.result = {'out': str(self.out), 'items': self.items, 'cores': self.cores}
        self.current = None

    def cancel(self):
        self.cancelled = True
        if self.current is not None:
            self.current.cancel()

    def info(self):
        cur = self.current
        if cur is not None and self.status == 'running':
            k = next((i for i, it in enumerate(self.items) if it['status'] == 'running'), None)
            if k is not None:
                done = sum(self.estimates[:k])
                self.progress = min(0.99, (done + self.estimates[k] * cur.progress) / sum(self.estimates))
                ratio = 1.0
                if k > 0 and self.items[k - 1].get('elapsed_s'):
                    spent = sum(it.get('elapsed_s', 0) for it in self.items[:k])
                    ratio = max(0.3, min(3.0, spent / max(1.0, done)))
                rest = sum(self.estimates[k + 1:]) * ratio
                self.eta_s = (cur.eta_s if cur.eta_s is not None else self.estimates[k]) + rest
                for s in cur.stages:
                    if s['state'] == 'active':
                        self.stage(f'item{k}', 'active', f'{s["label"]}' + (f' · {s["detail"]}' if s['detail'] else ''))
        return super().info()

    def run(self):
        self.out.mkdir(parents=True, exist_ok=True)
        self.note(f'Batch of {plural(len(self.items), "sequence")} started; each prediction uses '
                  f'{plural(self.cores, "CPU core")}. Results are saved in {self.out}.')
        for k, (rec, item) in enumerate(zip(self.records, self.items)):
            if self.cancelled:
                break
            params = dict(self.params, name=rec['name'], sequence=rec['sequence'], structure=rec.get('structure', ''),
                          cores=self.cores)
            job = PredictJob(self.state, params)
            job.out = self.out / rec['name']
            job.result['out'] = str(job.out)
            base_note = job.note

            def relay(text, level='info', _n=rec['name'], _base=base_note):
                _base(text, level)
                if level != 'info' or 'Saved' in str(text) or 'passed all-atom validation' in str(text):
                    self.note(f'{_n}: {text}', level)
            job.note = relay
            self.current = job
            item['status'] = 'running'
            self.stage(f'item{k}', 'active', 'starting')
            t0 = time.time()
            try:
                job.status = 'running'
                job.run()
                doc_summary = job.document
                entry = self.state.entries.pop(doc_summary['id'], None) if doc_summary else None
                self.state.invalidate_library()
                best = None
                if entry is not None:
                    scores = [c['values'].get('score') for c in entry.doc['conformations'] if c['values'].get('score') is not None]
                    best = min(scores) if scores else None
                    item['mfe_dbn'] = next((s.get('dbn') for s in (entry.doc.get('provenance') or {}).get('secondary_structures', [])
                                            if s.get('source') == 'vienna_mfe'), None)
                    item['dg_mfe'] = next((s.get('energy_kcal_mol') for s in (entry.doc.get('provenance') or {}).get('secondary_structures', [])
                                           if s.get('source') == 'vienna_mfe'), None)
                    item['distinct_2d'] = len({c['dbn'] for c in entry.doc['conformations']})
                    item['rmsd_max'] = round(max((v for row in entry.doc['ensemble']['rmsd_c4'] for v in row), default=0.0), 2)
                item.update(status='done', n_conformations=(doc_summary or {}).get('n_conformations'), best_score=best,
                            rcv=job.result.get('rcv'), out=str(job.out), valid=job.result.get('valid'),
                            trajectories=job.result.get('trajectories'))
                self.stage(f'item{k}', 'done', f'{plural(item["n_conformations"] or 0, "conformation")}')
                self.note(f'{rec["name"]}: {plural(item["n_conformations"] or 0, "conformation")} '
                          f'in {time.time() - t0:.0f} s.')
            except Exception as exc:
                if self.cancelled:
                    item['status'] = 'cancelled'
                    self.stage(f'item{k}', 'pending', 'cancelled')
                    break
                item.update(status='failed', error=str(exc))
                self.stage(f'item{k}', 'error', str(exc)[:120])
                self.note(f'{rec["name"]} failed: {exc}', 'error')
            finally:
                item['elapsed_s'] = round(time.time() - t0, 1)
                self.current = None
        self._write_summary()
        ok = [it for it in self.items if it['status'] == 'done']
        if self.cancelled:
            raise RuntimeError('cancelled')
        if not ok:
            raise RuntimeError('no sequence of the batch could be predicted')
        self.note(f'Batch finished: {len(ok)} of {plural(len(self.items), "sequence")} predicted. '
                  'Choose one to open in the viewer; the others stay in the batch folder.')

    def _write_summary(self):
        cols = ['name', 'length', 'status', 'n_conformations', 'best_score', 'dg_mfe', 'distinct_2d', 'rmsd_max',
                'elapsed_s', 'rcv', 'error', 'sequence', 'mfe_dbn']
        lines = [','.join(cols)]
        for it in self.items:
            row = []
            for c in cols:
                v = it.get(c)
                s = '' if v is None else str(v)
                row.append('"' + s.replace('"', '""') + '"' if any(ch in s for ch in ',"\n') else s)
            lines.append(','.join(row))
        try:
            (self.out / 'batch_summary.csv').write_text('\n'.join(lines) + '\n', encoding='utf-8')
            (self.out / 'batch.json').write_text(json.dumps({'items': self.items, 'parameters': self.params}, indent=1),
                                                 encoding='utf-8')
            self.result['summary_csv'] = str(self.out / 'batch_summary.csv')
        except OSError as exc:
            self.note(f'Could not write the batch summary: {exc}', 'warning')
