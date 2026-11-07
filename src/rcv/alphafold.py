"""AlphaFold 3 results in RCV: AlphaFold Server downloads (.zip), local AlphaFold 3 output folders, and running a
local AlphaFold 3 installation as a prediction engine.

AlphaFold Server zip:   fold_<job>_model_<k>.cif, fold_<job>_summary_confidences_<k>.json,
                        fold_<job>_full_data_<k>.json, fold_<job>_job_request.json (k = 0..4)
Local AlphaFold 3:      <out>/<job>/seed-<s>_sample-<k>/[<job>_seed-<s>_sample-<k>_]model.cif with
                        [..._]summary_confidences.json next to it, <job>_ranking_scores.csv, <job>_model.cif (best)

Every model becomes one conformation (best ranking score first). The per-atom pLDDT that AlphaFold 3 writes in
the B-factor column is averaged per nucleotide and shown as the confidence colouring; ranking score, pTM and
mean pLDDT are added as values. Only RNA chains are read (protein, DNA and ligands of the job are ignored).

AlphaFold 3 itself is not part of RCV: it needs its own installation (Linux, NVIDIA GPU) and model parameters
that Google DeepMind grants on request. RCV only writes the input file, starts the program and reads the result.
"""
from __future__ import annotations

import io
import json
import re
import shlex
import shutil
import subprocess
import zipfile
from pathlib import Path

from rnafold3d.structure_io import parse_cif, _mmcif_model, _mmcif_model_ids

MODEL_RE = re.compile(r'(?:^|[/\\_])model(?:_(\d+))?\.cif$', re.I)
SAMPLE_RE = re.compile(r'seed-(\d+)_sample-(\d+)', re.I)


class NotAlphaFold(ValueError):
    pass


# ----------------------------------------------------------------------------------------- reading
def _files_from_zip(data):
    zf = zipfile.ZipFile(io.BytesIO(data) if isinstance(data, (bytes, bytearray)) else data)
    out = {}
    for info in zf.infolist():
        if info.is_dir() or info.file_size > 200 * 1024 * 1024:
            continue
        name = info.filename.replace('\\', '/')
        if name.lower().endswith(('.cif', '.json', '.csv')):
            out[name] = zf.read(info)
    return out


def _files_from_dir(root):
    root = Path(root)
    out = {}
    for p in root.rglob('*'):
        if p.is_file() and p.suffix.lower() in ('.cif', '.json', '.csv') and p.stat().st_size < 200 * 1024 * 1024:
            out[p.relative_to(root).as_posix()] = p.read_bytes()
    return out


def looks_like_af3(names):
    names = list(names)
    has_model = any(MODEL_RE.search(n) for n in names)
    has_conf = any('summary_confidences' in n or 'ranking_scores' in n for n in names)
    return has_model and has_conf


def is_af3_dir(path):
    p = Path(path)
    if not p.is_dir() or (p / 'metadata' / 'manifest.json').exists():
        return False
    names = [q.name for q in p.glob('*')] + [q.relative_to(p).as_posix() for q in p.glob('*/*')]
    return looks_like_af3(names)


def _summary_for(model_name, files):
    """Name of the summary_confidences JSON that belongs to a model file."""
    cand = re.sub(r'model(_\d+)?\.cif$', lambda m: f'summary_confidences{m.group(1) or ""}.json', model_name, flags=re.I)
    if cand in files:
        return cand
    folder = model_name.rsplit('/', 1)[0] if '/' in model_name else ''
    for n in files:
        if n.startswith(folder) and n.endswith('summary_confidences.json') and n.count('/') == model_name.count('/'):
            return n
    return None


def _plddt_per_residue(cif, chain):
    """Mean pLDDT (B-factor column) of every residue of a chain, in residue order (label_seq_id)."""
    site = cif.get('_atom_site') or {}
    b = site.get('B_iso_or_equiv')
    if not b:
        return None
    vals = {}
    for asym, seq, v, mdl in zip(site.get('label_asym_id', []), site.get('label_seq_id', []), b,
                                 site.get('pdbx_PDB_model_num', ['1'] * len(b))):
        if asym != chain.label_asym_id or not seq.lstrip('-').isdigit():
            continue
        try:
            vals.setdefault(int(seq), []).append(float(v))
        except ValueError:
            continue
    return [round(sum(vals[r.index]) / len(vals[r.index]), 1) if r.index in vals else None for r in chain.residues]


def collect(files, source_label, chain_id=None):
    """files: {relative name: bytes}. Returns (name, confs, warnings, info)."""
    from .data import Conformation, _pick_chain, natural_key
    if not looks_like_af3(files):
        raise NotAlphaFold('no AlphaFold 3 models (model .cif with summary_confidences .json) found')
    warnings = set()
    models = [n for n in files if MODEL_RE.search(n)]
    sampled = [n for n in models if SAMPLE_RE.search(n)]
    if sampled:                        # local AlphaFold 3: the top-level <job>_model.cif repeats the best sample
        models = sampled
    job_name = None
    req = next((n for n in files if n.endswith('job_request.json')), None)
    if req:
        try:
            r = json.loads(files[req])
            r = r[0] if isinstance(r, list) else r
            job_name = r.get('name')
        except (ValueError, AttributeError, IndexError):
            pass
    data_json = next((n for n in files if n.endswith('_data.json') and not n.endswith('full_data.json')), None)
    if not job_name and data_json:
        try:
            job_name = json.loads(files[data_json]).get('name')
        except (ValueError, AttributeError):
            pass
    entries = []
    for n in models:
        summ = _summary_for(n, files)
        s = {}
        if summ:
            try:
                s = json.loads(files[summ])
            except ValueError:
                s = {}
        m = SAMPLE_RE.search(n)
        idx = MODEL_RE.search(n).group(1)
        label = (f'AF3 seed {m.group(1)} · sample {m.group(2)}' if m else
                 f'AF3 model {idx}' if idx is not None else 'AF3 model')
        entries.append({'name': n, 'summary': s, 'label': label, 'ranking': s.get('ranking_score')})
    entries.sort(key=lambda e: (-(e['ranking'] if e['ranking'] is not None else -1e9), natural_key(e['name'])))
    confs = []
    for rank, e in enumerate(entries, 1):
        text = files[e['name']].decode('utf-8', errors='replace')
        cif = parse_cif(text)
        chains, meta = _mmcif_model(cif, True, _mmcif_model_ids(cif)[0], e['name'])
        if not chains:
            raise ValueError(f'{e["name"]}: the AlphaFold 3 model contains no RNA chain')
        ch = _pick_chain(chains, chain_id, Path(e['name']).name, warnings)
        plddt = _plddt_per_residue(cif, ch)
        s = e['summary']
        known = [v for v in (plddt or []) if v is not None]
        extra = {'af3_ranking': s.get('ranking_score'), 'af3_ptm': s.get('ptm'), 'af3_iptm': s.get('iptm'),
                 'plddt': round(sum(known) / len(known), 1) if known else None,
                 'af3_disordered': s.get('fraction_disordered')}
        meta_c = {'id': f'af3_{rank:02d}', 'rank': rank, 'plddt': plddt,
                  'values_extra': {k: v for k, v in extra.items() if v is not None},
                  'af3': {k: s.get(k) for k in ('ranking_score', 'ptm', 'iptm', 'fraction_disordered', 'has_clash',
                                                'num_recycles') if k in s}}
        confs.append(Conformation(e['label'], e['name'], None, ch.auth_asym_id or ch.label_asym_id, ch.residues, meta_c,
                                  meta.get('title', '')))
    if not job_name:
        stem = Path(entries[0]['name']).name
        job_name = re.sub(r'(^fold_)|(_seed-\d+_sample-\d+)?_?model(_\d+)?\.cif$', '', stem, flags=re.I) or 'alphafold3'
    info = {'kind': 'alphafold3', 'label': source_label, 'models': len(confs), 'job': job_name}
    return job_name, confs, warnings, info


def collect_zip(path_or_bytes, chain_id=None):
    data = path_or_bytes if isinstance(path_or_bytes, (bytes, bytearray)) else Path(path_or_bytes).read_bytes()
    try:
        files = _files_from_zip(data)
    except zipfile.BadZipFile as exc:
        raise ValueError(f'not a valid zip file: {exc}') from exc
    return collect(files, 'AlphaFold Server (zip)', chain_id)


def collect_dir(path, chain_id=None):
    return collect(_files_from_dir(path), 'AlphaFold 3 (local)', chain_id)


# ----------------------------------------------------------------------------------- local install
DEFAULTS = {'python': '', 'script': '', 'model_dir': '', 'db_dir': '', 'command': '', 'seeds': 1, 'extra': ''}


def load_settings():
    from .paths import load_settings as _load
    return {**DEFAULTS, **(_load().get('alphafold3') or {})}


def save_settings(values):
    from .paths import save_settings as _save
    clean = {k: values.get(k, DEFAULTS[k]) for k in DEFAULTS}
    clean['seeds'] = max(1, min(20, int(clean.get('seeds') or 1)))
    _save({'alphafold3': clean})
    return clean


def input_json(name, sequence, seeds=1, msa_free=True):
    """AlphaFold 3 input (dialect alphafold3, version 1) for one RNA chain."""
    rna = {'id': 'A', 'sequence': sequence}
    if msa_free:
        rna['unpairedMsa'] = ''               # no MSA search: runs without the genetic databases
    return {'name': name, 'modelSeeds': list(range(1, int(seeds) + 1)), 'sequences': [{'rna': rna}],
            'dialect': 'alphafold3', 'version': 1}


def command(settings, json_path, out_dir):
    """Command line that runs AlphaFold 3 on json_path, writing into out_dir."""
    s = {**DEFAULTS, **settings}
    msa_free = not s.get('db_dir')
    if s.get('command'):
        fields = {'json': str(json_path), 'json_name': Path(json_path).name, 'input_dir': str(Path(json_path).parent),
                  'output_dir': str(out_dir), 'model_dir': s.get('model_dir', ''), 'db_dir': s.get('db_dir', '')}
        text = s['command']
        for k, v in fields.items():
            text = text.replace('{' + k + '}', v)
        return shlex.split(text, posix=True)
    if not s.get('script'):
        raise ValueError('AlphaFold 3 is not set up: choose run_alphafold.py (or enter a command) in the AlphaFold 3 settings')
    cmd = [s.get('python') or 'python3', s['script'], f'--json_path={json_path}', f'--output_dir={out_dir}']
    if s.get('model_dir'):
        cmd.append(f'--model_dir={s["model_dir"]}')
    cmd.append(f'--db_dir={s["db_dir"]}' if not msa_free else '--norun_data_pipeline')
    if s.get('extra'):
        cmd += shlex.split(s['extra'])
    return cmd


def check(settings):
    """Quick checks of an AlphaFold 3 setup (files exist, Python can import alphafold3). Returns a report."""
    s = {**DEFAULTS, **settings}
    items = []

    def add(ok, text):
        items.append({'ok': bool(ok), 'text': text})
    if s.get('command'):
        exe = shlex.split(s['command'])[0] if s['command'].strip() else ''
        add(bool(exe) and (shutil.which(exe) or Path(exe).exists()), f'program {exe or "(empty)"}')
        add('{json}' in s['command'] or '{json_name}' in s['command'], 'command contains {json} or {json_name}')
        add('{output_dir}' in s['command'], 'command contains {output_dir}')
    else:
        py = s.get('python') or 'python3'
        add(shutil.which(py) or Path(py).exists(), f'Python: {py}')
        add(s.get('script') and Path(s['script']).is_file(), f'run_alphafold.py: {s.get("script") or "(not set)"}')
        md = Path(s['model_dir']) if s.get('model_dir') else None
        add(md and md.is_dir() and any(md.glob('af3.bin*')), f'model parameters (af3.bin) in {s.get("model_dir") or "(not set)"}')
        if s.get('db_dir'):
            add(Path(s['db_dir']).is_dir(), f'databases: {s["db_dir"]}')
        else:
            items.append({'ok': None, 'text': 'no databases: AlphaFold 3 runs without MSA (faster, usually less accurate)'})
        if all(i['ok'] for i in items if i['ok'] is not None):
            try:
                p = subprocess.run([py, '-c', 'import alphafold3, jax; print(jax.devices())'], capture_output=True,
                                   text=True, timeout=120)
                gpu = 'cuda' in p.stdout.lower() or 'gpu' in p.stdout.lower()
                add(p.returncode == 0, 'alphafold3 and JAX import' + ('' if p.returncode == 0 else f': {p.stderr.strip()[-200:]}'))
                if p.returncode == 0:
                    add(gpu, f'GPU visible to JAX: {p.stdout.strip()[:120]}')
            except Exception as exc:
                add(False, f'could not start {py}: {exc}')
    return {'ok': all(i['ok'] for i in items if i['ok'] is not None), 'items': items}
