"""Build RCV ensemble documents (schema 'rcv.document/1').

Sources
    * rnafold3d ensemble packages (directory with metadata/manifest.json; engine runs and imports)
    * any set of PDB / mmCIF files, single- or multi-model (NMR ensembles, MD snapshots, other predictors)
    * in-memory uploads ({name, text}) from the GUI

All conformations must be the same RNA: residues are matched by author numbering and must carry
the same nucleotides. The document holds one shared atom table, per-conformation coordinates
(superposed on the first conformation), base pairs, helix families, energy terms (engine scores
when available plus the model-independent evaluator) and ensemble statistics.
"""
from __future__ import annotations

import base64
import datetime
import hashlib
import json
import re
from pathlib import Path

import numpy as np

import os

from rnafold3d.build import atom_order
from rnafold3d.evaluate import DESCRIPTORS, TERMS as EVAL_TERMS, TOTAL as EVAL_TOTAL, evaluate_many
from rnafold3d.structure_io import read_models

SCHEMA = 'rcv.document/1'
LEGACY_SCHEMAS = ('rnaconfswitch.document/1', 'rnafold3d.viewer/1')
STRUCT_EXT = ('.pdb', '.ent', '.cif', '.mmcif')
COVALENT = {'C': 0.76, 'N': 0.71, 'O': 0.66, 'P': 1.07, 'S': 1.05, 'F': 0.57, 'BR': 1.20, 'SE': 1.20, 'CL': 1.02, 'I': 1.39}
MAX_CONFORMATIONS = 500

ENGINE_TERMS = {
    'clash': ('Clash (coarse)', 'Coarse-grained steric overlap between rigid bodies and loops.'),
    'atom_clash': ('Clash (atoms)', 'All-atom overlap energy (Verlet list).'),
    'chainbreak': ('Chain closure', 'Residual chain-break penalty at loop closure points.'),
    'restraint': ('Restraints', 'Distance restraints (experimental pseudoknot option).'),
    'feasibility': ('Loop feasibility', 'Penalty when loop end-to-end distances exceed their reach.'),
    'compaction': ('Compaction', 'Radius-of-gyration restraint for branched RNAs.'),
    'stacking': ('Helix stacking', 'Reward for coaxially stacked helix ends.'),
    'hairpin': ('Hairpin geometry', 'Hairpin-loop geometry term.'),
    'electro': ('Electrostatics', 'Debye–Hückel phosphate repulsion (Debye length 10 Å).'),
    'torsion': ('Torsion prior', 'Knowledge-based backbone torsion prior of sampled residues.'),
    'loop_closure': ('Loop torsions + penalties', 'Torsion energy of closed loops plus 10 per severe clash.'),
    'ss_offset': ('2D offset', 'Free-energy gap of the conformer’s 2D hypothesis to the MFE (in kT).'),
}


# ----------------------------------------------------------------------------- input collection
class Conformation:
    __slots__ = ('label', 'source', 'model', 'chain', 'residues', 'meta', 'title')

    def __init__(self, label, source, model, chain, residues, meta=None, title=''):
        self.label, self.source, self.model, self.chain = label, source, model, chain
        self.residues, self.meta, self.title = residues, meta or {}, title


def _natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', str(s))]


def _pick_chain(chains, chain_id, where, warnings):
    if not chains:
        raise ValueError(f'{where}: no RNA residues found')
    if chain_id:
        for key, ch in chains.items():
            if chain_id in (key, ch.auth_asym_id, ch.label_asym_id):
                return ch
        raise ValueError(f'{where}: chain {chain_id!r} not found (available: {", ".join(sorted(chains))})')
    ch = max(chains.values(), key=lambda c: len(c.residues))
    others = [c for c in chains.values() if c is not ch and len(c.residues) >= 2]
    if others:
        warnings.add(f'{where}: several RNA chains; showing chain {ch.auth_asym_id or ch.label_asym_id} '
                     f'({len(ch.residues)} nt). Use the chain option to choose another.')
    return ch


def _models_to_confs(models, fname, chain_id, warnings, single_file=False):
    out = []
    multi = len(models) > 1
    stem = Path(fname).stem
    for mid, chains, meta in models:
        ch = _pick_chain(chains, chain_id, Path(fname).name, warnings)
        label = (f'Model {mid}' if single_file else f'{stem} · model {mid}') if multi else stem
        out.append(Conformation(label, Path(fname).name, mid if multi else None,
                                ch.auth_asym_id or ch.label_asym_id, ch.residues, {}, meta.get('title', '')))
    return out


def _package_confs(root, manifest, chain_id, warnings):
    confs = []
    for c in sorted(manifest['conformers'], key=lambda c: c.get('rank', 0)):
        files = c.get('files', {})
        rel = files.get('pdb') or files.get('mmcif')
        path = root / rel
        models = read_models(path)
        ch = _pick_chain(models[0][1], chain_id, path.name, warnings)
        confs.append(Conformation(c['id'], rel, None, ch.auth_asym_id or ch.label_asym_id, ch.residues, c))
    return confs


def collect_paths(paths, chain_id=None):
    """Read conformations from disk. Returns (name, confs, manifest, warnings, source)."""
    warnings = set()
    paths = [Path(p).expanduser() for p in paths]
    if len(paths) == 1 and paths[0].is_dir() and (paths[0] / 'metadata' / 'manifest.json').exists():
        root = paths[0]
        manifest = json.loads((root / 'metadata' / 'manifest.json').read_text())
        name = manifest.get('input', {}).get('name') or root.name
        return name, _package_confs(root, manifest, chain_id, warnings), manifest, warnings, \
            {'kind': 'package', 'path': str(root.resolve())}
    files = []
    for p in paths:
        if p.is_dir():
            files += sorted((f for f in p.iterdir() if f.suffix.lower() in STRUCT_EXT), key=_natural_key)
        elif p.exists():
            files.append(p)
        else:
            raise FileNotFoundError(f'not found: {p}')
    if not files:
        raise ValueError('no PDB/mmCIF files found')
    confs = []
    for f in files:
        confs += _models_to_confs(read_models(f), f.name, chain_id, warnings, single_file=len(files) == 1)
    name = files[0].stem if len(files) == 1 else (paths[0].name if paths[0].is_dir() else f'{files[0].stem} (+{len(files) - 1})')
    return name, confs, None, warnings, {'kind': 'files', 'paths': [str(f.resolve()) for f in files],
                                         'checksums': {f.name: _sha256(f) for f in files[:200]}}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def collect_uploads(files, chain_id=None):
    """files: [{'name': str, 'text': str}] (a package folder upload may include metadata/manifest.json)."""
    warnings = set()
    by_name = {f['name'].replace('\\', '/'): f['text'] for f in files}
    man_key = next((k for k in by_name if k.endswith('manifest.json')), None)
    if man_key:
        manifest = json.loads(by_name[man_key])
        if manifest.get('schema', '').startswith('rnafold3d.ensemble'):
            confs = []
            for c in sorted(manifest['conformers'], key=lambda c: c.get('rank', 0)):
                rel = c.get('files', {}).get('pdb') or c.get('files', {}).get('mmcif')
                key = next((k for k in by_name if k.endswith(rel)), None)
                if key is None:
                    continue
                models = read_models(text=by_name[key], filename=key)
                ch = _pick_chain(models[0][1], chain_id, rel, warnings)
                confs.append(Conformation(c['id'], rel, None, ch.auth_asym_id or ch.label_asym_id, ch.residues, c))
            if confs:
                name = manifest.get('input', {}).get('name') or 'ensemble'
                return name, confs, manifest, warnings, {'kind': 'upload-package'}
    confs = []
    names = sorted((k for k in by_name if k.lower().endswith(STRUCT_EXT)), key=_natural_key)
    if not names:
        raise ValueError('no PDB/mmCIF files among the uploaded files')
    for k in names:
        confs += _models_to_confs(read_models(text=by_name[k], filename=k), k, chain_id, warnings, single_file=len(names) == 1)
    stem = Path(names[0]).stem
    name = stem if len(names) == 1 else f'{stem} (+{len(names) - 1})'
    return name, confs, None, warnings, {'kind': 'upload', 'files': names}


# ----------------------------------------------------------------------------- geometry helpers
def _kabsch(mobile, target):
    """Rotation R and translation t minimising |mobile @ R + t - target|."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt)) or 1.0
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return r, tc - mc @ r


def _rmsd_fit(a, b):
    r, t = _kabsch(a, b)
    return float(np.sqrt(np.mean(np.sum((a @ r + t - b) ** 2, axis=1))))


def _element(name):
    n = name.upper()
    if n.startswith('OP'):
        return 'O'
    for two in ('BR', 'CL', 'SE'):
        if n.startswith(two):
            return two
    return n[0] if n[0] in 'CNOPSFI' else 'C'


def _bonds(xyz, elements, residue_of):
    """Covalent bonds from distances on one conformation (grid hashing)."""
    cell = 2.4
    keys = np.floor(xyz / cell).astype(int)
    grid = {}
    for i, k in enumerate(map(tuple, keys)):
        grid.setdefault(k, []).append(i)
    rad = np.array([COVALENT.get(e, 0.77) for e in elements])
    out = []
    for i, k in enumerate(map(tuple, keys)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((k[0] + dx, k[1] + dy, k[2] + dz), ()):
                        if j <= i or abs(residue_of[i] - residue_of[j]) > 1:
                            continue
                        d = float(np.linalg.norm(xyz[i] - xyz[j]))
                        if 0.4 < d < rad[i] + rad[j] + 0.45:
                            out.append((i, j))
    return out


def _helices(pairs):
    """Maximal stacked runs (i,j),(i+1,j-1),... of a pair list; returns list of pair-lists."""
    pset = {(i, j) for i, j in pairs}
    out, seen = [], set()
    for i, j in sorted(pset):
        if (i, j) in seen or (i - 1, j + 1) in pset:
            continue
        run = []
        a, b = i, j
        while (a, b) in pset:
            run.append((a, b))
            seen.add((a, b))
            a, b = a + 1, b - 1
        out.append(run)
    return out


def _nested_dbn(n, pairs):
    """Dot-bracket of the nested (non-pseudoknotted) pairs of a conformation."""
    s = ['.'] * n
    for i, j, _t, pk in pairs:
        if not pk and s[i] == '.' and s[j] == '.':
            s[i], s[j] = '(', ')'
    return ''.join(s)


def _layout2d(dbn):
    """Secondary-structure drawing coordinates (ViennaRNA NAView), centred, unit backbone step; None without ViennaRNA."""
    try:
        import RNA
        pts = RNA.naview_xy_coordinates(dbn)
    except Exception:
        return None
    xy = np.array([[p.X, p.Y] for p in list(pts)[:len(dbn)]], float)
    if len(xy) != len(dbn) or not np.isfinite(xy).all():
        return None
    step = np.median(np.linalg.norm(np.diff(xy, axis=0), axis=1)) if len(xy) > 1 else 1.0
    xy = (xy - xy.mean(axis=0)) / (step or 1.0)
    ends = (xy[0] + xy[-1]) / 2                 # rotate so the 5'/3' ends sit at the bottom and the fold grows upward
    v = -ends
    if np.linalg.norm(v) > 1e-6 and len(xy) > 2:
        a = np.arctan2(v[0], -v[1])             # angle of v from the screen 'up' direction (0, -1)
        c, s = np.cos(-a), np.sin(-a)
        xy = xy @ np.array([[c, s], [-s, c]])
    return [round(float(v), 3) for v in xy.reshape(-1)]


def _mds(dist):
    n = len(dist)
    if n < 3:
        return [[float(k), 0.0] for k in range(n)]
    d2 = np.asarray(dist) ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ d2 @ j
    w, v = np.linalg.eigh(b)
    order = np.argsort(w)[::-1][:2]
    xy = v[:, order] * np.sqrt(np.maximum(w[order], 0.0))
    if xy[0, 0] > 0:
        xy[:, 0] *= -1
    return [[round(float(a), 3), round(float(c), 3)] for a, c in xy]


# ----------------------------------------------------------------------------- engine scores
def _engine_values(meta):
    """Engine score and its additive decomposition (sums to total_score)."""
    if 'total_score' not in meta:
        return {}
    st = meta.get('score_terms', {})
    vals = {'score': meta['total_score'], 'score_3d': meta.get('score_3d')}
    if meta.get('ss_free_energy_kcal_mol') is not None:
        vals['ss_dg'] = meta['ss_free_energy_kcal_mol']
    acc = 0.0
    for k in ENGINE_TERMS:
        if k in ('loop_closure', 'ss_offset'):
            continue
        v = float(st.get(k, 0.0) or 0.0)
        vals[k] = round(v, 3)
        acc += v
    vals['loop_closure'] = round(float(meta.get('score_3d', 0.0)) - (float(st.get('total', acc)) - float(st.get('pull_bias', 0.0) or 0.0)), 3)
    vals['ss_offset'] = round(float(meta.get('ss_offset_kT', 0.0)), 3)
    return vals


# ----------------------------------------------------------------------------- main builder
def build_document(name, confs, manifest=None, warnings=(), source=None, plugins=(), progress=None):
    warnings = set(warnings)
    if not confs:
        raise ValueError('no conformations to show')
    if len(confs) > MAX_CONFORMATIONS:
        warnings.add(f'{len(confs)} conformations supplied; showing the first {MAX_CONFORMATIONS}.')
        confs = confs[:MAX_CONFORMATIONS]

    # --- residue correspondence: by position when every conformation has the same full sequence
    #     (e.g. NMR models numbered 17-45 next to engine models numbered 1-29), else by author numbering
    same_sequence = len({''.join(r.base for r in c.residues) for c in confs}) == 1

    def keys_of(conf):
        if same_sequence:
            return [(str(k + 1), '') for k in range(len(conf.residues))]
        if all(r.auth_seq_id for r in conf.residues):
            return [(r.auth_seq_id, r.ins_code) for r in conf.residues]
        return [(str(r.index), '') for r in conf.residues]

    ref = confs[0]
    ref_keys = keys_of(ref)
    common = set(ref_keys)
    per_conf_maps = []
    for c in confs:
        km = dict(zip(keys_of(c), c.residues))
        per_conf_maps.append(km)
        common &= set(km)
    order = [k for k in ref_keys if k in common]
    if len(order) < max(2, 0.5 * len(ref_keys)):
        raise ValueError('the structures do not share the same residues; all conformations must be the same RNA '
                         f'({len(order)} common of {len(ref_keys)} residues)')
    if len(order) < len(ref_keys) or any(len(m) > len(order) for m in per_conf_maps):
        warnings.add(f'{len(order)} residues are present in every conformation; residues missing from any '
                     'conformation are not shown.')
    seq = ''
    for k in order:
        bases = {m[k].base for m in per_conf_maps}
        if len(bases) > 1:
            raise ValueError(f'residue {k[0]}{k[1]} differs between conformations ({"/".join(sorted(bases))}); '
                             'all conformations must have the same sequence')
        seq += per_conf_maps[0][k].base

    # --- shared atom table (atoms present in every conformation)
    atom_names, atom_res = [], []
    residues_out = []
    for ri, k in enumerate(order):
        names = set(per_conf_maps[0][k].atoms)
        for m in per_conf_maps[1:]:
            names &= set(m[k].atoms)
        std = atom_order(seq[ri]) if seq[ri] in 'ACGU' else ()
        ordered = [a for a in std if a in names] + sorted(names - set(std))
        r0 = per_conf_maps[0][k]
        residues_out.append({'name': seq[ri], 'comp': r0.comp_id, 'num': f'{r0.auth_seq_id or k[0]}{r0.ins_code if r0.auth_seq_id else k[1]}',
                             'chain': ref.chain, 'first_atom': len(atom_names), 'n_atoms': len(ordered)})
        atom_names += ordered
        atom_res += [ri] * len(ordered)
    elements = [_element(a) for a in atom_names]
    n_atoms = len(atom_names)
    xyz = np.zeros((len(confs), n_atoms, 3))
    for ci, m in enumerate(per_conf_maps):
        for ri, k in enumerate(order):
            atoms = m[k].atoms
            s = residues_out[ri]['first_atom']
            for q in range(residues_out[ri]['n_atoms']):
                xyz[ci, s + q] = atoms[atom_names[s + q]]
    if progress:
        progress('superposing')

    # --- superposition on C4' (fallback: all atoms) of the first conformation, then centring
    c4 = [i for i, nm in enumerate(atom_names) if nm == "C4'"]
    fit_idx = c4 if len(c4) >= 3 else list(range(n_atoms))
    for ci in range(1, len(confs)):
        r, t = _kabsch(xyz[ci, fit_idx], xyz[0, fit_idx])
        xyz[ci] = xyz[ci] @ r + t
    mean = xyz.mean(axis=0)
    for ci in range(len(confs)):                         # one refinement round against the mean structure
        r, t = _kabsch(xyz[ci, fit_idx], mean[fit_idx])
        xyz[ci] = xyz[ci] @ r + t
    mean = xyz.mean(axis=0)
    centre = mean[fit_idx].mean(axis=0)
    xyz -= centre
    mean -= centre

    # --- ensemble statistics
    rmsf_atom = np.sqrt(np.mean(np.sum((xyz - mean) ** 2, axis=2), axis=0))
    rmsf = [round(float(rmsf_atom[r['first_atom']:r['first_atom'] + r['n_atoms']].mean()), 3) if r['n_atoms'] else 0.0
            for r in residues_out]
    nc = len(confs)
    rmsd = np.zeros((nc, nc))
    if len(c4) >= 3:
        for a in range(nc):
            for b in range(a + 1, nc):
                rmsd[a, b] = rmsd[b, a] = _rmsd_fit(xyz[a, c4], xyz[b, c4])
    embedding = _mds(rmsd)

    # --- per-conformation analysis
    engine_run = bool(manifest and manifest.get('software', {}).get('engine', '').startswith('rnafold3d-v2'))
    engine_mode = engine_run or any('total_score' in c.meta for c in confs)
    plugin_meta = {}
    conformations = []
    all_helices = []
    residue_sets = []
    for ci in range(nc):
        residues = []
        for r in residues_out:
            s = r['first_atom']
            residues.append({atom_names[s + q]: xyz[ci, s + q] for q in range(r['n_atoms'])})
        residue_sets.append(residues)
    if progress:
        progress(f'evaluating {nc} conformations')
    # engine packages carry the validator report already; imported structures are validated here
    known = [c.meta.get('validation') if c.meta.get('validation', {}).get('severe_clashes') is not None else None for c in confs]
    pending = sum(k is None for k in known)
    work = pending * len(seq)               # validation dominates: ~1 s per 200-nt conformation
    workers = min(max(1, (os.cpu_count() or 2) - 1), pending) if (pending >= 6 and work > 3000) else 1
    try:
        evaluations = evaluate_many(seq, residue_sets, known, workers)
    except (OSError, RuntimeError):              # no process pool available: fall back to serial
        evaluations = evaluate_many(seq, residue_sets, known, 1)
    for ci, c in enumerate(confs):
        residues = residue_sets[ci]
        ev = evaluations[ci]
        values = dict(ev['values'])
        values.update(_engine_values(c.meta))
        values['rmsd_first'] = round(float(rmsd[0, ci]), 3)
        values['rmsd_mean'] = round(float(np.sqrt(np.mean(np.sum((xyz[ci, c4 or fit_idx] - mean[c4 or fit_idx]) ** 2,
                                                                  axis=1)))), 3)
        for plug in plugins:
            try:
                extra = plug.compute(seq, residues) or {}
            except Exception as exc:                          # a plugin must never break loading
                warnings.add(f'plugin {plug.name}: {type(exc).__name__}: {exc}')
                extra = {}
            for key, val in extra.items():
                values[f'plugin.{plug.id}.{key}'] = val
            for m in getattr(plug, 'metrics', ()):
                plugin_meta[f'plugin.{plug.id}.{m["key"]}'] = {**m, 'group': 'plugin', 'plugin': plug.name}
        validation = c.meta.get('validation') or ev['validation']
        if validation and 'overlaps_per_1000' not in validation and 'heavy_atom_overlaps_per_1000_atoms' in validation:
            validation = {**validation, 'overlaps_per_1000': validation['heavy_atom_overlaps_per_1000_atoms']}
        pairs = [(i - 1, j - 1, t) for i, j, t in ev['pairs']]
        pk = {(i - 1, j - 1) for i, j in ev['pk_pairs']}
        hel = _helices([(i, j) for i, j, _ in pairs])
        all_helices.append(hel)
        conformations.append({
            'id': c.meta.get('id') or f'conf_{ci + 1:03d}',
            'label': c.label,
            'rank': c.meta.get('rank', ci + 1),
            'source': {'file': c.source, 'model': c.model, 'chain': c.chain, 'title': c.title},
            'coords': base64.b64encode(xyz[ci].astype('<f4').tobytes()).decode('ascii'),
            'pairs': [[i, j, t, (i, j) in pk] for i, j, t in pairs],
            'stacks': [[i - 1, j - 1] for i, j in ev['stacks']],
            'dbn': ev['dbn'],
            'input_dbn': c.meta.get('secondary_structure'),
            'secondary_structure_id': c.meta.get('secondary_structure_id'),
            'cluster': c.meta.get('cluster_id'),
            'values': values,
            'validation': validation,
            'explored': c.meta.get('explored'),
        })

    # --- helix families: helices (>= 2 bp) of all conformations linked by shared pairs
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    owner = {}
    nodes = []
    for ci, hel in enumerate(all_helices):
        for hi, h in enumerate(hel):
            if len(h) < 2:
                continue
            node = (ci, hi)
            parent[node] = node
            nodes.append(node)
            for p in h:
                if p in owner:
                    a, b = find(owner[p]), find(node)
                    if a != b:
                        parent[b] = a
                else:
                    owner[p] = node
    fam_pairs = {}
    for node in nodes:
        fam_pairs.setdefault(find(node), set()).update(all_helices[node[0]][node[1]])
    fam_order = sorted(fam_pairs, key=lambda f: (min(p[0] for p in fam_pairs[f]), max(p[1] for p in fam_pairs[f])))
    fam_index = {f: k for k, f in enumerate(fam_order)}
    families = []
    for k, f in enumerate(fam_order):
        ps = sorted(fam_pairs[f])
        families.append({'id': k, 'name': f'H{k + 1}', 'first': ps[0][0], 'last': ps[0][1],
                         'span': [min(p[0] for p in ps), max(p[1] for p in ps)]})
    for ci, hel in enumerate(all_helices):
        hres = [-1] * len(seq)
        for hi, h in enumerate(hel):
            fam = fam_index[find((ci, hi))] if len(h) >= 2 else -2
            for i, j in h:
                hres[i] = hres[j] = fam
        conformations[ci]['helix_of_residue'] = hres

    # --- 2D drawings: one layout per distinct nested secondary structure
    layouts = {}
    for c in conformations:
        key = _nested_dbn(len(seq), c['pairs'])
        if key not in layouts:
            layouts[key] = _layout2d(key)
        c['layout_key'] = key if layouts[key] is not None else None
    layouts = {k: v for k, v in layouts.items() if v is not None}

    # --- bonds (topology from the first conformation)
    bonds = _bonds(xyz[0], elements, atom_res)

    # --- metric catalogue
    metrics = []
    terms = {}
    if engine_mode:
        metrics.append({'key': 'score', 'label': 'Model score', 'unit': 'a.u.', 'group': 'engine',
                        'lower_is_better': True,
                        'terms': [k for k in ENGINE_TERMS if any(k in c['values'] for c in conformations)],
                        'description': 'rnafold3d engine score used for ranking (dimensionless, kT-like). '
                                       'Not a free energy; ranks are not populations.'})
        for k, (lab, desc) in ENGINE_TERMS.items():
            terms[k] = {'label': lab, 'description': desc, 'unit': 'a.u.', 'group': 'engine'}
        if any('ss_dg' in c['values'] for c in conformations):
            metrics.append({'key': 'ss_dg', 'label': 'ΔG of input 2D', 'unit': 'kcal/mol', 'group': 'engine',
                            'lower_is_better': True,
                            'description': 'ViennaRNA free energy of the secondary-structure hypothesis this '
                                           'conformer was built from.'})
    metrics.append({**EVAL_TOTAL, 'group': 'evaluator'})
    for k, t in EVAL_TERMS.items():
        terms[k] = {**t, 'group': 'evaluator'}
        metrics.append({'key': k, 'label': t['label'], 'unit': t['unit'], 'group': 'evaluator',
                        'lower_is_better': True, 'description': t['description']})
    for k, d in DESCRIPTORS.items():
        metrics.append({'key': k, **d, 'group': 'descriptor'})
    metrics.append({'key': 'rmsd_first', 'label': 'RMSD to first', 'unit': 'Å', 'group': 'descriptor',
                    'lower_is_better': None, 'description': "C4' RMSD to the first conformation."})
    metrics.append({'key': 'rmsd_mean', 'label': 'RMSD to mean', 'unit': 'Å', 'group': 'descriptor',
                    'lower_is_better': None, 'description': "C4' RMSD to the ensemble mean structure."})
    for k, m in plugin_meta.items():
        metrics.append({'key': k, 'label': m.get('label', k), 'unit': m.get('unit', ''), 'group': 'plugin',
                        'lower_is_better': m.get('lower_is_better'), 'description': m.get('description', '')})
    if all(c['values'].get('dg2d') is None for c in conformations):
        warnings.add('ViennaRNA is not available: ΔG 2D and the evaluator energy omit the 2D term.')

    digest = hashlib.sha256()
    for c in conformations:
        digest.update(c['coords'].encode())
    doc = {
        'schema': SCHEMA,
        'id': digest.hexdigest()[:12],
        'name': name,
        'created': datetime.datetime.now().isoformat(timespec='seconds'),
        'source': source or {},
        'origin': 'rnafold3d-engine' if engine_run else ('rnafold3d-package' if manifest else 'structures'),
        'sequence': seq,
        'residues': residues_out,
        'atoms': {'names': atom_names, 'elements': elements, 'residue': atom_res},
        'bonds': [v for b in bonds for v in b],
        'conformations': conformations,
        'default_metric': 'score' if engine_run else 'eval_total',
        'metrics': metrics,
        'terms': terms,
        'helix_families': families,
        'layouts2d': layouts,
        'ensemble': {'rmsd_c4': [[round(float(v), 3) for v in row] for row in rmsd],
                     'embedding': embedding, 'rmsf': rmsf,
                     'superposition': "C4' atoms, fitted to the first conformation then to the ensemble mean"},
        'provenance': _provenance(manifest),
        'plugins': [p.id for p in plugins],
        'warnings': sorted(warnings),
    }
    try:                                                # topology-aware ensemble map (Energy tab, Pathway tab)
        from .landscape import layout_doc
        doc['ensemble']['landscape'] = layout_doc(doc)
    except Exception as exc:                            # the map is optional; never block opening
        doc['warnings'].append(f'ensemble map not computed: {type(exc).__name__}: {exc}')
    return doc


def _provenance(manifest):
    if not manifest:
        return None
    inp = manifest.get('input', {})
    return {'software': manifest.get('software'),
            'secondary_structure_mode': inp.get('secondary_structure_mode'),
            'secondary_structures': [{k: v for k, v in s.items() if k != 'pairs'}
                                     for s in inp.get('secondary_structures', [])],
            'parameters': {k: manifest.get('parameters', {}).get(k) for k in
                           ('n_models', 'trajectories', 'trajectories_run', 'seed', 'temperature', 'ionic_conditions',
                            'cluster_rmsd_threshold', 'max_structures')
                           if k in manifest.get('parameters', {})},
            'statistics': manifest.get('statistics'),
            'notes': manifest.get('notes', [])}


def load_document(paths, chain_id=None, name=None, plugins=(), progress=None):
    nm, confs, manifest, warnings, source = collect_paths(paths, chain_id)
    return build_document(name or nm, confs, manifest, warnings, source, plugins, progress)


def load_uploads(files, chain_id=None, name=None, plugins=(), progress=None):
    nm, confs, manifest, warnings, source = collect_uploads(files, chain_id)
    return build_document(name or nm, confs, manifest, warnings, source, plugins, progress)
