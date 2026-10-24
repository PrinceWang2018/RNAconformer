"""Lightweight, dependency-free mmCIF / PDB reading and writing for RNA chains.

Reads the first model only, keeps the first alternate location (blank, 'A' or
'1'), maps modified nucleotides to their canonical parent base when the file
declares the parent (``_chem_comp.mon_nstd_parent_comp_id`` or
``_pdbx_struct_mod_residue``), and preserves both label and author numbering.

Residue ``index`` is the 1-based position in the polymer sequence
(``label_seq_id`` for mmCIF). For PDB files without SEQRES alignment the index is
the observed order, which is only valid for gap-free chains.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np

CANONICAL = ('A', 'C', 'G', 'U')
# Frequent modified nucleotides when the file omits parent annotations.
FALLBACK_PARENTS = {
    'PSU': 'U', 'H2U': 'U', '5MU': 'U', '4SU': 'U', 'OMU': 'U', 'UR3': 'U',
    '5MC': 'C', 'OMC': 'C', 'CBR': 'C', 'CCC': 'C',
    '1MA': 'A', 'MA6': 'A', '6IA': 'A', 'A2M': 'A', 'AET': 'A',
    '2MG': 'G', 'M2G': 'G', '7MG': 'G', 'OMG': 'G', 'YYG': 'G', '1MG': 'G', 'GTP': 'G', 'GDP': 'G',
    'RA': 'A', 'RC': 'C', 'RG': 'G', 'RU': 'U', 'ADE': 'A', 'CYT': 'C', 'GUA': 'G', 'URA': 'U', 'URI': 'U',
}
_TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'(?=\s|$)|"(?:[^"]|"(?=\S))*"(?=\s|$)|\S+""")


@dataclass
class Residue:
    index: int
    base: str
    comp_id: str
    atoms: dict = field(default_factory=dict)
    auth_seq_id: str = ''
    ins_code: str = ''
    auth_asym_id: str = ''

    @property
    def modified(self) -> bool:
        return self.comp_id not in CANONICAL

    def xyz(self, name):
        return self.atoms[name]


@dataclass
class Chain:
    label_asym_id: str
    auth_asym_id: str
    entity_id: str
    sequence: str            # canonical one-letter sequence of the full polymer (may include unobserved)
    residues: list           # observed residues, ordered by index

    def residue_map(self):
        return {r.index: r for r in self.residues}


def _strip(v):
    if len(v) >= 2 and v[0] == v[-1] and v[0] in '\'"':
        return v[1:-1]
    return v


def _tokenize(text):
    """Yield (token, is_text_field) pairs; handles quotes, comments and ';' fields."""
    lines = text.splitlines()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(';'):
            buf = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(';'):
                buf.append(lines[i])
                i += 1
            i += 1
            yield '\n'.join(buf).strip(), True
            continue
        for tok in _TOKEN.findall(line):
            if tok.startswith('#'):
                break
            yield tok, False
        i += 1


def parse_cif(text):
    """Return {category: {item: [values]}} for loops and single items (first data block)."""
    data = {}
    toks = list(_tokenize(text))
    k = 0
    n = len(toks)
    reserved = ('loop_', 'data_', 'save_', 'global_', 'stop_')
    while k < n:
        tok, is_text = toks[k]
        if not is_text and tok == 'loop_':
            k += 1
            names = []
            while k < n and not toks[k][1] and toks[k][0].startswith('_'):
                names.append(toks[k][0])
                k += 1
            values = []
            while k < n:
                t, txt = toks[k]
                if not txt and (t.startswith('_') or t.startswith(reserved)):
                    break
                values.append(t if txt else _strip(t))
                k += 1
            if names:
                width = len(names)
                for col, name in enumerate(names):
                    cat, _, item = name.partition('.')
                    data.setdefault(cat, {})[item] = values[col::width]
            continue
        if not is_text and tok.startswith('_') and k + 1 < n:
            value, txt = toks[k + 1]
            cat, _, item = tok.partition('.')
            data.setdefault(cat, {}).setdefault(item, []).append(value if txt else _strip(value))
            k += 2
            continue
        k += 1
    return data


def _parents(cif):
    parents = dict(FALLBACK_PARENTS)
    cc = cif.get('_chem_comp', {})
    for cid, parent in zip(cc.get('id', []), cc.get('mon_nstd_parent_comp_id', [])):
        if parent not in ('?', '.', '') and parent.upper() in CANONICAL:
            parents[cid] = parent.upper()
    mod = cif.get('_pdbx_struct_mod_residue', {})
    for cid, parent in zip(mod.get('label_comp_id', []), mod.get('parent_comp_id', [])):
        if parent.upper() in CANONICAL:
            parents[cid] = parent.upper()
    for b in CANONICAL:
        parents[b] = b
    return parents


def read_mmcif(path, rna_only=True):
    """Return (chains, metadata) for the first model. Chains keyed by label_asym_id."""
    text = Path(path).read_text(errors='replace')
    return _mmcif_model(parse_cif(text), rna_only, None, str(path))


def _mmcif_model_ids(cif):
    site = cif.get('_atom_site') or {}
    ids = []
    for m in site.get('pdbx_PDB_model_num', ['1'] if site else []):
        if m not in '?.' and (not ids or ids[-1] != m) and m not in ids:
            ids.append(m)
    return ids or ['1']


def _mmcif_model(cif, rna_only=True, model_id=None, label='mmCIF'):
    parents = _parents(cif)
    site = cif.get('_atom_site')
    if not site:
        raise ValueError(f'{label}: no _atom_site records')
    get = lambda k, default='?': site.get(k, [default] * len(site['Cartn_x']))
    model = get('pdbx_PDB_model_num', '1')
    first_model = model_id if model_id is not None else next((m for m in model if m not in '?.'), '1')
    seq_by_entity = {}
    ep = cif.get('_entity_poly', {})
    for eid, typ, can in zip(ep.get('entity_id', []), ep.get('type', []), ep.get('pdbx_seq_one_letter_code_can', [])):
        seq_by_entity[eid] = (typ, re.sub(r'\s+', '', can).upper())
    chains = {}
    seen_atoms = set()
    for comp, asym, seq, atom, alt, x, y, z, mdl, eid, aasym, aseq, ins, elem in zip(
            get('label_comp_id'), get('label_asym_id'), get('label_seq_id'), get('label_atom_id'),
            get('label_alt_id', '.'), site['Cartn_x'], site['Cartn_y'], site['Cartn_z'], model,
            get('label_entity_id'), get('auth_asym_id'), get('auth_seq_id'), get('pdbx_PDB_ins_code'),
            get('type_symbol')):
        if mdl != first_model or alt not in ('.', '?', 'A', '1') or not seq.lstrip('-').isdigit():
            continue
        if elem in ('H', 'D'):
            continue
        base = parents.get(comp)
        if base is None:
            continue
        typ = seq_by_entity.get(eid, ('polyribonucleotide', ''))[0]
        if rna_only and 'ribonucleotide' not in typ:
            continue
        atom = atom.replace('*', "'")
        if atom == 'O1P':
            atom = 'OP1'
        elif atom == 'O2P':
            atom = 'OP2'
        key = (asym, seq, atom)
        if key in seen_atoms:
            continue
        seen_atoms.add(key)
        if asym not in chains:
            chains[asym] = Chain(asym, aasym, eid, seq_by_entity.get(eid, ('', ''))[1], {})
        res = chains[asym].residues.setdefault(int(seq), Residue(int(seq), base, comp, {}, aseq, '' if ins in '?.' else ins, aasym))
        res.atoms[atom] = np.array([float(x), float(y), float(z)])
    for ch in chains.values():
        ch.residues = [ch.residues[k] for k in sorted(ch.residues)]
        if not ch.sequence:
            ch.sequence = ''.join(r.base for r in ch.residues)
    meta = {'title': ' '.join(cif.get('_struct', {}).get('title', [''])),
            'method': ' '.join(cif.get('_exptl', {}).get('method', [''])),
            'resolution': None}
    for cat, item in (('_refine', 'ls_d_res_high'), ('_reflns', 'd_resolution_high'), ('_em_3d_reconstruction', 'resolution')):
        vals = cif.get(cat, {}).get(item, [])
        for v in vals:
            try:
                meta['resolution'] = float(v)
                break
            except ValueError:
                continue
        if meta['resolution'] is not None:
            break
    return chains, meta


def read_pdb(path):
    """Read RNA residues from a PDB file (first model). Index = order of appearance per chain."""
    lines = []
    for line in Path(path).read_text(errors='replace').splitlines():
        if line.startswith('ENDMDL'):
            break
        lines.append(line)
    return _pdb_chains(lines)


def _pdb_model_blocks(text):
    """Split PDB text into per-model ATOM/HETATM line blocks: [(model_id, lines)]."""
    blocks, cur, mid = [], [], None
    for line in text.splitlines():
        if line.startswith('MODEL'):
            if cur:
                blocks.append((mid or str(len(blocks) + 1), cur))
            cur, mid = [], (line[10:14].strip() or line[5:].strip() or None)
        elif line.startswith('ENDMDL'):
            if cur:
                blocks.append((mid or str(len(blocks) + 1), cur))
            cur, mid = [], None
        elif line.startswith(('ATOM  ', 'HETATM')):
            cur.append(line)
    if cur:
        blocks.append((mid or str(len(blocks) + 1), cur))
    return blocks


def _pdb_chains(lines):
    chains = {}
    order = {}
    for line in lines:
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        alt = line[16].strip()
        if alt not in ('', 'A', '1'):
            continue
        comp = line[17:20].strip().upper()
        base = FALLBACK_PARENTS.get(comp, comp if comp in CANONICAL else None)
        if base is None:
            continue
        name = line[12:16].strip().replace('*', "'")
        name = {'O1P': 'OP1', 'O2P': 'OP2'}.get(name, name)
        element = line[76:78].strip().upper() if len(line) >= 78 else ''
        if element in ('H', 'D') or (not element and name.startswith('H')):
            continue
        ch = line[21].strip() or '_'
        key = (line[22:26].strip(), line[26].strip())
        chains.setdefault(ch, {})
        order.setdefault(ch, [])
        if key not in chains[ch]:
            chains[ch][key] = Residue(len(order[ch]) + 1, base, comp, {}, key[0], key[1], ch)
            order[ch].append(key)
        res = chains[ch][key]
        if name not in res.atoms:
            res.atoms[name] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    out = {}
    for ch, residues in chains.items():
        rs = [residues[k] for k in order[ch]]
        out[ch] = Chain(ch, ch, '', ''.join(r.base for r in rs), rs)
    return out, {}


def read_structure(path):
    p = str(path).lower()
    if p.endswith(('.cif', '.mmcif', '.cif.gz')):
        return read_mmcif(path)
    return read_pdb(path)


def _looks_like_cif(text):
    head = text.lstrip()[:2000]
    return head.startswith('data_') or '_atom_site.' in head or '\nloop_' in head


def read_models(path=None, text=None, filename=None, rna_only=True):
    """Read every model of a PDB or mmCIF file (path or in-memory text).

    Returns a list of (model_id, chains, meta) with the same chain structure as read_structure.
    NMR ensembles and multi-model PDB files yield one entry per model.
    """
    if text is None:
        text = Path(path).read_text(errors='replace')
    name = str(filename or path or '')
    lower = name.lower()
    is_cif = lower.endswith(('.cif', '.mmcif')) or (not lower.endswith(('.pdb', '.ent')) and _looks_like_cif(text))
    if is_cif:
        cif = parse_cif(text)
        out = []
        for mid in _mmcif_model_ids(cif):
            chains, meta = _mmcif_model(cif, rna_only, mid, name or 'mmCIF')
            out.append((mid, chains, meta))
        return out
    blocks = _pdb_model_blocks(text)
    title = ' '.join(l[10:].strip() for l in text.splitlines()[:200] if l.startswith('TITLE'))
    meta = {'title': title, 'method': '', 'resolution': None}
    return [(mid, *(_pdb_chains(lines)[0:1]), dict(meta)) for mid, lines in blocks]


# ----------------------------------------------------------------------------- writing
def _pdb_atom_name(name, element):
    if len(name) >= 4:
        return name[:4]
    return (' ' + name).ljust(4) if len(element) == 1 else name.ljust(4)


def write_pdb(path, residues, sequence, chain_id='A', remarks=(), model_number=None):
    """residues: list (per sequence position) of {atom_name: xyz}; numbering 1..N."""
    lines = [f'REMARK 250 {r}'[:80] for r in remarks]
    if model_number is not None:
        lines.append(f'MODEL     {model_number:4d}')
    serial = 1
    for i, (base, atoms) in enumerate(zip(sequence, residues), 1):
        for name, xyz in atoms.items():
            element = 'P' if name.startswith('P') and name in ('P',) else name[0]
            if name.startswith('OP'):
                element = 'O'
            x, y, z = (float(v) for v in xyz)
            lines.append(f'ATOM  {serial:5d} {_pdb_atom_name(name, element)} {base:>3s} {chain_id}{i:4d}    '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}')
            serial += 1
    lines.append(f'TER   {serial:5d}      {sequence[-1]:>3s} {chain_id}{len(sequence):4d}')
    if model_number is not None:
        lines.append('ENDMDL')
    lines.append('END')
    Path(path).write_text('\n'.join(lines) + '\n')


def write_mmcif(path, residues, sequence, chain_id='A', entry_id='RNAFOLD3D', title=''):
    """Minimal standards-conformant mmCIF (atom_site + entity_poly) for one RNA chain."""
    out = [f'data_{entry_id}', '#', f'_entry.id {entry_id}', '#']
    if title:
        out += [f"_struct.title '{title.replace(chr(39), ' ')}'", '#']
    out += ['_entity_poly.entity_id 1', "_entity_poly.type polyribonucleotide",
            f'_entity_poly.pdbx_seq_one_letter_code_can {sequence}', '#',
            'loop_', '_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id', '_entity_poly_seq.hetero']
    out += [f'1 {i} {b} n' for i, b in enumerate(sequence, 1)]
    out += ['#', 'loop_'] + [f'_atom_site.{k}' for k in (
        'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
        'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
        'B_iso_or_equiv', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num')]
    serial = 1
    for i, (base, atoms) in enumerate(zip(sequence, residues), 1):
        for name, xyz in atoms.items():
            element = 'O' if name.startswith('OP') else name[0]
            q = f'"{name}"' if "'" in name else name
            x, y, z = (float(v) for v in xyz)
            out.append(f'ATOM {serial} {element} {q} . {base} {chain_id} 1 {i} ? {x:.3f} {y:.3f} {z:.3f} 1.00 0.00 {i} {base} {chain_id} {q} 1')
            serial += 1
    out.append('#')
    Path(path).write_text('\n'.join(out) + '\n')
