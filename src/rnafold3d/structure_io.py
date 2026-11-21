"""Lightweight, dependency-free mmCIF / PDB reading and writing for RNA chains.

Reads the first model only (read_models: every model), keeps the first alternate
conformer of each atom, maps modified nucleotides to their canonical parent base
(declared in the file: ``_chem_comp.mon_nstd_parent_comp_id``,
``_pdbx_struct_mod_residue`` or MODRES; else the table below; else the atoms of
the base), and preserves both label and author numbering. What a model holds
besides RNA (ligands, ions, protein and DNA chains) is kept in a reduced form
in ``meta['others']``.

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
# Checked against the wwPDB chemical component dictionary (2026-11): mirror-image (L-) RNA, 2'-O-methoxyethyl,
# locked and 2'-substituted nucleotides that declare no parent, or a DNA parent, in the files.
FALLBACK_PARENTS.update({
    '0A': 'A', '0C': 'C', '0G': 'G', '0U': 'U',
    'A44': 'A', 'C43': 'C', 'G48': 'G', 'U36': 'U', 'T39': 'U', 'C5L': 'C',
    'A2L': 'A', 'C2L': 'C', 'G2L': 'G', 'U2L': 'U',
    'LKC': 'C', 'UMS': 'U', 'CSL': 'C', 'A5M': 'C', '10C': 'C',
})
L_RNA = ('0A', '0C', '0G', '0U')
# A modified nucleotide of an RNA chain whose declared parent is a DNA nucleotide (2'-fluoro, 2'-methylseleno ...)
DNA_PARENTS = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DU': 'U', 'DT': 'U', 'DI': 'G'}
WATER = frozenset(('HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'TIP3', 'SOL'))
AMINO = frozenset((
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL', 'SEC', 'PYL', 'MSE', 'ASX', 'GLX', 'UNK', 'HID', 'HIE', 'HIP', 'CYX', 'SEP', 'TPO',
    'PTR', 'CSO', 'HYP', 'MLY', 'KCX', 'LLP', 'CME', 'OCS', 'CSD', 'M3L', 'ALY', 'PCA', 'FME'))
# 2' substituents: a modified residue of a PDB-format file that declares a DNA parent is RNA-like when it has one
_SUGAR_2PRIME = frozenset(("O2'", "F2'", 'F', "N2'", "S2'", "SE2'", "C2M", "CM2", "O2*"))
try:                                                     # modified nucleotides written by the engine (0.7.6)
    from .modifications import CCD_PARENTS as _MOD_PARENTS
    for _k, _v in _MOD_PARENTS.items():
        FALLBACK_PARENTS.setdefault(_k, _v)
except ImportError:                                      # pragma: no cover
    pass
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
    """({component: parent base}, {component: base of a declared DNA parent}) of an mmCIF file."""
    parents, dna = dict(FALLBACK_PARENTS), {}
    cc = cif.get('_chem_comp', {})
    mod = cif.get('_pdbx_struct_mod_residue', {})
    for cid, parent in list(zip(cc.get('id', []), cc.get('mon_nstd_parent_comp_id', []))) + \
            list(zip(mod.get('label_comp_id', []), mod.get('parent_comp_id', []))):
        parent = parent.upper()
        if parent in CANONICAL:
            parents[cid] = parent
        elif parent in DNA_PARENTS and cid not in DNA_PARENTS:
            dna[cid] = DNA_PARENTS[parent]
    for b in CANONICAL:
        parents[b] = b
    return parents, dna


def base_from_atoms(names):
    """Parent base of a nucleotide that declares none, from the atom names of its base: a purine with O6 (or S6) is
    read as G, another purine as A, a pyrimidine with N4 as C, another pyrimidine as U. None when it has no sugar
    or no recognisable base."""
    if "C1'" not in names:
        return None
    if 'N9' in names and 'C8' in names and 'C4' in names:
        return 'G' if ('O6' in names or 'S6' in names) else 'A'
    if 'N1' in names and 'C6' in names and 'C5' in names and 'C4' in names:
        return 'C' if 'N4' in names else 'U'
    return None


def new_others():
    """What a model holds besides its RNA chains, in the reduced form the viewer shows as context."""
    return {'groups': {},        # (chain, number, insertion code, component) -> {name: (element, xyz)}: ligands, ions
            'protein': {},       # chain -> {(number, insertion code): (component, xyz of CA)}
            'dna': {},           # chain -> {(number, insertion code): (component, xyz of C4')}
            'water': 0}


def _element_guess(name, comp):
    """Element from a PDB atom name when the element columns are empty."""
    n = name.strip().upper()
    if comp.upper() == n and len(n) <= 2:               # an ion: residue and atom carry the same name
        return n
    n = n.lstrip('0123456789')
    return n[:1] or 'C'


def _note_other(others, kind, chain, num, ins, comp, atom, elem, xyz):
    """File an atom that is not part of an RNA chain. Proteins and DNA keep one atom per residue."""
    if kind == 'water':
        others['water'] += atom == 'O' or atom == 'OW'
    elif kind == 'protein':
        if atom == 'CA':
            others['protein'].setdefault(chain, {}).setdefault((num, ins), (comp, xyz))
    elif kind == 'dna':
        if atom == "C4'":
            others['dna'].setdefault(chain, {}).setdefault((num, ins), (comp, xyz))
    else:
        others['groups'].setdefault((chain, num, ins, comp), {}).setdefault(atom, (elem, xyz))


def _settle_others(others):
    """Residues with a peptide backbone that came in as hetero groups (modified amino acids) join the protein."""
    for key in [k for k, atoms in others['groups'].items() if {'N', 'CA', 'C', 'O'} <= set(atoms) and len(atoms) > 3]:
        chain, num, ins, comp = key
        atoms = others['groups'].pop(key)
        others['protein'].setdefault(chain, {}).setdefault((num, ins), (comp, atoms['CA'][1]))
    return others


def describe_others(others):
    """One phrase on what else a model holds, e.g. 'protein (129 residues), DNA (24 nucleotides) and 2 ligands or
    ions (MG, SAM)'; '' when there is nothing."""
    if not others:
        return ''
    parts = []
    n = sum(len(v) for v in others['protein'].values())
    if n:
        parts.append(f'protein ({n} residue{"s" if n != 1 else ""})')
    n = sum(len(v) for v in others['dna'].values())
    if n:
        parts.append(f'DNA ({n} nucleotide{"s" if n != 1 else ""})')
    if others['groups']:
        comps = sorted({k[3] for k in others['groups']})
        n = len(others['groups'])
        parts.append(f'{n} ligand{"s" if n != 1 else ""} or ion{"s" if n != 1 else ""} '
                     f'({", ".join(comps[:6])}{" …" if len(comps) > 6 else ""})')
    if others['water'] and not parts:
        parts.append(f'{others["water"]} water molecules')
    return ' and '.join([', '.join(parts[:-1]), parts[-1]] if len(parts) > 1 else parts)


def read_mmcif(path, rna_only=True):
    """Return (chains, metadata) for the first model. Chains keyed by label_asym_id."""
    text = Path(path).read_text(errors='replace')
    return _mmcif_model(parse_cif(text), rna_only, None, Path(path).name)


def _mmcif_model_ids(cif):
    site = cif.get('_atom_site') or {}
    ids = []
    for m in site.get('pdbx_PDB_model_num', ['1'] if site else []):
        if m not in '?.' and (not ids or ids[-1] != m) and m not in ids:
            ids.append(m)
    return ids or ['1']


def _mmcif_model(cif, rna_only=True, model_id=None, label='mmCIF'):
    parents, dna_parents = _parents(cif)
    site = cif.get('_atom_site')
    if not site:
        raise ValueError(f'{label}: this mmCIF file has no atom table (_atom_site), so it holds no coordinates. '
                         'Download the structure again, or check that it is a coordinate file and not a '
                         'structure-factor or dictionary file')
    missing = [k for k in ('Cartn_x', 'Cartn_y', 'Cartn_z', 'label_atom_id', 'label_comp_id') if k not in site
               and not (k.startswith('label_') and 'auth' + k[5:] in site)]
    if missing:
        raise ValueError(f'{label}: the atom table of this mmCIF file lacks the column{"s" if len(missing) > 1 else ""} '
                         f'{", ".join("_atom_site." + k for k in missing)}; the file is incomplete or was cut off. '
                         'Download or export it again')
    n_rows = len(site['Cartn_x'])

    def get(k, default='?', alt=None):
        col = site.get(k)
        if col is None and alt:
            col = site.get(alt)
        return col if col is not None else [default] * n_rows
    model = get('pdbx_PDB_model_num', '1')
    first_model = model_id if model_id is not None else next((m for m in model if m not in '?.'), '1')
    seq_by_entity = {}
    ep = cif.get('_entity_poly', {})
    for eid, typ, can in zip(ep.get('entity_id', []), ep.get('type', []),
                             ep.get('pdbx_seq_one_letter_code_can', [''] * len(ep.get('entity_id', [])))):
        seq_by_entity[eid] = (typ, re.sub(r'\s+', '', can).upper())
    comp_col, atom_col = get('label_comp_id', alt='auth_comp_id'), get('label_atom_id', alt='auth_atom_id')
    # nucleotides that declare no parent: the base is read from their atoms (inside RNA chains only, see below)
    unknown = {c for c in set(comp_col) if c not in parents and c not in dna_parents and c not in DNA_PARENTS
               and c not in AMINO and c not in WATER}
    inferred = {}
    if unknown:
        names = {c: set() for c in unknown}
        for c, a in zip(comp_col, atom_col):
            if c in names:
                names[c].add(a.replace('*', "'"))
        inferred = {c: b for c, b in ((c, base_from_atoms(v)) for c, v in names.items()) if b}
    chains = {}
    others = new_others()
    seen_atoms = set()
    row = 0
    try:
        for row, (comp, asym, seq, atom, x, y, z, mdl, eid, aasym, aseq, ins, elem) in enumerate(zip(
                comp_col, get('label_asym_id', alt='auth_asym_id'), get('label_seq_id', alt='auth_seq_id'), atom_col,
                site['Cartn_x'], site['Cartn_y'], site['Cartn_z'], model,
                get('label_entity_id'), get('auth_asym_id', alt='label_asym_id'), get('auth_seq_id', alt='label_seq_id'),
                get('pdbx_PDB_ins_code'), get('type_symbol'))):
            if mdl != first_model or elem in ('H', 'D'):
                continue
            # an entity that _entity_poly does not list is taken as RNA when its residues are numbered (files
            # written by modelling programs often have no entity table)
            typ = seq_by_entity[eid][0] if eid in seq_by_entity else ''
            polymer = seq.lstrip('-').isdigit()
            rna = eid not in seq_by_entity or 'ribonucleotide' in typ.replace('deoxyribonucleotide', '')
            atom = atom.replace('*', "'")
            base = None
            if polymer and (rna or not rna_only):
                base = parents.get(comp)
                if base is None and 'deoxy' not in typ:      # a DNA parent or none declared: RNA chains only
                    base = dna_parents.get(comp) or inferred.get(comp)
            if base is None:
                if comp in WATER:
                    kind = 'water'
                elif 'polypeptide' in typ or (polymer and comp in AMINO):
                    kind = 'protein'
                elif polymer and (comp in DNA_PARENTS or 'deoxyribonucleotide' in typ):
                    kind = 'dna'
                else:
                    kind = 'group'
                if kind in ('water', 'group') or atom in ('CA', "C4'"):
                    num = aseq if aseq not in '?.' else seq
                    _note_other(others, kind, aasym if aasym not in '?.' else asym, num, '' if ins in '?.' else ins,
                                comp, atom, elem.upper() if elem not in '?.' else _element_guess(atom, comp),
                                np.array([float(x), float(y), float(z)]))
                continue
            if atom == 'O1P':
                atom = 'OP1'
            elif atom == 'O2P':
                atom = 'OP2'
            key = (asym, seq, atom)
            if key in seen_atoms:                            # the first alternate conformer of each atom
                continue
            if asym not in chains:
                chains[asym] = Chain(asym, aasym, eid, seq_by_entity.get(eid, ('', ''))[1], {})
            res = chains[asym].residues.get(int(seq))
            if res is None:
                res = chains[asym].residues[int(seq)] = Residue(int(seq), base, comp, {}, aseq, '' if ins in '?.' else ins, aasym)
            elif res.comp_id != comp:                        # two components at one position: the first one
                continue
            seen_atoms.add(key)
            xyz = np.array([float(x), float(y), float(z)])
            if not (abs(xyz) < 1e9).all():                   # NaN or infinite
                raise ValueError('a coordinate is not a number (NaN or infinite)')
            res.atoms[atom] = xyz
    except ValueError as exc:
        raise ValueError(f'{label}: atom {row + 1} of the atom table has no readable coordinates ({exc}); the file '
                         'is damaged or was cut off. Download or export it again') from None
    for ch in chains.values():
        ch.residues = [ch.residues[k] for k in sorted(ch.residues)]
        if not ch.sequence:
            ch.sequence = ''.join(r.base for r in ch.residues)
    meta = {'title': ' '.join(cif.get('_struct', {}).get('title', [''])),
            'method': ' '.join(cif.get('_exptl', {}).get('method', [''])),
            'resolution': None, 'others': _settle_others(others)}
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


def _pdb_header(lines):
    """({component: parent} from MODRES, set of components in SEQRES) of the header lines of a PDB file."""
    modres, seqres = {}, set()
    for line in lines:
        if line.startswith('MODRES'):
            comp, std = line[12:15].strip().upper(), line[24:27].strip().upper()
            if comp and std:
                modres[comp] = std
        elif line.startswith('SEQRES'):
            seqres.update(line[19:70].upper().split())
        elif line.startswith(('ATOM  ', 'HETATM', 'MODEL')):
            break
    return modres, seqres


def read_pdb(path):
    """Read RNA residues from a PDB file (first model). Index = order of appearance per chain."""
    lines = []
    for line in Path(path).read_text(errors='replace').splitlines():
        if line.startswith('ENDMDL'):
            break
        lines.append(line)
    chains, others = _pdb_chains([(k, l) for k, l in enumerate(lines, 1)], *_pdb_header(lines), label=Path(path).name)
    return chains, {'others': others}


def _pdb_model_blocks(text):
    """Split PDB text into per-model ATOM/HETATM blocks: [(model_id, [(line number, line)])]."""
    blocks, cur, mid = [], [], None
    for k, line in enumerate(text.splitlines(), 1):
        if line.startswith('MODEL'):
            if cur:
                blocks.append((mid or str(len(blocks) + 1), cur))
            cur, mid = [], (line[10:14].strip() or line[5:].strip() or None)
        elif line.startswith('ENDMDL'):
            if cur:
                blocks.append((mid or str(len(blocks) + 1), cur))
            cur, mid = [], None
        elif line.startswith(('ATOM  ', 'HETATM')):
            cur.append((k, line))
    if cur:
        blocks.append((mid or str(len(blocks) + 1), cur))
    return blocks


def _pdb_chains(lines, modres=None, seqres=(), label='PDB file'):
    """RNA chains and the other contents of one model. lines: [(line number, line)] or plain lines."""
    modres = modres or {}
    chains = {}
    order = {}
    others = new_others()
    het = set()                    # residues that came from HETATM records
    from_dna = set()               # residues read as RNA through a DNA parent
    unknown = {}                   # SEQRES components without a parent: (chain, key) -> atoms, decided at the end
    for item in lines:
        k, line = item if isinstance(item, tuple) else (0, item)
        if not line.startswith(('ATOM  ', 'HETATM')):
            continue
        comp = line[17:20].strip().upper()
        name = line[12:16].strip().replace('*', "'")
        name = {'O1P': 'OP1', 'O2P': 'OP2'}.get(name, name)
        element = line[76:78].strip().upper() if len(line) >= 78 else ''
        if element in ('H', 'D') or (not element and name.lstrip('0123456789').startswith(('H', 'D'))
                                     and comp.upper() != name.upper()):
            continue
        try:
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        except ValueError:
            where = f'line {k}' if k else f'the record of atom {line[6:11].strip() or "?"}'
            raise ValueError(f'{label}: {where} is not a complete {line[:6].strip()} record (no coordinates in columns '
                             f'31-54): "{line.rstrip()[:60]}". The file is damaged, was cut off, or is not in PDB '
                             'format; export it again or use the mmCIF file') from None
        if not np.isfinite(xyz).all():
            raise ValueError(f'{label}: {"line " + str(k) if k else "an atom record"} has a coordinate that is not a '
                             'number (NaN or infinite); the program that wrote this file did not finish the model')
        ch = line[21].strip() or '_'
        key = (line[22:26].strip(), line[26].strip())
        base = FALLBACK_PARENTS.get(comp, comp if comp in CANONICAL else None)
        dna = False
        if base is None and comp in modres:
            std = modres[comp]
            base = std if std in CANONICAL else DNA_PARENTS.get(std)
            dna = std in DNA_PARENTS
        if base is None:
            if comp in seqres and comp not in DNA_PARENTS and comp not in AMINO and comp not in WATER:
                unknown.setdefault((ch, key, comp), {})[name] = (element or _element_guess(name, comp), xyz)
                if key not in chains.get(ch, {}):            # keep its place in the chain
                    chains.setdefault(ch, {})[key] = None
                    order.setdefault(ch, []).append(key)
                continue
            kind = 'water' if comp in WATER else 'protein' if comp in AMINO else 'dna' if comp in DNA_PARENTS else 'group'
            _note_other(others, kind, ch, key[0], key[1], comp, name, element or _element_guess(name, comp), xyz)
            continue
        chains.setdefault(ch, {})
        order.setdefault(ch, [])
        if key not in chains[ch]:
            chains[ch][key] = Residue(len(order[ch]) + 1, base, comp, {}, key[0], key[1], ch)
            order[ch].append(key)
            if line.startswith('HETATM'):
                het.add((ch, key))
            if dna:
                from_dna.add((ch, key))
        res = chains[ch][key]
        if res is None or res.comp_id != comp:               # two components at one position: the first one
            continue
        if name not in res.atoms:                            # the first alternate conformer of each atom
            res.atoms[name] = xyz
    for (ch, key, comp), atoms in unknown.items():           # no parent declared: the base is read from the atoms
        base = base_from_atoms(set(atoms))
        if base and chains[ch].get(key) is None:
            chains[ch][key] = Residue(0, base, comp, {a: x for a, (_e, x) in atoms.items()}, key[0], key[1], ch)
        elif chains[ch].get(key) is None:
            for a, (e, x) in atoms.items():
                _note_other(others, 'group', ch, key[0], key[1], comp, a, e, x)
    out = {}
    for ch, residues in chains.items():
        rs = []
        for key in order[ch]:
            r = residues[key]
            if r is None:
                continue
            # a hetero group without a sugar (a free base such as ADE) or a DNA-parent residue without a 2'
            # substituent is not part of the RNA chain: it goes to the other contents of the model
            if ((ch, key) in het and "C1'" not in r.atoms and len(r.atoms) > 3) or \
                    ((ch, key) in from_dna and not (_SUGAR_2PRIME & set(r.atoms))):
                for a, x in r.atoms.items():
                    _note_other(others, 'dna' if (ch, key) in from_dna else 'group', ch, key[0], key[1], r.comp_id, a,
                                _element_guess(a, r.comp_id), x)
                continue
            rs.append(r)
        for i, r in enumerate(rs, 1):
            r.index = i
        if rs:
            out[ch] = Chain(ch, ch, '', ''.join(r.base for r in rs), rs)
    return out, _settle_others(others)


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
    NMR ensembles and multi-model PDB files yield one entry per model. Models without RNA are left out when
    another model of the file has RNA (meta['others'] tells what a model holds besides RNA).
    """
    if text is None:
        text = Path(path).read_text(errors='replace')
    name = str(filename or path or '')
    lower = name.lower()
    label = Path(name).name if name else 'structure file'
    is_cif = lower.endswith(('.cif', '.mmcif')) or (not lower.endswith(('.pdb', '.ent')) and _looks_like_cif(text))
    if is_cif:
        cif = parse_cif(text)
        out = []
        for mid in _mmcif_model_ids(cif):
            chains, meta = _mmcif_model(cif, rna_only, mid, label)
            out.append((mid, chains, meta))
    else:
        blocks = _pdb_model_blocks(text)
        all_lines = text.splitlines()
        title = ' '.join(l[10:].strip() for l in all_lines[:200] if l.startswith('TITLE'))
        modres, seqres = _pdb_header(all_lines)
        out = []
        for mid, lines in blocks:
            chains, others = _pdb_chains(lines, modres, seqres, label)
            out.append((mid, chains, {'title': title, 'method': '', 'resolution': None, 'others': others}))
    if any(chains for _m, chains, _meta in out):
        out = [m for m in out if m[1]]
    return out


# ----------------------------------------------------------------------------- writing
def _pdb_atom_name(name, element):
    if len(name) >= 4:
        return name[:4]
    return (' ' + name).ljust(4) if len(element) == 1 else name.ljust(4)


def element_of(name):
    """Element of an atom name of a (modified) nucleotide."""
    n = name.upper()
    if n.startswith('OP'):
        return 'O'
    for two in ('BR', 'CL', 'SE'):
        if n.startswith(two):
            return two
    return n[0] if n[0] in 'CNOPSFI' else 'C'


def write_pdb(path, residues, sequence, chain_id='A', remarks=(), model_number=None, comps=None):
    """residues: list (per sequence position) of {atom_name: xyz}; numbering 1..N. comps: residue names (modified
    nucleotides, written as HETATM with MODRES records); default the sequence letters."""
    comps = list(comps) if comps else list(sequence)
    lines = [f'REMARK 250 {r}'[:80] for r in remarks]
    for i, (base, comp) in enumerate(zip(sequence, comps), 1):
        if comp != base:
            lines.append(f'MODRES XXXX {comp:>3s} {chain_id} {i:4d}  {base:>3s}  MODIFIED NUCLEOTIDE'[:80])
    if model_number is not None:
        lines.append(f'MODEL     {model_number:4d}')
    serial = 1
    for i, (base, comp, atoms) in enumerate(zip(sequence, comps, residues), 1):
        rec = 'ATOM  ' if comp == base else 'HETATM'
        for name, xyz in atoms.items():
            element = element_of(name)
            x, y, z = (float(v) for v in xyz)
            lines.append(f'{rec}{serial:5d} {_pdb_atom_name(name, element)} {comp:>3s} {chain_id}{i:4d}    '
                         f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}')
            serial += 1
    lines.append(f'TER   {serial:5d}      {comps[-1]:>3s} {chain_id}{len(sequence):4d}')
    if model_number is not None:
        lines.append('ENDMDL')
    lines.append('END')
    Path(path).write_text('\n'.join(lines) + '\n')


def write_mmcif(path, residues, sequence, chain_id='A', entry_id='RNAFOLD3D', title='', comps=None):
    """Minimal standards-conformant mmCIF (atom_site + entity_poly) for one RNA chain. comps: residue names of
    modified nucleotides (declared in _chem_comp with their parent); default the sequence letters."""
    comps = list(comps) if comps else list(sequence)
    out = [f'data_{entry_id}', '#', f'_entry.id {entry_id}', '#']
    if title:
        out += [f"_struct.title '{title.replace(chr(39), ' ')}'", '#']
    one = ''.join(b if c == b else f'({c})' for b, c in zip(sequence, comps))
    out += ['_entity_poly.entity_id 1', "_entity_poly.type polyribonucleotide",
            f'_entity_poly.pdbx_seq_one_letter_code {one}', f'_entity_poly.pdbx_seq_one_letter_code_can {sequence}', '#',
            'loop_', '_entity_poly_seq.entity_id', '_entity_poly_seq.num', '_entity_poly_seq.mon_id', '_entity_poly_seq.hetero']
    out += [f'1 {i} {c} n' for i, c in enumerate(comps, 1)]
    modified = sorted({(c, b) for b, c in zip(sequence, comps) if c != b})
    if modified:
        out += ['#', 'loop_', '_chem_comp.id', '_chem_comp.type', '_chem_comp.mon_nstd_parent_comp_id']
        out += [f"{c} 'RNA linking' {b}" for c, b in modified]
    out += ['#', 'loop_'] + [f'_atom_site.{k}' for k in (
        'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
        'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
        'B_iso_or_equiv', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num')]
    serial = 1
    for i, (base, comp, atoms) in enumerate(zip(sequence, comps, residues), 1):
        group = 'ATOM' if comp == base else 'HETATM'
        for name, xyz in atoms.items():
            element = element_of(name)
            q = f'"{name}"' if "'" in name else name
            x, y, z = (float(v) for v in xyz)
            out.append(f'{group} {serial} {element} {q} . {comp} {chain_id} 1 {i} ? {x:.3f} {y:.3f} {z:.3f} 1.00 0.00 {i} {comp} {chain_id} {q} 1')
            serial += 1
    out.append('#')
    Path(path).write_text('\n'.join(out) + '\n')
