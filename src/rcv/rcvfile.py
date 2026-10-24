"""The .rcv file: every conformation of one RNA in one compact, self-describing file.

An .rcv file is a ZIP archive (open it with any unzip tool) containing

    mimetype         'application/vnd.rcv+zip' (stored uncompressed, first entry, for file-type detection)
    rcv.json         manifest + the ensemble document without coordinates: sequence, residues, atom names,
                     per-conformation base pairs, dot-brackets, energies and energy terms, validation reports,
                     helix families, ensemble statistics (RMSD matrix, similarity map, RMSF), 2D layouts,
                     provenance (engine version, parameters, seeds, source-file checksums) and warnings
    coords.bin       all coordinates, 0.001 Å integers (see encode_coords), SHA-256 checked on reading
    bonds.bin        covalent bonds as little-endian uint32 atom-index pairs
    session.json     optional viewer state: camera, colours, representation, selections, measurements,
                     scenes, notes, command history
    thumbnail.jpg    optional preview image (or thumbnail.png)
    methods.txt      optional methods paragraph describing how the ensemble was obtained
    provenance/      optional verbatim engine manifest (manifest.json)
    README.txt       this description

Coordinates ('rcv-coords/1'): int32 little-endian milli-ångström values, written per conformation as three
planes (all x, then all y, then all z). Conformation 0 is stored as differences along the atom order
(v[i] - v[i-1]); every other conformation as differences to conformation 0 (same atom, same axis). The int32
stream is byte-shuffled (all first bytes, then all second bytes, ...) and deflated by the ZIP container.
Precision equals the PDB format (0.001 Å); files are typically 10-15x smaller than the same models as PDB.
"""
from __future__ import annotations

import base64
import copy
import datetime
import hashlib
import io
import json
import re
import zipfile
from pathlib import Path

import numpy as np

FORMAT = 'rcv/1'
MIME = 'application/vnd.rcv+zip'
COORDS_ENCODING = 'rcv-coords/1'
MAX_UNPACKED = 1 << 30            # refuse archives that inflate beyond 1 GiB (damaged or hostile files)
MAX_MEMBERS = 2000

README = __doc__


# ----------------------------------------------------------------------------- coordinates
def _shuffle(a):
    return np.ascontiguousarray(a.astype('<i4').view(np.uint8).reshape(-1, 4).T).tobytes()


def _unshuffle(data, n):
    b = np.frombuffer(data, np.uint8)
    if b.size != 4 * n:
        raise ValueError('coordinate block has the wrong size')
    return np.ascontiguousarray(b.reshape(4, n).T).view('<i4').reshape(-1)


def encode_coords(xyz):
    """xyz: (n_conf, n_atoms, 3) in Å -> bytes (see module docstring)."""
    q = np.rint(np.asarray(xyz, np.float64) * 1000.0)
    if np.abs(q).max(initial=0) > 2 ** 31 - 1:
        raise ValueError('coordinates out of range')
    q = q.astype(np.int64).transpose(0, 2, 1)           # (conf, axis, atom)
    d = q.copy()
    d[0, :, 1:] = q[0, :, 1:] - q[0, :, :-1]
    d[1:] = q[1:] - q[0]
    if np.abs(d).max(initial=0) > 2 ** 31 - 1:
        raise ValueError('coordinates out of range')
    return _shuffle(d.reshape(-1))


def decode_coords(data, n_conf, n_atoms):
    d = _unshuffle(data, n_conf * n_atoms * 3).astype(np.int64).reshape(n_conf, 3, n_atoms)
    q = d.copy()
    q[0] = np.cumsum(d[0], axis=1)
    q[1:] = d[1:] + q[0]
    return q.transpose(0, 2, 1).astype(np.float64) / 1000.0


def doc_coords(doc):
    return np.stack([np.frombuffer(base64.b64decode(c['coords']), '<f4').reshape(-1, 3) for c in doc['conformations']])


# ----------------------------------------------------------------------------- packing
def pdb_size_estimate(n_conf, n_atoms):
    """Bytes of the same ensemble as PDB files (81-byte ATOM records plus MODEL/TER/END lines)."""
    return n_conf * (n_atoms * 81 + 3 * 81)


def pack(doc, session=None, thumbnail=None, methods=None, provenance_manifest=None, app_version=''):
    """Ensemble document (+ optional viewer session, PNG thumbnail bytes, methods text) -> .rcv bytes."""
    xyz = doc_coords(doc)
    n_conf, n_atoms = xyz.shape[0], xyz.shape[1]
    coords = encode_coords(xyz)
    bonds = np.asarray(doc.get('bonds', []), '<u4').tobytes()
    meta = {k: v for k, v in doc.items() if k not in ('conformations', 'bonds', 'session')}
    meta['conformations'] = [{k: v for k, v in c.items() if k != 'coords'} for c in doc['conformations']]
    manifest = {
        'format': FORMAT, 'mimetype': MIME, 'app': 'RCV (RNA Conformation Viewer)', 'app_version': app_version,
        'saved': datetime.datetime.now().astimezone().isoformat(timespec='seconds'),
        'name': doc.get('name'), 'sequence': doc.get('sequence'), 'length': len(doc.get('sequence', '')),
        'n_conformations': n_conf, 'n_atoms': n_atoms,
        'coordinates': {'file': 'coords.bin', 'encoding': COORDS_ENCODING, 'unit': 'angstrom', 'scale': 1000,
                        'shape': [n_conf, n_atoms, 3], 'sha256': hashlib.sha256(coords).hexdigest()},
        'bonds': {'file': 'bonds.bin', 'dtype': 'uint32-le', 'count': len(doc.get('bonds', [])) // 2},
        'contents': ['rcv.json', 'coords.bin', 'bonds.bin'],
        'document': meta,
    }
    extra = []
    if session:
        extra.append(('session.json', json.dumps(session, separators=(',', ':'), allow_nan=False).encode('utf-8')))
    if thumbnail:
        thumbnail = bytes(thumbnail)
        extra.append(('thumbnail.jpg' if thumbnail[:2] == b'\xff\xd8' else 'thumbnail.png', thumbnail))
    if methods:
        extra.append(('methods.txt', methods.encode('utf-8')))
    if provenance_manifest:
        extra.append(('provenance/manifest.json', json.dumps(provenance_manifest, indent=1).encode('utf-8')))
    manifest['contents'] += [n for n, _ in extra]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w') as z:
        z.writestr(zipfile.ZipInfo('mimetype', (1980, 1, 1, 0, 0, 0)), MIME, compress_type=zipfile.ZIP_STORED)
        def put(name, data):
            z.writestr(name, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
        put('rcv.json', json.dumps(manifest, separators=(',', ':'), allow_nan=False).encode('utf-8'))
        put('coords.bin', coords)
        put('bonds.bin', bonds)
        for name, data in extra:
            put(name, data)
        put('README.txt', README.encode('utf-8'))
    return buf.getvalue()


def write(path, doc, **kw):
    data = pack(doc, **kw)
    path = Path(path)
    tmp = path.with_name(path.name + '.part')
    tmp.write_bytes(data)
    tmp.replace(path)
    return len(data)


# ----------------------------------------------------------------------------- unpacking
def is_rcv(data_or_path):
    try:
        if isinstance(data_or_path, (bytes, bytearray)):
            head = bytes(data_or_path[:64])
        else:
            with open(data_or_path, 'rb') as f:
                head = f.read(64)
    except OSError:
        return False
    return head.startswith(b'PK\x03\x04') and MIME.encode() in head


def unpack(data):
    """.rcv bytes -> (document, session or None, info dict). Raises ValueError for damaged or foreign files."""
    try:
        z = zipfile.ZipFile(io.BytesIO(data))
    except zipfile.BadZipFile:
        raise ValueError('not an .rcv file (not a ZIP archive)')
    with z:
        infos = z.infolist()
        if len(infos) > MAX_MEMBERS or sum(i.file_size for i in infos) > MAX_UNPACKED:
            raise ValueError('the .rcv file is too large or damaged')
        names = {i.filename for i in infos}
        if 'rcv.json' not in names or 'coords.bin' not in names:
            raise ValueError('not an .rcv file (rcv.json or coords.bin missing)')
        manifest = json.loads(z.read('rcv.json').decode('utf-8'))
        fmt = str(manifest.get('format', ''))
        if not fmt.startswith('rcv/'):
            raise ValueError('not an .rcv file')
        if fmt.split('/', 1)[1].split('.')[0] != '1':
            raise ValueError(f'this .rcv file uses format {fmt}; this version of RCV reads rcv/1 (please update RCV)')
        cinfo = manifest['coordinates']
        if cinfo.get('encoding') != COORDS_ENCODING:
            raise ValueError(f'unknown coordinate encoding {cinfo.get("encoding")}')
        raw = z.read('coords.bin')
        if cinfo.get('sha256') and hashlib.sha256(raw).hexdigest() != cinfo['sha256']:
            raise ValueError('the .rcv file is damaged (coordinate checksum mismatch)')
        n_conf, n_atoms = int(manifest['n_conformations']), int(manifest['n_atoms'])
        xyz = decode_coords(raw, n_conf, n_atoms)
        bonds = np.frombuffer(z.read('bonds.bin'), '<u4').astype(int).tolist() if 'bonds.bin' in names else []
        session = json.loads(z.read('session.json').decode('utf-8')) if 'session.json' in names else None
        thumb_name = next((n for n in ('thumbnail.jpg', 'thumbnail.png') if n in names), None)
        thumb = z.read(thumb_name) if thumb_name else None
        methods = z.read('methods.txt').decode('utf-8') if 'methods.txt' in names else None
        prov = json.loads(z.read('provenance/manifest.json')) if 'provenance/manifest.json' in names else None
    doc = copy.deepcopy(manifest['document'])
    if len(doc.get('conformations', [])) != n_conf or len(doc.get('atoms', {}).get('names', [])) != n_atoms:
        raise ValueError('the .rcv file is damaged (conformation or atom count mismatch)')
    for c, x in zip(doc['conformations'], xyz):
        c['coords'] = base64.b64encode(x.astype('<f4').tobytes()).decode('ascii')
    doc['bonds'] = bonds
    info = {'format': fmt, 'saved': manifest.get('saved'), 'app_version': manifest.get('app_version'),
            'bytes': len(data), 'pdb_bytes_estimate': pdb_size_estimate(n_conf, n_atoms),
            'thumbnail': base64.b64encode(thumb).decode('ascii') if thumb else None,
            'methods': methods, 'provenance_manifest': prov, 'contents': sorted(names)}
    return doc, session, info


def read(path):
    return unpack(Path(path).read_bytes())


# ----------------------------------------------------------------------------- back to structures
def conformations_from_doc(doc):
    """Rebuild data.Conformation objects from a document, so an ensemble opened from an .rcv file can be extended
    (pull exploration) exactly like one read from structure files. Engine scores are restored into the metadata."""
    from rnafold3d.structure_io import Residue
    from .data import ENGINE_TERMS, Conformation
    xyz = doc_coords(doc).astype(np.float64)
    names = doc['atoms']['names']
    out = []
    for c, x in zip(doc['conformations'], xyz):
        residues = []
        for ri, r in enumerate(doc['residues']):
            s = r['first_atom']
            atoms = {names[s + q]: x[s + q].copy() for q in range(r['n_atoms'])}
            m = re.match(r'(-?\d+)(.*)', str(r.get('num', ri + 1)))
            seq_id, ins = (m.group(1), m.group(2)[:1]) if m else (str(ri + 1), '')
            residues.append(Residue(index=ri, base=r['name'], comp_id=r.get('comp') or r['name'], atoms=atoms,
                                    auth_seq_id=seq_id, ins_code=ins, auth_asym_id=r.get('chain') or 'A'))
        meta = {'id': c.get('id'), 'rank': c.get('rank'), 'cluster_id': c.get('cluster'),
                'secondary_structure': c.get('input_dbn'), 'secondary_structure_id': c.get('secondary_structure_id'),
                'validation': c.get('validation'), 'explored': c.get('explored')}
        v = c.get('values', {})
        if v.get('score') is not None:
            terms = {k: v.get(k, 0.0) for k in ENGINE_TERMS if k not in ('loop_closure', 'ss_offset')}
            terms['total'] = (v.get('score_3d') or 0.0) - (v.get('loop_closure') or 0.0)
            meta.update(total_score=v['score'], score_3d=v.get('score_3d'), score_terms=terms,
                        ss_free_energy_kcal_mol=v.get('ss_dg'), ss_offset_kT=v.get('ss_offset', 0.0))
        meta = {k: val for k, val in meta.items() if val is not None}
        src = c.get('source') or {}
        out.append(Conformation(c.get('label') or c.get('id'), src.get('file') or 'rcv', src.get('model'),
                                src.get('chain') or 'A', residues, meta, src.get('title') or ''))
    return out


def pdb_text(doc, k=None):
    """PDB text of conformation k (or all conformations as models) of a document."""
    xyz = doc_coords(doc)
    names, elems = doc['atoms']['names'], doc['atoms']['elements']
    ks = range(len(xyz)) if k is None else [k]
    lines = [f'REMARK 250 RCV EXPORT OF {doc.get("name", "")}'[:80]]
    for m, ci in enumerate(ks, 1):
        if k is None:
            lines.append(f'MODEL     {m:4d}')
        serial = 1
        for r in doc['residues']:
            num = str(r.get('num', ''))
            seq = ''.join(ch for ch in num if ch.isdigit() or ch == '-') or '1'
            ins = num[len(seq):][:1] or ' '
            for q in range(r['n_atoms']):
                a = r['first_atom'] + q
                nm, el = names[a], elems[a]
                nm4 = nm[:4] if len(nm) >= 4 else ((' ' + nm) if len(el) == 1 else nm).ljust(4)
                x, y, z = xyz[ci, a]
                lines.append(f'ATOM  {serial % 100000:5d} {nm4} {(r.get("comp") or r["name"])[:3]:>3s} '
                             f'{(r.get("chain") or "A")[:1]}{int(seq):4d}{ins}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00'
                             f'          {el:>2s}')
                serial += 1
        lines.append('TER')
        if k is None:
            lines.append('ENDMDL')
    lines.append('END')
    return '\n'.join(lines) + '\n'
