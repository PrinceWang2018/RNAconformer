"""Sequence and 2D structure of the conformations of an ensemble document as text in the formats other programs
read (the same text as Export > Secondary structure in the window, static/js/export/secondary.js):

    dbn     dot-bracket, the ViennaRNA / FASTA-like form (RNAfold, RNAplot, forna, VARNA, RNAstructure dot2ct ...)
    ct      connect table of mfold and RNAstructure
    bpseq   BPSEQ

The 2D structure of a conformation is the one RCV shows for it: the base pairs found in its 3D coordinates.
"""
from __future__ import annotations

FORMATS = ('dbn', 'ct', 'bpseq')
_OPEN, _CLOSE = '([{<', ')]}>'


def partners(dbn):
    """Partner of every position (0-based, -1 = unpaired) of a dot-bracket with ( ) [ ] { } < >."""
    out = [-1] * len(dbn)
    stacks = [[], [], [], []]
    for k, ch in enumerate(dbn):
        if ch in _OPEN:
            stacks[_OPEN.index(ch)].append(k)
        elif ch in _CLOSE and stacks[_CLOSE.index(ch)]:
            i = stacks[_CLOSE.index(ch)].pop()
            out[i], out[k] = k, i
    return out


def _energy(e):
    return f'{e:.2f}' if isinstance(e, (int, float)) and not isinstance(e, bool) and e == e and abs(e) != float('inf') else None


def _whole(v):
    try:
        return int(str(v).strip()) if str(v).strip().lstrip('-').isdigit() else None
    except ValueError:
        return None


def _title(r):
    return ' '.join(' '.join(str(x) for x in (r.get('name'), r.get('id')) if x).split()) or 'RNA'


def short_id(doc, k):
    """The name the window shows for a conformation: imported models keep their label, predictions their id."""
    c = doc['conformations'][k]
    label = c.get('label')
    if doc.get('origin') != 'rnafold3d-engine' and label and len(label) <= 28:
        return label
    return c['id']


def records(doc, which='all', current=0):
    """[{'name', 'id', 'seq', 'dbn', 'energy', 'numbers', 'also'}]. which: 'current', 'all', or 'distinct' (every
    different 2D once, with the ids of the conformations that share it in 'also')."""
    numbers = [_whole(r.get('num')) for r in doc.get('residues') or []]       # numbering of the file, when it is plain
    plain = bool(numbers) and all(v is not None and (i == 0 or v > numbers[i - 1]) for i, v in enumerate(numbers))

    def rec(k):
        vals = doc['conformations'][k].get('values') or {}
        e = vals.get('ss_dg')
        return {'name': doc.get('name'), 'id': short_id(doc, k), 'seq': doc['sequence'], 'dbn': doc['conformations'][k]['dbn'],
                'energy': e if isinstance(e, (int, float)) else vals.get('dg2d'), 'numbers': numbers if plain else None}
    if which == 'current':
        if not 0 <= current < len(doc['conformations']):
            raise ValueError(f'conformation {current + 1} does not exist (the ensemble has {len(doc["conformations"])})')
        return [rec(current)]
    out = [rec(k) for k in range(len(doc['conformations']))]
    if which != 'distinct':
        return out
    seen = {}
    for r in out:
        if r['dbn'] in seen:
            seen[r['dbn']]['also'].append(r['id'])
        else:
            seen[r['dbn']] = dict(r, also=[])
    return list(seen.values())


def to_dot_bracket(recs):
    out = []
    for r in recs:
        e = _energy(r.get('energy'))
        extra = ' '.join(x for x in (f'dG={e} kcal/mol' if e is not None else '',
                                     'also ' + ','.join(r['also']) if r.get('also') else '') if x)
        out.append(f'>{_title(r)}{" " + extra if extra else ""}\n{r["seq"]}\n{r["dbn"]}\n')
    return ''.join(out)


def to_ct(recs):
    out = []
    for r in recs:
        p, n, e = partners(r['dbn']), len(r['seq']), _energy(r.get('energy'))
        out.append(f'{n:5d}  {"ENERGY = " + e + "  " if e is not None else ""}{_title(r)}')
        nums = r.get('numbers')
        for i, b in enumerate(r['seq']):
            out.append(f'{i + 1:5d} {b} {i:7d} {i + 2 if i + 1 < n else 0:4d} {p[i] + 1:4d} {nums[i] if nums else i + 1:4d}')
    return '\n'.join(out) + '\n' if out else ''


def to_bpseq(recs):
    blocks = []
    for r in recs:
        p, e = partners(r['dbn']), _energy(r.get('energy'))
        lines = [f'# {_title(r)}{" dG=" + e + " kcal/mol" if e is not None else ""}']
        lines += [f'{i + 1} {b} {p[i] + 1}' for i, b in enumerate(r['seq'])]
        blocks.append('\n'.join(lines) + '\n')
    return '\n'.join(blocks)


def text(doc, fmt='dbn', which='all', current=0):
    if fmt not in FORMATS:
        raise ValueError(f'unknown format {fmt!r}: choose dbn (dot-bracket), ct or bpseq')
    return {'dbn': to_dot_bracket, 'ct': to_ct, 'bpseq': to_bpseq}[fmt](records(doc, which, current))
