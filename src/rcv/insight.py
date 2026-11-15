"""Advanced analyses of an ensemble (0.7.4). They appear as closed sections inside the existing tabs:

Switch tab    experiments()        experiments that tell state A from state B: disrupt-and-rescue double mutants,
                                   nucleotides whose chemical-probing signal differs, label pairs for FRET / PELDOR
              variants()           what a list of variants (G34U, A12C+G30C ...) does to the ensemble
Pathway tab   cotranscriptional()  the structure of the growing transcript, length by length
Details tab   motifs()             named loops, junctions with coaxial stacks, bends, non-canonical contacts
              pockets()            where a compound could bind, and how the sites differ between two conformations
              target_sites()       accessible windows for antisense oligonucleotides / probes
              compare()            two ensembles of one RNA (wild type / mutant, two predictors, +/- ligand)

Secondary-structure energies are ViennaRNA nearest-neighbour free energies under the chosen conditions; 3D
quantities come from the coordinates of the conformations. Nothing here uses the model scores.

The pocket geometry itself is `rnafold3d.screening.pockets`, which has no dependency on RCV.
"""
from __future__ import annotations

import math
import re

import numpy as np

from .switch import CANON, compatible, nested_canonical, pairs_of

FRET_R0 = 50.0                  # Å; Förster radius assumed for the label-pair ranking (e.g. Cy3/Cy5: 50-60 Å)
POLAR = {'A': ('N1', 'N3', 'N6', 'N7'), 'G': ('N1', 'N2', 'N3', 'O6', 'N7'), 'C': ('O2', 'N3', 'N4'),
         'U': ('O2', 'N3', 'O4')}
EDGE = {'W': {'A': {'N1', 'N6'}, 'G': {'N1', 'N2', 'O6'}, 'C': {'N3', 'N4', 'O2'}, 'U': {'N3', 'O4', 'O2'}},
        'H': {'A': {'N6', 'N7'}, 'G': {'O6', 'N7'}, 'C': {'N4'}, 'U': {'O4'}},
        'S': {'A': {'N3'}, 'G': {'N2', 'N3'}, 'C': {'O2'}, 'U': {'O2'}}}
EDGE_NAME = {'W': 'Watson-Crick', 'H': 'Hoogsteen', 'S': 'sugar'}


# ------------------------------------------------------------------------------------------- helpers
def _atom_index(doc, name):
    """Atom index of `name` in every residue (-1 if absent)."""
    names, out = doc['atoms']['names'], []
    for r in doc['residues']:
        idx = -1
        for q in range(r['n_atoms']):
            if names[r['first_atom'] + q] == name:
                idx = r['first_atom'] + q
                break
        out.append(idx)
    return np.array(out)


def trace(doc, coords, atom="C4'"):
    """Coordinates of one atom per residue (conformations x residues x 3); residues without it use their first atom."""
    idx = _atom_index(doc, atom)
    first = np.array([r['first_atom'] for r in doc['residues']])
    return np.asarray(coords)[:, np.where(idx >= 0, idx, first), :].astype(float)


def conf_dbn(doc, k):
    """Nested canonical 2D of conformation k (the engine's 2D when present, else the pairs found in 3D)."""
    seq, c = doc['sequence'], doc['conformations'][k]
    dbn = c.get('input_dbn')
    if dbn and len(dbn) == len(seq):
        return nested_canonical(seq, pairs_of(dbn))
    return nested_canonical(seq, [(p[0], p[1]) for p in c.get('pairs', []) if not (len(p) > 3 and p[3])])


def _kabsch(P, Q, w=None):
    """Rotation R and translation t that put P on Q (weights w): P @ R + t."""
    w = np.ones(len(P)) if w is None else w
    p0 = (P * w[:, None]).sum(0) / w.sum()
    q0 = (Q * w[:, None]).sum(0) / w.sum()
    u, _s, vt = np.linalg.svd(((P - p0) * w[:, None]).T @ (Q - q0))
    R = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
    return R, q0 - p0 @ R


def _rmsd(P, Q):
    R, t = _kabsch(P, Q)
    return float(np.sqrt(np.mean(np.sum((P @ R + t - Q) ** 2, axis=1))))


def _runs(flags):
    """[start, end] (inclusive) of the runs of True in a boolean sequence."""
    out, s = [], None
    for i, f in enumerate(list(flags) + [False]):
        if f and s is None:
            s = i
        elif not f and s is not None:
            out.append([s, i - 1])
            s = None
    return out


def _two_state(cond, seq, A, B):
    """E(A) - E(B) for seq; pairs that the sequence cannot form are left out of their state."""
    fc = cond.fold_compound(seq, eval_only=True)
    sites = cond.motif_sites(seq)
    return cond.energy(fc, seq, compatible(seq, A), sites) - cond.energy(fc, seq, compatible(seq, B), sites)


def _unpaired_profile(cond, seq):
    """(probability unpaired per nucleotide, MFE structure, MFE, ensemble free energy)."""
    fc = cond.fold_compound(seq)
    s, e = fc.mfe()
    fc.exp_params_rescale(e)
    _, G = fc.pf()
    n = len(seq)
    P = np.array(fc.bpp())[1:n + 1, 1:n + 1]
    return np.clip(1.0 - P.sum(0) - P.sum(1), 0.0, 1.0), s, e, G


# ------------------------------------------------------------------------------- Switch tab: experiments
def rescue_designs(seq, A, B, cond, top=6):
    """Disrupt-and-rescue double mutants: a first mutation breaks a base pair that only one state has (the
    equilibrium moves to the other state), a second mutation at the partner restores the pair (it moves back).
    The classic proof that a helix exists and matters. Ranked by effect of the first mutation x rescue."""
    kT = cond.kT
    d0 = _two_state(cond, seq, A, B)
    frac = lambda d: 1.0 / (1.0 + math.exp(max(-60.0, min(60.0, d / kT))))
    pa, pb = set(pairs_of(A)), set(pairs_of(B))
    partner = {}
    for name, ps in (('A', pa), ('B', pb)):
        partner[name] = {x: y for i, j in ps for x, y in ((i, j), (j, i))}
    out = []
    for name, own, other in (('A', pa - pb, 'B'), ('B', pb - pa, 'A')):
        sign = 1.0 if name == 'A' else -1.0           # breaking an A pair raises E(A): ddG > 0
        best = {}
        for i, j in sorted(own):
            for p, q in ((i, j), (j, i)):
                o = partner[other].get(p)
                for x in 'ACGU':
                    if x == seq[p] or (x, seq[q]) in CANON:
                        continue
                    if o is not None and (x, seq[o]) not in CANON:
                        continue                       # would also break a pair of the other state
                    m1 = seq[:p] + x + seq[p + 1:]
                    dd1 = _two_state(cond, m1, A, B) - d0
                    if sign * dd1 < 1.0:
                        continue
                    o2 = partner[other].get(q)
                    for y in 'ACGU':
                        if y == seq[q] or (x, y) not in CANON:
                            continue
                        if o2 is not None and (y, m1[o2]) not in CANON:
                            continue
                        m2 = m1[:q] + y + m1[q + 1:]
                        dd2 = _two_state(cond, m2, A, B) - d0
                        rescue = 1.0 - abs(dd2) / abs(dd1)
                        if rescue < 0.5:
                            continue
                        row = {'state': name, 'pair': [min(i, j), max(i, j)],
                               'first': {'pos': p, 'from': seq[p], 'to': x}, 'second': {'pos': q, 'from': seq[q], 'to': y},
                               'ddg_first': round(dd1, 2), 'ddg_double': round(dd2, 2),
                               'frac_a_wt': round(frac(d0), 4), 'frac_a_first': round(frac(d0 + dd1), 4),
                               'frac_a_double': round(frac(d0 + dd2), 4), 'rescue': round(rescue, 3),
                               'sequence_first': m1, 'sequence_double': m2, 'score': abs(dd1) * rescue}
                        key = (min(i, j), max(i, j))
                        if key not in best or row['score'] > best[key]['score']:
                            best[key] = row
        rows = sorted(best.values(), key=lambda r: -r['score'])[:top]
        for r in rows:
            r['score'] = round(r['score'], 3)
        out += rows
    return out


def probing_segments(A, B, top=8):
    """Runs of nucleotides paired in one state and unpaired in the other: where SHAPE/DMS reactivity differs."""
    a = [ch != '.' for ch in A]
    b = [ch != '.' for ch in B]
    segs = []
    for name, flags in (('B', [x and not y for x, y in zip(a, b)]), ('A', [y and not x for x, y in zip(a, b)])):
        segs += [{'start': s, 'end': e, 'length': e - s + 1, 'reactive_in': name} for s, e in _runs(flags)]
    segs.sort(key=lambda s: (-s['length'], s['start']))
    return {'segments': segs[:top], 'informative': sum(x != y for x, y in zip(a, b)), 'length': len(A)}


def label_pairs(doc, coords, confs_a, confs_b, top=5, r0=FRET_R0, min_sep=8):
    """Pairs of nucleotides whose distance differs most between the conformations of two states, scored by the
    difference in FRET efficiency E = 1 / (1 + (d / R0)^6) minus the spread inside each state."""
    if not confs_a or not confs_b:
        return []
    X = trace(doc, coords)
    n = X.shape[1]

    def stats(ks):
        D = np.stack([np.linalg.norm(X[k][:, None] - X[k][None], axis=2) for k in ks[:12]])
        E = 1.0 / (1.0 + (D / r0) ** 6)
        return D.mean(0), E.mean(0), E.std(0)
    da, ea, sa = stats(confs_a)
    db, eb, sb = stats(confs_b)
    score = np.abs(ea - eb) - (sa + sb)
    i, j = np.triu_indices(n, min_sep)
    order = np.argsort(-score[i, j])
    out = []
    for o in order:
        a, b = int(i[o]), int(j[o])
        if score[a, b] < 0.15 or len(out) >= top:
            break
        if any(abs(a - r['i']) < 5 and abs(b - r['j']) < 5 for r in out):
            continue
        out.append({'i': a, 'j': b, 'd_a': round(float(da[a, b]), 1), 'd_b': round(float(db[a, b]), 1),
                    'fret_a': round(float(ea[a, b]), 2), 'fret_b': round(float(eb[a, b]), 2),
                    'score': round(float(score[a, b]), 3)})
    return out


def experiments(doc, states, a, b, cond, coords=None):
    """Experiments that tell state a from state b (indices into `states`, see switch.ensemble_states)."""
    if a is None or b is None or a == b or not (0 <= a < len(states)) or not (0 <= b < len(states)):
        raise ValueError('choose two different states A and B')
    seq, A, B = doc['sequence'], states[a]['dbn'], states[b]['dbn']
    out = {'a': a, 'b': b, 'dbn_a': A, 'dbn_b': B, 'rescue': rescue_designs(seq, A, B, cond),
           'probing': probing_segments(A, B), 'labels': [], 'r0': FRET_R0, 'conditions_text': cond.describe()}
    if coords is not None:
        out['labels'] = label_pairs(doc, coords, states[a]['conformations'], states[b]['conformations'])
    return out


# --------------------------------------------------------------------------------- Switch tab: variants
_VARIANT = re.compile(r'^([ACGUT])?(\d+)([ACGUT])$')


def parse_variants(text, seq, limit=60):
    """'G34U, A12C+G30C' -> [(label, mutant sequence, [positions])]. Substitutions only; T is read as U."""
    out = []
    for token in re.split(r'[\s,;]+', (text or '').strip().upper()):
        if not token:
            continue
        m_seq, pos = list(seq), []
        for part in token.split('+'):
            m = _VARIANT.match(part.replace('>', ''))
            if not m:
                raise ValueError(f'cannot read variant "{part}" (write G34U or 34U; several changes: A12C+G30C)')
            frm, p, to = m.group(1), int(m.group(2)), m.group(3).replace('T', 'U')
            if not 1 <= p <= len(seq):
                raise ValueError(f'variant "{part}": position {p} is outside the sequence (1-{len(seq)})')
            if frm and frm.replace('T', 'U') != seq[p - 1]:
                raise ValueError(f'variant "{part}": position {p} is {seq[p - 1]}, not {frm}')
            if to == seq[p - 1]:
                raise ValueError(f'variant "{part}": position {p} is already {to}')
            m_seq[p - 1] = to
            pos.append(p - 1)
        out.append(('+'.join(f'{seq[p]}{p + 1}{m_seq[p]}' for p in pos), ''.join(m_seq), pos))
        if len(out) > limit:
            raise ValueError(f'at most {limit} variants at a time')
    if not out:
        raise ValueError('enter at least one variant, e.g. G34U')
    return out


def variants(seq, states, cond, text, a=None, b=None):
    """Effect of sequence variants on the ensemble: change of the pairing profile (1 - Pearson correlation of
    the per-nucleotide unpaired probabilities, as in riboSNitch studies), where it changes most, how many base
    pairs of the most stable structure change, and the shift of the A/B equilibrium when two states are chosen."""
    wt_u, wt_s, wt_e, _G = _unpaired_profile(cond, seq)
    two = a is not None and b is not None and a != b and 0 <= a < len(states) and 0 <= b < len(states)
    kT = cond.kT
    d0 = _two_state(cond, seq, states[a]['dbn'], states[b]['dbn']) if two else None
    frac = lambda d: 1.0 / (1.0 + math.exp(max(-60.0, min(60.0, d / kT))))
    wt_pairs = set(pairs_of(wt_s))
    rows = []
    for label, mseq, pos in parse_variants(text, seq):
        u, s, e, _g = _unpaired_profile(cond, mseq)
        if wt_u.std() > 1e-9 and u.std() > 1e-9:
            change = 1.0 - float(np.corrcoef(wt_u, u)[0, 1])
        else:
            change = float(np.abs(wt_u - u).mean())
        diff = np.abs(u - wt_u)
        w = 5
        local = np.convolve(diff, np.ones(w) / w, mode='same')
        row = {'variant': label, 'positions': pos, 'sequence': mseq, 'change': round(max(0.0, change), 4),
               'hotspot': int(np.argmax(local)), 'changed_nt': int((diff > 0.3).sum()),
               'mfe_dbn': s, 'mfe': round(e, 2), 'dmfe': round(e - wt_e, 2),
               'pairs_changed': len(wt_pairs ^ set(pairs_of(s))),
               'profile_diff': [round(float(x), 3) for x in (u - wt_u)]}
        if two:
            d = _two_state(cond, mseq, states[a]['dbn'], states[b]['dbn'])
            row.update({'ddg': round(d - d0, 2), 'frac_a': round(frac(d), 4)})
        row['effect'] = 'large' if row['change'] >= 0.3 or row['pairs_changed'] >= 8 else (
            'moderate' if row['change'] >= 0.1 or row['pairs_changed'] >= 3 else 'small')
        rows.append(row)
    return {'rows': rows, 'wt': {'mfe_dbn': wt_s, 'mfe': round(wt_e, 2), 'frac_a': round(frac(d0), 4) if two else None},
            'two_state': two, 'conditions_text': cond.describe()}


# ------------------------------------------------------------------------ Pathway tab: growing transcript
def cotranscriptional(seq, cond, states=(), a=None, b=None, min_len=8):
    """Most stable structure of every nascent transcript 1..L (5' to 3', as the polymerase makes it).

    For each length: the structure, which base pairs are kept in the full-length structure and which are
    transient, and how much of states A and B has already formed. Rearrangements (several pairs breaking when
    one more nucleotide appears) mark decision points. Each length is at equilibrium: this shows what is
    possible at that length, not the kinetics of getting there."""
    import RNA
    n = len(seq)
    md = cond.model()
    rows, prev = [], set()
    final = None
    structs = {}
    for L in range(min(min_len, n), n + 1):
        s, e = RNA.fold_compound(seq[:L], md).mfe()
        structs[L] = (s, e)
    final = set(pairs_of(structs[n][0]))
    refs = {}
    for name, k in (('a', a), ('b', b)):
        if k is not None and 0 <= k < len(states):
            refs[name] = set(pairs_of(states[k]['dbn']))
    events = []
    for L in sorted(structs):
        s, e = structs[L]
        ps = set(pairs_of(s))
        broken, formed = prev - ps, ps - prev
        code = ['.'] * L
        for i, j in ps:
            code[i] = code[j] = 'f' if (i, j) in final else 't'
        row = {'length': L, 'dbn': s, 'energy': round(e, 2), 'code': ''.join(code),
               'transient': sum((p not in final) for p in ps)}
        for name, ref in refs.items():
            possible = [p for p in ref if p[1] < L]
            row[name] = round(sum(p in ps for p in possible) / len(ref), 3) if ref else 0.0
        rows.append(row)
        if len(broken) >= 3 and rows[:-1]:
            events.append({'length': L, 'broken': len(broken), 'formed': len(formed), 'nucleotide': seq[L - 1]})
        prev = ps
    out = {'rows': rows, 'final': structs[n][0], 'events': events, 'length': n,
           'max_transient': max((r['transient'] for r in rows), default=0),
           'conditions_text': f'{cond.temperature:g} °C, {cond.sodium:g} M Na+'}
    for name in refs:
        first = next((r['length'] for r in rows if r[name] >= 0.9), None)
        peak = max(rows, key=lambda r: r[name])
        out[f'first_{name}'] = first
        out[f'peak_{name}'] = {'length': peak['length'], 'fraction': peak[name]}
    return out


# ------------------------------------------------------------------------------- Details tab: motifs
def _loops(n, pairs):
    """Loops of a nested pair set: [(closing pair or None, [enclosed pairs])] with the exterior loop first."""
    pt = [-1] * n
    for i, j in pairs:
        pt[i], pt[j] = j, i

    def children(lo, hi):
        out, k = [], lo
        while k <= hi:
            if pt[k] > k:
                out.append((k, pt[k]))
                k = pt[k] + 1
            else:
                k += 1
        return out
    loops = [(None, children(0, n - 1))]
    for i, j in sorted(pairs):
        loops.append(((i, j), children(i + 1, j - 1)))
    return loops


def _stems(pairs):
    """Runs of stacked pairs: [[(i, j), (i+1, j-1), ...]]."""
    ps, out = set(pairs), []
    for i, j in sorted(pairs):
        if (i - 1, j + 1) in ps:
            continue
        run = [(i, j)]
        while (run[-1][0] + 1, run[-1][1] - 1) in ps:
            run.append((run[-1][0] + 1, run[-1][1] - 1))
        out.append(run)
    return out


def _axis(X, stem):
    """Unit direction of a helix from its first to its last pair, or None if shorter than 3 pairs.

    The C1'-C1' vectors of successive base pairs rotate about the helix axis, so the axis is the direction
    in which they vary least (a line through the pair midpoints is a poor axis for less than one turn)."""
    if len(stem) < 3:
        return None
    b = np.array([X[j] - X[i] for i, j in stem])
    b /= np.linalg.norm(b, axis=1)[:, None]
    v = np.linalg.svd(b - b.mean(0))[2][2] if len(stem) > 3 else np.cross(b[1] - b[0], b[2] - b[1])
    v = v / (np.linalg.norm(v) or 1.0)
    mid = np.array([(X[i] + X[j]) / 2 for i, j in stem])
    if np.dot(v, mid[-1] - mid[0]) < 0:
        v = -v
    return v


def _angle(u, v):
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))))


def _hairpin_name(loop_seq):
    k = len(loop_seq)
    if k == 4:
        if re.fullmatch('G[ACGU][AG]A', loop_seq):
            return 'GNRA tetraloop'
        if re.fullmatch('U[ACGU]CG', loop_seq):
            return 'UNCG tetraloop'
        if loop_seq == 'CUUG':
            return 'CUUG tetraloop'
        return 'tetraloop'
    return {3: 'triloop', 5: 'pentaloop', 6: 'hexaloop', 7: 'heptaloop'}.get(k, f'{k}-nt hairpin loop')


def _contacts(doc, xyz, seq, canon_partner):
    """Base-base hydrogen-bond contacts outside canonical pairs, A-minor contacts and base triples."""
    from scipy.spatial import cKDTree
    names = doc['atoms']['names']
    pts, owner, aname = [], [], []
    o2, plane = {}, {}
    for r, res in enumerate(doc['residues']):
        base = seq[r]
        ring = [xyz[res['first_atom'] + q] for q in range(res['n_atoms'])
                if names[res['first_atom'] + q] in ('N1', 'C2', 'N3', 'C4', 'C5', 'C6')]
        if len(ring) == 6:
            ring = np.array(ring)
            plane[r] = (ring.mean(0), np.linalg.svd(ring - ring.mean(0))[2][2])
        for q in range(res['n_atoms']):
            nm = names[res['first_atom'] + q]
            if nm in POLAR.get(base, ()):
                pts.append(xyz[res['first_atom'] + q])
                owner.append(r)
                aname.append(nm)
            elif nm == "O2'":
                o2[r] = xyz[res['first_atom'] + q]
    out = []
    if len(pts) < 2:
        return out
    pts = np.array(pts)
    found = {}
    for x, y in cKDTree(pts).query_pairs(3.3):
        i, j = owner[x], owner[y]
        if abs(i - j) < 2 or canon_partner.get(i) == j:
            continue
        key = (min(i, j), max(i, j))
        found.setdefault(key, []).append((aname[x], aname[y]) if i < j else (aname[y], aname[x]))
    for (i, j), bonds in sorted(found.items()):
        if len(bonds) < 2 or i not in plane or j not in plane:
            continue
        # the two bases must lie in one plane (stacked neighbours across a helix also have close N/O atoms)
        (ci, ni), (cj, nj) = plane[i], plane[j]
        if abs(float(np.dot(ni, nj))) < 0.80 or abs(float(np.dot(cj - ci, ni))) > 1.6 or abs(float(np.dot(cj - ci, nj))) > 1.6:
            continue

        def edge(base, atoms):
            score = {e: len(EDGE[e][base] & atoms) for e in 'WHS'}
            return max('WHS', key=lambda e: (score[e], e == 'W'))
        ei, ej = edge(seq[i], {b[0] for b in bonds}), edge(seq[j], {b[1] for b in bonds})
        triple = canon_partner.get(i) is not None or canon_partner.get(j) is not None
        if not triple and (seq[i], seq[j]) in CANON and ei == ej == 'W':
            continue                                   # a Watson-Crick or wobble pair, already part of the 2D
        third = [canon_partner[x] for x in (i, j) if canon_partner.get(x) is not None]
        out.append({'type': 'base triple' if triple else 'non-canonical pair',
                    'label': f'{seq[i]}{i + 1}·{seq[j]}{j + 1}', 'residues': sorted({i, j, *third}),
                    'detail': f'{EDGE_NAME[ei]} / {EDGE_NAME[ej]} edges, {len(bonds)} hydrogen bonds'})
    # A-minor: an unpaired adenine reaching into the minor groove of a helix (N1 or N3 to a paired residue's O2')
    if o2:
        keys = sorted(o2)
        tree = cKDTree(np.array([o2[k] for k in keys]))
        seen = set()
        for x, r in enumerate(owner):
            if seq[r] != 'A' or aname[x] not in ('N1', 'N3') or canon_partner.get(r) is not None:
                continue
            for h in tree.query_ball_point(pts[x], 3.5):
                t = keys[h]
                if abs(t - r) < 3 or canon_partner.get(t) is None or (r, t) in seen:
                    continue
                seen.add((r, t))
                out.append({'type': 'A-minor contact', 'label': f'A{r + 1} → {seq[t]}{t + 1}',
                            'residues': sorted({r, t, canon_partner[t]}),
                            'detail': f'adenine {aname[x]} to the O2′ of a paired nucleotide (minor groove)'})
    return out


def motifs(doc, k, coords):
    """Structural motifs of conformation k: loops by name and size, bends at internal loops and bulges,
    junctions with their coaxial stacks, pseudoknots, non-canonical pairs, base triples and A-minor contacts."""
    seq = doc['sequence']
    n = len(seq)
    c = doc['conformations'][k]
    all_pairs = [(p[0], p[1]) for p in c.get('pairs', [])]
    crossing = [(p[0], p[1]) for p in c.get('pairs', []) if len(p) > 3 and p[3]]
    dbn = conf_dbn(doc, k)
    pairs = pairs_of(dbn)
    xyz = np.asarray(coords)[k].astype(float)
    X = trace(doc, coords, "C1'")[k]
    stems = _stems(pairs)
    stem_of = {p: s for s, st in enumerate(stems) for p in st}
    axes = [_axis(X, st) for st in stems]
    out = []
    for close, kids in _loops(n, pairs):
        if close is None:
            continue
        i, j = close
        if not kids:
            ls = seq[i + 1:j]
            out.append({'type': 'hairpin loop', 'label': _hairpin_name(ls), 'residues': list(range(i, j + 1)),
                        'detail': f"{ls} closed by {seq[i]}{i + 1}-{seq[j]}{j + 1}"})
        elif len(kids) == 1:
            p, q = kids[0]
            l5, l3 = p - i - 1, j - q - 1
            if l5 == 0 and l3 == 0:
                continue
            name = 'bulge' if 0 in (l5, l3) else 'internal loop'
            size = f'{max(l5, l3)} nt' if 0 in (l5, l3) else f'{l5} × {l3} nt'
            u, v = axes[stem_of[(i, j)]], axes[stem_of[(p, q)]]
            bend = _angle(u, v) if u is not None and v is not None else None
            label = f'{name}, {size}'
            if name == 'bulge' and max(l5, l3) == 3 and bend is not None and bend >= 50:
                label = f'kink-turn-like bulge, {size}'
            out.append({'type': name, 'label': label,
                        'residues': list(range(i, p + 1)) + list(range(q, j + 1)),
                        'detail': (f'helix axes bend by {bend:.0f}°' if bend is not None else
                                   'flanking helices too short to measure the bend'),
                        **({'bend': round(bend, 1)} if bend is not None else {})})
        else:
            helices = [(close, True)] + [(kid, False) for kid in kids]
            dirs, ends = [], []
            for (p, q), outer in helices:
                ax = axes[stem_of[(p, q)]]
                dirs.append(None if ax is None else (-ax if outer else ax))    # pointing away from the junction
                ends.append((X[p] + X[q]) / 2)
            stacks = []
            for x in range(len(helices)):
                for y in range(x + 1, len(helices)):
                    if dirs[x] is None or dirs[y] is None:
                        continue
                    if _angle(dirs[x], dirs[y]) >= 135 and np.linalg.norm(ends[x] - ends[y]) <= 10.0:
                        stacks.append((x, y))
            res = sorted({r for (p, q), _o in helices for r in (p, q)} |
                         {r for a, b in zip([i] + [q for _p, q in kids], [p for p, _q in kids] + [j]) for r in range(a + 1, b)})
            name = lambda x: f'{seq[helices[x][0][0]]}{helices[x][0][0] + 1}-{seq[helices[x][0][1]]}{helices[x][0][1] + 1}'
            out.append({'type': 'junction', 'label': f'{len(helices)}-way junction', 'residues': res,
                        'detail': ('coaxial stack: ' + '; '.join(f'helix {name(x)} on helix {name(y)}' for x, y in stacks))
                        if stacks else 'no coaxial stack found', 'coaxial': len(stacks)})
    nested = set(pairs)
    knots = [p for p in crossing if p not in nested] or [p for p in all_pairs if p not in nested and (seq[p[0]], seq[p[1]]) in CANON]
    for st in _stems(sorted(set(knots))):
        if len(st) >= 2:
            out.append({'type': 'pseudoknot', 'label': f'pseudoknot helix, {len(st)} bp',
                        'residues': sorted({r for p in st for r in p}),
                        'detail': f'{seq[st[0][0]]}{st[0][0] + 1}-{seq[st[0][1]]}{st[0][1] + 1} to '
                                  f'{seq[st[-1][0]]}{st[-1][0] + 1}-{seq[st[-1][1]]}{st[-1][1] + 1} crosses the main fold'})
    canon_partner = {}
    for i, j in all_pairs:
        if (seq[i], seq[j]) in CANON:
            canon_partner[i], canon_partner[j] = j, i
    out += _contacts(doc, xyz, seq, canon_partner)
    out.sort(key=lambda m: (m['residues'][0] if m['residues'] else 0, m['type']))
    counts = {}
    for m in out:
        counts[m['type']] = counts.get(m['type'], 0) + 1
    return {'conformation': k, 'dbn': dbn, 'helices': len(stems), 'motifs': out, 'counts': counts}


# --------------------------------------------------------------------------- Details tab: target sites
def target_sites(doc, cond, width=18, top=5):
    """Windows an antisense oligonucleotide or hybridisation probe can reach: the free energy needed to open
    every window of `width` nucleotides (ensemble free energy with the window unpaired minus without that
    condition), and in what fraction of the conformations the window is unpaired."""
    import RNA
    seq = doc['sequence']
    n = len(seq)
    width = int(width)
    if not 8 <= width <= 30:
        raise ValueError('window width must be between 8 and 30 nucleotides')
    if n < width + 2:
        raise ValueError('the RNA is shorter than the window')
    md = cond.model()
    fc = RNA.fold_compound(seq, md)
    _s, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, G0 = fc.pf()
    P = np.array(fc.bpp())[1:n + 1, 1:n + 1]
    unpaired = np.clip(1.0 - P.sum(0) - P.sum(1), 0.0, 1.0)
    conf_unp = np.array([[ch == '.' for ch in conf_dbn(doc, k)] for k in range(len(doc['conformations']))], float)
    step = max(1, n // 400)
    rows = []
    for s in range(0, n - width + 1, step):
        f2 = RNA.fold_compound(seq, md)
        for t in range(s, s + width):
            f2.hc_add_up(t + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        f2.exp_params_rescale(mfe)
        _, G1 = f2.pf()
        win = seq[s:s + width]
        rows.append({'start': s, 'end': s + width - 1, 'dg_open': round(max(0.0, G1 - G0), 2),
                     'gc': round(sum(ch in 'GC' for ch in win) / width, 2),
                     'open_conformations': round(float((conf_unp[:, s:s + width].mean(1) >= 0.8).mean()), 3),
                     'unpaired': round(float(unpaired[s:s + width].mean()), 3), 'sequence': win})
    best = []
    for r in sorted(rows, key=lambda r: (r['dg_open'], -r['open_conformations'])):
        if all(abs(r['start'] - x['start']) >= width // 2 for x in best):
            comp = {'A': 'T', 'C': 'G', 'G': 'C', 'U': 'A'}
            best.append({**r, 'antisense_dna': ''.join(comp[ch] for ch in reversed(r['sequence']))})
        if len(best) >= top:
            break
    return {'width': width, 'profile': [r['dg_open'] for r in rows], 'starts': [r['start'] for r in rows],
            'unpaired': [round(float(x), 3) for x in unpaired], 'best': best,
            'conditions_text': f'{cond.temperature:g} °C, {cond.sodium:g} M Na+'}


# ------------------------------------------------------------------------------ Details tab: compare
def _sequence_map(sa, sb):
    """Matched positions of two sequences of one RNA (identity when the lengths agree)."""
    if len(sa) == len(sb):
        return list(range(len(sa))), list(range(len(sb)))
    from difflib import SequenceMatcher
    ia, ib = [], []
    for blk in SequenceMatcher(None, sa, sb, autojunk=False).get_matching_blocks():
        if blk.size < 4:                              # short chance matches do not make two RNAs the same
            continue
        ia += range(blk.a, blk.a + blk.size)
        ib += range(blk.b, blk.b + blk.size)
    return ia, ib


def compare(doc_a, coords_a, doc_b, coords_b, ka=0, kb=0, max_conf=40):
    """Two ensembles of one RNA, residue by residue: where conformation ka of A and kb of B differ after a
    superposition that follows the part that stays the same, how far the local environment agrees (a local
    distance difference score, insensitive to domain movements), where the base pairing differs across the
    ensembles, and which conformations of one ensemble have a counterpart in the other."""
    sa, sb = doc_a['sequence'], doc_b['sequence']
    ia, ib = _sequence_map(sa, sb)
    if len(ia) < max(10, 0.7 * min(len(sa), len(sb))):
        raise ValueError('the two ensembles do not look like the same RNA (their sequences do not match)')
    XA, XB = trace(doc_a, coords_a), trace(doc_b, coords_b)
    if not (0 <= ka < len(XA) and 0 <= kb < len(XB)):
        raise ValueError('conformation index out of range')
    A, B = XA[ka][ia], XB[kb][ib]
    L = len(ia)
    d0 = max(0.5, 0.6 * (max(L, 13) - 0.5) ** 0.5 - 2.5)
    w = np.ones(L)
    for _ in range(15):                              # superposition weighted towards the common core
        R, t = _kabsch(A, B, w)
        d = np.linalg.norm(A @ R + t - B, axis=1)
        w = 1.0 / (1.0 + (d / d0) ** 2)
    DA = np.linalg.norm(A[:, None] - A[None], axis=2)
    DB = np.linalg.norm(B[:, None] - B[None], axis=2)
    near = (DA < 15.0) & ~np.eye(L, dtype=bool)
    diff = np.abs(DA - DB)
    local = np.zeros(L)
    for thr in (0.5, 1.0, 2.0, 4.0):
        local += ((diff < thr) & near).sum(1) / np.maximum(1, near.sum(1)) / 4.0

    def paired(doc, idx):
        f = np.zeros(len(doc['sequence']))
        for k in range(len(doc['conformations'])):
            f += np.array([ch != '.' for ch in conf_dbn(doc, k)], float)
        return (f / max(1, len(doc['conformations'])))[idx]
    pair_diff = paired(doc_b, ib) - paired(doc_a, ia)
    na, nb = min(len(XA), max_conf), min(len(XB), max_conf)
    M = np.array([[_rmsd(XA[x][ia], XB[y][ib]) for y in range(nb)] for x in range(na)])
    moved = d > max(3.0, 2.0 * float(np.median(d)))
    seg = [{'start': ia[s], 'end': ia[e], 'max': round(float(d[s:e + 1].max()), 1)} for s, e in _runs(moved) if e - s >= 1]
    per_a = lambda v, nd: [None if x is None else round(float(x), nd) for x in
                           (dict(zip(ia, v)).get(r) for r in range(len(sa)))]
    return {'matched': L, 'length_a': len(sa), 'length_b': len(sb), 'ka': ka, 'kb': kb,
            'differences': [{'pos_a': x, 'a': sa[x], 'pos_b': y, 'b': sb[y]} for x, y in zip(ia, ib) if sa[x] != sb[y]],
            'rmsd': round(_rmsd(A, B), 2), 'tm': round(float(w.sum() / L), 3), 'core': int((d < 2.0 * d0).sum()),
            'local_agreement_mean': round(float(local.mean()), 3),
            'displacement': per_a(d, 2), 'local_agreement': per_a(local, 3), 'pair_difference': per_a(pair_diff, 3),
            'moved': sorted(seg, key=lambda s: -s['max'])[:8],
            'pairing_changed': [{'start': ia[s], 'end': ia[e]} for s, e in _runs(np.abs(pair_diff) >= 0.5)],
            'nearest_a': [{'conformation': x, 'other': int(M[x].argmin()), 'rmsd': round(float(M[x].min()), 2)} for x in range(na)],
            'nearest_b': [{'conformation': y, 'other': int(M[:, y].argmin()), 'rmsd': round(float(M[:, y].min()), 2)} for y in range(nb)],
            'cross_rmsd_min': round(float(M.min()), 2), 'cross_rmsd_mean': round(float(M.mean()), 2)}


# ------------------------------------------------------------------------------- Details tab: pockets
def _flexible(doc, k):
    """Per nucleotide, 1 when it is not in a Watson-Crick helix of conformation k, 0 when it is.

    Small molecules bind loops, bulges and junctions; the grooves of an A-form helix are where docking
    puts its false positives, so the pocket score needs to know the difference."""
    paired = np.array([ch != '.' for ch in conf_dbn(doc, k)], float)
    return 1.0 - paired


def pockets(doc, k, coords, other=None, top=10, spacing=None, burial_min=None):
    """Pockets of conformation k — where a small molecule could bind — and how they change in conformation
    `other`.

    The geometry is `rnafold3d.screening.pockets`: a buriedness scan of the free space around the atoms,
    clustered into pockets, each described by its volume, how enclosed it is, the nucleotides lining it and
    how much of that lining lies outside a helix (small molecules bind loops, bulges and junctions). With a
    second conformation, every site — the nucleotides lining a pocket of either conformation — is measured
    in both, and the difference in open volume says which site opens in which conformation: the sites a
    conformation-selective compound would use.

    `burial_min` is how enclosed a point must be to count: 5 of 7 directions by default, which is what
    finds the sites of known RNA-ligand complexes (`benchmark/screening/pocket_sites.py`); 6 reports only
    the deep pockets, which are tidier and miss the shallow sites most RNA binders use.

    Volumes are grid volumes (±10% at the default 0.8 Å spacing); nothing here is an affinity, and a pocket
    in a predicted model is only as reliable as the model."""
    from rnafold3d.screening import pockets as P
    confs = doc['conformations']
    if not 0 <= k < len(confs):
        raise ValueError('conformation index out of range')
    if other is not None and not 0 <= other < len(confs):
        raise ValueError('conformation index out of range')
    if other == k:
        other = None
    xyz = np.asarray(coords).astype(float)
    n = len(doc['sequence'])
    elems, res_of = doc['atoms']['elements'], doc['atoms']['residue']
    kw = dict(elements=elems, residue_of=res_of, **({} if spacing is None else {'spacing': spacing}),
              **({} if burial_min is None else {'burial_min': int(burial_min)}))
    found = P.find(xyz[k], flexible=_flexible(doc, k), **kw)
    segs = lambda res: [{'start': s, 'end': e} for s, e in _runs(np.isin(np.arange(n), res))]
    rows = [{'pocket': p['id'], 'volume': p['volume'], 'burial': p['burial'], 'enclosure': p['enclosure'],
             'extent': p['extent'], 'loop_fraction': p['loop_fraction'], 'score': p['score'],
             'wide': p['wide'], 'center': p['center'], 'residues': p['residues'], 'segments': segs(p['residues'])}
            for p in found[:top]]
    out = {'conformation': k, 'other': other, 'spacing': found[0]['spacing'] if found else P.SPACING,
           'n_pockets': len(found), 'pockets': rows,
           'params': {'probe': P.PROBE, 'burial_min': int(burial_min or P.BURIAL_MIN), 'min_volume': P.MIN_VOLUME,
                      'contact': P.CONTACT, 'wide_volume': P.WIDE_VOLUME}}
    if other is None:
        return out
    found_b = P.find(xyz[other], flexible=_flexible(doc, other), **kw)
    changes, by_a = [], {}
    for r in P.compare(found, found_b, xyz[k], xyz[other], res_of):
        a, b = r['a'], r['b']
        entry = {'side': r['side'], 'pocket': None if a is None else found[a]['id'],
                 'other_pocket': None if b is None else found_b[b]['id'], 'overlap': r['overlap'],
                 'volume': r['volume_a'], 'volume_other': r['volume_b'], 'difference': r['difference'],
                 'open_in': r['open_in'], 'residues': r['residues'], 'segments': segs(r['residues']),
                 'score': found[a]['score'] if a is not None else None,
                 'score_other': found_b[b]['score'] if b is not None else None}
        changes.append(entry)
        if r['side'] == 'A':
            by_a[entry['pocket']] = entry
    for row in rows:                                  # the table of this conformation's pockets carries the comparison
        c = by_a.get(row['pocket'])
        # volume_site / volume_other are the open volume at this pocket's nucleotides in the two conformations
        # (comparable between them, see screening.pockets.volume_near); 'volume' stays the pocket's own volume
        row.update({'volume_site': c['volume'] if c else row['volume'], 'volume_other': c['volume_other'] if c else 0.0,
                    'difference': c['difference'] if c else 1.0, 'open_in': c['open_in'] if c else 'A',
                    'other_pocket': c['other_pocket'] if c else None})
    out['n_pockets_other'] = len(found_b)
    out['changes'] = changes[:2 * top]
    return out
