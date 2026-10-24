"""Switch analysis at the secondary-structure level (ViennaRNA): populations of the ensemble's 2D states and a
single-mutation scan that shows which mutations shift the equilibrium between two states.

The 2D states are the canonical, nested base-pair sets of the ensemble's conformations (the engine's 2D hypothesis
when present, otherwise the pairs annotated in 3D). Energies are Turner nearest-neighbour free energies from
ViennaRNA; the 3D model scores are not used. For a mutation, a pair that is no longer canonical is removed from the
state it belonged to, so the state is evaluated without it.
"""
from __future__ import annotations

import math

CANON = {('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G'), ('G', 'U'), ('U', 'G')}
R_KCAL = 0.0019872036


def pairs_of(dbn):
    st, out = {'(': [], '[': [], '{': [], '<': []}, []
    close = {')': '(', ']': '[', '}': '{', '>': '<'}
    for i, ch in enumerate(dbn):
        if ch in st:
            st[ch].append(i)
        elif ch in close and st[close[ch]]:
            out.append((st[close[ch]].pop(), i))
    return sorted(out)


def nested_canonical(seq, pairs):
    """Dot-bracket of the canonical pairs that form a nested (pseudoknot-free) set, greedily by helix order."""
    s = ['.'] * len(seq)
    taken = []
    for i, j in sorted(pairs, key=lambda p: (p[1] - p[0], p[0])):
        if (seq[i], seq[j]) not in CANON or j - i < 4 or s[i] != '.' or s[j] != '.':
            continue
        if any(a < i < b < j or i < a < j < b for a, b in taken):
            continue
        taken.append((i, j))
        s[i], s[j] = '(', ')'
    return ''.join(s)


def compatible(seq, dbn):
    """Remove pairs that are not canonical in seq."""
    s = list(dbn)
    for i, j in pairs_of(dbn):
        if (seq[i], seq[j]) not in CANON:
            s[i] = s[j] = '.'
    return ''.join(s)


def ensemble_states(doc):
    """Distinct 2D states of a document: [{'dbn', 'conformations': [indices], 'source'}], most frequent first."""
    seq = doc['sequence']
    states = {}
    for k, c in enumerate(doc['conformations']):
        src = 'engine 2D hypothesis'
        dbn = c.get('input_dbn')
        if dbn and len(dbn) == len(seq):
            dbn = nested_canonical(seq, pairs_of(dbn))
        else:
            src = 'pairs annotated in 3D'
            dbn = nested_canonical(seq, [(p[0], p[1]) for p in c.get('pairs', []) if not (len(p) > 3 and p[3])])
        st = states.setdefault(dbn, {'dbn': dbn, 'conformations': [], 'source': src})
        st['conformations'].append(k)
    return sorted(states.values(), key=lambda s: (-len(s['conformations']), s['conformations'][0]))


def analyse(seq, states, temperature=37.0, a=None, b=None, top=8):
    """states: list of dicts with 'dbn'. Returns energies, populations and (for a != b) the mutation scan."""
    import RNA
    md = RNA.md()
    md.temperature = float(temperature)
    kT = R_KCAL * (float(temperature) + 273.15)
    fc = RNA.fold_compound(seq, md)
    mfe_dbn, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, G = fc.pf()
    out_states = []
    for st in states:
        e = fc.eval_structure(st['dbn'])
        out_states.append({**st, 'energy': round(e, 2), 'p_ensemble': math.exp(-(e - G) / kT)})
    emin = min(s['energy'] for s in out_states)
    w = [math.exp(-(s['energy'] - emin) / kT) for s in out_states]
    for s, x in zip(out_states, w):
        s['relative'] = x / sum(w)
    result = {'temperature': temperature, 'kT': kT, 'mfe': {'dbn': mfe_dbn, 'energy': round(mfe, 2)},
              'ensemble_energy': round(G, 2), 'mfe_probability': math.exp(-(mfe - G) / kT), 'states': out_states,
              'vienna': getattr(RNA, '__version__', '')}
    if a is None or b is None or a == b or not (0 <= a < len(states)) or not (0 <= b < len(states)):
        return result
    A, B = states[a]['dbn'], states[b]['dbn']
    eA, eB = out_states[a]['energy'], out_states[b]['energy']
    d0 = eA - eB
    frac = lambda d: 1.0 / (1.0 + math.exp(d / kT))            # fraction of A in the two-state A/B equilibrium
    scan = []
    for i, wt in enumerate(seq):
        row = {}
        for base in 'ACGU':
            if base == wt:
                continue
            mseq = seq[:i] + base + seq[i + 1:]
            mfc = RNA.fold_compound(mseq, md, RNA.OPTION_EVAL_ONLY)
            sa, sb = compatible(mseq, A), compatible(mseq, B)
            ea, eb = mfc.eval_structure(sa), mfc.eval_structure(sb)
            ddg = (ea - eb) - d0
            row[base] = {'ddg': round(ddg, 2), 'frac_a': round(frac(ea - eb), 4), 'breaks_a': sa != A, 'breaks_b': sb != B}
        scan.append(row)
    flat = [(i, base, v['ddg']) for i, row in enumerate(scan) for base, v in row.items()]
    to_a = sorted((x for x in flat if x[2] < -0.05), key=lambda x: x[2])[:top]
    to_b = sorted((x for x in flat if x[2] > 0.05), key=lambda x: -x[2])[:top]

    def detail(i, base):
        mseq = seq[:i] + base + seq[i + 1:]
        mfc = RNA.fold_compound(mseq, md)
        ms, me = mfc.mfe()
        mfc.exp_params_rescale(me)
        _, mg = mfc.pf()
        pa = math.exp(-(mfc.eval_structure(compatible(mseq, A)) - mg) / kT)
        pb = math.exp(-(mfc.eval_structure(compatible(mseq, B)) - mg) / kT)
        return {'pos': i, 'from': seq[i], 'to': base, 'ddg': scan[i][base]['ddg'], 'frac_a': scan[i][base]['frac_a'],
                'p_a': pa, 'p_b': pb, 'mfe_dbn': ms, 'mfe': round(me, 2), 'sequence': mseq}
    result['scan'] = {'a': a, 'b': b, 'dbn_a': A, 'dbn_b': B, 'ddg_wt': round(d0, 2), 'frac_a_wt': round(frac(d0), 4),
                      'p_a_wt': out_states[a]['p_ensemble'], 'p_b_wt': out_states[b]['p_ensemble'], 'matrix': scan,
                      'towards_a': [detail(i, base) for i, base, _ in to_a],
                      'towards_b': [detail(i, base) for i, base, _ in to_b]}
    return result
