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


def analyse(seq, states, temperature=37.0, a=None, b=None, top=8, cond=None):
    """states: list of dicts with 'dbn'. Returns energies, populations and (for a != b) the mutation scan.

    cond (conditions.Conditions) sets temperature, Na+, ligand and probing data; without it only the temperature."""
    import RNA
    from .conditions import Conditions
    cond = cond or Conditions(temperature=float(temperature))
    temperature = cond.temperature
    kT = cond.kT
    use_probe = cond.probing is not None
    fc = cond.fold_compound(seq, probing=use_probe)
    sites = cond.motif_sites(seq)
    mfe_dbn, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, G = fc.pf()
    out_states = []
    for st in states:
        e = cond.energy(fc, seq, st['dbn'], sites)
        out_states.append({**st, 'energy': round(e, 2), 'p_ensemble': math.exp(-(e - G) / kT)})
    emin = min(s['energy'] for s in out_states)
    w = [math.exp(-(s['energy'] - emin) / kT) for s in out_states]
    for s, x in zip(out_states, w):
        s['relative'] = x / sum(w)
    result = {'temperature': temperature, 'kT': kT, 'mfe': {'dbn': mfe_dbn, 'energy': round(mfe, 2)},
              'ensemble_energy': round(G, 2), 'mfe_probability': math.exp(-(mfe - G) / kT), 'states': out_states,
              'vienna': getattr(RNA, '__version__', ''), 'conditions': cond.to_json(), 'conditions_text': cond.describe(),
              'probing_used': use_probe}
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
            mfc = cond.fold_compound(mseq, eval_only=True)
            msites = cond.motif_sites(mseq)
            sa, sb = compatible(mseq, A), compatible(mseq, B)
            ea, eb = cond.energy(mfc, mseq, sa, msites), cond.energy(mfc, mseq, sb, msites)
            ddg = (ea - eb) - d0
            row[base] = {'ddg': round(ddg, 2), 'frac_a': round(frac(ea - eb), 4), 'breaks_a': sa != A, 'breaks_b': sb != B}
        scan.append(row)
    flat = [(i, base, v['ddg']) for i, row in enumerate(scan) for base, v in row.items()]
    to_a = sorted((x for x in flat if x[2] < -0.05), key=lambda x: x[2])[:top]
    to_b = sorted((x for x in flat if x[2] > 0.05), key=lambda x: -x[2])[:top]

    def detail(i, base):
        mseq = seq[:i] + base + seq[i + 1:]
        mfc = cond.fold_compound(mseq)
        msites = cond.motif_sites(mseq)
        ms, me = mfc.mfe()
        mfc.exp_params_rescale(me)
        _, mg = mfc.pf()
        pa = math.exp(-(cond.energy(mfc, mseq, compatible(mseq, A), msites) - mg) / kT)
        pb = math.exp(-(cond.energy(mfc, mseq, compatible(mseq, B), msites) - mg) / kT)
        return {'pos': i, 'from': seq[i], 'to': base, 'ddg': scan[i][base]['ddg'], 'frac_a': scan[i][base]['frac_a'],
                'p_a': pa, 'p_b': pb, 'mfe_dbn': ms, 'mfe': round(me, 2), 'sequence': mseq}
    result['scan'] = {'a': a, 'b': b, 'dbn_a': A, 'dbn_b': B, 'ddg_wt': round(d0, 2), 'frac_a_wt': round(frac(d0), 4),
                      'p_a_wt': out_states[a]['p_ensemble'], 'p_b_wt': out_states[b]['p_ensemble'], 'matrix': scan,
                      'towards_a': [detail(i, base) for i, base, _ in to_a],
                      'towards_b': [detail(i, base) for i, base, _ in to_b]}
    return result


def dose_response(seq, states, cond, concentrations):
    """Populations of the states (share among them and in the full ensemble) over ligand concentrations (µM)."""
    from .conditions import Conditions
    if not cond.ligand:
        raise ValueError('choose a ligand first')
    rows = []
    for c in concentrations:
        cc = Conditions(temperature=cond.temperature, sodium=cond.sodium, ligand={**cond.ligand, 'conc_uM': float(c)})
        fc = cc.fold_compound(seq)
        sites = cc.motif_sites(seq)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        _, G = fc.pf()
        e = [cc.energy(fc, seq, st['dbn'], sites) for st in states]
        emin = min(e)
        w = [math.exp(-(x - emin) / cc.kT) for x in e]
        rows.append({'conc_uM': float(c), 'binding_energy': round(cc.binding_energy(), 3),
                     'relative': [x / sum(w) for x in w], 'p_ensemble': [math.exp(-(x - G) / cc.kT) for x in e],
                     'bound': [cc.motif_bonus(st['dbn'], sites) < 0 for st in states]})
    return {'rows': rows, 'ligand': cond.ligand, 'sites': len(cond.motif_sites(seq))}
