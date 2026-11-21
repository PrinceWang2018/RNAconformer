"""Z-RNA: identification of left-handed double-helical regions (0.7.6).

Z-RNA (Hall et al. 1984; Popenda et al. 2004; Placido et al. 2007, ADAR1 Zalpha with Z-RNA) is a left-handed double
helix of alternating purine-pyrimidine repeats (CG best, also CA/UG). Its signature in 3D:

* left-handed twist: about -60 degrees per dinucleotide repeat (two base-pair steps of unequal twist; A-form RNA is
  about +31 degrees per step);
* the purines syn, the pyrimidines anti (each pair has its syn purine);
* consecutive Watson-Crick pairs (G-C, A-U, G-U) stacked into a helix.

A region is reported only with high confidence: at least MIN_PAIRS consecutive stacked pairs, every two-step
window left-handed (sum of twists <= WINDOW_MAX), and the syn purine / anti pyrimidine pattern in nearly every pair.
An isolated distorted step, a single syn purine or a kink never qualifies. The same rules are implemented in
rcv/static/js/zrna.js for frames that move (morph, pull); tests keep the two identical.

Twist of a step between stacked pairs (i, j) and (i+1, j-1): the angle from the C1'(i)->C1'(j) vector of the first
pair to that of the second, about the local helix axis (mean base-pair normal, oriented 5'->3' along the strand of
i); positive = right-handed.
"""
from __future__ import annotations

import math

import numpy as np

MIN_PAIRS = 3
WINDOW_MAX = -30.0         # degrees, sum of two consecutive step twists
STEP_MAX = 25.0            # no single step may be clearly right-handed
SYN_FRACTION = 0.75        # pairs whose purine is syn
ANTI_FRACTION = 0.9        # pyrimidines anti
SYN_RANGE = (-20.0, 120.0)  # chi O4'-C1'-N9-C4 / O4'-C1'-N1-C2
RING = ('C2', 'C4', 'C6')
C1_RANGE = (8.5, 12.5)     # C1'-C1' of a pair in the analysed frame


def chi_of(atoms, base):
    g, ref = ('N9', 'C4') if base in 'AG' else ('N1', 'C2')
    if not all(k in atoms for k in ("O4'", "C1'", g, ref)):
        return None
    from .geomutil import dihedral
    return float(dihedral(*[np.asarray(atoms[k], float) for k in ("O4'", "C1'", g, ref)]))


def is_syn(chi):
    return chi is not None and SYN_RANGE[0] <= chi <= SYN_RANGE[1]


def _normal(atoms):
    if not all(k in atoms for k in RING):
        return None
    a, b, c = (np.asarray(atoms[k], float) for k in RING)
    n = np.cross(b - a, c - a)
    ln = np.linalg.norm(n)
    return n / ln if ln > 1e-6 else None


def step_twist(ai, aj, ak, al):
    """Twist (degrees) of the step from pair (i, j) to pair (k, l) = (i+1, j-1); atoms dicts; None if incomplete."""
    need = ("C1'",)
    if not all(all(n in x for n in need) for x in (ai, aj, ak, al)):
        return None
    ns = [_normal(x) for x in (ai, aj, ak, al)]
    if any(n is None for n in ns):
        return None
    c = {k: np.asarray(x["C1'"], float) for k, x in zip('ijkl', (ai, aj, ak, al))}
    down = (c['k'] - c['i']) + (c['l'] - c['j'])         # 5'->3' along the strand of i (strand 2 runs the other way)
    h = np.zeros(3)
    for n in ns:
        h += n if np.dot(n, down) >= 0 else -n
    ln = np.linalg.norm(h)
    if ln < 1e-6:
        return None
    h /= ln
    a1, a2 = c['j'] - c['i'], c['l'] - c['k']
    a1 = a1 - np.dot(a1, h) * h
    a2 = a2 - np.dot(a2, h) * h
    if np.linalg.norm(a1) < 1e-6 or np.linalg.norm(a2) < 1e-6:
        return None
    return math.degrees(math.atan2(float(np.dot(np.cross(a1, a2), h)), float(np.dot(a1, a2))))


def helix_runs(pairs):
    """Runs of stacked pairs [(i, j), (i+1, j-1), ...] (0-based, i < j)."""
    ps = sorted((min(i, j), max(i, j)) for i, j in pairs)
    have = set(ps)
    runs, used = [], set()
    for p in ps:
        if p in used or (p[0] - 1, p[1] + 1) in have:
            continue
        run = [p]
        used.add(p)
        while (run[-1][0] + 1, run[-1][1] - 1) in have and run[-1][1] - run[-1][0] > 2:
            run.append((run[-1][0] + 1, run[-1][1] - 1))
            used.add(run[-1])
        runs.append(run)
    return runs


def detect(seq, residues, pairs, chi=None, min_pairs=MIN_PAIRS):
    """High-confidence Z-RNA regions.

    seq: parent sequence; residues: list of {atom: xyz} (needs C1', C2, C4, C6 and, without `chi`, O4' and the
    glycosidic atoms); pairs: (i, j) 0-based base pairs; chi: optional list of chi angles (degrees) per residue.
    Returns [{'pairs': [(i, j) ...] 1-based, 'residues': [...] 1-based, 'mean_twist', 'twists', 'syn_purines',
    'anti_pyrimidines', 'confidence': 'high'}]."""
    n = len(seq)
    chis = list(chi) if chi is not None else [chi_of(residues[r], seq[r]) if r < len(residues) else None
                                               for r in range(n)]
    def paired(i, j):
        a, b = residues[i].get("C1'"), residues[j].get("C1'")
        return a is not None and b is not None and C1_RANGE[0] <= float(np.linalg.norm(np.asarray(a) - np.asarray(b))) <= C1_RANGE[1]
    out = []
    for run in helix_runs([(i, j) for i, j in pairs if paired(min(i, j), max(i, j))]):
        if len(run) < min_pairs:
            continue
        tw = [step_twist(residues[i], residues[j], residues[k], residues[l])
              for (i, j), (k, l) in zip(run, run[1:])]
        z_pair = []
        for i, j in run:
            pur = [r for r in (i, j) if seq[r] in 'AG']
            pyr = [r for r in (i, j) if seq[r] in 'CU']
            syn_ok = len(pur) == 1 and is_syn(chis[pur[0]])
            anti_ok = all(chis[r] is not None and not is_syn(chis[r]) for r in pyr)
            z_pair.append((syn_ok, anti_ok, len(pyr)))
        # maximal windows of consecutive pairs meeting the rules
        a = 0
        while a < len(run):
            b = a
            while b + 1 < len(run):
                t = tw[b]
                if t is None or t > STEP_MAX:
                    break
                if b > a and (tw[b - 1] is None or tw[b - 1] + t > WINDOW_MAX):
                    break
                b += 1
            if b - a + 1 >= min_pairs:
                seg = list(range(a, b + 1))
                # trim pairs at either end whose purine is not syn (a helix may run on in A-form)
                while seg and not z_pair[seg[0]][0]:
                    seg = seg[1:]
                while seg and not z_pair[seg[-1]][0]:
                    seg = seg[:-1]
                if len(seg) >= min_pairs:
                    steps = [tw[k] for k in seg[:-1]]
                    syn_f = sum(z_pair[k][0] for k in seg) / len(seg)
                    n_pyr = sum(z_pair[k][2] for k in seg)
                    anti_f = (sum(z_pair[k][1] * z_pair[k][2] for k in seg) / n_pyr) if n_pyr else 0.0
                    windows_ok = all(s1 + s2 <= WINDOW_MAX for s1, s2 in zip(steps, steps[1:])) and len(steps) >= 2
                    if syn_f >= SYN_FRACTION and anti_f >= ANTI_FRACTION and windows_ok:
                        out.append(_region(run, seg, steps, syn_f, anti_f))
            a = b + 1
    return out


def _region(run, seg, steps, syn_f, anti_f):
    prs = [run[k] for k in seg]
    res = sorted({r for p in prs for r in p})
    return {'pairs': [(i + 1, j + 1) for i, j in prs], 'residues': [r + 1 for r in res],
            'mean_twist': round(float(np.mean(steps)), 1), 'twists': [round(float(t), 1) for t in steps],
            'syn_purines': round(syn_f, 2), 'anti_pyrimidines': round(anti_f, 2), 'confidence': 'high'}


# ---------------------------------------------------------------------------------------------- prediction (engine)
# Free energy of the A -> Z transition of a helical segment (kcal/mol, estimates). Z-RNA needs alternating
# purine-pyrimidine pairs; CG steps convert most easily, CA/UG steps less, AU steps hardly (as for Z-DNA, Ho et al.
# 1986, but RNA needs more salt than DNA). Scaled so that r(CG)n duplexes turn left-handed near 3-4 M Na+ (Hall et
# al. 1984: 6 M NaClO4; Popenda et al. 2004) and the r(CGC[m8G]CG) duplex near physiological salt (Balasubramaniyam
# et al. 2018). Use as a hypothesis generator, not a quantitative prediction.
Z_STEP = {'GC': 1.0, 'MIX': 1.7, 'AU': 2.8}      # per base-pair step at 0.15 M Na+
Z_SALT_SLOPE = 0.75                               # kcal/mol per step per decade of Na+
Z_JUNCTION = 3.5                                  # each A-Z junction inside a helix (one pair is broken, Ha 2005)
Z_END = 0.5                                       # Z segment ending at a loop or a helix end
Z_SYN_MOD = {'m8G': -3.0, '8BrG': -3.0, 'm8A': -2.0, '8oxoG': -1.0, '8oxoA': -0.8}
Z_BINDER = -1.5                                   # per step, with a Z-binding protein (ADAR1 Zalpha, ZBP1, PKZ; ~2 domains per 6 bp, Kd ~0.1 uM)
Z_MIN_PAIRS = 4
Z_PROPOSE_MAX = 3.0                               # propose a Z hypothesis up to this dG (about 1 % at 37 C)
RT37 = 0.6163


def _pair_ok(a, b):
    return a + b in ('GC', 'CG', 'AU', 'UA')


def _step_kind(p, q, seq):
    gc = sum(seq[i] in 'GC' for i in (p[0], q[0]))
    return 'GC' if gc == 2 else ('AU' if gc == 0 else 'MIX')


def segment_energy(seq, run, a, b, salt=0.15, mods=(), binder=False):
    """dG(A -> Z) of pairs run[a..b] of a helix `run` (list of stacked pairs), with its terms."""
    seg = run[a:b + 1]
    steps = [_step_kind(p, q, seq) for p, q in zip(seg, seg[1:])]
    salt_corr = -Z_SALT_SLOPE * math.log10(max(salt, 1e-3) / 0.15)
    e_steps = sum(Z_STEP[k] + salt_corr for k in steps)
    n_junction = int(a > 0) + int(b < len(run) - 1)
    e_ends = Z_JUNCTION * n_junction + Z_END * (2 - n_junction)
    res = {r for p in seg for r in p}
    by_pos = {int(m['pos']) - 1: m['code'] for m in mods or ()}
    e_mod = sum(Z_SYN_MOD.get(by_pos[r], 0.0) for r in res if r in by_pos and seq[r] in 'AG')
    e_bind = Z_BINDER * len(steps) if binder else 0.0
    dg = e_steps + e_ends + e_mod + e_bind
    return {'dG_kcal_mol': round(dg, 2), 'steps': round(e_steps, 2), 'junctions': round(e_ends, 2),
            'modifications': round(e_mod, 2), 'binder': round(e_bind, 2), 'n_junctions': n_junction,
            'z_fraction_estimate': round(1.0 / (1.0 + math.exp(min(50.0, dg / RT37))), 3)}


def candidates(seq, pairs, salt=0.15, mods=(), binder=False, min_pairs=Z_MIN_PAIRS):
    """Segments of helices that could turn left-handed, best first.

    A segment is a run of >= min_pairs stacked Watson-Crick pairs (G-C, A-U) whose strand alternates purine /
    pyrimidine; blocked or non-canonical pairs end it. Every window is scored with segment_energy, the best window
    of each helix is kept."""
    out = []
    for run in helix_runs(pairs):
        ok = [_pair_ok(seq[i], seq[j]) for i, j in run]
        best = None
        for a in range(len(run)):
            for b in range(a + min_pairs - 1, len(run)):
                if not all(ok[a:b + 1]):
                    break
                if any((seq[run[k][0]] in 'AG') == (seq[run[k + 1][0]] in 'AG') for k in range(a, b)):
                    break
                e = segment_energy(seq, run, a, b, salt, mods, binder)
                if best is None or e['dG_kcal_mol'] < best[2]['dG_kcal_mol']:
                    best = (a, b, e)
        if best:
            a, b, e = best
            out.append({'run': run, 'a': a, 'b': b, 'pairs': run[a:b + 1], **e})
    return sorted(out, key=lambda c: c['dG_kcal_mol'])


def hypothesis(seq, hyp, cand):
    """A 2D hypothesis with one helix segment in Z form: the pair next to each A-Z junction is opened (the junction
    extrudes a base pair) and flanking A-form fragments shorter than two pairs are dropped."""
    run, a, b = cand['run'], cand['a'], cand['b']
    drop = set()
    if a > 0:
        drop.add(run[a - 1])
        if a - 1 < 2:
            drop.update(run[:a - 1])
    if b < len(run) - 1:
        drop.add(run[b + 1])
        if len(run) - b - 2 < 2:
            drop.update(run[b + 2:])
    pairs = [p for p in hyp['pairs'] if tuple(p) not in drop]
    n = len(seq)
    db = ['.'] * n
    for i, j in pairs:
        db[i], db[j] = '(', ')'
    h = {k: v for k, v in hyp.items() if k not in ('torsions', 'template_atoms', 'distance_restraints', 'restraints')}
    seg = [(i + 1, j + 1) for i, j in cand['pairs']]
    e0 = hyp.get('energy_kcal_mol')
    h.update({'pairs': pairs, 'dbn': ''.join(db), 'z_pairs': [tuple(p) for p in cand['pairs']],
              'source': f"Z-RNA form of {hyp['id']} ({hyp.get('source', '')})", 'parent': hyp['id'],
              'zrna': {'pairs': seg, 'residues': sorted({r for p in seg for r in p}),
                       'opened_pairs': [(i + 1, j + 1) for i, j in sorted(drop)],
                       **{k: cand[k] for k in ('dG_kcal_mol', 'steps', 'junctions', 'modifications', 'binder',
                                               'n_junctions', 'z_fraction_estimate')}}})
    if e0 is not None:
        h['energy_kcal_mol'] = round(e0 + cand['dG_kcal_mol'], 2)
    re = hyp.get('rank_energy_kcal_mol')
    if re is not None:
        h['rank_energy_kcal_mol'] = round(re + cand['dG_kcal_mol'], 2)
    h.pop('template', None)
    return h
