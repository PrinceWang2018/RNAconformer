"""Long RNAs: fold independently folding domains separately, then assemble the whole chain.

The fold-tree search scales badly with length: more loops to close and more helices to arrange per trajectory,
so above a few hundred nucleotides most trajectories fail or take very long. Large RNAs are built from domains
that the secondary structure makes obvious: the helices that leave the exterior loop, and below a long-range
helix, the branches of the multiloop it closes. This module

1. finds the domains (contiguous stretches closed by one helix, or several adjacent ones together) of at most
   MAX_DOMAIN nt in the 2D of the whole chain (split());
2. predicts each domain as a chain of its own, with its own template search, 2D and trajectories
   (predict.predict on the subsequence; written to <run>/domains/);
3. turns the best model of every domain into one part of a multi-part template for the whole chain (assembly()):
   C4'-C4' distances and backbone torsions inside each domain, atom positions for the Cartesian refinement. The
   parts are not restrained to each other: their relative placement, the linkers and the multiloops between them
   are sampled by the full-chain search, exactly as for a template that covers only part of a chain.

The whole-chain hypotheses of the normal pipeline keep a share of the trajectories, so the assembly competes
with them in the ranking instead of replacing them.

Status: experimental and off by default (PredictConfig.domains = 'auto' or --domains auto). On six library chains
of 255-393 nt without templates (benchmark/long_rna.py) the assembly was worse than the whole-chain search (top
model TM-score 0.091 against 0.104, valid trajectories 2.5 against 3.3 of 16, 1.8x the time); at that length both
are near the level of unrelated structures.
"""
from __future__ import annotations

import math

import numpy as np

MIN_LENGTH = 240            # chains at least this long are split (PredictConfig.domains = 'auto')
MAX_DOMAIN = 160            # a domain longer than this is split further at its multiloop
MIN_DOMAIN = 25             # shorter stretches are left to the assembly search
TEMPLATE_ATOMS = ['P', "C4'", "C1'", 'NG', 'C2', 'C4', 'C6']


def _children(lo, hi, partner):
    """Helices directly enclosed in [lo, hi): [(i, j)]."""
    out, k = [], lo
    while k < hi:
        l = partner.get(k)
        if l is not None and l > k and l < hi:
            out.append((k, l))
            k = l + 1
        else:
            k += 1
    return out


def split(n, pairs, max_len=MAX_DOMAIN, min_len=MIN_DOMAIN):
    """Domains [(start, end)] (0-based, end exclusive) of a nested 2D, in sequence order."""
    partner = {}
    for i, j in pairs:
        partner[i], partner[j] = j, i
    out = []

    def place(lo, hi):
        kids = _children(lo, hi, partner)
        group = []

        def flush():
            if group:
                s, e = group[0][0], group[-1][1] + 1
                if e - s >= min_len:
                    out.append((s, e))
                group.clear()
        for i, j in kids:
            span = j - i + 1
            if span > max_len:
                flush()
                # descend through the helix and its internal loops to the first multiloop, split its branches
                a, b = i, j
                while True:
                    inner = _children(a + 1, b, partner)
                    if len(inner) != 1:
                        break
                    a, b = inner[0]
                if len(inner) >= 2:
                    place(a + 1, b)
                elif b - a + 1 >= min_len:
                    out.append((i, j + 1))              # a long rod without a multiloop: one domain
                continue
            if group and j + 1 - group[0][0] > max_len:
                flush()
            group.append((i, j))
        flush()
    place(0, n)
    return sorted(out)


def torsions_of(residues):
    """Backbone and glycosidic torsions [alpha, beta, gamma, delta, epsilon, zeta, chi] per residue (degrees,
    None where an atom is missing) from all-atom residues."""
    from .geomutil import dihedral
    from .build import GLY
    out = []
    n = len(residues)
    for r in range(n):
        R = residues[r]
        P_ = residues[r - 1] if r > 0 else {}
        N_ = residues[r + 1] if r + 1 < n else {}
        base = 'G' if 'N9' in R and 'O6' in R else 'A' if 'N9' in R else 'U' if 'O4' in R else 'C'
        quads = [(P_, "O3'", R, 'P', R, "O5'", R, "C5'"), (R, 'P', R, "O5'", R, "C5'", R, "C4'"),
                 (R, "O5'", R, "C5'", R, "C4'", R, "C3'"), (R, "C5'", R, "C4'", R, "C3'", R, "O3'"),
                 (R, "C4'", R, "C3'", R, "O3'", N_, 'P'), (R, "C3'", R, "O3'", N_, 'P', N_, "O5'")]
        row = []
        for q in quads:
            pts = [q[k].get(q[k + 1]) for k in range(0, 8, 2)]
            row.append(None if any(p is None for p in pts) else round(float(dihedral(*map(np.asarray, pts))), 1))
        g = GLY[base]
        pts = [R.get("O4'"), R.get("C1'"), R.get(g), R.get('C4' if g == 'N9' else 'C2')]
        row.append(None if any(p is None for p in pts) else round(float(dihedral(*map(np.asarray, pts))), 1))
        out.append(row)
    return out


def assembly(seq, domains, models, dmax=30.0):
    """A multi-part template hypothesis for the whole chain from domain models.

    domains: [(start, end)]; models: [(residues, nested pairs in domain numbering)] (best model of each domain).
    Returns {'pairs', 'distance_restraints', 'torsions', 'template_atoms', 'covered'}."""
    from .build import GLY
    pairs, torsions, atoms, group = [], {}, {}, {}
    rest = {'i': [], 'j': [], 'd': [], 'tol': []}
    covered = 0
    for g, ((s, e), (res, dpairs)) in enumerate(zip(domains, models)):
        covered += e - s
        pairs += [(s + i, s + j) for i, j in dpairs]
        X = np.array([r["C4'"] for r in res], float)
        D = np.linalg.norm(X[:, None] - X[None], axis=2)
        a, b = np.triu_indices(len(X), 3)
        keep = D[a, b] <= dmax
        a, b = a[keep], b[keep]
        rest['i'] += (a + s).tolist()
        rest['j'] += (b + s).tolist()
        rest['d'] += np.round(D[a, b], 2).tolist()
        rest['tol'] += np.round(0.6 + 0.04 * D[a, b], 2).tolist()
        for k, row in enumerate(torsions_of(res)):
            torsions[s + k] = row
        for k, r in enumerate(res):
            q = s + k
            vals = []
            for nm in TEMPLATE_ATOMS:
                v = r.get(GLY[seq[q]] if nm == 'NG' else nm)
                vals.append(None if v is None else [round(float(c), 3) for c in v])
            atoms[q] = vals
            group[q] = g
    return {'pairs': sorted(pairs), 'distance_restraints': rest, 'torsions': torsions,
            'template_atoms': {'names': TEMPLATE_ATOMS, 'atoms': atoms, 'group': group},
            'covered': covered / len(seq)}


def outer_pairs(n, pairs, domains):
    """Pairs of the whole-chain 2D outside every domain (they hold the domains together)."""
    inside = np.zeros(n, bool)
    for s, e in domains:
        inside[s:e] = True
    return [(i, j) for i, j in pairs if not inside[i] and not inside[j]]


def read_model(path):
    """Residues [{atom: xyz}] of a single-chain model written by structure_io.write_pdb."""
    res = {}
    with open(path, encoding='utf-8') as fh:
        lines = fh.readlines()
    for line in lines:
        if line.startswith('ATOM'):
            k = int(line[22:26])
            res.setdefault(k, {})[line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]),
                                                                   float(line[46:54])])
    return [res[k] for k in sorted(res)]


def fold_domains(seq, base_pairs, config, out_dir, predict, say=None):
    """Predict every domain as its own chain; returns (domains kept, models, report) for assembly()."""
    from dataclasses import replace
    from .annotate import dotbracket_to_pairs
    doms = split(len(seq), base_pairs)
    kept, models, report = [], [], []
    n_traj = config.trajectories or max(8, 4 * config.n_models)
    for k, (s, e) in enumerate(doms):
        cfg = replace(config, n_models=1, trajectories=max(4, n_traj // 2), domains='off', save_candidates=False,
                      export_rejected=False, reactivities=tuple(config.reactivities[s:e]) if config.reactivities else (),
                      alignment=())
        d_out = out_dir / 'domains' / f'domain{k + 1}_{s + 1}-{e}'
        if say:
            say(f'domain {k + 1}/{len(doms)}: residues {s + 1}-{e} ({e - s} nt)')
        man = predict(seq[s:e], d_out, name=f'domain{k + 1}', config=cfg, log=None)
        rep = {'domain': k + 1, 'residues': (s + 1, e), 'folder': str(d_out.relative_to(out_dir)),
               'valid': man['statistics']['valid'], 'trajectories': man['statistics']['trajectories']}
        if man['conformers']:
            c = man['conformers'][0]
            dbn = ''.join(ch if ch in '().' else '.' for ch in c['secondary_structure'])
            kept.append((s, e))
            models.append((read_model(d_out / c['files']['pdb']), [(i - 1, j - 1) for i, j in dotbracket_to_pairs(dbn)]))
            rep.update({'template': c.get('template'), 'score': c['total_score']})
        report.append(rep)
    return kept, models, report
