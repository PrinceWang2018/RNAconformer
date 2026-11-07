"""Template-based modelling from a compact library of experimental RNA structures (no deep learning).

The library (data/pdb_rna_v1.npz, built by tools/build_template_library.py) holds one chain per
non-redundant RNA class of the PDB (BGSU RNA 3D Hub list) with its release date, sequence, canonical base
pairs, backbone torsions and a few atoms per residue.

search()          Smith-Waterman local alignment of the query against every library chain (vectorised
                  over the library, affine gaps), Karlin-Altschul E-values, one hit per class. Hits can be
                  limited to structures released before a date (blind benchmarks) or exclude PDB entries.
hypothesis()      turns a hit into a secondary structure (template pairs that are complementary in the
                  query, nested part as the fold tree, crossing pairs as restraints; unaligned regions
                  folded by ViennaRNA with the template pairs enforced), C4'-C4' distance restraints
                  between aligned residues and starting torsions for aligned unpaired residues.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent / 'data'
LIBRARY_FILE = DATA / 'pdb_rna_v2.npz'
CODE = np.full(256, 4, np.int8)
for _k, _c in enumerate('ACGU'):
    CODE[ord(_c)] = _k
SEP = 5
# scoring as BLASTN (reward 2, penalty -3, gap open 5, extend 2; cost of a gap of length L = 5 + 2 L)
MATCH, MISMATCH, GAP_OPEN, GAP_EXT = 2, -3, 5, 2
LAMBDA, K_KA = 0.625, 0.41
CANON = {'AU', 'UA', 'GC', 'CG', 'GU', 'UG'}


class Library:
    def __init__(self, path=LIBRARY_FILE):
        z = np.load(path)
        self.path = Path(path)
        self.info = json.loads(bytes(z['info']).decode())
        self.meta = self.info['chains']
        self.seq = bytes(z['seq']).decode()
        self.off = z['offsets']
        self.compact = 'c4' in z.files                  # schema 2: C4' int16, other atoms int8 offsets
        if self.compact:
            self.c4_q, self.rel_q, self.tor_q = z['c4'], z['rel'], z['tor8']
            self.rel_scale = float(self.info['rel_scale_per_angstrom'])
            self.tor_scale = float(self.info['tor_scale_per_degree'])
        else:
            self.xyz_q, self.tor_q = z['xyz'], z['tor']
        self.origin = z['origin']
        self.scale = float(self.info['scale_per_angstrom'])
        self.missing = int(self.info['missing'])
        self.atoms = self.info['atoms']
        p = z['pairs']
        self.pair_chain = np.searchsorted(self.off, p[:, 0], side='right') - 1 if len(p) else np.zeros(0, int)
        self.pairs_global = p
        self.codes = CODE[np.frombuffer(self.seq.encode(), np.uint8)]
        self.n_chains = len(self.meta)
        self.class_of = [m['class'] for m in self.meta]

    @property
    def version(self):
        return f"{self.path.name} ({self.info.get('created')}, {self.n_chains} chains)"

    def chain_seq(self, k):
        return self.seq[self.off[k]:self.off[k + 1]]

    def chain_xyz(self, k, atom="C4'"):
        a = self.atoms.index(atom)
        s = slice(self.off[k], self.off[k + 1])
        if not self.compact:
            q = self.xyz_q[s, a].astype(np.float64)
            bad = q[:, 0] == self.missing
            x = q / self.scale + self.origin[k]
            x[bad] = np.nan
            return x
        c = self.c4_q[s].astype(np.float64)
        bad = c[:, 0] == self.missing
        x = c / self.scale + self.origin[k]
        if a != 1:
            r = self.rel_q[s, a if a < 1 else a - 1].astype(np.float64)
            bad = bad | (r[:, 0] == -128)
            x = x + r / self.rel_scale
        x[bad] = np.nan
        return x

    def chain_tor(self, k):
        t = self.tor_q[self.off[k]:self.off[k + 1]].astype(np.float64)
        if self.compact:
            bad = t == -128
            t = t / self.tor_scale
        else:
            bad = t == self.missing
            t = t / 100.0
        t[bad] = np.nan
        return t

    def chain_pairs(self, k):
        a, b = np.searchsorted(self.pair_chain, [k, k + 1])
        o = int(self.off[k])
        return [(int(i) - o, int(j) - o) for i, j in self.pairs_global[a:b]]

    def allowed(self, cutoff=None, exclude=()):
        """Chains usable as templates: released before `cutoff` (YYYY-MM-DD) and not in `exclude` (PDB ids)."""
        ex = {e.upper()[:4] for e in exclude}
        ok = np.ones(self.n_chains, bool)
        for k, m in enumerate(self.meta):
            if ex and m['pdb'].upper() in ex:
                ok[k] = False
            elif cutoff and (not m['release'] or m['release'] >= cutoff):
                ok[k] = False
        return ok


@lru_cache(maxsize=2)
def load_library(path=None):
    p = Path(path) if path else LIBRARY_FILE
    if not p.exists():
        return None
    return Library(p)


# ----------------------------------------------------------------------------------------- alignment
def _score_rows(q):
    """Substitution scores for each query letter against codes 0..5 (A C G U N separator)."""
    m = np.full((5, 6), MISMATCH, np.int64)
    for k in range(4):
        m[k, k] = MATCH
    m[4, :5] = -1
    m[:, 4] = -1
    m[:, SEP] = -10 ** 6
    return m[q]


def sw_scan(query, lib: Library, allowed=None):
    """Best local alignment score of the query against every library chain (int array, -1 = not allowed)."""
    q = CODE[np.frombuffer(query.encode(), np.uint8)]
    n_ch = lib.n_chains
    # database: chains separated by one separator cell; segment index for gap resets
    lens = np.diff(lib.off)
    db = np.full(int(lib.off[-1]) + n_ch, SEP, np.int8)
    seg = np.zeros(len(db), np.int64)
    pos = lib.off[:-1] + np.arange(n_ch)            # start of each chain in db (after its separator)
    for k in range(n_ch):
        db[pos[k] + 1:pos[k] + 1 + lens[k]] = lib.codes[lib.off[k]:lib.off[k + 1]]
        seg[pos[k]:pos[k] + 1 + lens[k]] = k
    M = len(db)
    rows = _score_rows(q)
    H = np.zeros(M + 1, np.int64)
    F = np.full(M + 1, -10 ** 9, np.int64)
    best = np.zeros(M + 1, np.int64)
    ar = np.arange(M + 1, dtype=np.int64) * GAP_EXT
    big = np.concatenate([[0], seg]) * np.int64(10 ** 10)
    for i in range(len(q)):
        s = rows[i][db]
        diag = H[:-1] + s
        F = np.maximum(H - GAP_OPEN - GAP_EXT, F - GAP_EXT)
        H0 = np.zeros(M + 1, np.int64)
        H0[1:] = np.maximum(diag, F[1:])
        np.maximum(H0, 0, out=H0)
        A = H0 + ar + big
        Pm = np.maximum.accumulate(A)
        E = np.full(M + 1, -10 ** 9, np.int64)
        E[1:] = Pm[:-1] - big[1:] - ar[1:] - GAP_OPEN
        H = np.maximum(H0, E)
        H[1:][db == SEP] = 0
        np.maximum(best, H, out=best)
    per = np.maximum.reduceat(best[1:], pos)
    if allowed is not None:
        per = np.where(allowed, per, -1)
    return per


def align(query, target):
    """Smith-Waterman (Gotoh) alignment with traceback. Returns (score, [(qi, ti), ...] aligned positions)."""
    q = CODE[np.frombuffer(query.encode(), np.uint8)]
    t = CODE[np.frombuffer(target.encode(), np.uint8)]
    n, m = len(q), len(t)
    S = _score_rows(q)[:, t]                         # n x m
    NEG = -10 ** 9
    H = np.zeros((n + 1, m + 1), np.int64)
    E = np.full((n + 1, m + 1), NEG, np.int64)       # gap in query (consume target)
    F = np.full((n + 1, m + 1), NEG, np.int64)       # gap in target (consume query)
    ar = np.arange(m + 1, dtype=np.int64) * GAP_EXT
    for i in range(1, n + 1):
        F[i] = np.maximum(H[i - 1] - GAP_OPEN - GAP_EXT, F[i - 1] - GAP_EXT)
        H0 = np.zeros(m + 1, np.int64)
        H0[1:] = np.maximum(H[i - 1, :-1] + S[i - 1], F[i, 1:])
        np.maximum(H0, 0, out=H0)
        Pm = np.maximum.accumulate(H0 + ar)
        E[i, 1:] = Pm[:-1] - ar[1:] - GAP_OPEN
        H[i] = np.maximum(H0, E[i])
    i, j = np.unravel_index(int(np.argmax(H)), H.shape)
    score = int(H[i, j])
    out = []
    state = 'H'
    while i > 0 and j > 0:
        if state == 'H':
            h = H[i, j]
            if h == 0:
                break
            if h == H[i - 1, j - 1] + S[i - 1, j - 1]:
                out.append((i - 1, j - 1))
                i, j = i - 1, j - 1
            elif h == E[i, j]:
                state = 'E'
            elif h == F[i, j]:
                state = 'F'
            else:
                break
        elif state == 'E':                           # E[i, j] = max_k H[i, k] - open - ext (j - k)
            e = E[i, j]
            if e == H[i, j - 1] - GAP_OPEN - GAP_EXT:
                state = 'H'
            j -= 1
        else:
            f = F[i, j]
            if f == H[i - 1, j] - GAP_OPEN - GAP_EXT:
                state = 'H'
            i -= 1
    return score, out[::-1]


def align_segments(query, target, min_score=24, depth=3):
    """Collinear local alignments (like BLAST's several HSPs): the best local alignment plus further alignments
    of the remaining query/target regions before and after it (e.g. both halves of a scaffold around an insert).
    Returns (best score, aligned positions sorted by query index)."""
    score, aln = align(query, target)
    if not aln or depth <= 0:
        return score, aln
    (q0, t0), (q1, t1) = aln[0], aln[-1]
    out = list(aln)
    if q0 >= 8 and t0 >= 8:
        s, left = align_segments(query[:q0], target[:t0], min_score, depth - 1)
        if s >= min_score:
            out = left + out
    if len(query) - q1 - 1 >= 8 and len(target) - t1 - 1 >= 8:
        s, right = align_segments(query[q1 + 1:], target[t1 + 1:], min_score, depth - 1)
        if s >= min_score:
            out += [(a + q1 + 1, b + t1 + 1) for a, b in right]
    return score, out


def evalue(score, m, n_db):
    return K_KA * m * n_db * math.exp(-LAMBDA * score)


def bits(score):
    return (LAMBDA * score - math.log(K_KA)) / math.log(2)


def homologs(query, lib: Library, scan=None, max_evalue=1e-2):
    """Boolean mask of library chains with a significant local sequence match to the query."""
    per = sw_scan(query, lib) if scan is None else scan
    n_db = int(lib.off[-1])
    smin = (math.log(K_KA * len(query) * n_db) - math.log(max_evalue)) / LAMBDA
    return per >= smin


FULL_LENGTH_EVALUE, FULL_LENGTH_COVERAGE = 10.0, 0.6


def search(query, lib: Library | None = None, cutoff=None, exclude=(), max_hits=5, max_evalue=1e-3,
           min_coverage=0.45, scan=None):
    """Significant template hits (best first), one per non-redundant class.

    A hit needs E <= max_evalue, or E <= 10 x max_evalue when it covers >= 60% of the query (a full-length
    homolog with scattered differences, e.g. a tRNA of another species), and must cover >= 45% of the query
    (a template that covers less does not define the fold; parts are added by segment_hits)."""
    lib = lib or load_library()
    if lib is None:
        return []
    allowed = lib.allowed(cutoff, exclude)
    per = np.where(allowed, sw_scan(query, lib) if scan is None else scan, -1)
    n_db = int(lib.off[-1])
    m = len(query)
    order = np.argsort(-per)
    hits, seen = [], set()
    for k in order[:200]:
        k = int(k)
        if per[k] <= 0:
            break
        e = evalue(int(per[k]), m, n_db)
        if e > FULL_LENGTH_EVALUE * max_evalue:
            break
        cls = lib.class_of[k]
        if cls in seen:
            continue
        seen.add(cls)
        tseq = lib.chain_seq(k)
        score, aln = align_segments(query, tseq)
        if not aln:
            continue
        xyz = lib.chain_xyz(k)
        aln = [(qi, ti) for qi, ti in aln if not np.isnan(xyz[ti, 0])]
        cov = len(aln) / m
        if cov < min_coverage or (e > max_evalue and cov < FULL_LENGTH_COVERAGE):
            continue
        ident = sum(query[qi] == tseq[ti] for qi, ti in aln) / max(1, len(aln))
        tpairs = lib.chain_pairs(k)
        tmap = {ti: qi for qi, ti in aln}
        mapped = [(tmap[a], tmap[b]) for a, b in tpairs if a in tmap and b in tmap]
        compat = sum(query[i] + query[j] in CANON for i, j in mapped)
        meta = lib.meta[k]
        hits.append({'chain_index': k, 'pdb': meta['pdb'], 'chain': meta['chain'], 'model': meta['model'],
                     'class': cls, 'release': meta['release'], 'method': meta['method'],
                     'resolution': meta['resolution'], 'title': meta['title'], 'template_length': meta['length'],
                     'score': score, 'bits': round(bits(score), 1), 'evalue': float(f'{evalue(score, m, n_db):.3g}'),
                     'identity': round(ident, 3), 'coverage': round(cov, 3), 'aligned': aln,
                     'pairs_mapped': len(mapped), 'pairs_compatible': compat})
        if len(hits) >= max_hits:
            break
    hits.sort(key=lambda h: (h['evalue'], -h['coverage']))
    return hits


def _similarity(lib, h1, h2, L):
    """TM-score-like agreement (0-1) of two hits on the query residues both align (C4' atoms), or None."""
    X1, X2 = lib.chain_xyz(h1['chain_index']), lib.chain_xyz(h2['chain_index'])
    m1, m2 = dict(h1['aligned']), dict(h2['aligned'])
    com = [q for q in m1 if q in m2 and not np.isnan(X1[m1[q], 0]) and not np.isnan(X2[m2[q], 0])]
    if len(com) < 15:
        return None
    A = np.array([X1[m1[q]] for q in com])
    B = np.array([X2[m2[q]] for q in com])
    d0 = max(0.5, 1.24 * (max(L, 19) - 15) ** (1 / 3) - 1.8)
    w = np.ones(len(com))
    for _ in range(8):                      # weighted superposition, weights as in the TM-score
        a0 = (A * w[:, None]).sum(0) / w.sum()
        b0 = (B * w[:, None]).sum(0) / w.sum()
        u, _s, vt = np.linalg.svd(((A - a0) * w[:, None]).T @ (B - b0))
        R = u @ np.diag([1.0, 1.0, np.sign(np.linalg.det(u @ vt))]) @ vt
        d = np.linalg.norm((A - a0) @ R - (B - b0), axis=1)
        w = 1.0 / (1.0 + (d / d0) ** 2)
    return float(w.mean())


def rank_templates(query, hits, lib, n=3):
    """Order hits by identity x coverage x agreement with the other hits (3D consensus) and keep n.

    A template in an unusual conformation (e.g. a tRNA bent by a bound enzyme) disagrees with the other
    hits of the same family and moves down; with a single hit the agreement is 1."""
    L = len(query)
    for h in hits:
        sims = [x for x in (_similarity(lib, h, o, L) for o in hits if o is not h) if x is not None]
        h['consensus'] = round(float(np.mean(sims)), 3) if sims else 1.0
        h['rank_score'] = round(h['identity'] * h['coverage'] * (0.5 + 0.5 * h['consensus']), 3)
    return sorted(hits, key=lambda h: -h['rank_score'])[:n]


# ------------------------------------------------------------------------------- template -> model
def _drop_short_hairpins(pairs, min_loop=3):
    pairs = sorted(pairs)
    out = []
    for i, j in pairs:
        inner = [p for p in pairs if i < p[0] and p[1] < j]
        if not inner and j - i - 1 < min_loop:
            continue
        out.append((i, j))
    return out


def hypothesis(query, hit, lib: Library | None = None, temperature=37.0, dmax=30.0, extra=()):
    """Secondary structure, crossing-pair restraints, distance restraints and torsions from one hit.

    extra: further hits for parts of the query the first hit does not cover (alignments already in query
    numbering, see segment_hits). Each part contributes its own base pairs, restraints, torsions and atoms;
    parts are not restrained to each other (their relative placement is left to the physics model)."""
    from .annotate import nested_subset, pairs_to_dotbracket
    lib = lib or load_library()
    parts = []
    covered = set()
    for h in (hit, *extra):
        aln = [(int(a), int(b)) for a, b in h['aligned'] if int(a) not in covered]
        covered.update(a for a, _ in aln)
        if aln:
            parts.append((int(h['chain_index']), aln))
    pairs = []
    used = set()
    for k, aln in parts:
        tmap = {ti: qi for qi, ti in aln}
        for a, b in lib.chain_pairs(k):
            a, b = int(a), int(b)
            if a in tmap and b in tmap:
                i, j = sorted((tmap[a], tmap[b]))
                if query[i] + query[j] in CANON and i not in used and j not in used and j - i > 3:
                    pairs.append((i, j))
                    used.update((i, j))
    nested, crossing = nested_subset(pairs)
    nested = _drop_short_hairpins(nested)
    n = len(query)
    energy = None
    try:                                   # fold the residues the templates do not cover, template pairs enforced
        import RNA
        md = RNA.md()
        md.temperature = temperature
        md.noLP = 0
        fc = RNA.fold_compound(query, md)
        for i, j in nested:
            fc.hc_add_bp(i + 1, j + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS | RNA.CONSTRAINT_CONTEXT_ENFORCE)
        for qi in covered - {x for p in nested for x in p}:
            fc.hc_add_up(qi + 1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
        s, e = fc.mfe()
        from .annotate import dotbracket_to_pairs
        folded = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(s)]
        if set(nested) <= set(folded):
            nested = _drop_short_hairpins(sorted(folded))
            energy = round(e, 2)
    except Exception:
        pass
    taken = {x for p in nested for x in p}
    crossing = [(i, j) for i, j in crossing if i not in taken and j not in taken]
    restraints = {'i': [], 'j': [], 'd': [], 'tol': []}
    torsions, atoms, group = {}, {}, {}
    base_atoms = {'NG', 'C2', 'C4', 'C6'}
    for g, (k, aln) in enumerate(parts):
        # distance restraints between aligned residues of this part (C4' atoms)
        xyz = lib.chain_xyz(k)
        qi = np.array([a for a, _ in aln])
        ti = np.array([b for _, b in aln])
        X = xyz[ti]
        D = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        a_idx, b_idx = np.triu_indices(len(qi), 1)
        sep = np.abs(qi[b_idx] - qi[a_idx])
        keep = (sep >= 3) & (D[a_idx, b_idx] <= dmax)
        a_idx, b_idx = a_idx[keep], b_idx[keep]
        d = D[a_idx, b_idx]
        restraints['i'] += qi[a_idx].tolist()
        restraints['j'] += qi[b_idx].tolist()
        restraints['d'] += np.round(d, 2).tolist()
        restraints['tol'] += np.round(0.6 + 0.04 * d, 2).tolist()
        tor = lib.chain_tor(k)
        for a, b in aln:
            torsions[a] = [None if np.isnan(v) else round(float(v), 1) for v in tor[b]]
        # template atom positions for the Cartesian refinement (base atoms only if purine/pyrimidine match)
        tseq = lib.chain_seq(k)
        allx = {nm: lib.chain_xyz(k, nm) for nm in lib.atoms}
        for a, b in aln:
            same_class = (query[a] in 'AG') == (tseq[b] in 'AG') and tseq[b] in 'ACGU'
            vals = []
            for nm in lib.atoms:
                v = allx[nm][b]
                ok = not np.isnan(v[0]) and (same_class or nm not in base_atoms)
                vals.append([round(float(c), 3) for c in v] if ok else None)
            if any(v is not None for v in vals):
                atoms[a] = vals
                group[a] = g
    dbn = pairs_to_dotbracket(n, nested, offset=0)
    return {'dbn': dbn, 'pairs': nested, 'restraints': crossing, 'energy_kcal_mol': energy,
            'distance_restraints': restraints, 'torsions': torsions,
            'template_atoms': {'names': list(lib.atoms), 'atoms': atoms, 'group': group}}


def uncovered_segments(n, aligned, min_len=20, max_gap=3):
    """Maximal runs of query residues (0-based [start, end)) not covered by an alignment; runs separated by
    at most max_gap covered residues are merged."""
    cov = np.zeros(n, bool)
    for a, _ in aligned:
        cov[int(a)] = True
    segs, k = [], 0
    while k < n:
        if cov[k]:
            k += 1
            continue
        s = k
        while k < n and not cov[k]:
            k += 1
        segs.append([s, k])
    merged = []
    for s, e in segs:
        if merged and s - merged[-1][1] <= max_gap:
            merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged if e - s >= min_len]


def segment_hits(query, hit, lib, cutoff=None, exclude=(), min_len=20, max_evalue=0.05, min_identity=0.85,
                 min_aligned=20):
    """Templates for the parts of the query that `hit` leaves uncovered (e.g. an insert in a scaffold):
    one ranked hit per uncovered segment of >= min_len nt, alignments shifted to query numbering.

    A short segment cannot reach a strict E-value against the whole library, so a segment hit needs
    E <= 0.05 together with >= 85% identity over >= 20 aligned residues."""
    out = []
    for s, e in uncovered_segments(len(query), hit['aligned'], min_len):
        sub = query[s:e]
        found = search(sub, lib, cutoff=cutoff, exclude=exclude, max_hits=8, max_evalue=max_evalue, min_coverage=0.4)
        found = [h for h in found if h['class'] != hit['class'] and h['identity'] >= min_identity
                 and len(h['aligned']) >= min_aligned]
        if not found:
            continue
        best = rank_templates(sub, found, lib, 1)[0]
        best = {**best, 'aligned': [(int(a) + s, int(b)) for a, b in best['aligned']], 'segment': (s + 1, e)}
        out.append(best)
    return out
