"""Derive local RNA covalent geometry from a high-resolution calibration set.

Output: src/rnafold3d/data/geometry_v1.json (versioned, with provenance).

The calibration set must be disjoint from every benchmark target and family. It
only provides *single-nucleotide* covalent geometry, canonical base-pair
geometry and single-residue backbone torsion statistics - never tertiary
fragments or motifs.

Usage:
    PYTHONPATH=src python tools/calibrate_geometry.py --raw benchmark/data/raw \
        --candidates benchmark/data/calibration_candidates.txt \
        --exclude benchmark/manifest_v1.json --out src/rnafold3d/data/geometry_v1.json
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

from rnafold3d.structure_io import read_mmcif
from rnafold3d.annotate import canonical_pairs
from rnafold3d.geomutil import dihedral, angle, kabsch, pseudorotation, circular_mean

EXCLUDE_TITLE_WORDS = ('riboswitch', 'aptamer', 'ribozyme', 'sarcin', 'ricin', 'trna', 'transfer rna',
                       'guanidine', 'fluoride', 'pseudoknot', 'tar ', 'hairpin', 'intron', 'ribonuclease p',
                       'srl', 'xrrna', 'tetraloop')
SUGAR = ("C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", 'N')
RING = ("C4'", "O4'", "C1'", "C2'", "C3'")
BASE = {
    'A': ("C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'),
    'G': ("C1'", 'N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'),
    'C': ("C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'),
    'U': ("C1'", 'N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'),
}
BASE_FIT = {b: tuple(a for a in names if a != "C1'") for b, names in BASE.items()}
GLY = {'A': 'N9', 'G': 'N9', 'C': 'N1', 'U': 'N1'}
CHI_REF = {'A': 'C4', 'G': 'C4', 'C': 'C2', 'U': 'C2'}
BACKBONE = ('P', 'OP1', 'OP2', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'")


def complete(res):
    return all(a in res.atoms for a in BACKBONE[3:]) and all(a in res.atoms for a in BASE[res.base])


def sugar_coords(res):
    at = dict(res.atoms)
    at['N'] = at[GLY[res.base]]
    return np.array([at[a] for a in SUGAR])


def iterative_average(blocks, fit_idx, rounds=4):
    ref = blocks[0]
    for _ in range(rounds):
        aligned = []
        for x in blocks:
            r, t = kabsch(x[fit_idx], ref[fit_idx])
            aligned.append(x @ r.T + t)
        aligned = np.array(aligned)
        ref = aligned.mean(axis=0)
    spread = np.sqrt(((aligned - ref) ** 2).sum(axis=2).mean(axis=0))
    return ref, spread


def planarize(coords):
    c = coords.mean(axis=0)
    _, _, vt = np.linalg.svd(coords - c)
    n = vt[2]
    return coords - np.outer((coords - c) @ n, n)


def residue_frame(res, base_tpl):
    """Rigid transform (R, t) mapping the base template into the observed base."""
    names = BASE_FIT[res.base]
    tpl = np.array([base_tpl[res.base][a] for a in names])
    obs = np.array([res.atoms[a] for a in names])
    return kabsch(tpl, obs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--raw', default='benchmark/data/raw')
    ap.add_argument('--candidates', default='benchmark/data/calibration_candidates.txt')
    ap.add_argument('--exclude', default='benchmark/manifest_v1.json')
    ap.add_argument('--max-chain', type=int, default=50)
    ap.add_argument('--out', default='src/rnafold3d/data/geometry_v1.json')
    a = ap.parse_args()

    excluded = set()
    if Path(a.exclude).exists():
        excluded = {t['pdb_id'].upper() for t in json.loads(Path(a.exclude).read_text())['targets']}
    ids = [x.strip().upper() for x in Path(a.candidates).read_text().split() if x.strip()]
    used, skipped = [], {}
    chains_all = []
    entries_all = []
    for pid in ids:
        path = Path(a.raw) / f'{pid}.cif'
        if pid in excluded:
            skipped[pid] = 'benchmark_target'
            continue
        if not path.exists():
            skipped[pid] = 'missing_file'
            continue
        chains, meta = read_mmcif(path)
        title = (meta.get('title') or '').lower() + ' '
        word = next((w for w in EXCLUDE_TITLE_WORDS if w in title), None)
        if word:
            skipped[pid] = f'family_keyword:{word.strip()}'
            continue
        if not chains or max(len(c.sequence) for c in chains.values()) > a.max_chain:
            skipped[pid] = 'chain_too_long_or_absent'
            continue
        used.append({'pdb_id': pid, 'resolution': meta.get('resolution'), 'title': meta.get('title', '')[:120],
                     'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        chains_all.extend(chains.values())
        entries_all.append(list(chains.values()))
    print(f'calibration entries used: {len(used)}; skipped: {len(skipped)}')

    # ---------------------------------------------------------------- per-residue measurements
    records = []        # (res, prev, nxt) for complete canonical residues
    for ch in chains_all:
        rmap = ch.residue_map()
        for r in ch.residues:
            if r.modified or not complete(r):
                continue
            prev = rmap.get(r.index - 1)
            nxt = rmap.get(r.index + 1)
            if prev is not None and ("O3'" not in prev.atoms or "C3'" not in prev.atoms or 'P' not in r.atoms or
                                     not 1.45 < np.linalg.norm(prev.atoms["O3'"] - r.atoms['P']) < 1.75):
                prev = None
            if nxt is not None and ('P' not in nxt.atoms or "O5'" not in nxt.atoms or
                                    not 1.45 < np.linalg.norm(r.atoms["O3'"] - nxt.atoms['P']) < 1.75):
                nxt = None
            records.append((r, prev, nxt))
    print('complete residues:', len(records))

    def tors(r, prev, nxt):
        A = r.atoms
        t = {}
        if prev is not None:
            t['alpha'] = dihedral(prev.atoms["O3'"], A['P'], A["O5'"], A["C5'"])
        if 'P' in A:
            t['beta'] = dihedral(A['P'], A["O5'"], A["C5'"], A["C4'"])
        t['gamma'] = dihedral(A["O5'"], A["C5'"], A["C4'"], A["C3'"])
        t['delta'] = dihedral(A["C5'"], A["C4'"], A["C3'"], A["O3'"])
        if nxt is not None:
            t['epsilon'] = dihedral(A["C4'"], A["C3'"], A["O3'"], nxt.atoms['P'])
            t['zeta'] = dihedral(A["C3'"], A["O3'"], nxt.atoms['P'], nxt.atoms["O5'"])
        t['chi'] = dihedral(A["O4'"], A["C1'"], A[GLY[r.base]], A[CHI_REF[r.base]])
        ph, amp = pseudorotation([A[x] for x in ("C4'", "O4'", "C1'", "C2'", "C3'")])
        t['phase'] = ph
        t['amplitude'] = amp
        return t

    meas = [(r, p, n, tors(r, p, n)) for r, p, n in records]

    def pucker_class(t):
        ph = t['phase'] % 360
        if ph < 50 or ph > 340:
            return 'C3endo'
        if 130 <= ph <= 200:
            return 'C2endo'
        return None

    # ---------------------------------------------------------------- bonds and angles
    bond_defs = {
        'P-OP1': ('P', 'OP1'), 'P-OP2': ('P', 'OP2'), "P-O5'": ('P', "O5'"), "O5'-C5'": ("O5'", "C5'"),
        "C5'-C4'": ("C5'", "C4'"), "C4'-O4'": ("C4'", "O4'"), "O4'-C1'": ("O4'", "C1'"), "C1'-C2'": ("C1'", "C2'"),
        "C2'-C3'": ("C2'", "C3'"), "C3'-C4'": ("C3'", "C4'"), "C3'-O3'": ("C3'", "O3'"), "C2'-O2'": ("C2'", "O2'"),
    }
    angle_defs = {
        "OP1-P-OP2": ('OP1', 'P', 'OP2'), "OP1-P-O5'": ('OP1', 'P', "O5'"), "OP2-P-O5'": ('OP2', 'P', "O5'"),
        "P-O5'-C5'": ('P', "O5'", "C5'"), "O5'-C5'-C4'": ("O5'", "C5'", "C4'"), "C5'-C4'-C3'": ("C5'", "C4'", "C3'"),
        "C5'-C4'-O4'": ("C5'", "C4'", "O4'"), "C4'-C3'-O3'": ("C4'", "C3'", "O3'"), "C2'-C3'-O3'": ("C2'", "C3'", "O3'"),
        "C4'-O4'-C1'": ("C4'", "O4'", "C1'"), "O4'-C1'-C2'": ("O4'", "C1'", "C2'"), "C1'-C2'-C3'": ("C1'", "C2'", "C3'"),
        "C2'-C3'-C4'": ("C2'", "C3'", "C4'"), "C3'-C4'-O4'": ("C3'", "C4'", "O4'"), "C1'-C2'-O2'": ("C1'", "C2'", "O2'"),
        "C3'-C2'-O2'": ("C3'", "C2'", "O2'"),
    }
    bonds = {k: [] for k in bond_defs}
    bonds["O3'-P"] = []
    bonds["C1'-N"] = []
    angles = {k: [] for k in angle_defs}
    angles["C3'-O3'-P"] = []
    angles["O3'-P-O5'"] = []
    angles["O3'-P-OP1"] = []
    angles["O3'-P-OP2"] = []
    angles["O4'-C1'-N"] = []
    angles["C2'-C1'-N"] = []
    for r, prev, nxt, t in meas:
        A = r.atoms
        for k, (x, y) in bond_defs.items():
            if x in A and y in A:
                bonds[k].append(np.linalg.norm(A[x] - A[y]))
        for k, (x, y, z) in angle_defs.items():
            if all(q in A for q in (x, y, z)):
                angles[k].append(angle(A[x], A[y], A[z]))
        bonds["C1'-N"].append(np.linalg.norm(A["C1'"] - A[GLY[r.base]]))
        angles["O4'-C1'-N"].append(angle(A["O4'"], A["C1'"], A[GLY[r.base]]))
        angles["C2'-C1'-N"].append(angle(A["C2'"], A["C1'"], A[GLY[r.base]]))
        if prev is not None:
            bonds["O3'-P"].append(np.linalg.norm(prev.atoms["O3'"] - A['P']))
            angles["C3'-O3'-P"].append(angle(prev.atoms["C3'"], prev.atoms["O3'"], A['P']))
            angles["O3'-P-O5'"].append(angle(prev.atoms["O3'"], A['P'], A["O5'"]))
            if 'OP1' in A and 'OP2' in A:
                angles["O3'-P-OP1"].append(angle(prev.atoms["O3'"], A['P'], A['OP1']))
                angles["O3'-P-OP2"].append(angle(prev.atoms["O3'"], A['P'], A['OP2']))

    def stat(v):
        v = np.asarray(v, float)
        med = np.median(v)
        mad = np.median(np.abs(v - med)) * 1.4826
        core = v[np.abs(v - med) <= 5 * max(mad, 1e-3)]   # drop gross outliers / errors
        return {'mean': round(float(core.mean()), 4), 'sd': round(float(core.std()), 4), 'n': int(len(core)),
                'p001': round(float(np.quantile(core, 0.001)), 4), 'p999': round(float(np.quantile(core, 0.999)), 4)}

    bond_stats = {k: stat(v) for k, v in bonds.items() if v}
    angle_stats = {k: stat(v) for k, v in angles.items() if v}

    # ---------------------------------------------------------------- templates
    by_pucker = {'C3endo': [], 'C2endo': []}
    for r, p, n, t in meas:
        c = pucker_class(t)
        if c:
            by_pucker[c].append(sugar_coords(r))
    fit_idx = [SUGAR.index(x) for x in RING]
    sugar_templates = {}
    for k, blocks in by_pucker.items():
        avg, spread = iterative_average(np.array(blocks), fit_idx)
        avg = avg - avg[SUGAR.index("C4'")]
        sugar_templates[k] = {'atoms': {nm: [round(float(v), 4) for v in xyz] for nm, xyz in zip(SUGAR, avg)},
                              'rms_spread': {nm: round(float(s), 3) for nm, s in zip(SUGAR, spread)},
                              'n': len(blocks)}
    base_templates = {}
    for b in 'ACGU':
        blocks = np.array([np.array([r.atoms[x] for x in BASE[b]]) for r, p, n, t in meas if r.base == b])
        fidx = list(range(1, len(BASE[b])))
        avg, spread = iterative_average(blocks, fidx)
        base_part = planarize(avg[1:])
        avg = np.vstack([avg[:1], base_part])
        avg = avg - avg[1]   # glycosidic N at origin
        base_templates[b] = {'atoms': {nm: [round(float(v), 4) for v in xyz] for nm, xyz in zip(BASE[b], avg)},
                             'rms_spread': {nm: round(float(s), 3) for nm, s in zip(BASE[b], spread)}, 'n': len(blocks)}

    # phosphate template in frame (O3'(i-1), P, O5')
    ph_blocks = []
    for r, prev, nxt, t in meas:
        if prev is not None and 'OP1' in r.atoms and 'OP2' in r.atoms:
            ph_blocks.append(np.array([prev.atoms["O3'"], r.atoms['P'], r.atoms["O5'"], r.atoms['OP1'], r.atoms['OP2']]))
    ph_avg, ph_spread = iterative_average(np.array(ph_blocks), [0, 1, 2])
    ph_avg = ph_avg - ph_avg[1]
    phosphate_template = {'atoms': {nm: [round(float(v), 4) for v in xyz] for nm, xyz in
                                    zip(("O3'prev", 'P', "O5'", 'OP1', 'OP2'), ph_avg)},
                          'rms_spread': [round(float(s), 3) for s in ph_spread], 'n': len(ph_blocks)}

    # ---------------------------------------------------------------- pair templates + helix torsions
    btpl = {b: {k: np.array(v) for k, v in base_templates[b]['atoms'].items()} for b in 'ACGU'}
    pair_tf = {}
    helix_torsions = []
    loop_conformers = []
    # Annotate per entry across chains (calibration duplexes are often two strands).
    class _R:
        __slots__ = ('index', 'base', 'atoms', 'src')

        def __init__(self, index, src):
            self.index = index
            self.base = src.base
            self.atoms = src.atoms
            self.src = src

    tmap = {id(r): t for r, p, n, t in meas}
    for entry_chains in entries_all:
        glob = []
        for k, ch in enumerate(entry_chains):
            for r in ch.residues:
                if not r.modified and complete(r):
                    glob.append(_R(k * 100000 + r.index, r))
        gmap = {g.index: g for g in glob}
        pairs = canonical_pairs(glob)
        partner = {}
        for i, j, typ in pairs:
            partner[i] = j
            partner[j] = i
            ri, rj = gmap[i], gmap[j]
            Ri, ti = residue_frame(ri, btpl)
            Rj, tj = residue_frame(rj, btpl)
            pair_tf.setdefault(ri.base + rj.base, []).append((Ri.T @ Rj, Ri.T @ (tj - ti)))
            pair_tf.setdefault(rj.base + ri.base, []).append((Rj.T @ Ri, Rj.T @ (ti - tj)))
        for g in glob:
            t = tmap.get(id(g.src))
            if t is None or not all(k in t for k in ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi')):
                continue
            i = g.index
            j = partner.get(i)
            stacked_helix = (j is not None and partner.get(i + 1) == j - 1 and partner.get(i - 1) == j + 1)
            vec = [t[k] for k in ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi')] + [t['phase']]
            if stacked_helix and pucker_class(t) == 'C3endo':
                helix_torsions.append(vec)
            elif j is None:
                loop_conformers.append(vec)

    def avg_tf(items):
        Rs = np.array([x[0] for x in items])
        ts = np.array([x[1] for x in items])
        u, _, vt = np.linalg.svd(Rs.sum(axis=0))
        d = np.sign(np.linalg.det(u @ vt))
        R = u @ np.diag([1, 1, d]) @ vt
        tmed = np.median(ts, axis=0)
        keep = np.linalg.norm(ts - tmed, axis=1) < 1.5
        return R, ts[keep].mean(axis=0), int(keep.sum())

    pair_templates = {}
    for k, items in sorted(pair_tf.items()):
        if len(items) < 5:
            continue
        R, t, n = avg_tf(items)
        pair_templates[k] = {'R': [[round(float(v), 5) for v in row] for row in R],
                             't': [round(float(v), 4) for v in t], 'n': n}

    hel = np.array(helix_torsions)
    names = ('alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi')
    aform = {nm: round(float(circular_mean(hel[:, k])), 2) for k, nm in enumerate(names)}
    aform_sd = {nm: round(float(np.degrees(np.sqrt(-2 * np.log(np.abs(np.mean(np.exp(1j * np.radians(hel[:, k])))))))), 2)
                for k, nm in enumerate(names)}

    # torsion histograms (all residues with full torsions), 10-degree bins
    hist = {}
    for k, nm in enumerate(names):
        vals = [m[3][nm] for m in meas if nm in m[3]]
        h, _ = np.histogram(np.asarray(vals), bins=36, range=(-180, 180))
        hist[nm] = [int(x) for x in h]

    rng = np.random.default_rng(0)
    loop_conformers = np.array(loop_conformers)
    lib = loop_conformers
    if len(lib) > 1500:
        lib = lib[rng.choice(len(lib), 1500, replace=False)]

    out = {
        'schema': 'rnafold3d.geometry/1',
        'version': '1.0.0',
        'created': datetime.date.today().isoformat(),
        'description': 'Local nucleotide covalent geometry, canonical pair geometry and single-residue torsion '
                       'statistics derived from high-resolution RNA-only X-ray structures disjoint from the '
                       'benchmark set. Contains no tertiary fragments.',
        'provenance': {'entries': used, 'n_entries': len(used), 'skipped': skipped,
                       'filters': {'max_chain_length': a.max_chain, 'resolution_max': 1.6,
                                   'excluded_title_keywords': EXCLUDE_TITLE_WORDS,
                                   'excluded_benchmark_ids': sorted(excluded)},
                       'n_complete_residues': len(records)},
        'bonds': bond_stats,
        'angles': angle_stats,
        'sugar_templates': sugar_templates,
        'base_templates': base_templates,
        'phosphate_template': phosphate_template,
        'pair_templates': pair_templates,
        'aform_torsions': aform,
        'aform_torsion_sd': aform_sd,
        'n_helix_residues': int(len(hel)),
        'torsion_histograms_10deg': hist,
        'loop_conformer_library': {'columns': list(names) + ['phase'],
                                   'rows': [[round(float(v), 1) for v in row] for row in lib],
                                   'n_total_observed': int(len(loop_conformers))},
    }
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(out, indent=1))
    print('A-form torsions:', aform, 'sd', aform_sd)
    print('pair templates:', {k: v['n'] for k, v in pair_templates.items()})
    print('sugar templates n:', {k: v['n'] for k, v in sugar_templates.items()})
    print('loop conformers:', len(loop_conformers))
    for k in ("O3'-P", "C3'-O3'-P", "O3'-P-O5'"):
        print(k, bond_stats.get(k) or angle_stats.get(k))


if __name__ == '__main__':
    main()
