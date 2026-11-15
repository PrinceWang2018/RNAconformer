"""Secondary-structure hypotheses from sequence alone (thermodynamic, no deep learning).

Uses the ViennaRNA package (Lorenz et al. 2011; Turner 2004 nearest-neighbour energies):
the minimum-free-energy structure plus structurally distinct low-energy suboptimal
structures from deterministic Wuchty enumeration. Energies are kcal/mol at the given
temperature and are genuine nearest-neighbour model free energies (unlike the 3D score).
"""
from __future__ import annotations

from .annotate import dotbracket_to_pairs


def vienna_version():
    try:
        import RNA
        return RNA.__version__
    except ImportError:
        return None


SHAPE_SLOPE, SHAPE_INTERCEPT = 1.8, -0.6           # kcal/mol, Deigan et al. 2009


def read_reactivities(path, n):
    """Chemical-probing reactivities (SHAPE, DMS ...) for n nucleotides from a text file: one value per line, or
    "position value" lines (1-based; other columns ignored). Missing: blank, NA, nan or a negative value."""
    from pathlib import Path
    return parse_reactivities(Path(path).read_text(), n)


def parse_reactivities(text, n):
    """read_reactivities for text."""
    rows = [l.replace(',', ' ').split() for l in text.splitlines()
            if l.strip() and not l.lstrip().startswith(('#', '>'))]

    def num(tok):
        try:
            v = float(tok)
        except ValueError:
            return None
        return None if (v != v or v < 0) else v
    rows = [r for r in rows if r and (r[0].lstrip('-').replace('.', '', 1).isdigit() or r[0].lower() in ('na', 'nan'))]
    indexed = bool(rows) and all(len(r) >= 2 and r[0].isdigit() for r in rows)
    vals = [None] * n
    if indexed:
        for r in rows:
            k = int(r[0])
            if not 1 <= k <= n:
                raise ValueError(f'probing data: position {k} is outside the sequence (1-{n})')
            vals[k - 1] = num(r[1])
    else:
        if len(rows) != n:
            raise ValueError(f'probing data: {len(rows)} values for {n} nucleotides')
        vals = [num(r[0]) for r in rows]
    if sum(v is not None for v in vals) < max(5, n // 4):
        raise ValueError('probing data: fewer than a quarter of the nucleotides have a value')
    return vals


def read_alignment(path, sequence):
    """Aligned homologous sequences (FASTA with gaps, Stockholm or Clustal) with the query in the first row.
    The file must contain the query sequence; it is moved to the first row. Returns rows of equal length."""
    from pathlib import Path
    return parse_alignment(Path(path).read_text(), sequence)


def parse_alignment(text, sequence):
    """read_alignment for text."""
    names, rows = [], {}
    cur = None
    for line in text.splitlines():
        line = line.rstrip()
        if not line or line.startswith(('#', '//', 'CLUSTAL', 'MUSCLE')) or set(line) <= set(' *:.'):
            continue
        if line.startswith('>'):
            cur = line[1:].strip() or f'seq{len(names) + 1}'
            if cur not in rows:
                names.append(cur)
                rows[cur] = ''
            continue
        if cur is not None:                       # FASTA
            rows[cur] += ''.join(line.split())
        else:                                     # Stockholm / Clustal: name, aligned block
            parts = line.split()
            if len(parts) < 2:
                continue
            if parts[0] not in rows:
                names.append(parts[0])
                rows[parts[0]] = ''
            rows[parts[0]] += parts[1]
    aln = [rows[k].upper().replace('T', 'U').replace('.', '-').replace('~', '-') for k in names]
    aln = [a for a in aln if a]
    if len(aln) < 2 or len({len(a) for a in aln}) != 1:
        raise ValueError('alignment: needs at least two aligned sequences of the same length')
    q = [k for k, a in enumerate(aln) if a.replace('-', '') == sequence]
    if not q:
        raise ValueError('alignment: none of the aligned sequences is the query sequence')
    aln.insert(0, aln.pop(q[0]))
    return [''.join(ch if ch in 'ACGU-' else 'N' for ch in a) for a in aln]


def alignment_consensus(sequence, alignment, temperature=37.0):
    """Consensus secondary structure of an alignment (ViennaRNA RNAalifold: thermodynamics plus covariation),
    mapped onto the query (first row); only pairs that are canonical in the query are kept."""
    import RNA
    md = RNA.md()
    md.temperature = float(temperature)
    cons, _e = RNA.fold_compound(list(alignment), md).mfe()
    col_to_q, k = {}, 0
    for c, ch in enumerate(alignment[0]):
        if ch != '-':
            col_to_q[c] = k
            k += 1
    pairs = []
    for i, j in dotbracket_to_pairs(cons):
        a, b = col_to_q.get(i - 1), col_to_q.get(j - 1)
        if a is not None and b is not None and b - a > 3 and sequence[a] + sequence[b] in ('AU', 'UA', 'GC', 'CG', 'GU', 'UG'):
            pairs.append((a, b))
    s = ['.'] * len(sequence)
    for a, b in pairs:
        s[a], s[b] = '(', ')'
    return ''.join(s), pairs


def predict_structures(sequence, max_structures=3, temperature=37.0, no_lonely_pairs=True,
                       max_delta=4.0, min_distance=None, reactivities=None, alignment=None):
    """Return a list of hypotheses [{id, dbn, pairs (0-based), energy, source, ...}] (MFE first).

    reactivities: SHAPE-like values per nucleotide (None = missing); they enter the folding as pseudo-energies
    (Deigan et al. 2009), so the hypotheses agree with the experiment. alignment: aligned homologs (query
    first); the covariation-supported consensus structure becomes the first hypothesis."""
    try:
        import RNA
    except ImportError as exc:
        raise RuntimeError('ViennaRNA (python module "RNA") is required to predict secondary structure '
                           'from sequence alone; install it (e.g. conda install -c bioconda viennarna) '
                           'or provide --dbn') from exc
    md = RNA.md()
    md.noLP = 1 if no_lonely_pairs else 0
    md.temperature = float(temperature)
    fc = RNA.fold_compound(sequence, md)
    probing = reactivities is not None and len(reactivities) == len(sequence)
    if probing:
        fc.sc_add_SHAPE_deigan([-999.0] + [(-999.0 if v is None else max(0.0, float(v))) for v in reactivities],
                               SHAPE_SLOPE, SHAPE_INTERCEPT)
    plain = RNA.fold_compound(sequence, md, RNA.OPTION_EVAL_ONLY) if probing else None

    def energies(s, e):
        """Turner free energy; with probing data also the energy that includes the pseudo-energies (ranking)."""
        if not probing:
            return {'energy_kcal_mol': round(e, 2)}
        return {'energy_kcal_mol': round(plain.eval_structure(s), 2), 'rank_energy_kcal_mol': round(e, 2)}
    mfe_s, mfe_e = fc.mfe()
    fc.exp_params_rescale(mfe_e)
    _, ens_g = fc.pf()
    n_mfe_pairs = mfe_s.count('(')
    min_distance = min_distance or max(3, int(round(0.15 * max(n_mfe_pairs, 1))))
    out = [{'id': 'ss1', 'dbn': mfe_s, **energies(mfe_s, mfe_e), 'source': 'vienna_mfe_probing' if probing else 'vienna_mfe',
            'probability_in_ensemble': round(fc.pr_structure(mfe_s), 5), 'bp_distance_to_mfe': 0}]
    if max_structures > 1:
        delta = 1.0
        subs = []
        while delta <= max_delta:
            subs = sorted(((s.energy, s.structure) for s in fc.subopt(int(delta * 100))), key=lambda x: (x[0], x[1]))
            if len(subs) >= 60 or len(subs) > 20000:
                break
            delta += 1.0
        for e, s in subs:
            if len(out) >= max_structures:
                break
            if all(RNA.bp_distance(s, o['dbn']) >= min_distance for o in out):
                out.append({'id': f'ss{len(out) + 1}', 'dbn': s, **energies(s, e),
                            'source': 'vienna_suboptimal_probing' if probing else 'vienna_suboptimal',
                            'probability_in_ensemble': round(fc.pr_structure(s), 5),
                            'bp_distance_to_mfe': RNA.bp_distance(s, mfe_s)})
    meta = {'program': 'ViennaRNA', 'version': RNA.__version__, 'temperature_C': temperature,
            'no_lonely_pairs': bool(no_lonely_pairs), 'dangles': md.dangles,
            'ensemble_free_energy_kcal_mol': round(ens_g, 2), 'min_bp_distance_between_hypotheses': min_distance,
            'energy_parameters': 'ViennaRNA default (Turner 2004)'}
    if probing:
        meta['probing'] = {'nucleotides_with_data': sum(v is not None for v in reactivities),
                           'method': f'pseudo-energies, slope {SHAPE_SLOPE} intercept {SHAPE_INTERCEPT} kcal/mol (Deigan et al. 2009)'}
    if alignment:
        cons, _pairs = alignment_consensus(sequence, alignment, temperature)
        meta['alignment'] = {'sequences': len(alignment), 'columns': len(alignment[0]), 'consensus_pairs': cons.count('(')}
        if cons.count('(') >= 3 and all(o['dbn'] != cons for o in out):
            e = fc.eval_structure(cons)
            out.insert(0, {'id': '', 'dbn': cons, **energies(cons, e), 'source': 'alignment_consensus',
                           'probability_in_ensemble': round(fc.pr_structure(cons), 5),
                           'bp_distance_to_mfe': RNA.bp_distance(cons, mfe_s)})
            out = out[:max(max_structures, 2)]
        elif cons.count('(') >= 3:
            k = next(k for k, o in enumerate(out) if o['dbn'] == cons)
            out.insert(0, out.pop(k))
            out[0]['source'] += '+alignment_consensus'
    for k, h in enumerate(out, 1):
        h['id'] = f'ss{k}'
        h['pairs'] = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(h['dbn'])]
    return out, meta
