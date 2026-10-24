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


def predict_structures(sequence, max_structures=3, temperature=37.0, no_lonely_pairs=True,
                       max_delta=4.0, min_distance=None):
    """Return a list of hypotheses [{id, dbn, pairs (0-based), energy, source, ...}] (MFE first)."""
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
    mfe_s, mfe_e = fc.mfe()
    fc.exp_params_rescale(mfe_e)
    _, ens_g = fc.pf()
    n_mfe_pairs = mfe_s.count('(')
    min_distance = min_distance or max(3, int(round(0.15 * max(n_mfe_pairs, 1))))
    out = [{'id': 'ss1', 'dbn': mfe_s, 'energy_kcal_mol': round(mfe_e, 2), 'source': 'vienna_mfe',
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
                out.append({'id': f'ss{len(out) + 1}', 'dbn': s, 'energy_kcal_mol': round(e, 2),
                            'source': 'vienna_suboptimal', 'probability_in_ensemble': round(fc.pr_structure(s), 5),
                            'bp_distance_to_mfe': RNA.bp_distance(s, mfe_s)})
    for h in out:
        h['pairs'] = [(i - 1, j - 1) for i, j in dotbracket_to_pairs(h['dbn'])]
    meta = {'program': 'ViennaRNA', 'version': RNA.__version__, 'temperature_C': temperature,
            'no_lonely_pairs': bool(no_lonely_pairs), 'dangles': md.dangles,
            'ensemble_free_energy_kcal_mol': round(ens_g, 2), 'min_bp_distance_between_hypotheses': min_distance,
            'energy_parameters': 'ViennaRNA default (Turner 2004)'}
    return out, meta
