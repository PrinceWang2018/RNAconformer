"""Methods paragraph and citation text for an ensemble document (copied into papers, stored in .rcv files)."""
from __future__ import annotations

from . import __author__, __version__


def _range(vals, unit=''):
    vals = [v for v in vals if isinstance(v, (int, float))]
    if not vals:
        return None
    lo, hi = min(vals), max(vals)
    return f'{lo:.1f}{unit}' if abs(hi - lo) < 0.05 else f'{lo:.1f}–{hi:.1f}{unit}'


def methods_text(doc, vienna_version=None):
    """A factual methods paragraph built only from what the document records."""
    seq = doc.get('sequence', '')
    n = len(doc.get('conformations', []))
    name = doc.get('name', 'the RNA')
    prov = doc.get('provenance') or {}
    sw = prov.get('software') or {}
    par = prov.get('parameters') or {}
    stats = prov.get('statistics') or {}
    origin = doc.get('origin')
    parts = []
    short = seq if len(seq) <= 60 else f'{seq[:30]}…{seq[-20:]}'
    parts.append(f'Conformational ensemble of {name} ({len(seq)} nt; 5′-{short}-3′).')
    if origin == 'rnafold3d-engine':
        ss = prov.get('secondary_structures') or []
        mfe = [s for s in ss if s.get('source') == 'vienna_mfe']
        sub = [s for s in ss if s.get('source') == 'vienna_suboptimal']
        user = [s for s in ss if s.get('source') == 'user']
        v2d = f'ViennaRNA {vienna_version}' if vienna_version else 'ViennaRNA'
        if user:
            hyp = 'the secondary structure supplied by the user'
        else:
            hyp = (f'{len(ss)} secondary-structure hypotheses from {v2d} (minimum free energy'
                   + (f' and {len(sub)} suboptimal structure{"s" if len(sub) != 1 else ""}' if sub else '') + ')') if ss \
                else f'secondary structures from {v2d}'
        temp = par.get('temperature')
        parts.append(
            f'Three-dimensional conformations were predicted with the rnafold3d engine '
            f'{sw.get("engine_version", "")}'.rstrip() + f' bundled with RCV {__version__} (RNA Conformation Viewer), a '
            'physics-inspired method that does not use deep learning: '
            f'{hyp}{f" at {temp:g} °C" if isinstance(temp, (int, float)) else ""} were assembled from calibrated A-form '
            'helices by fold-tree simulated annealing with a coarse-grained score (steric clash, loop feasibility, '
            'helix stacking, Debye–Hückel electrostatics and a knowledge-based torsion prior), loops were closed by '
            'cyclic coordinate descent, and every model was checked by an all-atom stereochemical validator.')
        tr = stats.get('trajectories') or par.get('trajectories_run') or par.get('trajectories')
        if tr:
            kept = f'{n} structurally distinct conformation{"s" if n != 1 else ""} were kept'
            cl = par.get('cluster_rmsd_threshold')
            if cl:
                kept += f' after clustering at {cl:.1f} Å C4′ RMSD'
            parts.append(f'{tr} independent trajectories were sampled (random seed {par.get("seed", "–")}); '
                         f'{stats.get("valid", "–")} passed validation and {kept}.')
        if par.get('ionic_conditions'):
            parts.append(f'Ionic conditions: {par["ionic_conditions"].replace(" A ", " Å ")}.')
        explored = [c for c in doc.get('conformations', []) if c.get('explored')]
        if explored:
            parts.append(f'{len(explored)} conformation{"s were" if len(explored) != 1 else " was"} added by '
                         'pull-guided exploration in RCV: temporary distance restraints derived from an interactive '
                         'pull biased the sampling, and each model was then relaxed without restraints.')
    else:
        src = doc.get('source') or {}
        files = [p.replace('\\', '/').rsplit('/', 1)[-1] for p in src.get('paths', [])] or src.get('files') or []
        what = f' from {", ".join(files[:4])}{" and others" if len(files) > 4 else ""}' if files else ''
        parts.append(f'{n} conformation{"s were" if n != 1 else " was"} read{what} and analysed with RCV {__version__} '
                     '(RNA Conformation Viewer).')
    ens = doc.get('ensemble') or {}
    rms = [v for row in ens.get('rmsd_c4', []) for v in row if v > 0]
    sup = ens.get('superposition', "C4' atoms")
    parts.append(f'Conformations were superposed on {sup.replace(chr(39), "′")}' +
                 (f'; pairwise C4′ RMSD ranged from {min(rms):.1f} to {max(rms):.1f} Å.' if rms else '.'))
    parts.append('Base pairs, stacking and pseudoknots were annotated from the 3D coordinates; conformations were '
                 'compared with RCV’s model-independent evaluator (Turner free energy of the annotated 2D structure'
                 + (f' computed with ViennaRNA {vienna_version}' if vienna_version else '') +
                 ', loop/tertiary stacking, electrostatic and steric-clash terms). '
                 'Model scores rank conformations of one sequence and are not experimental free energies; '
                 'morph animations are linear interpolations between superposed conformations, not simulated pathways.')
    return ' '.join(p.strip() for p in parts if p)


def citation_text():
    return (f'{__author__}. RCV — RNA Conformation Viewer, version {__version__} (software). Open-source software '
            '(GPL-3.0) for predicting, exploring and presenting RNA conformational ensembles with the physics-inspired '
            'rnafold3d engine.')
