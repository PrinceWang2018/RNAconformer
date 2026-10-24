"""M1.5: conservative validation/filtering of *guide-coordinate* assemblies.

This does not certify atomic geometry, stereochemistry, or thermodynamic stability.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .closure import close_assembly

@dataclass(frozen=True)
class ValidationConfig:
    guide_step: float = 6.0
    closure_tolerance: float = 0.05
    clash_distance: float = 4.0
    rmsd_threshold: float = 3.0
    require_complete: bool = True

@dataclass(frozen=True)
class CandidateReport:
    candidate_index: int
    valid: bool
    reasons: tuple
    missing_residues: tuple
    n_clashes: int
    min_nonlocal_distance: float | None
    n_closures: int
    failed_closures: int
    max_endpoint_error: float
    max_link_error: float
    placement_score: float
    cluster_id: int | None = None
    nearest_selected_rmsd: float | None = None


def _coords(model, assembly, config):
    coords = {}
    for stem in assembly.stem_sites:
        for site in stem:
            if site.residue in coords:
                raise ValueError('Duplicate residue in stem sites')
            coords[site.residue] = np.asarray(site.C4, dtype=float)
    closures = close_assembly(model, assembly, step=config.guide_step,
                              tolerance=config.closure_tolerance,
                              seed=assembly.seed + assembly.proposal)
    for closure in closures:
        if closure.converged:
            for residue, xyz in zip(closure.residues, closure.xyz):
                if residue in coords:
                    raise ValueError('Duplicate residue in closure')
                coords[residue] = np.asarray(xyz, dtype=float)
    return coords, closures


def validate_candidate(model, assembly, config=ValidationConfig()):
    if config.guide_step <= 0 or config.closure_tolerance <= 0 or config.clash_distance <= 0 or config.rmsd_threshold <= 0:
        raise ValueError('Validation thresholds must be positive')
    coords, closures = _coords(model, assembly, config)
    missing = tuple(i for i in range(len(model.sequence)) if i not in coords)
    reasons = []
    if config.require_complete and missing:
        reasons.append('missing_residues')
    failures = [c for c in closures if not c.converged]
    if failures:
        reasons.append('loop_closure_failed')
    if any(p.shape != (3,) or not np.isfinite(p).all() for p in coords.values()):
        reasons.append('nonfinite_coordinates')
    # Adjacent modeled residues are checked as a coarse guide-chain diagnostic.
    link_errors = [abs(np.linalg.norm(coords[i+1]-coords[i])-config.guide_step)
                   for i in range(len(model.sequence)-1) if i in coords and i+1 in coords]
    # Idealized stem geometry is not parameterized to a fixed 6A guide link.
    # Only non-stem adjacent links are checked; stem-internal geometry is provisional.
    owner = {r:k for k,s in enumerate(model.stems) for pair in s.pairs for r in pair}
    external_link_errors = [abs(np.linalg.norm(coords[i+1]-coords[i])-config.guide_step)
                            for i in range(len(model.sequence)-1) if i in coords and i+1 in coords
                            and not (i in owner and i+1 in owner and owner[i]==owner[i+1])]
    if external_link_errors and max(external_link_errors) > config.closure_tolerance:
        reasons.append('guide_link_mismatch')
    # Nonlocal pairs only; adjacent/nearby residues are not steric clash candidates.
    # Same-stem sites excluded because analytic helix is an uncalibrated scaffold.
    clashes = 0
    minimum = float('inf')
    for i in sorted(coords):
        for j in sorted(coords):
            if j <= i+2 or (i in owner and j in owner and owner[i]==owner[j]):
                continue
            d = float(np.linalg.norm(coords[i]-coords[j]))
            minimum = min(minimum, d)
            if d < config.clash_distance:
                clashes += 1
    if clashes:
        reasons.append('nonlocal_clashes')
    report = CandidateReport(assembly.proposal, not reasons, tuple(reasons), missing, clashes,
                             None if minimum==float('inf') else minimum, len(closures), len(failures),
                             max((c.endpoint_error for c in closures), default=0.),
                             max((c.max_link_error for c in closures), default=0.),
                             float(assembly.placement_score))
    return report, coords, closures


def aligned_rmsd(a, b):
    """Kabsch-aligned RMSD of matching residue-coordinate dictionaries."""
    keys = sorted(set(a) & set(b))
    if not keys or set(a) != set(b):
        raise ValueError('RMSD requires matching, nonempty residue sets')
    x = np.stack([a[i] for i in keys]); y = np.stack([b[i] for i in keys])
    x = x-x.mean(axis=0); y = y-y.mean(axis=0)
    u, _, vt = np.linalg.svd(x.T @ y)
    sign = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1., 1., np.sign(sign)]) @ vt
    return float(np.sqrt(np.mean(np.sum((x@rotation-y)**2, axis=1))))


def filter_candidates(model, assemblies, config=ValidationConfig(), max_candidates=20):
    """Return (accepted assemblies, all reports, accepted coordinate dictionaries).

    Greedy quality-first RMSD clustering; invalid candidates never enter clusters.
    """
    if max_candidates < 1:
        raise ValueError('max_candidates must be >= 1')
    inspected = [(a, *validate_candidate(model,a,config)) for a in assemblies]
    inspected.sort(key=lambda item:(not item[1].valid, item[1].placement_score, item[0].proposal))
    selected=[]; selected_coords=[]; reports={}
    for assembly, report, coords, _ in inspected:
        if not report.valid:
            reports[assembly.proposal]=report
            continue
        distances = [aligned_rmsd(coords, prev) for prev in selected_coords]
        nearest = min(distances, default=None)
        if nearest is not None and nearest < config.rmsd_threshold:
            reports[assembly.proposal] = CandidateReport(**{**asdict(report), 'valid':False,'reasons':('duplicate_cluster',),
                                                           'cluster_id':int(np.argmin(distances)),
                                                           'nearest_selected_rmsd':nearest})
        elif len(selected) < max_candidates:
            cid=len(selected)
            selected.append(assembly);selected_coords.append(coords)
            reports[assembly.proposal] = CandidateReport(**{**asdict(report),'cluster_id':cid,
                                                           'nearest_selected_rmsd':nearest})
        else:
            reports[assembly.proposal] = CandidateReport(**{**asdict(report),'valid':False,'reasons':('selection_limit',),
                                                           'nearest_selected_rmsd':nearest})
    return tuple(selected), tuple(reports[a.proposal] for a in assemblies), tuple(selected_coords)
