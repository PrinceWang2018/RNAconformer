"""M2.1 bridge: reconstruct P/C4'/N guide sites for closed internal loops.

All positions are synthetic guides, NOT calibrated atom coordinates. Stem sites
remain fixed. Terminal tails and failed closures remain explicitly unmodeled.
"""
from dataclasses import dataclass
import numpy as np
from .geometry import NucleotideSites
from .closure import close_assembly

@dataclass(frozen=True)
class FullChainGuides:
    sites: tuple
    closures: tuple
    missing_residues: tuple
    complete: bool
    max_junction_link_error: float
    n_loop_sites: int
    geometry_status: str = 'synthetic_three_site_guides_not_atomistic'


def _unit(v, fallback=None):
    v=np.asarray(v,dtype=float); n=float(np.linalg.norm(v))
    if n>1e-9:return v/n
    if fallback is not None:return _unit(fallback)
    return np.array([1.,0.,0.])


def _side_vector(tangent, preferred):
    """Project a nearby stem C4->P direction perpendicular to loop tangent."""
    tangent=_unit(tangent)
    projected=np.asarray(preferred)-tangent*np.dot(preferred,tangent)
    if np.linalg.norm(projected)<1e-8:
        axis=np.eye(3)[int(np.argmin(np.abs(tangent)))]
        projected=np.cross(tangent,axis)
    return _unit(projected)


def build_full_chain_guides(model, assembly, *, step=6.0, tolerance=0.05, seed=0):
    """Create guide sites for all *successfully closed* internal unpaired runs.

    At each loop point, C4 is the FABRIK position, P is an offset in a locally
    transported perpendicular direction, and N is a second orthogonal offset.
    The construction preserves all stem sites exactly and checks every link
    crossing between a loop and a stem. No atomic P-O or stereochemical
    continuity is asserted.
    """
    if step<=0 or tolerance<=0:raise ValueError('step and tolerance must be positive')
    sites={}
    for stem in assembly.stem_sites:
        for s in stem:
            if s.residue in sites:raise ValueError('Duplicate stem residue')
            sites[s.residue]=s
    closures=close_assembly(model,assembly,step=step,tolerance=tolerance,seed=seed)
    loop_count=0
    for closure in closures:
        if not closure.converged or not closure.residues:continue
        ids=closure.residues
        a=ids[0]-1;b=ids[-1]+1
        if a not in sites or b not in sites:continue
        left=sites[a];right=sites[b]
        points=[np.asarray(left.C4),*[np.asarray(x) for x in closure.xyz],np.asarray(right.C4)]
        prev_side=_side_vector(points[1]-points[0],np.asarray(left.P)-np.asarray(left.C4))
        for k,idx in enumerate(ids,1):
            tangent=_unit(points[k+1]-points[k-1])
            side=_side_vector(tangent,prev_side)
            # Parallel transport by reprojecting the previous side vector.
            normal=_unit(np.cross(tangent,side))
            c4=points[k]
            p=c4+2.5*side-0.8*tangent
            base=c4-3.6*normal
            sites[idx]=NucleotideSites(idx,tuple(map(float,p)),tuple(map(float,c4)),tuple(map(float,base)))
            prev_side=side;loop_count+=1
    missing=tuple(i for i in range(len(model.sequence)) if i not in sites)
    owner={r:k for k,stem in enumerate(model.stems) for pair in stem.pairs for r in pair}
    link_errors=[]
    for i in range(len(model.sequence)-1):
        if i not in sites or i+1 not in sites:continue
        if i in owner and i+1 in owner and owner[i]==owner[i+1]:continue
        link_errors.append(abs(np.linalg.norm(np.asarray(sites[i].C4)-sites[i+1].C4)-step))
    return FullChainGuides(tuple(sites[i] for i in sorted(sites)),closures,missing,
                           bool(not missing and all(c.converged for c in closures) and
                           (not link_errors or max(link_errors)<=tolerance)),
                           float(max(link_errors,default=0.)),loop_count)
