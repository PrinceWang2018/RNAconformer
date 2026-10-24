"""M2.1 geometric consistency diagnostics and local guide-frame regularization.

These are MODEL targets for three-site guides, not experimental atomistic force-field
parameters. This module does not claim to reconstruct sugar puckers, torsions or bonds.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .geometry import NucleotideSites

@dataclass(frozen=True)
class GuideGeometryConfig:
    p_radius: float = 3.0
    n_radius: float = 3.6
    min_frame_sine: float = 0.12
    c4_neighbor_min: float = 3.0
    c4_neighbor_max: float = 8.5
    interface_target: float = 6.0
    interface_tolerance: float = 1.0

@dataclass(frozen=True)
class GeometryReport:
    n_residues: int
    n_interfaces: int
    n_degenerate_frames: int
    n_bad_c4_links: int
    n_bad_interfaces: int
    max_interface_deviation: float
    max_p_radius_deviation: float
    max_n_radius_deviation: float
    n_base_flip_events: int
    passed_guide_checks: bool
    status: str = 'three_site_guide_consistency_only_not_atomistic_validation'
    def as_dict(self): return asdict(self)

def _unit(x):
    x=np.asarray(x,dtype=float)
    n=np.linalg.norm(x)
    if n<1e-9: return None
    return x/n

def regularize_loop_frames(sites, stem_residues, config=GuideGeometryConfig()):
    """Project loop P/N guide offsets onto a stable local frame.

    Stem coordinates are left untouched; C4 coordinates are NEVER modified.
    Uses parallel-transport-like projection of the previous transverse axis.
    """
    if config.p_radius<=0 or config.n_radius<=0: raise ValueError('Guide radii must be positive')
    ordered=sorted(sites,key=lambda s:s.residue)
    if len({s.residue for s in ordered})!=len(ordered): raise ValueError('Duplicate residue')
    lookup={s.residue:s for s in ordered}
    fixed=set(stem_residues)
    out=[];prev_side=None
    for s in ordered:
        i=s.residue
        if i in fixed:
            out.append(s);prev_side=None;continue
        left=lookup.get(i-1);right=lookup.get(i+1)
        if left is None or right is None:
            out.append(s);prev_side=None;continue
        tangent=_unit(np.asarray(right.C4)-np.asarray(left.C4))
        if tangent is None:
            out.append(s);continue
        pref=prev_side if prev_side is not None else np.asarray(s.P)-np.asarray(s.C4)
        side=pref-np.dot(pref,tangent)*tangent
        side=_unit(side)
        if side is None:
            axis=np.eye(3)[np.argmin(np.abs(tangent))]
            side=_unit(np.cross(tangent,axis))
        normal=_unit(np.cross(tangent,side))
        c4=np.asarray(s.C4)
        p=c4+config.p_radius*side
        n=c4-config.n_radius*normal
        out.append(NucleotideSites(i,tuple(map(float,p)),s.C4,tuple(map(float,n))))
        prev_side=side
    return tuple(out)

def validate_guide_geometry(sites, stem_residues=(), config=GuideGeometryConfig()):
    ordered=sorted(sites,key=lambda s:s.residue)
    if len({s.residue for s in ordered})!=len(ordered): raise ValueError('Duplicate residue')
    lookup={s.residue:s for s in ordered};fixed=set(stem_residues)
    degenerate=bad_links=bad_interfaces=flips=0
    max_iface=max_p=max_n=0.0
    n_iface=0;last_normal=None
    for s in ordered:
        coords=np.asarray([s.P,s.C4,s.N],dtype=float)
        if not np.all(np.isfinite(coords)):raise ValueError('Non-finite coordinate')
        p=coords[0]-coords[1];n=coords[2]-coords[1]
        if s.residue not in fixed:
            max_p=max(max_p,abs(np.linalg.norm(p)-config.p_radius))
            max_n=max(max_n,abs(np.linalg.norm(n)-config.n_radius))
        pu=_unit(p);nu=_unit(n)
        if pu is None or nu is None or np.linalg.norm(np.cross(pu,nu))<config.min_frame_sine:
            degenerate+=1
        if s.residue not in fixed and nu is not None:
            if last_normal is not None and np.dot(last_normal,nu)<-0.8:flips+=1
            last_normal=nu
        else:last_normal=None
        other=lookup.get(s.residue+1)
        if other is None:continue
        d=np.linalg.norm(np.asarray(s.C4)-np.asarray(other.C4))
        if d<config.c4_neighbor_min or d>config.c4_neighbor_max:bad_links+=1
        if (s.residue in fixed)!=(other.residue in fixed):
            n_iface+=1;err=abs(d-config.interface_target)
            max_iface=max(max_iface,err)
            if err>config.interface_tolerance:bad_interfaces+=1
    passed=degenerate==0 and bad_links==0 and bad_interfaces==0
    return GeometryReport(len(ordered),n_iface,degenerate,bad_links,bad_interfaces,
                          float(max_iface),float(max_p),float(max_n),flips,passed)
