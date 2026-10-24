"""Analytic A-form-inspired three-site nucleotide scaffold.

Coordinates are idealized P/C4'/N base-anchor guide sites, not validated atom positions.
They must not be exported as an all-atom PDB. The P and C4' sites are constructed
by offsets from a helical local frame; geometry is for rigid-body initialization.
"""
from dataclasses import dataclass
import numpy as np
from .core import transform

@dataclass(frozen=True)
class NucleotideSites:
    residue: int
    P: tuple
    C4: tuple
    N: tuple

def aform_three_site(stem, rise=2.8, twist_deg=32.7, radius=9.0, phase_deg=0.0):
    if not stem.pairs:raise ValueError('Empty stem')
    output=[]
    for k,(i,j) in enumerate(stem.pairs):
        for idx,offset,orientation in ((i,0.,1),(j,np.pi,-1)):
            a=np.deg2rad(phase_deg+k*twist_deg)+offset
            radial=np.array([np.cos(a),np.sin(a),0.])
            tang=np.array([-np.sin(a),np.cos(a),0.])
            axis=np.array([0.,0.,1.])
            # Synthetic geometry: offsets are intentionally not force-field parameters.
            c4=radius*radial + k*rise*axis
            p=c4+1.9*radial-2.1*orientation*tang-1.0*orientation*axis
            base=c4-3.6*radial+0.6*orientation*axis
            output.append(NucleotideSites(idx,tuple(p),tuple(c4),tuple(base)))
    return tuple(sorted(output,key=lambda n:n.residue))

def rigid_transform_sites(sites,rotation,translation):
    out=[]
    for site in sites:
        xyz=transform(np.array([site.P,site.C4,site.N]),rotation,translation)
        out.append(NucleotideSites(site.residue,*map(tuple,xyz)))
    return tuple(out)

def stem_endpoint_frames(sites,stem):
    lookup={s.residue:s for s in sites}
    frames=[]
    for pair in (stem.pairs[0],stem.pairs[-1]):
        p1=np.asarray(lookup[pair[0]].C4);p2=np.asarray(lookup[pair[1]].C4)
        center=(p1+p2)/2
        radial=(p1-p2)/np.linalg.norm(p1-p2)
        z=np.array([0.,0.,1.])
        y=np.cross(z,radial);y/=np.linalg.norm(y)
        frames.append({'pair':pair,'center':tuple(center),'axes':tuple(map(tuple,np.column_stack((radial,y,z)).T))})
    return tuple(frames)
