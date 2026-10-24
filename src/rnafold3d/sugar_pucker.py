"""Experimental restrained five-membered ribose ring pucker optimizer.

Fits five ring atoms to a target pseudorotation using Altona-Sundaralingam
five torsions while restraining bond lengths, angles, and local substituents.
Not a force field or chemically validated RNA reconstruction.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares
from .atomistic import Atom

RING = ("O4'", "C1'", "C2'", "C3'", "C4'")
# ν0 C4'-O4'-C1'-C2'; ν1 O4'-C1'-C2'-C3'; ...
TORSION_QUADS = ((4,0,1,2),(0,1,2,3),(1,2,3,4),(2,3,4,0),(3,4,0,1))

@dataclass(frozen=True)
class PuckerConfig:
    phase_deg: float = 18.0  # nominal C3'-endo region, convention approximate
    amplitude_deg: float = 35.0
    torsion_weight: float = 1.0
    bond_weight: float = 20.0
    angle_weight: float = 8.0
    anchor_weight: float = 1.5
    chirality_weight: float = 20.0
    max_displacement: float = 0.65
    max_nfev: int = 180


def dihedral(a,b,c,d):
    a,b,c,d=(np.asarray(v,float) for v in (a,b,c,d))
    v1=a-b; axis=c-b; norm=np.linalg.norm(axis)
    if norm<1e-10:raise ValueError('degenerate torsion axis')
    axis=axis/norm
    v2=d-c
    x=v1-np.dot(v1,axis)*axis; y=v2-np.dot(v2,axis)*axis
    if np.linalg.norm(x)<1e-10 or np.linalg.norm(y)<1e-10:raise ValueError('degenerate dihedral')
    return float(np.arctan2(np.dot(np.cross(axis,x),y),np.dot(x,y)))


def _angle(a,b,c):
    u=a-b;v=c-b
    return float(np.arccos(np.clip(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)+1e-12),-1,1)))


def _wrap(x):return np.arctan2(np.sin(x),np.cos(x))


def ring_torsions(residue):
    x=np.array([residue[k].xyz for k in RING],float)
    return np.array([dihedral(*(x[j] for j in q)) for q in TORSION_QUADS])


def pseudorotation(residue):
    """Return approximate phase/amplitude fitted to ν_j=A cos(P+4πj/5)."""
    nu=ring_torsions(residue)
    basis=np.array([[np.cos(4*np.pi*j/5),-np.sin(4*np.pi*j/5)] for j in range(5)])
    a,b=np.linalg.lstsq(basis,nu,rcond=None)[0]
    return {'phase_deg':float(np.rad2deg(np.arctan2(b,a))%360),
            'amplitude_deg':float(np.rad2deg(np.hypot(a,b))),
            'fit_rms_deg':float(np.rad2deg(np.sqrt(np.mean((basis@np.array([a,b])-nu)**2))))}


def _signed_volume(a,b,c,d):return float(np.dot(np.cross(a-d,b-d),c-d))


def refine_sugar_pucker(residue, config=None):
    """Locally adjust ring coordinates while preserving covalent geometry approximately.

    Adjacent backbone and base atoms remain fixed; ring substituents are not
    rebuilt, so results require independent stereochemical checks.
    """
    cfg=config or PuckerConfig()
    required=set(RING)|{"O2'","O3'","C5'"}
    if not required.issubset(residue):raise ValueError('incomplete ribose residue')
    if not 0<cfg.max_displacement<=2 or cfg.max_nfev<1:raise ValueError('invalid pucker settings')
    start=np.array([residue[k].xyz for k in RING],float)
    if not np.isfinite(start).all():raise ValueError('nonfinite coordinates')
    bonds=[(i,(i+1)%5) for i in range(5)]
    angles=[((i-1)%5,i,(i+1)%5) for i in range(5)]
    bond_targets=[np.linalg.norm(start[i]-start[j]) for i,j in bonds]
    angle_targets=[_angle(start[i],start[j],start[k]) for i,j,k in angles]
    chirality0=_signed_volume(start[1],start[2],start[3],start[4])
    # Phase sign depends on dihedral convention; optimize relative to measured
    # initial phase to avoid an unjustified stereochemical inversion.
    target=np.deg2rad(cfg.amplitude_deg)*np.cos(np.deg2rad(cfg.phase_deg)+4*np.pi*np.arange(5)/5)
    def residual(v):
        x=v.reshape(5,3)
        nu=np.array([dihedral(*(x[j] for j in q)) for q in TORSION_QUADS])
        out=list(cfg.torsion_weight*_wrap(nu-target))
        out.extend(cfg.bond_weight*(np.linalg.norm(x[i]-x[j])-t) for (i,j),t in zip(bonds,bond_targets))
        out.extend(cfg.angle_weight*(_angle(x[i],x[j],x[k])-t) for (i,j,k),t in zip(angles,angle_targets))
        out.extend((cfg.anchor_weight*(x-start)).ravel())
        volume=_signed_volume(x[1],x[2],x[3],x[4])
        out.append(cfg.chirality_weight*max(0.0,0.4*abs(chirality0)-np.sign(chirality0)*volume))
        return np.asarray(out,float)
    low=(start-cfg.max_displacement).ravel();high=(start+cfg.max_displacement).ravel()
    res=least_squares(residual,start.ravel(),bounds=(low,high),max_nfev=cfg.max_nfev)
    xyz=res.x.reshape(5,3)
    out=dict(residue)
    for i,k in enumerate(RING):
        a=residue[k];out[k]=Atom(a.name,a.element,tuple(map(float,xyz[i])))
    return out,{'success':bool(res.success),'objective_before':float(np.sum(residual(start.ravel())**2)),
                'objective_after':float(np.sum(res.fun**2)),
                'max_ring_displacement':float(np.max(np.linalg.norm(xyz-start,axis=1))),
                'initial_pucker':pseudorotation(residue),'final_pucker':pseudorotation(out),
                'experimental':True,'warning':'Ring-only restrained optimization; substituent bonds and global chain stereochemistry not guaranteed'}


def refine_chain_puckers(residues,config=None):
    out=[];reports=[]
    for r in residues:
        updated,report=refine_sugar_pucker(r,config)
        out.append(updated);reports.append(report)
    return out,reports
