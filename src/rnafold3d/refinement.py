"""M3.2 constrained rigid-nucleotide chain refinement (experimental).

This optimizer repairs inter-residue placement without changing intra-residue
chemistry. It is NOT torsional RNA backbone refinement or force-field MD.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from .atomistic import Atom, validate_atoms

@dataclass(frozen=True)
class RefinementConfig:
    bond_target: float = 1.61
    bond_weight: float = 12.0
    anchor_weight: float = 0.8
    angle_weight: float = 2.0
    clash_weight: float = 3.0
    clash_distance: float = 2.0
    max_nfev: int = 120
    max_residues: int = 40  # finite-difference prototype guard


def _angle(a, b, c):
    v=a-b; w=c-b
    return np.arccos(np.clip(np.dot(v,w)/(np.linalg.norm(v)*np.linalg.norm(w)+1e-12),-1,1))


def refine_rigid_chain(residues, config=None):
    """Constrained placement of intact residues; return (residues, diagnostics).

    The inter-residue linkage and two flanking bond angles are restrained.
    All residue atoms move together as rigid bodies, preserving internal chirality.
    """
    cfg=config or RefinementConfig()
    n=len(residues)
    if not 1 <= n <= cfg.max_residues:
        raise ValueError(f'prototype supports 1..{cfg.max_residues} residues')
    mandatory={"P","O3'","O5'","C3'"}
    for i,res in enumerate(residues):
        if not mandatory.issubset(res):
            raise ValueError(f'residue {i+1} missing linkage atoms')
    names=[list(r) for r in residues]
    original=[np.array([r[k].xyz for k in ns],float) for r,ns in zip(residues,names)]
    centers=[x.mean(axis=0) for x in original]
    local=[x-c for x,c in zip(original,centers)]
    idx=[{name:j for j,name in enumerate(ns)} for ns in names]
    initial=np.concatenate([np.zeros((n,3)),np.asarray(centers)],axis=1).ravel()
    # Use initial inter-residue angles as a soft fallback only if sensible.
    def xyz(x):
        poses=x.reshape(n,6)
        return [Rotation.from_rotvec(p[:3]).apply(loc)+p[3:] for p,loc in zip(poses,local)]
    def residual(x):
        pts=xyz(x); out=[]
        for i in range(n):
            # Preserve input guide placement, but allow modest movement.
            for key in ('P',"C4'",'N9','N1'):
                if key in idx[i]:
                    j=idx[i][key]
                    out.extend((cfg.anchor_weight*(pts[i][j]-original[i][j])).tolist())
        for i in range(n-1):
            left,right=pts[i],pts[i+1]
            o3=left[idx[i]["O3'"]]; p=right[idx[i+1]['P']]
            c3=left[idx[i]["C3'"]]; o5=right[idx[i+1]["O5'"]]
            out.append(cfg.bond_weight*(np.linalg.norm(o3-p)-cfg.bond_target))
            out.append(cfg.angle_weight*(_angle(c3,o3,p)-np.deg2rad(119.0)))
            out.append(cfg.angle_weight*(_angle(o3,p,o5)-np.deg2rad(104.0)))
        # Nonlocal clashes: exclude adjacent residues, avoid O(N^2) within-residue.
        for i in range(n):
            for j in range(i+2,n):
                for ai,k in enumerate(names[i]):
                    if residues[i][k].element.upper()=='H':continue
                    for aj,l in enumerate(names[j]):
                        if residues[j][l].element.upper()=='H':continue
                        d=np.linalg.norm(pts[i][ai]-pts[j][aj])
                        out.append(cfg.clash_weight*min(d-cfg.clash_distance,0.0))
        return np.asarray(out)
    before=validate_atoms(residues)
    if n==1:
        return residues,{'success':True,'message':'single residue; no inter-residue refinement','before':before,'after':before,'optimizer_evaluations':0,'experimental':True}
    result=least_squares(residual,initial,max_nfev=cfg.max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-7)
    coords=xyz(result.x)
    refined=[{key:Atom(key,residues[i][key].element,tuple(map(float,coords[i][j]))) for j,key in enumerate(names[i])} for i in range(n)]
    after=validate_atoms(refined)
    return refined,{'success':bool(result.success),'message':result.message,'before':before,'after':after,'optimizer_evaluations':int(result.nfev),'objective':float(np.sum(result.fun**2)),'experimental':True,'method':'rigid-residue restrained least squares; no torsion optimization'}
