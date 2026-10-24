"""M1.4: C4' guide-chain closure using FABRIK inverse kinematics.

All endpoints are fixed. No atomic torsion, chirality or collision guarantees.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Closure:
    residues: tuple
    xyz: tuple
    converged: bool
    endpoint_error: float
    max_link_error: float
    iterations: int
    reason: str


def close_chain(start,end,residues,step=6.0,tolerance=1e-3,max_iterations=300,seed=0):
    """Return guide coordinates for consecutive unpaired residues.

    Uses n+1 equal-length virtual C4'-C4' links. Reports infeasible chains
    instead of silently stretching the backbone. Zero-residue joins checked.
    """
    a=np.asarray(start,dtype=float);b=np.asarray(end,dtype=float)
    if a.shape!=(3,) or b.shape!=(3,) or not np.isfinite([*a,*b]).all():raise ValueError('Invalid endpoints')
    n=len(residues);links=n+1
    if step<=0 or tolerance<=0 or max_iterations<1:raise ValueError('Invalid solver parameters')
    d=float(np.linalg.norm(b-a))
    if n==0:
        error=abs(d-step)
        return Closure((),(),error<=tolerance,0.,error,0,'closed' if error<=tolerance else 'fixed_endpoint_mismatch')
    if d>links*step+tolerance:
        return Closure(tuple(residues),(),False,d-links*step,d-links*step,0,'unreachable')
    rng=np.random.default_rng(seed)
    direction=(b-a)/d if d>1e-10 else np.array([1.,0.,0.])
    perpendicular=np.cross(direction,[0.,0.,1.])
    if np.linalg.norm(perpendicular)<1e-8:perpendicular=np.cross(direction,[0.,1.,0.])
    perpendicular/=np.linalg.norm(perpendicular)
    # Noncollinear initialization helps reach compact endpoints.
    points=np.array([a+(b-a)*i/links+perpendicular*np.sin(np.pi*i/links)*max(1.,step)*0.8 for i in range(links+1)])
    points[0]=a;points[-1]=b
    iterations=0
    for iteration in range(max_iterations):
        iterations=iteration+1
        points[-1]=b
        for i in range(links-1,-1,-1):
            v=points[i]-points[i+1];norm=np.linalg.norm(v)
            if norm<1e-12:v=rng.normal(size=3);norm=np.linalg.norm(v)
            points[i]=points[i+1]+step*v/norm
        points[0]=a
        for i in range(links):
            v=points[i+1]-points[i];norm=np.linalg.norm(v)
            if norm<1e-12:v=rng.normal(size=3);norm=np.linalg.norm(v)
            points[i+1]=points[i]+step*v/norm
        if np.linalg.norm(points[-1]-b)<tolerance:break
    endpoint_error=float(np.linalg.norm(points[-1]-b))
    max_link_error=float(np.max(np.abs(np.linalg.norm(np.diff(points,axis=0),axis=1)-step)))
    converged=endpoint_error<=tolerance and max_link_error<=tolerance
    return Closure(tuple(residues),tuple(tuple(float(v) for v in p) for p in points[1:-1]),converged,endpoint_error,max_link_error,iterations,'closed' if converged else 'not_converged')


def close_assembly(model,assembly,step=6.0,tolerance=1e-2,seed=0):
    """Close unpaired runs flanked by fixed stem residues; return all reports.

    Terminal tails lack two fixed endpoints and remain unmodeled.
    """
    known={s.residue:np.asarray(s.C4) for stem in assembly.stem_sites for s in stem}
    results=[];i=0;n=len(model.sequence)
    while i<n:
        if i in known:
            i+=1;continue
        j=i
        while j<n and j not in known:j+=1
        if i>0 and j<n and i-1 in known:
            results.append(close_chain(known[i-1],known[j],tuple(range(i,j)),step,tolerance,seed=seed+i))
        i=j
    # Check adjacent fixed residues belonging to different stems.
    for i in range(n-1):
        if i in known and i+1 in known:
            owner={r:k for k,stem in enumerate(model.stems) for pair in stem.pairs for r in pair}
            if owner[i]!=owner[i+1]:results.append(close_chain(known[i],known[i+1],(),step,tolerance,seed=seed+i))
    return tuple(results)
