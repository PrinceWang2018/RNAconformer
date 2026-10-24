"""M1.3: seeded, topology-aware rigid-body placement of idealized RNA stems.

This is a geometric hypothesis generator, not a thermodynamic predictor.
"""
from dataclasses import dataclass
import numpy as np
from .topology import build_topology
from .geometry import aform_three_site, rigid_transform_sites

@dataclass(frozen=True)
class Assembly:
    stem_sites: tuple
    rotations: tuple
    translations: tuple
    seed: int
    proposal: int
    placement_score: float


def _rotation(axis, angle):
    axis=np.asarray(axis,dtype=float)
    axis=axis/np.linalg.norm(axis)
    x,y,z=axis; c=np.cos(angle); s=np.sin(angle); q=1-c
    return np.array([[c+x*x*q,x*y*q-z*s,x*z*q+y*s],
                     [y*x*q+z*s,c+y*y*q,y*z*q-x*s],
                     [z*x*q-y*s,z*y*q+x*s,c+z*z*q]])


def _connections(model):
    """All sequence-adjacent residues belonging to distinct stems.

    A run of unpaired residues between two stem residues contributes a single
    endpoint-to-endpoint connection with len(run)+1 C4'-C4' links.
    """
    owner={r:k for k,stem in enumerate(model.stems) for pair in stem.pairs for r in pair}
    nodes=sorted(owner)
    result=[]
    for a,b in zip(nodes,nodes[1:]):
        if owner[a]!=owner[b] and all(i not in owner for i in range(a+1,b)):
            result.append((owner[a],a,owner[b],b,b-a))
    return tuple(result)


def _score(sites, model, step=6.0):
    coords={s.residue:np.asarray(s.C4) for stem in sites for s in stem}
    score=0.
    for _,a,_,b,nlinks in _connections(model):
        dist=np.linalg.norm(coords[a]-coords[b])
        score+=max(0.,dist-nlinks*step)**2
    # Coarse clash penalty between nonlocal residues on different stems.
    owner={r:k for k,stem in enumerate(model.stems) for p in stem.pairs for r in p}
    residues=sorted(coords)
    for i,a in enumerate(residues):
        for b in residues[i+1:]:
            if owner[a]==owner[b] or abs(a-b)<3:continue
            d=np.linalg.norm(coords[a]-coords[b]);score+=0.5*max(0.,5.0-d)**2
    return float(score)


def assemble(model, n_candidates=8, seed=0, proposals=100, step=6.0):
    """Generate distinct seeded stem arrangements using a planar topology tree.

    Each child stem starts near its parent's inner pair. Whole-stem rigid
    perturbations are accepted via greedy stochastic search against coarse
    reachability/clash penalties. Root stems are placed independently.
    """
    if n_candidates<1 or proposals<0:raise ValueError('Invalid sampling settings')
    top=build_topology(model,strict=True)
    templates=[aform_three_site(s) for s in model.stems]
    parent={child:p for p,child in top.edges}
    order=[]
    def visit(k):
        order.append(k)
        for c in [j for p,j in top.edges if p==k]:visit(c)
    for root in top.junctions[0].child_stems:visit(root)
    rng=np.random.default_rng(seed);out=[]
    for ci in range(n_candidates):
        rotations=[np.eye(3) for _ in templates]
        translations=[np.zeros(3) for _ in templates]
        placed=[None]*len(templates)
        for k in order:
            if k in parent:
                p=parent[k]
                last=model.stems[p].pairs[-1]
                pcoords={s.residue:np.asarray(s.C4) for s in placed[p]}
                anchor=(pcoords[last[0]]+pcoords[last[1]])/2
                direction=rotations[p]@np.array([0.,0.,1.])
                t=anchor+direction*(7.+rng.uniform(-3.,3.))+rng.normal(0,3,3)
            else:
                t=np.array([k*25.,0.,0.])+rng.normal(0,2,3)
            axis=rng.normal(size=3); r=_rotation(axis,rng.uniform(-np.pi,np.pi))
            rotations[k]=r;translations[k]=t
            placed[k]=rigid_transform_sites(templates[k],r,t)
        best=_score(placed,model,step)
        # Coordinate descent with simulated-annealing acceptance to escape traps.
        for it in range(proposals):
            if not order:break
            k=int(rng.choice(order));old_r=rotations[k];old_t=translations[k];old_sites=placed[k]
            temperature=max(0.2, 4.0*(1-it/max(1,proposals)))
            axis=rng.normal(size=3)
            trial_r=_rotation(axis,rng.normal(0,0.25))@old_r
            trial_t=old_t+rng.normal(0,2.0,3)
            placed[k]=rigid_transform_sites(templates[k],trial_r,trial_t)
            trial=_score(placed,model,step)
            if trial<best or rng.random()<np.exp(np.clip((best-trial)/temperature,-700,0)):
                rotations[k]=trial_r;translations[k]=trial_t;best=trial
            else:placed[k]=old_sites
        out.append(Assembly(tuple(tuple(s) for s in placed),tuple(tuple(map(tuple,r)) for r in rotations),tuple(tuple(t) for t in translations),seed,ci,_score(placed,model,step)))
    return tuple(out)
