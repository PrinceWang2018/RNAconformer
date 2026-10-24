"""Experimental RNA backbone torsion-space refinement.

Rotates chemically defined *acyclic* bonds while keeping bond lengths and
intra-fragment geometry exactly unchanged. This is not a complete RNA force
field and does not optimize sugar pucker or ring torsions.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares
from .atomistic import Atom, validate_atoms

# Edges for canonical ribose/phosphate and glycosidic attachments.
BACKBONE_EDGES = [("P","OP1"),("P","OP2"),("P","O5'"),("O5'","C5'"),
 ("C5'","C4'"),("C4'","O4'"),("O4'","C1'"),("C1'","C2'"),
 ("C2'","C3'"),("C3'","C4'"),("C2'","O2'"),("C3'","O3'")]
# Base ring bonds; chi rotates entire base around C1'-N9/N1.
BASE_EDGES = {
 'A': [('N9','C8'),('C8','N7'),('N7','C5'),('C5','C4'),('C4','N9'),('C4','N3'),('N3','C2'),('C2','N1'),('N1','C6'),('C6','C5'),('C6','N6')],
 'G': [('N9','C8'),('C8','N7'),('N7','C5'),('C5','C4'),('C4','N9'),('C4','N3'),('N3','C2'),('C2','N1'),('N1','C6'),('C6','C5'),('C6','O6'),('C2','N2')],
 'C': [('N1','C2'),('C2','N3'),('N3','C4'),('C4','C5'),('C5','C6'),('C6','N1'),('C2','O2'),('C4','N4')],
 'U': [('N1','C2'),('C2','N3'),('N3','C4'),('C4','C5'),('C5','C6'),('C6','N1'),('C2','O2'),('C4','O4')]
}

@dataclass(frozen=True)
class TorsionConfig:
    max_residues: int = 20
    max_nfev: int = 35
    max_rotation_deg: float = 35.0
    bond_target: float = 1.61
    bond_weight: float = 8.0
    angle_weight: float = 2.0
    anchor_weight: float = 0.35
    clash_weight: float = 1.5
    clash_distance: float = 2.0
    include_chi: bool = False


def _graph(residues, sequence):
    adjacency = {(i,k):set() for i,r in enumerate(residues) for k in r}
    def edge(a,b):
        if a in adjacency and b in adjacency:
            adjacency[a].add(b);adjacency[b].add(a)
    for i,r in enumerate(residues):
        for a,b in BACKBONE_EDGES + BASE_EDGES[sequence[i]]:
            edge((i,a),(i,b))
        edge((i,"C1'"),(i,'N9' if sequence[i] in 'AG' else 'N1'))
        if i+1<len(residues):edge((i,"O3'"),(i+1,'P'))
    return adjacency


def _rotatable_edges(residues,sequence,include_chi=False):
    graph=_graph(residues,sequence)
    selected=[]
    for i,r in enumerate(residues):
        # beta: P-O5'; gamma: O5'-C5'; epsilon: C3'-O3';
        # zeta: O3'(i)-P(i+1). Exclude bonds inside the ribose ring.
        proposals=[('beta',(i,'P'),(i,"O5'")),('gamma',(i,"O5'"),(i,"C5'")),
                   ('epsilon',(i,"C3'"),(i,"O3'"))]
        if i+1<len(residues):proposals.append(('zeta',(i,"O3'"),(i+1,'P')))
        if include_chi:proposals.append(('chi',(i,"C1'"),(i,'N9' if sequence[i] in 'AG' else 'N1')))
        for label,a,b in proposals:
            if b not in graph.get(a,()):continue
            # Remove the axis edge and find the b-side component. Cyclic bonds
            # are excluded automatically (alternative graph path exists).
            seen={b};stack=[b]
            while stack:
                u=stack.pop()
                for v in graph[u]:
                    if (u==b and v==a) or (u==a and v==b):continue
                    if v not in seen:seen.add(v);stack.append(v)
            if a in seen:continue
            selected.append((label,a,b,tuple(sorted(seen))))
    return selected


def _angle(a,b,c):
    u=a-b;v=c-b
    return np.arccos(np.clip(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)+1e-12),-1,1))


def refine_torsions(residues, sequence, config=None):
    """Return (refined residues, diagnostics), preserving every rotatable bond length.

    Input residues must be chemically complete canonical RNA templates. The
    method is a local torsion-space repair; it cannot bridge large breaks.
    """
    cfg=config or TorsionConfig()
    n=len(residues)
    if len(sequence)!=n or not 1<=n<=cfg.max_residues:
        raise ValueError('sequence length mismatch or prototype residue limit exceeded')
    if any(b not in 'ACGU' for b in sequence):raise ValueError('canonical RNA only')
    required={"P","O5'","C5'","C4'","C3'","O3'","C1'"}
    for i,r in enumerate(residues):
        if not required.issubset(r):raise ValueError(f'incomplete residue {i+1}')
        if ('N9' if sequence[i] in 'AG' else 'N1') not in r:raise ValueError('missing glycosidic atom')
    keys=[(i,k) for i,r in enumerate(residues) for k in r]
    index={k:j for j,k in enumerate(keys)}
    start=np.asarray([residues[i][k].xyz for i,k in keys],float)
    if not np.isfinite(start).all():raise ValueError('nonfinite coordinates')
    moves=_rotatable_edges(residues,sequence,cfg.include_chi)
    # Optimize only backbone rotations; chi is optional and off by default.
    def apply(angles):
        xyz=start.copy()
        for angle,(_,a,b,affected) in zip(angles,moves):
            origin=xyz[index[a]];axis=xyz[index[b]]-origin
            length=np.linalg.norm(axis)
            if length<1e-9:raise ValueError('degenerate torsion axis')
            unit=axis/length
            ids=[index[k] for k in affected]
            v=xyz[ids]-origin
            xyz[ids]=origin+v*np.cos(angle)+np.cross(unit,v)*np.sin(angle)+np.outer(v@unit,unit)*(1-np.cos(angle))
        return xyz
    links=[(index[(i,"O3'")],index[(i+1,'P')],index[(i,"C3'")],index[(i+1,"O5'")]) for i in range(n-1)]
    anchors=[index[(i,k)] for i,r in enumerate(residues) for k in ("P","C4'",'N9','N1') if k in r]
    # Clash checks only between non-neighboring residues.
    heavy=[j for j,(i,k) in enumerate(keys) if residues[i][k].element.upper()!='H']
    pairs=[(a,b) for ia,a in enumerate(heavy) for b in heavy[ia+1:] if abs(keys[a][0]-keys[b][0])>1]
    def residual(x):
        xyz=apply(x);out=[]
        for o,p,c,o5 in links:
            out.append(cfg.bond_weight*(np.linalg.norm(xyz[o]-xyz[p])-cfg.bond_target))
            out.append(cfg.angle_weight*(_angle(xyz[c],xyz[o],xyz[p])-np.deg2rad(119)))
            out.append(cfg.angle_weight*(_angle(xyz[o],xyz[p],xyz[o5])-np.deg2rad(104)))
        out.extend((cfg.anchor_weight*(xyz[anchors]-start[anchors])).ravel())
        for a,b in pairs:
            d=np.linalg.norm(xyz[a]-xyz[b])
            out.append(cfg.clash_weight*min(d-cfg.clash_distance,0.0))
        return np.asarray(out,float)
    before=validate_atoms(residues)
    if not moves:
        return residues,{'success':True,'message':'no rotatable bonds','before':before,'after':before,'torsions':0,'experimental':True}
    bound=np.deg2rad(cfg.max_rotation_deg)
    result=least_squares(residual,np.zeros(len(moves)),bounds=(-bound,bound),max_nfev=cfg.max_nfev)
    coords=apply(result.x)
    out=[{} for _ in residues]
    for j,(i,k) in enumerate(keys):
        a=residues[i][k]
        out[i][k]=Atom(a.name,a.element,tuple(map(float,coords[j])))
    after=validate_atoms(out)
    return out,{'success':bool(result.success),'message':str(result.message),'before':before,'after':after,
       'torsions':len(moves),'torsion_labels':[f'{label}:{a[0]+1}' for label,a,b,_ in moves],
       'max_rotation_deg':float(np.rad2deg(np.max(np.abs(result.x)))),
       'optimizer_evaluations':int(result.nfev),'objective_before':float(np.sum(residual(np.zeros(len(moves)))**2)),
       'objective_after':float(np.sum(result.fun**2)),
       'experimental':True,'method':'acyclic bond torsion rotations; no sugar pucker optimization or force field'}
