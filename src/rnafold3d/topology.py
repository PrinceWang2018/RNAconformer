"""Secondary-structure tree and junction graph (0-based indices).

Crossing pairs are detected and explicitly rejected for planar junction decomposition.
"""
from dataclasses import dataclass, asdict
from .core import parse

@dataclass(frozen=True)
class Junction:
    kind: str
    parent_stem: int | None
    child_stems: tuple
    segments: tuple
    incident_stems: tuple

@dataclass(frozen=True)
class Topology:
    junctions: tuple
    edges: tuple
    pseudoknot_pairs: tuple
    supported_planar: bool

def build_topology(model, strict=False):
    pairs=model.pairs
    crossing=tuple((a,b) for a in pairs for b in pairs if a[0]<b[0]<a[1]<b[1])
    if crossing:
        if strict: raise ValueError('Crossing base pairs require pseudoknot-aware topology; planar junction graph unavailable')
        return Topology((),(),tuple(sorted(set(p for ab in crossing for p in ab))),False)
    stems=model.stems
    if not stems:
        return Topology((Junction('external',None,(),(tuple(range(len(model.sequence))),),()),),(),(),True)
    bounds=[(s.pairs[0][0],s.pairs[0][1],s.pairs[-1][0],s.pairs[-1][1]) for s in stems]
    parent={}
    for i,(l,r,_,_) in enumerate(bounds):
        containers=[j for j,(L,R,_,_) in enumerate(bounds) if j!=i and L<l<r<R]
        parent[i]=min(containers,key=lambda j:bounds[j][1]-bounds[j][0]) if containers else None
    children={i:[] for i in range(len(stems))}
    roots=[]
    for i,p in parent.items():
        (roots if p is None else children[p]).append(i)
    for v in children.values():v.sort(key=lambda i:bounds[i][0])
    roots.sort(key=lambda i:bounds[i][0])
    junctions=[]; edges=[]
    def segments_between(lo,hi,kids):
        out=[]; cursor=lo
        for k in kids:
            left,right,_,_=bounds[k]
            out.append(tuple(range(cursor,left)))
            cursor=right+1
        out.append(tuple(range(cursor,hi+1)))
        return tuple(out)
    junctions.append(Junction('external',None,tuple(roots),segments_between(0,len(model.sequence)-1,roots),tuple(roots)))
    for i,(left,right,inner_left,inner_right) in enumerate(bounds):
        kids=tuple(children[i]); segs=segments_between(inner_left+1,inner_right-1,kids)
        if len(kids)==0: kind='hairpin'
        elif len(kids)==1:
            a,b=segs
            kind='stack_continuation' if not a and not b else ('bulge' if not a or not b else 'internal_loop')
        else:kind='multibranch_junction' if len(kids)>=2 else 'unknown'
        junctions.append(Junction(kind,i,kids,segs,(i,)+kids))
        edges.extend((i,k) for k in kids)
    return Topology(tuple(junctions),tuple(edges),(),True)

def topology_dict(t):return asdict(t)
