from dataclasses import dataclass
import math
import numpy as np

@dataclass(frozen=True)
class Stem:
    pairs: tuple
    @property
    def length(self): return len(self.pairs)

@dataclass(frozen=True)
class Element:
    kind: str
    residues: tuple
    incident_stems: tuple

@dataclass(frozen=True)
class SecondaryModel:
    sequence: str
    structure: str
    pairs: tuple
    stems: tuple
    elements: tuple
    stem_adjacency: tuple


def parse(sequence, structure):
    seq = sequence.upper().replace('T','U').replace(' ','').replace('\n','')
    db = ''.join(structure.split())
    if not 1 <= len(seq) <= 200: raise ValueError('RNA length must be 1..200 nt')
    if len(seq) != len(db): raise ValueError('Sequence/structure length mismatch')
    if any(c not in 'ACGU' for c in seq): raise ValueError('Sequence must contain A/C/G/U only')
    brackets = {'(':')','[':']','{':'}','<':'>'}
    stacks = {k:[] for k in brackets}
    pairs=[]
    reverse={v:k for k,v in brackets.items()}
    for idx,ch in enumerate(db):
        if ch in stacks: stacks[ch].append(idx)
        elif ch in reverse:
            op=reverse[ch]
            if not stacks[op]: raise ValueError(f'Unmatched bracket at {idx+1}')
            pairs.append((stacks[op].pop(),idx))
        elif ch!='.': raise ValueError(f'Unsupported character {ch!r}')
    if any(stacks.values()): raise ValueError('Unclosed brackets')
    pairs=tuple(sorted(pairs))
    pairset=set(pairs)
    starts=sorted(p for p in pairs if (p[0]-1,p[1]+1) not in pairset)
    stems=[]
    for i,j in starts:
        run=[]
        while (i,j) in pairset:
            run.append((i,j)); i+=1;j-=1
        stems.append(Stem(tuple(run)))
    stems=tuple(stems)
    ownership={idx:k for k,stem in enumerate(stems) for p in stem.pairs for idx in p}
    unpaired=[i for i in range(len(seq)) if i not in ownership]
    # Unpaired runs in primary-sequence order; classification is based on incident stems.
    elements=[]; pos=0
    while pos<len(unpaired):
        run=[unpaired[pos]]; pos+=1
        while pos<len(unpaired) and unpaired[pos]==run[-1]+1:
            run.append(unpaired[pos]);pos+=1
        a,b=run[0],run[-1]
        left=ownership.get(a-1);right=ownership.get(b+1)
        adjacent=tuple(sorted({x for x in (left,right) if x is not None}))
        kind=('external' if len(adjacent)<2 else 'connector')
        if left is not None and left==right: kind='hairpin_or_internal'
        elements.append(Element(kind,tuple(run),adjacent))
    edges=set()
    for i in range(len(seq)-1):
        a=ownership.get(i);b=ownership.get(i+1)
        if a is not None and b is not None and a!=b: edges.add(tuple(sorted((a,b))))
    for el in elements:
        if len(el.incident_stems)==2: edges.add(el.incident_stems)
    return SecondaryModel(seq,db,pairs,stems,tuple(elements),tuple(sorted(edges)))


def aform_stem(stem, rise=2.81, twist_deg=32.7, radius=9.0, phase_deg=0.0):
    """Idealized phosphate-free C1'-like paired backbone guide points, not atomistic coordinates.
    Returns (residue_indices, xyz) in matching order. Geometry is a tunable initialization.
    """
    indices=[]; coords=[]
    twist=math.radians(twist_deg); phase=math.radians(phase_deg)
    for k,(i,j) in enumerate(stem.pairs):
        theta=k*twist+phase
        for index,offset in ((i,0.0),(j,math.pi)):
            t=theta+offset
            indices.append(index)
            coords.append([radius*math.cos(t),radius*math.sin(t),k*rise])
    return np.asarray(indices,dtype=int),np.asarray(coords,dtype=float)


def transform(xyz, rotation, translation):
    r=np.asarray(rotation,dtype=float);t=np.asarray(translation,dtype=float)
    if r.shape!=(3,3) or t.shape!=(3,): raise ValueError('Expected 3x3 rotation and 3-vector')
    if not np.allclose(r.T@r,np.eye(3),atol=1e-6) or not np.isclose(np.linalg.det(r),1,atol=1e-6):
        raise ValueError('Rotation must be proper orthogonal')
    return np.asarray(xyz)@r.T+t


def endpoint_reachable(start, end, n_links, max_step=7.5):
    if n_links<0: return False
    return np.linalg.norm(np.asarray(end)-np.asarray(start)) <= n_links*max_step+1e-8
