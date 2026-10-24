"""Independent RNA stereochemistry diagnostics; conservative, configurable and non-calibrated.

No dependence on reconstruction or refinement objectives. Uses atom names and coordinates.
Checks do not establish experimentally accurate RNA geometry.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np

# Approximate reference distances (angstrom), deliberately generous acceptance ranges.
BONDS = [("P","OP1",1.35,1.75),("P","OP2",1.35,1.75),("P","O5'",1.35,1.85),
 ("O5'","C5'",1.25,1.65),("C5'","C4'",1.35,1.70),("C4'","O4'",1.25,1.65),
 ("O4'","C1'",1.25,1.65),("C1'","C2'",1.35,1.70),("C2'","C3'",1.35,1.70),
 ("C3'","C4'",1.35,1.70),("C3'","O3'",1.25,1.65),("C2'","O2'",1.25,1.65)]
ANGLES = [("O3'","P","O5'",85,135), ("P","O5'","C5'",100,145),
 ("O5'","C5'","C4'",95,135), ("C5'","C4'","C3'",90,135),
 ("C4'","C3'","O3'",95,135)]
RING = ("O4'","C1'","C2'","C3'","C4'")
BACKBONE = {
 'alpha': ((-1,"O3'"),(0,'P'),(0,"O5'"),(0,"C5'")),
 'beta': ((0,'P'),(0,"O5'"),(0,"C5'"),(0,"C4'")),
 'gamma': ((0,"O5'"),(0,"C5'"),(0,"C4'"),(0,"C3'")),
 'delta': ((0,"C5'"),(0,"C4'"),(0,"C3'"),(0,"O3'")),
 'epsilon': ((0,"C4'"),(0,"C3'"),(0,"O3'"),(1,'P')),
 'zeta': ((0,"C3'"),(0,"O3'"),(1,'P'),(1,"O5'")),
}

def xyz(atom):
    return np.asarray(atom.xyz,dtype=float)

def distance(a,b):
    return float(np.linalg.norm(xyz(a)-xyz(b)))

def angle(a,b,c):
    v=xyz(a)-xyz(b); w=xyz(c)-xyz(b)
    n=np.linalg.norm(v)*np.linalg.norm(w)
    if n<1e-10:return float('nan')
    return math.degrees(math.acos(float(np.clip(np.dot(v,w)/n,-1,1))))

def dihedral(a,b,c,d):
    p=[xyz(t) for t in (a,b,c,d)]
    b0=-(p[1]-p[0]); b1=p[2]-p[1]; b2=p[3]-p[2]
    norm=np.linalg.norm(b1)
    if norm<1e-10:return float('nan')
    b1/=norm; v=b0-np.dot(b0,b1)*b1; w=b2-np.dot(b2,b1)*b1
    if np.linalg.norm(v)*np.linalg.norm(w)<1e-10:return float('nan')
    return math.degrees(math.atan2(np.dot(np.cross(b1,v),w),np.dot(v,w)))

def validate_stereochemistry(residues, sequence=None, *, clash_cutoff=1.8):
    """Return JSON-serializable diagnostics; any detected critical issue fails.

    Geometry reference intervals are screening tolerances, not a validated
    distribution of crystallographic RNA structures. Missing atoms fail closed.
    """
    issues=[]; torsions=[]; puckers=[]; chirality=[]; lengths=[]; angles=[]
    if sequence is not None and len(sequence)!=len(residues):
        raise ValueError('sequence and residue count differ')
    def issue(kind,i,detail,**kw):issues.append({'kind':kind,'residue':i+1,'detail':detail,**kw})
    for i,res in enumerate(residues):
        required={x for pair in BONDS for x in pair[:2]}|set(RING)|{"C1'","N9" if sequence and sequence[i] in 'AG' else 'N1'}
        # Do not demand both non-bridging oxygen naming conventions (OP1/OP2 vs O1P/O2P).
        for name in required:
            if name not in res:issue('missing_atom',i,name)
        for name,atom in res.items():
            if not np.isfinite(xyz(atom)).all():issue('nonfinite',i,name)
        for a,b,lo,hi in BONDS:
            if a in res and b in res:
                d=distance(res[a],res[b]);lengths.append({'residue':i+1,'atoms':[a,b],'distance':d})
                if not lo<=d<=hi:issue('bond_length',i,f'{a}-{b}',observed=d,range=[lo,hi])
        for a,b,c,lo,hi in ANGLES:
            if all(n in res for n in (a,b,c)):
                v=angle(res[a],res[b],res[c]);angles.append({'residue':i+1,'atoms':[a,b,c],'degrees':v})
                if not lo<=v<=hi:issue('bond_angle',i,f'{a}-{b}-{c}',observed=v,range=[lo,hi])
        if all(n in res for n in RING):
            vals=[dihedral(*(res[RING[(j+k)%5]] for k in range(4))) for j in range(5)]
            if not all(math.isfinite(x) for x in vals):issue('ring_degenerate',i,'ring torsions undefined')
            puckers.append({'residue':i+1,'ring_torsions_deg':vals})
        if all(n in res for n in ("C1'","C2'","C3'","O4'")):
            v=[xyz(res[n])-xyz(res["C1'"]) for n in ("C2'","C3'","O4'")]
            signed=float(np.linalg.det(np.stack(v)))
            chirality.append({'residue':i+1,'signed_volume':signed})
            if abs(signed)<0.2:issue('chirality_degenerate',i,'near-planar sugar stereocenter',observed=signed)
        for name,defs in BACKBONE.items():
            atoms=[]
            for offset,n in defs:
                j=i+offset
                if j<0 or j>=len(residues) or n not in residues[j]:break
                atoms.append(residues[j][n])
            if len(atoms)==4:
                v=dihedral(*atoms);torsions.append({'residue':i+1,'name':name,'degrees':v})
                if not math.isfinite(v):issue('torsion_degenerate',i,name)
        if i+1<len(residues):
            nxt=residues[i+1]
            if "O3'" in res and 'P' in nxt:
                d=distance(res["O3'"],nxt['P'])
                if not 1.35<=d<=1.85:issue('phosphodiester_distance',i,'O3\'-P(next)',observed=d)
            else:issue('missing_link',i,'O3\'-P(next)')
            if all(n in res for n in ("C3'","O3'")) and 'P' in nxt:
                v=angle(res["C3'"],res["O3'"],nxt['P'])
                if not 100<=v<=145:issue('phosphodiester_angle',i,'C3\'-O3\'-P(next)',observed=v)
    # Non-neighbor heavy-atom clashes; skip bonded neighbors and adjacent residues
    for i in range(len(residues)):
        for j in range(i+2,len(residues)):
            for an,a in residues[i].items():
                if a.element.upper() in ('H','D'):continue
                for bn,b in residues[j].items():
                    if b.element.upper() in ('H','D'):continue
                    d=distance(a,b)
                    if d<clash_cutoff:issue('nonlocal_clash',i,f'{an} vs residue {j+1} {bn}',observed=d)
    # Report chirality sign, but cannot certify absolute stereochemistry without calibrated reference.
    return {'passed':len(issues)==0,'issue_count':len(issues),'issues':issues[:250],
            'issues_truncated':len(issues)>250,'bond_lengths':lengths,'bond_angles':angles,
            'torsions':torsions,'sugar_pucker':puckers,'chirality':chirality,
            'limitations':['No calibrated sugar-pucker classification','No absolute stereochemical assignment',
                           'No full force-field validation','No complete bonded-graph clash exclusions']}
