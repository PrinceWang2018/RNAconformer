"""M3 experimental template-fitted atomistic export.

Requires user-provided chemically complete RNA residue templates. Never claims
that independent per-residue fitting yields a stereochemically valid RNA chain.
"""
from dataclasses import dataclass, asdict
from pathlib import Path
import json
import numpy as np

ANCHORS = ("P", "C4'", "N")
PURINES = frozenset('AG')
BASE_N = {'A':'N9','G':'N9','C':'N1','U':'N1'}

@dataclass(frozen=True)
class Atom:
    name: str
    element: str
    xyz: tuple


def parse_template_pdb(path, base):
    """Read a single complete residue from a supplied PDB, rejecting ambiguity."""
    atoms={}
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM  '):continue
        name=line[12:16].strip(); residue=line[17:20].strip().upper()
        if residue != base:continue
        if line[16:17] not in (' ','A'):continue
        if name in atoms:raise ValueError(f'duplicate atom {name} in template {base}')
        try: xyz=tuple(float(line[i:i+8]) for i in (30,38,46))
        except ValueError as e:raise ValueError('malformed template coordinates') from e
        element=line[76:78].strip().upper() or ''.join(x for x in name if x.isalpha())[0]
        atoms[name]=Atom(name,element,xyz)
    needed={"P","C4'",BASE_N[base],"O3'","O5'","C3'","C5'","C1'"}
    if not needed.issubset(atoms):raise ValueError(f'template {base} missing {sorted(needed-set(atoms))}')
    if len(atoms)<18:raise ValueError('incomplete RNA template')
    return atoms


def fit_residue(atoms, site, base):
    names=['P',"C4'",BASE_N[base]]
    a=np.array([atoms[x].xyz for x in names],float)
    b=np.array([site.P,site.C4,site.N],float)
    if not np.isfinite(b).all() or np.linalg.matrix_rank(a-a.mean(axis=0))<2 or np.linalg.matrix_rank(b-b.mean(axis=0))<2:
        raise ValueError('degenerate anchor geometry')
    ac=a.mean(axis=0);bc=b.mean(axis=0)
    u,_,vt=np.linalg.svd((a-ac).T@(b-bc))
    d=np.eye(3);d[2,2]=np.linalg.det(u@vt)
    rot=u@d@vt
    result={k:Atom(v.name,v.element,tuple((np.array(v.xyz)-ac)@rot+bc)) for k,v in atoms.items()}
    residual=float(np.sqrt(np.mean(np.sum(((a-ac)@rot+bc-b)**2,axis=1))))
    return result,residual


def reconstruct(sequence, sites, templates, max_anchor_rmsd=1.5):
    """Template-fit each residue; fail if anchors cannot be fitted closely."""
    lookup={s.residue:s for s in sites}
    if set(lookup)!=set(range(len(sequence))) or len(sites)!=len(sequence):
        raise ValueError('complete, unique guide sites required')
    out=[]; residuals=[]
    for i,base in enumerate(sequence):
        if base not in 'ACGU':raise ValueError('only canonical RNA bases supported')
        if base not in templates:raise ValueError(f'missing template for {base}')
        atoms,err=fit_residue(templates[base],lookup[i],base)
        if err>max_anchor_rmsd:raise ValueError(f'residue {i+1} anchor RMSD {err:.2f} exceeds threshold')
        out.append(atoms);residuals.append(err)
    return out,residuals


def validate_atoms(residues, clash_distance=1.7, linkage_min=1.3, linkage_max=1.9):
    """Conservative diagnostics, not complete stereochemical validation."""
    broken=[];links=[]
    for i in range(len(residues)-1):
        a=np.array(residues[i]["O3'"].xyz);b=np.array(residues[i+1]['P'].xyz)
        d=float(np.linalg.norm(a-b));links.append(d)
        if not linkage_min<=d<=linkage_max:broken.append({'after_residue':i+1,'distance':d})
    clashes=[]
    for i,atoms1 in enumerate(residues):
        for j in range(i+2,len(residues)):
            for a in atoms1.values():
                for b in residues[j].values():
                    if a.element=='H' or b.element=='H':continue
                    d=float(np.linalg.norm(np.array(a.xyz)-b.xyz))
                    if d<clash_distance:clashes.append((i+1,j+1,a.name,b.name,round(d,3)))
    return {'link_distances':links,'broken_links':broken,'nonlocal_clashes':clashes[:100],
            'n_nonlocal_clashes':len(clashes),'passed_basic_checks':not broken and not clashes,
            'note':'Additional bond-angle, torsion, sugar-pucker and force-field checks required.'}


def write_pdb(path,sequence,residues,report,allow_invalid=False):
    if not report['passed_basic_checks'] and not allow_invalid:
        raise ValueError('invalid atomistic model; export blocked (use allow_invalid only for debugging)')
    lines=['REMARK 900 RNAFOLD3D EXPERIMENTAL TEMPLATE-FITTED MODEL',
           'REMARK 900 NOT VALIDATED AS AN RNA TERTIARY STRUCTURE']
    idx=1
    for i,(base,atoms) in enumerate(zip(sequence,residues),1):
        for atom in atoms.values():
            x,y,z=atom.xyz
            if i>9999 or idx>99999:raise ValueError('PDB numbering limit; use mmCIF')
            lines.append(f'ATOM  {idx:5d} {atom.name:4s} {base:>3s} A{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{0.:6.2f}          {atom.element:>2s}')
            idx+=1
    lines+=['TER','END']
    Path(path).write_text('\n'.join(lines)+'\n')


def export_ensemble(sequence, conformers, template_paths, outdir, allow_invalid=False, refine=False, refine_torsion=False, refine_pucker=False):
    """Conformers are lists of NucleotideSites. Save manifest and validated PDBs."""
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    templates={base:parse_template_pdb(path,base) for base,path in template_paths.items()}
    manifest={'format_version':'0.1','sequence':sequence,'method':'independent residue template fitting',
              'experimental':True,'conformers':[]}
    for k,sites in enumerate(conformers,1):
        entry={'id':k,'status':'rejected'}
        try:
            residues,residuals=reconstruct(sequence,sites,templates)
            if refine:
                from .refinement import refine_rigid_chain
                residues,refinement=refine_rigid_chain(residues)
                entry['refinement']=refinement
            if refine_torsion:
                from .torsion_refinement import refine_torsions
                residues,torsion_report=refine_torsions(residues,sequence)
                entry['torsion_refinement']=torsion_report
            if refine_pucker:
                from .sugar_pucker import refine_chain_puckers
                residues,pucker_reports=refine_chain_puckers(residues)
                entry['sugar_pucker_refinement']=pucker_reports
            report=validate_atoms(residues)
            from .stereochemistry import validate_stereochemistry
            stereo=validate_stereochemistry(residues, sequence)
            entry["stereochemistry"]=stereo
            report["passed_basic_checks"]=report["passed_basic_checks"] and stereo["passed"]
            entry.update({'anchor_rmsd_max':max(residuals),'validation':report})
            if report['passed_basic_checks'] or allow_invalid:
                name=f'conf_{k:03d}.pdb';write_pdb(out/name,sequence,residues,report,allow_invalid)
                entry.update({'file':name,'status':'passed_basic_checks' if report['passed_basic_checks'] else 'invalid_debug_only'})
        except (ValueError,KeyError) as e:entry['error']=str(e)
        manifest['conformers'].append(entry)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
    return manifest
