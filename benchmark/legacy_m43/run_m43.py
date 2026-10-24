"""M4.3: conservative reference preparation and honest benchmark execution.
No annotation is inferred from target 3D coordinates.
"""
import argparse,csv,hashlib,json,subprocess,sys
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from benchmark_dataset import curate

BASE={'A','C','G','U'}
def get(d,k):
    v=d.get(k,[])
    return v if isinstance(v,list) else [v]

def prepare_reference(cif,chain,out):
    d=MMCIF2Dict(str(cif))
    names=['label_asym_id','label_seq_id','label_comp_id','label_atom_id','Cartn_x','Cartn_y','Cartn_z','pdbx_PDB_model_num','label_alt_id','auth_asym_id','auth_seq_id','pdbx_PDB_ins_code']
    cols=[get(d,'_atom_site.'+n) for n in names]
    if not cols[0] or len(set(map(len,cols)))!=1:raise ValueError('Incomplete atom_site table')
    residues={}; mapping={}; first_model=next((m for m in cols[7] if m not in ('.','?')), '1')
    for label,pos,base,atom,x,y,z,model,alt,auth_chain,auth_pos,ins in zip(*cols):
        if label!=chain or model!=first_model or alt not in ('.','?','A','1') or not pos.isdigit():continue
        if base not in BASE:continue
        idx=int(pos); residues.setdefault(idx,{})[atom.replace('*',"'")]=(base,(x,y,z))
        mapping[idx]={'label_seq_id':idx,'auth_asym_id':auth_chain,'auth_seq_id':auth_pos,'insertion_code':ins}
    if not residues:raise ValueError('No canonical RNA atoms for selected chain')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    lines=[];serial=1; coverage=0
    for idx,atoms in sorted(residues.items()):
        if "C4'" in atoms:coverage+=1
        for atom,(base,xyz) in sorted(atoms.items()):
            if len(atom)>4:continue
            try:xx,yy,zz=map(float,xyz)
            except ValueError:continue
            element=atom[0] if atom[0] in 'CNOPSH' else 'C'
            lines.append(f'ATOM  {serial%100000:5d} {atom:4s} {base:>3s} A{idx:4d}    {xx:8.3f}{yy:8.3f}{zz:8.3f}{1.0:6.2f}{0.0:6.2f}          {element:>2s}')
            serial+=1
    (out/'reference.pdb').write_text('\n'.join(lines)+'\nEND\n')
    (out/'residue_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
    sha=hashlib.sha256(Path(cif).read_bytes()).hexdigest()
    report={'source_cif':str(cif),'source_sha256':sha,'label_asym_id':chain,'observed_canonical_residues':len(residues),'c4_coverage':coverage/len(residues),'contiguous_label_ids':sorted(residues)==list(range(min(residues),max(residues)+1)),'excluded_modified_nucleotides':True,'reference_pdb':'reference.pdb','residue_mapping':'residue_mapping.json','review_status':'unreviewed'}
    (out/'reference_qc.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

def run(root,manifest,offline=True,predictions=None):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    curate(manifest,root,download=not offline)
    rows=list(csv.DictReader((root/'inventory.csv').open()))
    output=[]
    for row in rows:
        item={'pdb_id':row['pdb_id'],'sample_id':row['sample_id'],'dataset_status':row['status'],'reference_ready':False,'benchmark_status':'not_run'}
        if row['status']=='needs_2d_and_review':
            try:
                sid=row['sample_id'];pdb=row['pdb_id'];ch=row['chain']
                rep=prepare_reference(root/'raw'/f'{pdb}.cif',ch,root/'samples'/sid/'reference')
                item.update(reference_ready=True,reference_qc=rep)
                ann=root/'samples'/sid/'annotations'/'independent_2d'/'annotation_qc.json'
                ora=root/'samples'/sid/'annotations'/'oracle_2d'/'annotation_qc.json'
                item['annotation_available']={'independent_2d':ann.exists(),'oracle_2d':ora.exists()}
                pred_dir=Path(predictions)/sid if predictions else None
                candidates=sorted(pred_dir.glob('*.pdb')) if pred_dir and pred_dir.exists() else []
                if candidates:
                    from rnafold3d.benchmark import benchmark
                    result=benchmark(root/'samples'/sid/'reference'/'reference.pdb',candidates,chain='A')
                    target=root/'samples'/sid/'benchmark.json';target.write_text(json.dumps(result,indent=2)+'\n')
                    item['benchmark_status']='evaluated' if result['n_evaluable'] else 'no_evaluable_predictions'
                    item['best_of_k_rmsd']=result['best_of_k_rmsd']
                else:item['benchmark_status']='no_prediction_pdbs'
            except Exception as e:item['benchmark_status']='reference_preparation_failed';item['error']=str(e)
        output.append(item)
    summary={'candidate_entries':len(rows),'references_prepared':sum(r['reference_ready'] for r in output),'benchmarks_evaluated':sum(r['benchmark_status']=='evaluated' for r in output),'mean_best_of_k_rmsd':None,'records':output,'notes':['No synthetic data are counted as experimental evidence.','No secondary structure annotations are inferred from target 3D structures.','References require manual review before scientific inclusion.']}
    vals=[r['best_of_k_rmsd'] for r in output if r.get('best_of_k_rmsd') is not None]
    if vals:summary['mean_best_of_k_rmsd']=sum(vals)/len(vals)
    (root/'m43_execution_report.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dataset',default='dataset');ap.add_argument('--manifest',default='candidate_manifest.json');ap.add_argument('--online',action='store_true');ap.add_argument('--predictions',default=None)
    a=ap.parse_args();print(json.dumps(run(a.dataset,a.manifest,not a.online,a.predictions),indent=2))
