"""Independent C4-prime RNA ensemble benchmark. No training or structure lookup."""
import argparse
import json
from pathlib import Path
import numpy as np


def read_pdb(path, chain=None):
    """Return ordered (chain, resseq, insertion, base, C4' coordinate) records.

    First model only. Alternate locations A/blank only. No residue renumbering.
    """
    seen = set()
    records = []
    mapping = {'A':'A','C':'C','G':'G','U':'U','RA':'A','RC':'C','RG':'G','RU':'U'}
    for line in Path(path).read_text().splitlines():
        if line.startswith('ENDMDL'): break
        if not line.startswith(('ATOM  ', 'HETATM')): continue
        if line[12:16].strip().replace('*', "'") != "C4'": continue
        if line[16].strip() not in ('', 'A'): continue
        ch = line[21].strip() or '_'
        if chain is not None and ch != chain: continue
        base = mapping.get(line[17:20].strip().upper())
        if base is None: continue
        try:
            rid=(ch,int(line[22:26]),line[26].strip())
            xyz=np.array([float(line[i:i+8]) for i in (30,38,46)],dtype=float)
        except ValueError: continue
        if not np.isfinite(xyz).all() or rid in seen: continue
        seen.add(rid)
        records.append((rid,base,xyz))
    if not records: raise ValueError(f'No RNA C4-prime coordinates in {path}')
    return records


def aligned_rmsd(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if a.shape != b.shape or a.ndim!=2 or a.shape[1]!=3 or not len(a):
        raise ValueError('Coordinates must be matching nonempty Nx3 arrays')
    ac=a-a.mean(0); bc=b-b.mean(0)
    u,s,vt=np.linalg.svd(ac.T@bc)
    d=np.linalg.det(u@vt)
    rot=u@np.diag([1.,1.,1. if d>=0 else -1.])@vt
    return float(np.sqrt(np.mean(np.sum((ac@rot-bc)**2,axis=1))))


def compare(ref,pred,min_coverage=0.8):
    r={k:(base,xyz) for k,base,xyz in ref}
    p={k:(base,xyz) for k,base,xyz in pred}
    common=[k for k in r if k in p and r[k][0]==p[k][0]]
    coverage=len(common)/len(r)
    if coverage < min_coverage or len(common)<3:
        return {'status':'insufficient_coverage','matched':len(common),'reference_nt':len(r),'coverage':coverage,'c4_rmsd':None}
    return {'status':'ok','matched':len(common),'reference_nt':len(r),'coverage':coverage,
            'c4_rmsd':aligned_rmsd([r[k][1] for k in common],[p[k][1] for k in common])}


def benchmark(reference,predictions,chain=None,min_coverage=0.8,success_threshold=5.0):
    if not predictions: raise ValueError('At least one prediction required')
    ref=read_pdb(reference,chain)
    results=[]; usable=[]
    for path in predictions:
        try:
            pred=read_pdb(path,chain)
            result=compare(ref,pred,min_coverage)
        except (OSError,ValueError) as exc:
            result={'status':'read_error','error':str(exc),'c4_rmsd':None}
        result['file']=str(path)
        results.append(result)
        if result['status']=='ok': usable.append(result)
    rmsds=[x['c4_rmsd'] for x in usable]
    return {'schema_version':'m4.1','metric':'C4-prime aligned RMSD (Angstrom)',
            'reference':str(reference),'reference_residues':len(ref),
            'n_predictions':len(predictions),'n_evaluable':len(usable),
            'best_of_k_rmsd':min(rmsds) if rmsds else None,
            'median_rmsd':float(np.median(rmsds)) if rmsds else None,
            'success_threshold_angstrom':success_threshold,
            'success_at_k':bool(rmsds and min(rmsds)<=success_threshold),
            'per_conformer':results,
            'limitations':['C4-prime RMSD is not eRMSD, INF, TM-score, or stereochemical validity',
                           'Matching uses exact chain/residue IDs and matching bases; no sequence alignment',
                           'Best-of-K uses all evaluable predictions, not top-ranked accuracy']}


def main():
    ap=argparse.ArgumentParser(description='Benchmark RNA conformer PDBs against one experimental reference')
    ap.add_argument('--reference',required=True)
    ap.add_argument('--predictions',nargs='+',required=True)
    ap.add_argument('--chain',default=None)
    ap.add_argument('--min-coverage',type=float,default=0.8)
    ap.add_argument('--success-threshold',type=float,default=5.0)
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    if not 0<args.min_coverage<=1: ap.error('min-coverage must be within (0,1]')
    result=benchmark(args.reference,args.predictions,args.chain,args.min_coverage,args.success_threshold)
    Path(args.output).write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(f"Evaluable {result['n_evaluable']}/{result['n_predictions']}; Best-of-K C4 RMSD: {result['best_of_k_rmsd']}")

if __name__=='__main__': main()
