"""CLI for optional, explicitly template-based experimental atomistic reconstruction."""
import argparse,json
from pathlib import Path
from .geometry import NucleotideSites
from .atomistic import export_ensemble

def main():
    p=argparse.ArgumentParser(description='Experimental RNA atomistic reconstruction from supplied residue templates')
    p.add_argument('--input',required=True,help='M2 JSON with selected ensemble or full-chain guides')
    p.add_argument('--template-a');p.add_argument('--template-c');p.add_argument('--template-g');p.add_argument('--template-u')
    p.add_argument('--outdir',required=True)
    p.add_argument('--allow-invalid-debug',action='store_true')
    p.add_argument('--refine-rigid-chain',action='store_true',help='Experimental chain-aware rigid residue refinement')
    p.add_argument('--refine-torsions',action='store_true',help='Experimental local acyclic backbone torsion refinement (max 20 nt)')
    p.add_argument('--refine-sugar-pucker',action='store_true',help='Experimental restrained ribose ring fitting; does not ensure stereochemical validity')
    a=p.parse_args()
    obj=json.loads(Path(a.input).read_text());seq=obj['sequence']
    raw=[]
    for key in ('m2_4_pool','m2_3_ensemble'):
        if key in obj:
            raw.extend(r['sites'] for r in obj[key].get('records',[]) if r.get('selected') and r.get('sites'))
    if not raw:
        raw.extend(x['full_chain_guides']['sites'] for x in obj.get('assemblies',[]) if x.get('full_chain_guides',{}).get('complete'))
    if not raw:raise SystemExit('No complete guide conformations found')
    paths={k:getattr(a,'template_'+k.lower()) for k in 'ACGU' if getattr(a,'template_'+k.lower())}
    conformers=[[NucleotideSites(**s) for s in conf] for conf in raw]
    result=export_ensemble(seq,conformers,paths,a.outdir,a.allow_invalid_debug,refine=a.refine_rigid_chain,refine_torsion=a.refine_torsions,refine_pucker=a.refine_sugar_pucker)
    print(f"{len(result['conformers'])} candidates; {sum(c['status']=='passed_basic_checks' for c in result['conformers'])} passed basic atomistic checks")
if __name__=='__main__':main()
