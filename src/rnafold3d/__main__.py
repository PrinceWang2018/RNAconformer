import argparse,json
from dataclasses import asdict
from .core import parse
from .topology import build_topology
from .geometry import aform_three_site,stem_endpoint_frames
from .assembly import assemble
from .closure import close_assembly
from .loop_sites import build_full_chain_guides
from .validation import filter_candidates, ValidationConfig
from .energy import score_sites, EnergyConfig
from .calibration import regularize_loop_frames, validate_guide_geometry
from .sampling import SamplingConfig, run_annealing, run_replica_exchange, summarize, evaluate
from .ensemble import EnsembleConfig, select_ensemble
from .ensemble_sampling import PoolConfig, run_pool

def main():
    import sys
    if len(sys.argv) > 1 and sys.argv[1] in ('predict', 'validate', 'ss', 'import', 'view', 'export-view'):
        from .cli import main as v2_main
        raise SystemExit(v2_main(sys.argv[1:]))
    ap=argparse.ArgumentParser(description='RNAFold3D M1 geometry kernel (not an atomistic RNA predictor)')
    ap.add_argument('--sequence',required=True);ap.add_argument('--structure',required=True)
    ap.add_argument('--output',default='m1_output.json');ap.add_argument('--strict',action='store_true')
    ap.add_argument('--assemble',type=int,default=0,metavar='N',help='Generate N rigid-body topology hypotheses')
    ap.add_argument('--seed',type=int,default=0);ap.add_argument('--proposals',type=int,default=100)
    ap.add_argument('--score-energy',action='store_true',help='Report experimental M2.1 guide-site energies for complete chains only')
    ap.add_argument('--calibrate-guides',action='store_true',help='Regularize loop local frames and report three-site guide consistency')
    ap.add_argument('--filter',action='store_true',help='Run M1.5 conservative guide-coordinate filtering')
    ap.add_argument('--max-selected',type=int,default=20)
    ap.add_argument('--rmsd-threshold',type=float,default=3.0)
    ap.add_argument('--allow-incomplete',action='store_true')
    ap.add_argument('--sample',choices=['anneal','remc'])
    ap.add_argument('--sample-steps',type=int,default=100)
    ap.add_argument('--sample-temperature',type=float,default=3.0)
    ap.add_argument('--select-ensemble',action='store_true',help='Rank sampled guide states and select diverse valid candidates')
    ap.add_argument('--ensemble-max-clashes',type=int,default=0)
    ap.add_argument('--pool-sampling',choices=['anneal','remc'],help='M2.4: save cross-run snapshots and sample loop axial rotations')
    ap.add_argument('--snapshot-interval',type=int,default=10)
    ap.add_argument('--loop-move-probability',type=float,default=0.5)
    ap.add_argument('--loop-rotation-sigma',type=float,default=0.4)
    a=ap.parse_args();m=parse(a.sequence,a.structure);top=build_topology(m,strict=a.strict)
    obj=asdict(m);obj['topology']=asdict(top);obj['stem_geometry']=[]
    for s in m.stems:
        sites=aform_three_site(s)
        obj['stem_geometry'].append({'sites':[asdict(x) for x in sites], 'end_frames':stem_endpoint_frames(sites,s), 'geometry_status':'idealized_unvalidated'})
    if a.assemble:
        candidates=assemble(m,a.assemble,a.seed,a.proposals)
        obj['assemblies']=[]
        if a.sample: obj['m2_2_sampling']=[]
        ensemble_entries=[]
        if a.filter:
            cfg=ValidationConfig(rmsd_threshold=a.rmsd_threshold,require_complete=not a.allow_incomplete)
            selected,reports,_=filter_candidates(m,candidates,cfg,max_candidates=a.max_selected)
            obj['m1_5_filter']={'config':asdict(cfg),'reports':[asdict(r) for r in reports],
                                'selected_candidate_indices':[c.proposal for c in selected],
                                'coordinate_status':'unvalidated_guide_coordinates'}
        for candidate in candidates:
            if a.sample:
                cfg_s=SamplingConfig(steps=a.sample_steps,temperature=a.sample_temperature,seed=a.seed+candidate.proposal)
                try:
                    sampled=(run_annealing(m,candidate,cfg_s) if a.sample=='anneal' else run_replica_exchange(m,candidate,cfg_s))
                    obj['m2_2_sampling'].append({'candidate_index':candidate.proposal,'result':summarize(sampled,a.sample)})
                    for label in ('best','final'):
                        ensemble_entries.append({'source':f'{a.sample}:{candidate.proposal}:{label}', 'state':sampled[label]})
                except ValueError as exc:
                    obj['m2_2_sampling'].append({'candidate_index':candidate.proposal,'error':str(exc)})
            if a.select_ensemble and not a.sample:
                state=evaluate(m,candidate,seed=a.seed+candidate.proposal)
                if state is not None:ensemble_entries.append({'source':f'initial:{candidate.proposal}','state':state})
            full=build_full_chain_guides(m,candidate,seed=a.seed+candidate.proposal)
            closures=full.closures
            stem_ids={r for stem in m.stems for pair in stem.pairs for r in pair}
            calibrated=regularize_loop_frames(full.sites,stem_ids) if a.calibrate_guides else full.sites
            calibration=validate_guide_geometry(calibrated,stem_ids).as_dict() if a.calibrate_guides else None
            energy = score_sites(m,calibrated,EnergyConfig(),require_complete=True).as_dict() if a.score_energy and full.complete else None
            obj['assemblies'].append({'m2_1_energy':energy,'m2_1_geometry_calibration':calibration,'seed':candidate.seed,'candidate_index':candidate.proposal,
                'placement_score':candidate.placement_score,
                'rotations':candidate.rotations,'translations':candidate.translations,
                'stem_sites':[[asdict(s) for s in sites] for sites in candidate.stem_sites],
                'full_chain_guides':{'sites':[asdict(s) for s in calibrated], 'missing_residues':full.missing_residues, 'complete':full.complete, 'n_loop_sites':full.n_loop_sites, 'max_junction_link_error':full.max_junction_link_error, 'geometry_status':full.geometry_status},
                'loop_closures':[asdict(c) for c in closures],
                'all_reported_loops_closed':all(c.converged for c in closures),
                'geometry_status':'unvalidated_guide_coordinates'})
        if a.pool_sampling:
            pc=PoolConfig(steps=a.sample_steps,snapshot_interval=a.snapshot_interval,loop_move_probability=a.loop_move_probability,loop_rotation_sigma=a.loop_rotation_sigma,mode=a.pool_sampling,seed=a.seed)
            ec_pool=EnsembleConfig(max_selected=a.max_selected,rmsd_threshold=a.rmsd_threshold,max_clashes=a.ensemble_max_clashes)
            obj['m2_4_pool']=run_pool(m,candidates,pc,ensemble_config=ec_pool)
        if a.select_ensemble:
            ec=EnsembleConfig(max_selected=a.max_selected,rmsd_threshold=a.rmsd_threshold,max_clashes=a.ensemble_max_clashes)
            obj['m2_3_ensemble']=select_ensemble(m,ensemble_entries,ec)
    with open(a.output,'w',encoding='utf-8') as f:json.dump(obj,f,indent=2)
    print(f'{len(m.sequence)} nt; {len(m.stems)} stems; planar={top.supported_planar}; output={a.output}')
if __name__=='__main__':main()
