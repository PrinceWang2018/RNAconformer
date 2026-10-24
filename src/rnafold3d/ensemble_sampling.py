"""M2.4: loop axial-rotation moves, snapshots and pooled ensemble sampling.

Loop moves preserve guide C4 link lengths and fixed stem anchors. These are
geometric moves, not RNA backbone torsion moves. Energies are arbitrary units.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from .sampling import SamplingConfig, evaluate, perturb, _accept, State
from .geometry import NucleotideSites
from .energy import score_sites, EnergyConfig
from .ensemble import EnsembleConfig, select_ensemble, _geometry

@dataclass(frozen=True)
class PoolConfig:
    steps: int = 100
    snapshot_interval: int = 10
    loop_move_probability: float = 0.5
    loop_rotation_sigma: float = 0.4
    temperature: float = 2.0
    final_temperature: float = 0.6
    mode: str = 'anneal'
    temperatures: tuple = (0.7, 1.5, 3.0, 6.0)
    exchange_interval: int = 10
    seed: int = 42
    def __post_init__(self):
        if self.steps < 0 or self.snapshot_interval < 1 or self.exchange_interval < 1 or self.mode not in ('anneal','remc'):
            raise ValueError('Invalid sampling configuration')
        if not 0 <= self.loop_move_probability <= 1 or any(not math.isfinite(x) or x<=0 for x in (self.loop_rotation_sigma,self.temperature,self.final_temperature,*self.temperatures)):
            raise ValueError('Invalid sampling parameters')
        if sorted(self.temperatures)!=list(self.temperatures):raise ValueError('Replica temperatures must be ordered')

def loop_runs(model):
    """Return internal unpaired runs bounded by stem residues."""
    stem={r for s in model.stems for pair in s.pairs for r in pair}
    runs=[]; i=0; n=len(model.sequence)
    while i<n:
        if i in stem:i+=1;continue
        j=i
        while j<n and j not in stem:j+=1
        if i>0 and j<n and i-1 in stem and j in stem:runs.append(tuple(range(i,j)))
        i=j
    return tuple(runs)

def rotate_loop(model,state,rng,sigma,energy_config=None):
    """Rigidly rotate one internal loop about its fixed stem-to-stem endpoint axis.

    This is a nontrivial loop-only conformational move that exactly preserves
    distances to both anchor C4 sites and within-loop virtual C4 bonds.
    """
    runs=loop_runs(model)
    if not runs:return None
    run=runs[int(rng.integers(len(runs)))]
    lookup={s.residue:s for s in state.sites}
    if any(i not in lookup for i in (*run,run[0]-1,run[-1]+1)):return None
    a=np.asarray(lookup[run[0]-1].C4);b=np.asarray(lookup[run[-1]+1].C4)
    axis=b-a;norm=np.linalg.norm(axis)
    if norm<1e-9:return None
    axis/=norm;angle=float(rng.normal(0,sigma))
    v=axis;K=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
    R=np.eye(3)+math.sin(angle)*K+(1-math.cos(angle))*(K@K)
    sites=[]
    for s in state.sites:
        if s.residue in run:
            xyz=[tuple(float(x) for x in (a+R@(np.asarray(pt)-a))) for pt in (s.P,s.C4,s.N)]
            sites.append(NucleotideSites(s.residue,*xyz))
        else:sites.append(s)
    sites=tuple(sites)
    # Re-score complete chain. No stochastic reclosure: move is reversible.
    try:energy=score_sites(model,sites,energy_config or EnergyConfig(),require_complete=True)
    except ValueError:return None
    if not math.isfinite(energy.total):return None
    return State(state.assembly,energy.total,energy.as_dict(),sites)

def _propose(model,state,rng,cfg,step,energy_config):
    if rng.random()<cfg.loop_move_probability:
        trial=rotate_loop(model,state,rng,cfg.loop_rotation_sigma,energy_config)
        return trial,'loop'
    sc=SamplingConfig(steps=1,seed=cfg.seed,translation_sigma=0.6,rotation_sigma=0.12)
    trial=evaluate(model,perturb(state.assembly,rng,sc),seed=cfg.seed+step+1,energy_config=energy_config)
    return trial,'stem'

def run_pool(model,assemblies,config=PoolConfig(),*,ensemble_config=EnsembleConfig(),energy_config=None):
    """Aggregate snapshots from multiple independent runs, then select diverse states."""
    rng=np.random.default_rng(config.seed)
    entries=[];stats=[]
    for run_id,assembly in enumerate(assemblies):
        initial=evaluate(model,assembly,seed=config.seed+run_id,energy_config=energy_config)
        if initial is None:
            stats.append({'run':run_id,'error':'invalid_initial_chain'});continue
        temperatures=config.temperatures if config.mode=='remc' else (config.temperature,)
        states=[initial for _ in temperatures]; best=initial
        accepted={'stem':0,'loop':0}; attempted={'stem':0,'loop':0}; invalid=0; swaps=0; swap_attempts=0
        def capture(step):
            for k,state in enumerate(states):entries.append({'source':f'run{run_id}:step{step}:replica{k}','state':state})
        capture(0)
        for step in range(config.steps):
            for k,t in enumerate(temperatures):
                trial,kind=_propose(model,states[k],rng,config,step+run_id*100000+k*1000,energy_config)
                attempted[kind]+=1
                if trial is None:invalid+=1;continue
                # M1.5 geometry gate for candidate snapshots / state transitions.
                reasons,_,_=_geometry(model,trial.sites,ensemble_config)
                if reasons:invalid+=1;continue
                temp=(config.temperature*(config.final_temperature/config.temperature)**(step/max(1,config.steps-1))) if config.mode=='anneal' else t
                if _accept(states[k].energy,trial.energy,temp,rng):
                    states[k]=trial;accepted[kind]+=1
                    if trial.energy<best.energy:best=trial
            if config.mode=='remc' and (step+1)%config.exchange_interval==0:
                parity=((step+1)//config.exchange_interval-1)%2
                for k in range(parity,len(states)-1,2):
                    swap_attempts+=1
                    logp=(1/temperatures[k]-1/temperatures[k+1])*(states[k].energy-states[k+1].energy)
                    if math.log(max(float(rng.random()),1e-300))<min(0.,logp):
                        states[k],states[k+1]=states[k+1],states[k];swaps+=1
            if (step+1)%config.snapshot_interval==0:capture(step+1)
        entries.append({'source':f'run{run_id}:best','state':best})
        stats.append({'run':run_id,'attempted':attempted,'accepted':accepted,'invalid':invalid,'swaps':swaps,'swap_attempts':swap_attempts,'best_energy':best.energy})
    report=select_ensemble(model,entries,ensemble_config)
    report['m2_4']={'config':asdict(config),'runs':stats,'snapshot_count':len(entries),'sampling_note':'Guide-space geometric search; snapshots are not equilibrium populations'}
    return report
