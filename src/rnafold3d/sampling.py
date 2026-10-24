"""M2.2 experimental Metropolis / annealing / replica-exchange guide sampling.

Energy is in arbitrary units. Temperature is a dimensionless search parameter;
these trajectories are NOT thermodynamic ensembles. Every topology proposal
rebuilds loops and rejects incomplete or non-finite candidates.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from .assembly import Assembly, _rotation
from .geometry import rigid_transform_sites
from .loop_sites import build_full_chain_guides
from .calibration import regularize_loop_frames, validate_guide_geometry
from .energy import score_sites, EnergyConfig

@dataclass(frozen=True)
class SamplingConfig:
    steps: int = 100
    temperature: float = 3.0
    final_temperature: float = 0.5
    translation_sigma: float = 0.6
    rotation_sigma: float = 0.12
    exchange_interval: int = 10
    temperatures: tuple = (0.7, 1.5, 3.0, 6.0)
    seed: int = 42
    def __post_init__(self):
        if self.steps < 0 or self.exchange_interval < 1 or not self.temperatures:
            raise ValueError('Invalid sampling iterations')
        for x in (self.temperature,self.final_temperature,self.translation_sigma,self.rotation_sigma,*self.temperatures):
            if not math.isfinite(x) or x <= 0: raise ValueError('Sampling parameters must be positive finite')

@dataclass
class State:
    assembly: Assembly
    energy: float
    breakdown: dict
    sites: tuple


def evaluate(model, assembly, *, seed=0, energy_config=None):
    """Return a valid complete three-site guide state, or None.

    All stems retain rigid geometry; loops are regenerated with closure and
    calibrated frames. Reject incomplete chains and severe guide discontinuities.
    """
    full=build_full_chain_guides(model,assembly,seed=seed)
    if not full.complete or len(full.sites)!=len(model.sequence):return None
    stem_ids={r for stem in model.stems for pair in stem.pairs for r in pair}
    sites=regularize_loop_frames(full.sites,stem_ids)
    report=validate_guide_geometry(sites,stem_ids)
    # Keep rejection conservative: continuity is the only hard requirement here.
    if full.max_junction_link_error > 0.051:return None
    try:
        result=score_sites(model,sites,energy_config or EnergyConfig(),require_complete=True)
    except ValueError:
        return None
    if not math.isfinite(result.total):return None
    return State(assembly,result.total,result.as_dict(),sites)


def perturb(assembly, rng, config):
    """Perturb one complete stem as a rigid body, preserving its local shape."""
    count=len(assembly.stem_sites)
    if count==0: return assembly
    k=int(rng.integers(count))
    rot=np.asarray(assembly.rotations[k],dtype=float)
    trans=np.asarray(assembly.translations[k],dtype=float)
    delta_r=_rotation(rng.normal(size=3),rng.normal(0,config.rotation_sigma))
    # Rotate about the stem centroid, not the global origin.
    center=np.mean([s.C4 for s in assembly.stem_sites[k]],axis=0)
    displacement=rng.normal(0,config.translation_sigma,3)
    new_sites=[]
    for s in assembly.stem_sites[k]:
        pts=np.asarray([s.P,s.C4,s.N]); new_sites.append(type(s)(s.residue,*[tuple(p) for p in ((pts-center)@delta_r.T+center+displacement)]))
    rotations=list(assembly.rotations); translations=list(assembly.translations); stems=list(assembly.stem_sites)
    rotations[k]=tuple(map(tuple,delta_r@rot))
    translations[k]=tuple(delta_r@(trans-center)+center+displacement)
    stems[k]=tuple(new_sites)
    return Assembly(tuple(stems),tuple(rotations),tuple(translations),assembly.seed,assembly.proposal,assembly.placement_score)


def _accept(old,new,temp,rng):
    return new <= old or rng.random() < math.exp(max(-700.0,(old-new)/temp))


def run_annealing(model, assembly, config=SamplingConfig(), *, energy_config=None):
    rng=np.random.default_rng(config.seed)
    state=evaluate(model,assembly,seed=config.seed,energy_config=energy_config)
    if state is None: raise ValueError('Initial assembly cannot form a complete valid guide chain')
    best=state; accepted=0; invalid=0; trace=[state.energy]
    for step in range(config.steps):
        frac=step/max(1,config.steps-1)
        temp=config.temperature*(config.final_temperature/config.temperature)**frac
        trial=evaluate(model,perturb(state.assembly,rng,config),seed=config.seed+step+1,energy_config=energy_config)
        if trial is None: invalid+=1
        elif _accept(state.energy,trial.energy,temp,rng):
            state=trial;accepted+=1
            if state.energy<best.energy:best=state
        trace.append(state.energy)
    return {'best':best,'final':state,'trace':trace,'accepted':accepted,'invalid':invalid,'attempted':config.steps}


def run_replica_exchange(model, assembly, config=SamplingConfig(), *, energy_config=None):
    """Fixed-temperature Metropolis chains with adjacent-temperature swaps."""
    rng=np.random.default_rng(config.seed)
    temps=tuple(config.temperatures)
    initial=evaluate(model,assembly,seed=config.seed,energy_config=energy_config)
    if initial is None:raise ValueError('Initial assembly cannot form a complete valid guide chain')
    states=[initial for _ in temps];best=initial
    accepted=[0]*len(temps);invalid=[0]*len(temps);swaps=0;swap_attempts=0
    traces=[[initial.energy] for _ in temps]
    for step in range(config.steps):
        for i,temp in enumerate(temps):
            trial=evaluate(model,perturb(states[i].assembly,rng,config),seed=config.seed+step*len(temps)+i+1,energy_config=energy_config)
            if trial is None:invalid[i]+=1
            elif _accept(states[i].energy,trial.energy,temp,rng):
                states[i]=trial;accepted[i]+=1
                if trial.energy<best.energy:best=trial
        if (step+1)%config.exchange_interval==0:
            parity=((step+1)//config.exchange_interval-1)%2
            for i in range(parity,len(temps)-1,2):
                j=i+1;swap_attempts+=1
                logp=(1/temps[i]-1/temps[j])*(states[i].energy-states[j].energy)
                if math.log(max(float(rng.random()),1e-300))<min(0.,logp):
                    states[i],states[j]=states[j],states[i];swaps+=1
        for i,state in enumerate(states):traces[i].append(state.energy)
    return {'best':best,'final':states[0],'traces':traces,'accepted':accepted,'invalid':invalid,
            'swaps':swaps,'swap_attempts':swap_attempts,'attempted_per_replica':config.steps}


def summarize(result, mode):
    keys=('accepted','invalid','attempted','swaps','swap_attempts','attempted_per_replica')
    return {'mode':mode,'best_energy':result['best'].energy,'final_energy':result['final'].energy,
            'best_breakdown':result['best'].breakdown,
            'best_sites':[asdict(s) for s in result['best'].sites],
            'best_assembly':{'rotations':result['best'].assembly.rotations,'translations':result['best'].assembly.translations},
            **{k:result[k] for k in keys if k in result},
            **({'energy_trace':result['trace']} if 'trace' in result else {'energy_traces':result['traces']})}
