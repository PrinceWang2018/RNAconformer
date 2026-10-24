"""M2.3: geometry-gated ranking and diverse ensemble selection of guide-site states.

Scores are arbitrary model units; clusters are structural, not thermodynamic states.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from .validation import aligned_rmsd

@dataclass(frozen=True)
class EnsembleConfig:
    max_selected: int = 10
    rmsd_threshold: float = 3.0
    clash_distance: float = 4.0
    max_clashes: int = 0
    max_link_error: float = 0.1
    def __post_init__(self):
        if self.max_selected < 1 or self.rmsd_threshold <= 0 or self.clash_distance <= 0 or self.max_clashes < 0 or self.max_link_error < 0:
            raise ValueError('Invalid ensemble parameters')


def _geometry(model, sites, cfg):
    n=len(model.sequence)
    coords={}
    for s in sites:
        if s.residue in coords: return ['duplicate_residue'], None, None
        v=np.asarray(s.C4,dtype=float)
        if v.shape!=(3,) or not np.isfinite(v).all():return ['nonfinite_coordinates'],None,None
        coords[s.residue]=v
    if set(coords)!=set(range(n)):return ['incomplete_chain'],None,None
    owner={r:k for k,stem in enumerate(model.stems) for pair in stem.pairs for r in pair}
    # Same-stem local guide lengths are provisional and not treated as calibrated bonds.
    errors=[abs(np.linalg.norm(coords[i+1]-coords[i])-6.) for i in range(n-1)
            if not (i in owner and i+1 in owner and owner[i]==owner[i+1])]
    clashes=sum(np.linalg.norm(coords[i]-coords[j]) < cfg.clash_distance
                for i in range(n) for j in range(i+3,n)
                if not (i in owner and j in owner and owner[i]==owner[j]))
    reasons=[]
    if errors and max(errors)>cfg.max_link_error: reasons.append('guide_link_mismatch')
    if clashes>cfg.max_clashes:reasons.append('nonlocal_clashes')
    return reasons, coords, {'max_link_error':max(errors,default=0.),'n_clashes':int(clashes)}


def select_ensemble(model, entries, config=EnsembleConfig()):
    """entries: iterable of {source, state}; returns JSON-serializable report.

    Every entry is inspected. Ranking is energy-first; distinct conformers
    are chosen by Kabsch-aligned C4' RMSD. Duplicate conformers remain in audit.
    """
    records=[]
    for index,entry in enumerate(entries):
        state=entry['state']
        reasons,coords,geometry=_geometry(model,state.sites,config)
        energy=float(state.energy)
        if not math.isfinite(energy):reasons.append('nonfinite_energy')
        records.append({'candidate_id':index,'source':str(entry['source']), 'energy':energy if math.isfinite(energy) else None,
                        'energy_breakdown':state.breakdown,'valid':not reasons,'reasons':reasons,
                        'geometry':geometry,'_coords':coords,'_sites':state.sites})
    ordered=sorted((r for r in records if r['valid']),key=lambda r:(r['energy'],r['candidate_id']))
    representatives=[]
    for r in ordered:
        d=[aligned_rmsd(r['_coords'],x['_coords']) for x in representatives]
        if d and min(d)<config.rmsd_threshold:
            r['cluster_id']=int(np.argmin(d))
            r['nearest_representative_rmsd']=float(min(d))
            r['reasons']=['duplicate_cluster']
        elif len(representatives)<config.max_selected:
            r['cluster_id']=len(representatives)
            r['nearest_representative_rmsd']=float(min(d)) if d else None
            representatives.append(r)
        else:
            r['reasons']=['selection_limit']
    for r in records:
        r.pop('_coords',None)
        sites=r.pop('_sites')
        if r in representatives:
            r['selected']=True
            r['sites']=[asdict(s) for s in sites]
        else:r['selected']=False
    return {'config':asdict(config),'input_count':len(records),'geometry_valid_count':len(ordered),
            'selected_count':len(representatives),'selected_ids':[r['candidate_id'] for r in representatives],
            'records':records,'coordinate_status':'unvalidated_three_site_guides',
            'score_units':'arbitrary; not a thermodynamic free energy',
            'selection_note':'Energy-first greedy RMSD clustering; not Boltzmann sampling'}
