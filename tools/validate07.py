"""Validation of the RCV 0.7 analyses on the benchmark ensembles (docs/V07_VALIDATION.md).

For every benchmark target the engine's exported models of both modes (benchmark/runs/v2_0_1_pilot/{independent,
oracle}/<target>/conformations, 5 + 5 models of one sequence) form one ensemble. On these ensembles:

  A. ensemble map vs classical MDS: stress, Spearman correlation of distances, trustworthiness, continuity and
     topology purity of map neighbours;
  B. pathways between the two most populated 2D states: direct vs best-route vs effective (network) barrier, number
     of intermediates, flux carried by the top route, reproducibility, run time;
  C. phi-value scan: are the kinetic hot spots (largest |ddG_TS|) at the cooperative unit?
  D. mechanism morph: chain distortion (largest relative change of consecutive C4'-C4' distances) of linear
     interpolation, the rigid-module morph and the staged mechanism morph (Node.js);
  E. probing data: synthetic SHAPE-like profiles from a known mixture of conformations (noise 0.1-0.4) - does the
     maximum-entropy reweighting recover the mixture?

    PYTHONPATH=src python tools/validate07.py --out docs/benchmark_v1/v07_validation.json
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from rcv import landscape, probing                       # noqa: E402
from rcv.landscape import classical_mds, quality         # noqa: E402
from rcv.conditions import Conditions                    # noqa: E402
from rcv.data import load_document                       # noqa: E402
from rcv.rcvfile import doc_coords                       # noqa: E402
from rcv.switch import ensemble_states                   # noqa: E402
from rcv import transition as T                          # noqa: E402


def med(v):
    v = [x for x in v if x is not None]
    return round(float(np.median(v)), 4) if v else None


def map_quality(doc):
    L = doc['ensemble']['landscape']
    out = {'classes': len(L['classes']), 'states': len(L['states']), 'map': L['quality']['map'], 'mds': L['quality']['mds']}
    # do the 2D states stay together inside the topology disks? (purity of map neighbours by 2D state)
    st = [None] * len(doc['conformations'])
    for i, s in enumerate(ensemble_states(doc)):
        for k in s['conformations']:
            st[k] = i
    if len(set(st)) > 1 and L['quality']['map']:
        D = np.asarray(doc['ensemble']['rmsd_c4'], float)
        out['state_purity'] = {'map': quality(D, np.asarray(L['positions'], float), st)['topology_purity'],
                               'mds': quality(D, classical_mds(D), st)['topology_purity']}
    return out


def pathway_case(doc):
    st = ensemble_states(doc)
    if len(st) < 2:
        return None
    # the two most different 2D states (largest base-pair distance; the more populated one is A)
    i, j = max(((i, j) for i in range(len(st)) for j in range(i + 1, len(st))), key=lambda ij: (T.bp_distance(st[ij[0]]['dbn'], st[ij[1]]['dbn']), len(st[ij[0]]['conformations']) + len(st[ij[1]]['conformations'])))
    a, b = st[i]['conformations'][0], st[j]['conformations'][0]
    t0 = time.time()
    r = T.analyse(doc, a, b, Conditions.from_json({}), coords=doc_coords(doc))
    dt = time.time() - t0
    again = T.analyse(doc, a, b, Conditions.from_json({}))
    paths = (r['tpt'] or {}).get('pathways', [])
    ab = next((e['saddle'] for e in r['edges'] if {e['a'], e['b']} == {0, 1}), None)
    out = {'a': a, 'b': b, 'bp_distance': r['bp_distance'], 'dG': r['dG'], 'direct': r['direct']['barrier'],
           'same_basin': r['same_basin'], 'above_basin_a': r['A']['above_basin'], 'above_basin_b': r['B']['above_basin'],
           'direct_between_basins': round(ab - r['A']['basin_energy'], 2) if ab is not None else None,
           'best': r['best']['route_barrier'], 'effective': r['effective_barrier'], 'nodes': len(r['nodes']),
           'route_legs': len(r['best']['route']) - 1, 'top_route_flux': paths[0]['flux_share'] if paths else None,
           'trap': bool(r['trap']), 'modules': len(r['motion']['modules']) if r.get('motion') else 0,
           'unit_bp': len(r['bottleneck']['unit']) if r['bottleneck'] else 0, 'seconds': round(dt, 2),
           'reproducible': again['effective_barrier'] == r['effective_barrier']}
    # phi scan: kinetic hot spots at the cooperative unit?
    if r['bottleneck'] and r['bottleneck']['unit']:
        t0 = time.time()
        sc = T.phi_scan(doc['sequence'], [r['nodes'][k]['dbn'] for k in r['best']['route']], Conditions.from_json({}))
        out['scan_seconds'] = round(time.time() - t0, 2)
        unit = {x for p in r['bottleneck']['unit'] for x in p}
        sens = np.array([p['max_abs_ddg_ts'] for p in sc['positions']])
        top = set(np.argsort(-sens)[:max(1, len(unit))].tolist())
        out['hotspots_in_unit'] = round(len(top & unit) / max(1, len(top)), 3)
        out['hotspots_expected_by_chance'] = round(len(unit) / len(sens), 3)
        order = np.argsort(np.argsort(sens)) / max(1, len(sens) - 1)      # 0 = least, 1 = most sensitive position
        out['unit_sensitivity_percentile'] = round(float(order[sorted(unit)].mean()), 3)
        out['mean_sens_unit'] = round(float(sens[sorted(unit)].mean()), 3)
        out['mean_sens_other'] = round(float(np.delete(sens, sorted(unit)).mean()), 3)
    out['_mech'] = {'motion': r.get('motion'), 'timing': r['timing']}
    return out


def probing_recovery(doc, rng, noises=(0.1, 0.2, 0.4)):
    U = probing.unpaired_matrix(doc)
    K = len(U)
    st = ensemble_states(doc)
    if K < 3 or len(st) < 2:
        return None
    rows = []
    for noise in noises:
        # true mixture: 70 % on the conformations of one 2D state, 30 % spread over the others
        major = st[1]['conformations']
        w = np.full(K, 0.3 / max(1, K - len(major)))
        w[major] = 0.7 / len(major)
        y = 0.2 + 1.2 * (w @ U) + rng.normal(0, noise, U.shape[1])
        r = probing.reweight(doc, list(np.clip(y, 0, None)))
        wf = np.array(r['weights'])
        rows.append({'noise': noise, 'true_major_share': 0.7, 'prior_major_share': round(len(major) / K, 3),
                     'fitted_major_share': round(float(wf[major].sum()), 3),
                     'weight_correlation': round(float(np.corrcoef(w, wf)[0, 1]), 3) if wf.std() > 0 else None,
                     'chi2_prior': r['chi2_prior'], 'chi2': r['chi2']})
    return rows


def morph_strain(cases, docs):
    """Chain distortion of linear, rigid and staged morphs (Node.js, rcv mechmorph.js)."""
    if not shutil.which('node'):
        return None
    script = r"""
import { readFileSync } from 'node:fs';
const [mm, ensMod, file] = process.argv.slice(2);
const { buildMechanism, chainStrain } = await import(mm);
const { Ensemble } = await import(ensMod);
const out = [];
for (const c of JSON.parse(readFileSync(file, 'utf8'))) {
  const ens = new Ensemble(c.doc);
  const A = ens.coords[c.a], B = ens.coords[c.b], trace = [...ens.trace];
  const lin = [], rig = [], stg = [];
  const mr = buildMechanism(A, B, ens.atomRes, c.motion, c.timing, { staged: false, names: ens.atomNames });
  const ms = buildMechanism(A, B, ens.atomRes, c.motion, c.timing, { staged: true, names: ens.atomNames });
  for (let k = 0; k <= 24; k++) {
    const t = k / 24, L = new Float32Array(A.length);
    for (let i = 0; i < A.length; i++) L[i] = A[i] + (B[i] - A[i]) * t;
    lin.push(L); rig.push(mr.at(t)); stg.push(ms.at(t));
  }
  out.push({ target: c.target, linear: chainStrain(lin, trace), rigid: chainStrain(rig, trace), staged: chainStrain(stg, trace) });
}
console.log(JSON.stringify(out));
"""
    with tempfile.TemporaryDirectory() as tmp:
        f = Path(tmp) / 'cases.json'
        f.write_text(json.dumps([{'target': t, 'doc': docs[t], 'a': c['a'], 'b': c['b'], **c['_mech']} for t, c in cases.items() if c and c['_mech']['motion']]))
        s = Path(tmp) / 'strain.mjs'
        s.write_text(script)
        js = ROOT / 'src' / 'rcv' / 'static' / 'js'
        r = subprocess.run(['node', str(s), (js / 'mechmorph.js').as_uri(), (js / 'ensemble.js').as_uri(), str(f)], capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            print(r.stderr[-2000:])
            return None
        return json.loads(r.stdout)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', default=str(ROOT / 'benchmark' / 'runs' / 'v2_0_1_pilot'))
    ap.add_argument('--data', default=str(ROOT / 'benchmark' / 'data'))
    ap.add_argument('--out', default=str(ROOT / 'docs' / 'benchmark_v1' / 'v07_validation.json'))
    ap.add_argument('--targets', default='')
    a = ap.parse_args()
    run, data = Path(a.run), Path(a.data)
    ds = json.loads((data / 'dataset_v1.json').read_text())
    targets = [s['sample_id'] for s in sorted(ds['samples'], key=lambda s: s['length'])]
    if a.targets:
        targets = [t for t in targets if t in a.targets.split(',')]
    rng = np.random.default_rng(11)
    res = {'maps': {}, 'pathways': {}, 'probing': {}}
    docs = {}
    for t in targets:
        dirs = [run / m / t for m in ('independent', 'oracle') if (run / m / t / 'metadata' / 'manifest.json').exists()]
        if not dirs:
            continue
        t0 = time.time()
        files = [str(f) for d in dirs for f in sorted((d / 'conformations').glob('*.pdb'))]     # .pdb only (each model also has a .cif)
        doc = load_document(files, name=t)
        docs[t] = doc
        res['maps'][t] = map_quality(doc)
        res['pathways'][t] = pathway_case(doc)
        res['probing'][t] = probing_recovery(doc, rng)
        p = res['pathways'][t]
        print(f"{t:8s} n={len(doc['sequence']):3d} conf={len(doc['conformations'])} topologies={res['maps'][t]['classes']} "
              f"states={res['maps'][t]['states']} " + (('same basin' if p['same_basin'] else f"barrier {p['direct']} (basins {p['direct_between_basins']})->{p['effective']} nodes {p['nodes']} hot@unit {p.get('hotspots_in_unit')}") if p else 'one 2D state')
              + f" ({time.time() - t0:.0f} s)", flush=True)
    res['morph'] = morph_strain(res['pathways'], docs)
    for c in res['pathways'].values():
        if c:
            c.pop('_mech', None)
    # summary
    M = [m for m in res['maps'].values()]
    multi = [m for m in M if m['classes'] > 1]
    P0 = [p for p in res['pathways'].values() if p]
    P = [p for p in P0 if not p['same_basin']]                 # A and B in different basins: barriers and rates
    S = {
        'targets': len(M), 'ensembles_with_several_topologies': len(multi),
        'map_vs_mds': {key: {'map': med([m['map'][key] for m in M if m['map']]), 'mds': med([m['mds'][key] for m in M if m['mds']]),
                             'map_better_or_equal': sum(1 for m in M if m['map'] and m['mds'] and ((m['map'][key] <= m['mds'][key]) if key == 'stress' else (m['map'][key] >= m['mds'][key])))}
                       for key in ('stress', 'spearman', 'trustworthiness', 'continuity', 'topology_purity')},
        'state_purity': {'map': med([m['state_purity']['map'] for m in M if m.get('state_purity')]),
                         'mds': med([m['state_purity']['mds'] for m in M if m.get('state_purity')]),
                         'cases': sum(1 for m in M if m.get('state_purity'))},
        'map_vs_mds_several_topologies': {key: {'map': med([m['map'][key] for m in multi if m['map']]), 'mds': med([m['mds'][key] for m in multi if m['mds']])}
                                          for key in ('trustworthiness', 'continuity', 'topology_purity', 'stress')},
        'pathways': {'pairs': len(P0), 'same_basin': len(P0) - len(P),
                     'endpoints_above_basin': sum(1 for p in P0 for x in (p['above_basin_a'], p['above_basin_b']) if x > 0.05),
                     'median_above_basin': med([x for p in P0 for x in (p['above_basin_a'], p['above_basin_b'])]),
                     'cases': len(P), 'median_direct': med([p['direct'] for p in P]), 'median_best': med([p['best'] for p in P]),
                     'median_direct_between_basins': med([p['direct_between_basins'] for p in P]),
                     'effective_below_direct_between_basins': sum(1 for p in P if p['effective'] < p['direct_between_basins'] - 0.05),
                     'median_effective': med([p['effective'] for p in P]),
                     'best_below_direct': sum(1 for p in P if p['best'] < p['direct'] - 0.05),
                     'effective_below_direct': sum(1 for p in P if p['effective'] < p['direct'] - 0.05),
                     'traps': sum(1 for p in P if p['trap']), 'reproducible': sum(1 for p in P0 if p['reproducible']),
                     'median_seconds': med([p['seconds'] for p in P]), 'max_seconds': max([p['seconds'] for p in P], default=None),
                     'median_scan_seconds': med([p.get('scan_seconds') for p in P]),
                     'hotspots_in_unit': med([p.get('hotspots_in_unit') for p in P]),
                     'unit_sensitivity_percentile': med([p.get('unit_sensitivity_percentile') for p in P]),
                     'unit_more_sensitive': sum(1 for p in P if p.get('mean_sens_unit') is not None and p['mean_sens_unit'] > p['mean_sens_other']),
                     'with_unit': sum(1 for p in P if p.get('mean_sens_unit') is not None),
                     'hotspots_by_chance': med([p.get('hotspots_expected_by_chance') for p in P]),
                     'sensitivity_unit_vs_other': [med([p.get('mean_sens_unit') for p in P]), med([p.get('mean_sens_other') for p in P])]},
    }
    if res['morph']:
        S['morph_strain'] = {k: med([x[k] for x in res['morph']]) for k in ('linear', 'rigid', 'staged')}
        S['morph_strain']['cases'] = len(res['morph'])
        S['morph_strain']['rigid_below_linear'] = sum(1 for x in res['morph'] if x['rigid'] < x['linear'])
        S['morph_strain']['staged_below_linear'] = sum(1 for x in res['morph'] if x['staged'] < x['linear'])
    rec = [r for rows in res['probing'].values() if rows for r in rows]
    S['probing'] = {str(nz): {'median_fitted_major_share': med([r['fitted_major_share'] for r in rec if r['noise'] == nz]),
                              'median_prior_major_share': med([r['prior_major_share'] for r in rec if r['noise'] == nz]),
                              'median_weight_correlation': med([r['weight_correlation'] for r in rec if r['noise'] == nz]),
                              'cases': sum(1 for r in rec if r['noise'] == nz)} for nz in (0.1, 0.2, 0.4)}
    res['summary'] = S
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps(S, indent=1))


if __name__ == '__main__':
    main()
