"""Command-line interface for the v2 engine.

    python -m rnafold3d predict  --fasta x.fa [--dbn x.dbn] --models 5 --seed 42 --out run_001 [--workers 8]
    python -m rnafold3d validate run_001/conformations/conf_001.pdb [more.pdb ...]
    python -m rnafold3d ss       --fasta x.fa            # secondary-structure hypotheses only
    python -m rnafold3d import   model1.pdb model2.cif --out ens/   # package external conformers
    python -m rnafold3d view     [run_001 | models.pdb ...] [--demo]  # opens the RCV app
"""
import os as _os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    _os.environ.setdefault(_v, '1')   # one BLAS thread per process; parallelism is across trajectories

import argparse
import glob
import json
import os
import re
import sys

from .predict import PredictConfig, predict, read_fasta, read_dbn, normalize_sequence


def expand_files(files):
    """Wildcards expanded here too (the Windows shells pass model_*.pdb on unexpanded); model_2 before model_10."""
    out = []
    for f in files:
        if any(c in f for c in '*?[') and not os.path.exists(f):
            hits = glob.glob(os.path.expanduser(f))
            if not hits:
                raise FileNotFoundError(f'no files match {f}')
            out += sorted(hits, key=lambda h: [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', h)])
        else:
            out.append(f)
    return out


def _predict(a):
    if a.fasta:
        name, seq = read_fasta(a.fasta)
    else:
        name, seq = a.name, a.sequence
    seq = normalize_sequence(seq)
    dbn = read_dbn(a.dbn, len(seq)) if a.dbn else (a.structure or None)
    extra = {}
    if getattr(a, 'job_json', None):
        with open(a.job_json, encoding='utf-8') as f:
            job = json.load(f)
        extra = {k: (tuple(job[k]) if isinstance(job[k], list) else job[k])
                 for k in ('seed_structures', 'pull_restraints', 'relax_steps', 'max_structures', 'mfe_share') if k in job}
        if 'save_candidates' in job:
            extra['save_candidates'] = bool(job['save_candidates'])
    cfg = PredictConfig(n_models=a.models, trajectories=a.trajectories, seed=a.seed,
                        max_structures=extra.pop('max_structures', a.max_structures),
                        workers=a.workers, anneal_steps=a.steps, temperature=a.temperature,
                        cluster_rmsd=a.cluster_rmsd, **extra)
    m = predict(seq, a.out, name=name, dbn=dbn, config=cfg, log=None if a.quiet else print)
    print(json.dumps({'out': a.out, 'conformers': len(m['conformers']), **m['statistics']}, indent=1))
    return 0 if m['conformers'] else 2


def _validate(a):
    from .structure_io import read_structure
    from .checks import validate_model
    from .build import atom_order
    status = 0
    for path in a.files:
        chains, _ = read_structure(path)
        for cid, ch in chains.items():
            if a.chain and cid != a.chain:
                continue
            rm = ch.residue_map()
            n = len(ch.sequence)
            res = [{k: v for k, v in (rm[i + 1].atoms.items() if i + 1 in rm else []) if k in atom_order(ch.sequence[i])}
                   for i in range(n)]
            rep = validate_model(ch.sequence, res)
            print(json.dumps({'file': str(path), 'chain': cid, 'length': n, 'passed': rep['passed'],
                              'issue_counts': rep['issue_counts'], 'severe_clashes': rep['severe_clashes'],
                              'sugar_puckers': rep['sugar_puckers']}))
            status |= 0 if rep['passed'] else 1
    return status


def _ss(a):
    from .ss import predict_structures
    name, seq = read_fasta(a.fasta) if a.fasta else (a.name, a.sequence)
    hyps, meta = predict_structures(normalize_sequence(seq), a.max_structures, a.temperature)
    print(json.dumps({'name': name, 'method': meta,
                      'hypotheses': [{k: v for k, v in h.items() if k != 'pairs'} for h in hyps]}, indent=1))
    return 0


def _import(a):
    from .ensemble_io import import_ensemble
    m = import_ensemble(a.files, a.out, name=a.name, chain_id=a.chain)
    print(json.dumps({'out': a.out, 'conformers': len(m['conformers']),
                      'valid': sum(bool(c['validation'] and c['validation']['passed']) for c in m['conformers'])}, indent=1))
    return 0


def _view(a):
    """Compatibility: `rnafold3d view` opens the RCV application."""
    from rcv.app import parser, run
    argv = list(a.paths) + ['--port', str(a.port)]
    if a.no_browser:
        argv += ['--shell', 'none']
    elif a.browser_tab:
        argv += ['--shell', 'browser']
    if a.demo:
        argv.append('--demo')
    if a.runs_dir:
        argv += ['--projects-dir', a.runs_dir]
    if a.chain:
        argv += ['--chain', a.chain]
    if a.verbose:
        argv.append('--verbose')
    return run(parser().parse_args(argv))


def _export_view(a):
    from rcv.app import _export_document
    argv = list(a.paths) + ['--out', a.out] + (['--chain', a.chain] if a.chain else []) + (['--name', a.name] if a.name else [])
    return _export_document(argv)


def main(argv=None):
    ap = argparse.ArgumentParser(prog='rnafold3d', description='Physics-inspired RNA 3D ensemble prediction (v2 engine)')
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('predict', help='predict an all-atom conformer ensemble from sequence')
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument('--fasta')
    src.add_argument('--sequence')
    p.add_argument('--name', default='query')
    p.add_argument('--dbn', help='optional dot-bracket file (otherwise ViennaRNA predicts 2D)')
    p.add_argument('--structure', help='optional dot-bracket string')
    p.add_argument('--models', type=int, default=5)
    p.add_argument('--trajectories', type=int, default=0)
    p.add_argument('--max-structures', type=int, default=3)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--steps', type=int, default=0, help='annealing steps per trajectory (0 = automatic)')
    p.add_argument('--temperature', type=float, default=37.0)
    p.add_argument('--cluster-rmsd', type=float, default=0.0)
    p.add_argument('--workers', type=int, default=max(1, min(8, (os.cpu_count() or 2) - 1)))
    p.add_argument('--out', required=True)
    p.add_argument('--quiet', action='store_true')
    p.add_argument('--job-json', help='JSON with guided-exploration options (seed_structures, pull_restraints, relax_steps)')
    p.set_defaults(func=_predict)
    v = sub.add_parser('validate', help='independent all-atom validation of PDB/mmCIF files')
    v.add_argument('files', nargs='+')
    v.add_argument('--chain')
    v.set_defaults(func=_validate)
    s = sub.add_parser('ss', help='secondary-structure hypotheses from sequence (ViennaRNA)')
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument('--fasta')
    g.add_argument('--sequence')
    s.add_argument('--name', default='query')
    s.add_argument('--max-structures', type=int, default=3)
    s.add_argument('--temperature', type=float, default=37.0)
    s.set_defaults(func=_ss)
    im = sub.add_parser('import', help='package external PDB/mmCIF conformers as an ensemble (viewer contract)')
    im.add_argument('files', nargs='+')
    im.add_argument('--out', required=True)
    im.add_argument('--name', default='imported')
    im.add_argument('--chain')
    im.set_defaults(func=_import)
    vw = sub.add_parser('view', help='open the RCV application (same as the rcv command)')
    vw.add_argument('paths', nargs='*', help='ensemble package folders, PDB/mmCIF files or folders of them '
                                             '(structure files given together become one ensemble)')
    vw.add_argument('--port', type=int, default=8765)
    vw.add_argument('--no-browser', action='store_true', help='only print the URL')
    vw.add_argument('--browser-tab', action='store_true', help='open in a normal browser tab instead of an app window')
    vw.add_argument('--demo', action='store_true', help='also load the demo ensembles')
    vw.add_argument('--runs-dir', help='where predictions started from the GUI are written (default ./runs)')
    vw.add_argument('--chain', help='RNA chain to show when files contain several')
    vw.add_argument('--verbose', action='store_true')
    vw.set_defaults(func=_view)
    ev = sub.add_parser('export-view', help='write the viewer document (JSON) for an ensemble, e.g. for sessions')
    ev.add_argument('paths', nargs='+')
    ev.add_argument('--out', required=True, help="output JSON file ('-' for stdout)")
    ev.add_argument('--chain')
    ev.add_argument('--name')
    ev.set_defaults(func=_export_view)
    a = ap.parse_args(argv)
    try:
        if a.cmd in ('validate', 'import'):
            a.files = expand_files(a.files)
        return a.func(a)
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        print(f'rnafold3d {a.cmd}: error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
