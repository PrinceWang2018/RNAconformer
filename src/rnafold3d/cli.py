"""Command-line interface for the v2 engine.

    python -m rnafold3d predict  --fasta x.fa [--dbn x.dbn] --models 5 --seed 42 --out run_001 [--workers 8]
    python -m rnafold3d validate run_001/conformations/conf_001.pdb [more.pdb ...]
    python -m rnafold3d ss       --fasta x.fa            # secondary-structure hypotheses only
    python -m rnafold3d import   model1.pdb model2.cif --out ens/   # package external conformers
    python -m rnafold3d view     [run_001 | models.pdb ...] [--demo]  # opens the RCV app
    python -m rnafold3d pockets  structure.pdb               # where a compound could bind
    python -m rnafold3d dock     --receptor rna.pdb --ligand compound.sdf --pocket 1 --out poses.sdf
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
                 for k in ('seed_structures', 'pull_restraints', 'relax_steps', 'max_structures', 'mfe_share',
                           'reactivities', 'alignment') if k in job}
        if 'save_candidates' in job:
            extra['save_candidates'] = bool(job['save_candidates'])
    if a.shape:
        from .ss import read_reactivities
        extra['reactivities'] = tuple(read_reactivities(a.shape, len(seq)))
    if a.alignment:
        from .ss import read_alignment
        extra['alignment'] = tuple(read_alignment(a.alignment, seq))
    for key in ('reactivities', 'alignment'):
        if key in extra:
            extra[key] = tuple(extra[key])
    cfg = PredictConfig(n_models=a.models, trajectories=a.trajectories, seed=a.seed,
                        max_structures=extra.pop('max_structures', a.max_structures),
                        workers=a.workers, anneal_steps=a.steps, temperature=a.temperature,
                        cluster_rmsd=a.cluster_rmsd, templates=a.templates, template_cutoff=a.template_cutoff,
                        template_exclude=tuple(x for x in a.template_exclude.replace(' ', '').split(',') if x),
                        pseudoknots=a.pseudoknots, fragments=a.fragments, domains=a.domains, **extra)
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


def _pockets(a):
    """Pockets of a structure: where a compound could bind (rnafold3d.screening.pockets)."""
    import numpy as np

    from rnafold3d.screening import pockets as P
    from rnafold3d.screening.receptor import Receptor
    rec = Receptor.from_pdb(a.structure, chain=a.chain)
    found = P.find(rec.xyz, elements=rec.element, residue_of=rec.residue, burial_min=a.burial_min,
                   keep_points=bool(a.out_map), top=a.top)
    if not found:
        print('no pocket of at least %.0f A^3 found' % P.MIN_VOLUME)
        return 0
    names = {int(r): f'{rec.resname[k]}{rec.name[k] and ""}' for k, r in enumerate(rec.residue)}
    print(f'{len(found)} pockets in {a.structure}')
    for p in found:
        res = ','.join(str(r + 1) for r in p['residues'][:24]) + ('…' if len(p['residues']) > 24 else '')
        print(f"  {p['id']:2d}  {p['volume']:7.0f} A^3  enclosed {100 * p['burial']:3.0f}%  score {p['score']:.2f}"
              f"  centre {p['center'][0]:7.1f}{p['center'][1]:8.1f}{p['center'][2]:8.1f}"
              f"{'  (a groove rather than a pocket)' if p['wide'] else ''}\n      lining residues: {res}")
    if a.json:
        import json as _json
        with open(a.json, 'w', encoding='utf-8') as fh:
            _json.dump([{k: v for k, v in p.items() if k != 'points'} for p in found], fh, indent=1)
        print(f'written to {a.json}')
    if a.out_map:
        with open(a.out_map, 'w', encoding='utf-8') as fh:     # the pocket points as a PDB of pseudo-atoms
            for p in found:
                for n, (x, y, z) in enumerate(np.asarray(p['points'])[::4]):
                    fh.write(f'HETATM{n % 99999 + 1:5d}  C   PKT P{p["id"]:3d}    '
                             f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00{p["score"]:6.2f}           C\n')
        print(f'pocket points written to {a.out_map}')
    return 0


def _dock(a):
    """Dock one or more compounds into a pocket of an RNA (rnafold3d.screening)."""
    import numpy as np

    from rnafold3d.screening import chem, dock as D, grid as G, pockets as P
    from rnafold3d.screening.receptor import Receptor, read_complex
    if a.from_complex:
        rec, found_ligands = read_complex(a.receptor, chain=a.chain)
    else:
        rec, found_ligands = Receptor.from_pdb(a.receptor, chain=a.chain), []
    compounds = chem.read(a.ligand) if a.ligand else found_ligands
    compounds = [m for m in compounds if len(m) and 'error' not in m.properties]
    if not compounds:
        raise ValueError('no compound could be read')
    if a.centre:
        box = {'centre': [float(x) for x in a.centre.split(',')], 'size': a.size}
        where = f'box of {a.size:.0f} A at {a.centre}'
    elif a.site_of:
        ref = [m for m in found_ligands if m.properties.get('residue') == a.site_of]
        if not ref:
            raise ValueError(f'{a.site_of} is not a compound in this file')
        span = float((ref[0].xyz.max(0) - ref[0].xyz.min(0)).max())
        box = {'centre': ref[0].xyz.mean(0), 'size': span + 2 * a.pad}
        where = f'the site of {a.site_of}'
    else:
        found = P.find(rec.xyz, elements=rec.element, residue_of=rec.residue)
        if not found:
            raise ValueError('no pocket found: give --centre x,y,z or --site-of RES')
        if not 1 <= a.pocket <= len(found):
            raise ValueError(f'this structure has {len(found)} pockets, not {a.pocket}')
        box = {'pocket': found[a.pocket - 1]}
        where = (f'pocket {a.pocket} of {len(found)} ({found[a.pocket - 1]["volume"]:.0f} A^3, score '
                 f'{found[a.pocket - 1]["score"]:.2f})')
    grid = G.Grid.for_pocket(rec, pad=a.pad, **box)
    print(f'receptor {len(rec)} atoms, {len(rec.ring_centre)} bases; docking into {where}; {grid!r}')
    rows, out_mols, out_props = [], [], []
    for m in compounds[:a.limit]:
        poses = D.dock(m, grid, receptor=rec, exhaustiveness=a.exhaustiveness, seed=a.seed, n_poses=a.poses)
        if not poses:
            continue
        best = poses[0]
        rows.append((m, poses))
        print(f'  {m.name or "compound":24s} {len(m):3d} atoms {m.n_rotatable} rotatable  '
              f'score {best.score:7.2f}  ({len(poses)} poses kept)')
        for n, p in enumerate(poses):
            out_mols.append(m)
            out_props.append({'pose': n + 1, 'score': round(p.score, 3), 'compound': m.name,
                              'grid_score': round(p.grid_score, 3)})
    if a.out and rows:
        chem.write_sdf([m for m in out_mols], a.out, xyz=[p.xyz for _m, ps in rows for p in ps],
                       properties=out_props)
        print(f'{len(out_props)} poses written to {a.out}')
    print('Scores rank poses and compounds; they are not affinities and not kcal/mol.')
    return 0


def _screen(a):
    """Screen a library against one or two conformations of an RNA (rnafold3d.screening.screen)."""
    import json

    import numpy as np

    from rnafold3d.screening import chem, pockets as P, screen as SC
    from rnafold3d.screening.receptor import Receptor, read_complex

    def receptor_of(path, chain=None):
        rec, ligands = read_complex(path, chain=chain)
        return rec, ligands

    def site_of(rec, ligands, which):
        if a.centre:
            return SC.Site(centre=tuple(float(x) for x in a.centre.split(',')), size=a.size,
                           label=f'a box at {a.centre}')
        if a.site_of:
            ref = [m for m in ligands if a.site_of in m.properties.get('residues', [])]
            if not ref:
                raise ValueError(f'{a.site_of} is not a compound in the {which} structure')
            return SC.Site.from_atoms(ref[0].xyz, pad=a.pad, label=f'the site of {a.site_of}')
        found = P.find(rec.xyz, elements=rec.element, residue_of=rec.residue, burial_min=a.burial_min)
        if not found:
            raise ValueError('no pocket found: give --centre x,y,z or --site-of RES')
        if not 1 <= a.pocket <= len(found):
            raise ValueError(f'the {which} structure has {len(found)} pockets, not {a.pocket}')
        p = found[a.pocket - 1]
        print(f'  {which}: pocket {a.pocket} of {len(found)}, {p["volume"]:.0f} A^3, '
              f'{len(p["residues"])} nucleotides lining it')
        return SC.Site.from_pocket(p, pad=a.pad)

    rec_a, lig_a = receptor_of(a.state_a, a.chain)
    site_a = site_of(rec_a, lig_a, 'state A')
    rec_b = site_b = None
    if a.state_b:
        rec_b, lig_b = receptor_of(a.state_b, a.chain)
        if not a.same_site:
            site_b = site_of(rec_b, lig_b, 'state B')
        elif site_a.residues:               # the same nucleotides, wherever they have moved to in B
            res = np.isin(rec_b.residue, np.asarray(site_a.residues, int))
            site_b = SC.Site.from_atoms(rec_b.xyz[res], pad=a.pad, label='the same nucleotides',
                                        residues=site_a.residues)
        else:                               # a box, not a pocket: the same box, since A and B are aligned
            site_b = SC.Site(centre=site_a.centre, size=site_a.size, label='the same box')
    compounds = [m for m in chem.read(a.library) if len(m) and 'error' not in (m.properties or {})]
    if not compounds:
        raise ValueError('no compound could be read from the library')
    compounds = compounds[:a.limit]
    screen = SC.Screen(rec_a, site_a, rec_b, site_b, exhaustiveness=a.exhaustiveness, seed=a.seed,
                       n_poses=a.poses)
    jobs = SC.workers(a.jobs, per_worker_mb=120)
    screen.prepare()
    print(f'{len(compounds)} compounds, {"two states" if rec_b is not None else "one state"}, '
          f'{a.exhaustiveness} search, {jobs} worker{"s" if jobs > 1 else ""}, '
          f'{screen.nbytes / 1e6:.0f} MB of grids')
    done = [0]

    def progress(n, total, row):
        done[0] = n
        tag = '' if row.get('selectivity') is None else f'  selectivity {row["selectivity"]:+6.2f}'
        print(f'  [{n:5d}/{total}] {str(row.get("name"))[:28]:28s} '
              f'A {_fmt(row.get("score_a"))}  B {_fmt(row.get("score_b"))}{tag}'
              f'{"  " + row["error"] if row.get("error") else ""}')

    rows = SC.run(screen, compounds, jobs=jobs, keep_poses=a.keep_poses, progress=progress)
    s = SC.summarise(rows, top=a.top)
    print(f'\n{s["compounds"]} compounds, {s["scored"]} scored, {s["failed"]} could not be placed, '
          f'{s["seconds"] / 60:.1f} min of docking')
    if 'selectivity_median' in s:
        print(f'  selectivity: median {s["selectivity_median"]:+.2f}, '
              f'from {s["selectivity_min"]:+.2f} to {s["selectivity_max"]:+.2f}')
    print(f'  the {len(s["top"])} best:')
    for i, r in enumerate(s['top']):
        print(f'   {i + 1:3d}. {str(r["name"])[:30]:30s} A {_fmt(r["score_a"])}  B {_fmt(r["score_b"])}  '
              f'selectivity {_fmt(r["selectivity"])}  strain {r["strain_a"]:.2f}')
    if a.out:
        from pathlib import Path
        out = Path(a.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({'states': {'a': site_a.as_dict(),
                                              'b': None if site_b is None else site_b.as_dict()},
                                   'parameters': {'exhaustiveness': a.exhaustiveness, 'seed': a.seed,
                                                  'poses': a.poses, 'library': a.library},
                                   'summary': {k: v for k, v in s.items() if k != 'top'},
                                   'rows': sorted(rows, key=SC.sort_key)}, indent=1), encoding='utf-8')
        print(f'written to {out}')
    print('\nSelectivity ranks compounds; it is not an affinity, not kcal/mol, and says nothing about\n'
          'whether a compound will block the transition. See docs/SCREENING_STATUS.md for what is measured.')
    return 0


def _fmt(x):
    return '    -' if x is None else f'{x:6.2f}'


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
    p.add_argument('--templates', choices=('auto', 'off'), default='auto',
                   help='use related experimental structures from the built-in library as templates (default auto)')
    p.add_argument('--template-cutoff', default='', metavar='YYYY-MM-DD',
                   help='only use templates released before this date (blind benchmarks)')
    p.add_argument('--template-exclude', default='', metavar='PDB,PDB', help='PDB ids never used as templates')
    p.add_argument('--pseudoknots', choices=('auto', 'off'), default='off',
                   help='also try a 2D with one crossing helix when it is more stable than the nested MFE (default off)')
    p.add_argument('--fragments', choices=('auto', 'off'), default='auto',
                   help='graft hairpin-loop and junction fragments of library structures with a matching sequence (default auto)')
    p.add_argument('--domains', choices=('auto', 'off'), default='off',
                   help='experimental: chains of 240 nt or more fold their domains separately, then assemble them (default off)')
    p.add_argument('--shape', metavar='FILE', help='chemical-probing reactivities (SHAPE, DMS ...): one value per line or '
                   '"position value" lines; they guide the 2D hypotheses as pseudo-energies')
    p.add_argument('--alignment', metavar='FILE', help='aligned homologous sequences (FASTA with gaps, Stockholm or Clustal) '
                   'that contain the query: the covariation-supported consensus 2D becomes the first hypothesis')
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
    pk = sub.add_parser('pockets', help='pockets of a structure: where a compound could bind')
    pk.add_argument('structure', help='PDB file of the RNA')
    pk.add_argument('--chain')
    pk.add_argument('--top', type=int, default=10)
    pk.add_argument('--burial-min', type=int, default=5, choices=(4, 5, 6, 7),
                    help='of seven directions, how enclosed a site must be (5 = the sites known RNA '
                         'binders use, 6 = only deep pockets; see benchmark/screening/pocket_sites.py)')
    pk.add_argument('--json', help='write the pockets as JSON')
    pk.add_argument('--out-map', help='write the pocket points as a PDB of pseudo-atoms, to look at in RCV')
    pk.set_defaults(func=_pockets)
    dk = sub.add_parser('dock', help='dock compounds into a pocket (experimental: no affinity is predicted)')
    dk.add_argument('--receptor', required=True, help='PDB file of the RNA')
    dk.add_argument('--ligand', help='compounds to dock (SDF, MOL, MOL2, PDB); without it, the compounds in '
                                     'the receptor file are redocked')
    dk.add_argument('--chain')
    dk.add_argument('--pocket', type=int, default=1, help='which pocket of the structure (default the best)')
    dk.add_argument('--centre', help='x,y,z of the docking box instead of a pocket')
    dk.add_argument('--size', type=float, default=22.0, help='edge of the box when --centre is given')
    dk.add_argument('--site-of', metavar='RES', help='dock into the site of a compound already in the file '
                                                     '(a redocking test)')
    dk.add_argument('--pad', type=float, default=6.0, help='Å around the pocket or reference compound')
    dk.add_argument('--exhaustiveness', default='standard', choices=('quick', 'standard', 'thorough'))
    dk.add_argument('--seed', type=int, default=1)
    dk.add_argument('--poses', type=int, default=5)
    dk.add_argument('--limit', type=int, default=1000, help='at most this many compounds from the file')
    dk.add_argument('--from-complex', action='store_true', default=True,
                    help='read compounds in the receptor file as compounds, not as part of the RNA')
    dk.add_argument('--out', help='write the poses as an SDF')
    dk.set_defaults(func=_dock)
    sn = sub.add_parser('screen', help='screen a compound library against one or two conformations '
                                       '(experimental: a ranking, not an affinity)')
    sn.add_argument('--state-a', required=True, help='PDB file of conformation A')
    sn.add_argument('--state-b', help='PDB file of conformation B; without it, a one-state ranking')
    sn.add_argument('--library', required=True, help='compounds to screen (SDF, MOL, MOL2)')
    sn.add_argument('--chain')
    sn.add_argument('--pocket', type=int, default=1, help='which pocket of the structure (default the best)')
    sn.add_argument('--centre', help='x,y,z of the docking box instead of a pocket')
    sn.add_argument('--size', type=float, default=22.0, help='edge of the box when --centre is given')
    sn.add_argument('--site-of', metavar='RES', help='the site of a compound already in the file')
    sn.add_argument('--pad', type=float, default=6.0, help='Å around the pocket')
    sn.add_argument('--same-site', action='store_true', default=True,
                    help='measure state B at the nucleotides that line the site in state A (the default: '
                         'what a two-state screen asks is what happens to this site, and a pocket that has '
                         'closed in B has no pocket to match)')
    sn.add_argument('--own-pocket', dest='same_site', action='store_false',
                    help='instead, use the pocket ranking of state B in its own right')
    sn.add_argument('--burial-min', type=int, default=5, choices=(4, 5, 6, 7))
    sn.add_argument('--exhaustiveness', default='quick', choices=('quick', 'standard', 'thorough'))
    sn.add_argument('--seed', type=int, default=1)
    sn.add_argument('--poses', type=int, default=3)
    sn.add_argument('--keep-poses', type=int, default=0, help='keep the coordinates of this many compounds')
    sn.add_argument('--jobs', type=int, default=0, help='worker processes (0 = every core)')
    sn.add_argument('--limit', type=int, default=200000)
    sn.add_argument('--top', type=int, default=20)
    sn.add_argument('--out', help='write the ranked table as JSON')
    sn.set_defaults(func=_screen)
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
