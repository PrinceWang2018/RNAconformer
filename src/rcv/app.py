"""RCV (RNA Conformation Viewer) entry point.

    rcv [PATH ...]                 open the app (start page, or the given ensembles / .rcv files)
    rcv engine predict ...         run the bundled prediction engine (rnafold3d CLI)
    rcv export-document PATH --out doc.json
    rcv convert PATH ... -o ensemble.rcv     structures / prediction folder -> one .rcv file
    rcv info FILE.rcv                        what an .rcv file contains
    rcv extract FILE.rcv -o DIR              one PDB per conformation (or --multi for one multi-model PDB)
    rcv batch SEQS.fasta [-o DIR]            predict every sequence of a FASTA file (one .rcv each)
    rcv html INPUT -o page.html              interactive single-file web page (.rcv or structures as input)
    rcv pathway INPUT --a A --b B [--scan]   transition A -> B: barriers, routes, bottleneck, rigid modules, phi values
    rcv probing INPUT --data FILE [--dms]    agreement with SHAPE/DMS reactivities, maximum-entropy reweighting
    rcv confidence INPUT                     predicted per-nucleotide accuracy (confidence model)

Options: --shell auto|webview|window|browser|none, --port, --demo, --projects-dir, --chain, --verbose.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import sys
import threading
import time

from . import APP_NAME, __version__


def _engine(argv):
    from rnafold3d.cli import main as engine_main
    return engine_main(argv)


def _export_document(argv):
    ap = argparse.ArgumentParser(prog='rcv export-document')
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--out', required=True, help="output JSON file ('-' for stdout)")
    ap.add_argument('--chain')
    ap.add_argument('--name')
    a = ap.parse_args(argv)
    from .data import load_document
    doc = load_document(a.paths, chain_id=a.chain, name=a.name)
    text = json.dumps(doc, separators=(',', ':'))
    if a.out == '-':
        print(text)
    else:
        with open(a.out, 'w', encoding='utf-8') as f:
            f.write(text)
        print(json.dumps({'out': a.out, 'conformations': len(doc['conformations']), 'length': len(doc['sequence'])}))
    return 0


def _convert(argv):
    ap = argparse.ArgumentParser(prog='rcv convert', description='Write PDB/mmCIF files, an NMR multi-model file or a '
                                 'prediction folder as one .rcv file (all conformations, energies, 2D structures).')
    ap.add_argument('paths', nargs='+')
    ap.add_argument('-o', '--out', required=True, help='output .rcv file')
    ap.add_argument('--chain')
    ap.add_argument('--name')
    a = ap.parse_args(argv)
    from . import rcvfile
    from .data import collect_paths, build_document
    from .methods import methods_text
    from .server import vienna_version
    nm, confs, manifest, warnings, source = collect_paths(a.paths, a.chain)
    doc = build_document(a.name or nm, confs, manifest, warnings, source)
    out = a.out if a.out.lower().endswith('.rcv') else a.out + '.rcv'
    size = rcvfile.write(out, doc, methods=methods_text(doc, vienna_version()), provenance_manifest=manifest,
                         app_version=__version__)
    pdb = rcvfile.pdb_size_estimate(len(doc['conformations']), len(doc['atoms']['names']))
    print(json.dumps({'out': out, 'conformations': len(doc['conformations']), 'length': len(doc['sequence']),
                      'bytes': size, 'pdb_bytes_estimate': pdb, 'ratio': round(pdb / max(1, size), 1)}))
    return 0


def _info(argv):
    ap = argparse.ArgumentParser(prog='rcv info', description='Describe an .rcv file.')
    ap.add_argument('file')
    ap.add_argument('--methods', action='store_true', help='print the stored methods paragraph')
    a = ap.parse_args(argv)
    from . import rcvfile
    doc, session, info = rcvfile.read(a.file)
    confs = doc['conformations']
    print(f"{doc['name']}: {len(doc['sequence'])} nt, {len(confs)} conformations, {len(doc['atoms']['names'])} atoms")
    print(f"format {info['format']}, saved {info['saved']} by RCV {info['app_version'] or '?'}; "
          f"{info['bytes'] / 1e3:.0f} kB (as PDB about {info['pdb_bytes_estimate'] / 1e3:.0f} kB)")
    print(f"sequence  {doc['sequence']}")
    key = doc.get('default_metric')
    for c in confs[:50]:
        v = c['values'].get(key)
        print(f"  {c['id']:<16s} {c['dbn']}  {key}={v if v is not None else '-'}")
    if len(confs) > 50:
        print(f'  ... {len(confs) - 50} more')
    print('contents: ' + ', '.join(info['contents']) + (' (with viewer session)' if session else ''))
    if a.methods and info.get('methods'):
        print('\n' + info['methods'])
    return 0


def _extract(argv):
    ap = argparse.ArgumentParser(prog='rcv extract', description='Write the conformations of an .rcv file as PDB.')
    ap.add_argument('file')
    ap.add_argument('-o', '--out', required=True, help='output folder (or file with --multi)')
    ap.add_argument('--multi', action='store_true', help='one multi-model PDB file instead of one file per conformation')
    a = ap.parse_args(argv)
    from pathlib import Path
    from . import rcvfile
    doc, _, _ = rcvfile.read(a.file)
    if a.multi:
        Path(a.out).write_text(rcvfile.pdb_text(doc))
        print(a.out)
        return 0
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for k, c in enumerate(doc['conformations']):
        f = out / f"{c['id']}.pdb"
        f.write_text(rcvfile.pdb_text(doc, k))
        print(f)
    return 0


def _batch(argv):
    ap = argparse.ArgumentParser(prog='rcv batch', description='Predict every sequence of a FASTA file on this computer '
                                 '(one after another, each with all chosen cores) and save one .rcv per sequence.')
    ap.add_argument('fasta')
    ap.add_argument('-o', '--out', help='projects folder (default: the app projects folder)')
    ap.add_argument('--models', type=int, default=8)
    ap.add_argument('--sampling', choices=['fast', 'standard', 'thorough'], default='standard')
    ap.add_argument('--cores', type=int)
    ap.add_argument('--seed', type=int, default=42)
    a = ap.parse_args(argv)
    from pathlib import Path
    from . import paths as P
    from .jobs import BatchJob
    from .server import AppState
    projects = Path(a.out) if a.out else P.default_projects_dir()
    state = AppState(projects, [projects], [], False)
    job = BatchJob(state, {'fasta': Path(a.fasta).read_text(encoding='utf-8'), 'models': a.models,
                           'sampling': a.sampling, 'cores': a.cores, 'seed': a.seed})
    seen = 0
    job.start()
    while True:
        time.sleep(1.0)
        info = job.info()
        for n in info['notes'][seen:]:
            print(f"[{n['t']:7.0f}s] {n['text']}", flush=True)
        seen = len(info['notes'])
        if info['status'] in ('done', 'failed', 'cancelled'):
            break
    print(json.dumps({'status': info['status'], 'out': job.result['out'],
                      'items': [{k: it.get(k) for k in ('name', 'status', 'n_conformations', 'rcv')} for it in job.items]}))
    return 0 if info['status'] == 'done' else 1


def _html(argv):
    ap = argparse.ArgumentParser(prog='rcv html', description='Write an interactive, self-contained HTML page of an '
                                 'ensemble (for supplementary material, reviewers and collaborators).')
    ap.add_argument('paths', nargs='+', help='an .rcv file, or structures / a prediction folder')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--title')
    a = ap.parse_args(argv)
    from pathlib import Path
    from . import rcvfile
    from .htmlexport import build_html
    from .methods import methods_text
    from .server import vienna_version
    if len(a.paths) == 1 and a.paths[0].lower().endswith('.rcv'):
        doc, session, _ = rcvfile.read(a.paths[0])
        view = (session or {}).get('view')
    else:
        from .data import load_document
        doc, view = load_document(a.paths), None
    Path(a.out).write_text(build_html(doc, view, a.title, methods_text(doc, vienna_version())), encoding='utf-8')
    print(json.dumps({'out': a.out, 'bytes': Path(a.out).stat().st_size, 'conformations': len(doc['conformations'])}))
    return 0


def _load_any(paths):
    """Document (and saved view) of an .rcv file or of structures / a prediction folder."""
    from . import rcvfile
    if len(paths) == 1 and paths[0].lower().endswith('.rcv'):
        doc, session, _ = rcvfile.read(paths[0])
        return doc, (session or {}).get('view')
    from .data import load_document
    return load_document(paths), None


def _conf_index(doc, key):
    """Conformation by id (conf_003), label or 1-based number."""
    for k, c in enumerate(doc['conformations']):
        if key in (c.get('id'), c.get('label')):
            return k
    if str(key).isdigit() and 1 <= int(key) <= len(doc['conformations']):
        return int(key) - 1
    raise SystemExit(f'unknown conformation {key!r} (use an id such as conf_001 or a number 1..{len(doc["conformations"])})')


def _conditions_args(ap):
    ap.add_argument('--temperature', type=float, default=37.0, help='°C (default 37)')
    ap.add_argument('--sodium', type=float, default=1.021, help='Na+ in mol/L (default 1.021, Turner standard)')
    ap.add_argument('--ligand', help='aptamer preset (theophylline) or MOTIF_SEQ:MOTIF_STRUCT')
    ap.add_argument('--kd', type=float, help='ligand Kd in µM (presets have their own)')
    ap.add_argument('--conc', type=float, default=0.0, help='ligand concentration in µM')


def _conditions(a):
    from .conditions import Conditions
    lig = None
    if a.ligand:
        if ':' in a.ligand:
            seq, struct = a.ligand.split(':', 1)
            lig = {'name': 'custom aptamer', 'motif_seq': seq, 'motif_struct': struct, 'kd_uM': a.kd or 1.0, 'conc_uM': a.conc}
        else:
            lig = {'preset': a.ligand, 'conc_uM': a.conc, **({'kd_uM': a.kd} if a.kd else {})}
    return Conditions.from_json({'temperature': a.temperature, 'sodium': a.sodium, 'ligand': lig})


def _pathway(argv):
    ap = argparse.ArgumentParser(prog='rcv pathway', description='Transition analysis from conformation A to B: '
                                 'refolding paths and barriers (ViennaRNA findpath), routes through intermediates with '
                                 'transition-path-theory flux, the cooperative unit of base pairs, rigid modules and, '
                                 'with --scan, the in-silico phi-value scan of all point mutations.')
    ap.add_argument('paths', nargs='+', help='an .rcv file, or structures / a prediction folder')
    ap.add_argument('--a', required=True, help='conformation A (id such as conf_001, or number)')
    ap.add_argument('--b', required=True, help='conformation B')
    _conditions_args(ap)
    ap.add_argument('--k0', type=float, default=1e6, help='Arrhenius prefactor in 1/s (default 1e6)')
    ap.add_argument('--no-intermediates', action='store_true', help='direct path only')
    ap.add_argument('--scan', action='store_true', help='also the phi-value scan (3 x length mutants)')
    ap.add_argument('-o', '--out', help='write the full result as JSON')
    a = ap.parse_args(argv)
    from .rcvfile import doc_coords
    from .transition import analyse, phi_scan
    doc, _ = _load_any(a.paths)
    ia, ib = _conf_index(doc, a.a), _conf_index(doc, a.b)
    cond = _conditions(a)
    r = analyse(doc, ia, ib, cond, k0=a.k0, intermediates=not a.no_intermediates, coords=doc_coords(doc))
    if a.scan and r['A']['dbn'] != r['B']['dbn'] and not r['same_basin']:
        r['phi_scan'] = phi_scan(doc['sequence'], [r['nodes'][k]['dbn'] for k in r['best']['route']], cond)
    ids = [c.get('id') for c in doc['conformations']]
    lines = [f"{doc.get('name')}: {ids[ia]} -> {ids[ib]} ({cond.describe()})",
             f"  A {r['A']['dbn']}  {r['A']['energy']:.2f} kcal/mol", f"  B {r['B']['dbn']}  {r['B']['energy']:.2f} kcal/mol"]
    lines.append(f"  basins: A {r['A']['above_basin']:+.2f} and B {r['B']['above_basin']:+.2f} kcal/mol above their minima; dG between basins {r['dG_basins']:+.2f} kcal/mol")
    if r['A']['dbn'] == r['B']['dbn']:
        lines.append('  same secondary structure: tertiary rearrangement without a 2D barrier')
    elif r['same_basin']:
        lines.append('  A and B relax into the same 2D basin: no barrier between them')
    else:
        lines.append(f"  barrier: direct {r['direct']['barrier']:.2f}, best route {r['best']['route_barrier']:.2f}, effective {r['effective_barrier']:.2f} kcal/mol;"
                     f" dG(A->B) {r['dG']:+.2f} kcal/mol; half-life of A ~ {r['half_life_s']:.3g} s (k0 = {a.k0:.0e}/s)")
        for i, p in enumerate((r['tpt'] or {}).get('pathways', [])[:5]):
            lines.append(f"  route {i + 1}: {' -> '.join(str(k) for k in p['nodes'])}  flux {p['flux_share']:.1%}")
        b = r['bottleneck']
        if b:
            lines.append(f"  cooperative unit: {len(b['unit'])} bp {[(i + 1, j + 1) for i, j in b['unit']]}; nucleus of B: {len(b['nucleus'])} bp")
        if r.get('trap'):
            where = {'a-side': 'on the A side', 'deep': 'on the way, deeper than B'}.get(r['trap']['kind'], 'on the way to B')
            share = f"; {r['trap']['time_share']:.0%} of the switching time is spent there" if r['trap'].get('time_share') is not None else ''
            lines.append(f"  trap: node {r['trap']['node']} lies {r['trap']['depth']:.2f} kcal/mol below the basin of A, {where}{share}")
    mo = r.get('motion') or {}
    for m in mo.get('modules', []):
        lines.append(f"  module M{m['id'] + 1}: {m['size']} nt, turn {m['rotation']:.0f} deg, shift {m['translation']:.1f} A")
    if r.get('phi_scan'):
        sc = r['phi_scan']
        def mutation(x):
            return f"{x['from']}{x['pos'] + 1}{x['to']} ({x['ddg_ts']:+.1f})"
        lines.append('  fastest: ' + ', '.join(mutation(x) for x in sc['faster'][:5]))
        lines.append('  slowest: ' + ', '.join(mutation(x) for x in sc['slower'][:5]))
    print('\n'.join(lines))
    if a.out:
        with open(a.out, 'w', encoding='utf-8') as f:
            json.dump(r, f, indent=1)
    return 0


def _probing(argv):
    ap = argparse.ArgumentParser(prog='rcv probing', description='Agreement of every conformation with a SHAPE or DMS '
                                 'reactivity profile and maximum-entropy reweighting of the ensemble to fit it.')
    ap.add_argument('paths', nargs='+', help='an .rcv file, or structures / a prediction folder')
    ap.add_argument('--data', required=True, help='reactivity file (position and value, .shape, .map, CSV or one value per line)')
    ap.add_argument('--dms', action='store_true', help='DMS data (A and C only); default SHAPE')
    ap.add_argument('--raw', action='store_true', help='do not normalise (values are already normalised)')
    ap.add_argument('--theta', type=float, help='reweighting strength (default 0.1 x number of data points)')
    ap.add_argument('-o', '--out', help='write the result as JSON')
    a = ap.parse_args(argv)
    from pathlib import Path
    from . import probing
    doc, _ = _load_any(a.paths)
    kind = 'DMS' if a.dms else 'SHAPE'
    vals, fmt, notes = probing.parse(Path(a.data).read_text(encoding='utf-8', errors='replace'), len(doc['sequence']))
    if not a.raw:
        vals, _ = probing.normalise(vals, doc['sequence'], kind)
    agr = probing.agreement(doc, vals, kind)
    rw = probing.reweight(doc, vals, kind, theta=a.theta) if len(doc['conformations']) > 1 else None
    ids = [c.get('id') for c in doc['conformations']]
    print(f"{doc.get('name')}: {kind}, {sum(v is not None for v in vals)} of {len(vals)} nucleotides ({fmt}{'; ' + '; '.join(notes) if notes else ''})")
    for k, cid in enumerate(ids):
        print(f"  {cid:12s} AUC {agr[k] if agr[k] is not None else float('nan'):.3f}" + (f"   weight {rw['weights'][k]:.3f}" if rw else ''))
    if rw:
        print(f"  chi2/n {rw['chi2_prior']:.3f} -> {rw['chi2']:.3f}; effective conformations {rw['n_eff_conformations']:.1f} of {len(ids)}")
    if a.out:
        with open(a.out, 'w', encoding='utf-8') as f:
            json.dump({'kind': kind, 'values': vals, 'agreement': dict(zip(ids, agr)), 'reweighting': rw}, f, indent=1)
    return 0


def _confidence(argv):
    ap = argparse.ArgumentParser(prog='rcv confidence', description='Predicted per-nucleotide accuracy (lDDT) of every '
                                 'conformation from the RCV confidence model (a small neural network; coordinates only).')
    ap.add_argument('paths', nargs='+', help='an .rcv file, or structures / a prediction folder')
    ap.add_argument('-o', '--out', help='write per-nucleotide values as JSON')
    a = ap.parse_args(argv)
    from rnafold3d import qa
    from .server import confidence_of
    if not qa.available():
        raise SystemExit('the confidence model is not installed (rnafold3d/data/qa_v1.npz)')
    doc, _ = _load_any(a.paths)
    r = confidence_of(doc)
    for c, v in zip(doc['conformations'], r['conformations']):
        print(f"  {c.get('id'):12s} {v['global']:.3f}")
    if a.out:
        with open(a.out, 'w', encoding='utf-8') as f:
            json.dump({c.get('id'): v for c, v in zip(doc['conformations'], r['conformations'])}, f, indent=1)
    return 0


def parser():
    ap = argparse.ArgumentParser(prog='rcv', description=f'RCV ({APP_NAME}): predict, explore and present RNA '
                                 'conformations on this computer (prediction engine + interactive viewer).',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='tools (no window):\n'
                                        '  rcv convert PATH... -o x.rcv     structures / prediction folder -> one .rcv file\n'
                                        '  rcv info x.rcv [--methods]       what an .rcv file contains\n'
                                        '  rcv extract x.rcv -o DIR         PDB files again (--multi: one multi-model file)\n'
                                        '  rcv batch seqs.fasta [-o DIR]    predict every sequence of a FASTA file\n'
                                        '  rcv html x.rcv -o page.html      interactive single-file web page\n'
                                        '  rcv pathway x.rcv --a 1 --b 5    barriers, routes, bottleneck, modules (--scan: phi values)\n'
                                        '  rcv probing x.rcv --data r.shape agreement with SHAPE/DMS data, reweighting\n'
                                        '  rcv confidence x.rcv             predicted per-nucleotide accuracy\n'
                                        '  rcv engine predict|ss|validate   the rnafold3d prediction engine\n'
                                        '  rcv export-document PATH --out doc.json')
    ap.add_argument('paths', nargs='*', help='rnafold3d run folders, PDB/mmCIF files (multi-model files give one '
                                             'conformation per model) or folders of them; each argument opens as one ensemble')
    ap.add_argument('--shell', choices=['auto', 'webview', 'window', 'browser', 'none'], default='auto',
                    help='window type (default: native window if pywebview is installed, else an Edge/Chrome app window)')
    ap.add_argument('--port', type=int, default=8765, help='local port (default 8765; the next free one is used if busy)')
    ap.add_argument('--demo', action='store_true', help='open the demo ensembles')
    ap.add_argument('--projects-dir', help='folder for predictions (default: Documents/RCV/projects)')
    ap.add_argument('--chain', help='RNA chain to show when files contain several')
    ap.add_argument('--keep-running', action='store_true', help='do not quit when the window is closed')
    ap.add_argument('--verbose', action='store_true')
    ap.add_argument('--version', action='version', version=f'{APP_NAME} {__version__}')
    return ap


def run(a, ready=None):
    from .desktop import open_external, run_webview, webview_available
    from .server import make_server
    from . import paths as P
    server, state = make_server(a.port, a.verbose, a.projects_dir)
    for p in a.paths:
        try:
            s = state.open_paths([p], a.chain)
            print(f'opened {s["name"]}: {s["n_conformations"]} conformations, {s["length"]} nt', flush=True)
        except Exception as exc:
            print(f'could not open {p}: {exc}', file=sys.stderr, flush=True)
    if a.demo and P.demo_dir():
        for name in ('HIV1_TAR_NMR_1ANR.pdb', 'HIV1_TAR_prediction', '1EHZ_tRNA_prediction'):
            if (P.demo_dir() / name).exists():
                s = state.open_paths([P.demo_dir() / name])
                print(f'opened {s["name"]}: {s["n_conformations"]} conformations, {s["length"]} nt', flush=True)
    threading.Thread(target=server.serve_forever, kwargs={'poll_interval': 0.3}, daemon=True).start()
    url = f'http://localhost:{server.server_address[1]}/#token={state.token}'
    print(f'{APP_NAME} {__version__} running at {url}', flush=True)
    print(f'projects folder: {state.projects_dir}', flush=True)
    if ready:
        ready(server, state, url)

    def shutdown():
        for job in state.jobs.values():
            job.cancel()
        server.shutdown()
        server.server_close()

    try:
        if a.shell in ('auto', 'webview') and webview_available():
            run_webview(url, state, on_closed=lambda: None)        # blocks until the window closes
            return 0
        if a.shell == 'webview':
            print('pywebview is not available here; using an app window instead', flush=True)
        used = 'none' if a.shell == 'none' else open_external(url, 'browser' if a.shell == 'browser' else 'window')
        state.shell = used
        if used != 'none':
            print('opened in an app window' if used == 'window' else 'opened in the default browser', flush=True)
            print('Close the window (or press Ctrl+C here) to quit.', flush=True)
        while True:
            time.sleep(0.5)
            if (used != 'none' and not a.keep_running and state.closing_at and time.time() > state.closing_at):
                print('window closed; quitting', flush=True)
                return 0
    except KeyboardInterrupt:
        print(f'\nstopping {APP_NAME}', flush=True)
        return 0
    finally:
        shutdown()


def _ensure_streams():
    """Windowed installer builds have no console: send output to a log file in the user data folder."""
    if sys.stdout is not None and sys.stderr is not None:
        return
    from . import paths as P
    try:
        log_dir = P.user_data_dir() / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)
        stream = open(log_dir / 'rcv.log', 'a', encoding='utf-8', buffering=1)
    except OSError:
        import io
        stream = io.StringIO()
    if sys.stdout is None:
        sys.stdout = stream
    if sys.stderr is None:
        sys.stderr = stream


def main(argv=None):
    multiprocessing.freeze_support()              # engine worker processes in frozen builds
    _ensure_streams()
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == 'engine':
        return _engine(argv[1:])
    sub = {'export-document': _export_document, 'convert': _convert, 'info': _info, 'extract': _extract,
           'batch': _batch, 'html': _html, 'pathway': _pathway, 'probing': _probing, 'confidence': _confidence}
    if argv and argv[0] in sub:
        return sub[argv[0]](argv[1:])
    return run(parser().parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
