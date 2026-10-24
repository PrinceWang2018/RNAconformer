"""RCV (RNA Conformation Viewer) entry point.

    rcv [PATH ...]                 open the app (start page, or the given ensembles / .rcv files)
    rcv engine predict ...         run the bundled prediction engine (rnafold3d CLI)
    rcv export-document PATH --out doc.json
    rcv convert PATH ... -o ensemble.rcv     structures / prediction folder -> one .rcv file
    rcv info FILE.rcv                        what an .rcv file contains
    rcv extract FILE.rcv -o DIR              one PDB per conformation (or --multi for one multi-model PDB)
    rcv batch SEQS.fasta [-o DIR]            predict every sequence of a FASTA file (one .rcv each)
    rcv html INPUT -o page.html              interactive single-file web page (.rcv or structures as input)

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
           'batch': _batch, 'html': _html}
    if argv and argv[0] in sub:
        return sub[argv[0]](argv[1:])
    return run(parser().parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
