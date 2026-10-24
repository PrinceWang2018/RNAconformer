"""GUI smoke test for RCV in a headless Chromium browser (Microsoft Edge or Chrome).

    python tools/viewer_smoketest.py [--out viewer_smoketest] [--browser PATH] [--only morph,video]

Starts RCV (no window) with the demo ensembles and the test harness (tools/viewer_harness), drives the
real application through the scenarios in tools/viewer_harness/steps.js, saves screenshots of the visual
scenarios and prints a pass/fail report. Works on Linux, Windows and WSL (using the Windows browser).
A temporary browser profile is used, so your own browser profile is never touched.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from rcv.desktop import find_browser, _is_wsl  # noqa: E402

SCREENSHOT = ['overview', 'morph', 'dark', 'overlay', 'welcome', 'welcomebatch', 'ss2d', 'console1', 'switchtab', 'panels', 'transitionframe',
              'tools2', 'plots', 'plotsall', 'surface1', 'plugindev', 'panelsize', 'winmax', 'langzh', 'langde', 'langja', 'guide1']
REALTIME = ['pull', 'plugins', 'exports', 'video', 'transition', 'rcvsave', 'htmlexport', 'ctxmenu', 'script1', 'tools1', 'pdffont', 'native1', 'perf',
            'welcomedemo', 'batchui', 'windows1', 'lang1',
            'predict', 'explore']   # pull needs real time
LAST = ['empty']            # closes every ensemble; must run last


def host_path(p, browser):
    """Path as seen by the browser process (Windows path for a Windows browser under WSL)."""
    if _is_wsl() and str(browser).startswith('/mnt/'):
        p = Path(p)
        parent = subprocess.check_output(['wslpath', '-w', str(p.parent)], text=True).strip()
        return parent.rstrip('\\') + '\\' + p.name
    return str(p)


class Server:
    def __init__(self, projects):
        env = dict(os.environ)
        env['PYTHONPATH'] = os.pathsep.join([str(ROOT / 'src')] + ([env['PYTHONPATH']] if env.get('PYTHONPATH') else []))
        env['RCV_DEVDIR'] = str(ROOT / 'tools' / 'viewer_harness')
        env['RCV_HOME'] = str(Path(projects).parent / 'apphome')
        env['PYTHONUNBUFFERED'] = '1'
        self.proc = subprocess.Popen([sys.executable, '-m', 'rcv', '--shell', 'none', '--demo', '--port', '0',
                                      '--projects-dir', str(projects)],
                                     cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        self.lines = []
        self.lock = threading.Lock()
        threading.Thread(target=self._pump, args=(self.proc.stderr,), daemon=True).start()
        self.url = None
        for line in self.proc.stdout:
            print('  server:', line.rstrip())
            m = re.search(r'(http://localhost:(\d+))/#token=(\S+)', line)
            if m:
                self.base, self.port, self.token = m.group(1), m.group(2), m.group(3)
                break
        if not getattr(self, 'token', None):
            raise RuntimeError('RCV did not start')
        threading.Thread(target=self._pump, args=(self.proc.stdout,), daemon=True).start()

    def _pump(self, stream):
        for line in stream:
            with self.lock:
                self.lines.append(line.rstrip())

    def logs(self, scenario):
        tag = f'[client test] [{scenario}] '
        with self.lock:
            return [l[len(tag):] for l in self.lines if l.startswith(tag)]

    def stop(self):
        self.proc.terminate()
        try:
            self.proc.wait(10)
        except subprocess.TimeoutExpired:
            self.proc.kill()


def kill_browser(proc, profile, browser):
    if _is_wsl() and str(browser).startswith('/mnt/'):
        ps = '/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe'
        name = Path(profile).name
        # stdin must not be inherited: PowerShell started from WSL can otherwise wait on it forever.
        try:
            subprocess.run([ps, '-NoProfile', '-NonInteractive', '-Command',
                            "Get-CimInstance Win32_Process -Filter \"name='msedge.exe' or name='chrome.exe'\" | "
                            f"Where-Object {{ $_.CommandLine -like '*{name}*' }} | "
                            "ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }"],
                           stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=60)
        except subprocess.TimeoutExpired:
            print(f'warning: could not stop the test browser ({name}) within 60 s')
    elif sys.platform == 'win32':
        subprocess.run(['taskkill', '/T', '/F', '/PID', str(proc.pid)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        proc.kill()
    except OSError:
        pass


def run(args):
    browser = args.browser or find_browser()
    if not browser:
        print('No Edge/Chrome browser found; pass --browser PATH.')
        return 2
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix='rnaview_smoke_', dir=out))
    print(f'browser: {browser}\noutput:  {out}')
    server = Server(work / 'projects')
    results = []
    only = set(args.only.split(',')) if args.only else None
    try:
        for name in SCREENSHOT + REALTIME + LAST:
            if only and name not in only:
                continue
            profile = work / f'profile_{name}'
            profile.mkdir()
            url = f'{server.base}/_dev/drive.html?token={server.token}&steps={name}'
            if name == 'script1':
                from urllib.parse import quote
                url += '&script=' + quote(str(ROOT / 'examples' / 'scripts' / 'tar_switch_figure.rcvs'))
            common = [browser, '--headless=new', '--hide-scrollbars', f'--window-size={args.size}', f'--force-device-scale-factor={args.scale}', '--no-first-run',
                      '--lang=en-US',            # the scenarios read English interface text (other languages: lang1, langzh, ...)
                      f'--user-data-dir={host_path(profile, browser)}']
            t0 = time.time()
            if name in SCREENSHOT or name in LAST:
                shot = out / f'{name}.png'
                cmd = common + ['--virtual-time-budget=15000', f'--screenshot={host_path(shot, browser)}', url]
                p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                try:
                    p.wait(args.timeout)
                except subprocess.TimeoutExpired:
                    pass
                kill_browser(p, profile, browser)
            else:
                p = subprocess.Popen(common + ['--remote-debugging-port=0', url], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                limit = args.timeout * (6 if name in ('predict', 'explore', 'batchui') else 1)
                while time.time() - t0 < limit:
                    logs = server.logs(name)
                    if any(l.startswith(('DONE', 'ERROR')) for l in logs):
                        break
                    time.sleep(0.5)
                kill_browser(p, profile, browser)
            logs = server.logs(name)
            ok = any(l.startswith('DONE') for l in logs) and not any(l.startswith('ERROR') for l in logs)
            detail = '; '.join(l for l in logs if not l.startswith('DONE'))[:300]
            results.append((name, ok, round(time.time() - t0, 1), detail))
            print(f"  {'PASS' if ok else 'FAIL'}  {name:9s} {results[-1][2]:6.1f}s  {detail}")
    finally:
        server.stop()
        shutil.rmtree(work, ignore_errors=True)
    failed = [r for r in results if not r[1]]
    print(f'\n{len(results) - len(failed)}/{len(results)} scenarios passed; screenshots in {out}')
    return 1 if failed else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='viewer_smoketest', help='folder for screenshots (default ./viewer_smoketest)')
    ap.add_argument('--browser', help='path to msedge/chrome executable (auto-detected by default)')
    ap.add_argument('--only', help='comma-separated scenario names')
    ap.add_argument('--timeout', type=float, default=90, help='seconds per scenario (predict/explore get 6x)')
    ap.add_argument('--size', default='1600,1000', help='browser window size W,H (default 1600,1000)')
    ap.add_argument('--scale', type=float, default=1, help='device pixel ratio of the screenshots (2 = sharp on HiDPI screens)')
    sys.exit(run(ap.parse_args()))


if __name__ == '__main__':
    main()
