"""Wall time and peak memory of RCV predictions with a given number of CPU cores (Linux).

Runs CASP16 targets with the benchmark settings (5 models, 40 trajectories, templates before 2024-05-01) and
samples the memory of the whole process tree twice a second (proportional set size, so pages shared between
the worker processes count once). Results are added to results/timing.json (used by figure.py, panel b).

    python benchmark/casp16/timing.py --workers 4 --only R1263,R1203
"""
from __future__ import annotations

import argparse
import json
import os
import threading
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / 'src'))
import compare  # noqa: E402


def _children(pid):
    out = []
    for d in os.listdir('/proc'):
        if d.isdigit():
            try:
                with open(f'/proc/{d}/stat') as f:
                    if int(f.read().rsplit(')', 1)[1].split()[1]) == pid:
                        out.append(int(d))
            except (OSError, ValueError, IndexError):
                pass
    return out


def tree_pss_mb(root):
    todo, seen, total = [root], set(), 0
    while todo:
        pid = todo.pop()
        if pid in seen:
            continue
        seen.add(pid)
        try:
            with open(f'/proc/{pid}/smaps_rollup') as f:
                for line in f:
                    if line.startswith('Pss:'):
                        total += int(line.split()[1])
                        break
        except OSError:
            continue
        todo += _children(pid)
    return total / 1024


class Peak(threading.Thread):
    def __init__(self):
        super().__init__(daemon=True)
        self.peak, self.stop = 0.0, threading.Event()

    def run(self):
        while not self.stop.is_set():
            self.peak = max(self.peak, tree_pss_mb(os.getpid()))
            self.stop.wait(0.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--only', default='R1263,R1203')
    ap.add_argument('--trajectories', type=int, default=40)
    ap.add_argument('--runs', default=str(HERE / 'runs_timing'))
    a = ap.parse_args()
    from rnafold3d.predict import predict, PredictConfig
    res_p = HERE / 'results' / 'timing.json'
    res = json.loads(res_p.read_text()) if res_p.exists() else {}
    key = 'runs' if a.workers != 16 else 'runs16'
    res['cores'] = a.workers if a.workers != 16 else res.get('cores', 4)
    res['logical_cpus'] = os.cpu_count()
    res['cpu'] = next((l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')), '')
    res.setdefault(key, [])
    for t in compare.TARGETS:
        if t['id'] not in a.only.split(','):
            continue
        out = Path(a.runs) / f"{t['id']}_w{a.workers}"
        cfg = PredictConfig(n_models=5, trajectories=a.trajectories, seed=42, workers=a.workers,
                            templates='auto', template_cutoff='2024-05-01')
        mon = Peak()
        mon.start()
        t0 = time.time()
        man = predict(t['sequence'], out, name=t['id'], config=cfg, log=None)
        wall = time.time() - t0
        mon.stop.set()
        mon.join()
        row = {'target': t['id'], 'length': len(t['sequence']), 'workers': a.workers, 'wall_time_s': round(wall, 1),
               'peak_memory_mb': round(mon.peak), 'valid': man['statistics']['valid']}
        res[key] = [r for r in res[key] if r['target'] != t['id']] + [row]
        print(row, flush=True)
    allrows = res.get('runs', []) + res.get('runs16', [])
    if allrows:
        top = max(allrows, key=lambda r: r['peak_memory_mb'])
        res['peak_memory_gb'] = round(top['peak_memory_mb'] / 1024, 2)
        res['peak_workers'] = top['workers']
    res_p.write_text(json.dumps(res, indent=1))


if __name__ == '__main__':
    main()
