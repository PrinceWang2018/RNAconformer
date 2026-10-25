"""Choose the decoy targets: drop tRNAs hidden in ribosome entries (benchmark family tRNA-Phe), then draw a
length-stratified subset that fits the compute budget (fixed seed). Rewrites <data>/dataset_v1.json (the full list
stays in dataset_all.json).

    python tools/qa_model/select_targets.py --data benchmark/data_qa
"""
import argparse
import json
import random
from pathlib import Path

RIBOSOME_WORDS = ('ribosom', 'disome', 'polysome', 'rnc', 'translat', 'stalled', 'nascent', 'decoding', 'frameshift complex',
                  '70s', '80s', '50s', '30s', '40s', '60s', 'elongation')
QUOTA = ((15, 40, 60), (41, 80, 50), (81, 120, 22))      # (min length, max length, targets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default='benchmark/data_qa')
    ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args()
    data = Path(a.data)
    full = data / 'dataset_all.json'
    if not full.exists():
        full.write_text((data / 'dataset_v1.json').read_text())
    d = json.loads(full.read_text())
    keep, dropped = [], []
    for s in d['samples']:
        seq = (data / 'samples' / s['sample_id'] / 'sequence.fasta').read_text().split('\n')[1].strip()
        title = s['title'].lower()
        trna_like = 70 <= s['length'] <= 95 and (seq.endswith('CCA') or any(w in title for w in RIBOSOME_WORDS))
        (dropped if trna_like else keep).append(s)
    rng = random.Random(a.seed)
    chosen = []
    for lo, hi, k in QUOTA:
        pool = [s for s in keep if lo <= s['length'] <= hi]
        rng.shuffle(pool)
        chosen += pool[:k]
    d2 = {**d, 'samples': sorted(chosen, key=lambda s: s['length']), 'selection': {
        'seed': a.seed, 'quota': QUOTA, 'dropped_trna_like': [s['sample_id'] for s in dropped],
        'eligible': len(keep)}}
    (data / 'dataset_v1.json').write_text(json.dumps(d2, indent=1))
    print(f"{len(dropped)} tRNA-like chains dropped, {len(keep)} eligible, {len(chosen)} chosen")


if __name__ == '__main__':
    main()
