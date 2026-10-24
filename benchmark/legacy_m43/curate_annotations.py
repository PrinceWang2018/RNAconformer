"""RNAFold3D-Bench M4.2: conservative, auditable secondary-structure curation.

No 2D labels are inferred from 3D coordinates here. Input annotation provenance is mandatory.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

OPEN = '([{<'
CLOSE = ')]}>'
PAIR = dict(zip(OPEN, CLOSE))
REVERSE = {v:k for k,v in PAIR.items()}
CANON = {frozenset(x) for x in ('AU','GC','GU')}


def parse_dotbracket(sequence, structure):
    sequence = sequence.upper().replace('T','U')
    if not sequence or any(x not in 'ACGU' for x in sequence):
        raise ValueError('Sequence must contain only A/C/G/U; resolve modifications explicitly')
    if len(sequence) != len(structure):
        raise ValueError('Sequence/dot-bracket length mismatch')
    stacks = {x:[] for x in OPEN}
    pairs=[]
    for i,ch in enumerate(structure, 1):
        if ch in OPEN: stacks[ch].append(i)
        elif ch in REVERSE:
            opening = REVERSE[ch]
            if not stacks[opening]: raise ValueError(f'Unmatched closing bracket at {i}')
            j=stacks[opening].pop()
            pairs.append((j,i))
        elif ch != '.': raise ValueError(f'Unsupported dot-bracket symbol {ch!r}')
    if any(stacks[x] for x in OPEN): raise ValueError('Unclosed bracket')
    pairs.sort()
    crossing = any(i<k<j<l or k<i<l<j for n,(i,j) in enumerate(pairs) for k,l in pairs[n+1:])
    noncanonical=[(i,j) for i,j in pairs if frozenset((sequence[i-1],sequence[j-1])) not in CANON]
    return {'pairs':pairs,'pair_count':len(pairs),'crossing_pairs':crossing,'noncanonical_pairs':noncanonical}


def parse_bpseq(text, sequence):
    rows=[]
    for line in text.splitlines():
        line=line.strip()
        if not line or line.startswith('#'):continue
        fields=line.split()
        if len(fields)!=3:raise ValueError('BPSEQ requires 3 columns')
        i,base,j=int(fields[0]),fields[1].upper().replace('T','U'),int(fields[2])
        rows.append((i,base,j))
    n=len(sequence)
    if len(rows)!=n or [r[0] for r in rows]!=list(range(1,n+1)):
        raise ValueError('BPSEQ positions must cover 1..N')
    for i,base,j in rows:
        if base!=sequence[i-1]:raise ValueError(f'BPSEQ sequence mismatch at {i}')
        if not 0<=j<=n or j==i:raise ValueError(f'Invalid BPSEQ partner at {i}')
        if j and rows[j-1][2]!=i:raise ValueError(f'Nonreciprocal BPSEQ pairing at {i}')
    return sorted((i,j) for i,_,j in rows if i<j)


def bpseq_to_dotbracket(sequence,pairs):
    # Four bracket levels; fail loudly for structures requiring more levels.
    n=len(sequence)
    chars=['.']*n
    levels=[]
    for i,j in sorted(pairs):
        if not 1<=i<j<=n:raise ValueError('Pair out of range')
        level=next((k for k,assigned in enumerate(levels)
                    if all(not (a<i<b<j or i<a<j<b) for a,b in assigned)),None)
        if level is None:
            if len(levels)==len(OPEN):raise ValueError('Pseudoknot requires >4 bracket levels')
            levels.append([]);level=len(levels)-1
        if chars[i-1]!='.' or chars[j-1]!='.':raise ValueError('Residue paired twice')
        levels[level].append((i,j))
        chars[i-1]=OPEN[level];chars[j-1]=CLOSE[level]
    return ''.join(chars)


def annotate_sample(sample_dir, annotation_path, mode, source, method, reviewed=False):
    sample_dir=Path(sample_dir)
    metadata=json.loads((sample_dir/'metadata.json').read_text())
    sequence=metadata['sequence'].upper().replace('T','U')
    if mode not in ('oracle_2d','independent_2d'):raise ValueError('Invalid annotation mode')
    if not source.strip() or not method.strip():raise ValueError('Source and method are required')
    raw=Path(annotation_path).read_text().strip()
    if annotation_path.suffix.lower()=='.bpseq':
        pairs=parse_bpseq(raw,sequence)
        dbn=bpseq_to_dotbracket(sequence,pairs)
    else:
        lines=[s.strip() for s in raw.splitlines() if s.strip() and not s.startswith('>')]
        dbn=lines[-1]
    qc=parse_dotbracket(sequence,dbn)
    c4_missing=int(metadata.get('missing_c4_count',len(sequence)))
    contiguous=bool(metadata.get('contiguous',False))
    coverage=(len(sequence)-c4_missing)/len(sequence)
    reasons=[]
    if not 30<=len(sequence)<=200:reasons.append('length_outside_30_200')
    if coverage<0.95:reasons.append('c4_coverage_below_0.95')
    if not contiguous:reasons.append('noncontiguous_label_seq_ids')
    if qc['noncanonical_pairs']:reasons.append('noncanonical_pairs_require_review')
    if not reviewed:reasons.append('manual_review_pending')
    # Modified residues and protein/RNA complexes require manual provenance review.
    status='accepted' if not reasons else 'needs_review'
    ann_dir=sample_dir/'annotations'/mode
    ann_dir.mkdir(parents=True,exist_ok=True)
    (ann_dir/'secondary_structure.dbn').write_text(f'>{metadata["sample_id"]}|{mode}\n{sequence}\n{dbn}\n')
    report={'sample_id':metadata['sample_id'],'mode':mode,'status':status,'reasons':reasons,
            'source':source,'method':method,'manual_reviewed':reviewed,'sequence_length':len(sequence),
            'c4_coverage':coverage,'pair_count':qc['pair_count'],'pseudoknot':qc['crossing_pairs'],
            'noncanonical_pair_count':len(qc['noncanonical_pairs']),
            'annotation_sha256':hashlib.sha256(Path(annotation_path).read_bytes()).hexdigest(),
            'reference_cif':metadata.get('reference_cif'),
            'note':'Oracle-2D labels must not be reported as independent sequence-to-structure predictions.'}
    (ann_dir/'annotation_qc.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


def inventory(dataset):
    root=Path(dataset)
    rows=[]
    for metadata_path in sorted(root.glob('samples/*/metadata.json')):
        meta=json.loads(metadata_path.read_text())
        for mode in ('oracle_2d','independent_2d'):
            report_path=metadata_path.parent/'annotations'/mode/'annotation_qc.json'
            if report_path.exists():
                report=json.loads(report_path.read_text())
                status=report['status']; reasons=';'.join(report['reasons'])
            else:status='missing_annotation';reasons='annotation_not_supplied'
            rows.append({'sample_id':meta['sample_id'],'mode':mode,'length':len(meta['sequence']),
                         'status':status,'reasons':reasons})
    out=root/'curation_inventory.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sample_id','mode','length','status','reasons'])
        w.writeheader();w.writerows(rows)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('annotate')
    a.add_argument('--sample',required=True)
    a.add_argument('--annotation',required=True,type=Path)
    a.add_argument('--mode',choices=['oracle_2d','independent_2d'],required=True)
    a.add_argument('--source',required=True)
    a.add_argument('--method',required=True)
    a.add_argument('--reviewed',action='store_true')
    b=sub.add_parser('inventory');b.add_argument('--dataset',required=True)
    args=p.parse_args()
    result=annotate_sample(args.sample,args.annotation,args.mode,args.source,args.method,args.reviewed) if args.command=='annotate' else inventory(args.dataset)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
