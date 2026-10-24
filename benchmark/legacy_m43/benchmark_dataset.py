"""Curate experimental RNA chains from RCSB mmCIF; no tertiary-derived 2D labels."""
import argparse, csv, json, re, urllib.request
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

BASES = {'A':'A','C':'C','G':'G','U':'U','I':'I'}
URL = 'https://files.rcsb.org/download/{pdb}.cif'

def _list(d, key):
    x=d.get(key,[])
    return x if isinstance(x,list) else [x]

def extract_chains(cif, low=30, high=200):
    """Return chain candidates with distinct polymer residue positions, not atom counts."""
    d=MMCIF2Dict(str(cif))
    cols=['_atom_site.label_asym_id','_atom_site.label_seq_id','_atom_site.label_comp_id','_atom_site.label_atom_id','_atom_site.label_alt_id','_atom_site.pdbx_PDB_model_num']
    data=[_list(d,c) for c in cols]
    if not data[0]: raise ValueError('No atom_site records')
    if len(set(map(len,data)))!=1: raise ValueError('Inconsistent atom_site columns')
    chain_res={}
    for chain,pos,res,atom,alt,model in zip(*data):
        if model!='1' or alt not in ('.','?','A','1') or not pos.isdigit(): continue
        key=(chain,int(pos))
        if res in BASES:
            chain_res.setdefault(chain,{}).setdefault(int(pos),{'base':BASES[res], 'atoms':set()})['atoms'].add(atom)
    result=[]
    for chain,rs in sorted(chain_res.items()):
        positions=sorted(rs)
        n=len(positions)
        if not low<=n<=high:continue
        missing_c4=sum("C4'" not in rs[p]['atoms'] and 'C4*' not in rs[p]['atoms'] for p in positions)
        result.append({'label_asym_id':chain,'length':n,'sequence':''.join(rs[p]['base'] for p in positions),
                       'label_seq_ids':positions,'missing_c4_count':missing_c4,
                       'contiguous':positions==list(range(positions[0],positions[-1]+1)),
                       'modified_residues_excluded':True})
    return result

def fetch(pdb, out, timeout=30):
    out.parent.mkdir(parents=True,exist_ok=True)
    req=urllib.request.Request(URL.format(pdb=pdb.upper()),headers={'User-Agent':'RNAFold3D-Bench/0.1 (research)'})
    with urllib.request.urlopen(req,timeout=timeout) as r: data=r.read()
    if not data or b'_atom_site.' not in data:raise ValueError('Not a valid atom_site mmCIF response')
    out.write_bytes(data)
    return out

def curate(manifest, root, download=True):
    root=Path(root); root.mkdir(parents=True,exist_ok=True)
    entries=json.loads(Path(manifest).read_text())['entries']; rows=[]
    for e in entries:
        pdb=e['pdb_id'].upper(); raw=root/'raw'/f'{pdb}.cif'
        try:
            if not raw.exists():
                if not download:raise FileNotFoundError('raw file absent (offline mode)')
                fetch(pdb,raw)
            chains=extract_chains(raw)
            if not chains: raise ValueError('No canonical-only 30-200 nt RNA chain detected')
            for c in chains:
                sid=f'{pdb}_{c["label_asym_id"]}'
                dest=root/'samples'/sid;dest.mkdir(parents=True,exist_ok=True)
                (dest/'sequence.fasta').write_text(f'>{sid}\n{c["sequence"]}\n')
                (dest/'metadata.json').write_text(json.dumps({'sample_id':sid,'pdb_id':pdb,'source':e['source'],
                    'reference_cif':f'../../raw/{pdb}.cif','secondary_structure_status':'NOT_PROVIDED',
                    'secondary_structure_provenance':None,**c},indent=2)+'\n')
                rows.append({'sample_id':sid,'pdb_id':pdb,'chain':c['label_asym_id'],'length':c['length'],
                             'missing_c4':c['missing_c4_count'],'contiguous':c['contiguous'],
                             'status':'needs_2d_and_review','error':''})
        except Exception as ex:
            rows.append({'sample_id':'','pdb_id':pdb,'chain':'','length':'','missing_c4':'',
                         'contiguous':'','status':'fetch_or_filter_failed','error':str(ex)[:240]})
    with (root/'inventory.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['sample_id','pdb_id','chain','length','missing_c4','contiguous','status','error']);writer.writeheader();writer.writerows(rows)
    print(json.dumps({'candidate_pdbs':len(entries),'rna_chains':sum(x['status']=='needs_2d_and_review' for x in rows),
                      'failures':sum(x['status']=='fetch_or_filter_failed' for x in rows),'inventory':str(root/'inventory.csv')},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',default='candidate_manifest.json');p.add_argument('--out',default='dataset');p.add_argument('--offline',action='store_true');a=p.parse_args()
    curate(a.manifest,a.out,not a.offline)
