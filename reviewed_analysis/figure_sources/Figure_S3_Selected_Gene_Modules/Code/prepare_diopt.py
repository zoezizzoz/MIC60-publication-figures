"""Archive AMPK pathway membership and all DIOPT mapping responses.

Version 9.1 is returned by the working v9 API. The v10 API currently reports
'Version not suppported'; never describe this output as DIOPT v10.
Membership selection is independent of RNA-seq values.
"""
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json, re, csv, hashlib

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'Supporting_Data/AMPK_sources'
RAW=SRC/'diopt_v9_1'
RAW.mkdir(exist_ok=True)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def csvout(name,rows):
    with (SRC/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

humans={};section=''
for line in (SRC/'hsa04152.txt').read_text().splitlines():
    if line[:12].strip():section=line[:12].strip()
    if section=='GENE':
        m=re.match(r'(\d+)\s+([^;]+);',line[12:])
        if m:humans[m.group(1)]=m.group(2)
assert len(humans)==122

def query(gene):
    dest=RAW/(gene+'.json')
    url=f'https://www.flyrnai.org/tools/diopt/web/diopt_api/v9/get_orthologs_from_entrez/9606/{gene}/7227/none'
    if not dest.exists():
        raise FileNotFoundError('Frozen DIOPT response is missing: '+str(dest))
    return gene,json.loads(dest.read_text()),{'file':str(dest.relative_to(ROOT)),'url':url,'sha256':sha(dest)}

rows=[];audits=[];manifest=[]
with ThreadPoolExecutor(max_workers=2) as pool:
    for i,(gene,obj,meta) in enumerate(pool.map(query,humans),1):
        manifest.append(meta)
        results=obj.get('results') or {}
        matches=results.get(gene,{}) if isinstance(results,dict) else {}
        keep=[]
        for fly_entrez,m in matches.items():
            assert int(m['species_id'])==7227
            score=int(m['score']);best=m['best_score']=='Yes';rev=m['best_score_rev']=='Yes'
            rank='high' if best and rev and score>=2 else 'moderate' if score>=4 or (score>=2 and (best or rev)) else 'low'
            if m.get('confidence'): assert m['confidence'].lower()==rank,(gene,m,rank)
            retained=rank in ['high','moderate']
            fbgn=m['species_specific_geneid']
            assert m['species_specific_geneid_type']=='FlyBase' and fbgn.startswith('FBgn')
            if retained:keep.append(fbgn)
            rows.append({'human_entrez':gene,'human_symbol':humans[gene],
                'fly_entrez':fly_entrez,'flybase_id':fbgn,'fly_symbol':m['symbol'],
                'diopt_score':score,'max_score':m['max_score'],'best_forward':best,
                'best_reverse':rev,'confidence':rank,'retained':retained,
                'methods':';'.join(m['methods'])})
        audits.append({'human_entrez':gene,'human_symbol':humans[gene],
            'retained_flybase_ids':';'.join(sorted(set(keep))),
            'mapping_status':'mapped' if keep else 'no_high_or_moderate_confidence_fly_match',
            'all_candidate_match_count':len(matches)})
        if i%20==0 or i==len(humans):print(f'Archived {i}/{len(humans)} human-gene responses',flush=True)

csvout('diopt_all_mapping_edges.csv',rows)
csvout('diopt_human_membership_audit.csv',audits)
kept=[r for r in rows if r['retained']]
csvout('diopt_retained_mapping_edges.csv',kept)
ids=sorted({r['flybase_id'] for r in kept})
(SRC/'diopt_flybase_ids.txt').write_text('\n'.join(ids)+'\n')
summary={'retrieved_utc':json.loads((SRC/'diopt_manifest.json').read_text())['retrieved_utc'],
    'pathway':'KEGG hsa04152 AMPK signaling pathway','human_gene_count':len(humans),
    'diopt_version':9.1,'selection':'All high- and moderate-confidence DIOPT matches; deduplicate by FlyBase ID',
    'mapped_human_count':sum(a['mapping_status']=='mapped' for a in audits),
    'unique_fly_count':len(ids),'retained_edge_count':len(kept),
    'low_confidence_edges_excluded':len(rows)-len(kept),
    'rank_rule':'high: best both directions and score >= 2; moderate: score >= 4 or best either direction and score >= 2; otherwise low',
    'source_pathway_sha256':sha(SRC/'hsa04152.txt'),'api_responses':manifest}
(SRC/'diopt_manifest.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='api_responses'},indent=2))
