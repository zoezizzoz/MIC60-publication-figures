"""Join independently defined membership to the unchanged archived contrast."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,json,math,hashlib

ROOT=Path(__file__).resolve().parents[1]
def read(p):return list(csv.DictReader(p.open()))
def write(p,rows,fields=None):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def number(x):
    try:return float(x)
    except (ValueError,TypeError):return math.nan
def same(a,b):
    a,b=number(a),number(b)
    return (math.isnan(a) and math.isnan(b)) or math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)

de=read(ROOT/'Supporting_Data/female_DE_source.csv'); by_symbol={r['gene']:r for r in de}
assert len(de)==len(by_symbol)
edges=read(ROOT/'Supporting_Data/AMPK_sources/diopt_retained_mapping_edges.csv')
group=defaultdict(list)
for edge in edges:group[edge['flybase_id']].append(edge)
annot=read(ROOT/'Supporting_Data/AMPK_sources/diopt_fly_identifiers.csv')
canonical=defaultdict(set)
for a in annot:canonical[a['FLYBASE']].add(a['SYMBOL'])
numeric=['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
ampk=[]
for fbgn,ee in group.items():
    symbols={e['fly_symbol'] for e in ee}
    assert len(symbols)==1
    symbol=next(iter(symbols))
    # Exact canonical symbol match; never use the DE outcome to choose a mapping.
    assert symbol in canonical[fbgn],(fbgn,symbol,canonical[fbgn])
    assert symbol in by_symbol,(fbgn,symbol)
    source=by_symbol[symbol];fc,p=number(source['log2FoldChange']),number(source['padj'])
    state='No fold-change estimate' if not math.isfinite(fc) else 'Adjusted P unavailable' if not math.isfinite(p) else 'Up in CS' if p<.05 and fc>=.58 else 'Down in CS' if p<.05 and fc<=-.58 else 'Below DE cutoffs'
    ampk.append({'module':'Fly homologs of human AMPK-pathway genes','gene':symbol,**{k:source[k] for k in numeric},
        'status':state,'flybase_id':fbgn,'human_entrez_ids':';'.join(sorted({e['human_entrez'] for e in ee})),
        'human_symbols':';'.join(sorted({e['human_symbol'] for e in ee})),
        'diopt_confidence_levels':';'.join(sorted({e['confidence'] for e in ee})),
        'mapping_rule':'Exact canonical symbol in DIOPT, org.Dm.eg.db 3.19.1 and archived DE table',
        'source':'KEGG hsa04152; DIOPT 9.1 high/moderate confidence'})
ampk.sort(key=lambda r:(r['gene'].casefold(),r['gene']))
before=read(ROOT/'Supporting_Data/selected_gene_modules_previous.csv')
unchanged=[r for r in before if r['module']!='AMPK']
fields=list(ampk[0])
candidate=[]
for module in ['FOXO','DNA replication','Spargel','Fly homologs of human AMPK-pathway genes','Chromatin','Checkpoint']:
    candidate.extend(ampk if module=='Fly homologs of human AMPK-pathway genes' else [{k:r.get(k,'') for k in fields} for r in unchanged if r['module']==module])
assert len(ampk)==131 and len(candidate)==208
assert len({r['flybase_id'] for r in ampk})==len(ampk)
for row in candidate:
    for k in numeric:assert same(row[k],by_symbol[row['gene']][k]),(row['gene'],k)
for row in unchanged:
    match=next(c for c in candidate if c['module']==row['module'] and c['gene']==row['gene'])
    assert all(match[k]==v for k,v in row.items()),row['gene']
write(ROOT/'Rebuilt_Output/AMPK_gene_set_and_expression.csv',ampk)
write(ROOT/'Rebuilt_Output/FigS3_supporting_data.csv',candidate)
old_ampk={r['gene'] for r in before if r['module']=='AMPK'}
new_ampk={r['gene'] for r in ampk}
delta=[]
for gene in sorted(old_ampk|new_ampk,key=str.casefold):
    delta.append({'gene':gene,'in_previous_AMPK':gene in old_ampk,'in_documented_AMPK':gene in new_ampk,
        'change':'retained' if gene in old_ampk&new_ampk else 'added' if gene in new_ampk else 'removed',
        'reason':'Membership determined solely by KEGG hsa04152 and high/moderate DIOPT mapping'})
write(ROOT/'Rebuilt_Output/AMPK_membership_changes.csv',delta)

summary={'human_pathway_genes':122,'human_genes_with_retained_fly_match':115,'unique_fly_members':131,
    'finite_fold_change':sum(math.isfinite(number(r['log2FoldChange'])) for r in ampk),
    'adjusted_P_available':sum(math.isfinite(number(r['padj'])) for r in ampk),
    'plot_categories':dict(Counter(r['status'] for r in ampk)),
    'all_symbols_exact_canonical_matches':True,'other_five_modules_unchanged_rows':len(unchanged),
    'source_value_checks':'All six numeric columns match archived female DE; old module rows remain byte-identical by field',
    'old_AMPK_removed':sorted(old_ampk-new_ampk),'old_AMPK_retained':sorted(old_ampk&new_ampk),
    'new_AMPK_added_count':len(new_ampk-old_ampk)}
(ROOT/'QA/membership_validation.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='original_files'},indent=2))
