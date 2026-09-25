"""Add transparent display groups without changing membership or expression.

Groups describe the human source components in KEGG hsa04152. They are display
annotations created for this figure, not official KEGG subpathway gene sets or
experimental validation of a fly gene's function. All mappings are retained.
"""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,json,xml.etree.ElementTree as ET,hashlib

ROOT=Path(__file__).resolve().parents[1]
URL='https://www.kegg.jp/pathway/hsa04152'
def read(name):return list(csv.DictReader((ROOT/name).open()))
def write(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

groups=[
 (1,1,'AMPK subunits','PRKAA1 PRKAA2 PRKAB1 PRKAB2 PRKAG1 PRKAG2 PRKAG3',
  'Catalytic alpha and regulatory beta/gamma subunits of the human AMPK complex.'),
 (1,2,'Upstream signalling','CAMKK2 STK11 CAB39 CAB39L STRADA STRADB MAP3K7 AKT1 AKT2 AKT3 PDPK1 INS INSR IGF1R IRS1 IRS2 IRS4 PIK3CA PIK3CB PIK3CD PIK3R1 PIK3R2 PIK3R3 P3R3URF-PIK3R3 ADIPOR1 ADIPOR2 ADRA1A',
  'Upstream kinases and cofactors, hormone/receptor signalling and the insulin/PI3K-Akt branch shown in the human map; not all are direct AMPK activators.'),
 (1,3,'PP2A-related proteins','PPP2CA PPP2CB PPP2R1A PPP2R1B PPP2R2A PPP2R2B PPP2R2C PPP2R2D PPP2R3A PPP2R3B PPP2R5A PPP2R5B PPP2R5C PPP2R5D PPP2R5E PPP2R3C',
  'Human PP2A catalytic, scaffold and regulatory subunits. DIOPT matches include related fly phosphatases; membership does not establish PP2A-complex membership in flies.'),
 (2,1,'Glucose metabolism','FBP1 FBP2 G6PC1 G6PC2 G6PC3 GYS1 GYS2 PCK1 PCK2 PFKFB1 PFKFB2 PFKFB3 PFKFB4 PFKL PFKM PFKP',
  'Human enzymes involved in glycolysis, gluconeogenesis or glycogen synthesis.'),
 (2,2,'GLUT / TBC1D1 / Rab-related','SLC2A4 TBC1D1 RAB10 RAB8A RAB14 RAB2A RAB11B',
  'Human GLUT4, TBC1D1 and Rab components associated with the glucose-transporter trafficking branch.'),
 (2,3,'CFTR-related proteins','CFTR',
  'Fly proteins mapped to the human CFTR component; this family label does not imply that every fly match is a CFTR chloride channel.'),
 (3,1,'Lipid metabolic enzymes','ACACA ACACB FASN SCD SCD5 CPT1A CPT1B CPT1C LIPE HMGCR',
  'Human enzymes involved in fatty-acid synthesis, desaturation or oxidation, lipolysis, and cholesterol synthesis.'),
 (3,2,'CD36-related proteins','CD36',
  'Fly proteins mapped to the human CD36 lipid-transport component; family members can have additional or different fly functions.'),
 (3,3,'CIDEA-related proteins','CIDEA',
  'Fly proteins mapped to human CIDEA, shown in the map as an AMPK-linked regulator; fly Drep family membership does not establish a lipid-storage function.'),
 (4,1,'mTOR / translation / autophagy','MTOR RPTOR AKT1S1 TSC1 TSC2 RHEB RPS6KB1 RPS6KB2 EIF4EBP1 EEF2 ULK1',
  'Human mTOR complex/regulators, translation effectors and ULK1 autophagy initiator.'),
 (4,2,'Transcriptional regulators','PPARGC1A SIRT1 FOXO1 FOXO3 HNF4A PPARG SREBF1 CREB1 CREB3 CREB3L1 CREB3L2 CREB3L3 CREB3L4 CREB5 CRTC2',
  'Human transcription factors, coactivators and the SIRT1 regulator represented in the pathway.'),
 (4,3,'RNA-binding proteins','ELAVL1',
  'Fly proteins mapped to the human ELAVL1/HuR RNA-binding component.'),
 (4,4,'Cyclins','CCND1 CCNA1 CCNA2',
  'Fly proteins mapped to the human cyclin A and D cell-cycle components.')
]
by_human={}
for col,rank,label,symbols,rationale in groups:
    for sym in symbols.split():
        assert sym not in by_human,sym
        by_human[sym]=(col,rank,label,rationale)

kgml=ROOT/'Supporting_Data/AMPK_sources/hsa04152_grouping_reference.kgml'
if not kgml.exists():
    raise FileNotFoundError('Frozen KEGG grouping graph is missing: '+str(kgml))
nodes=defaultdict(list)
for e in ET.fromstring(kgml.read_bytes()).findall('entry'):
    if e.attrib.get('type')=='gene':
        for gene in e.attrib['name'].split():nodes[gene.split(':')[1]].append(e.attrib['id'])

edges=read('Supporting_Data/AMPK_sources/diopt_retained_mapping_edges.csv')
source_audit=read('Supporting_Data/AMPK_sources/diopt_human_membership_audit.csv')
assert {a['human_entrez'] for a in source_audit}==set(nodes),'KEGG source membership changed; use the frozen entry before revising groups'
mapped_by_fly=defaultdict(list);trace=[]
for e in edges:
    col,rank,label,why=by_human[e['human_symbol']]
    mapped_by_fly[e['flybase_id']].append((col,rank,label,why))
    trace.append({**e,'display_column':col,'component_group':label,
        'kegg_node_ids':';'.join(nodes[e['human_entrez']]),'grouping_rationale':why,
        'grouping_reference':URL,'annotation_scope':'Human-source role or protein family; fly role inferred by orthology'})

ampk=read('Rebuilt_Output/AMPK_gene_set_and_expression.csv')
out=[]
for row in ampk:
    labels=set(mapped_by_fly[row['flybase_id']])
    assert len(labels)==1,(row['gene'],labels)
    col,rank,label,why=next(iter(labels))
    core_order={'AMPKalpha':1,'alc':2,'SNF4Agamma':3}
    out.append({**row,'display_column':col,'component_group_order':rank,'component_group':label,
        'subunit_display_order':core_order.get(row['gene'],0) if label=='AMPK subunits' else 0,
        'grouping_rationale':why,'grouping_reference':URL})
out.sort(key=lambda x:(x['display_column'],x['component_group_order'],x['subunit_display_order'],x['gene'].casefold()))
assert len(out)==131 and len({r['flybase_id'] for r in out})==131
original={r['flybase_id']:r for r in ampk}
assert all(all(row[k]==v for k,v in original[row['flybase_id']].items()) for row in out)
write('Rebuilt_Output/AMPK_grouped_gene_set_and_expression.csv',out)
write('Rebuilt_Output/AMPK_component_group_mapping.csv',trace)
defs=[]
for col,rank,label,symbols,why in groups:
    actual=[r for r in out if r['component_group']==label]
    defs.append({'display_column':col,'component_group_order':rank,'component_group':label,
        'fly_gene_count':len(actual),'source_human_symbols':';'.join(sorted({e['human_symbol'] for e in trace if e['component_group']==label})),
        'rationale':why,'reference':URL})
write('Rebuilt_Output/AMPK_component_group_definitions.csv',defs)
summary={'created_utc':datetime.now(timezone.utc).isoformat(),'unchanged_unique_fly_genes':131,
    'membership_and_all_existing_fields_unchanged':True,'multiple_group_assignment_count':0,
    'groups':{r['component_group']:r['fly_gene_count'] for r in defs},
    'column_gene_counts':dict(Counter(r['display_column'] for r in out)),
    'kgml_sha256':hashlib.sha256(kgml.read_bytes()).hexdigest(),
    'definition':'Display groups based on human-source pathway roles/families, manually specified and traced to source KEGG nodes; not new gene sets.'}
(ROOT/'QA/grouping_validation.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
