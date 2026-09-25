#!/usr/bin/env python3
"""Rebuild S3 offline from the frozen KEGG/DIOPT records and archived female DE."""
from pathlib import Path
import argparse,csv,hashlib,json,math,os,shutil,subprocess,sys
from collections import Counter
from export_png import export_png
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'Rebuilt_Output'
LABEL='Fly homologs of human AMPK-pathway genes'
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def number(s):
    try:return float(s)
    except (TypeError,ValueError):return math.nan
def same(a,b):
    a,b=number(a),number(b)
    return math.isnan(a) and math.isnan(b) or math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'));a=p.parse_args()
    rscript=shutil.which(a.rscript)
    if not rscript:p.error('Rscript not found; supply --rscript /path/to/Rscript')
    (ROOT/'QA').mkdir(exist_ok=True);OUT.mkdir(exist_ok=True)
    sources=json.loads((ROOT/'Supporting_Data/FROZEN_SOURCE_MANIFEST.json').read_text())
    for item in sources:
        path=ROOT/item['path']
        if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('Frozen input changed: '+str(path))
    for script in ['prepare_diopt.py','assemble_ampk_data.py','group_ampk_components.py']:
        subprocess.run([sys.executable,str(ROOT/'Code'/script)],check=True,cwd=ROOT)
    data=read(OUT/'FigS3_supporting_data.csv');ampk=read(OUT/'AMPK_grouped_gene_set_and_expression.csv')
    by_gene={r['gene']:r for r in ampk};group_cols=[c for c in ampk[0] if c not in data[0]]
    for row in data:
        for col in group_cols:row[col]=by_gene[row['gene']][col] if row['module']==LABEL else ''
    write(OUT/'FigS3_supporting_data.csv',data)
    de={r['gene']:r for r in read(ROOT/'Supporting_Data/female_DE_source.csv')}
    old=read(ROOT/'Supporting_Data/selected_gene_modules_previous.csv')
    numeric=['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
    assert len(data)==208 and len(ampk)==131
    assert len({r['flybase_id'] for r in ampk})==131
    assert len({r['component_group'] for r in ampk})==13
    for row in data:
        assert all(same(row[col],de[row['gene']][col]) for col in numeric),row['gene']
    current={(r['module'],r['gene']):r for r in data}
    for row in old:
        if row['module']=='AMPK':continue
        assert all(current[(row['module'],row['gene'])][k]==v for k,v in row.items())
    states=Counter(r['status'] for r in ampk)
    assert states=={'Below DE cutoffs':80,'Adjusted P unavailable':31,'No fold-change estimate':18,'Down in CS':2},states
    assert {r['gene'] for r in ampk if r['status']=='Down in CS'}=={'Takl1','ninaD'}
    # Preserve the comprehensive table before assembling the literature-selected figure.
    shutil.copy2(OUT/'FigS3_supporting_data.csv',OUT/'FigS3_full_orthology_survey.csv')
    selection=read(ROOT/'Supporting_Data/AMPK_literature_selection.csv')
    assert len(selection)==11 and len({r['gene'] for r in selection})==11
    selected=[]
    for evidence in selection:
        gene=evidence['gene']
        # Publication-supported fly genes need not belong to the human-orthology survey.
        row=dict(by_gene[gene]) if gene in by_gene else dict(next(r for r in old if r['module']=='AMPK' and r['gene']==gene))
        row['in_full_orthology_survey']=str(gene in by_gene).lower()
        row['module']='AMPK-associated genes'
        row.update({k:v for k,v in evidence.items() if k!='gene'})
        selected.append(row)
    write(OUT/'AMPK_literature_selected_gene_set_and_expression.csv',selected)
    keep={r['gene'] for r in selection}
    flagged=[dict(row,shown_in_focused_panel=str(row['gene'] in keep).lower()) for row in ampk]
    write(OUT/'AMPK_full_survey_with_display_selection.csv',flagged)
    decisions=read(ROOT/'Supporting_Data/S3_original_gene_display_decisions.csv')
    assert [(r['original_module'],r['gene']) for r in decisions]==[(r['module'],r['gene']) for r in old]
    old_lookup={(r['module'],r['gene']):r for r in old}
    figure=[]
    for module in ['FOXO','DNA replication','Spargel','AMPK','Chromatin','Checkpoint']:
        if module=='AMPK':
            figure.extend(selected)
            continue
        for decision in decisions:
            if decision['original_module']!=module or decision['shown_in_revised_figure']!='true':continue
            row=dict(old_lookup[(module,decision['gene'])])
            row['original_module']=module
            row['module']=decision['display_module']
            row.update({k:v for k,v in decision.items() if k not in ['gene','original_module','display_module']})
            figure.append(row)
    columns=list(dict.fromkeys(k for row in figure for k in row))
    figure=[{k:row.get(k,'') for k in columns} for row in figure]
    assert len(figure)==81
    for row in figure:
        assert all(same(row[col],de[row['gene']][col]) for col in numeric),row['gene']
    write(OUT/'FigS3_supporting_data.csv',figure)
    shutil.copy2(OUT/'FigS3_supporting_data.csv',ROOT/'Supporting_Data/FigS3_supporting_data.csv')
    shutil.copy2(ROOT/'Documentation/FigS3_figure_legend.txt',OUT/'FigS3_figure_legend.txt')
    subprocess.run([rscript,str(ROOT/'Code/generate_AMPK_full_survey.R')],check=True,cwd=ROOT)
    subprocess.run([rscript,str(ROOT/'Code/generate_FigS3_selected_modules.R')],check=True,cwd=ROOT)
    dimensions={}
    for stem in ['FigS3_selected_gene_modules','AMPK_signaling_components_grouped_panel','AMPK_literature_selected_panel']:
        dimensions[stem]=export_png(OUT/(stem+'.pdf'),OUT/(stem+'.png'))
    report={'frozen_inputs_checked':len(sources),'module_entries':len(figure),'full_survey_module_entries':len(data),'focused_AMPK_genes':len(selected),'focused_AMPK_categories':dict(Counter(r['status'] for r in selected)),'unique_fly_genes_across_modules':len({r['gene'] for r in figure}),'full_survey_unique_fly_genes_across_modules':len({r['gene'] for r in data}),'AMPK_unique_fly_genes':131,'AMPK_display_groups':13,'other_modules_retained_entries':70,'original_entries_audited':89,'excluded_from_previous_figure':7,'original_AMPK_members_restored':['Sesn','Atg8a'],'focused_AMPK_outside_orthology_survey':[r['gene'] for r in selected if r['in_full_orthology_survey']=='false'],'display_module_counts':dict(Counter(r['module'] for r in figure)),'all_six_numeric_fields_match_frozen_female_contrast':True,'AMPK_categories':dict(states),'PNG_dpi':600,'PNG_dimensions':dimensions,'DESeq2_rerun':False,'network_access':False}
    (ROOT/'QA/rebuild_validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
