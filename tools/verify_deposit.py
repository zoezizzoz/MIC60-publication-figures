#!/usr/bin/env python3
"""Verify the frozen graph tables and key source relationships without third-party packages."""
from pathlib import Path
import argparse, csv, hashlib, json, math, statistics, collections
ROOT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def verify():
    idx=json.loads((ROOT/'plotted_data/INDEX.json').read_text())
    assert len(idx)==43
    for item in idx:
        p=ROOT/item['file']
        assert p.is_file(),item['file']
        assert sha(p)==item['sha256'],'Modified plotted table: '+item['file']
        assert len(rows(p))==item['rows'],item['file']
        assert (ROOT/item['code']).is_file(),'Missing generator: '+item['code']
        assert (ROOT/item['source']).is_file(),'Missing source: '+item['source']
    assert len(rows(ROOT/'plotted_data/Fig1C/volcano_plotted_genes.csv'))==7456
    assert len(rows(ROOT/'plotted_data/FigS1C/volcano_plotted_genes.csv'))==9695
    assert len(rows(ROOT/'plotted_data/Fig4A/TIMELESS_DAPI_values.csv'))==64
    assert len(rows(ROOT/'plotted_data/Fig4B/mtt_plot_values.csv'))==78
    mtt=collections.defaultdict(list)
    for row in rows(ROOT/'plotted_data/Fig4B/mtt_plot_values.csv'):
        mtt[(int(row['dose_mM']),row['genotype'])].append(float(row['normalized_pct']))
    pooled=rows(ROOT/'plotted_data/Fig4B/mtt_pooled_well_summary.csv')
    assert len(pooled)==len(mtt)==10
    for row in pooled:
        dose=int(row['dose_mM']);values=mtt[(dose,row['genotype'])]
        n=len(values);sd=statistics.stdev(values)
        assert int(row['well_n'])==n==(6 if dose in (5,10) else 9)
        assert int(row['culture_preparations'])==(1 if dose in (5,10) else 2)
        assert math.isclose(float(row['mean']),statistics.mean(values),abs_tol=1e-10)
        assert math.isclose(float(row['sd']),sd,abs_tol=1e-10)
        assert math.isclose(float(row['sem']),sd/math.sqrt(n),abs_tol=1e-10)
    assert len(rows(ROOT/'plotted_data/FigS4C/transfection_efficiency.csv'))==20
    stress=rows(ROOT/'plotted_data/FigS4A/figS4_plotted_values.csv')
    de={r['gene']:r for r in rows(ROOT/'reviewed_analysis/data/rnaseq_female.csv')}
    assert len(stress)==64
    for r in stress:
        d=de[r['gene']]
        assert math.isclose(float(r['lfc_F']),float(d['log2FoldChange']),rel_tol=1e-10,abs_tol=1e-10)
        if r['padj_F'] not in ('','NA'):assert math.isclose(float(r['padj_F']),float(d['padj']),rel_tol=1e-10,abs_tol=1e-12)
        else:assert d['padj'] in ('','NA')
    s3=rows(ROOT/'plotted_data/FigS3/selected_gene_modules.csv')
    ampk=rows(ROOT/'plotted_data/FigS3/AMPK_grouped_gene_set_and_expression.csv')
    edges=rows(ROOT/'plotted_data/FigS3/AMPK_component_group_mapping.csv')
    definitions=rows(ROOT/'plotted_data/FigS3/AMPK_component_group_definitions.csv')
    assert len(s3)==86 and len({r['gene'] for r in s3})==86
    assert len(ampk)==131 and len({r['flybase_id'] for r in ampk})==131
    assert len(edges)==200 and len(definitions)==13
    assert {r['flybase_id'] for r in edges}=={r['flybase_id'] for r in ampk}
    group_counts=collections.Counter(r['component_group'] for r in ampk)
    assert group_counts=={r['component_group']:int(r['fly_gene_count']) for r in definitions}
    numeric=['baseMean','log2FoldChange','lfcSE','stat','pvalue','padj']
    for row in s3:
        for col in numeric:
            a,b=row[col],de[row['gene']][col]
            if a in ('','NA'):assert b in ('','NA'),(row['gene'],col)
            else:assert math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-14),(row['gene'],col)
    current={(row['module'],row['gene']):row for row in s3}
    old=rows(ROOT/'reviewed_analysis/figure_sources/Figure_S3_Selected_Gene_Modules/Supporting_Data/selected_gene_modules_previous.csv')
    unchanged=[row for row in old if row['module']!='AMPK']
    assert len(unchanged)==77
    for row in unchanged:assert all(current[(row['module'],row['gene'])][key]==value for key,value in row.items())
    package=ROOT/'reviewed_analysis/figure_sources/Figure_S3_Selected_Gene_Modules'
    evidence=rows(package/'Supporting_Data/AMPK_literature_selection.csv')
    focused=[r for r in s3 if r['module']=='AMPK-associated genes']
    assert [r['gene'] for r in focused]==[r['gene'] for r in evidence]
    assert len(focused)==9 and all(r['doi'] and r['published_fly_evidence'] for r in evidence)
    full={r['gene']:r for r in ampk}
    for row in focused:
        for col in numeric:assert row[col]==full[row['gene']][col]
    assert collections.Counter(r['status'] for r in focused)=={'Below DE cutoffs':5,'No fold-change estimate':3,'Adjusted P unavailable':1}
    flags=rows(package/'Rebuilt_Output/AMPK_full_survey_with_display_selection.csv')
    assert len(flags)==131
    assert {r['gene'] for r in flags if r['shown_in_focused_panel']=='true'}=={r['gene'] for r in focused}
    for row in ampk:
        for col in numeric:
            a,b=row[col],de[row['gene']][col]
            if a in ('','NA'):assert b in ('','NA')
            else:assert math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-14)
    assert {r['gene'] for r in focused} <= set(full)
    assert collections.Counter(row['status'] for row in ampk)=={'Below DE cutoffs':80,'Adjusted P unavailable':31,'No fold-change estimate':18,'Down in CS':2}
    assert {row['gene'] for row in ampk if row['status']=='Down in CS'}=={'Takl1','ninaD'}
    fields=rows(ROOT/'panels/Fig3N/Source_Data/TMRM_MTG_connected_object_image_summary.csv')
    objects=rows(ROOT/'panels/Fig3N/Source_Data/TMRM_MTG_connected_object_measurements.csv')
    groups=collections.defaultdict(list)
    for r in objects:groups[r['file']].append(float(r['tmrm_mtg_object_ratio']))
    assert len(objects)==22531 and len(fields)==32 and len(groups)==32
    wt=statistics.mean(float(r['image_mean_object_ratio']) for r in fields if r['condition']=='WT') if any(r['condition']=='WT' for r in fields) else statistics.mean(float(r['image_mean_object_ratio']) for r in fields if r['condition']=='WR')
    counts=collections.Counter(r['condition'] for r in fields)
    assert sorted(counts.values())==[7,25]
    errors=[]
    for r in fields:
        values=groups[r['file']];assert len(values)==int(r['objects_n'])
        mean=statistics.mean(values);error=abs(mean-float(r['image_mean_object_ratio']));errors.append(error)
        assert math.isclose(mean,float(r['image_mean_object_ratio']),rel_tol=1e-12,abs_tol=1e-12)
        assert math.isclose(100*mean/wt,float(r['image_mean_percent_of_WT']),rel_tol=1e-12,abs_tol=1e-10)
    for item in idx:
        if item['panel'] in ['Fig3B','Fig3C','Fig3E','Fig3F']:
            assert all(r['exclude']=='FALSE' for r in rows(ROOT/item['file']))
    return {'status':'passed','plotted_tables':len(idx),'graph_panels':len({r['panel'] for r in idx}),'S3_entries_checked':len(s3),'S3_AMPK_displayed_genes':len(focused),'S3_AMPK_full_survey_genes':len(ampk),'S3_mapping_edges':len(edges),'S3_display_groups':len(definitions),'S4A_DE_values_checked':64,'TMRM_objects':len(objects),'TMRM_fields':len(fields),'TMRM_max_mean_error':max(errors),'survival':'historical source retained; not certified against current curve'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--checksums',action='store_true');a=p.parse_args()
    report=verify()
    if a.checksums:
        items=json.loads((ROOT/'provenance/FILE_MANIFEST.json').read_text())['files']
        for item in items:
            path=ROOT/item['path'];assert path.is_file(),item['path'];assert sha(path)==item['sha256'],'Snapshot changed: '+item['path']
        report['snapshot_files_checked']=len(items)
    print(json.dumps(report,indent=2))
