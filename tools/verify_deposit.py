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
    assert len(idx)==40
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
    assert len(rows(ROOT/'plotted_data/FigS4C/transfection_efficiency.csv'))==20
    stress=rows(ROOT/'plotted_data/FigS4A/figS4_plotted_values.csv')
    de={r['gene']:r for r in rows(ROOT/'reviewed_analysis/data/rnaseq_female.csv')}
    assert len(stress)==64
    for r in stress:
        d=de[r['gene']]
        assert math.isclose(float(r['lfc_F']),float(d['log2FoldChange']),rel_tol=1e-10,abs_tol=1e-10)
        if r['padj_F'] not in ('','NA'):assert math.isclose(float(r['padj_F']),float(d['padj']),rel_tol=1e-10,abs_tol=1e-12)
        else:assert d['padj'] in ('','NA')
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
    return {'status':'passed','plotted_tables':len(idx),'graph_panels':len({r['panel'] for r in idx}),'S4A_DE_values_checked':64,'TMRM_objects':len(objects),'TMRM_fields':len(fields),'TMRM_max_mean_error':max(errors),'survival':'historical source retained; not certified against current curve'}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--checksums',action='store_true');a=p.parse_args()
    report=verify()
    if a.checksums:
        items=json.loads((ROOT/'provenance/FILE_MANIFEST.json').read_text())['files']
        for item in items:
            path=ROOT/item['path'];assert path.is_file(),item['path'];assert sha(path)==item['sha256'],'Snapshot changed: '+item['path']
        report['snapshot_files_checked']=len(items)
    print(json.dumps(report,indent=2))
