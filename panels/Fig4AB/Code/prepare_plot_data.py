"""Extract immutable input copies using the archived assay-audit mapping."""
from pathlib import Path
import csv, json, hashlib, statistics
from openpyxl import load_workbook
ROOT = Path(__file__).resolve().parents[1]
def write_csv(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
a = ROOT/'Sources/Figure_4A_TIMELESS';b = ROOT/'Sources/Figure_4B_MTT'
w = load_workbook(a/'Original_Data/Fig4A.xlsx', data_only=True, read_only=True)
sheet = w.active
assert tuple(next(sheet.values)) == ('WT','CS','WT-H2O2','CS-H2O2')
rows=[]
for col,group in enumerate(('WT','CS','WT +H2O2','CS +H2O2'),1):
    group_rows=[]
    for rr in range(2,sheet.max_row+1):
        value=sheet.cell(rr,col).value
        if value is None:continue
        assert isinstance(value,(int,float)) and value>=0
        group_rows.append({'group':group,'value':value,'source_sheet':sheet.title,'source_cell':sheet.cell(rr,col).coordinate})
    assert len(group_rows)==16
    rows.extend(group_rows)
write_csv(a/'Supporting_Data/TIMELESS_plotted_values.csv',rows)
old = ROOT.parents[1]/'reviewed_analysis/figure_sources/Figure_4C_TIMELESS/Supporting_Data/Fig4A_TIMELESS_DAPI_plotted_values.csv'
with old.open() as f: prior=list(csv.DictReader(f))
assert [(r['group'],r['value']) for r in rows]==[(r['group'],float(r['value'])) for r in prior]
wa=load_workbook(b/'Original_Data/MTT Data.xlsx',data_only=True,read_only=True)
wb=load_workbook(b/'Original_Data/7-28-26-MTT.xlsx',data_only=True,read_only=True)
for name in wa.sheetnames:
    assert list(wa[name].values)==list(wb[name].values),name
records=[]
for date,exp,sheetname,blankcells,wellrows,doses in [
 ('2026-07-15',1,'7.15.26 raw',[(6,3),(6,4)],range(3,6),[(0,(4,5)),(20,(6,7)),(40,(8,9))]),
 ('2026-07-28',2,'7.28.26 raw',[(35,2),(35,3),(35,4)],range(29,35),[(0,(4,5)),(5,(12,13)),(10,(10,11)),(20,(8,9)),(40,(6,7))])]:
    sh=wa[sheetname];blank=statistics.mean(sh.cell(r,c).value for r,c in blankcells)
    for dose,cols in doses:
        for genotype,col in zip(('WT','CS'),cols):
            for well,row in enumerate(wellrows,1):
                raw=sh.cell(row,col).value
                assert isinstance(raw,(int,float))
                records.append(dict(experiment=exp,acquisition_date=date,dose_mM=dose,genotype=genotype,well=well,source_sheet=sheetname,source_cell=sh.cell(row,col).coordinate,raw_absorbance=raw,blank=blank,corrected_absorbance=raw-blank))
for r in records:
    controls=[x['corrected_absorbance'] for x in records if x['experiment']==r['experiment'] and x['genotype']==r['genotype'] and x['dose_mM']==0]
    r['matched_control_mean']=statistics.mean(controls)
    assert r['matched_control_mean']>0
    r['normalized_pct']=100*r['corrected_absorbance']/r['matched_control_mean']
assert len(records)==78
write_csv(b/'Supporting_Data/mtt_plot_values.csv',records)
counts={str(d):sum(r['dose_mM']==d and r['genotype']=='WT' for r in records) for d in [0,5,10,20,40]}
assert list(counts.values())==[9,6,6,9,9]
report={'TIMELESS':{'observations':len(rows),'per_condition':16,'matches_recovered_values_exactly':True,'range':[min(r['value'] for r in rows),max(r['value'] for r in rows)],'inference':'Descriptive; cell/image/culture/experiment identifiers are absent. The later September audit removed legacy P=.0477/.0218.'},'MTT':{'wells':78,'wells_per_genotype_by_dose':counts,'culture_preparations':{'1':'2026-07-15','2':'2026-07-28'},'raw_workbook_numeric_and_text_contents_identical':True,'normalization':'100 * (raw absorbance - within-experiment mean blank) / mean corrected matched 0 mM genotype control','inference':'Mean and SEM across normalized wells pooled within genotype and dose; each well is treated as a biological replicate per author direction. No inferential tests.','source_note':'July 15 control is labelled DMSO. Preserved the matched-control normalization; vehicle equivalence is not established here.'},'input_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in ROOT.glob('Sources/*/Original_Data/*.xlsx')}}
(ROOT/'Documentation/DATA_VERIFICATION.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
