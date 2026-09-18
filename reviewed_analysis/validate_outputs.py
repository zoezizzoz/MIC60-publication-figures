#!/usr/bin/env python3
"""Validate delivered outputs and shared statistical sources without refitting analyses."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import pdfplumber
from publish_figures import SOURCES
R=Path(__file__).resolve().parent.parent
checks={}
manifest=json.loads((R/'Figures/manifest.json').read_text())
manifest_files={m['file'] for m in manifest}
expected_generated={f'{ext.upper()}/{name}.{ext}' for name in SOURCES for ext in ('pdf','png')}
assert expected_generated.issubset(manifest_files)
assert len(manifest)==len(manifest_files) and len(manifest) % 2 == 0
for m in manifest:
 p=R/'Figures'/m['file'];src=R/m['generated_source']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
 assert hashlib.sha256(src.read_bytes()).hexdigest()==m['generated_source_sha256']
 assert (p.read_bytes()==src.read_bytes())==m['source_matches_delivery']
 assert m['delivery_type'] in {'generated','manually_refined_export'}
checks['manifest_files']=len(manifest)
checks['single_page_PDFs']=0
for p in (R/'Figures/PDF').glob('*.pdf'):
 with pdfplumber.open(p) as d:
  assert len(d.pages)==1,p
  page=d.pages[0];text=page.extract_text() or ''
  assert 'italic(' not in text and 'plain(' not in text,p
  assert all(c['x0']>=-1 and c['x1']<=page.width+1 and c['top']>=-1 and c['bottom']<=page.height+1 for c in page.chars),p
  checks['single_page_PDFs']+=1
b=json.loads((R/'Review/go_reference_style_numeric_baseline.json').read_text())
assert all(hashlib.sha256((R/p).read_bytes()).hexdigest()==h for p,h in b.items())
checks['prior_numerical_tables_unchanged']=len(b)
A=R/'Analysis';out=A/'results/full_figures'
for sex in ['female','male']:
 de=pd.read_csv(A/f'data/rnaseq_{sex}.csv');sig=de[de.padj<.05]
 expected=pd.concat([sig[sig.log2FoldChange>0].sort_values(['padj','gene']).head(30),sig[sig.log2FoldChange<0].sort_values(['padj','gene']).head(30)])
 selected=pd.read_csv(out/f'heatmap_{sex}_selected_genes.csv')
 assert set(selected.gene)==set(expected.gene)  # R/Python locale order may differ within exact ties.
 z=pd.read_csv(out/f'heatmap_{sex}_row_z_scores.csv',index_col=0)
 assert z.shape==(60,12) and np.allclose(z.mean(axis=1),0,atol=1e-12) and np.allclose(z.std(axis=1,ddof=1),1,atol=1e-12)
 checks[sex+'_heatmap_selection_and_zscores']=True
coords=pd.read_csv(out/'Fig1B_ordination_coordinates.csv');assert len(coords)==12 and (coords.groupby('group').size()==3).all()
s=pd.read_csv(A/'figure_sources/Figure_S3_Selected_Gene_Modules/Rebuilt_Output/FigS3_supporting_data.csv')
de=pd.read_csv(A/'data/rnaseq_female.csv');j=s.merge(de,on='gene',suffixes=('_plot','_main'))
assert len(s)==len(j)==89 and s.gene.nunique()==86
assert np.allclose(j.padj_plot,j.padj_main,rtol=1e-12,atol=0) and np.allclose(j.log2FoldChange_plot,j.log2FoldChange_main,rtol=1e-12,atol=0)
checks['S3_matches_main_female_contrast']=True
l=pd.read_csv(out/'Fig3H_statistics.csv').iloc[0];assert l.n_WT==4 and l.n_CS==4 and l.allocations==70 and abs(l.exact_rank_permutation_P-2/70)<1e-12
assert len(pd.read_csv(out/'Fig4A_plotted_values.csv'))==64
checks['locomotion_exact_test_and_TIMELESS_count']=True
runs=json.loads((A/'results/logs/full_panels_run.json').read_text());assert len(runs)==8 and all(x['exit_code']==0 for x in runs)
checks['original_generators_completed']=8
assert len(json.loads((R/'Review/every_figure_audit.json').read_text()))==36
(R/'Review/output_validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
