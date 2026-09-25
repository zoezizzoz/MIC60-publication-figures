#!/usr/bin/env python3
"""Run the edited original generators that have sufficient accessible inputs."""
from pathlib import Path
import subprocess,os,json,time,sys
root=Path(__file__).resolve().parents[1]
rscript=os.environ.get('RSCRIPT','Rscript')
scripts={
 'Figure_1C_Female_Volcano':'generate_consistent_female_male_volcano_plots.R',
 'Figure_2E_Individual_Genes':'generate_figure2E.R',
 'Figure_3K_Mitochondrial_Area':'generate_tem_figures_2026_08_04.R',
 'Figure_S1B_Western_Blot':'generate_MIC60_WB_anti_myc_quantification.R',
 'Figure_S2_KEGG_GSEA':'plot_FigS2_KEGG_GSEA.R',
 'Figure_S3_Selected_Gene_Modules':'rebuild_FigS3.py',
 'Figure_S4_Stress':'generate_FigS4A_stress.R',
 'Figure_S4B_Transfection_Efficiency':'generate_transfection_efficiency.R',
}
records=[]
for folder,name in scripts.items():
 p=root/'figure_sources'/folder/'Code'/name;t=time.monotonic()
 command=([sys.executable,str(p),'--rscript',rscript] if p.suffix=='.py' else [rscript,str(p)])
 with (root/'results/logs'/f'{folder}.log').open('w') as log:
  result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 records.append({'folder':folder,'script':str(p.relative_to(root)),'exit_code':result.returncode,'elapsed_seconds':round(time.monotonic()-t,2)})
 (root/'results/logs/full_panels_run.json').write_text(json.dumps(records,indent=2)+'\n')
 print(folder,':',result.returncode,flush=True)
 if result.returncode:raise SystemExit(result.returncode)
