#!/usr/bin/env python3
"""Run the selected recovered analysis or current panel workflow from this checkout."""
from pathlib import Path
import argparse,datetime,json,os,shutil,subprocess,sys
ROOT=Path(__file__).resolve().parent
A='reviewed_analysis';P='panels'
TASKS={
'recent-panels':[
 ('R',P+'/Fig3N/Code/generate_Fig3N.R'),
 ('PY',P+'/Fig4AB/Code/prepare_plot_data.py'),
 ('R',P+'/Fig4AB/Sources/Figure_4A_TIMELESS/Code/generate_timeless_vector.R'),
 ('R',P+'/Fig4AB/Sources/Figure_4B_MTT/Code/generate_mtt_vector.R'),
 ('PY',P+'/FigS4A/restyle_s4a.py')],
'current-layout':[
 ('R',P+'/layout_2026-09-23/GO_Rebuild/code/05_go_pathways.R'),
 ('R',P+'/layout_2026-09-23/Fig3_Sleep_Rebuild/code/02_sleep_activity.R')],
'networks':[
 ('R',A+'/figure_sources/Figure_2C_Program_Ring/Code/plot_expanded_paper_style.R'),
 ('R',A+'/figure_sources/Figure_2D_STRING_Network/Code/plot_compact_landscape.R')],
'reviewed':[('PY',A+'/run_all.py')]
}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--step',choices=TASKS);p.add_argument('--list',action='store_true');p.add_argument('--rscript',default='Rscript');a=p.parse_args()
 if a.list:
  for name,jobs in TASKS.items():
   print(name+':');[print('  '+path) for _,path in jobs]
  return
 if not a.step:p.error('Choose --step or --list.')
 r=shutil.which(a.rscript)
 if not r:p.error('Rscript was not found. Install R or pass --rscript.')
 out=ROOT/'verification/local_runs'/datetime.datetime.now().strftime('%Y%m%d_%H%M%S');out.mkdir(parents=True)
 env=os.environ.copy();env['RSCRIPT']=r;records=[]
 for i,(kind,rel) in enumerate(TASKS[a.step]):
  cmd=[r if kind=='R' else sys.executable,str(ROOT/rel)]
  if rel.endswith('/run_all.py'):cmd+=['--rscript',r]
  logfile=out/(str(i+1)+'_'+Path(rel).stem+'.log')
  with logfile.open('w') as log:result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
  records.append({'script':rel,'exit_code':result.returncode,'log':str(logfile.relative_to(ROOT))})
  (out/'results.json').write_text(json.dumps(records,indent=2)+'\n')
  print(rel+': '+('passed' if result.returncode==0 else 'FAILED; see '+str(logfile)))
  if result.returncode:raise SystemExit(result.returncode)
if __name__=='__main__':main()
