#!/usr/bin/env python3
"""Run the reviewed analyses in dependency order, preserving inputs and source artwork."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--rscript', default='Rscript')
parser.add_argument('--step', choices=['rnaseq','sleep','assays','panels','go','full_panels','additional_assays','rnaseq_displays','schematics'],
                    help='Run one step; panels require the assay outputs.')
parser.add_argument('--publish', action='store_true',
                    help='Copy generated panels into ../Figures and replace its manifest. Omit to preserve Illustrator-refined final exports.')
args = parser.parse_args()
rscript = shutil.which(args.rscript)
if not rscript:
    parser.error('Rscript not found; specify --rscript /path/to/Rscript')
env = os.environ.copy()
env['RSCRIPT'] = rscript
steps = [
    ('rnaseq',[rscript,str(ROOT/'code/01_rnaseq.R')]),
    ('sleep',[rscript,str(ROOT/'code/02_sleep_activity.R')]),
    ('assays',[sys.executable,str(ROOT/'code/03_assay_audit.py')]),
    ('panels',[rscript,str(ROOT/'code/04_review_panels.R')]),
    ('go',[rscript,str(ROOT/'code/05_go_pathways.R')]),
    ('additional_assays',[rscript,str(ROOT/'code/06_additional_assays.R')]),
    ('rnaseq_displays',[rscript,str(ROOT/'code/07_rnaseq_displays.R')]),
    ('full_panels',[sys.executable,str(ROOT/'code/08_full_figure_panels.py')]),
    ('schematics',[rscript,str(ROOT/'code/09_schematics.R')]),
]
if args.step:
    steps = [s for s in steps if s[0]==args.step]
logs=ROOT/'results/logs';logs.mkdir(parents=True,exist_ok=True)
inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (ROOT/'data').rglob('*') if p.is_file()}
records=[]
for name,command in steps:
    started=time.monotonic()
    with (logs/(name+'.log')).open('w') as log:
        result=subprocess.run(command,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
    records.append({'step':name,'return_code':result.returncode,'elapsed_seconds':round(time.monotonic()-started,2)})
    (logs/'run_results.json').write_text(json.dumps(records,indent=2)+'\n')
    print(f'{name}: '+('passed' if result.returncode==0 else f'FAILED; see {logs/name}.log'),flush=True)
    if result.returncode:
        sys.exit(result.returncode)
for rel,digest in inputs.items():
    if hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=digest:
        raise RuntimeError(f'Input was modified: {rel}')
(logs/'input_integrity.json').write_text(json.dumps({'unchanged':True,'files_checked':len(inputs)},indent=2)+'\n')
if args.publish:
    from publish_figures import publish
    publish()
