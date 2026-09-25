#!/usr/bin/env python3
"""Render a one-page vector PDF (or PDF-compatible AI file) at 600 dpi."""
from pathlib import Path
import argparse
import pypdfium2 as pdfium

def export_png(source,destination,dpi=600):
    with pdfium.PdfDocument(str(source)) as pdf:
        if len(pdf)!=1:raise ValueError('Expected one artboard/page')
        image=pdf[0].render(scale=dpi/72).to_pil().convert('RGB')
        image.save(destination,dpi=(dpi,dpi),optimize=True)
        return image.size

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source',type=Path);p.add_argument('destination',type=Path);p.add_argument('--dpi',type=int,default=600);a=p.parse_args()
    print(export_png(a.source,a.destination,a.dpi))
