from pathlib import Path
import json, shutil, subprocess
import pdfplumber
ROOT=Path(__file__).resolve().parents[1]
report={}
renderer=shutil.which('pdftoppm')
if not renderer: raise SystemExit('Install Poppler and put pdftoppm on PATH.')
for pdf in ROOT.glob('Sources/*/Final_Graphs/*_Vector.pdf'):
    with pdfplumber.open(pdf) as doc:
        page=doc.pages[0]
        assert len(page.images)==0, 'Raster image found in '+str(pdf)
        assert page.chars, 'No editable PDF text'
        cropped=[c['text'] for c in page.chars if c['x0']<-.5 or c['x1']>page.width+.5 or c['top']<-.5 or c['bottom']>page.height+.5]
        assert not cropped, ('Text outside page',cropped)
        sizes=sorted({round(c['size'],2) for c in page.chars if c['upright']})
        r={'page_width_pt':page.width,'page_height_pt':page.height,'image_objects':len(page.images),'vector_rectangles':len(page.rects),'vector_curves':len(page.curves),'text_characters':len(page.chars),'fonts':sorted({c['fontname'] for c in page.chars}),'upright_text_sizes_pt':sizes,'text_outside_page':cropped}
        if 'TIMELESS' in pdf.stem:
            boxes=[d for d in page.rects if isinstance(d.get('non_stroking_color'),(list,tuple)) and max(d['non_stroking_color'])-min(d['non_stroking_color'])>.2 and d['width']>15 and d['height']>5]
            boxes=list({tuple(round(d[k],3) for k in ('x0','top','x1','bottom')):d for d in boxes}.values())
            r['box_visible_widths_pt']=[round(d['width']+d['linewidth'],4) for d in boxes]
            assert len(boxes)==4
            assert all(abs(v-28)<.01 for v in r['box_visible_widths_pt'])
        report[pdf.stem]=r
    subprocess.run([renderer,'-scale-to','1100','-png','-singlefile',str(pdf),str(pdf.with_suffix(''))],check=True,capture_output=True,timeout=60)
(ROOT/'Documentation/PDF_VECTOR_QA.json').write_text(json.dumps(report,indent=2))
print(json.dumps(report,indent=2))
