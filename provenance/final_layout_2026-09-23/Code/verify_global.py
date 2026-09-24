from pathlib import Path
import pypdfium2 as pdf,ctypes,math,json,hashlib,sys,collections
R=Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23/Layout_Data');O=Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23');plans={p['name']:p for p in json.loads((R/'PDF_CONTENT_BOUNDS.json').read_text())};out=[]
def read(f):
 d=pdf.PdfDocument(str(f));p=d[0];tp=p.get_textpage();chars=[];images=[]
 for i in range(tp.count_chars()):
  u=pdf.raw.FPDFText_GetUnicode(tp,i)
  if not u or chr(u).isspace():continue
  m=pdf.raw.FS_MATRIX();pdf.raw.FPDFText_GetMatrix(tp,i,ctypes.byref(m));s=pdf.raw.FPDFText_GetFontSize(tp,i)*math.hypot(m.c,m.d);chars.append((chr(u),round(s,2)))
 for o in p.get_objects(filter=[pdf.raw.FPDF_PAGEOBJ_IMAGE]):
  im=o.get_bitmap().to_pil();images.append((im.size,hashlib.sha256(im.tobytes()).hexdigest()))
 im=p.render(scale=2).to_pil().convert('RGB');b=im.convert('L').point(lambda v:255 if v<245 else 0).getbbox();return d,p,collections.Counter(chars),sorted(images),im,[v/2 for v in b]
for name in sys.argv[1:] or list(plans):
 b,bp,bc,bi,bim,bb=read(O/'Before_AI'/name);a,ap,ac,ai,aim,ab=read(O/'After_AI'/name);assert abs(ap.get_width()-612)<.001,(name,ap.get_width());expected=bc.copy();
 if name=='Fig3.ai':expected[('2',7.0)]+=1
 assert expected==ac,(name,'text/font-size differences',expected-ac,ac-expected);assert bi==ai,(name,'raster pixel changes');aim.save(O/'Previews'/name.replace('.ai','_normalized.png'));row={'name':name,'canvas_pt':[ap.get_width(),ap.get_height()],'paint_bounds_pdf':ab,'original_text_and_font_sizes_preserved':True,'intentional_text_additions':['Missing area-axis tick 2 restored at 7 pt'] if name=='Fig3.ai' else [],'embedded_images_identical':True,'after_sha256':hashlib.sha256((O/'After_AI'/name).read_bytes()).hexdigest()};out.append(row);print(name,'PASS',row['canvas_pt'],ab,flush=True);a.close();b.close()
(R/'GLOBAL_VERIFICATION.json').write_text(json.dumps(out,indent=2))
