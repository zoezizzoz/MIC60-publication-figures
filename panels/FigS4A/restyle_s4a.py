from pathlib import Path
import csv,json,math,statistics,xml.etree.ElementTree as ET
P=Path(__file__).parent
rows=list(csv.DictReader((P/'figS4_plotted_values.csv').open()));enrich=list(csv.DictReader((P/'figS4_enrichment.csv').open()))
W,H=612,310;L,R,T,B=58,602,38,274;ymin,ymax=-1.7,10.7
X=lambda x:L+(float(x)-.4)/6.2*(R-L)
Y=lambda y:B-(float(y)-ymin)/(ymax-ymin)*(B-T)
cols={'higher_cs':'#B2182B','higher_wt':'#5AB4E5','not_significant':'#A6A6A6','not_tested':'#D1D1D1'}
def status(r):
 if r['padj_F'] in ('NA',''):return 'not_tested'
 p=float(r['padj_F']);v=float(r['lfc_F'])
 return 'higher_cs' if p<.05 and v>=.58 else 'higher_wt' if p<.05 and v<=-.58 else 'not_significant'
S='http://www.w3.org/2000/svg';ET.register_namespace('',S)
svg=ET.Element('{'+S+'}svg',{'width':'612pt','height':'310pt','viewBox':'0 0 612 310','version':'1.1'})
def el(tag,attrs,parent=svg):return ET.SubElement(parent,'{'+S+'}'+tag,{k:str(v) for k,v in attrs.items()})
def text(x,y,s,size=8.5,color='#000000',anchor='start',italic=False,parent=svg,**attrs):
 t=el('text',{'x':round(x,4),'y':round(y,4),'font-family':'Arial','font-size':size,'fill':color,'text-anchor':anchor,**({'font-style':'italic'} if italic else {}),**attrs},parent);t.text=s;return t
def span(t,s,**attrs):e=el('tspan',attrs,t);e.text=s;return e
def line(x1,y1,x2,y2,color='#000000',width=.55,**attrs):return el('line',{'x1':round(x1,5),'y1':round(y1,5),'x2':round(x2,5),'y2':round(y2,5),'stroke':color,'stroke-width':width,**attrs})
# Title and legend follow the supplied visual reference.
t=text(58,16,'',9);span(t,'dMIC60',**{'font-style':'italic'});span(t,'-null female')
x=215
for label,key,advance in [('Higher in CS','higher_cs',70),('Higher in WT','higher_wt',71),('Not significant','not_significant',85),('Not tested','not_tested',64)]:
 el('circle',{'cx':x,'cy':13,'r':2.6,'fill':cols[key]});text(x+6,16,label,8.5);x+=advance
# Threshold region and axes.
el('rect',{'x':L,'y':Y(.58),'width':R-L,'height':Y(-.58)-Y(.58),'fill':'#EEEEEE'})
line(L,Y(0),R,Y(0),'#666666',.65)
for v in [-.58,.58]:line(L,Y(v),R,Y(v),'#A6A6A6',.55,**{'stroke-dasharray':'2.5 2.5'})
line(L,T,L,B);line(L,B,R,B)
for v in [0,2.5,5,7.5,10]:
 y=Y(v);line(L-2.1,y,L,y);text(L-4.5,y+2.9,f'{v:.1f}',8.5,anchor='end')
t=text(0,0,'',9.5,anchor='middle',transform='translate(21 154) rotate(-90)');span(t,'log');span(t,'2',**{'baseline-shift':'sub','font-size':'6.65'});span(t,' fold change (');span(t,'dMIC60',**{'font-style':'italic'});span(t,'-CS/');span(t,'dMIC60',**{'font-style':'italic'});span(t,'-WT)')
pathways=['Control','ISR','UPR (ATF6)','UPR (IRE1/XBP1s)','HSR','OSR'];xt=[['Non-stress','control'],['ISR'],['ER-UPR','(ATF6)'],['ER-UPR','(IRE1/XBP1s)'],['HSR'],['OSR']]
for i,parts in enumerate(xt,1):
 for j,s in enumerate(parts):text(X(i),287+j*10.5,s,8.7,anchor='middle')
# Source statistics, unchanged; reference-like muted group annotations.
for i,pw in enumerate(pathways,1):
 e=next(e for e in enrich if e['pathway']==pw);odds='∞' if e['odds_ratio']=='Inf' else f"{float(e['odds_ratio']):.2f}"
 text(X(i),51,'OR = '+odds,8.2,'#666666','middle')
 t=text(X(i),65,'',8.2,'#666666','middle');span(t,'p',**{'font-style':'italic'});span(t,'BH',**{'baseline-shift':'sub','font-size':'5.75'});span(t,' = '+f"{float(e['fisher_padj']):.2f}")
# Sixty-four numeric observations, with the source table's recorded jitter.
for r in rows:
 key=status(r);el('circle',{'id':'point_'+r['pathway'].replace(' ','_')+'_'+r['gene'],'cx':X(r['x_plot']),'cy':Y(r['lfc_F']),'r':1.1,'fill':cols[key]})
for i,pw in enumerate(pathways,1):
 m=statistics.median(float(r['lfc_F']) for r in rows if r['pathway']==pw);line(X(i-.175),Y(m),X(i+.175),Y(m),'#4D4D4D',1.3)
# Deliberate placements keep all seventeen selected gene labels legible.
placements={
 ('Control','mtSSB'):(.59,1.15),('Control','betaGlu'):(1.15,.50),('Control','mRpL9'):(.57,-.97),('Control','Prosbeta6'):(1.21,-.97),
 ('ISR','PPP1R15'):(1.73,1.15),('ISR','Gadd45'):(2.16,.28),('ISR','crc'):(1.88,-.97),
 ('UPR (ATF6)','Mis12'):(3.15,2.10),('UPR (ATF6)','Hsc70-3'):(3.00,-.77),
 ('UPR (IRE1/XBP1s)','Hsc70-4'):(3.53,.89),('UPR (IRE1/XBP1s)','l(1)G0320'):(4.08,-1.08),
 ('HSR','Hsc70-4'):(5.03,.87),('HSR','Hsp83'):(4.54,1.43),('HSR','Hsp60A'):(5.00,-1.06),
 ('OSR','Gclm'):(6.23,.88),('OSR','Gclc'):(5.72,-.82),('OSR','Sod1'):(6.24,-.92)}
label_plan=[]
for r in rows:
 key=(r['pathway'],r['gene'])
 if key not in placements:continue
 xx,yy=placements[key];x,y=X(xx),Y(yy);px,py=X(r['x_plot']),Y(r['lfc_F']);font=7.9;estimated_width=sum(.30 if c in 'ilt1()' else .57 for c in r['gene'])*font
 tx=min(max(px,x),x+estimated_width);ty=min(max(py,y-font*.75),y+.7);dx,dy=tx-px,ty-py;n=math.hypot(dx,dy)
 if n>3:line(px+dx*2/n,py+dy*2/n,tx-dx/n,ty-dy/n,'#999999',.43)
 color=cols[status(r)]
 # A white underlay avoids collisions with the threshold dashes.
 el('rect',{'x':x-1,'y':y-font*.82,'width':estimated_width+2,'height':font*1.05,'fill':'#FFFFFF','fill-opacity':.96})
 text(x,y,r['gene'],font,color,italic=True)
 if r['gene']=='Mis12':text(x,y+9.5,'adj. p = 0.01068',6.4,color)
 if r['gene']=='l(1)G0320':text(x,y+8.2,'adj. p = 0.001212',6.4,color)
 label_plan.append({'pathway':r['pathway'],'gene':r['gene'],'color':color,'point_xy':[px,py],'label_xy':[x,y]})
assert len(label_plan)==17 and len(rows)==64
ET.ElementTree(svg).write(P/'S4A_reference_style.svg',encoding='utf-8',xml_declaration=True)
(P/'RESTYLE_MANIFEST.json').write_text(json.dumps({'panel_size_pt':[W,H],'point_count':len(rows),'gene_label_count':len(label_plan),'colors':cols,'labels':label_plan,'data_source':'figS4_plotted_values.csv','enrichment_source':'figS4_enrichment.csv','jitter':'Recorded x_plot values from the recovered source table','retained_gene_p_annotations':True},indent=2))
print('Built editable vector S4A: 64 points, 17 gene labels, 2 gene P labels, 6 existing pathway summaries')
