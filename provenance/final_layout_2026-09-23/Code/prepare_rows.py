from pathlib import Path
import json,copy,math
R=Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23/Layout_Data'); O=Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23')
OLD=json.loads(Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23/Layout_Data/CATEGORY_PATH_ROLES.json').read_text())
def cx(b):return (b[0]+b[2])/2
def cy(b):return (b[1]+b[3])/2
def width(b):return b[2]-b[0]
def height(b):return b[1]-b[3]
def newdoc(name):
 d=json.loads((R/(name+'_NORMALIZED.json')).read_text());p={'name':name+'.ai','path':d['path'],'paths':{},'texts':{},'items':[],'layout':[]};return d,p

def path(p,q,fn,role):
 pts=[[[*fn(v)] for v in trip] for trip in q['points']]
 if max((abs(a-b) for u,v in zip(pts,q['points']) for aa,bb in zip(u,v) for a,b in zip(aa,bb)),default=0)<1e-8:return
 if q['id'] in p['paths']:raise ValueError(('duplicatepath',q['id']))
 p['paths'][q['id']]={'id':q['id'],'before':q['b'],'after':pts,'role':role,'visible':q['visible']}
def txt(p,t,dx=0,dy=0,**kw):
 assert t['id'] not in p['texts'],('duplicatetext',t['id'])
 p['texts'][t['id']]={'id':t['id'],'before':t['b'],'original':t['text'],'dx':dx,'dy':dy,**kw}
def item(p,d,id,dx=0,dy=0,**kw):
 t=next(x for x in d['top'] if x['id']==id);p['items'].append({'id':id,'before':t['b'],'dx':dx,'dy':dy,**kw})
def letter(p,d,id,left,top=None):
 t=next(t for t in d['texts'] if t['id']==id);txt(p,t,dx=left-t['b'][0],dy=0 if top is None else top-t['b'][1])
def boxplot(d,p,rt,newlimits=None,compact=False,single=False,bgid=None):
 root=next(r for r in d['roots'] if r['id']==rt);ps={q['id']:q for q in root['paths']};ts={q['id']:q for q in d['texts']};old=next(q for q in OLD if q['doc']==d['name'] and q['rootId']==rt)
 bg=ps[bgid or str(int(old['groupId'])+1)];b=bg['b'];lo,hi=b[0],b[2];nl,nh=newlimits or (lo,hi)
 def fx(x):return nl+(x-lo)*(nh-nl)/(hi-lo)
 boxes=sorted([ps[q['id']] for q in old['boxes']],key=lambda q:cx(q['b']));centers=[cx(q['b']) for q in boxes];targets=[]
 if compact:
  for i in range(0,len(boxes),2):m=sum(centers[i:i+2])/2;targets.extend([m-18.5,m+18.5])
 else:
  fractions=[.22,.78] if len(boxes)==2 else [.14,.36,.64,.86];targets=[nl+(nh-nl)*f for f in fractions]
 pairs=[(sum(centers[i:i+2])/2,sum(targets[i:i+2])/2) for i in range(0,len(centers),2)]
 movers={q['id']:q for q in old['pathMoves']};boxw=width(boxes[0]['b'])
 def nearest(x):return min(range(len(centers)),key=lambda i:abs(x-centers[i]))
 for q in root['paths']:
  if not q['visible']:continue
  bb=q['b'];w=width(bb);h=height(bb);x=cx(bb);y=cy(bb);m=movers.get(q['id'])
  if m:
   i=nearest(x);dx=targets[i]-centers[i]
   if 'xmap' in m:
    k=min(range(0,len(centers),2),key=lambda k:abs(x-(centers[k]+centers[k+1])/2))
    fn=lambda v,k=k:(targets[k]+(v[0]-centers[k])*(targets[k+1]-targets[k])/(centers[k+1]-centers[k]),v[1]);role='significance bracket'
   elif len(q['points'])>=4 and q['fill'] is not None and w<5 and h<5:
    fn=lambda v,x=x,y=y,w=w,h=h,dx=dx:(x+dx+(v[0]-x)*2.4/w,y+(v[1]-y)*2.4/h);role='point (same data center)'
   elif abs(w-boxw)<.01 and (abs(x-centers[i])<.01):
    fn=lambda v,x=x,w=w,dx=dx:(x+dx+(v[0]-x)*24/w,v[1]);role='24 pt box/median'
   elif 'whisker cap' in q['name']:
    fn=lambda v,x=x,w=w,dx=dx:(x+dx+(v[0]-x)*16/w,v[1]);role='16 pt whisker cap'
   else:fn=lambda v,dx=dx:(v[0]+dx,v[1]);role=m['kind']
   path(p,q,fn,role)
  elif newlimits:
   # Keep the axis/tick gutter physical. Extend plot and clipping boundaries.
   fn=lambda v:(v[0]+nl-lo if v[0]<lo else v[0]+nh-hi if v[0]>hi else fx(v[0]),v[1])
   path(p,q,fn,'plot/axis boundary')
 oldtexts={q['id']:q for q in old['textMoves']}
 for t in d['texts']:
  if rt not in t['chain'] and t['id'] not in oldtexts:continue
  x=cx(t['b']);s=t['text'];m=oldtexts.get(t['id']);kw={};dx=0
  if m and m['kind']=='category label':
   i=nearest(x);dx=targets[i]-x
   if single:kw={'content':s.replace('\r','').replace('\n',''),'targetCenter':targets[i],'targetTop':b[3]-5}
  elif m:
   pair=min(pairs,key=lambda q:abs(x-q[0]));dx=pair[1]-pair[0]
  elif newlimits:
   dx=nl-lo if x<lo else fx(x)-x
  if dx or kw:txt(p,t,dx=dx,**kw)
 p['layout'].append({'root':rt,'plot_before':[lo,hi],'plot_after':[nl,nh],'centers_before':centers,'centers_after':targets,'box_width':24,'cap_width':16,'point_diameter':2.4,'data_y_preserved':True})
 return (lo,hi,nl,nh)
plans=[]
# Figure3: narrower bodies within compact plots; wide TEM and TMRM rows.
d,p=newdoc('Fig3');a=d['artboard'];off=a[0]
for rt in ['1843','525','2078','1048','1464']:boxplot(d,p,rt,compact=True)
for rt,le in zip(['1534','1645','1753'],[55,247,439]):
 r=boxplot(d,p,rt,(off+le,off+le+155),single=True)
 if rt=='1645':
  t=next(t for t in d['texts'] if t['id']=='7408');txt(p,t,dx=r[2]-r[0])
boxplot(d,p,'7293',(off+390,off+594),single=True)
for id,x in [('7427',18),('7430',18),('7431',210),('7432',402),('7433',18),('7434',346)]:letter(p,d,id,off+x)
# M moves as one group, including its images, labels and scale bars.
item(p,d,'7263',dx=-32)
# Align the right edge of the TEM row; retain image proportions and attached scales.
iroot=next(q for q in d['top'] if q['id']=='2306');dx=off+594-iroot['b'][2]
for id in ['2308','2306','2304','2313','2312','2311','2310','7429']:item(p,d,id,dx=dx)
plans.append(p)
# Figure2D: extend three compact plots to the common right edge.
d,p=newdoc('Fig2');ps={q['id']:q for r in d['roots'] for q in r['paths']};ts={q['id']:q for q in d['texts']};prev=json.loads(Path('/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Common_Width_2026-09-23/Layout_Data/FIG2D_PATH_ROLES.json').read_text());oldlimits=[ps[i]['b'] for i in ['831','889','951']];newlefts=[58,249,440];newwidth=154
cps=[next(z for z in OLD if z['doc']==d['name'] and z['groupId']==rt) for rt in ['830','888','950']];cs=[[cx(ps[z['id']]['b']) for z in sorted(q['boxes'],key=lambda z:cx(ps[z['id']]['b']))] for q in cps];nc=[[l+newwidth*.22,l+newwidth*.78] for l in newlefts]
for z in prev['paths']:
 if z['panel'] is None:continue
 q=ps[z['id']];k=z['panel'];lo,hi=oldlimits[k][0],oldlimits[k][2];nl=newlefts[k];nh=nl+newwidth;x=cx(q['b']);w=width(q['b']);role=z['role'];j=min(range(2),key=lambda j:abs(x-cs[k][j]));dx=nc[k][j]-cs[k][j]
 if role in ['header','plot-boundary','x-axis']:fn=lambda v,lo=lo,hi=hi,nl=nl,nh=nh:(nl+(v[0]-lo)*(nh-nl)/(hi-lo),v[1])
 elif role in ['y-axis-tick','y-axis']:fn=lambda v,nl=nl,lo=lo:(v[0]+nl-lo,v[1])
 elif role=='significance-bracket':fn=lambda v,k=k:(nc[k][0]+(v[0]-cs[k][0])*(nc[k][1]-nc[k][0])/(cs[k][1]-cs[k][0]),v[1])
 elif role=='box-whisker-median':
  nw=24 if abs(w-24.2243)<.02 else 16 if abs(w-16.1495)<.02 else w
  fn=lambda v,x=x,dx=dx,nw=nw,w=w:(x+dx+(v[0]-x)*nw/w if w else v[0]+dx,v[1])
 elif role=='data-point':
  y=cy(q['b']);h=height(q['b']);fn=lambda v,x=x,y=y,w=w,h=h,dx=dx:(x+dx+(v[0]-x)*2.4/w,y+(v[1]-y)*2.4/h)
 else:fn=lambda v,dx=dx:(v[0]+dx,v[1])
 path(p,q,fn,role)
for z in prev['texts']:
 t=ts[z['id']];role=z['role'];x=cx(t['b']);k=min(range(3),key=lambda k:abs(x-cx(oldlimits[k])))
 if role=='panel-letter':txt(p,t,dx=18-t['b'][0]);continue
 if role in ['numeric-tick-label','centered-axis-title']:dx=newlefts[k]-oldlimits[k][0]
 elif role=='single-line-category-label':j=0 if t['text'].endswith('WT') else 1;dx=nc[k][j]-x
 else:dx=newlefts[k]+newwidth/2-x
 txt(p,t,dx=dx)
p['layout'].append({'Figure2D':{'lefts':newlefts,'width':newwidth,'right':594,'centers':nc}});plans.append(p)
# Figure4 TIMELESS and MTT plots.
d,p=newdoc('Fig4');boxplot(d,p,'383',(58,300),bgid='602')
r=next(r for r in d['roots'] if r['id']=='608');ps={q['id']:q for q in r['paths']};b=ps['891']['b'];lo,hi=b[0],b[2];nl,nh=395,594
fx=lambda x:nl+(x-lo)*(nh-nl)/(hi-lo)
for q in r['paths']:
 if not q['visible']:continue
 x=cx(q['b']);w=width(q['b']);h=height(q['b'])
 if w<22 and h<13 and len(q['points'])>=3:fn=lambda v,x=x:(v[0]+fx(x)-x,v[1]);role='marker retains aspect'
 elif w<22 and h<.005 and lo<x<hi:fn=lambda v,x=x:(v[0]+fx(x)-x,v[1]);role='error cap/legend'
 else:fn=lambda v:(v[0]+nl-lo if v[0]<lo else v[0]+nh-hi if v[0]>hi else fx(v[0]),v[1]);role='continuous x-coordinate mapping'
 path(p,q,fn,role)
for t in d['texts']:
 if '608' not in t['chain']:continue
 x=cx(t['b'])
 if t['id'] in ['654','655','656','657','658','659']:dx=(nl+nh-lo-hi)/2
 elif t['id'] in ['650','651','687','688','689','690','691','692']:dx=nl-lo
 else:dx=fx(x)-x
 txt(p,t,dx=dx)
letter(p,d,'41',nl-21)
p['layout'].append({'MTT_plot_before':[lo,hi],'MTT_plot_after':[nl,nh],'data_y_preserved':True})
# Enlarge the two-part explanatory model proportionally, leaving physical font sizes unchanged.
# Its source canvas/clip paths are included, to keep all vector components registered.
ids=[q['id'] for q in d['top'] if (44<=int(q['id'])<=269)]+['43']
modeltexts=[t for t in d['texts'] if any(i in t['chain'] for i in ids)]
src=[92,110.5];scale=1.23;dst=[18,110.5]
for id in ids:item(p,d,id,uniform={'scale':scale,'source':src,'target':dst})
p['restoreTexts']=[{'id':t['id'],'before':t['b'],'uniform':{'scale':scale,'source':src,'target':dst}} for t in modeltexts]
p['artboardBottom']=-204
plans.append(p)
# S4: balance validation blot and a wide C plot beneath A.
d,p=newdoc('FigS4');boxplot(d,p,'436',(410,594),single=True)
letter(p,d,'526',370,260)
rt='1014';src=[39.45,234.42];scale=1.45;dst=[18,234.42]
item(p,d,rt,uniform={'scale':scale,'source':src,'target':dst})
p['restoreTexts']=[{'id':t['id'],'before':t['b'],'uniform':{'scale':scale,'source':src,'target':dst}} for t in d['texts'] if rt in t['chain']]
p['postTexts']=[{'id':'517','targetLeft':18,'targetTop':260}]
plans.append(p)
for p in plans:
 p['paths']=list(p['paths'].values());p['texts']=list(p['texts'].values())
 (R/(p['name'].replace('.ai','')+'_ROW_PLAN.json')).write_text(json.dumps(p,indent=2));print(p['name'],len(p['paths']),len(p['texts']),len(p['items']))
