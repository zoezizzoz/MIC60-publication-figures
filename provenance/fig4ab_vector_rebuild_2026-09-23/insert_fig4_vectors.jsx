(function(){function q(s){return '"'+String(s).replace(/\\/g,'\\\\').replace(/"/g,'\\"').replace(/\r/g,'\\r').replace(/\n/g,'\\n')+'"';}
function json(v){if(v===null)return'null';if(typeof v=='string')return q(v);if(typeof v=='number'||typeof v=='boolean')return String(v);if(v instanceof Array){var a=[];for(var k=0;k<v.length;k++)a.push(json(v[k]));return '['+a.join(',')+']';}var a=[];for(var k in v)if(v.hasOwnProperty(k))a.push(q(k)+':'+json(v[k]));return '{'+a.join(',')+'}';}
var root="/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Fig4AB_Vector_Rebuild_2026-09-23",master="/Users/picklejuice/Desktop/MIC60 Final Figs/Fig4.ai",jobs=[{"label": "A - TIMELESS descriptive, 64 observations", "bounds": [113.607145985127, 433.231751112585, 324.856893258626, 277.011258993047], "ai": "/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Fig4AB_Vector_Rebuild_2026-09-23/Sources/Figure_4A_TIMELESS/Final_Graphs/Fig4A_TIMELESS_Vector.ai"}, {"label": "B - MTT descriptive, 78 technical wells", "bounds": [350.857623753298, 432.773616640276, 545.768368884039, 276.183286524436], "ai": "/Users/picklejuice/Desktop/dMIC60 Graphs Code and Data/Fig4AB_Vector_Rebuild_2026-09-23/Sources/Figure_4B_MTT/Final_Graphs/Fig4B_MTT_Vector.ai"}],old=app.userInteractionLevel,src=null,d=null,created=[],removed=false,log=[];
function rect(x){return [x[0],x[1],x[2],x[3]];}
function at(path){for(var i=0;i<app.documents.length;i++){var d=app.documents[i];if(d.fullName.fsName===path)return d;}throw Error('Open document missing: '+path);}
function match(b,c){for(var k=0;k<4;k++)if(Math.abs(b[k]-c[k])>.01)return false;return true;}
function findRaster(d,b){var a=[];for(var k=0;k<d.rasterItems.length;k++){var r=d.rasterItems[k];if(match(r.geometricBounds,b))a.push(r);}if(a.length!==1)throw Error('Ambiguous/missing target raster '+b);return a[0];}
function snapshot(d){var a=[];for(var i=0;i<d.rasterItems.length;i++){var r=d.rasterItems[i];a.push({name:r.name,bounds:rect(r.geometricBounds)});}return a;}
function exportPNG(d,path){var x=new ExportOptionsPNG24();x.artBoardClipping=true;x.antiAliasing=true;x.transparency=false;x.horizontalScale=150;x.verticalScale=150;d.exportFile(new File(path),ExportType.PNG24,x);}
function checkpoint(msg){var f=new File(root+'/Documentation/insertion_progress.txt');f.encoding='UTF-8';f.open('a');f.writeln(new Date().toString()+' '+msg);f.close();}
try{app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;d=at(master);d.activate();
for(var i=0;i<jobs.length;i++)findRaster(d,jobs[i].bounds);
var ar=rect(d.artboards[0].artboardRect),allBefore=snapshot(d),scientificBefore=[];
for(var i=0;i<allBefore.length;i++)if(!match(allBefore[i].bounds,jobs[0].bounds)&&!match(allBefore[i].bounds,jobs[1].bounds))scientificBefore.push(allBefore[i]);
if(scientificBefore.length!==20)throw Error('Expected 20 scientific photographs');
var so=new IllustratorSaveOptions();so.pdfCompatible=true;so.compressed=true;
d.saveAs(new File(root+'/Documentation/Before_Fig4_Retry.ai'),so);d.saveAs(new File(master),so);exportPNG(d,root+'/Previews/Fig4_Before.png');checkpoint('Fresh live backup and preview saved');
for(var i=0;i<jobs.length;i++){var j=jobs[i];src=app.open(new File(j.ai));src.activate();if(src.rasterItems.length||src.placedItems.length||src.pluginItems.length)throw Error('Source must be wholly vector');
src.selection=null;app.executeMenuCommand('selectall');app.copy();src.close(SaveOptions.DONOTSAVECHANGES);src=null;
d=at(master);d.activate();var layer=d.layers.add();layer.name='Figure 4'+j.label;created.push(layer);var g=layer.groupItems.add();g.name='Fig4'+j.label+' - editable vectors';
d.selection=null;app.paste();var pasted=d.selection;if(!pasted||!pasted.length)throw Error('No vector artwork pasted');for(var k=pasted.length-1;k>=0;k--)pasted[k].move(g,ElementPlacement.PLACEATBEGINNING);
var b=rect(g.geometricBounds),w=j.bounds[2]-j.bounds[0],h=j.bounds[1]-j.bounds[3];if(Math.abs((b[2]-b[0])-w)>.1||Math.abs((b[1]-b[3])-h)>.1)throw Error('Imported dimensions differ: '+b);
g.translate(j.bounds[0]-b[0],j.bounds[1]-b[1]);g.note='Regenerated from original observations and latest recoverable analytical code; placed at 100%. Source package: '+root;
log.push({panel:j.label,source:j.ai,group:g.name,bounds:rect(g.geometricBounds),raster_items:layer.rasterItems.length,text_frames:layer.textFrames.length,path_items:layer.pathItems.length});
d=at(master);d.activate();checkpoint('Imported vector '+j.label);}
for(var i=0;i<jobs.length;i++)findRaster(d,jobs[i].bounds).remove();removed=true;
var after=snapshot(d);if(after.length!==20)throw Error('Unexpected remaining raster count');
for(var i=0;i<scientificBefore.length;i++){var found=false;for(var k=0;k<after.length;k++)if(after[k].name===scientificBefore[i].name&&match(after[k].bounds,scientificBefore[i].bounds))found=true;if(!found)throw Error('Scientific photograph changed');}
if(!match(rect(d.artboards[0].artboardRect),ar))throw Error('Artboard changed');
d.selection=null;exportPNG(d,root+'/Previews/Fig4_After.png');d.saveAs(new File(root+'/Illustrator_Figures/Fig4.ai'),so);d.saveAs(new File(master),so);
var report={master:d.fullName.fsName,saved:d.saved,package_copy:root+'/Illustrator_Figures/Fig4.ai',artboard:ar,remaining_rasters:after.length,scientific_photographs_preserved:true,panels:log};var f=new File(root+'/Documentation/NATIVE_INSERTION_QA.json');f.encoding='UTF-8';f.open('w');f.write(json(report));f.close();checkpoint('Updated master and package copy saved');return json(report);
}catch(e){if(!removed){try{if(d){d.activate();for(var k=created.length-1;k>=0;k--)created[k].remove();}}catch(ee){}}throw e;}finally{if(src)src.close(SaveOptions.DONOTSAVECHANGES);app.userInteractionLevel=old;}})()
