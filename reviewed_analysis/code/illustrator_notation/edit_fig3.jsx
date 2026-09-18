(function(){
 var original='<SOURCE_FIGURE_DIRECTORY>/Assembled_Figures/Fig3.ai';
 var doc=null;for(var i=0;i<app.documents.length;i++){try{if(app.documents[i].fullName.fsName==original)doc=app.documents[i];}catch(e){}}
 if(!doc)throw new Error('Current Fig3 not found');
 var opts=new IllustratorSaveOptions();opts.pdfCompatible=true;opts.compressed=true;
 var before=new File('<PROJECT_WORKSPACE>/review_work/assembly_notation/before/Fig3.ai');
 if(before.exists)throw new Error('Backup already exists; refusing to replace');
 doc.saveAs(before,opts);
 var regular=app.textFonts.getByName('ArialMT'),italic=app.textFonts.getByName('Arial-ItalicMT');
 var log=new File('<PROJECT_WORKSPACE>/review_work/assembly_notation/Fig3_changes.tsv');log.encoding='UTF-8';log.open('w');
 for(var t=0;t<doc.textFrames.length;t++){
  var tf=doc.textFrames[t],s=tf.contents,trim=s.replace(/^\s+|\s+$/g,''),newText=trim,kind='',italStart=0,italEnd=0,size=null;
  if(/^UAS-dMIC60(?:-(WT|CS)| \((WT|CS)\))?$/.test(trim)){
   kind='construct';newText=trim.replace(/-(WT|CS)$/,' ($1)');italEnd=10;size=7;
  }else if(/^(d?MIC60)-[Nn]ull(?:\s+(?:Flies|flies|HeLa [Cc]ells))?$/.test(trim)){
   kind='null';newText=trim.replace('-Null','-null').replace(' Flies',' flies').replace(' Cells',' cells');italEnd=newText.indexOf('-');size=7;
  }else if(/^Day\s+\d+$/.test(trim)) {kind='day';size=7;
  }else if(/^[♀♂]$/.test(trim)){kind='sex';size=7;
  }else if(/^(?:dMIC60-)?(?:WT|CS)$/.test(trim)){kind='protein';
  }else if(/^p\s*[<=>]/.test(trim)){kind='stat';italEnd=1;
  }
  if(!kind)continue;
  var cx=tf.left+tf.width/2,top=tf.top;var oldFonts=[];
  for(var c=0;c<tf.characters.length;c++)oldFonts.push(tf.characters[c].characterAttributes.textFont.name);
  tf.contents=newText;
  for(var c=0;c<tf.characters.length;c++){tf.characters[c].characterAttributes.textFont=(c>=italStart&&c<italEnd)?italic:regular;if(size!==null)tf.characters[c].characterAttributes.size=size;}
  if(kind!='stat')tf.left=cx-tf.width/2;tf.top=top;
  log.writeln([t,kind,encodeURIComponent(s),encodeURIComponent(newText),tf.left,tf.top,size,oldFonts.join(',')].join('\t'));
 }
 log.close();doc.saveAs(new File(original),opts);
 if(!new File(original).copy('<PROJECT_WORKSPACE>/review_output/MIC60_review/Figures/Illustrator/Fig3.ai'))throw new Error('Copy failed');
 return 'Fig3 labels updated and backup saved';
})();
