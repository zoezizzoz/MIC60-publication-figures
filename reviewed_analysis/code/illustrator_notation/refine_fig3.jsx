(function(){
 var path='<SOURCE_FIGURE_DIRECTORY>/Assembled_Figures/Fig3.ai',doc=null;
 for(var i=0;i<app.documents.length;i++){try{if(app.documents[i].fullName.fsName==path)doc=app.documents[i];}catch(e){}}
 if(!doc)throw new Error('Fig3 missing');
 // Preserve the left-aligned survival legend beside its colored line keys.
 for(var i=0;i<doc.textFrames.length;i++){
  var tf=doc.textFrames[i];
  if(/^UAS-dMIC60 \((WT|CS)\)$/.test(tf.contents)&&tf.top>505&&tf.top<525&&tf.left<160)tf.left=105.8310546875;
 }
 // Pair each existing sex symbol with its adjacent null label at the same height.
 for(var i=0;i<doc.textFrames.length;i++){
  var sym=doc.textFrames[i];if(sym.contents!='♀')continue;
  var best=null,dist=999;
  for(var j=0;j<doc.textFrames.length;j++){
   var tf=doc.textFrames[j];if(tf.contents!='dMIC60-null')continue;
   var dy=Math.abs((sym.top-sym.height/2)-(tf.top-tf.height/2));
   var dx=Math.abs(sym.left-(tf.left+tf.width));
   if(dy<6&&dx<20&&dx+dy<dist){best=tf;dist=dx+dy;}
  }
  if(best){sym.left=best.left+best.width+3;sym.top=best.top-best.height/2+sym.height/2;}
 }
 doc.save();new File(path).copy('<PROJECT_WORKSPACE>/review_output/MIC60_review/Figures/Illustrator/Fig3.ai');
 return 'Fig3 legend and symbol spacing checked';
})();
