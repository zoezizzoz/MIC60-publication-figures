(function(){
 var targets=[
 '<SOURCE_FIGURE_DIRECTORY>/Assembled_Figures/Fig4.ai',
 '<SOURCE_FIGURE_DIRECTORY>/Assembled_Figures/FigS4.ai'];
 var reg=app.textFonts.getByName('ArialMT'),ital=app.textFonts.getByName('Arial-ItalicMT');
 for(var n=0;n<targets.length;n++){
  var doc=null;for(var j=0;j<app.documents.length;j++){try{if(app.documents[j].fullName.fsName==targets[n])doc=app.documents[j];}catch(e){}}
  if(!doc)throw new Error('Document missing');
  if(doc.name=='Fig4.ai'){
   // This legacy MTT graph is an embedded raster. Cover only its old title;
   // retain the complete original artwork and every plotted value underneath.
   var group=doc.groupItems.add();group.name='Notation correction - embedded MTT title';
   var white=new RGBColor();white.red=white.green=white.blue=255;
   var cover=group.pathItems.rectangle(425,363,82,10);cover.filled=true;cover.fillColor=white;cover.stroked=false;
   var label=group.textFrames.add();label.contents='MIC60-null HeLa cells';
   label.textRange.characterAttributes.textFont=reg;label.textRange.characterAttributes.size=7;
   var black=new RGBColor();black.red=black.green=black.blue=0;label.textRange.characterAttributes.fillColor=black;
   for(var c=0;c<5;c++)label.characters[c].characterAttributes.textFont=ital;
   label.left=404-label.width/2;label.top=424;
  }else{
   for(var t=0;t<doc.textFrames.length;t++){var tf=doc.textFrames[t];if(tf.contents=='WT'){
    var cx=tf.left+tf.width/2,top=tf.top;tf.textRange.characterAttributes.size=7;tf.textRange.characterAttributes.textFont=reg;tf.left=cx-tf.width/2;tf.top=top;
   }}
  }
  doc.save();new File(targets[n]).copy('<PROJECT_WORKSPACE>/review_output/MIC60_review/Figures/Illustrator/'+doc.name);
 }
 return 'Embedded MTT title and paired blot label finalized';
})();
