(function(){
 var f=new File('<PROJECT_WORKSPACE>/review_work/assembly_notation/inventory.tsv');f.encoding='UTF-8';f.open('w');
 for(var d=0;d<app.documents.length;d++){
  var doc=app.documents[d],path='';try{path=doc.fullName.fsName;}catch(e){}
  f.writeln(['DOC',d,encodeURIComponent(doc.name),encodeURIComponent(path),doc.saved,doc.textFrames.length].join('\t'));
  if(!/^Fig(?:[1-4]|S[1-4])\.ai$/.test(doc.name))continue;
  for(var t=0;t<doc.textFrames.length;t++){
   var tf=doc.textFrames[t],s=tf.contents;if(!/MIC60|UAS|Day\s*\d|WT|CS|[♀♂]|p\s*[=<>]/i.test(s))continue;
   var runs=[];for(var j=0;j<tf.characters.length;j++){var a=tf.characters[j].characterAttributes;runs.push(encodeURIComponent(tf.characters[j].contents)+':'+a.textFont.name+':'+a.size);}
   f.writeln(['TEXT',d,t,encodeURIComponent(s),tf.position.join(','),tf.width,tf.height,tf.locked,tf.hidden,runs.join('|')].join('\t'));
  }
 }
 f.close();return 'Inventory complete';
})();
