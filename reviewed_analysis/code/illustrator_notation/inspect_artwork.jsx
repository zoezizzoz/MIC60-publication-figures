(function(){var lines=[];for(var i=0;i<app.documents.length;i++){var doc=app.documents[i],path='';try{path=doc.fullName.fsName;}catch(e){}if(path.indexOf('/dMIC60 Final Graphs and Figures/Assembled_Figures/')<0||!/^Fig[34]\.ai$/.test(doc.name))continue;
 lines.push('DOC '+path);
 for(var j=0;j<doc.placedItems.length;j++){var p=doc.placedItems[j];var file='';try{file=p.file.fsName;}catch(e){}lines.push(['PLACED',j,file,p.position,p.width,p.height].join('\t'));}
 for(var j=0;j<doc.rasterItems.length;j++){var p=doc.rasterItems[j];lines.push(['RASTER',j,p.name,p.position,p.width,p.height].join('\t'));}
}var f=new File('<PROJECT_WORKSPACE>/review_work/assembly_notation/artwork.tsv');f.open('w');f.write(lines.join('\n'));f.close();return 'Artwork inspected';})();
