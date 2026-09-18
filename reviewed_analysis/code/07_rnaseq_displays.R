# Extend the count-matrix workflow in 01_rnaseq.R to explicitly reproducible displays.
# Original ordination/heatmap generators are unavailable. These are documented
# reconstructions, not assertions of identical original transformation settings.
suppressPackageStartupMessages({library(DESeq2);library(ggplot2);library(pheatmap);library(patchwork)})
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
root<-dirname(dirname(script));source(file.path(dirname(script),'figure_style.R'))
out<-file.path(root,'results/full_figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-as.matrix(read.csv(file.path(root,'data/rnaseq_counts.csv'),row.names=1,check.names=FALSE));storage.mode(x)<-'integer'
samples<-colnames(x);group<-factor(paste0(ifelse(grepl('^CS',samples),'CS','WT'),'_',ifelse(grepl('F[123]$',samples),'F','M')),levels=c('WT_F','CS_F','WT_M','CS_M'))
stopifnot(length(samples)==12,all(table(group)==3),!anyDuplicated(rownames(x)))
dds<-DESeqDataSetFromMatrix(x,data.frame(row.names=samples,group),design=~group)
# Design-blind rlog for sample exploration; it does not supply any DE P values.
rd<-rlog(dds,blind=TRUE);mat<-assay(rd)
write.csv(mat,file.path(out,'RNAseq_rlog_all_libraries.csv'))
label<-sub('^WR','WT',samples);label<-sub('([FCM])([123])$','\\1\\2',label)
label<-sub('^(WT|CS)(F|M)([123])$','\\1-\\2\\3',label)
variable<-apply(mat,1,var)>0
pc<-prcomp(t(mat[variable,,drop=FALSE]),center=TRUE,scale.=FALSE)
pct<-100*pc$sdev^2/sum(pc$sdev^2)
md<-cmdscale(dist(t(mat[variable,,drop=FALSE]),method='euclidean'),k=2,eig=TRUE)
coords<-data.frame(sample=samples,label,group,PC1=pc$x[,1],PC2=pc$x[,2],MDS1=md$points[,1],MDS2=md$points[,2])
write.csv(coords,file.path(out,'Fig1B_ordination_coordinates.csv'),row.names=FALSE)
ord<-function(ax,ay,lx,ly) ggplot(coords,aes(x=.data[[ax]],y=.data[[ay]],color=group))+
 geom_point(size=1.6)+ggrepel::geom_text_repel(aes(label=label),size=FIG_ANNOT_SIZE,family=FIG_FONT,seed=FIG_SEED,max.overlaps=Inf,box.padding=.35,show.legend=FALSE)+
 scale_color_manual(values=FIG_GROUP_COLORS,labels=FIG_SEX_GROUP_LABELS)+labs(x=lx,y=ly)+theme_fig(legend_position='top')+
 theme(legend.title=element_blank())+guides(color=guide_legend(nrow=2))
p<-(ord('PC1','PC2',sprintf('PC1 (%.1f%%)',pct[1]),sprintf('PC2 (%.1f%%)',pct[2])) |
 ord('MDS1','MDS2','MDS1','MDS2'))+plot_layout(guides='collect') & theme(legend.position='top')
fig_save(p,file.path(out,'Fig1B_PCA_MDS'),width=580/72,height=215/72)
zs<-list();selection<-list()
for(sex in c('female','male')) {
 de<-read.csv(file.path(root,'data',paste0('rnaseq_',sex,'.csv')))
 # Preserve the stated top-30-per-direction rule; rank by the archived sex-specific BH P.
 chosen<-do.call(rbind,lapply(c('up','down'),function(direction) {
  d<-de[!is.na(de$padj)&de$padj<.05 & if(direction=='up') de$log2FoldChange>=.58 else de$log2FoldChange<=-.58,]
  d<-d[order(d$padj,d$gene),];d<-head(d,30);d$direction<-direction;d
 }))
 stopifnot(nrow(chosen)==60,!anyDuplicated(chosen$gene),all(chosen$gene %in% rownames(mat)))
 z<-t(scale(t(mat[chosen$gene,,drop=FALSE])));stopifnot(all(is.finite(z)))
 zs[[sex]]<-z;selection[[sex]]<-chosen
 write.csv(chosen,file.path(out,paste0('heatmap_',sex,'_selected_genes.csv')),row.names=FALSE)
 write.csv(z,file.path(out,paste0('heatmap_',sex,'_row_z_scores.csv')))
}
lim<-ceiling(max(abs(unlist(zs))))
# pheatmap measures text on a temporary PostScript device. Register the
# standard sans metrics alias there; the final Quartz export embeds Arial.
grDevices::pdfFonts(Arial=grDevices::pdfFonts('Helvetica')[[1]])
annotation<-data.frame(Group=group,row.names=samples)
for(sex in names(zs)) {
 z<-zs[[sex]]
 grDevices::quartz(type='pdf',file=tempfile(fileext='.pdf'),width=330/72,height=540/72,family=FIG_FONT)
 hm<-pheatmap(z,color=colorRampPalette(unname(FIG_HEATMAP_SCALE))(101),
 breaks=seq(-lim,lim,length.out=102),cluster_rows=TRUE,cluster_cols=TRUE,clustering_distance_rows='euclidean',
 clustering_distance_cols='euclidean',clustering_method='complete',annotation_col=annotation,
 annotation_colors=list(Group=FIG_GROUP_COLORS),annotation_names_col=FALSE,annotation_legend=FALSE,
 labels_row=as.expression(lapply(rownames(z),function(g) bquote(italic(.(g))))),
 labels_col=label,fontsize=7,fontsize_row=7,fontsize_col=7,fontfamily=FIG_FONT,border_color=NA,
 treeheight_row=18,treeheight_col=18,angle_col='90',legend_breaks=c(-lim,0,lim),silent=TRUE)
 grDevices::dev.off()
 stem<-if(sex=='female')'Fig1D_Female_Heatmap' else 'FigS1D_Male_Heatmap'
 fig_save(hm$gtable,file.path(out,stem),width=330/72,height=540/72)
 write.csv(data.frame(row_order=rownames(z)[hm$tree_row$order]),file.path(out,paste0(stem,'_row_order.csv')),row.names=FALSE)
 write.csv(data.frame(sample_order=colnames(z)[hm$tree_col$order]),file.path(out,paste0(stem,'_column_order.csv')),row.names=FALSE)
}
write.csv(data.frame(group=names(FIG_GROUP_COLORS),label=FIG_SEX_GROUP_LABELS[names(FIG_GROUP_COLORS)],hex=unname(FIG_GROUP_COLORS)),file.path(out,'heatmap_group_key.csv'),row.names=FALSE)
jsonlite::write_json(list(transformation='DESeq2 rlog, blind=TRUE, all 12 libraries',PCA_genes=sum(variable),PCA_variance_percent=pct[1:2],
 MDS='classical MDS, Euclidean distances, same nonconstant rlog genes',heatmap_selection='30 smallest BH P values per direction among adjusted P < 0.05 and absolute log2FC >= 0.58; ties by gene name',
 scaling='row z scores over all 12 libraries',clustering='Euclidean, complete linkage',shared_z_limit=lim),file.path(out,'RNAseq_display_parameters.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,'RNAseq_displays_sessionInfo.txt'))

writeLines(capture.output(warnings()),file.path(out,"RNAseq_display_warnings.txt"))
