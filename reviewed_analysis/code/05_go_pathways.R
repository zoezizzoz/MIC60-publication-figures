# Reproduce GO over-representation analyses and render legible standalone panels.
suppressPackageStartupMessages({library(clusterProfiler);library(org.Dm.eg.db);library(ggplot2);library(patchwork)})
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
root <- dirname(dirname(script)); out <- file.path(root,'results/go')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(dirname(script),'figure_style.R'))
source(file.path(dirname(script),'go_analysis_config.R'))
source(file.path(dirname(script),'go_highlight_config.R'))
highlight_audit <- list()
palette <- FIG_GO_COLORS
all_results <- list(); selections <- list(); plots <- list(); audit <- list(); y_spans <- list()
for (sex in GO_CONFIG$sexes) {
  d <- read.csv(file.path(root,'data',paste0('rnaseq_',sex,'.csv')))
  d <- d[!is.na(d$padj),]; stopifnot(!anyDuplicated(d$gene))
  for (direction in GO_CONFIG$directions) {
    genes <- d$gene[d$padj < GO_CONFIG$de_adjusted_p_cutoff & if(direction=='up')
      d$log2FoldChange >= GO_CONFIG$absolute_log2_fold_change_cutoff else
      d$log2FoldChange <= -GO_CONFIG$absolute_log2_fold_change_cutoff]
    fit <- enrichGO(gene=genes,OrgDb=org.Dm.eg.db,keyType=GO_CONFIG$key_type,ont='ALL',universe=d$gene,
      pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod=GO_CONFIG$p_adjust_method,
      minGSSize=GO_CONFIG$minimum_gene_set_size,maxGSSize=GO_CONFIG$maximum_gene_set_size,
      pool=GO_CONFIG$pool_ontologies)
    result <- fit@result
    stopifnot(all(result$p.adjust>=0 & result$p.adjust<=1),!anyDuplicated(result$ID),
      all(lengths(strsplit(result$geneID,'/',fixed=TRUE))==result$Count))
    result$sex <- sex; result$direction <- direction
    key <- paste(sex,direction,sep='_'); all_results[[key]] <- result
    write.csv(result,file.path(out,paste0(key,'_all_results.csv')),row.names=FALSE)
    selected <- do.call(rbind,lapply(GO_CONFIG$ontologies,function(ont) {
      a<-result[result$ONTOLOGY==ont & result$p.adjust<GO_CONFIG$go_adjusted_p_cutoff,]
      a<-a[do.call(order,a[GO_CONFIG$selection_order]),];head(a,GO_CONFIG$top_terms_per_ontology)
    }))
    selections[[key]] <- selected
    write.csv(selected,file.path(out,paste0(key,'_displayed_terms.csv')),row.names=FALSE)
    audit[[key]] <- data.frame(sex,direction,testable_genes=nrow(d),DEG_input=length(genes),
      tested_terms=nrow(result),BH_significant_terms=sum(result$p.adjust<GO_CONFIG$go_adjusted_p_cutoff),displayed_terms=nrow(selected))
    if(sex=='female' && direction=='up' && file.exists(file.path(root,'data/go_female_up_archived.csv'))) {
      archived<-read.csv(file.path(root,'data/go_female_up_archived.csv'))
      cmp<-merge(result,archived,by='ID',suffixes=c('_rerun','_archived'))
      check<-data.frame(rerun_terms=nrow(result),archived_terms=nrow(archived),matched_terms=nrow(cmp),
        max_adjusted_P_difference=max(abs(cmp$p.adjust_rerun-cmp$p.adjust_archived)),
        count_mismatches=sum(cmp$Count_rerun!=cmp$Count_archived))
      stopifnot(nrow(cmp)==nrow(result),nrow(cmp)==nrow(archived),check$max_adjusted_P_difference<1e-12,check$count_mismatches==0)
      write.csv(check,file.path(out,'female_up_archive_validation.csv'),row.names=FALSE)
    }
    # Reference layout: labels share the bar area, with a separate count column.
    # Full gene memberships remain in the companion CSV, at readable table scale.
    selected$ONTOLOGY <- factor(selected$ONTOLOGY,levels=GO_CONFIG$ontologies)
    panel_width <- FIG_PANEL_PT[[paste0('GO_',sex)]][1] / if(sex=='female') 2 else 1
    label_width <- panel_width - 50
    selected$term_label <- vapply(selected$Description,fig_go_wrap,character(1),width_pt=label_width)
    line_n <- lengths(strsplit(selected$term_label,'\n',fixed=TRUE))
    highlight <- GO_HIGHLIGHTS[GO_HIGHLIGHTS$sex==sex & GO_HIGHLIGHTS$direction==direction,]
    stopifnot(all(highlight$ID %in% selected$ID))
    selected$gene_label <- ''
    for (h in seq_len(nrow(highlight))) {
      idx <- match(highlight$ID[h],selected$ID)
      genes_to_show <- strsplit(highlight$genes[h],'/',fixed=TRUE)[[1]]
      members <- strsplit(selected$geneID[idx],'/',fixed=TRUE)[[1]]
      stopifnot(length(genes_to_show)>=2,length(genes_to_show)<=3,
                !anyDuplicated(genes_to_show),all(genes_to_show %in% members))
      selected$gene_label[idx] <- paste(genes_to_show,collapse=', ')
      stat <- d[match(genes_to_show,d$gene),c('gene','log2FoldChange','padj')]
      stopifnot(all(stat$padj<GO_CONFIG$de_adjusted_p_cutoff),
        all(if(direction=='up') stat$log2FoldChange>=GO_CONFIG$absolute_log2_fold_change_cutoff
            else stat$log2FoldChange<=-GO_CONFIG$absolute_log2_fold_change_cutoff))
      highlight_audit[[paste(key,highlight$ID[h])]] <- cbind(
        data.frame(sex=sex,direction=direction,ID=highlight$ID[h],
          Description=selected$Description[idx],reason=highlight$reason[h]),stat)
    }
    row_height <- pmax(12,line_n*FIG_GO_LINE_PT+FIG_GO_ROW_GAP_PT) +
      ifelse(nzchar(selected$gene_label),FIG_GO_GENE_LINE_PT,0)
    # Include an extra blank interval when the ontology changes.
    group_gap <- c(0,ifelse(diff(as.integer(selected$ONTOLOGY))!=0,4,0))
    cursor <- sum(row_height)+sum(group_gap)+3
    selected$y_top <- selected$y_bottom <- selected$y <- NA_real_
    for(i in seq_len(nrow(selected))) {
      cursor <- cursor-group_gap[i]
      selected$y_top[i] <- cursor
      selected$y_bottom[i] <- cursor-row_height[i]
      selected$y[i] <- cursor-FIG_GO_LINE_PT/2
      cursor <- cursor-row_height[i]
    }
    # Each text line uses a centered anchor. The first line shares exactly
    # the bar/count center; wrapped continuation lines remain below it.
    label_lines <- do.call(rbind,lapply(seq_len(nrow(selected)),function(i) {
      lines <- strsplit(selected$term_label[i],"\n",fixed=TRUE)[[1]]
      data.frame(line_label=lines,y=selected$y[i]-(seq_along(lines)-1)*FIG_GO_LINE_PT)
    }))
    gene_lines <- selected[nzchar(selected$gene_label),]
    gene_lines$gene_y <- selected$y[nzchar(selected$gene_label)] -
      line_n[nzchar(selected$gene_label)]*FIG_GO_LINE_PT
    selected$evidence <- -log10(selected$p.adjust)
    y_spans[[key]] <- max(selected$y_top)+1
    x_max <- max(selected$evidence)*1.04
    # Slim ontology strips use a shared physical width with padding around 7 pt labels.
    strip <- do.call(rbind,lapply(split(selected,selected$ONTOLOGY,drop=TRUE),function(z)
      data.frame(ONTOLOGY=z$ONTOLOGY[1],ymin=min(z$y_bottom)+1,ymax=max(z$y_top)-1)))
    absent <- setdiff(GO_CONFIG$ontologies,as.character(unique(selected$ONTOLOGY)))
    subtitle <- paste0(length(genes),' DEGs; ',nrow(selected),' of ',
      sum(result$p.adjust<GO_CONFIG$go_adjusted_p_cutoff),' significant GO terms shown')
    if(length(absent)) subtitle <- paste0(subtitle,'\nNo significant ',paste(absent,collapse='/'),' terms')
    p <- ggplot(selected)+
      geom_go_roundrect(aes(xmin=0,xmax=evidence,ymin=y-FIG_GO_LINE_PT*.46,ymax=y+FIG_GO_LINE_PT*.46,fill=ONTOLOGY),
        alpha=FIG_GO_BAR_ALPHA,show.legend=FALSE)+
      geom_go_roundrect(data=strip,aes(ymin=ymin,ymax=ymax,fill=ONTOLOGY),xmin=-x_max*.18,xmax=-x_max*.095,
        show.legend=FALSE,alpha=FIG_GO_STRIP_ALPHA,width_pt=FIG_GO_STRIP_WIDTH_PT)+
      geom_text(data=strip,aes(y=(ymin+ymax)/2,label=ONTOLOGY),x=-x_max*.1375,
        family=FIG_FONT,size=FIG_ANNOT_SIZE,color='black')+
      geom_point(aes(y=y,fill=ONTOLOGY),x=-x_max*.052,shape=21,size=FIG_GO_CIRCLE_MM,
        stroke=FIG_LINE_WIDTH,color='black',show.legend=FALSE)+
      geom_text(aes(y=y,label=Count),x=-x_max*.052,family=FIG_FONT,size=FIG_GO_COUNT_SIZE,color='black')+
      geom_text(data=label_lines,aes(y=y,label=line_label),x=x_max*.012,hjust=0,vjust=.5,
        lineheight=1,family=FIG_FONT,size=FIG_ANNOT_SIZE,color='black')+
      geom_text(data=gene_lines,aes(y=gene_y,label=gene_label),
        x=x_max*.012,hjust=0,vjust=.5,family=FIG_FONT,fontface='italic',
        size=FIG_GO_GENE_SIZE,color='#4D4D4D')+
      # An explicit baseline starts at zero, excluding the count/ontology gutter.
      annotate('segment',x=0,xend=x_max,y=0,yend=0,linewidth=FIG_LINE_WIDTH)+
      scale_fill_manual(values=palette)+
      scale_x_continuous(limits=c(-x_max*.19,x_max),expand=c(0,0),
        breaks=local({b<-pretty(c(0,x_max),n=3);b[b>=0 & b<=x_max]}))+
      scale_y_continuous(limits=c(0,max(selected$y_top)+1),expand=c(0,0),breaks=NULL)+
      labs(title=if(direction=='up')'Upregulated genes' else 'Downregulated genes',subtitle=subtitle,
        x=expression(-log[10](italic(p)[plain(BH)])),y=NULL)+
      theme_fig()+
      theme(axis.line=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),
        axis.title.x=element_text(size=7,margin=margin(t=1)),
        plot.title=element_text(face='bold',size=FIG_SUBTITLE_SIZE + 3,margin=margin(b=6)),
        plot.subtitle=element_text(size=FIG_SUBTITLE_SIZE,lineheight=1.05,margin=margin(b=1)),
        plot.margin=margin(if(sex=="male" && direction=="up") 0 else 5,6,3,3))
    plots[[key]]<-p
  }
  fig <- if(sex=='female') 'Fig2A_GO_Female' else 'FigS1E_GO_Male'
  if(sex=='female') {
    layout <- plots$female_up | plots$female_down
  } else {
    # Stack directly, with plotting height proportional to each panel's rows.
    # This keeps male row/bar spacing consistent without a blank spacer panel.
    layout <- (plots$male_up / plots$male_down) +
      plot_layout(heights=c(y_spans$male_up,y_spans$male_down))
  }
  # Direction/term-count headings are needed to identify the separate analyses.
  # Detailed methods live in the companion legend, as for the other panels.
  combined <- layout + plot_annotation(
    title=if(sex=="male") "Male" else "Female",
    theme=theme(
      # Keep the sex heading exactly 2 pt larger than the
      # 10 pt Upregulated/Downregulated heading, with a clear gap below it.
      plot.title=element_text(family=FIG_FONT,size=FIG_SUBTITLE_SIZE + 5,
        face='bold',hjust=.5,margin=margin(b=8)),
      plot.margin=margin(3,3,3,3)))
  fig_save_panel(combined,file.path(out,fig),paste0('GO_',sex))

}
write.csv(do.call(rbind,audit),file.path(out,'analysis_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,selections),file.path(out,'all_displayed_terms_and_genes.csv'),row.names=FALSE)
# Keep the complete gene lists in a clearly named supplementary table.
# Existing raw/result CSVs above remain byte-for-byte compatible.
full_memberships <- do.call(rbind,selections)
keys <- paste(full_memberships$sex,full_memberships$direction,full_memberships$ID)
hkeys <- paste(GO_HIGHLIGHTS$sex,GO_HIGHLIGHTS$direction,GO_HIGHLIGHTS$ID)
idx <- match(keys,hkeys)
full_memberships$highlighted_genes <- ifelse(is.na(idx),'',GO_HIGHLIGHTS$genes[idx])
full_memberships$highlight_reason <- ifelse(is.na(idx),'',GO_HIGHLIGHTS$reason[idx])
write.csv(full_memberships,file.path(out,'GO_displayed_terms_full_gene_lists.csv'),row.names=FALSE)
write.csv(do.call(rbind,highlight_audit),file.path(out,'GO_highlighted_gene_audit.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
writeLines(capture.output(AnnotationDbi::metadata(org.Dm.eg.db)),file.path(out,'annotation_metadata.txt'))
jsonlite::write_json(list(parameters=GO_CONFIG,ontology_colors=as.list(FIG_GO_COLORS),
  package_versions=list(clusterProfiler=as.character(packageVersion('clusterProfiler')),
    org.Dm.eg.db=as.character(packageVersion('org.Dm.eg.db')))),
  file.path(out,'analysis_configuration.json'),pretty=TRUE,auto_unbox=TRUE)
print(do.call(rbind,audit))
