# Fit each sex once, validate archived results, and reuse normalized counts.
suppressPackageStartupMessages(library(DESeq2))
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), fixed=TRUE))
root<-dirname(dirname(script));out<-file.path(root,'results','rnaseq');dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-as.matrix(read.csv(file.path(root,'data','rnaseq_counts.csv'),row.names=1,check.names=FALSE))
stopifnot(all(is.finite(x)),all(x>=0),all(x==round(x)),!anyDuplicated(rownames(x)))
storage.mode(x)<-'integer'
checks<-list()
for(sex in c('female','male')) {
  samples<-colnames(x)[grepl(if(sex=='female')'F[123]$' else 'M[123]$',colnames(x))]
  genotype<-factor(ifelse(grepl('^CS',samples),'CS','WR'),levels=c('WR','CS'))
  stopifnot(length(samples)==6,all(table(genotype)==3))
  dds<-DESeqDataSetFromMatrix(x[,samples],data.frame(row.names=samples,genotype),design=~genotype)
  dds<-DESeq(dds,quiet=TRUE)
  # Default results alpha=0.1 reproduces the archived pipeline; DEG threshold remains 0.05.
  res<-as.data.frame(results(dds,contrast=c('genotype','CS','WR')));res$gene<-rownames(res)
  archived<-read.csv(file.path(root,'data',paste0('rnaseq_',sex,'.csv')))
  a<-merge(res,archived,by='gene',suffixes=c('_new','_archived'))
  checks[[sex]]<-data.frame(sex,genes=nrow(res),testable=sum(!is.na(res$padj)),up=sum(res$padj<.05 & res$log2FoldChange>=.58,na.rm=TRUE),down=sum(res$padj<.05 & res$log2FoldChange<=-.58,na.rm=TRUE),max_LFC_difference=max(abs(a$log2FoldChange_new-a$log2FoldChange_archived),na.rm=TRUE),max_padj_difference=max(abs(a$padj_new-a$padj_archived),na.rm=TRUE),NA_padj_disagreements=sum(is.na(a$padj_new)!=is.na(a$padj_archived)))
  write.csv(res,file.path(out,paste0(sex,'_DE.csv')),row.names=FALSE)
  write.csv(counts(dds,normalized=TRUE),file.path(out,paste0(sex,'_normalized_counts.csv')))
}
write.csv(do.call(rbind,checks),file.path(out,'archive_comparison.csv'),row.names=FALSE)
print(do.call(rbind,checks))
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
