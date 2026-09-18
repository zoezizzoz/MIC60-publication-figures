# Descriptive replacement MTT panel and validation of saved blot/efficiency data.
suppressPackageStartupMessages(library(ggplot2))
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), fixed=TRUE))
root<-dirname(dirname(script));out<-file.path(root,'results');dir.create(file.path(out,'panels'),showWarnings=FALSE)
source(file.path(dirname(script),'generate_mtt_figure.R'))
wb<-read.csv(file.path(root,'data','western_blot.csv'))
wt<-wb[wb$genotype=='dMIC60WT',];cs<-wb[wb$genotype=='dMIC60CS',]
a<-merge(wt,cs,by='blot_date',suffixes=c('_WT','_CS'))
test<-t.test(a$normalized_MIC60_Myc_over_ATP5B_WT,a$normalized_MIC60_Myc_over_ATP5B_CS,paired=TRUE)
write.csv(data.frame(n_pairs=nrow(a),t=unname(test$statistic),df=unname(test$parameter),P=test$p.value),file.path(out,'western_blot_audit.csv'),row.names=FALSE)
te<-read.csv(file.path(root,'data','transfection_efficiency.csv'))
write.csv(aggregate(efficiency~experiment+group,te,mean),file.path(out,'transfection_experiment_means.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'panels','sessionInfo.txt'))
