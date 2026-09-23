suppressPackageStartupMessages(library(data.table))
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
r <- file.path(dirname(dirname(script)), 'results/sleep')
w <- r
# A stratified rank test is sensitivity analysis, with recording start date as stratum.
# For each stratum: centered rank sum and exact null variance accounting for tied ranks.
# Weight 1/(N+1) prevents larger strata from dominating through the rank scale alone.
rows <- list();cohorts <- list()
for(metric in c('sleep','activity')) for(ph in c('Total','Day','Night')) {
 z <- fread(file=file.path(r,paste0(metric,if(ph=='Total') '_total' else '_phase','_values_and_exclusions.csv')))[exclude==FALSE]
 if(ph!='Total') z <- z[phase==ph]
 value <- if(metric=='sleep') 'percent_asleep' else 'activity'
 a <- z[,{
  ranks<-rank(get(value));N<-.N;n1<-sum(genotype=='CS');n0<-N-n1
  centered<-sum(ranks[genotype=='CS'])-n1*(N+1)/2
  variance<-n1*n0/(N*(N-1))*sum((ranks-mean(ranks))^2)
  list(WT_n=n0,CS_n=n1,WT_mean=mean(get(value)[genotype=='WR']),CS_mean=mean(get(value)[genotype=='CS']),centered=centered,variance=variance,weight=1/(N+1))
 },by=cohort_start]
 Z<-sum(a$weight*a$centered)/sqrt(sum(a$weight^2*a$variance))
 rows[[length(rows)+1]]<-data.table(metric,phase=ph,Z,P=2*pnorm(-abs(Z)))
 a[,`:=`(metric=metric,phase=ph)];cohorts[[length(cohorts)+1]]<-a
}
x<-rbindlist(rows);x[,Holm_P:=p.adjust(P,'holm')]
fwrite(x,file.path(w,'sleep_cohort_stratified_sensitivity.csv'));fwrite(rbindlist(cohorts),file.path(w,'sleep_per_cohort_summary.csv'));print(x)
for(metric in c('sleep','activity')){
 z<-fread(file=file.path(r,paste0(metric,'_profile_summary.csv')));print(z[,.(minimum_n=min(flies),maximum_n=max(flies)),by=genotype])
}
