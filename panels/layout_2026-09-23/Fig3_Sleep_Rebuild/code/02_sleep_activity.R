# Reproducible audit of analysis windows and genotype/phase filtering.
# Primary corrected window: 12 <= elapsed hours < 60, matching the methods.
# Light schedule confirmed by author: 10:00 lights on, 22:00 lights off.
# Preserve the documented 1.5-IQR rule; report unfiltered and cohort sensitivity results.
suppressPackageStartupMessages({library(damr); library(sleepr); library(data.table); library(ggplot2)})
script <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1]), fixed=TRUE))
root <- dirname(dirname(script)); input <- file.path(root,'data','sleep'); output <- file.path(root,'results','sleep')
dir.create(output, recursive=TRUE, showWarnings=FALSE)
metadata <- fread(file=file.path(input,'sleep_metadata.csv'))
linked <- link_dam_metadata(metadata, result_dir=input)
raw <- load_dam(linked, FUN=sleepr::sleep_dam_annotation)
raw_meta <- as.data.table(meta(raw))
fwrite(raw_meta[,names(raw_meta)[!vapply(raw_meta,is.list,logical(1))],with=FALSE],file.path(output,'raw_fly_metadata.csv'))
curated <- curate_dead_animals(raw)
female <- curated[xmv(sex)=='Female' & xmv(genotype) %in% c('WR','CS')]
d <- as.data.table(female); m <- as.data.table(meta(female))
d <- merge(d, m[, .(id, genotype, replicate)], by='id')
d[,cohort_start:=sub('\\|.*$','',as.character(id))]
cadence <- d[,.(interval_seconds=median(diff(sort(unique(t))))),by=id]
stopifnot(all(cadence$interval_seconds==300))
fwrite(cadence,file.path(output,'acquisition_intervals.csv'))

d[, hours:=t/3600]
d[, phase:=ifelse(((hours-12)%%24)>=12,'Night','Day')]
stopifnot(all(is.finite(d$activity)))
fwrite(as.data.table(metadata), file.path(output,'metadata_used.csv'))
fwrite(d[,.(n_timepoints=.N,first_hour=min(hours),last_hour=max(hours)),by=.(id,genotype,replicate)],file.path(output,'curated_fly_coverage.csv'))
# Correct the activity-profile scale as well as the summary-panel scales.
for(metric in c('sleep','activity')) {
  z<-copy(d[hours>=12 & hours<60]);z[,bin:=floor((hours-12)*2)/2]
  per_fly<-z[,.(value=if(metric=='sleep')mean(asleep)*30 else mean(activity)/5),by=.(id,genotype,bin)]
  profile<-per_fly[,.(mean=mean(value),SEM=sd(value)/sqrt(.N),flies=.N),by=.(genotype,bin)]
  fwrite(profile,file.path(output,paste0(metric,'_profile_summary.csv')))

}
filter_iqr <- function(z, value, paired=FALSE) {
  z <- copy(z)
  bycols <- if(paired)c('genotype','phase') else 'genotype'
  z[, c('lower','upper'):= {q<-quantile(get(value),c(.25,.75),na.rm=TRUE);list(q[1]-1.5*diff(q),q[2]+1.5*diff(q))},by=bycols]
  z[,exclude:=!is.finite(get(value))|get(value)<lower|get(value)>upper]
  if(paired)z[,exclude:=any(exclude)|uniqueN(phase)!=2,by=id]
  z
}
alltests <- list(); k<-0
for(window in c('legacy_0_60','corrected_12_60')) {
  z<-d[hours<60 & hours>=ifelse(window=='corrected_12_60',12,0)]
  for(metric in c('sleep','activity')) {
    value<-if(metric=='sleep')'percent_asleep' else 'activity'
    totals<-z[,.(percent_asleep=100*mean(asleep,na.rm=TRUE),activity=mean(activity,na.rm=TRUE)),by=.(id,genotype,replicate,cohort_start)]
    phases<-z[,.(percent_asleep=100*mean(asleep,na.rm=TRUE),activity=mean(activity,na.rm=TRUE)),by=.(id,genotype,replicate,cohort_start,phase)]
    for(mode in c('unfiltered','legacy_independent_phase_IQR','paired_phase_IQR')) {
      tt<-if(mode=='unfiltered')copy(totals)[,exclude:=FALSE] else filter_iqr(totals,value)
      pp<-if(mode=='unfiltered')copy(phases)[,exclude:=FALSE] else if(mode=='paired_phase_IQR')filter_iqr(phases,value,TRUE) else {
        a<-copy(phases);a[,c('lower','upper'):={q<-quantile(get(value),c(.25,.75),na.rm=TRUE);list(q[1]-1.5*diff(q),q[2]+1.5*diff(q))},by=.(genotype,phase)];a[,exclude:=!is.finite(get(value))|get(value)<lower|get(value)>upper];a
      }
      for(phase in c('Total','Day','Night')) {
        phase_name<-phase
        a<-if(phase=='Total')tt[exclude==FALSE] else pp[exclude==FALSE & get("phase")==phase_name]
        x<-a[genotype=='WR'][[value]];y<-a[genotype=='CS'][[value]]
        k<-k+1;alltests[[k]]<-data.table(window,metric,mode,phase,WT_n=length(x),CS_n=length(y),WT_mean=mean(x),CS_mean=mean(y),WT_median=median(x),CS_median=median(y),P=wilcox.test(x,y,exact=FALSE,correct=TRUE)$p.value)
      }
      if(window=='corrected_12_60' && mode=='paired_phase_IQR') {
        # Rescale only after rank tests to preserve identical floating-point tie handling.
        if(metric=='activity'){tt[,activity:=activity/5];pp[,activity:=activity/5]}
        fwrite(tt,file.path(output,paste0(metric,'_total_values_and_exclusions.csv')))
        fwrite(pp,file.path(output,paste0(metric,'_phase_values_and_exclusions.csv')))
        # Descriptive batch summaries retain experiment structure for review.
        fwrite(pp[exclude==FALSE,lapply(.SD,mean),by=.(cohort_start,genotype,phase),.SDcols=value],file.path(output,paste0(metric,'_cohort_means.csv')))

      }
    }
  }
}
tests<-rbindlist(alltests)
tests[metric=="activity",c("WT_mean","CS_mean","WT_median","CS_median"):=lapply(.SD,function(x)x/5),.SDcols=c("WT_mean","CS_mean","WT_median","CS_median")]
tests[,Holm_P_six_summary_tests:=p.adjust(P,'holm'),by=.(window,mode)]
fwrite(tests,file.path(output,'window_and_filter_sensitivity.csv'))
print(tests[window=='corrected_12_60' & mode=='paired_phase_IQR'])
# Primary analysis and per-bin n are exported explicitly for captions/reproduction.
primary <- tests[window=='corrected_12_60' & mode=='paired_phase_IQR']
fwrite(primary,file.path(output,'primary_statistics.csv'))
fwrite(d[hours>=12 & hours<60,.(flies=uniqueN(id)),by=.(cohort_start,genotype)],file.path(output,'cohort_sample_counts.csv'))
jsonlite::write_json(list(window_hours=c(12,60),upper_bound_exclusive=TRUE,
  lights_on='10:00',lights_off='22:00',light_schedule_author_confirmed=TRUE,
  acquisition_interval_seconds=300,sleep_minimum_immobility_seconds=300,
  death_curation='sleepr::curate_dead_animals; 24-hour windows, <=1% moving, hourly evaluation',
  IQR_multiplier=1.5,phase_exclusions_paired=TRUE,
  test='two-sided Wilcoxon rank-sum, normal approximation with tie and continuity correction',
  adjustment='Holm across six total/day/night sleep/activity genotype comparisons',
  cohort_analysis='Sensitivity only; four recording start dates; replicate column is not an independent-run ID',
  activity_units='beam crossings per minute',profiles='mean and SEM across available flies per 30-minute bin'),
  file.path(output,'analysis_parameters.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(output,'sessionInfo.txt'))

# Use the edited original plot builders and common style.
source(file.path(dirname(script),'generate_sleep_activity_figures.R'))

# Record the cohort-stratified sensitivity analysis alongside the primary tests.
source(file.path(dirname(script),'sleep_cohort_sensitivity.R'))
