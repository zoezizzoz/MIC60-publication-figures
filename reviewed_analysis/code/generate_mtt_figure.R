# Compatibility entrypoint: use the current pooled-well Figure 4B generator.
# The implementation is shared with the publication panel, avoiding stale plots.
args <- grep('^--file=',commandArgs(FALSE),value=TRUE)
base <- dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',args[1]),fixed=TRUE)))
while(!file.exists(file.path(base,'reproduce.py'))) {
 parent <- dirname(base)
 if(parent==base) stop('Repository root not found')
 base <- parent
}
script <- file.path(base,'panels/Fig4AB/Sources/Figure_4B_MTT/Code/generate_mtt_vector.R')
status <- system2(file.path(R.home('bin'),'Rscript'),shQuote(script))
if(status!=0) stop('Figure 4B rebuild failed')
source_dir <- file.path(dirname(dirname(script)),'Final_Graphs')
targets <- file.path(base,c('reviewed_analysis/results/panels/Fig4B_MTT_descriptive'))
for(stem in targets) {
 dir.create(dirname(stem),recursive=TRUE,showWarnings=FALSE)
 for(ext in c('pdf','png')) if(!file.copy(file.path(source_dir,paste0('Fig4B_MTT_Vector.',ext)),paste0(stem,'.',ext),overwrite=TRUE)) stop('Export copy failed')
}
