suppressPackageStartupMessages(library(ggplot2))
args <- gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE)
script <- normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1]))
root <- dirname(dirname(script))
source_file <- file.path(root,'figure_sources/Figure_3K_Mitochondrial_Area/Rebuilt_Output/TEM_2026-08-04_image_level_summary.csv')
out <- file.path(root,'results/full_figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(source_file,stringsAsFactors=FALSE)
specs <- list(
  J=list(metric='Mitochondrial perimeter',label=expression('Mean perimeter ('*mu*'m)'),limits=c(0,10),breaks=c(0,2.5,5,7.5,10)),
  K=list(metric='Mitochondrial area',label=expression('Mean area ('*mu*'m'^2*')'),limits=c(0,5.2),breaks=0:5),
  L=list(metric='Aspect ratio',label='Mean aspect ratio',limits=c(1.4,3.6),breaks=c(1.5,2,2.5,3,3.5))
)
for(panel in names(specs)){
  s <- specs[[panel]]
  x <- d[d$Panel==panel & d$Metric==s$metric,,drop=FALSE]
  stopifnot(nrow(x)>=20,!anyNA(x$Image_Mean),!anyDuplicated(paste(x$Image,x$Condition)))
  x$Condition <- factor(x$Condition,levels=c('WR','CS'))
  set.seed(3)
  p <- ggplot(x,aes(Condition,Image_Mean,fill=Condition))+
    geom_boxplot(width=.42,outlier.shape=NA,linewidth=.33)+
    geom_point(position=position_jitter(width=.08,height=0,seed=3),size=1.1,color='#222222')+
    scale_fill_manual(values=c(WR='#5AB4E5',CS='#CB78A8'),guide='none')+
    scale_x_discrete(labels=c('dMIC60-\nWT','dMIC60-\nCS'))+
    scale_y_continuous(limits=s$limits,breaks=s$breaks,expand=expansion(mult=c(0,.02)))+
    labs(x=NULL,y=s$label)+theme_classic(base_family='sans',base_size=7)+
    theme(axis.title.y=element_text(size=7),axis.text=element_text(size=6.5,color='black'),
          axis.text.x=element_text(lineheight=.8),plot.margin=margin(5,5,5,5))
  ggsave(file.path(out,paste0('Fig3',panel,'_TEM_field_means.pdf')),p,width=1.75,height=2.1,device=grDevices::pdf)
  cat(panel,table(x$Condition),'\n')
}
