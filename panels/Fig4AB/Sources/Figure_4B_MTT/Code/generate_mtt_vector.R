suppressPackageStartupMessages({library(ggplot2);library(jsonlite);library(grid);library(gtable)})
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
panel <- dirname(dirname(script));root <- dirname(dirname(panel))
s <- jsonlite::read_json(file.path(root,'Code/render_settings.json'),simplifyVector=TRUE)
font <- s$font_family
grDevices::pdfFonts(Arial=grDevices::pdfFonts("Helvetica")[[1]])
base_theme <- function() theme_classic(base_family=font,base_size=s$text_pt)+theme(
 text=element_text(family=font,size=s$text_pt,color='black'),
 axis.text=element_text(size=s$text_pt,color='black'),axis.title=element_text(size=s$text_pt),
 axis.line=element_line(linewidth=.20,color='black'),axis.ticks=element_line(linewidth=.20,color='black'),
 axis.ticks.length=unit(2,'pt'),axis.title.x=element_text(margin=margin(t=3)),
 axis.title.y=element_text(margin=margin(r=3)),plot.title=element_text(size=s$text_pt,hjust=.5,margin=margin(b=3)),
 legend.text=element_text(size=s$text_pt),plot.margin=margin(3,4,3,3))
save_plot <- function(plot,stem,width,height) {
 grDevices::quartz(type='pdf',file=file.path(panel,'Final_Graphs',paste0(stem,'.pdf')),width=width/72,height=height/72,family=font)
 if(inherits(plot,'ggplot'))print(plot) else {grid.newpage();grid.draw(plot)}
 grDevices::dev.off()
 writeLines(capture.output(sessionInfo()),file.path(panel,'Documentation/R_sessionInfo.txt'))
}

# Uses the archived generate_mtt_figure.R analysis and final September inside-legend placement.
d <- read.csv(file.path(panel,'Supporting_Data/mtt_plot_values.csv'))
d$genotype <- factor(d$genotype,levels=c('WT','CS'))
stopifnot(nrow(d)==78,all(is.finite(d$normalized_pct)))
# Independent check against the raw-sheet mapping in R, without using Python's normalized column.
suppressPackageStartupMessages(library(readxl))
source <- file.path(panel,'Original_Data/MTT Data.xlsx')
cell_value <- function(sheet,address) as.numeric(read_excel(source,sheet=sheet,range=address,col_names=FALSE,.name_repair='minimal')[[1]][1])
raw_check <- vapply(seq_len(nrow(d)),function(i) cell_value(d$source_sheet[i],d$source_cell[i]),numeric(1))
stopifnot(max(abs(raw_check-d$raw_absorbance))<1e-12)
corrected <- raw_check-d$blank
control <- ave(ifelse(d$dose_mM==0,corrected,NA_real_),d$experiment,d$genotype,FUN=function(x)mean(x,na.rm=TRUE))
stopifnot(max(abs(100*corrected/control-d$normalized_pct))<1e-10)
experiment_summ <- do.call(rbind,lapply(split(d,list(d$experiment,d$dose_mM,d$genotype),drop=TRUE),function(x)
 data.frame(experiment=x$experiment[1],dose_mM=x$dose_mM[1],genotype=x$genotype[1],experiment_mean=mean(x$normalized_pct),technical_sd=sd(x$normalized_pct),technical_n=nrow(x))))
summ <- do.call(rbind,lapply(split(experiment_summ,list(experiment_summ$dose_mM,experiment_summ$genotype),drop=TRUE),function(x)
 data.frame(dose_mM=x$dose_mM[1],genotype=x$genotype[1],mean=mean(x$experiment_mean),sd=if(nrow(x)>1)sd(x$experiment_mean) else NA_real_,biological_n=nrow(x),technical_n=sum(x$technical_n))))
experiment_summ$experiment <- factor(experiment_summ$experiment,levels=c(1,2));experiment_summ$genotype <- factor(experiment_summ$genotype,levels=c('WT','CS'));summ$genotype <- factor(summ$genotype,levels=c('WT','CS'))
stopifnot(nrow(experiment_summ)==16,nrow(summ)==10)
p <- ggplot(d,aes(dose_mM,normalized_pct))+
 geom_hline(yintercept=100,color='#A6A6A6',linetype='dashed',linewidth=.20)+
 geom_point(color='#242424',size=.9,alpha=.72,show.legend=FALSE)+
 geom_line(data=summ,aes(y=mean,color=genotype,group=genotype),linewidth=.30)+
 geom_errorbar(data=summ,aes(y=mean,ymin=mean-sd,ymax=mean+sd,group=genotype),width=1.8,linewidth=.20,na.rm=TRUE)+
 geom_point(data=experiment_summ,aes(y=experiment_mean,color=genotype,shape=experiment),size=1.8)+
 scale_x_continuous(breaks=c(0,5,10,20,40),expand=expansion(add=2.2))+
 scale_y_continuous(limits=c(0,150),breaks=seq(0,125,25),expand=expansion(mult=0))+
 scale_color_manual(name=NULL,values=c(WT=s$WT,CS=s$CS),breaks=c('WT','CS'),labels=c('dMIC60-WT','dMIC60-CS'))+
 scale_shape_manual(name=NULL,values=c('1'=16,'2'=17),labels=c('Experiment 1','Experiment 2'))+
 labs(x=expression(H[2]*O[2]~'concentration (mM)'),y='Metabolic activity (MTT)\n(% of matched 0 mM control)',title=expression(italic(MIC60)*'-null HeLa cells'))+
 guides(color=guide_legend(order=1,nrow=1,override.aes=list(shape=16)),shape=guide_legend(order=2,nrow=1,override.aes=list(color='black')))+
 base_theme()+theme(legend.position='inside',legend.position.inside=c(.5,.97),legend.justification.inside=c(.5,1),legend.direction='horizontal',legend.box='vertical',legend.key.size=unit(8,'pt'),legend.key.width=unit(9,'pt'),legend.key.height=unit(8,'pt'),legend.spacing.x=unit(3,'pt'),legend.spacing.y=unit(1,'pt'),legend.box.spacing=unit(0,'pt'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),legend.background=element_blank())
built <- ggplot_build(p);stopifnot(nrow(built$data[[2]])==78,nrow(built$data[[5]])==16)
write.csv(experiment_summ,file.path(panel,'Supporting_Data/mtt_experiment_means.csv'),row.names=FALSE)
write.csv(summ,file.path(panel,'Supporting_Data/mtt_biological_summary.csv'),row.names=FALSE)
save_plot(p,'Fig4B_MTT_Vector',s$B_width_pt,s$B_height_pt)
