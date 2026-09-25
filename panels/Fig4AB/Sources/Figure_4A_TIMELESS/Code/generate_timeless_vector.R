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

# Adapted from 06_additional_assays.R; continuous baseline follows the later recorded layout.
# All 64 observations are retained. Missing experiment IDs preclude biological-replicate tests.
d <- read.csv(file.path(panel,'Supporting_Data/TIMELESS_plotted_values.csv'),check.names=FALSE)
groups <- c('WT','CS','WT +H2O2','CS +H2O2');d$group <- factor(d$group,levels=groups)
stopifnot(nrow(d)==64,all(table(d$group)==16),all(is.finite(d$value)))
positions <- c(0,40,100,140);d$x <- positions[as.integer(d$group)]
fills <- setNames(c(s$WT,s$CS,s$WT_H2O2,s$CS_H2O2),groups)
p <- ggplot(d,aes(x=x,y=value,group=group,fill=group))+
 stat_boxplot(geom='errorbar',width=s$A_box_width_data*.5,linewidth=.20)+
 geom_boxplot(width=s$A_box_width_data,outlier.shape=NA,linewidth=.20)+
 geom_point(position=position_jitter(width=4,height=0,seed=s$point_seed),size=.85,alpha=.82,color='#242424')+
 scale_fill_manual(values=fills,guide='none')+
 scale_x_continuous(breaks=positions,labels=c('dMIC60-\nWT','dMIC60-\nCS','dMIC60-\nWT','dMIC60-\nCS'),limits=c(-24,164),expand=expansion(mult=0))+
 scale_y_continuous(limits=c(0,1.3),breaks=c(0,.5,1),expand=expansion(mult=0))+
 labs(x=NULL,y='TIMELESS / DAPI',title=expression(italic(MIC60)*'-null HeLa cells'))+base_theme()+theme(axis.text.x=element_text(lineheight=.9))
# Treatment signs align with the same four category centers; no faceted baseline.
grDevices::quartz(type='pdf',file=tempfile(fileext='.pdf'),width=s$A_width_pt/72,height=s$A_height_pt/72,family=font)
built <- ggplot_build(p);g <- ggplotGrob(p)
pr <- g$layout[g$layout$name=='panel',];after <- max(g$layout$b[grepl('^axis-b',g$layout$name)])
g <- gtable_add_rows(g,unit(10,'pt'),pos=after)
signs <- textGrob(c('-','-','+','+'),x=unit(built$layout$panel_params[[1]]$x$break_positions(),'npc'),y=.5,gp=gpar(fontfamily=font,fontsize=s$text_pt))
g <- gtable_add_grob(g,signs,t=after+1,l=pr$l,r=pr$r,clip='off')
treatment <- textGrob(expression(H[2]*O[2]),x=unit(1,'npc')-unit(3,'pt'),y=.5,just='right',gp=gpar(fontfamily=font,fontsize=s$text_pt))
g <- gtable_add_grob(g,treatment,t=after+1,l=1,r=pr$l-1,clip='off')
grDevices::dev.off()
stopifnot(nrow(built$data[[3]])==64)
write.csv(built$data[[2]][,c('group','x','ymin','lower','middle','upper','ymax')],file.path(panel,'Supporting_Data/boxplot_summary.csv'),row.names=FALSE)
save_plot(g,'Fig4A_TIMELESS_Vector',s$A_width_pt,s$A_height_pt)
