# Reconstructions where the supplied original generator is unreadable.
# Data and shared plotting helpers are unchanged; original source paths are in the audit.
suppressPackageStartupMessages({library(ggplot2);library(readxl);library(patchwork)})
script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
root <- dirname(dirname(script));source(file.path(dirname(script),'figure_style.R'))
out <- file.path(root,'results/full_figures');dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(FIG_SEED)
# Locomotion: exact rank permutation preserves the observed tie and all 8 values.
d <- read.csv(file.path(root,'data/locomotion.csv'));d$group <- factor(d$group,levels=c('WT','CS'))
ranks <- rank(d$performance_index);observed <- sum(ranks[d$group=='CS'])
null <- combn(ranks,4,FUN=sum);P <- mean(abs(null-mean(null))>=abs(observed-mean(null))-1e-12)
stopifnot(abs(P-2/70)<1e-12)
b <- do.call(rbind,lapply(levels(d$group),function(g) {
 x<-d$performance_index[d$group==g];z<-fivenum(x)
 data.frame(group=g,ymin=z[1],lower=z[2],middle=z[3],upper=z[4],ymax=z[5])
}))
p <- ggplot(d,aes(group,performance_index,fill=group))+
 geom_boxplot(data=b,aes(x=group,ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,fill=group),inherit.aes=FALSE,stat='identity',width=FIG_SUMMARY_WIDTH,linewidth=FIG_LINE_WIDTH)+
 geom_errorbar(data=b,aes(x=group,ymin=ymin,ymax=ymax),inherit.aes=FALSE,width=FIG_EB_CAP,linewidth=FIG_LINE_WIDTH)+
 fig_points()+fig_sig_bracket(P,d$performance_index)+fig_scale_fill(c('WT','CS'))+
 fig_scale_x_group(FIG_GROUP_LABELS[c('WT','CS')])+fig_scale_y('Performance index')+labs(x=NULL)+theme_fig()
fig_save(p,file.path(out,'Fig3H_Locomotion'),width=130/72,height=125/72)
write.csv(data.frame(n_WT=4,n_CS=4,allocations=length(null),exact_rank_permutation_P=P),file.path(out,'Fig3H_statistics.csv'),row.names=FALSE)
# TIMELESS/DAPI: preserve all supplied values; no culture-level inference without IDs.
w <- as.data.frame(read_excel(file.path(root,'data/TIMELESS_DAPI.xlsx')))
stopifnot(ncol(w)==4)
w <- as.data.frame(lapply(w,function(x) x[is.finite(x)]),check.names=FALSE)
stopifnot(nrow(w)==16)
print(names(w))
# The workbook order is explicitly checked to avoid assigning conditions by guesswork.
stopifnot(identical(names(w),c('WT','CS','WT-H2O2','CS-H2O2')))
d <- data.frame(value=unlist(w,use.names=FALSE),genotype=factor(rep(c('WT','CS','WT','CS'),each=nrow(w)),levels=c('WT','CS')),
 treatment=factor(rep(c('Vehicle','Vehicle','H2O2','H2O2'),each=nrow(w)),levels=c('Vehicle','H2O2')))
p <- ggplot(d,aes(genotype,value,fill=paste(genotype,treatment)))+fig_summary()+fig_points()+
 scale_fill_manual(values=FIG_TREATMENT_COLORS)+
 fig_scale_x_group(FIG_GROUP_LABELS[c('WT','CS')])+fig_scale_y('Nuclear TIMELESS / DAPI')+facet_wrap(~treatment,nrow=1,labeller=as_labeller(c(Vehicle="plain(Vehicle)",H2O2='"20 mM"~H[2]*O[2]'),default=label_parsed))+labs(x=NULL)+theme_fig()+
 theme(strip.text=element_blank(),strip.background=element_blank(),
       axis.text.x=element_text(size=8),
       plot.margin=margin(3,4,3,9))
# Put treatment indicators below the genotype axis, aligned with box centres.
# Concentration remains 20 mM for the positive treatment, as in the source data.
timeless_treatment_row <- function(plot, signs, text_size, treatment_size) {
  # Match the export device when measuring Arial label widths.
  metrics_file <- tempfile(fileext=".pdf")
  grDevices::quartz(type="pdf", file=metrics_file)
  on.exit({grDevices::dev.off(); unlink(metrics_file)}, add=TRUE)
  built <- ggplot_build(plot)
  g <- ggplotGrob(plot)
  panels <- g$layout[grepl("^panel($|-)", g$layout$name), ]
  panels <- panels[order(panels$l), ]
  row_after <- max(g$layout$b[grepl("^axis-b", g$layout$name)])
  g <- gtable::gtable_add_rows(g, grid::unit(max(text_size * 1.6, treatment_size * 1.4), "pt"), pos=row_after)
  gp <- grid::gpar(fontfamily=FIG_FONT, fontsize=text_size, col="black")
  for (i in seq_len(nrow(panels))) {
    xpos <- built$layout$panel_params[[i]]$x$break_positions()
    stopifnot(length(xpos) == length(signs[[i]]))
    marks <- grid::textGrob(signs[[i]], x=grid::unit(xpos,"npc"), y=.5, gp=gp)
    g <- gtable::gtable_add_grob(g, marks, t=row_after+1, l=panels$l[i], r=panels$r[i], clip="off")
  }
  treatment <- grid::textGrob(expression(H[2]*O[2]),
    x=grid::unit(1,"npc")-grid::unit(4,"pt"), y=.5, just="right",
    gp=grid::gpar(fontfamily=FIG_FONT, fontsize=treatment_size, col="black"))
  gtable::gtable_add_grob(g, treatment, t=row_after+1, l=1, r=min(panels$l)-1, clip="off")
}

p <- timeless_treatment_row(p, list(c("\u2212", "\u2212"), c("+", "+")), text_size=8, treatment_size=12)
fig_save(p,file.path(out,'Fig4A_TIMELESS_Descriptive'),width=260/72,height=160/72)
write.csv(d,file.path(out,'Fig4A_plotted_values.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'additional_assays_sessionInfo.txt'))
