# Edited from the archived generate_mtt_figure.R.
# Technical wells are descriptive; no significance tests/brackets are drawn.
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep('^--file=',commandArgs(FALSE),value=TRUE)
script_dir <- dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',script_arg[[1]]),fixed=TRUE)))
figure_dir <- dirname(script_dir)
source(file.path(script_dir,'figure_style.R'))
d <- read.csv(file.path(figure_dir, "results", "MTT_reconstructed_wells.csv"), stringsAsFactors = FALSE)
d$genotype <- factor(d$genotype, levels = c("WT", "CS"))
d$x_plot <- d$dose_mM
doses <- c(0, 5, 10, 20, 40)
experiment_summ <- do.call(rbind, lapply(split(d, list(d$experiment, d$dose_mM, d$genotype), drop = TRUE), function(x) {
    data.frame(experiment = x$experiment[1], dose_mM = x$dose_mM[1], genotype = x$genotype[1], experiment_mean = mean(x$normalized_pct), 
        technical_sd = sd(x$normalized_pct), technical_n = nrow(x))
}))
experiment_summ$genotype <- factor(experiment_summ$genotype, levels = c("WT", "CS"))
experiment_summ$x_plot <- experiment_summ$dose_mM
summ <- do.call(rbind, lapply(split(experiment_summ, list(experiment_summ$dose_mM, experiment_summ$genotype)), function(x) {
    data.frame(dose_mM = x$dose_mM[1], genotype = x$genotype[1], mean = mean(x$experiment_mean), sd = if (nrow(x) > 1) 
        sd(x$experiment_mean)
    else NA_real_, biological_n = nrow(x), technical_n = sum(x$technical_n))
}))
summ$genotype <- factor(summ$genotype, levels = c("WT", "CS"))
summ$x_plot <- summ$dose_mM
p <- ggplot(d, aes(x_plot, normalized_pct)) + fig_reference_line(100) + geom_point(color = FIG_PT_COLOR, position = "identity", 
    shape = FIG_PT_SHAPE, size = FIG_PT_SIZE_DENSE, alpha = 0.72, show.legend = FALSE) + geom_line(data = summ, aes(x = x_plot, y = mean, 
    group = genotype, color = genotype), linewidth = FIG_PROFILE_LINEWIDTH) + fig_errorbar(data = summ, mapping = aes(x = x_plot, 
    y = mean, ymin = mean - sd, ymax = mean + sd, group = genotype), width = 1.8, na.rm = TRUE) + geom_point(data = experiment_summ, aes(x = x_plot, 
    y = experiment_mean, shape = factor(experiment), color = genotype), size = FIG_MARKER_SIZE, alpha = 0.95) + scale_x_continuous(breaks = doses, labels = doses,
    expand = expansion(add = 2.2)) + fig_scale_color(c('WT','CS'),labels=FIG_GROUP_LABELS[c('WT','CS')]) +
    fig_scale_y("Metabolic activity (MTT; % of matched 0 mM control)", limits=c(0,130),breaks=seq(0,125,25)) +
    labs(x=expression(H[2]*O[2]~"concentration (mM)"),color=NULL,shape='Experiment') +
    scale_shape_manual(values=c(16,17),labels=c('1','2')) +
    guides(color=guide_legend(order=1,override.aes=list(linewidth=FIG_PROFILE_LINEWIDTH)),shape=guide_legend(order=2)) +
    theme_fig(legend_position='top') +
    theme(legend.title=element_text(size=FIG_LEGEND_TEXT_SIZE),legend.box='vertical',legend.spacing.y=grid::unit(1,'pt'),
      legend.box.spacing=grid::unit(2,'pt'),legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0),axis.title.y=element_text(margin=margin(r=3)),
      axis.title.x=element_text(margin=margin(t=3)))
fig_save_panel(p,file.path(figure_dir,'results/panels/Fig4B_MTT_descriptive'),'MTT')
