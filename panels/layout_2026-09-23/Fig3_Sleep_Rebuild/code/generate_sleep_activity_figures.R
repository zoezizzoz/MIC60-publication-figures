# Edited from the archived generate_sleep_activity_figures.R.
# Original plot builders/style retained; corrected numerical outputs are loaded once.
# Statistics are calculated before activity-unit rescaling in 02_sleep_activity.R.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
SCRIPT_DIR <- dirname(normalizePath(gsub('~+~',' ',sub('^--file=','',script_arg[[1]]),fixed=TRUE)))
source(file.path(SCRIPT_DIR, 'figure_style.R'))
FIG_DIR <- file.path(dirname(SCRIPT_DIR), 'results/sleep')
FIG_P_ADJUST_LABEL <- 'Holm'
tests <- fread(file=file.path(FIG_DIR, 'window_and_filter_sensitivity.csv'))
tests <- tests[window=='corrected_12_60' & mode=='paired_phase_IQR']
reviewed_p <- function(metric_name, phase_name) {
    p <- tests[metric==metric_name & phase==phase_name, Holm_P_six_summary_tests]
    stopifnot(length(p)==1L, is.finite(p)); p
}
genotype_labels <- FIG_GROUP_LABELS[c("WR","CS")]
geno_order <- c("WR", "CS")
fly_n_caption <- function(data, phase_col = NULL) {
    d <- unique(as.data.table(data)[genotype %in% geno_order, c("id", "genotype", phase_col), with = FALSE])
    if (is.null(phase_col)) {
        counts <- d[, .(n = uniqueN(id)), by = genotype]
        n_for <- function(g) counts[as.character(genotype) == g, n]
        return(stringr::str_wrap(paste0("Adult female flies monitored individually; WT n = ", n_for("WR"), ", CS n = ", n_for("CS"),
            "."), width = 72))
    }
    counts <- d[, .(n = uniqueN(id)), by = c(phase_col, "genotype")]
    phases <- unique(as.character(d[[phase_col]]))
    phase_text <- vapply(phases, function(ph) {
        n_for <- function(g) counts[get(phase_col) == ph & as.character(genotype) == g, n]
        paste0(ph, ": WT n = ", n_for("WR"), ", CS n = ", n_for("CS"))
    }, character(1))
    stringr::str_wrap(paste0("Adult female flies monitored individually; ", paste(phase_text, collapse = "; "), "."), width = 72)
}
phase_fill <- c(on = "#FFF6C2", off = "#D9D9D9")
phase_bar <- c(on = "#FFD700", off = "black")
two_group_boxplot <- function(data, yvar, title, ylab, metric, ylim = NULL) {
    plot_data <- as.data.table(data)[genotype %in% geno_order]
    plot_data[, `:=`(genotype, factor(genotype, levels = geno_order))]
    comparison_p <- reviewed_p(metric, "Total")
    set.seed(FIG_SEED)
    ggplot(plot_data, aes(x = genotype, y = .data[[yvar]], fill = genotype)) + fig_summary(p_value = comparison_p) + fig_points(color = FIG_PT_COLOR) +
        fig_sig_bracket(comparison_p, plot_data[[yvar]]) + fig_scale_x_group(genotype_labels) + fig_scale_fill(geno_order,
        labels = genotype_labels, guide = "none") + fig_scale_y(ylab, limits = NULL, breaks = if(yvar=="percent_asleep") seq(0,100,25) else waiver()) + labs(x = NULL, title = stringr::str_wrap(title,
        width = 60), caption = fly_n_caption(plot_data)) + theme_fig() + coord_cartesian(ylim = ylim, clip = "off")
}
phase_pvals <- function(metric, dodge_width = FIG_DODGE_WIDTH) {
    data.frame(phase = c("Day", "Night"), x1 = (1:2) - dodge_width/4,
               x2 = (1:2) + dodge_width/4,
               p_value = vapply(c("Day", "Night"), function(ph) reviewed_p(metric, ph), numeric(1)))
}
phase_boxplot <- function(summary_dt, yvar, title, ylab, metric, ylim = NULL) {
    plot_data <- as.data.table(summary_dt)[genotype %in% geno_order]
    plot_data[, `:=`(genotype, factor(genotype, levels = geno_order))]
    comparisons <- phase_pvals(metric)
    set.seed(FIG_SEED)
    ggplot(plot_data, aes(x = phase, y = .data[[yvar]], fill = genotype)) + fig_summary(width = 0.70, position = position_dodge(width = FIG_DODGE_WIDTH),
        show_legend = TRUE) + geom_point(shape = FIG_PT_SHAPE, size = FIG_PT_SIZE, alpha = FIG_PT_ALPHA, color = FIG_PT_COLOR,
        show.legend = FALSE, position = position_jitterdodge(jitter.width = FIG_JITTER_WIDTH, jitter.height = 0, dodge.width = FIG_DODGE_WIDTH,
            seed = FIG_SEED)) + fig_sig_brackets_at(comparisons, plot_data[[yvar]]) + fig_scale_fill(geno_order, labels = genotype_labels,
        guide = guide_legend(title = NULL)) + fig_scale_y(ylab, limits = NULL, breaks = if(yvar=="percent_asleep") seq(0,100,25) else waiver()) + scale_x_discrete(name = NULL, expand = expansion(add = FIG_X_EXPAND)) +
        labs(title = stringr::str_wrap(title, width = 40), fill = NULL, caption = fly_n_caption(plot_data, "phase")) + theme_fig(legend_position = "top") +
        theme(legend.direction = "horizontal", legend.justification = "center") + coord_cartesian(ylim = ylim, clip = "off")
}
build_profile <- function(prof, ylab, title, legend_position = "top") {
    prof[, `:=`(genotype, factor(genotype, levels = geno_order))]
    ymax <- max(prof$m + prof$sem, na.rm = TRUE) * 1.05
    bar_lo <- ymax * 1.02
    bar_hi <- ymax * 1.08
    ld <- data.frame(xmin = c(0, 12, 24, 36), xmax = c(12, 24, 36, 48), phase = c("on", "off", "on", "off"))
    ggplot() + geom_rect(data = ld, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax), fill = phase_fill[ld$phase], alpha = 0.5) +
        geom_ribbon(data = prof, aes(x = zbin, ymin = m - sem, ymax = m + sem, fill = genotype), alpha = FIG_RIBBON_ALPHA) +
        geom_rect(data = ld, aes(xmin = xmin, xmax = xmax, ymin = bar_lo, ymax = bar_hi), fill = phase_bar[ld$phase], colour = "black",
            linewidth = FIG_LINE_WIDTH) + geom_line(data = prof, aes(x = zbin, y = m, colour = genotype), linewidth = FIG_PROFILE_LINEWIDTH) +
        scale_fill_manual(values = fig_colors(geno_order), guide = "none") + fig_scale_color(geno_order, labels = genotype_labels) +
        scale_x_continuous(name = FIG_LABEL_TIME, breaks = seq(0, 48, 6), limits = c(0, 48), expand = c(0, 0)) + scale_y_continuous(name = ylab,
        expand = c(0, 0)) + coord_cartesian(ylim = c(0, bar_hi), clip = "off") + labs(title = stringr::str_wrap(title, width = 60),
        caption = "Mean and SEM across individual flies; elapsed hours 12–60.") + theme_fig(legend_position = legend_position) + theme(legend.direction = "horizontal",
        legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA))
}

for (metric in c('sleep','activity')) {
    prof <- fread(file=file.path(FIG_DIR,paste0(metric,'_profile_summary.csv')))
    setnames(prof,c('mean','SEM','bin'),c('m','sem','zbin'))
    ylab <- if(metric=='sleep') FIG_LABEL_SLEEP_PROFILE else FIG_LABEL_ACTIVITY
    p <- build_profile(prof,ylab=ylab,title='',legend_position='top')
    fig_save_panel(p,file.path(FIG_DIR,paste0(metric,'_profile_corrected')),paste0(metric,'_profile'))
    total <- fread(file=file.path(FIG_DIR,paste0(metric,'_total_values_and_exclusions.csv')))[exclude==FALSE]
    phase <- fread(file=file.path(FIG_DIR,paste0(metric,'_phase_values_and_exclusions.csv')))[exclude==FALSE]
    phase[,phase:=factor(phase,levels=c('Day','Night'),labels=FIG_LABEL_PHASES)]
    yvar <- if(metric=='sleep') 'percent_asleep' else 'activity'
    ylab <- if(metric=='sleep') FIG_LABEL_SLEEP else FIG_LABEL_ACTIVITY
    p <- two_group_boxplot(total,yvar,'',ylab,metric,ylim=if(metric=='sleep') c(0,100) else c(0,NA))
    fig_save_panel(p,file.path(FIG_DIR,paste0(metric,'_total_corrected')),paste0(metric,'_total'))
    # Leave room for the adjusted-P annotation above the full observed range.
    p <- phase_boxplot(phase,yvar,'',ylab,metric,ylim=if(metric=='sleep') c(0,115) else c(0,NA))
    fig_save_panel(p,file.path(FIG_DIR,paste0(metric,'_phase_corrected')),paste0(metric,'_phase'))
}
