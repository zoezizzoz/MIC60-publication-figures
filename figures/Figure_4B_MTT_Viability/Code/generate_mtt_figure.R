# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
    normalizePath(getwd())
}
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
dir.create(file.path(figure_dir, "Rebuilt_Output"), recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(figure_dir, "Supporting_Data", "mtt_plot_values.csv"), stringsAsFactors = FALSE)
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
sig_tests <- read.csv(file.path(figure_dir, "Supporting_Data", "mtt_welch_tests_rebuilt_TECHNICAL_WELLS_NOT_FOR_INFERENCE.csv"),
    stringsAsFactors = FALSE)
sig_tests <- sig_tests[is.finite(sig_tests$p_value) & sig_tests$p_value < FIG_ALPHA, ]
sig_tests$x1 <- sig_tests$dose_mM - 1.05
sig_tests$x2 <- sig_tests$dose_mM + 1.05
sig_tests$label_x <- ifelse(sig_tests$dose_mM == 5, 7, sig_tests$dose_mM)
sig_tests$y <- vapply(sig_tests$dose_mM, function(dose) max(d$normalized_pct[d$dose_mM == dose], na.rm = TRUE) + 5, numeric(1))
sig_tests$label <- vapply(sig_tests$p_value, function(p_value) {
    if (p_value < 0.001)
        "p < 0.001"
    else paste0("p = ", format.pval(p_value, digits = 2))
}, character(1))
p <- ggplot(d, aes(x_plot, normalized_pct)) + fig_reference_line(100) + geom_point(color = FIG_PT_COLOR, position = "identity",
    shape = FIG_PT_SHAPE, size = 1.6, alpha = 0.72, show.legend = FALSE) + geom_line(data = summ, aes(x = x_plot, y = mean,
    group = genotype, color = genotype), linewidth = FIG_PROFILE_LINEWIDTH) + fig_errorbar(data = summ, mapping = aes(x = x_plot,
    y = mean, ymin = mean - sd, ymax = mean + sd, group = genotype), width = 1.8) + geom_point(data = experiment_summ, aes(x = x_plot,
    y = experiment_mean), color = FIG_PT_COLOR, shape = FIG_PT_SHAPE, size = 2.5, alpha = 0.95, show.legend = FALSE) + geom_segment(data = sig_tests,
    aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH, color = "black") + geom_segment(data = sig_tests,
    aes(x = x1, xend = x1, y = y, yend = y - 1.2), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH, color = "black") +
    geom_segment(data = sig_tests, aes(x = x2, xend = x2, y = y, yend = y - 1.2), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH,
        color = "black") + geom_text(data = sig_tests, aes(x = label_x, y = y + 1.5, label = label), inherit.aes = FALSE,
    vjust = 0, size = FIG_ANNOT_SIZE, family = FIG_FONT, color = "black") + scale_x_continuous(breaks = doses, labels = doses,
    expand = expansion(add = 2.2)) + scale_color_manual(name = expression(italic("MIC60 Null") ~ "HeLa Cells"), values = c(WT = "#0072B2",
    CS = "#B83B6B"), limits = c("WT", "CS"), labels = c("dMIC60-WT", "dMIC60-CS")) + fig_scale_y("Cell viability (% of matched 0 mM control)",
    limits = c(0, 130), breaks = seq(0, 125, 25)) + labs(x = "Hydrogen peroxide concentration (mM)", title = "MTT dose response",
    caption = paste0("Technical wells per genotype: n = 9, 6, 6, 9, and 9 at 0, 5, 10, 20, and 40 mM H₂O₂, respectively.\n",
        "Data were pooled from two experiments (biological n = 2 at 0/20/40 mM; n = 1 at 5/10 mM).\n", "Small black points: wells; larger black points: experiment means.\n",
        "Colored lines: mean of experiment means; error bars: SD across experiments (n = 2 only).\n", "Brackets show significant two-sided Welch tests of technical wells; these p values are descriptive.")) +
    guides(color = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.35))) +
    theme_fig(legend_position = "inside") + theme(legend.position.inside = c(0.5, 0.84), legend.justification.inside = c(0.5,
    0), legend.title = element_text(size = 14.5, color = "black", hjust = 0.5, margin = margin(b = 1)), legend.text = element_text(size = 14,
    color = "black"), legend.key.width = grid::unit(24, "pt"), legend.key.height = grid::unit(10, "pt"), legend.spacing.x = grid::unit(6,
    "pt"), legend.spacing.y = grid::unit(0, "pt"), legend.margin = margin(t = 0, r = 0, b = 0, l = 0), legend.box.spacing = grid::unit(0,
    "pt"), plot.margin = margin(t = 0, r = 10, b = 8, l = 8))
fig_save(p, file.path(figure_dir, "Rebuilt_Output", "MTT_dose_response_combined"), width = FIG_W_WIDE, height = FIG_H_WIDE,
    formats = c("png", "pdf"))
write.csv(experiment_summ[, c("experiment", "dose_mM", "genotype", "experiment_mean", "technical_sd", "technical_n")], .rebuild_file(file.path(figure_dir,
    "Rebuilt_Output", "mtt_experiment_means.csv")), row.names = FALSE)
write.csv(summ[, c("dose_mM", "genotype", "mean", "sd", "biological_n", "technical_n")], .rebuild_file(file.path(figure_dir,
    "Rebuilt_Output", "mtt_biological_summary.csv")), row.names = FALSE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
