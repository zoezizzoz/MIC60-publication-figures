# Descriptive MTT: experiment means and between-experiment SD; no technical-well significance brackets.
.label_script <- normalizePath(gsub('~+~',' ',sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1]),fixed=TRUE))
source(file.path(dirname(dirname(dirname(.label_script))), 'Graph_Label_Conventions.R'))
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
# Fig. 4B is reproduced at a small size in the assembled figure, so use a
# substantially larger type scale for legibility at final placement size.
FIG_BASE_SIZE <- 11
FIG_AXIS_TITLE_SIZE <- 12
FIG_AXIS_TEXT_SIZE <- 11
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
p <- ggplot(d, aes(x_plot, normalized_pct)) + fig_reference_line(100) + geom_point(color = FIG_PT_COLOR, position = "identity", 
    shape = FIG_PT_SHAPE, size = 1.6, alpha = 0.72, show.legend = FALSE) + geom_line(data = summ, aes(x = x_plot, y = mean, 
    group = genotype, color = genotype), linewidth = FIG_PROFILE_LINEWIDTH) + fig_errorbar(data = summ, mapping = aes(x = x_plot, 
    y = mean, ymin = mean - sd, ymax = mean + sd, group = genotype), width = 1.8, na.rm = TRUE) + geom_point(data = experiment_summ, aes(x = x_plot, 
    y = experiment_mean, color = genotype, shape = factor(experiment)), size = 2.5, alpha = 0.95) + scale_x_continuous(breaks = doses, labels = doses, 
    expand = expansion(add = 2.2)) + scale_color_manual(name = NULL, values = c(WT = "#5AB4E5", 
    CS = "#CB78A8"), limits = c("WT", "CS"), labels = c("dMIC60-WT", "dMIC60-CS")) + fig_scale_y("Metabolic activity (MTT; % of matched 0 mM control)", 
    limits = c(0, 130), breaks = seq(0, 125, 25)) + labs(x = expression(H[2]*O[2]~"concentration (mM)"), title = "MTT dose response", 
    caption = paste0("Technical wells per genotype: n = 9, 6, 6, 9, and 9 at 0, 5, 10, 20, and 40 mM H₂O₂, respectively.\n", 
        "Data were pooled from two experiments (biological n = 2 at 0/20/40 mM; n = 1 at 5/10 mM).\n", "Small black points: wells; larger colored symbols: experiment means.\n", 
        "Colored lines: mean of experiment means; error bars: SD across experiments (n = 2 only).\n", "No inferential tests: technical wells are subsamples; 1–2 independent experiments per dose.")) + 
    scale_shape_manual(name = "Experiment", values = c(16,17), labels = c("Experiment 1","Experiment 2")) +
    guides(shape=guide_legend(order=2), color = guide_legend(order=1, title.position = "top", title.hjust = 0.5, nrow = 1, byrow = TRUE, override.aes = list(linewidth = 1.35))) + 
    theme_fig(legend_position = "top") + theme(legend.direction="horizontal", legend.justification="center",
      legend.text=element_text(size=FIG_AXIS_TEXT_SIZE), legend.box="vertical",
      # Compact legend rows and axis-title gaps without changing label sizes.
      legend.spacing.y=grid::unit(1,"pt"), legend.box.spacing=grid::unit(2,"pt"),
      legend.margin=margin(0,0,0,0), legend.box.margin=margin(30,0,-30,0),
      axis.title.y=element_text(margin=margin(r=3)),
      axis.title.x=element_text(margin=margin(t=3)),
      plot.margin=margin(4,10,8,8))

fig_save(p, file.path(figure_dir, "Rebuilt_Output", "MTT_dose_response_combined"), width = FIG_W_WIDE, height = FIG_H_WIDE, 
    formats = c("png", "pdf"))
write.csv(experiment_summ[, c("experiment", "dose_mM", "genotype", "experiment_mean", "technical_sd", "technical_n")], .rebuild_file(file.path(figure_dir, 
    "Rebuilt_Output", "mtt_experiment_means.csv")), row.names = FALSE)
write.csv(summ[, c("dose_mM", "genotype", "mean", "sd", "biological_n", "technical_n")], .rebuild_file(file.path(figure_dir, 
    "Rebuilt_Output", "mtt_biological_summary.csv")), row.names = FALSE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
