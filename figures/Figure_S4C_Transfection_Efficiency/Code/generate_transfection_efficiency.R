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
package_dir <- normalizePath(file.path(figure_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
dir.create(file.path(figure_dir, "Rebuilt_Output"), recursive = TRUE, showWarnings = FALSE)
TE_AXIS_TITLE_SIZE <- 16
TE_AXIS_TEXT_SIZE <- 14
TE_GROUP_LABEL_SIZE <- 16
FIG_COLORS[c("control", "mutant")] <- c("#5AB3E5", "#CA79A6")
FIG_FILL_LIGHTEN <- 0
transfection <- read.csv(file.path(figure_dir, "Supporting_Data", "transfection_efficiency_plot_data.csv"), stringsAsFactors = FALSE)
transfection$group <- factor(transfection$group, levels = c("WT", "CS"))
p <- ggplot(transfection, aes(group, efficiency, fill = group, color = group)) + fig_reference_line(50) + fig_summary(p_value = NA_real_) +
    fig_points(color = "black") + fig_scale_fill(c("WT", "CS")) + fig_scale_color(c("WT", "CS"), guide = "none") + fig_scale_x_group(c(WT = "dMIC60-WT-FLAG",
    CS = "dMIC60-CS-FLAG")) + fig_scale_y("Transfection efficiency (%)", breaks = seq(0, 100, 25)) + labs(x = NULL, title = "Transfection efficiency",
    subtitle = "Each point is one image; 2 independent experiments", caption = paste0("n = 10 images/genotype; 2 independent experiments.\n",
        "Points: images; boxes: median and IQR; dashed line: 50%.")) + theme_fig() + theme(axis.title = element_text(size = TE_AXIS_TITLE_SIZE),
    axis.text = element_text(size = TE_AXIS_TEXT_SIZE), axis.text.x = element_text(size = TE_GROUP_LABEL_SIZE))
fig_save(p, file.path(figure_dir, "Rebuilt_Output", "FigS4C_transfection_efficiency"), width = FIG_W_WIDE, height = FIG_H_1COL,
    formats = c("png", "pdf"))
cat("WT mean:", mean(transfection$efficiency[transfection$group == "WT"]), "\n")
cat("CS mean:", mean(transfection$efficiency[transfection$group == "CS"]), "\n")
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
