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

transfection <- read.csv(
  file.path(figure_dir, "Supporting_Data", "transfection_efficiency_plot_data.csv"),
  stringsAsFactors = FALSE
)
transfection$group <- factor(transfection$group, levels = c("WT", "CS"))

p <- ggplot(transfection, aes(group, efficiency, fill = group, color = group)) +
  fig_reference_line(50) +
  fig_summary(p_value = NA_real_) +
  fig_points() +
  fig_scale_fill(c("WT", "CS")) +
  fig_scale_color(c("WT", "CS"), guide = "none") +
  fig_scale_x_group(c(WT = "MIC60 WT", CS = "MIC60 CS")) +
  fig_scale_y(
    "Transfection efficiency (%)",
    breaks = seq(0, 100, 25)
  ) +
  labs(
    x = NULL,
    title = "Transfection efficiency",
    subtitle = "Each point is one image; 2 independent experiments",
    caption = paste0(
      "n = 10 images/genotype; 2 independent experiments.\n",
      "Points: images; boxes: median and IQR; dashed line: 50%."
    )
  ) +
  theme_fig()

fig_save(
  p,
  file.path(figure_dir, "Rebuilt_Output", "Figure_S3C_Transfection_Efficiency"),
  width = FIG_W_1COL,
  height = FIG_H_1COL,
  formats = c("png", "pdf")
)

cat("WT mean:", mean(transfection$efficiency[transfection$group == "WT"]), "\n")
cat("CS mean:", mean(transfection$efficiency[transfection$group == "CS"]), "\n")
