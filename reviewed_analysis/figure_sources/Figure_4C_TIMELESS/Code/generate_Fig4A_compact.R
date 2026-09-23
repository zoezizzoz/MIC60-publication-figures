#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- grep("^--file=",gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE),value=TRUE)
project_dir <- dirname(dirname(normalizePath(sub("^--file=","",args[1]))))
output_dir <- file.path(project_dir,"Rebuilt_Output")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
source(file.path(project_dir, "Code", "figure_style.R"))
FIG_ROLE_OF <- c(
  FIG_ROLE_OF,
  "WT +H2O2" = "control_male",
  "CS +H2O2" = "mutant_male"
)

# Compact, presentation-style proportions based on the supplied example.
FIG_BASE_SIZE <- 15
FIG_AXIS_TITLE_SIZE <- 15
FIG_AXIS_TEXT_SIZE <- 14
FIG_ANNOT_SIZE <- 4.5
FIG_SUMMARY_WIDTH <- 0.65
# One shared discrete interval is used across all four groups; the larger,
# symmetric outer expansion keeps those equal intervals compact on the canvas.
FIG_X_EXPAND <- 1.65
FIG_MARGIN_PT <- c(t = 5, r = 8, b = 5, l = 6)

groups <- c("WT", "CS", "WT +H2O2", "CS +H2O2")
display_labels <- c(
  "WT\n−",
  "CS\n−",
  "WT\n+",
  "CS\n+"
)

dat <- read.csv(
  file.path(project_dir, "Supporting_Data", "Fig4A_TIMELESS_DAPI_plotted_values.csv"),
  stringsAsFactors = FALSE
)
dat$group <- factor(dat$group, levels = groups)

p_reference <- c(untreated = 0.0477, h2o2 = 0.0218)

p <- ggplot(dat, aes(x = group, y = value, fill = group)) +
  fig_summary(style = "box") +
  fig_points(width = 0.08, color = FIG_PT_COLOR) +
  fig_sig_bracket(
    p_value = p_reference[["untreated"]],
    values = dat$value[dat$group %in% c("WT", "CS")],
    x1 = 1,
    x2 = 2,
    label = expression(italic(p) == 0.0477)
  ) +
  fig_sig_bracket(
    p_value = p_reference[["h2o2"]],
    values = dat$value[dat$group %in% c("WT +H2O2", "CS +H2O2")],
    x1 = 3,
    x2 = 4,
    label = expression(italic(p) == 0.0218)
  ) +
  fig_scale_fill(groups) +
  fig_scale_x_group(labels = display_labels) +
  fig_scale_y(name = "TIMELESS / DAPI", limits = c(0, NA)) +
  labs(x = NULL) +
  annotation_custom(
    grid::textGrob("dMIC60", gp = grid::gpar(fontsize = 12.5)),
    xmin = -1.20, xmax = -1.20, ymin = -0.070, ymax = -0.070
  ) +
  annotation_custom(
    grid::textGrob(expression(H[2]*O[2]), gp = grid::gpar(fontsize = 12.5)),
    xmin = -1.20, xmax = -1.20, ymin = -0.155, ymax = -0.155
  ) +
  guides(fill = "none") +
  theme_fig() +
  theme(
    axis.text.x = element_text(size = 12.5, lineheight = 0.92),
    plot.margin = margin(
      FIG_MARGIN_PT[["t"]], FIG_MARGIN_PT[["r"]],
      FIG_MARGIN_PT[["b"]], FIG_MARGIN_PT[["l"]], unit = "pt"
    )
  ) +
  coord_cartesian(clip = "off")

fig_save(
  p,
  file.path(output_dir, "Fig4A_TIMELESS_DAPI_provided_data_compact"),
  width = 5.4,
  height = 3.8
)

message("Created compact source-backed Fig4A TIMELESS/DAPI plot.")

writeLines(capture.output(sessionInfo()),file.path(output_dir,"R_sessionInfo.txt"))
