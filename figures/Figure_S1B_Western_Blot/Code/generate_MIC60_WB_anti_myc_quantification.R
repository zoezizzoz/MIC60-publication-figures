#!/usr/bin/env Rscript

# Rebuild the publication-ready MIC60-Myc protein abundance figure.
# Run from any working directory:
#   Rscript Code/generate_MIC60_WB_anti_myc_quantification.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
})

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Run this file with Rscript so its package-relative paths can be resolved.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
package_dir <- dirname(script_dir)
data_file <- file.path(package_dir, "Original_Data", "myc WB quantifications.xlsx")
output_dir <- file.path(package_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "figure_style.R"))
FIG_SHOW_NS <- TRUE

# This panel is reduced substantially when placed beside the western blot in
# Figure S1B, so use larger plot-specific typography for readability at the
# final assembled-figure size. These overrides do not affect other figures.
WB_AXIS_TITLE_SIZE <- 21
WB_AXIS_TEXT_SIZE <- 19
WB_ANNOT_SIZE <- 5.8
FIG_ANNOT_SIZE <- WB_ANNOT_SIZE

# Plot-specific fills sampled from the approved Figure S1B reference image.
FIG_COLORS[c("control", "mutant")] <- c("#5AB3E5", "#CA79A6")
FIG_FILL_LIGHTEN <- 0

wide <- read_excel(data_file, sheet = "Sheet1", .name_repair = "minimal")
date_columns <- names(wide)[2:4]
blot_dates <- as.Date(as.numeric(date_columns), origin = "1899-12-30")

extract_group <- function(group_name) {
  row_index <- which(wide[[1]] == group_name)
  if (length(row_index) != 1L) {
    stop("Expected exactly one row for ", group_name, ".")
  }
  as.numeric(wide[row_index, 2:4])
}

wt <- extract_group("WT Rescue")
cs <- extract_group("CS Rescue")

plot_data <- data.frame(
  blot_date = rep(format(blot_dates, "%Y-%m-%d"), each = 2),
  genotype = factor(
    rep(c("dMIC60WT", "dMIC60CS"), times = length(blot_dates)),
    levels = c("dMIC60WT", "dMIC60CS")
  ),
  normalized_MIC60_Myc_over_ATP5B = as.vector(rbind(wt, cs)),
  larvae_per_replicate = 5L
)

# Compact x positions keep the two genotypes close together while preserving
# the original group order and all plotted values.
x_positions <- c(dMIC60WT = 1.00, dMIC60CS = 1.48)
plot_data$x_position <- unname(x_positions[as.character(plot_data$genotype)])

paired_test <- t.test(wt, cs, paired = TRUE, alternative = "two.sided")
p_label <- sprintf("n.s.\np = %.3f", paired_test$p.value)

figure <- ggplot(
  plot_data,
  aes(
    x = x_position,
    y = normalized_MIC60_Myc_over_ATP5B,
    group = genotype,
    fill = genotype
  )
) +
  fig_summary(p_value = paired_test$p.value, width = 0.28) +
  fig_points(color = "black", width = 0.035) +
  fig_sig_bracket(
    p_value = paired_test$p.value,
    values = plot_data$normalized_MIC60_Myc_over_ATP5B,
    x1 = x_positions[["dMIC60WT"]],
    x2 = x_positions[["dMIC60CS"]],
    label = p_label
  ) +
  fig_scale_fill(c("dMIC60WT", "dMIC60CS")) +
  scale_x_continuous(
    breaks = unname(x_positions),
    labels = c("WT", "CS"),
    limits = c(0.68, 1.80),
    expand = expansion(mult = 0)
  ) +
  fig_scale_y(
    "MIC60-Myc / ATP5B\n(normalized intensity)",
    limits = c(0, NA),
    breaks = seq(0, 2.0, 0.5)
  ) +
  labs(
    x = NULL,
    title = "MIC60-Myc protein abundance",
    subtitle = "Three matched blots; 5 larvae pooled per replicate",
    caption = paste0(
      "n = 3 biological replicates per genotype; each point represents one blot.\n",
      "Each replicate contains lysate pooled from 5 larvae.\n",
      "Boxes: median/IQR.\n",
      "Whiskers: 1.5 x IQR; two-sided paired t-test."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_fig() +
  theme(
    axis.title = element_text(size = WB_AXIS_TITLE_SIZE),
    axis.text = element_text(size = WB_AXIS_TEXT_SIZE),
    plot.margin = margin(t = 10, r = 9, b = 9, l = 12)
  )

output_base <- file.path(output_dir, "MIC60_WB_anti_myc_quantification")
fig_save(
  figure,
  output_base,
  width = FIG_W_1COL,
  height = FIG_H_1COL + 0.35,
  formats = c("png", "pdf")
)

write.csv(
  plot_data,
  file.path(output_dir, "MIC60_WB_graph_ready_data.csv"),
  row.names = FALSE
)

summary_statistics <- data.frame(
  comparison = "dMIC60WT vs dMIC60CS",
  n_biological_replicates_per_genotype = length(wt),
  larvae_pooled_per_replicate = 5L,
  dMIC60WT_mean = mean(wt),
  dMIC60WT_median = median(wt),
  dMIC60WT_sd = sd(wt),
  dMIC60CS_mean = mean(cs),
  dMIC60CS_median = median(cs),
  dMIC60CS_sd = sd(cs),
  paired_t_statistic = unname(paired_test$statistic),
  paired_t_df = unname(paired_test$parameter),
  paired_t_p_value = paired_test$p.value,
  stringsAsFactors = FALSE
)

write.csv(
  summary_statistics,
  file.path(output_dir, "MIC60_WB_summary_statistics.csv"),
  row.names = FALSE
)

message("Wrote publication outputs to: ", output_dir)
