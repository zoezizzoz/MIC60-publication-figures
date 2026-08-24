# Publication analysis of TMRM/MTG in dMIC60CS versus WT flight muscle.
# Uses the single acquisition setting shared by both genotypes and the
# project-wide formatting functions in figure_style.R.
#
# Analysis unit: image field (not assumed here to be a biological replicate).
# Display: image fields as dots, median and IQR as boxes, group means as
# white diamonds. Statistical comparison: two-sided Mann-Whitney test.

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required package(s) first: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
project_dir <- normalizePath(file.path(figure_dir, ".."), mustWork = TRUE)
style_file <- file.path(script_dir, "figure_style.R")
source(style_file)

measurement_path <- file.path(figure_dir, "Original_Data", "tmrm_raw_measurements.csv")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

matched_setting <- "905_-99 1143_-51"
repeat_field_exclusions <- c(
  "TMRM MTG CS WR ZZ.lif - 7.16.26 CS_Series024_2.tif"
)

raw <- read.csv(measurement_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "file", "condition", "settings_group", "tmrm_mtg_ratio_raw",
  "primary_analysis_include"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing required column(s): ", paste(missing_columns, collapse = ", "))
}

dat <- raw[
  raw$primary_analysis_include == "Yes" &
    raw$condition %in% c("WR", "CS") &
    raw$settings_group == matched_setting &
    !raw$file %in% repeat_field_exclusions &
    is.finite(raw$tmrm_mtg_ratio_raw),
  required_columns
]

if (nrow(dat) == 0) {
  stop("No primary-analysis images were found.")
}

if (!identical(unique(dat$settings_group), matched_setting)) {
  stop("The retained data do not match the specified acquisition setting.")
}

dat$group <- ifelse(dat$condition == "CS", "CSF", "WR")
dat$group <- factor(dat$group, levels = c("WR", "CSF"))
dat$image <- dat$file

wr_mean <- mean(
  dat$tmrm_mtg_ratio_raw[dat$group == "WR"],
  na.rm = TRUE
)
if (!is.finite(wr_mean) || wr_mean <= 0) {
  stop("The WR reference mean is invalid.")
}

dat$percent_wr <- dat$tmrm_mtg_ratio_raw / wr_mean * 100

group_means <- aggregate(
  percent_wr ~ group,
  data = dat,
  FUN = mean
)

if (sum(dat$group == "WR") < 2 || sum(dat$group == "CSF") < 2) {
  stop("At least two image-level replicates are required in both groups.")
}

has_ties <- anyDuplicated(dat$percent_wr) > 0
mw <- wilcox.test(
  percent_wr ~ group,
  data = dat,
  alternative = "two.sided",
  exact = !has_ties,
  correct = has_ties,
  conf.int = FALSE
)

group_order <- c("WR", "CSF")
group_labels <- c(
  "WR" = "dMIC60 WT\n\u2640",
  "CSF" = "dMIC60 CS\n\u2640"
)

set.seed(FIG_SEED)

p <- ggplot(
  dat,
  aes(x = group, y = percent_wr, fill = group, color = group)
) +
  fig_reference_line(100) +
  fig_summary(p_value = mw$p.value) +
  fig_points(dense = FALSE) +
  fig_markers(
    data = group_means,
    mapping = aes(x = group, y = percent_wr),
    width = 0
  ) +
  fig_sig_bracket(mw$p.value, dat$percent_wr) +
  fig_scale_fill(group_order) +
  fig_scale_color(group_order, guide = "none") +
  fig_scale_x_group(group_labels) +
  fig_scale_y("TMRM / MTG (% of WT image mean)") +
  labs(
    x = NULL,
    title = paste0(
      "Mitochondrial membrane potential\n",
      "in adult female Drosophila indirect flight muscle"
    ),
    subtitle = "Image-level measurements",
    caption = paste0(
      "n = ", sum(dat$group == "WR"), " WT and ",
      sum(dat$group == "CSF"), " CS images (n denotes images).\n",
      "Normalized to the mean WT TMRM/MTG ratio.\n",
      "Boxes: median/IQR; diamonds: means.\n",
      "Two-sided Mann-Whitney test.\n",
      "Settings: TMRM 905/-99; MTG 1143/-51.\n",
      "Excluded: 3 saturated and 1 repeated CS image."
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_fig()

output_base <- file.path(
  output_dir,
  "TMRM_MTG_CSF_vs_WR_matched_setting_image_level"
)
fig_save(
  p,
  output_base,
  width = FIG_W_1COL,
  height = FIG_H_1COL + 0.7
)

write.csv(
  dat,
  paste0(output_base, "_processed_values.csv"),
  row.names = FALSE
)

statistics <- data.frame(
  comparison = "dMIC60CS vs WT",
  analysis_unit = "image",
  n_WT_images = sum(dat$group == "WR"),
  n_CS_images = sum(dat$group == "CSF"),
  WT_mean_percent = mean(dat$percent_wr[dat$group == "WR"]),
  CS_mean_percent = mean(dat$percent_wr[dat$group == "CSF"]),
  test = "two-sided Mann-Whitney",
  p_value = mw$p.value,
  matched_setting = matched_setting,
  stringsAsFactors = FALSE
)
write.csv(
  statistics,
  paste0(output_base, "_statistics.csv"),
  row.names = FALSE
)

figure_legend <- paste0(
  "Fig. 3O. Mitochondrial membrane potential in adult female Drosophila ",
  "indirect flight muscle. TMRM intensity was normalized to MitoTracker ",
  "Green (MTG) intensity for each image and expressed relative to the mean ",
  "WT TMRM/MTG ratio. Each point represents one image; boxes show the median ",
  "and interquartile range, and white diamonds show group means. n = ",
  sum(dat$group == "WR"), " WT and ", sum(dat$group == "CSF"),
  " dMIC60CS images (n denotes images). P = ",
  format.pval(mw$p.value, digits = 2, eps = 0.0001),
  ", two-sided Mann-Whitney test. Images were acquired using identical ",
  "settings (TMRM gain/offset, 905/-99; MTG gain/offset, 1143/-51). Three ",
  "saturated dMIC60CS images and one repeated dMIC60CS field were excluded ",
  "before analysis."
)
writeLines(
  figure_legend,
  paste0(output_base, "_figure_legend.txt")
)

cat("\nMatched setting: ", unique(dat$settings_group), "\n", sep = "")
cat("WR images: ", sum(dat$group == "WR"), "\n", sep = "")
cat("CS images: ", sum(dat$group == "CSF"), "\n", sep = "")
cat("WR reference mean TMRM/MTG: ", sprintf("%.6f", wr_mean), "\n", sep = "")
cat("CS mean (% WR): ", sprintf("%.2f", mean(dat$percent_wr[dat$group == "CSF"])), "\n", sep = "")
cat("Mann-Whitney, two-sided p = ", format.pval(mw$p.value, digits = 3), "\n", sep = "")
cat("Saved: ", paste0(output_base, ".pdf"), "\n", sep = "")
cat("Saved: ", paste0(output_base, ".png"), "\n", sep = "")
