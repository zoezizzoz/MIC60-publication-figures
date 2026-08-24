# Audited TEM mitochondrial morphology analysis
# Female dMIC60WT (WR) versus dMIC60CS (CS), 5,000x images
#
# Primary panels use one point per TEM field (image-level mean).
# Mitochondria within a field are subsamples. Because fly/specimen IDs are not
# available, inferential p-values are exploratory and must not be described as
# biological-replicate statistics.

required_packages <- c("ggplot2", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

library(ggplot2)
library(readxl)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
base_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
# Never annotate non-significant comparisons in these TEM graphs.
FIG_SHOW_NS <- FALSE
# Use the Unicode-capable project font and Quartz PDF device so the female
# symbol is preserved in both raster and vector exports on macOS.
FIG_FONT <- "Arial Unicode MS"
quartz_pdf_unicode <- function(filename, width, height, ...) {
  grDevices::quartz(
    type = "pdf",
    file = filename,
    width = width,
    height = height,
    family = FIG_FONT,
    ...
  )
}
FIG_PDF_DEVICE <- quartz_pdf_unicode
input_file <- file.path(base_dir, "Original_Data", "7.28.26 TEM mito analysis raw data.xlsx")
output_dir <- file.path(base_dir, "Rebuilt_Output")
inventory_file <- file.path(base_dir, "Supporting_Data", "TEM_image_inventory.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- as.data.frame(read_excel(input_file, sheet = "Raw Data"))
names(raw) <- trimws(names(raw))

required_columns <- c(
  "Record ID", "Condition", "Image", "Area (µm²)", "Perimeter (µm)",
  "Major (µm)", "Minor (µm)", "Aspect ratio", "Circularity",
  "Roundness", "Solidity", "Feret (µm)"
)
if (length(setdiff(required_columns, names(raw)))) {
  stop("Input workbook does not contain all required columns.")
}

numeric_columns <- setdiff(required_columns, c("Condition", "Image"))
for (column_name in numeric_columns) {
  raw[[column_name]] <- suppressWarnings(as.numeric(raw[[column_name]]))
}

raw$Condition <- toupper(trimws(raw$Condition))
raw$Image <- trimws(raw$Image)
# Cnt is the WR/control genotype. Enforce the alias from either the condition
# field or the legacy Cnt filename before any grouping or statistics.
raw$Condition[raw$Condition == "CNT" | grepl("^Cnt ", raw$Image)] <- "WR"
raw$Source_Image <- raw$Image
raw$Calibration_Status <- "Embedded/calibrated µm values"
raw$Pixel_Size_um <- NA_real_
raw$QC_Status <- "INCLUDED"
raw$Exclusion_Reason <- ""

# Six repeated ImageJ measurements were identified by near-identical centroid,
# area, and major-axis values in the full 40-column source export. The later
# member of each pair is excluded; source rows remain in the QC table.
probable_duplicate_ids <- c(76, 99, 201, 223, 251, 262)
is_duplicate <- raw$`Record ID` %in% probable_duplicate_ids
raw$QC_Status[is_duplicate] <- "EXCLUDED"
raw$Exclusion_Reason[is_duplicate] <- "Probable repeated ImageJ measurement"

# The four older WR-control TIFFs lost spatial calibration on export. Their
# pixel dimensions are converted using the visible 2 µm scale-bar lengths in
# the original TIFFs. The two pixel-unit rows for image 0008 duplicate the two
# already-calibrated rows and are excluded in favor of the calibrated rows.
scale_bar_pixels <- c(
  "Cnt fmale_a_5000X_0005-1.tif" = 360,
  "Cnt fmale_a_5000X_0008-1.tif" = 358,
  "Cnt fmale_a_5000X_0010-1.tif" = 359,
  "Cnt fmale_a_5000X_0012-1.tif" = 360
)

is_cnt_pixel <- raw$Image %in% names(scale_bar_pixels) & raw$`Area (µm²)` > 100
raw$Pixel_Size_um[is_cnt_pixel] <- 2 / scale_bar_pixels[raw$Image[is_cnt_pixel]]
raw$Calibration_Status[is_cnt_pixel] <- "Converted from pixels using visible 2 µm scale bar"

raw$`Area (µm²)`[is_cnt_pixel] <- (
  raw$`Area (µm²)`[is_cnt_pixel] * raw$Pixel_Size_um[is_cnt_pixel]^2
)
for (column_name in c("Perimeter (µm)", "Major (µm)", "Minor (µm)", "Feret (µm)")) {
  raw[[column_name]][is_cnt_pixel] <- (
    raw[[column_name]][is_cnt_pixel] * raw$Pixel_Size_um[is_cnt_pixel]
  )
}

cnt_0008_pixel_ids <- raw$`Record ID` %in% c(309, 310)
raw$QC_Status[cnt_0008_pixel_ids] <- "EXCLUDED"
raw$Exclusion_Reason[cnt_0008_pixel_ids] <- (
  "Same two objects remeasured in calibrated rows 311-312"
)

# Canonical source filenames for coverage matching.
raw$Source_Image <- sub("-1\\.tif$", ".tif", raw$Source_Image)
raw$Source_Image[raw$Image == "Cnt fmale_a_5000X_0008.dm3"] <- (
  "Cnt fmale_a_5000X_0008.tif"
)

finite_required <- is.finite(raw$`Area (µm²)`) &
  is.finite(raw$`Perimeter (µm)`) &
  raw$`Area (µm²)` > 0 & raw$`Perimeter (µm)` > 0
bad_measurement <- !finite_required & raw$QC_Status == "INCLUDED"
raw$QC_Status[bad_measurement] <- "EXCLUDED"
raw$Exclusion_Reason[bad_measurement] <- "Missing, non-finite, or non-physical morphology value"

analysis_dat <- raw[
  raw$QC_Status == "INCLUDED" & raw$Condition %in% c("WR", "CS"),
]

# Physical consistency checks (rounded ImageJ outputs allow small differences).
analysis_dat$AR_Recomputed <- analysis_dat$`Major (µm)` / analysis_dat$`Minor (µm)`
analysis_dat$AR_Absolute_Error <- abs(
  analysis_dat$`Aspect ratio` - analysis_dat$AR_Recomputed
)
analysis_dat$AR_QC <- ifelse(
  is.finite(analysis_dat$`Aspect ratio`) & analysis_dat$`Aspect ratio` >= 1,
  "AVAILABLE",
  "NOT AVAILABLE"
)
analysis_dat$Circularity_Recomputed <- (
  4 * pi * analysis_dat$`Area (µm²)` / analysis_dat$`Perimeter (µm)`^2
)
analysis_dat$Circularity_Absolute_Error <- abs(
  analysis_dat$Circularity - analysis_dat$Circularity_Recomputed
)

aggregate_mean <- function(value, condition, image) {
  aggregate(
    value,
    by = list(Condition = condition, Image = image),
    FUN = function(z) {
      z <- z[is.finite(z)]
      if (length(z)) mean(z) else NA_real_
    }
  )
}

area_summary <- aggregate_mean(
  analysis_dat$`Area (µm²)`, analysis_dat$Condition, analysis_dat$Source_Image
)
names(area_summary)[3] <- "Mean_Area_um2"

perimeter_summary <- aggregate_mean(
  analysis_dat$`Perimeter (µm)`, analysis_dat$Condition, analysis_dat$Source_Image
)
names(perimeter_summary)[3] <- "Mean_Perimeter_um"

ar_summary <- aggregate_mean(
  analysis_dat$`Aspect ratio`, analysis_dat$Condition, analysis_dat$Source_Image
)
names(ar_summary)[3] <- "Mean_Aspect_Ratio"

image_counts <- aggregate(
  rep.int(1L, nrow(analysis_dat)),
  by = list(Condition = analysis_dat$Condition, Image = analysis_dat$Source_Image),
  FUN = sum
)
names(image_counts)[3] <- "Measured_Mitochondria_n"

image_summary <- Reduce(
  function(x, y) merge(x, y, by = c("Condition", "Image"), all = TRUE),
  list(image_counts, area_summary, ar_summary, perimeter_summary)
)
image_summary$Acquisition_Subset <- ifelse(
  grepl("^Cnt ", image_summary$Image),
  "WR control TIFF set (2023-12-01)",
  ifelse(
    grepl("^WR-F", image_summary$Image),
    "WR-F DM3 set (2024-04-29)",
    ifelse(
      grepl("-00(27|28)\\.dm3$", image_summary$Image),
      "CS grid 1 (2024-04-11)",
      "CS grid 2 (2024-04-11)"
    )
  )
)
image_summary$QC_Note <- ifelse(
  image_summary$Measured_Mitochondria_n < 3,
  "LOW OBJECT COUNT (<3)",
  "OK"
)
image_summary <- image_summary[order(image_summary$Condition, image_summary$Image),]

metric_specs <- list(
  list(panel = "K", metric = "Mitochondrial area", column = "Mean_Area_um2", unit = "µm²"),
  list(panel = "L", metric = "Aspect ratio", column = "Mean_Aspect_Ratio", unit = "unitless"),
  list(panel = "M", metric = "Mitochondrial perimeter", column = "Mean_Perimeter_um", unit = "µm")
)

cliffs_delta <- function(x, y) {
  comparisons <- outer(x, y, FUN = "-")
  (sum(comparisons > 0) - sum(comparisons < 0)) / length(comparisons)
}

summary_rows <- list()
for (spec in metric_specs) {
  x <- image_summary[image_summary$Condition == "WR", spec$column]
  y <- image_summary[image_summary$Condition == "CS", spec$column]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  test <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    Panel = spec$panel,
    Metric = spec$metric,
    Unit = spec$unit,
    WR_Images_n = length(x),
    CS_Images_n = length(y),
    WR_Image_Mean = mean(x),
    CS_Image_Mean = mean(y),
    WR_Image_Median = median(x),
    CS_Image_Median = median(y),
    Mean_Difference_CS_minus_WR = mean(y) - mean(x),
    Fold_Change_CS_over_WR = mean(y) / mean(x),
    Cliffs_Delta_CS_vs_WR = cliffs_delta(y, x),
    Exploratory_Mann_Whitney_p = unname(test$p.value),
    stringsAsFactors = FALSE
  )
}
statistics <- do.call(rbind, summary_rows)

# Sensitivity analysis omitting the older Cnt TIFF set from WR.
sensitivity <- image_summary[!grepl("^Cnt ", image_summary$Image),]
sensitivity_rows <- list()
for (spec in metric_specs) {
  x <- sensitivity[sensitivity$Condition == "WR", spec$column]
  y <- sensitivity[sensitivity$Condition == "CS", spec$column]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  test <- wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)
  sensitivity_rows[[length(sensitivity_rows) + 1]] <- data.frame(
    Panel = spec$panel,
    Metric = spec$metric,
    WR_Images_n = length(x),
    CS_Images_n = length(y),
    WR_Image_Mean = mean(x),
    CS_Image_Mean = mean(y),
    Mean_Difference_CS_minus_WR = mean(y) - mean(x),
    Fold_Change_CS_over_WR = mean(y) / mean(x),
    Exploratory_Mann_Whitney_p = unname(test$p.value),
    stringsAsFactors = FALSE
  )
}
sensitivity_statistics <- do.call(rbind, sensitivity_rows)

# Coverage audit of all images supplied in the requested folder.
inventory <- read.csv(inventory_file, check.names = FALSE, stringsAsFactors = FALSE)
coverage <- merge(
  inventory,
  image_counts,
  by.x = c("Condition", "Source file"),
  by.y = c("Condition", "Image"),
  all.x = TRUE
)
coverage$Measured_Mitochondria_n[is.na(coverage$Measured_Mitochondria_n)] <- 0
coverage$Morphology_Coverage <- ifelse(
  coverage$Measured_Mitochondria_n > 0,
  "MEASURED",
  "NO MEASUREMENTS IN SOURCE WORKBOOK"
)
coverage <- coverage[order(coverage$Condition, coverage$`Source file`),]

# Export auditable flat tables.
write.csv(raw, file.path(output_dir, "TEM_QC_all_source_rows.csv"), row.names = FALSE)
write.csv(analysis_dat, file.path(output_dir, "TEM_clean_mitochondria.csv"), row.names = FALSE)
write.csv(image_summary, file.path(output_dir, "TEM_image_level_summary.csv"), row.names = FALSE)
write.csv(statistics, file.path(output_dir, "TEM_primary_statistics.csv"), row.names = FALSE)
write.csv(
  sensitivity_statistics,
  file.path(output_dir, "TEM_sensitivity_statistics_WRF_only.csv"),
  row.names = FALSE
)
write.csv(coverage, file.path(output_dir, "TEM_image_coverage_QC.csv"), row.names = FALSE)

# Publication plots: each dot is one TEM field; boxes show median and IQR,
# with Tukey whiskers extending to the most extreme value within 1.5 × IQR.
# Appearance is inherited from the project-wide figure_style.R settings.
geno_order <- c("WR", "CS")
group_labels <- c("WR" = "dMIC60WT \u2640", "CS" = "dMIC60CS \u2640")

make_panel <- function(data, value_column, y_label, plot_title, p_value) {
  plot_dat <- data.frame(
    Condition = factor(data$Condition, levels = geno_order),
    Value = data[[value_column]]
  )
  plot_dat <- plot_dat[is.finite(plot_dat$Value),]
  ggplot(plot_dat, aes(x = Condition, y = Value, fill = Condition)) +
    fig_summary(style = "box") +
    geom_point(
      aes(color = Condition),
      position = position_jitter(
        width = FIG_JITTER_WIDTH,
        height = 0,
        seed = FIG_SEED
      ),
      shape = FIG_PT_SHAPE,
      size = FIG_PT_SIZE,
      alpha = FIG_PT_ALPHA,
      show.legend = FALSE
    ) +
    fig_scale_fill(geno_order) +
    fig_scale_color(geno_order) +
    fig_scale_x_group(group_labels) +
    fig_scale_y(y_label) +
    fig_sig_bracket(p_value, plot_dat$Value) +
    labs(x = NULL, title = plot_title) +
    theme_fig() +
    theme(
      plot.margin = margin(
        t = FIG_MARGIN_PT[["t"]],
        r = FIG_MARGIN_PT[["r"]],
        b = FIG_MARGIN_PT[["b"]],
        l = FIG_MARGIN_PT[["l"]]
      )
    )
}

area_p <- statistics$Exploratory_Mann_Whitney_p[statistics$Metric == "Mitochondrial area"]
aspect_ratio_p <- statistics$Exploratory_Mann_Whitney_p[statistics$Metric == "Aspect ratio"]
perimeter_p <- statistics$Exploratory_Mann_Whitney_p[statistics$Metric == "Mitochondrial perimeter"]

area_plot <- make_panel(
  image_summary,
  "Mean_Area_um2",
  "Mean mitochondrial area (µm²)",
  "Mitochondrial area",
  area_p
)
aspect_ratio_plot <- make_panel(
  image_summary,
  "Mean_Aspect_Ratio",
  "Mean mitochondrial aspect ratio",
  "Mitochondrial aspect ratio",
  aspect_ratio_p
)
perimeter_plot <- make_panel(
  image_summary,
  "Mean_Perimeter_um",
  "Mean mitochondrial perimeter (µm)",
  "Mitochondrial perimeter",
  perimeter_p
)

fig_save(area_plot, file.path(output_dir, "TEM_mitochondrial_area_audited"))
fig_save(
  aspect_ratio_plot,
  file.path(output_dir, "TEM_mitochondrial_aspect_ratio_audited")
)
fig_save(
  perimeter_plot,
  file.path(output_dir, "TEM_mitochondrial_perimeter_audited")
)

cat("Included mitochondria:", nrow(analysis_dat), "\n")
cat("WR images:", sum(image_summary$Condition == "WR"), "\n")
cat("CS images:", sum(image_summary$Condition == "CS"), "\n")
print(statistics)
