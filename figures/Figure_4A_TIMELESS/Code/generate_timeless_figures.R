# ============================================================
# HeLa TIMELESS immunofluorescence analysis
# Source: 7-27-26 HeLa Timeless Quantification.xlsx
#
# Statistical design:
#   - MIC60-KO HeLa cells complemented with dMIC60-WT or dMIC60-CS
#   - Analyze the 20X untreated and 63X 20 mM H2O2 acquisitions separately
#   - Average cell-level CTCF values within each image for visualization
#   - One coverglass per condition: biological n = 1 per group
#   - Images and cells are technical subsamples; no inferential test is valid
#   - No outlier removal
# ============================================================

library(readxl)
library(ggplot2)
library(patchwork)

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_arg) > 0) {
  normalizePath(gsub("~\\+~", " ", sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath("HeLa Timeless Data Analysis.R")
}

script_dir <- dirname(script_path)
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
publication_dir <- normalizePath(file.path(figure_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))

input_file <- file.path(
  figure_dir,
  "Original_Data",
  "7-27-26 HeLa Timeless Quantification.xlsx"
)

output_dir <- file.path(
  figure_dir,
  "Rebuilt_Output"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Helpers
# ============================================================

trim_character <- function(x) {
  ifelse(is.na(x), NA_character_, trimws(as.character(x)))
}

first_non_missing <- function(x) {
  x <- trim_character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else x[[1]]
}

image_summary <- function(cell_data, acquisition, treatment) {
  means <- aggregate(
    cell_data$ctcf,
    by = list(group = cell_data$group, image = cell_data$image),
    FUN = mean,
    na.rm = TRUE
  )
  names(means)[3] <- "image_mean_ctcf"

  counts <- aggregate(
    cell_data$ctcf,
    by = list(group = cell_data$group, image = cell_data$image),
    FUN = length
  )
  names(counts)[3] <- "cell_count"

  means <- merge(means, counts, by = c("group", "image"), sort = FALSE)
  means$acquisition <- acquisition
  means$treatment <- treatment
  means <- means[, c(
    "acquisition", "treatment", "group", "image",
    "cell_count", "image_mean_ctcf"
  )]

  wt_reference <- mean(
    means$image_mean_ctcf[means$group == "WT"],
    na.rm = TRUE
  )
  means$wt_reference_mean_ctcf <- wt_reference
  means$percent_wt <- means$image_mean_ctcf / wt_reference * 100
  means$group <- factor(means$group, levels = c("WT", "CS"))
  means[order(means$group, means$image), ]
}

analyze_acquisition <- function(image_data) {
  group_summary <- do.call(
    rbind,
    lapply(c("WT", "CS"), function(g) {
      x <- image_data$image_mean_ctcf[image_data$group == g]
      y <- image_data$percent_wt[image_data$group == g]
      data.frame(
        acquisition = unique(image_data$acquisition),
        treatment = unique(image_data$treatment),
        group = g,
        coverglasses = 1,
        images = length(x),
        cells = sum(image_data$cell_count[image_data$group == g]),
        mean_ctcf = mean(x),
        median_ctcf = median(x),
        sd_ctcf = sd(x),
        sem_ctcf = sd(x) / sqrt(length(x)),
        mean_percent_wt = mean(y),
        median_percent_wt = median(y),
        stringsAsFactors = FALSE
      )
    })
  )

  statistics <- data.frame(
    acquisition = unique(image_data$acquisition),
    treatment = unique(image_data$treatment),
    comparison = "dMIC60-CS versus dMIC60-WT",
    biological_unit = "Coverglass",
    wt_coverglasses = 1,
    cs_coverglasses = 1,
    inferential_test = "Not performed: n = 1 coverglass per group",
    wt_images = sum(image_data$group == "WT"),
    cs_images = sum(image_data$group == "CS"),
    wt_cells = sum(image_data$cell_count[image_data$group == "WT"]),
    cs_cells = sum(image_data$cell_count[image_data$group == "CS"]),
    wt_mean_ctcf = mean(image_data$image_mean_ctcf[image_data$group == "WT"]),
    cs_mean_ctcf = mean(image_data$image_mean_ctcf[image_data$group == "CS"]),
    cs_percent_of_wt_mean = mean(image_data$percent_wt[image_data$group == "CS"]),
    direction = if (
      mean(image_data$image_mean_ctcf[image_data$group == "CS"]) >
        mean(image_data$image_mean_ctcf[image_data$group == "WT"])
    ) "CS higher" else "CS lower",
    interpretation = paste0(
      "Descriptive only; images and cells are technical subsamples from ",
      "one coverglass per group"
    ),
    outlier_removal = "None",
    stringsAsFactors = FALSE
  )

  list(summary = group_summary, statistics = statistics)
}


# ============================================================
# Parse 20X untreated sheet
# ============================================================

nt <- read_excel(
  input_file,
  sheet = "20X NT",
  skip = 1,
  .name_repair = "unique"
)
names(nt)[1:8] <- c(
  "cell_index", "image", "area", "mean_fluorescence",
  "integrated_density", "raw_integrated_density", "ctcf", "condition"
)
nt$image <- trim_character(nt$image)

nt_condition_map <- tapply(nt$condition, nt$image, first_non_missing)
nt$group <- unname(nt_condition_map[nt$image])

nt_cells <- nt[
  !is.na(nt$ctcf) &
    is.finite(nt$ctcf) &
    nt$group %in% c("WT", "CS"),
  c(
    "image", "group", "cell_index", "area", "mean_fluorescence",
    "integrated_density", "raw_integrated_density", "ctcf"
  )
]
nt_cells$acquisition <- "20X NT"
nt_cells$treatment <- "Untreated"
nt_cells <- nt_cells[, c(
  "acquisition", "treatment", "group", "image", "cell_index",
  "area", "mean_fluorescence", "integrated_density",
  "raw_integrated_density", "ctcf"
)]


# ============================================================
# Parse 63X 20 mM H2O2 sheet
# ============================================================

top <- read_excel(
  input_file,
  sheet = "63X 20mM",
  range = "A2:I23",
  .name_repair = "unique"
)
names(top)[1:9] <- c(
  "image_condition", "cell", "area", "mean_fluorescence", "min",
  "max", "integrated_density", "raw_integrated_density", "ctcf"
)

top_rows <- list()
current_image <- NA_character_
current_group <- NA_character_
for (i in seq_len(nrow(top))) {
  label <- trim_character(top$image_condition[[i]])
  cell_label <- trim_character(top$cell[[i]])
  if (!is.na(label) && identical(cell_label, "blank")) {
    current_image <- label
  }
  if (!is.na(label) && grepl("^(WT|CS)", label)) {
    current_group <- sub(" .*$", "", label)
  }
  if (!is.na(top$ctcf[[i]]) && is.finite(top$ctcf[[i]])) {
    top_rows[[length(top_rows) + 1]] <- data.frame(
      acquisition = "63X 20mM",
      treatment = "20 mM H2O2",
      group = current_group,
      image = current_image,
      cell_index = cell_label,
      area = top$area[[i]],
      mean_fluorescence = top$mean_fluorescence[[i]],
      integrated_density = top$integrated_density[[i]],
      raw_integrated_density = top$raw_integrated_density[[i]],
      ctcf = top$ctcf[[i]],
      stringsAsFactors = FALSE
    )
  }
}
top_cells <- do.call(rbind, top_rows)

lower <- read_excel(
  input_file,
  sheet = "63X 20mM",
  range = "A26:I49",
  .name_repair = "unique"
)
names(lower)[1:9] <- c(
  "condition", "cell_index", "image", "area", "mean_fluorescence",
  "integrated_density", "raw_integrated_density", "ctcf", "image_average"
)
lower$image <- trim_character(lower$image)
lower_condition_map <- tapply(lower$condition, lower$image, first_non_missing)
lower$group <- sub(" .*$", "", unname(lower_condition_map[lower$image]))

lower_cells <- lower[
  !is.na(lower$ctcf) &
    is.finite(lower$ctcf) &
    lower$group %in% c("WT", "CS"),
  c(
    "group", "image", "cell_index", "area", "mean_fluorescence",
    "integrated_density", "raw_integrated_density", "ctcf"
  )
]
lower_cells$acquisition <- "63X 20mM"
lower_cells$treatment <- "20 mM H2O2"
lower_cells <- lower_cells[, names(top_cells)]

h2o2_cells <- rbind(top_cells, lower_cells)


# ============================================================
# Image-level analysis and verification
# ============================================================

nt_images <- image_summary(nt_cells, "20X NT", "Untreated")
h2o2_images <- image_summary(h2o2_cells, "63X 20mM", "20 mM H2O2")
image_data <- rbind(nt_images, h2o2_images)
image_data$group <- factor(image_data$group, levels = c("WT", "CS"))
image_data$coverglass <- paste(
  image_data$acquisition,
  image_data$group,
  "coverglass 1",
  sep = " | "
)

nt_result <- analyze_acquisition(nt_images)
h2o2_result <- analyze_acquisition(h2o2_images)
group_summary <- rbind(nt_result$summary, h2o2_result$summary)
statistics <- rbind(nt_result$statistics, h2o2_result$statistics)

nt_cells$percent_wt <- (
  nt_cells$ctcf / unique(nt_images$wt_reference_mean_ctcf) * 100
)
h2o2_cells$percent_wt <- (
  h2o2_cells$ctcf / unique(h2o2_images$wt_reference_mean_ctcf) * 100
)
nt_cells$group <- factor(nt_cells$group, levels = c("WT", "CS"))
h2o2_cells$group <- factor(h2o2_cells$group, levels = c("WT", "CS"))
nt_cells$coverglass <- paste("20X NT", nt_cells$group, "coverglass 1", sep = " | ")
h2o2_cells$coverglass <- paste("63X 20mM", h2o2_cells$group, "coverglass 1", sep = " | ")

# Check that parsing recovered the workbook's intended sample hierarchy.
stopifnot(
  sum(nt_images$group == "WT") == 5,
  sum(nt_images$group == "CS") == 5,
  sum(nt_images$cell_count[nt_images$group == "WT"]) == 41,
  sum(nt_images$cell_count[nt_images$group == "CS"]) == 42,
  sum(h2o2_images$group == "WT") == 7,
  sum(h2o2_images$group == "CS") == 6,
  sum(h2o2_images$cell_count[h2o2_images$group == "WT"]) == 13,
  sum(h2o2_images$cell_count[h2o2_images$group == "CS"]) == 15
)


# ============================================================
# Publication-style figures
# ============================================================

make_panel <- function(data, result, subtitle_text) {
  n_wt <- sum(data$group == "WT")
  n_cs <- sum(data$group == "CS")

  set.seed(FIG_SEED)
  ggplot(data, aes(x = group, y = percent_wt, fill = group, color = group)) +
    fig_reference_line(100) +
    fig_summary(p_value = NA_real_) +
    fig_points() +
    fig_scale_fill(c("WT", "CS")) +
    fig_scale_color(c("WT", "CS"), guide = "none") +
    fig_scale_x_group(c("WT" = "dMIC60 WT", "CS" = "dMIC60 CS")) +
    fig_scale_y("TIMELESS CTCF (% of WT image mean)") +
    labs(
      x = NULL,
      title = "TIMELESS abundance in MIC60-KO HeLa cells",
      subtitle = paste0(
        subtitle_text,
        "; each dot is the mean CTCF from one image"
      ),
      caption = paste0(
        "Technical subsamples (not independent replicates):\n",
        "WT: ", n_wt, " images (",
        sum(data$cell_count[data$group == "WT"]),
        " cells); CS: ", n_cs, " images (",
        sum(data$cell_count[data$group == "CS"]),
        " cells). No outliers removed.\n",
        "Boxes: median and IQR of image means.\n",
        "Descriptive only; no condition-level test with n = 1 coverglass/group."
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_fig()
}

p_nt <- make_panel(nt_images, nt_result, "Untreated; 20x acquisition")
p_h2o2 <- make_panel(h2o2_images, h2o2_result, "20 mM H2O2; 63x acquisition")

fig_save(
  p_nt,
  file.path(output_dir, "HeLa_TIMELESS_untreated_WT_vs_CS"),
  width = 4.8,
  height = 5.4
)

fig_save(
  p_h2o2,
  file.path(output_dir, "HeLa_TIMELESS_20mM_H2O2_WT_vs_CS"),
  width = 4.8,
  height = 5.4
)

p_combined <- (
  (p_nt + labs(title = "Untreated") +
     theme(plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0.5))) |
  (p_h2o2 + labs(title = "20 mM H2O2") +
     theme(plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0.5)))
)
fig_save(
  p_combined,
  file.path(output_dir, "HeLa_TIMELESS_WT_vs_CS_all_acquisitions"),
  width = 9.6,
  height = 5.4
)


# ============================================================
# Versions showing every cell-level measurement
# ============================================================

make_all_points_panel <- function(image_data, cell_data, result, subtitle_text) {
  n_wt <- sum(image_data$group == "WT")
  n_cs <- sum(image_data$group == "CS")

  set.seed(FIG_SEED)
  ggplot(image_data, aes(x = group, y = percent_wt, fill = group)) +
    fig_reference_line(100) +
    fig_summary(p_value = NA_real_) +
    geom_point(
      data = cell_data,
      mapping = aes(x = group, y = percent_wt, color = group),
      inherit.aes = FALSE,
      position = position_jitter(
        width = FIG_JITTER_WIDTH,
        height = 0,
        seed = FIG_SEED
      ),
      shape = FIG_PT_SHAPE,
      size = FIG_PT_SIZE_DENSE,
      alpha = FIG_PT_ALPHA_DENSE,
      show.legend = FALSE
    ) +
    fig_markers(
      image_data,
      aes(x = group, y = percent_wt)
    ) +
    fig_scale_fill(c("WT", "CS")) +
    fig_scale_color(c("WT", "CS"), guide = "none") +
    fig_scale_x_group(c("WT" = "dMIC60 WT", "CS" = "dMIC60 CS")) +
    fig_scale_y("TIMELESS CTCF (% of WT image mean)") +
    labs(
      x = NULL,
      title = "TIMELESS abundance in MIC60-KO HeLa cells",
      subtitle = paste0(
        subtitle_text,
        "; all quantified cells are shown"
      ),
      caption = paste0(
        "Small dots: individual cells; white diamonds: image means.\n",
        "Boxes: median and IQR of image means; no outliers removed.\n",
        "WT: ", n_wt, " images (",
        sum(cell_data$group == "WT"), " cells); CS: ",
        n_cs, " images (", sum(cell_data$group == "CS"),
        " cells).\n",
        "Descriptive only; no condition-level test with n = 1 coverglass/group."
      )
    ) +
    coord_cartesian(clip = "off") +
    theme_fig()
}

p_nt_all <- make_all_points_panel(
  nt_images,
  nt_cells,
  nt_result,
  "Untreated; 20x acquisition"
)
p_h2o2_all <- make_all_points_panel(
  h2o2_images,
  h2o2_cells,
  h2o2_result,
  "20 mM H2O2; 63x acquisition"
)

fig_save(
  p_nt_all,
  file.path(output_dir, "HeLa_TIMELESS_untreated_WT_vs_CS_all_points"),
  width = 4.8,
  height = 5.4
)
fig_save(
  p_h2o2_all,
  file.path(output_dir, "HeLa_TIMELESS_20mM_H2O2_WT_vs_CS_all_points"),
  width = 4.8,
  height = 5.4
)

p_combined_all <- (
  (p_nt_all + labs(title = "Untreated") +
     theme(plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0.5))) |
  (p_h2o2_all + labs(title = "20 mM H2O2") +
     theme(plot.title = element_text(size = FIG_AXIS_TITLE_SIZE, hjust = 0.5)))
)

fig_save(
  p_combined_all,
  file.path(output_dir, "HeLa_TIMELESS_WT_vs_CS_all_acquisitions_all_points"),
  width = 9.6,
  height = 5.4
)


# ============================================================
# Save auditable data and results
# ============================================================

cell_data <- rbind(nt_cells, h2o2_cells)
cell_data$group <- factor(cell_data$group, levels = c("WT", "CS"))
cell_data <- cell_data[order(
  cell_data$acquisition,
  cell_data$group,
  cell_data$image,
  cell_data$cell_index
), ]

write.csv(
  cell_data,
  file.path(output_dir, "HeLa_TIMELESS_cell_level_CTCF.csv"),
  row.names = FALSE
)
write.csv(
  image_data,
  file.path(output_dir, "HeLa_TIMELESS_image_level_means.csv"),
  row.names = FALSE
)
write.csv(
  group_summary,
  file.path(output_dir, "HeLa_TIMELESS_group_summary.csv"),
  row.names = FALSE
)
write.csv(
  statistics,
  file.path(output_dir, "HeLa_TIMELESS_statistical_results.csv"),
  row.names = FALSE
)

cat("\nImage-level values:\n")
print(image_data)
cat("\nStatistical results:\n")
print(statistics)
