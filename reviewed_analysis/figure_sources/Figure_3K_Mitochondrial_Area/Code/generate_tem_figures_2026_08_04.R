# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
required <- c("ggplot2", "readxl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
    library(ggplot2)
    library(readxl)
})
args <- grep("^--file=", gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE), value = TRUE)
script_dir <- if (length(args)) dirname(normalizePath(sub("^--file=", "", args[[1]]))) else normalizePath(getwd())
package_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
FIG_SHOW_NS <- FALSE
input_file <- file.path(package_dir, "Original_Data", "8.4.26 TEM Retrace.xlsx")
output_dir <- file.path(package_dir, "Rebuilt_Output")
support_dir <- file.path(package_dir, "Supporting_Data")
qc_dir <- file.path(package_dir, "QC")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
raw <- as.data.frame(read_excel(input_file, sheet = "Sheet1", .name_repair = "minimal"))
names(raw) <- make.unique(trimws(names(raw)))
if (names(raw)[1] == "") names(raw)[1] <- "ImageJ_Row"
raw$Source_Row <- seq_len(nrow(raw)) + 1L
required_columns <- c("Label", "Area", "Perim.", "Major", "Minor", "AR")
if (length(setdiff(required_columns, names(raw)))) {
    stop("Missing required columns: ", paste(setdiff(required_columns, names(raw)), collapse = ", "))
}
raw$Label <- trimws(as.character(raw$Label))
raw$Condition <- ifelse(grepl("^(WR|Cnt)", raw$Label, ignore.case = TRUE), "WR", ifelse(grepl("^CS", raw$Label, ignore.case = TRUE), 
    "CS", NA_character_))
raw$Image <- raw$Label
numeric_columns <- intersect(c("Area", "Perim.", "Major", "Minor", "Feret", "MinFeret", "Width", "Height", "AR", "Round", 
    "Solidity"), names(raw))
for (column in numeric_columns) raw[[column]] <- suppressWarnings(as.numeric(raw[[column]]))
scale_bar_pixels <- c(`Cnt fmale_a_5000X_0005.tif` = 360, `Cnt fmale_a_5000X_0008.tif` = 358, `Cnt fmale_a_5000X_0010.tif` = 359, 
    `Cnt fmale_a_5000X_0012.tif` = 360, `WR-F -5000X -0009.tif` = 868)
raw$Scale_Bar_Pixels <- unname(scale_bar_pixels[raw$Label])
raw$Pixel_Size_um <- ifelse(is.na(raw$Scale_Bar_Pixels), NA_real_, 2/raw$Scale_Bar_Pixels)
raw$Calibration_Status <- ifelse(is.na(raw$Pixel_Size_um), "DM3 spatial calibration retained in workbook", "Converted from pixels using visible 2-um scale bar")
pixel_rows <- is.finite(raw$Pixel_Size_um)
raw$Area_um2 <- raw$Area
raw$Perimeter_um <- raw$Perim.
raw$Major_um <- raw$Major
raw$Minor_um <- raw$Minor
raw$Area_um2[pixel_rows] <- raw$Area[pixel_rows] * raw$Pixel_Size_um[pixel_rows]^2
for (pair in list(c("Perim.", "Perimeter_um"), c("Major", "Major_um"), c("Minor", "Minor_um"))) {
    source_col <- pair[[1]]
    target_col <- pair[[2]]
    raw[[target_col]][pixel_rows] <- raw[[source_col]][pixel_rows] * raw$Pixel_Size_um[pixel_rows]
}
raw$Aspect_Ratio <- raw$AR
raw$Area_Include <- is.finite(raw$Area_um2) & raw$Area_um2 > 0 & !is.na(raw$Condition)
raw$Perimeter_Include <- is.finite(raw$Perimeter_um) & raw$Perimeter_um > 0 & !is.na(raw$Condition)
raw$Aspect_Ratio_Include <- is.finite(raw$Aspect_Ratio) & raw$Aspect_Ratio > 0 & !is.na(raw$Condition)
raw$Aspect_Ratio_Exclusion <- ifelse(raw$Aspect_Ratio_Include, "", ifelse(is.na(raw$Condition), "Unrecognized condition", 
    "Ellipse major/minor axes absent; AR unavailable"))
duplicate_key <- paste(raw$Label, raw$Area, raw$Perim., raw$Major, raw$Minor, raw$AR, sep = "|")
raw$Exact_Duplicate <- duplicated(duplicate_key)
if (any(raw$Exact_Duplicate)) stop("Exact duplicate measurement rows detected in the August 4 workbook.")
write.csv(raw, .rebuild_file(file.path(qc_dir, "TEM_2026-08-04_all_source_rows_QC.csv")), row.names = FALSE, na = "")
calibration_audit <- unique(raw[, c("Image", "Condition", "Scale_Bar_Pixels", "Pixel_Size_um", "Calibration_Status")])
calibration_audit <- calibration_audit[order(calibration_audit$Condition, calibration_audit$Image), ]
write.csv(calibration_audit, .rebuild_file(file.path(qc_dir, "TEM_2026-08-04_calibration_audit.csv")), row.names = FALSE, 
    na = "")
metric_definitions <- data.frame(Panel = c("K", "L", "J"), Metric = c("Mitochondrial area", "Aspect ratio", "Mitochondrial perimeter"), 
    Value_Column = c("Area_um2", "Aspect_Ratio", "Perimeter_um"), Include_Column = c("Area_Include", "Aspect_Ratio_Include", 
        "Perimeter_Include"), Y_Label = c("Mitochondrial area (µm²)", "Mitochondrial aspect ratio", "Mitochondrial perimeter (µm)"), 
    File_Stem = c("TEM_mitochondrial_area_2026-08-04", "TEM_mitochondrial_aspect_ratio_2026-08-04", "TEM_mitochondrial_perimeter_2026-08-04"), 
    stringsAsFactors = FALSE)
clean_parts <- list()
image_parts <- list()
stats_parts <- list()
cliffs_delta <- function(wr, cs) {
    comparisons <- outer(cs, wr, FUN = "-")
    (sum(comparisons > 0) - sum(comparisons < 0))/length(comparisons)
}
for (i in seq_len(nrow(metric_definitions))) {
    spec <- metric_definitions[i, ]
    keep <- raw[[spec$Include_Column]]
    dat <- raw[keep, c("Source_Row", "Image", "Condition", spec$Value_Column, "Calibration_Status")]
    names(dat)[4] <- "Value"
    dat$Panel <- spec$Panel
    dat$Metric <- spec$Metric
    clean_parts[[i]] <- dat
    image_summary <- aggregate(Value ~ Image + Condition, data = dat, FUN = mean)
    names(image_summary)[3] <- "Image_Mean"
    image_summary$Panel <- spec$Panel
    image_summary$Metric <- spec$Metric
    image_summary$Objects_n <- as.integer(table(interaction(dat$Image, dat$Condition, drop = TRUE))[interaction(image_summary$Image, 
        image_summary$Condition, drop = TRUE)])
    image_summary$Condition <- factor(image_summary$Condition, levels = c("WR", "CS"))
    image_parts[[i]] <- image_summary
    wr <- image_summary$Image_Mean[image_summary$Condition == "WR"]
    cs <- image_summary$Image_Mean[image_summary$Condition == "CS"]
    has_ties <- anyDuplicated(image_summary$Image_Mean) > 0
    test <- wilcox.test(cs, wr, alternative = "two.sided", exact = !has_ties, correct = has_ties)
    stats_parts[[i]] <- data.frame(Panel = spec$Panel, Metric = spec$Metric, Analysis_Unit = "TEM image field mean", WR_Images_n = length(wr), 
        CS_Images_n = length(cs), WR_Objects_n = sum(dat$Condition == "WR"), CS_Objects_n = sum(dat$Condition == "CS"), WR_Image_Mean = mean(wr), 
        CS_Image_Mean = mean(cs), WR_Image_Median = median(wr), CS_Image_Median = median(cs), Mean_Difference_CS_minus_WR = mean(cs) - 
            mean(wr), Fold_Change_CS_over_WR = mean(cs)/mean(wr), Cliffs_Delta_CS_vs_WR = cliffs_delta(wr, cs), Exploratory_Mann_Whitney_p = test$p.value, 
        stringsAsFactors = FALSE)
    plot_data <- image_summary
    plot_data$Group <- factor(plot_data$Condition, levels = c("WR", "CS"))
    raw_plot_data <- dat
    raw_plot_data$Group <- factor(raw_plot_data$Condition, levels = c("WR", "CS"))
    group_order <- c("WR", "CS")
    group_labels <- FIG_GROUP_LABELS
    p <- ggplot(plot_data, aes(x = Group, y = Image_Mean, fill = Group, color = Group)) + fig_summary(p_value = test$p.value) + 
        fig_points(dense = TRUE, data = raw_plot_data, mapping = aes(x = Group, y = Value), inherit.aes = FALSE, width = 0.14, 
            color = "black") + fig_markers(plot_data, aes(x = Group, y = Image_Mean), width = 0.08) + fig_sig_bracket(test$p.value, raw_plot_data$Value, label=paste0('plain("Field-level")~italic(p)=="',format.pval(test$p.value,digits=2),'"'),parse_label=TRUE) + fig_scale_fill(group_order) + fig_scale_color(group_order, guide = "none") + fig_scale_x_group(group_labels) + 
        fig_scale_y(spec$Y_Label) + labs(x = NULL) + theme_fig()
    fig_save(p, file.path(output_dir, spec$File_Stem), width = 168/72, height = 179/72, formats = c("png", "pdf"))
}
clean_all <- do.call(rbind, clean_parts)
image_all <- do.call(rbind, image_parts)
stats_all <- do.call(rbind, stats_parts)
write.csv(clean_all, .rebuild_file(file.path(support_dir, "TEM_2026-08-04_clean_measurements_long.csv")), row.names = FALSE, 
    na = "")
write.csv(image_all, .rebuild_file(file.path(support_dir, "TEM_2026-08-04_image_level_summary.csv")), row.names = FALSE, 
    na = "")
write.csv(stats_all, .rebuild_file(file.path(support_dir, "TEM_2026-08-04_primary_statistics.csv")), row.names = FALSE, na = "")
cat("August 4 TEM retrace analysis complete.\n")
print(stats_all)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
