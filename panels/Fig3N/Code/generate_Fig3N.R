# Descriptive TMRM fields from three flies/genotype; fly-to-field mapping is unavailable.
# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", .rebuild_args[1]), fixed=TRUE))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
required <- c("ggplot2", "scales", "ggtext")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
})
args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(args)) dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", args[[1]]), fixed=TRUE))) else normalizePath(getwd())
package_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
FIG_SHOW_NS <- FALSE
# Fit the existing panel without the historical reserved blank footer.
ggsave <- ggplot2::ggsave
if (identical(Sys.info()[["sysname"]], "Darwin") && capabilities("aqua")) {
    FIG_FONT <- "Arial"
    quartz_pdf_unicode <- function(filename, width, height, ...) {
        grDevices::quartz(type = "pdf", file = filename, width = width, height = height, family = FIG_FONT, ...)
    }
    FIG_PDF_DEVICE <- quartz_pdf_unicode
}
support_dir <- file.path(package_dir, "Source_Data")
output_dir <- file.path(package_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
objects <- read.csv(file.path(support_dir, "TMRM_MTG_connected_object_measurements.csv"), stringsAsFactors = FALSE, check.names = FALSE)
images <- read.csv(file.path(support_dir, "TMRM_MTG_connected_object_image_summary.csv"), stringsAsFactors = FALSE, check.names = FALSE)
required_object_columns <- c("file", "condition", "object_ratio_percent_of_WT_image_mean")
required_image_columns <- c("file", "condition", "objects_n", "image_mean_percent_of_WT")
if (length(setdiff(required_object_columns, names(objects)))) stop("Object table is missing required columns")
if (length(setdiff(required_image_columns, names(images)))) stop("Image table is missing required columns")
if (any(!is.finite(objects$object_ratio_percent_of_WT_image_mean))) stop("Non-finite object ratios detected")
if (any(!is.finite(images$image_mean_percent_of_WT))) stop("Non-finite image means detected")
if (anyDuplicated(paste(objects$file, objects$object_id_within_image, sep = "|"))) stop("Duplicate object IDs detected")
if (anyDuplicated(images$file)) stop("Duplicate image summaries detected")
objects$Group <- factor(objects$condition, levels = c("WR", "CS"))
images$Group <- factor(images$condition, levels = c("WR", "CS"))
# Shared reference is the pooled mean of WT field means, not a per-experiment control.
pooled_WT_mean <- mean(images$image_mean_object_ratio[images$condition=="WR"])
recomputed_pct <- 100 * images$image_mean_object_ratio / pooled_WT_mean
stopifnot(max(abs(recomputed_pct-images$image_mean_percent_of_WT))<1e-10)
images$image_mean_percent_of_WT <- recomputed_pct

group_order <- c("WR", "CS")
# Use real Arial italic text runs; Unicode-only/math fonts can fall back to serif.
group_labels <- c(WR = "<i>dMIC60</i>-WT", CS = "<i>dMIC60</i>-CS")
wr <- images$image_mean_percent_of_WT[images$Group == "WR"]
cs <- images$image_mean_percent_of_WT[images$Group == "CS"]
has_ties <- anyDuplicated(images$image_mean_percent_of_WT) > 0
test <- wilcox.test(cs, wr, alternative = "two.sided", exact = !has_ties, correct = has_ties)
cliffs_delta <- function(wr_values, cs_values) {
    comparisons <- outer(cs_values, wr_values, FUN = "-")
    (sum(comparisons > 0) - sum(comparisons < 0))/length(comparisons)
}
statistics <- data.frame(Analysis = "Connected MTG-positive object TMRM/MTG ratio", Mask = "3x3 median-filtered MTG; per-image Otsu; 8-connected; minimum 0.05 um2",
    Analysis_Unit = "image mean of connected-object ratios", WT_Images_n = length(wr), CS_Images_n = length(cs), WT_Objects_n = sum(objects$Group ==
        "WR"), CS_Objects_n = sum(objects$Group == "CS"), WT_Image_Mean_Percent = mean(wr), CS_Image_Mean_Percent = mean(cs),
    WT_Image_Median_Percent = median(wr), CS_Image_Median_Percent = median(cs), CS_over_WT_Image_Mean = mean(cs)/mean(wr),
    Cliffs_Delta_CS_vs_WT = cliffs_delta(wr, cs), Exploratory_Mann_Whitney_p = test$p.value, stringsAsFactors = FALSE)
write.csv(statistics, .rebuild_file(file.path(support_dir, "TMRM_MTG_connected_object_primary_statistics.csv")), row.names = FALSE)
set.seed(FIG_SEED)
p <- ggplot(images, aes(x = Group, y = image_mean_percent_of_WT, fill = Group)) + fig_reference_line(100) + fig_summary() +
    fig_points(width = 0.09, color = "black") + fig_scale_fill(group_order) +
    fig_scale_x_group(group_labels) + scale_y_continuous(breaks = seq(0, 400, 100), labels = scales::label_number(accuracy = 1),
    name = "Field mean TMRM / MTG (% of WT)", expand = expansion(mult = c(0.02, 0.08))) + labs(x = NULL) + theme_fig() +
    theme(axis.text.x = ggtext::element_markdown(family = "Arial"))
# Presentation adapted to the existing Figure 3N footprint; analysis unchanged.
p <- p + theme(axis.text = element_text(size=6, color='black'),
  axis.text.x = ggtext::element_markdown(family='Arial', size=6, lineheight=0.9),
  axis.title.y = element_text(size=7), plot.margin = margin(3, 3, 3, 3)) +
  scale_x_discrete(labels=c(WR='<i>dMIC60</i>-<br>WT', CS='<i>dMIC60</i>-<br>CS'),
    expand=expansion(add=FIG_X_EXPAND))
output_base <- file.path(output_dir, "Fig3N_TMRM_regenerated")
fig_save(p, output_base, width = 125/72, height = 140/72, formats = c("png", "pdf"))
cat("Connected-object Figure 3N complete.\n")
print(statistics)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))

write.csv(images, file.path(output_dir, "Fig3N_plotted_values.csv"), row.names=FALSE)
write.csv(ggplot_build(p)$data[[4]], file.path(output_dir, "Fig3N_point_geometry.csv"), row.names=FALSE)
