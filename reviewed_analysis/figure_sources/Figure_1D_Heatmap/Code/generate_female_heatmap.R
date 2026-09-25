#!/usr/bin/env Rscript
# Reproduce the corrected 20 September 2026 Fig. 1D panel from its
# verified row-z-score matrix. This script does not rerun DESeq2.
# Run: Rscript Code/generate_female_heatmap.R [output_directory]
# Default output: Rebuilt_Output/Fig1D_Female_Heatmap.pdf
suppressPackageStartupMessages({library(pheatmap); library(grid)})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this script with Rscript.")
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE))
root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else file.path(root, "Rebuilt_Output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

z <- as.matrix(read.csv(file.path(root, "Supporting_Data", "heatmap_female_row_z_scores.csv"),
                       row.names = 1, check.names = FALSE))
selected <- read.csv(file.path(root, "Supporting_Data", "heatmap_female_selected_genes.csv"))
stopifnot(nrow(z) == 60L, ncol(z) == 12L, all(is.finite(z)),
          !anyDuplicated(rownames(z)), setequal(rownames(z), selected$gene),
          sum(selected$direction == "up") == 30L, sum(selected$direction == "down") == 30L)

group_colors <- c(WT_F = "#5AB4E5", WT_M = "#455A9F", CS_F = "#CB78A8", CS_M = "#881840")
samples <- colnames(z)
group <- factor(paste0(ifelse(grepl("^CS", samples), "CS", "WT"), "_",
                       ifelse(grepl("F[123]$", samples), "F", "M")), levels = names(group_colors))
stopifnot(!anyNA(group), all(table(group) == 3L))
annotation <- data.frame(Group = group, row.names = samples)
labels <- sub("^WR", "WT", samples)
labels <- sub("^(WT|CS)(F|M)([123])$", "\\1-\\2\\3", labels)

# Match the archived panel: 300 x 550 points, Arial 7.5 pt,
# Euclidean distance, complete linkage, and a shared -3 to +3 color scale.
grDevices::pdfFonts(Arial = grDevices::pdfFonts("Helvetica")[[1]])
output_pdf <- file.path(out, "Fig1D_Female_Heatmap.pdf")
if (capabilities("aqua")) {
  grDevices::quartz(type = "pdf", file = output_pdf, width = 300/72, height = 550/72, family = "Arial")
} else if (capabilities("cairo")) {
  grDevices::cairo_pdf(output_pdf, width = 300/72, height = 550/72, family = "Arial")
} else {
  stop("Quartz (macOS) or Cairo is required for PDF export.")
}
hm <- pheatmap(z,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-3, 3, length.out = 102),
  cluster_rows = TRUE, cluster_cols = TRUE,
  clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
  clustering_method = "complete", annotation_col = annotation,
  annotation_colors = list(Group = group_colors),
  annotation_names_col = FALSE, annotation_legend = FALSE,
  labels_row = as.expression(lapply(rownames(z), function(g) bquote(italic(.(g))))),
  labels_col = labels, fontsize = 7.5, fontsize_row = 7.5, fontsize_col = 7.5,
  fontfamily = "Arial", border_color = NA, treeheight_row = 17, treeheight_col = 17,
  angle_col = "90", legend_breaks = c(-3, 0, 3), silent = TRUE)
grid.newpage()
grid.draw(hm$gtable)
dev.off()
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo.txt"))
message("Saved ", normalizePath(output_pdf))
