#!/usr/bin/env Rscript
# A compact, lightly gridded PCA with nearby sample labels.
# Plots archived, verified scores without recomputing or moving observations.
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript.")
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg), fixed = TRUE))
root <- dirname(dirname(script_path))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else file.path(root, "Rebuilt_Output")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
co <- read.csv(file.path(root, "Supporting_Data", "Fig1B_ordination_coordinates.csv"))
checks <- read.csv(file.path(root, "Supporting_Data", "fresh_ordination_check.csv"))
colors <- c(WT_F = "#5AB4E5", WT_M = "#455A9F", CS_F = "#CB78A8", CS_M = "#881840")
co$group <- factor(co$group, levels = names(colors))
stopifnot(nrow(co) == 12L, !anyDuplicated(co$sample), !anyNA(co$group),
          all(table(co$group) == 3L), all(is.finite(co$PC1)), all(is.finite(co$PC2)))

p <- ggplot(co, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.1) +
  geom_text_repel(aes(label = label), size = 2.4, family = "Arial",
    seed = 12, box.padding = 0.22, point.padding = 0.18,
    force = 0.6, force_pull = 2, max.overlaps = Inf, max.iter = 20000,
    min.segment.length = 0.8, segment.size = 0.25, show.legend = FALSE) +
  scale_color_manual(values = colors,
    labels = c(WT_F = "WT-F", WT_M = "WT-M", CS_F = "CS-F", CS_M = "CS-M")) +
  scale_x_continuous(breaks = seq(-100, 100, 50), expand = expansion(mult = 0.13)) +
  scale_y_continuous(breaks = seq(-50, 50, 25), expand = expansion(mult = 0.14)) +
  labs(x = sprintf("PC1 (%.1f%%)", checks$PC1_variance_percent[1]),
       y = sprintf("PC2 (%.1f%%)", checks$PC2_variance_percent[1])) +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(color = "#EBEBEB", linewidth = 0.22),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.margin = margin(4, 5, 4, 4),
        legend.text = element_text(size = 7.5), legend.position = "top",
        legend.key.width = grid::unit(10, "pt"),
        legend.key.height = grid::unit(10, "pt"),
        legend.spacing.x = grid::unit(4, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = grid::unit(3, "pt")) +
  guides(color = guide_legend(nrow = 1))

stem <- file.path(out, "Fig1B_PCA_Compact_Labels")
if (capabilities("aqua")) {
  device <- function(filename, ...) grDevices::quartz(type = "pdf", file = filename, ...)
} else if (capabilities("cairo")) {
  device <- grDevices::cairo_pdf
} else stop("Quartz (macOS) or Cairo is required.")
ggsave(paste0(stem, ".pdf"), p, width = 260/72, height = 210/72,
       device = device, bg = "white", limitsize = FALSE)
ggsave(paste0(stem, ".png"), p, width = 260/72, height = 210/72,
       dpi = 300, bg = "white", limitsize = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "R_sessionInfo_compact.txt"))
message("Saved compact PCA PDF and PNG to ", normalizePath(out))
