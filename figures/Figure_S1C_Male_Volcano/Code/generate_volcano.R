#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[[1]])))
panel_dir <- dirname(script_dir)
panel_name <- basename(panel_dir)
is_male <- grepl("Male", panel_name)
sex <- if (is_male) "Male" else "Female"
input_name <- if (is_male) "DEG_CSM_vs_WRM.csv" else "DEG_CSF_vs_WRF.csv"
output_stem <- if (is_male) "Figure_S1C_Male_CS_vs_WT_volcano" else "Figure_1C_Female_CS_vs_WT_volcano"

d <- read.csv(file.path(panel_dir, "Original_Data", input_name), check.names = FALSE)
d <- d[, nzchar(names(d)), drop = FALSE]
if (!"gene" %in% names(d)) d$gene <- rownames(d)
d <- d |>
  filter(!is.na(padj)) |>
  mutate(
    neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
    direction = case_when(
      padj < 0.05 & log2FoldChange >= 0.58 ~ "Up in CS",
      padj < 0.05 & log2FoldChange <= -0.58 ~ "Down in CS",
      TRUE ~ "Not significant"
    )
  )

labels <- bind_rows(
  d |> filter(direction == "Up in CS") |> arrange(padj) |> slice_head(n = 6),
  d |> filter(direction == "Down in CS") |> arrange(padj) |> slice_head(n = 6)
)

p <- ggplot(d, aes(log2FoldChange, neg_log10_padj)) +
  geom_point(aes(color = direction), size = 0.85, alpha = 0.8) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dotdash", color = "#9E9E9E") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotdash", color = "#9E9E9E") +
  geom_label_repel(
    data = labels, aes(label = gene), color = "black", fill = "white",
    fontface = "bold", size = 3.4, box.padding = 0.45, point.padding = 0.2,
    min.segment.length = 0, max.overlaps = Inf, seed = 60, show.legend = FALSE
  ) +
  scale_color_manual(values = c("Up in CS" = "#D62728", "Down in CS" = "#4D4398", "Not significant" = "#B5B5B5")) +
  labs(
    title = paste0(sex, ": dMIC60-CS vs dMIC60-WT"),
    subtitle = sprintf(
      "|log2FC| >= 0.58; adjusted P < 0.05; down: %d; up: %d",
      sum(d$direction == "Down in CS"), sum(d$direction == "Up in CS")
    ),
    x = "log2 fold change (CS/WT)", y = "-log10 adjusted P value", color = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "bottom")

for (out_dir in c("Rebuilt_Output", "Final_Graphs")) {
  dir.create(file.path(panel_dir, out_dir), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(panel_dir, out_dir, paste0(output_stem, ".png")), p, width = 7, height = 8, dpi = 600, bg = "white")
  ggsave(file.path(panel_dir, out_dir, paste0(output_stem, ".pdf")), p, width = 7, height = 8, device = "pdf", useDingbats = FALSE)
}

write.csv(
  data.frame(sex, testable_genes = nrow(d), up_in_CS = sum(d$direction == "Up in CS"), down_in_CS = sum(d$direction == "Down in CS")),
  file.path(panel_dir, "Supporting_Data", "volcano_DEG_count_summary.csv"), row.names = FALSE
)
