#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
panel_dir <- dirname(dirname(normalizePath(sub("^--file=", "", args[[1]]))))
d <- read_csv(file.path(panel_dir, "Supporting_Data", "FigS2_plotted_pathways.csv"), show_col_types = FALSE) |>
  mutate(
    sex = factor(sex, levels = c("Female", "Male")),
    direction = factor(direction, levels = c("Higher_in_WT", "Higher_in_CS")),
    pathway = str_wrap(Description, 34),
    pathway_key = factor(paste(sex, pathway, sep = "___"), levels = unique(paste(sex, pathway, sep = "___")))
  )

p <- ggplot(d, aes(NES, pathway_key, color = direction)) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_key), linewidth = 0.7) +
  geom_point(aes(size = -log10(p.adjust))) +
  facet_wrap(~sex, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_color_manual(values = c(Higher_in_WT = "#4C72B0", Higher_in_CS = "#A8322B"), labels = c("Higher in WT", "Higher in CS")) +
  scale_size_continuous(name = expression(-log[10] ~ "adjusted P"), range = c(2, 5.5)) +
  labs(x = "Normalized enrichment score (NES)", y = NULL, color = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.text = element_text(face = "bold"), legend.position = "top")

for (out_dir in c("Rebuilt_Output", "Final_Graphs")) {
  dir.create(file.path(panel_dir, out_dir), recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(panel_dir, out_dir, "Figure_S2_KEGG_GSEA.png"), p, width = 10, height = 6.2, dpi = 600, bg = "white")
  ggsave(file.path(panel_dir, out_dir, "Figure_S2_KEGG_GSEA.pdf"), p, width = 10, height = 6.2, device = "pdf", useDingbats = FALSE)
}
