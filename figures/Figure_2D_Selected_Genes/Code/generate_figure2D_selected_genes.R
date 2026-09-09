#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(patchwork)
  library(tidyverse)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
set.seed(FIG_SEED)

input_dir <- file.path(figure_dir, "Original_Data")
support_dir <- file.path(figure_dir, "Supporting_Data")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
final_dir <- file.path(figure_dir, "Final_Graphs")
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

count_mat <- read.csv(file.path(input_dir, "count_matrix_symbol.csv"), row.names = 1, check.names = FALSE)
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"
female_de <- read.csv(file.path(input_dir, "DEG_CSF_vs_WRF.csv"), check.names = FALSE)
male_de <- read.csv(file.path(input_dir, "DEG_CSM_vs_WRM.csv"), check.names = FALSE)

fit_sex <- function(sample_names) {
  x <- count_mat[, sample_names, drop = FALSE]
  genotype <- factor(ifelse(grepl("^CS", sample_names), "CS", "WR"), levels = c("WR", "CS"))
  coldata <- data.frame(row.names = sample_names, genotype = genotype)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ genotype)
  dds <- DESeq(dds, quiet = TRUE)
  list(dds = dds, normalized = counts(dds, normalized = TRUE))
}

female_samples <- colnames(count_mat)[grepl("F[123]$", colnames(count_mat))]
male_samples <- colnames(count_mat)[grepl("M[123]$", colnames(count_mat))]
female_fit <- fit_sex(female_samples)
male_fit <- fit_sex(male_samples)

genes_show <- c("timeout", "ImpL2", "DNAlig3")
stopifnot(all(genes_show %in% rownames(count_mat)))

make_long <- function(norm, samples, sex_label) {
  bind_rows(lapply(genes_show, function(g) data.frame(
    gene = g,
    Sample = samples,
    count = as.numeric(norm[g, samples]),
    Genotype = ifelse(grepl("^CS", samples), "dMIC60-CS", "dMIC60-WT"),
    Sex = sex_label,
    stringsAsFactors = FALSE
  )))
}
long <- bind_rows(
  make_long(female_fit$normalized, female_samples, "Female"),
  make_long(male_fit$normalized, male_samples, "Male")
)
long$Genotype <- factor(long$Genotype, levels = c("dMIC60-WT", "dMIC60-CS"))
long$Sex <- factor(long$Sex, levels = c("Female", "Male"))
long$gene <- factor(long$gene, levels = genes_show)
long$PaletteGroup <- factor(
  paste(long$Sex, long$Genotype, sep = " | "),
  levels = c(
    "Female | dMIC60-WT", "Female | dMIC60-CS",
    "Male | dMIC60-WT", "Male | dMIC60-CS"
  )
)

lookup_padj <- function(g, sex_label) {
  d <- if (sex_label == "Female") female_de else male_de
  d$padj[match(g, d$gene)]
}
sig_ann <- bind_rows(lapply(genes_show, function(g) bind_rows(lapply(c("Female", "Male"), function(s) {
  vals <- long$count[long$gene == g & long$Sex == s]
  p <- lookup_padj(g, s)
  span <- max(diff(range(vals)), abs(max(vals)) * 0.15, 1e-9)
  data.frame(gene = g, Sex = s, padj = p,
             y = max(vals) + span * FIG_BRACKET_OFFSET,
             tick = span * FIG_BRACKET_TICK)
})))) %>% filter(is.finite(padj), padj < FIG_ALPHA)
sig_ann$label <- vapply(sig_ann$padj, function(p) {
  if (p < 0.001) "adjusted~italic(p) < 0.001" else if (p < 0.01) sprintf("adjusted~italic(p) == %.3f", p) else sprintf("adjusted~italic(p) == %.2f", p)
}, character(1))
sig_ann$gene <- factor(sig_ann$gene, levels = genes_show)
sig_ann$Sex <- factor(sig_ann$Sex, levels = c("Female", "Male"))

facet_brackets <- function(d, x1 = 1, x2 = 2) {
  if (nrow(d) == 0) return(NULL)
  d$x1 <- x1; d$x2 <- x2; d$xmid <- (x1 + x2) / 2
  list(
    geom_segment(data = d, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_segment(data = d, aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_segment(data = d, aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_text(data = d, aes(x = xmid, y = y + tick * 0.6, label = label), inherit.aes = FALSE,
              vjust = 0, size = FIG_ANNOT_SIZE, family = FIG_FONT, parse = TRUE)
  )
}

p <- ggplot(long, aes(Genotype, count, fill = PaletteGroup)) +
  fig_summary() + fig_points(color = FIG_PT_COLOR) + facet_brackets(sig_ann) +
  scale_fill_manual(values = c(
    "Female | dMIC60-WT" = "#5BB4E5",
    "Female | dMIC60-CS" = "#CC79A7",
    "Male | dMIC60-WT" = "#47599E",
    "Male | dMIC60-CS" = "#891740"
  ), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(FIG_Y_EXPAND[1], FIG_Y_EXPAND[2] + 0.08))) +
  facet_grid(gene ~ Sex, scales = "free_y", switch = "y",
             labeller = labeller(Sex = c("Female" = "\u2640 Female", "Male" = "\u2642 Male"))) +
  labs(
    x = NULL, y = "DESeq2 normalized counts", title = "Individual gene expression",
    caption = paste0(
      "Female and male comparisons were fitted separately with DESeq2; dMIC60-WT n = 3 and dMIC60-CS n = 3 per sex.\n",
      "All libraries, including CSF3, were analyzed without manual gene prefiltering. ",
      "Wald-test P values were BH-adjusted genome-wide; adjusted P values are shown where padj < 0.05."
    )
  ) +
  theme_fig() +
  theme(
    panel.border = element_rect(fill = NA, color = FIG_FACET_BORDER_COLOR, linewidth = FIG_LINE_WIDTH),
    panel.spacing.y = unit(FIG_FACET_GAP_PT, "pt"),
    strip.background = element_rect(fill = FIG_FACET_STRIP_FILL, color = FIG_FACET_BORDER_COLOR, linewidth = FIG_LINE_WIDTH),
    strip.text.y.left = element_text(angle = 0, face = "bold.italic", size = FIG_AXIS_TEXT_SIZE + 2),
    strip.text.x = element_text(size = FIG_AXIS_TEXT_SIZE + 1),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

fig_save(p, file.path(output_dir, "Figure2D_individual_genes"), width = 7.2, height = 6.7)

# Female-only horizontal layout used in the assembled manuscript figure.
# Build it from the same normalized counts, annotations, palette, and figure
# helpers as the complete female/male panel so the two exports stay synchronized.
horizontal_gene_order <- c("timeout", "DNAlig3", "ImpL2")

make_female_horizontal_panel <- function(gene_name, show_y_title = FALSE) {
  panel_data <- droplevels(subset(long, Sex == "Female" & gene == gene_name))
  panel_data$gene <- factor(as.character(panel_data$gene), levels = gene_name)
  panel_ann <- droplevels(subset(sig_ann, Sex == "Female" & gene == gene_name))
  panel_ann$gene <- factor(as.character(panel_ann$gene), levels = gene_name)

  ggplot(panel_data, aes(Genotype, count, fill = PaletteGroup)) +
    fig_summary() +
    fig_points(color = FIG_PT_COLOR) +
    facet_brackets(panel_ann) +
    scale_fill_manual(values = c(
      "Female | dMIC60-WT" = "#5BB4E5",
      "Female | dMIC60-CS" = "#CC79A7"
    ), guide = "none") +
    scale_x_discrete(
      labels = c("dMIC60-WT", "dMIC60-CS"),
      expand = expansion(add = FIG_X_EXPAND)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(FIG_Y_EXPAND[1], FIG_Y_EXPAND[2] + 0.08))
    ) +
    facet_grid(gene ~ ., scales = "free_y", switch = "y") +
    labs(
      x = expression("♀" ~ italic("dMIC60-Null")),
      y = if (show_y_title) "DESeq2 normalized counts" else NULL
    ) +
    theme_fig() +
    theme(
      panel.border = element_rect(
        fill = NA,
        color = FIG_FACET_BORDER_COLOR,
        linewidth = FIG_LINE_WIDTH
      ),
      strip.background = element_rect(
        fill = FIG_FACET_STRIP_FILL,
        color = FIG_FACET_BORDER_COLOR,
        linewidth = FIG_LINE_WIDTH
      ),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold.italic",
        size = FIG_AXIS_TEXT_SIZE + 1
      ),
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        face = "italic",
        lineheight = 0.95
      ),
      axis.title.x = element_text(
        face = "plain",
        size = FIG_AXIS_TEXT_SIZE,
        margin = margin(t = 1)
      ),
      plot.margin = margin(6, 6, 10, 6)
    )
}

horizontal_panels <- lapply(seq_along(horizontal_gene_order), function(i) {
  make_female_horizontal_panel(
    horizontal_gene_order[[i]],
    show_y_title = i == 1L
  )
})

p_horizontal <- wrap_plots(horizontal_panels, nrow = 1)
fig_save(
  p_horizontal,
  file.path(output_dir, "Figure2D_individual_genes_female_horizontal"),
  width = 9.0,
  height = 2.55
)

write.csv(long, file.path(support_dir, "Figure2D_normalized_counts_plotted.csv"), row.names = FALSE)
write.csv(sig_ann, file.path(support_dir, "Figure2D_adjusted_p_annotations.csv"), row.names = FALSE)
write.csv(female_de, file.path(support_dir, "DEG_CSF_vs_WRF.csv"), row.names = FALSE)
write.csv(male_de, file.path(support_dir, "DEG_CSM_vs_WRM.csv"), row.names = FALSE)
file.copy(file.path(output_dir, "Figure2D_individual_genes.png"), file.path(final_dir, "Figure_2D_Selected_Genes_full.png"), overwrite = TRUE)
file.copy(file.path(output_dir, "Figure2D_individual_genes.pdf"), file.path(final_dir, "Figure_2D_Selected_Genes_full.pdf"), overwrite = TRUE)
file.copy(
  file.path(output_dir, "Figure2D_individual_genes_female_horizontal.png"),
  file.path(final_dir, "Figure_2D_Selected_Genes.png"),
  overwrite = TRUE
)
file.copy(
  file.path(output_dir, "Figure2D_individual_genes_female_horizontal.pdf"),
  file.path(final_dir, "Figure_2D_Selected_Genes.pdf"),
  overwrite = TRUE
)
cat("Figure 2D rebuilt from the corrected sex-specific DESeq2 analyses.\n")
