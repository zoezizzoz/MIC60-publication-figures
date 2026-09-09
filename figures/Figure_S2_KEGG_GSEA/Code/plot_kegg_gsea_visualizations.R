#!/usr/bin/env Rscript

# Standard publication visualizations for the completed KEGG GSEA analysis:
# NES overview, running-enrichment curves, and leading-edge fold-change maps.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(enrichplot)
  library(ggplot2)
  library(org.Dm.eg.db)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

set.seed(1)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
panel_dir <- dirname(script_dir)
out_dir <- file.path(panel_dir, "Rebuilt_Output")
plot_dir <- file.path(out_dir, "GSEA_visualizations")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

term2gene <- read_csv(file.path(panel_dir, "Original_Data", "KEGG_Drosophila_TERM2GENE_snapshot.csv"), show_col_types = FALSE)
term2name <- read_csv(file.path(panel_dir, "Original_Data", "KEGG_Drosophila_TERM2NAME_snapshot.csv"), show_col_types = FALSE)
mapping <- read_csv(file.path(out_dir, "gene_ID_mapping_and_DE_results.csv"), show_col_types = FALSE) %>%
  mutate(ENTREZID = as.character(ENTREZID))
gsea_results <- read_csv(file.path(out_dir, "KEGG_GSEA_all_results.csv"), show_col_types = FALSE)

red <- "#B2182B"
blue <- "#2166AC"
navy <- "#17365D"
theme_pub <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, colour = navy),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#EAF2F8", colour = navy),
    legend.position = "top"
  )

run_gsea <- function(comparison) {
  ranked <- mapping %>%
    filter(comparison == !!comparison, is.finite(stat), !is.na(ENTREZID)) %>%
    arrange(desc(abs(stat))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  ranks <- sort(setNames(ranked$stat, ranked$ENTREZID), decreasing = TRUE)
  GSEA(
    geneList = ranks,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )
}

gsea_objects <- list(
  CSF_vs_WRF = run_gsea("CSF_vs_WRF"),
  CSM_vs_WRM = run_gsea("CSM_vs_WRM")
)

# -----------------------------------------------------------------------------
# 1. Global NES overview
# -----------------------------------------------------------------------------
nes_plot_data <- gsea_results %>%
  filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  group_by(comparison) %>%
  slice_max(order_by = abs(NES), n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Direction = if_else(NES > 0, "Higher in CS", "Higher in WT"),
    facet_label = recode(comparison, CSF_vs_WRF = "Female: CS-F vs WT-F", CSM_vs_WRM = "Male: CS-M vs WT-M"),
    pathway_key = paste(comparison, str_wrap(Description, 38), sep = "___")
  ) %>%
  arrange(facet_label, NES) %>%
  mutate(pathway_key = factor(pathway_key, levels = unique(pathway_key)))

p_nes <- ggplot(nes_plot_data, aes(NES, pathway_key)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.35) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_key, colour = Direction), linewidth = 0.8) +
  geom_point(aes(size = -log10(p.adjust), colour = Direction)) +
  facet_wrap(~facet_label, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_colour_manual(values = c("Higher in CS" = red, "Higher in WT" = blue)) +
  scale_size_continuous(name = expression(-log[10]~BH~italic(P)), range = c(2.5, 6)) +
  labs(
    title = "KEGG GSEA pathway overview",
    subtitle = "Top significant pathways by |NES|; all pathways shown have BH-adjusted P < 0.05",
    x = "Normalized enrichment score (NES)", y = NULL, colour = NULL
  ) + theme_pub +
  theme(axis.text.y = element_text(size = 9.5), panel.spacing.x = grid::unit(1.2, "lines"))

ggsave(file.path(plot_dir, "KEGG_GSEA_NES_overview.pdf"), p_nes, width = 11, height = 5.25)
ggsave(file.path(plot_dir, "KEGG_GSEA_NES_overview.png"), p_nes, width = 11, height = 5.25, dpi = 300)

if ("--overview-only" %in% commandArgs(trailingOnly = TRUE)) quit(save = "no")

# -----------------------------------------------------------------------------
# 2. Running-enrichment plots
# -----------------------------------------------------------------------------
selected_curves <- tribble(
  ~comparison,   ~pathway_id, ~file_stub,                   ~display_name,
  "CSF_vs_WRF", "dme00020", "female_TCA_cycle",           "Female: citrate cycle (TCA cycle)",
  "CSF_vs_WRF", "dme03010", "female_ribosome",            "Female: ribosome",
  "CSF_vs_WRF", "dme04624", "female_Toll_Imd",            "Female: Toll and Imd signaling",
  "CSM_vs_WRM", "dme04146", "male_peroxisome",            "Male: peroxisome",
  "CSF_vs_WRF", "dme04711", "female_circadian_rhythm",    "Female: circadian rhythm - fly",
  "CSM_vs_WRM", "dme04711", "male_circadian_rhythm",      "Male: circadian rhythm - fly"
)

curve_plots <- list()
curve_summary <- list()

for (i in seq_len(nrow(selected_curves))) {
  sel <- selected_curves[i, ]
  obj <- gsea_objects[[sel$comparison]]
  result_row <- as.data.frame(obj) %>% filter(ID == sel$pathway_id)
  if (nrow(result_row) != 1) stop("Pathway missing or duplicated: ", sel$comparison, " / ", sel$pathway_id)
  direction <- ifelse(result_row$NES > 0, "Higher in CS", "Higher in WT")
  line_colour <- ifelse(result_row$NES > 0, red, blue)
  title <- sprintf(
    "%s | NES %.2f | BH P = %.3g",
    sel$display_name, result_row$NES, result_row$p.adjust
  )
  p <- gseaplot2(
    obj,
    geneSetID = sel$pathway_id,
    title = title,
    color = line_colour,
    base_size = 10.5,
    rel_heights = c(1.5, 0.45, 0.9),
    subplots = 1:3,
    pvalue_table = FALSE,
    ES_geom = "line"
  )
  ggsave(file.path(plot_dir, paste0("GSEA_running_", sel$file_stub, ".pdf")), p, width = 7.4, height = 5.3)
  ggsave(file.path(plot_dir, paste0("GSEA_running_", sel$file_stub, ".png")), p, width = 7.4, height = 5.3, dpi = 300)
  curve_plots[[sel$file_stub]] <- p
  curve_summary[[sel$file_stub]] <- tibble(
    comparison = sel$comparison,
    pathway_id = sel$pathway_id,
    pathway = result_row$Description,
    NES = result_row$NES,
    pvalue = result_row$pvalue,
    p_adjust = result_row$p.adjust,
    direction = direction,
    leading_edge_genes = result_row$gene_symbols
  )
}

write_csv(bind_rows(curve_summary), file.path(plot_dir, "selected_running_plot_statistics.csv"))

# A multi-page PDF preserves the full running-score, hit-index, and rank panels.
pdf(file.path(plot_dir, "KEGG_GSEA_selected_running_plots.pdf"), width = 11, height = 8.5, onefile = TRUE)
for (p in curve_plots) print(p)
dev.off()

# -----------------------------------------------------------------------------
# 3. Leading-edge log2FC heatmap for major significant pathways
# -----------------------------------------------------------------------------
leading_specs <- tribble(
  ~comparison,   ~pathway_id, ~pathway_label,
  "CSF_vs_WRF", "dme00020", "TCA cycle\n(female leading edge)",
  "CSF_vs_WRF", "dme03010", "Ribosome\n(female leading edge)",
  "CSF_vs_WRF", "dme04624", "Toll/Imd signaling\n(female leading edge)",
  "CSM_vs_WRM", "dme04146", "Peroxisome\n(male leading edge)"
)

leading_gene_sets <- list()
for (i in seq_len(nrow(leading_specs))) {
  spec <- leading_specs[i, ]
  row <- gsea_results %>% filter(comparison == spec$comparison, ID == spec$pathway_id)
  ids <- strsplit(row$core_enrichment, "/", fixed = TRUE)[[1]]
  top_ids <- mapping %>%
    filter(comparison == spec$comparison, ENTREZID %in% ids, is.finite(stat)) %>%
    arrange(desc(abs(stat))) %>%
    distinct(ENTREZID, .keep_all = TRUE) %>%
    slice_head(n = 12) %>%
    pull(ENTREZID)
  leading_gene_sets[[spec$pathway_label]] <- top_ids
}

leading_long <- bind_rows(lapply(names(leading_gene_sets), function(label) {
  ids <- leading_gene_sets[[label]]
  gene_order <- mapping %>%
    filter(ENTREZID %in% ids) %>%
    distinct(ENTREZID, gene) %>%
    mutate(order = match(ENTREZID, ids)) %>%
    arrange(order)
  tidyr::expand_grid(pathway = label, ENTREZID = ids, comparison = c("CSF_vs_WRF", "CSM_vs_WRM")) %>%
    left_join(mapping %>% dplyr::select(comparison, ENTREZID, gene, log2FoldChange), by = c("comparison", "ENTREZID")) %>%
    mutate(
      gene = coalesce(gene, gene_order$gene[match(ENTREZID, gene_order$ENTREZID)]),
      sex_label = recode(comparison, CSF_vs_WRF = "Female", CSM_vs_WRM = "Male"),
      gene_key = paste(pathway, gene, sep = "___"),
      gene_order = match(ENTREZID, ids)
    )
}))

leading_levels <- leading_long %>%
  distinct(pathway, gene_key, gene_order) %>%
  arrange(pathway, desc(gene_order)) %>%
  pull(gene_key)
leading_long <- leading_long %>% mutate(gene_key = factor(gene_key, levels = unique(leading_levels)))
heat_limit <- max(2, quantile(abs(leading_long$log2FoldChange), 0.95, na.rm = TRUE))

p_leading <- ggplot(leading_long, aes(sex_label, gene_key, fill = log2FoldChange)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(log2FoldChange), "NA", sprintf("%.2f", log2FoldChange))), size = 2.6) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_fill_gradient2(
    low = blue, mid = "white", high = red, midpoint = 0,
    limits = c(-heat_limit, heat_limit), oob = squish,
    na.value = "grey85", name = "log2FC\n(CS / WT)"
  ) +
  labs(
    title = "Leading-edge genes driving major KEGG GSEA signals",
    subtitle = "Top 12 leading-edge genes by |DESeq2 Wald statistic|; values are log2 fold changes",
    x = NULL, y = NULL
  ) + theme_pub +
  theme(axis.text.y = element_text(size = 7.5), axis.text.x = element_text(face = "bold"))

ggsave(file.path(plot_dir, "KEGG_GSEA_leading_edge_log2FC_heatmap.pdf"), p_leading, width = 9.5, height = 11)
ggsave(file.path(plot_dir, "KEGG_GSEA_leading_edge_log2FC_heatmap.png"), p_leading, width = 9.5, height = 11, dpi = 300)
write_csv(leading_long %>% dplyr::select(pathway, ENTREZID, gene, comparison, sex_label, log2FoldChange),
          file.path(plot_dir, "leading_edge_heatmap_values.csv"))

# -----------------------------------------------------------------------------
# 4. Circadian KEGG-member fold-change heatmap
# -----------------------------------------------------------------------------
circadian_ids <- term2gene %>% filter(term == "dme04711") %>% pull(ENTREZID) %>% as.character()
circadian_symbols <- suppressMessages(suppressWarnings(mapIds(
  org.Dm.eg.db, keys = circadian_ids, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)))
circadian <- expand_grid(
  ENTREZID = circadian_ids,
  comparison = c("CSF_vs_WRF", "CSM_vs_WRM")
) %>%
  left_join(mapping %>% dplyr::select(comparison, ENTREZID, log2FoldChange, stat, padj), by = c("comparison", "ENTREZID")) %>%
  mutate(
    gene = unname(circadian_symbols[ENTREZID]),
    sex_label = recode(comparison, CSF_vs_WRF = "Female", CSM_vs_WRM = "Male")
  )
circadian_order <- circadian %>%
  group_by(gene) %>%
  summarise(max_abs_stat = if (all(is.na(stat))) -Inf else max(abs(stat), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_stat) %>% pull(gene)
circadian <- circadian %>% mutate(gene = factor(gene, levels = unique(circadian_order)))

p_circadian <- ggplot(circadian, aes(sex_label, gene, fill = log2FoldChange)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(log2FoldChange), "NA", sprintf("%.2f", log2FoldChange))), size = 3) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0,
                       limits = c(-5.5, 5.5), oob = squish, na.value = "grey85",
                       name = "log2FC\n(CS / WT)") +
  labs(
    title = "KEGG circadian rhythm - fly: member-gene fold changes",
    subtitle = "Not significant in either sex (BH P > 0.75); grey = no estimable fold change",
    x = NULL, y = NULL
  ) + theme_pub +
  theme(axis.text.x = element_text(face = "bold"), axis.text.y = element_text(size = 9))

ggsave(file.path(plot_dir, "KEGG_circadian_member_log2FC_heatmap.pdf"), p_circadian, width = 6.5, height = 7.5)
ggsave(file.path(plot_dir, "KEGG_circadian_member_log2FC_heatmap.png"), p_circadian, width = 6.5, height = 7.5, dpi = 300)
write_csv(circadian, file.path(plot_dir, "KEGG_circadian_member_values.csv"))

cat("Visualization directory: ", plot_dir, "\n", sep = "")
cat("Created 1 NES overview, 6 running-enrichment plots, 2 heatmaps, and supporting CSVs.\n")
