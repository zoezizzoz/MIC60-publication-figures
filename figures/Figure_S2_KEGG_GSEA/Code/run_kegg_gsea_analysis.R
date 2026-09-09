#!/usr/bin/env Rscript

# KEGG over-representation analysis (ORA) and preranked KEGG GSEA for the
# female and male CS-versus-WR DESeq2 result tables.

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(org.Dm.eg.db)
  library(readr)
  library(stringr)
  library(tidyr)
})

set.seed(1)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
panel_dir <- dirname(script_dir)
input_dir <- file.path(panel_dir, "Original_Data")
output_dir <- file.path(panel_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Download only public KEGG annotation tables when no local snapshot exists.
# All study gene lists/statistics remain local; enricher()/GSEA() test locally.
term2gene_file <- file.path(input_dir, "KEGG_Drosophila_TERM2GENE_snapshot.csv")
term2name_file <- file.path(input_dir, "KEGG_Drosophila_TERM2NAME_snapshot.csv")
if (file.exists(term2gene_file) && file.exists(term2name_file)) {
  kegg_term2gene <- readr::read_csv(term2gene_file, show_col_types = FALSE)
  kegg_term2name <- readr::read_csv(term2name_file, show_col_types = FALSE)
} else {
  kegg_links <- KEGGREST::keggLink("pathway", "dme")
  kegg_names <- KEGGREST::keggList("pathway", "dme")
  kegg_to_entrez <- KEGGREST::keggConv("ncbi-geneid", "dme")

  kegg_gene_map <- tibble(
    kegg_gene = names(kegg_to_entrez),
    ENTREZID = sub("^ncbi-geneid:", "", unname(kegg_to_entrez))
  )
  kegg_term2gene <- tibble(
    kegg_gene = names(kegg_links),
    term = sub("^path:", "", unname(kegg_links))
  ) %>%
    inner_join(kegg_gene_map, by = "kegg_gene") %>%
    distinct(term, ENTREZID) %>%
    select(term, ENTREZID)
  kegg_term2name <- tibble(
    term = names(kegg_names),
    name = sub(" - Drosophila melanogaster \\(fruit fly\\)$", "", unname(kegg_names))
  )

  readr::write_csv(kegg_term2gene, term2gene_file)
  readr::write_csv(kegg_term2name, term2name_file)
}

comparisons <- tibble::tribble(
  ~comparison,   ~sex,      ~input_file,
  "CSF_vs_WRF", "Female",  file.path(input_dir, "DEG_CSF_vs_WRF.csv"),
  "CSM_vs_WRM", "Male",    file.path(input_dir, "DEG_CSM_vs_WRM.csv")
)

padj_cutoff <- 0.05
lfc_cutoff <- 1
min_gs_size <- 10
max_gs_size <- 500

theme_pub <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )

map_symbols <- function(symbols) {
  suppressWarnings(AnnotationDbi::mapIds(
    org.Dm.eg.db,
    keys = unique(symbols),
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  ))
}

entrez_to_symbol_string <- function(x) {
  ids <- unique(unlist(strsplit(x, "/", fixed = TRUE)))
  syms <- suppressMessages(suppressWarnings(AnnotationDbi::mapIds(
    org.Dm.eg.db, keys = ids, column = "SYMBOL",
    keytype = "ENTREZID", multiVals = "first"
  )))
  paste(unique(na.omit(unname(syms))), collapse = "/")
}

tidy_kegg <- function(x, comparison, sex, analysis, direction = NA_character_) {
  out <- as.data.frame(x)
  if (nrow(out) == 0) {
    return(tibble(
      comparison = character(), sex = character(), analysis = character(),
      direction = character(), ID = character(), Description = character()
    ))
  }
  gene_column <- if (identical(analysis, "ORA")) "geneID" else "core_enrichment"
  out <- as_tibble(out) %>%
    mutate(
      comparison = comparison,
      sex = sex,
      analysis = analysis,
      direction = direction,
      gene_symbols = vapply(.data[[gene_column]], entrez_to_symbol_string, character(1)),
      significant_BH_0.05 = !is.na(p.adjust) & p.adjust < padj_cutoff
    ) %>%
    relocate(comparison, sex, analysis, direction)
  out
}

all_ora <- list()
all_gsea <- list()
all_mapping <- list()
qc_rows <- list()

for (i in seq_len(nrow(comparisons))) {
  cmp <- comparisons[i, ]
  message("Analyzing ", cmp$comparison, " ...")

  deg <- readr::read_csv(cmp$input_file, show_col_types = FALSE) %>%
    transmute(
      gene = as.character(gene),
      baseMean = as.numeric(baseMean),
      log2FoldChange = as.numeric(log2FoldChange),
      lfcSE = as.numeric(lfcSE),
      stat = as.numeric(stat),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj)
    )

  symbol_map <- map_symbols(deg$gene)
  mapped <- deg %>%
    mutate(ENTREZID = unname(symbol_map[gene]))

  # Preserve the strongest absolute Wald statistic where multiple symbols map
  # to the same Entrez ID.
  mapped_unique <- mapped %>%
    filter(!is.na(ENTREZID), is.finite(stat)) %>%
    arrange(desc(abs(stat))) %>%
    distinct(ENTREZID, .keep_all = TRUE)

  universe <- unique(mapped$ENTREZID[!is.na(mapped$ENTREZID) & !is.na(mapped$padj)])
  up_ids <- mapped %>%
    filter(!is.na(ENTREZID), !is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(ENTREZID) %>% unique()
  down_ids <- mapped %>%
    filter(!is.na(ENTREZID), !is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(ENTREZID) %>% unique()

  all_mapping[[cmp$comparison]] <- mapped %>%
    mutate(comparison = cmp$comparison, sex = cmp$sex) %>%
    relocate(comparison, sex)

  for (direction_name in c("Up_in_CS", "Down_in_CS")) {
    ids <- if (direction_name == "Up_in_CS") up_ids else down_ids
    ek <- enricher(
      gene = ids,
      universe = universe,
      TERM2GENE = kegg_term2gene,
      TERM2NAME = kegg_term2name,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size
    )
    all_ora[[paste(cmp$comparison, direction_name, sep = "__")]] <-
      tidy_kegg(ek, cmp$comparison, cmp$sex, "ORA", direction_name)
  }

  ranks <- mapped_unique$stat
  names(ranks) <- mapped_unique$ENTREZID
  ranks <- sort(ranks, decreasing = TRUE)

  gk <- GSEA(
    geneList = ranks,
    TERM2GENE = kegg_term2gene,
    TERM2NAME = kegg_term2name,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )
  all_gsea[[cmp$comparison]] <- tidy_kegg(gk, cmp$comparison, cmp$sex, "GSEA") %>%
    mutate(direction = if_else(NES > 0, "Higher_in_CS", "Higher_in_WR"))

  qc_rows[[cmp$comparison]] <- tibble(
    comparison = cmp$comparison,
    sex = cmp$sex,
    input_file = cmp$input_file,
    input_rows = nrow(deg),
    genes_with_padj = sum(!is.na(deg$padj)),
    mapped_rows = sum(!is.na(mapped$ENTREZID)),
    unique_ranked_entrez = length(ranks),
    background_entrez = length(universe),
    up_DEG_entrez = length(up_ids),
    down_DEG_entrez = length(down_ids)
  )
}

ora <- bind_rows(all_ora) %>% arrange(comparison, direction, p.adjust, pvalue)
gsea <- bind_rows(all_gsea) %>% arrange(comparison, p.adjust, desc(abs(NES)))
mapping <- bind_rows(all_mapping)
qc <- bind_rows(qc_rows)

readr::write_csv(ora, file.path(output_dir, "KEGG_ORA_all_results.csv"), na = "")
readr::write_csv(filter(ora, significant_BH_0.05), file.path(output_dir, "KEGG_ORA_BH_significant.csv"), na = "")
readr::write_csv(gsea, file.path(output_dir, "KEGG_GSEA_all_results.csv"), na = "")
readr::write_csv(filter(gsea, significant_BH_0.05), file.path(output_dir, "KEGG_GSEA_BH_significant.csv"), na = "")
readr::write_csv(mapping, file.path(output_dir, "gene_ID_mapping_and_DE_results.csv"), na = "")
readr::write_csv(qc, file.path(output_dir, "analysis_QC_summary.csv"), na = "")

# ORA figure: up to eight BH-significant pathways per comparison and direction.
ora_plot_data <- ora %>%
  filter(significant_BH_0.05) %>%
  mutate(GeneRatio_numeric = Count / as.numeric(sub(".*/", "", GeneRatio))) %>%
  group_by(comparison, direction) %>%
  slice_min(order_by = p.adjust, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    comparison = factor(comparison, levels = comparisons$comparison),
    direction = factor(direction, levels = c("Up_in_CS", "Down_in_CS")),
    pathway = str_wrap(Description, width = 38)
  )

if (nrow(ora_plot_data) > 0) {
  p_ora <- ggplot(ora_plot_data, aes(GeneRatio_numeric, reorder(pathway, GeneRatio_numeric))) +
    geom_point(aes(size = Count, colour = -log10(p.adjust))) +
    facet_grid(direction ~ comparison, scales = "free_y", space = "free_y") +
    scale_colour_viridis_c(option = "magma", direction = -1, name = expression(-log[10]~BH~italic(P))) +
    labs(
      title = "KEGG over-representation among differentially expressed genes",
      subtitle = "DEG: BH-adjusted P < 0.05 and |log2 fold change| > 1; tested-gene background",
      x = "Gene ratio", y = NULL, size = "DEG count"
    ) + theme_pub
  ggsave(file.path(output_dir, "KEGG_ORA_significant_dotplot.pdf"), p_ora, width = 10, height = 7.5)
  ggsave(file.path(output_dir, "KEGG_ORA_significant_dotplot.png"), p_ora, width = 10, height = 7.5, dpi = 300)
}

# GSEA figure: up to eight BH-significant pathways in each direction/comparison.
gsea_plot_data <- gsea %>%
  filter(significant_BH_0.05) %>%
  group_by(comparison, direction) %>%
  slice_max(order_by = abs(NES), n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    comparison = factor(comparison, levels = comparisons$comparison),
    pathway = str_wrap(Description, width = 42)
  )

if (nrow(gsea_plot_data) > 0) {
  p_gsea <- ggplot(gsea_plot_data, aes(NES, reorder(pathway, NES), fill = direction)) +
    geom_col(width = 0.72) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    facet_wrap(~comparison, scales = "free_y") +
    scale_fill_manual(values = c(Higher_in_CS = "#B2182B", Higher_in_WR = "#2166AC")) +
    labs(
      title = "Preranked KEGG gene-set enrichment analysis",
      subtitle = "All mapped genes ranked by the DESeq2 Wald statistic; pathways shown at BH-adjusted P < 0.05",
      x = "Normalized enrichment score (NES)", y = NULL, fill = NULL
    ) + theme_pub
  ggsave(file.path(output_dir, "KEGG_GSEA_significant_NES.pdf"), p_gsea, width = 10, height = 7.5)
  ggsave(file.path(output_dir, "KEGG_GSEA_significant_NES.png"), p_gsea, width = 10, height = 7.5, dpi = 300)
}

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "sessionInfo.txt"))

cat("\nQC summary:\n")
print(qc)
cat("\nBH-significant ORA pathways:\n")
print(ora %>% count(comparison, direction, wt = significant_BH_0.05, name = "n_significant"))
cat("\nBH-significant GSEA pathways:\n")
print(gsea %>% count(comparison, direction, wt = significant_BH_0.05, name = "n_significant"))
cat("\nOutput directory: ", output_dir, "\n", sep = "")
