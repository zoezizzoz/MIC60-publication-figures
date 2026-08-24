#!/usr/bin/env Rscript

# Rebuild all RNA-seq-derived panels used in Fig1, Fig2, and FigS1 with CSF3
# excluded before every statistical and visualization step.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(AnnotationDbi)
  library(ggrepel)
  library(pheatmap)
  library(patchwork)
  library(Rtsne)
  library(uwot)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
package_dir <- normalizePath(file.path(figure_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
set.seed(FIG_SEED)

count_path <- file.path(figure_dir, "Original_Data", "count_matrix_symbol.csv")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(output_dir, name)

excluded_samples <- "CSF3"
padj_cutoff <- 0.05
lfc_cutoff <- 1

# -----------------------------------------------------------------------------
# One shared CSF3-excluded DESeq2 object drives every panel below.
# -----------------------------------------------------------------------------
count_mat <- read.csv(count_path, row.names = 1, check.names = FALSE)
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"
if (!all(excluded_samples %in% colnames(count_mat))) {
  stop("Excluded sample not found: ", excluded_samples)
}
count_mat <- count_mat[, !colnames(count_mat) %in% excluded_samples, drop = FALSE]

samples <- colnames(count_mat)
geno <- ifelse(grepl("^CS", samples), "CS", "WR")
sex <- ifelse(substr(samples, 3, 3) == "F", "F", "M")
coldata <- data.frame(
  row.names = samples,
  genotype = factor(geno, levels = c("WR", "CS")),
  sex = factor(sex, levels = c("F", "M")),
  group = factor(paste0(geno, sex), levels = c("WRF", "CSF", "WRM", "CSM"))
)

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)
res_f <- results(dds, contrast = c("group", "CSF", "WRF"))
res_m <- results(dds, contrast = c("group", "CSM", "WRM"))
norm_counts <- counts(dds, normalized = TRUE)
vsd <- vst(dds, blind = TRUE)

format_result <- function(x) {
  d <- as.data.frame(x)
  d$gene <- rownames(d)
  d[order(d$padj), ]
}
df_f <- format_result(res_f)
df_m <- format_result(res_m)

write.csv(df_f, out("deseq2_female_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)
write.csv(df_m, out("deseq2_male_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)
write.csv(norm_counts, out("normalized_counts_excluding_CSF3.csv"))
saveRDS(
  list(dds = dds, vsd = vsd, norm = norm_counts, coldata = coldata,
       df_f = df_f, df_m = df_m, excluded_samples = excluded_samples),
  out("deseq2_objects_excluding_CSF3.rds")
)

group_n <- table(coldata$group)
n_group <- function(x) unname(as.integer(group_n[x]))
all_caption <- paste0(
  "Drosophila bulk RNA-seq biological libraries: female WT n = ", n_group("WRF"),
  ", female CS n = ", n_group("CSF"), ", male WT n = ", n_group("WRM"),
  ", male CS n = ", n_group("CSM"), ". CSF3 excluded before all analysis steps."
)
female_caption <- paste0(
  "Female Drosophila bulk RNA-seq: WT n = ", n_group("WRF"),
  ", CS n = ", n_group("CSF"),
  " biological libraries. CSF3 excluded before all analysis steps."
)
male_caption <- paste0(
  "Male Drosophila bulk RNA-seq: WT n = ", n_group("WRM"),
  ", CS n = ", n_group("CSM"),
  " biological libraries; CSF3 excluded from the shared DESeq2 model."
)

# -----------------------------------------------------------------------------
# Fig1B: PCA, MDS, UMAP, and t-SNE from the same 500 most-variable VST genes.
# -----------------------------------------------------------------------------
vst_mat <- assay(vsd)
top_n <- min(500, nrow(vst_mat))
top_genes <- names(sort(matrixStats::rowVars(vst_mat), decreasing = TRUE))[seq_len(top_n)]
x <- t(vst_mat[top_genes, , drop = FALSE])

pca <- prcomp(x, center = TRUE, scale. = FALSE)
pca_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
mds <- cmdscale(dist(x), k = 2, eig = TRUE)
set.seed(FIG_SEED)
umap_xy <- uwot::umap(
  x, n_neighbors = 5, min_dist = 0.30, n_components = 2,
  metric = "euclidean", n_epochs = 500, verbose = FALSE,
  ret_model = FALSE, n_threads = 1
)
set.seed(FIG_SEED)
tsne_xy <- Rtsne::Rtsne(
  x, dims = 2, perplexity = 3, check_duplicates = FALSE,
  pca = TRUE, theta = 0.5, max_iter = 1000, verbose = FALSE
)$Y

qc_meta <- tibble(Sample = rownames(x)) %>%
  mutate(
    Genotype = factor(ifelse(grepl("^CS", Sample), "CS", "WT"),
                      levels = c("WT", "CS")),
    Sex = factor(ifelse(substr(Sample, 3, 3) == "F", "Female", "Male"),
                 levels = c("Female", "Male"))
  )

embedding_plot <- function(coords, xlab, ylab, title, equal_units = FALSE) {
  d <- bind_cols(qc_meta, tibble(x = coords[, 1], y = coords[, 2]))
  p <- ggplot(d, aes(x, y)) +
    geom_point(
      aes(fill = Genotype, shape = Sex),
      color = FIG_PT_COLOR, stroke = FIG_MARKER_STROKE,
      size = FIG_PT_SIZE + 0.8, alpha = FIG_PT_ALPHA
    ) +
    ggrepel::geom_text_repel(
      aes(label = Sample), color = FIG_PT_COLOR, size = FIG_ANNOT_SIZE,
      box.padding = 0.35, point.padding = 0.25, min.segment.length = 0,
      segment.color = FIG_REF_LINE_COLOR,
      segment.size = FIG_REF_LINE_LINEWIDTH, seed = FIG_SEED,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    fig_scale_fill(c("WT", "CS"), labels = c("WT", "CS")) +
    scale_shape_manual(
      values = c("Female" = 21, "Male" = 22),
      labels = c("Female" = "\u2640 Female", "Male" = "\u2642 Male")
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, size = FIG_PT_SIZE + 0.8)
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(fill = "white", size = FIG_PT_SIZE + 0.8)
      )
    ) +
    labs(x = xlab, y = ylab, title = title) +
    theme_fig(legend_position = "bottom") +
    theme(
      plot.title = element_text(size = FIG_SUBTITLE_SIZE, face = "bold"),
      legend.title = element_blank()
    )
  if (equal_units) {
    axis_limit <- max(abs(coords), na.rm = TRUE) * 1.08
    p <- p + coord_fixed(
      ratio = 1,
      xlim = c(-axis_limit, axis_limit),
      ylim = c(-axis_limit, axis_limit)
    )
  }
  p
}

p_pca <- embedding_plot(
  pca$x[, 1:2], paste0("PC1 (", pca_var[1], "%)"),
  paste0("PC2 (", pca_var[2], "%)"), "PCA", equal_units = TRUE
)
p_mds <- embedding_plot(mds$points, "MDS1", "MDS2", "MDS", equal_units = TRUE)
p_umap <- embedding_plot(umap_xy, "UMAP1", "UMAP2", "UMAP")
p_tsne <- embedding_plot(tsne_xy, "t-SNE1", "t-SNE2", "t-SNE")

qc_four <- (p_pca + p_mds + p_umap + p_tsne) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "RNA-seq sample-level quality control",
    caption = paste0(
      "Drosophila bulk RNA-seq: female WT n = 3, CS n = 2; ",
      "male WT n = 3, CS n = 3 biological libraries.\n",
      "CSF3 excluded before all analysis steps. DESeq2 VST counts; ",
      "500 most-variable genes."
    ),
    theme = theme(
      plot.title = element_text(
        family = FIG_FONT, size = FIG_TITLE_SIZE, face = "bold", hjust = 0
      ),
      plot.caption = element_text(
        family = FIG_FONT, size = FIG_CAPTION_SIZE, color = "grey35", hjust = 0
      )
    )
  ) & theme(legend.position = "bottom")
fig_save(qc_four, out("fig1b_qc_pca"), width = FIG_W_2COL, height = FIG_H_2COL)

# -----------------------------------------------------------------------------
# Volcano plots for Fig1 and FigS1.
# -----------------------------------------------------------------------------
make_volcano <- function(d, title, caption) {
  plot_df <- d %>% mutate(
    neglog10 = -log10(pmax(padj, .Machine$double.xmin)),
    status = case_when(
      !is.na(padj) & padj < padj_cutoff & log2FoldChange > lfc_cutoff ~ "Up",
      !is.na(padj) & padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
      TRUE ~ "Not significant"
    )
  )
  labels <- bind_rows(
    plot_df %>% filter(status == "Up") %>% slice_min(padj, n = 8),
    plot_df %>% filter(status == "Down") %>% slice_min(padj, n = 8)
  ) %>% distinct(gene, .keep_all = TRUE)
  counts <- table(factor(plot_df$status, levels = c("Down", "Not significant", "Up")))

  ggplot(plot_df, aes(log2FoldChange, neglog10, color = status)) +
    geom_point(size = 1.6, alpha = 0.70) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", color = "grey55", linewidth = 0.45) +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "grey55", linewidth = 0.45) +
    ggrepel::geom_label_repel(
      data = labels, aes(label = gene), size = FIG_ANNOT_SIZE - 1.0,
      color = "black", fill = "white", label.size = 0.20,
      box.padding = 0.35, point.padding = 0.20, min.segment.length = 0,
      segment.size = 0.22, seed = FIG_SEED, max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("Down" = unname(FIG_COLORS["control"]),
                 "Not significant" = "grey75",
                 "Up" = unname(FIG_COLORS["mutant"])),
      breaks = c("Down", "Not significant", "Up")
    ) +
    labs(
      x = expression(log[2] ~ fold ~ change ~ (CS/WT)),
      y = expression(-log[10] ~ adjusted ~ italic(p)),
      title = title,
      subtitle = paste0(
        "Down: ", counts[["Down"]], "   Up: ", counts[["Up"]],
        "   thresholds: |log2FC| > ", lfc_cutoff, ", padj < ", padj_cutoff
      ),
      caption = caption,
      color = NULL
    ) +
    theme_fig(legend_position = "top")
}

fig_save(
  make_volcano(df_f, "Female CS versus WT differential expression", female_caption),
  out("fig1c_volcano_female"), width = 6.6, height = 5.3
)
fig_save(
  make_volcano(df_m, "Male CS versus WT differential expression", male_caption),
  out("figS1_volcano_male"), width = 6.6, height = 5.3
)

# -----------------------------------------------------------------------------
# Heatmaps for Fig1 and FigS1.
# -----------------------------------------------------------------------------
heatmap_genes <- function(d, n_each = 30) {
  bind_rows(
    d %>% filter(!is.na(padj), log2FoldChange > 0) %>% slice_min(padj, n = n_each),
    d %>% filter(!is.na(padj), log2FoldChange < 0) %>% slice_min(padj, n = n_each)
  ) %>% distinct(gene, .keep_all = TRUE) %>% pull(gene)
}

save_heatmap <- function(d, base_name, title) {
  genes <- heatmap_genes(d)
  z <- assay(vsd)[genes, , drop = FALSE]
  z <- t(scale(t(z)))
  z[!is.finite(z)] <- 0
  ann <- data.frame(Group = factor(coldata[colnames(z), "group"],
                                   levels = c("CSF", "CSM", "WRF", "WRM")))
  rownames(ann) <- colnames(z)
  ann_cols <- list(Group = c(CSF = "#F8766D", CSM = "#7CAE00",
                             WRF = "#00BFC4", WRM = "#C77CFF"))
  heat_cols <- colorRampPalette(c(
    unname(FIG_COLORS["control"]), "white", unname(FIG_COLORS["mutant"])
  ))(101)
  for (ext in c("png", "pdf")) {
    pheatmap::pheatmap(
      z, color = heat_cols, cluster_rows = TRUE, cluster_cols = TRUE,
      annotation_col = ann, annotation_colors = ann_cols,
      show_colnames = TRUE, show_rownames = TRUE,
      fontsize = FIG_LEGEND_TEXT_SIZE,
      fontsize_row = max(5, FIG_LEGEND_TEXT_SIZE - 1),
      border_color = NA,
      main = paste0(title, "\nCSF3 excluded before all analysis steps"),
      filename = paste0(out(base_name), ".", ext),
      width = 7.0, height = 8.4
    )
  }
}
save_heatmap(df_f, "fig1d_heatmap_female", "Top female CS-versus-WT genes")
save_heatmap(df_m, "figS1_heatmap_male", "Top male CS-versus-WT genes")

# -----------------------------------------------------------------------------
# GSEA/GO panels for Fig2 and FigS1.
# -----------------------------------------------------------------------------
run_gsea <- function(d) {
  ranks <- d %>% filter(is.finite(stat)) %>% dplyr::select(gene, stat)
  ranks$ENTREZID <- mapIds(
    org.Dm.eg.db, keys = ranks$gene, column = "ENTREZID",
    keytype = "SYMBOL", multiVals = "first"
  )
  ranks <- ranks %>%
    filter(!is.na(ENTREZID)) %>% arrange(desc(abs(stat))) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  gene_list <- ranks$stat
  names(gene_list) <- ranks$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  set.seed(FIG_SEED)
  gseGO(
    geneList = gene_list, OrgDb = org.Dm.eg.db, keyType = "ENTREZID",
    ont = "ALL", minGSSize = 10, maxGSSize = 500,
    pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE
  )
}

gsea_f <- run_gsea(df_f)
gsea_m <- run_gsea(df_m)
write.csv(as.data.frame(gsea_f), out("GSEA_GO_female_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_m), out("GSEA_GO_male_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)

gsea_dot <- function(gsea_obj, title, caption, n_each = 7) {
  d <- as.data.frame(gsea_obj) %>%
    filter(!is.na(p.adjust), p.adjust < padj_cutoff) %>%
    mutate(Direction = ifelse(NES > 0, "Higher in CS", "Higher in WT"))
  if (nrow(d) == 0) stop("No significant GSEA terms for: ", title)
  selected <- d %>%
    group_by(Direction) %>%
    arrange(p.adjust, desc(abs(NES)), .by_group = TRUE) %>%
    slice_head(n = n_each) %>% ungroup() %>%
    arrange(NES) %>%
    mutate(Display = stringr::str_wrap(Description, width = 38))
  selected$Display <- factor(selected$Display, levels = unique(selected$Display))
  ggplot(selected, aes(NES, Display)) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.45) +
    geom_segment(aes(x = 0, xend = NES, yend = Display),
                 color = "grey80", linewidth = 0.6) +
    geom_point(aes(size = setSize, color = -log10(p.adjust))) +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "#92C5DE", high = "#B2182B",
                         name = expression(-log[10] ~ adjusted ~ italic(p))) +
    scale_size_continuous(range = c(3, 8), name = "Gene-set size") +
    labs(
      x = "Normalized enrichment score (NES)", y = NULL,
      title = title,
      subtitle = "Rank-based GSEA across all tested genes; BH-adjusted p < 0.05",
      caption = caption
    ) +
    theme_fig(legend_position = "right") +
    theme(
      strip.background = element_rect(fill = FIG_FACET_STRIP_FILL,
                                      color = FIG_FACET_BORDER_COLOR),
      strip.text = element_text(face = "bold")
    )
}

fig_save(
  gsea_dot(gsea_f, "Female CS-versus-WT biological programs", female_caption),
  out("fig2a_GO_female"), width = 9.0, height = 7.2
)
fig_save(
  gsea_dot(gsea_m, "Male CS-versus-WT biological programs", male_caption),
  out("figS1_GO_male"), width = 9.0, height = 7.2
)
rm(gsea_f, gsea_m)
invisible(gc())

# -----------------------------------------------------------------------------
# Fig2B: selected individual genes using CSF3-excluded normalized counts.
# -----------------------------------------------------------------------------
genes_show <- intersect(c("timeout", "ImpL2", "DNAlig3"), rownames(norm_counts))
long <- bind_rows(lapply(genes_show, function(g) {
  data.frame(
    gene = g, count = norm_counts[g, ],
    genotype = coldata$genotype, sex = coldata$sex
  )
}))
long$Genotype <- factor(
  ifelse(long$genotype == "CS", "dMIC60-CS", "dMIC60-WT"),
  levels = c("dMIC60-WT", "dMIC60-CS")
)
long$Sex <- factor(ifelse(long$sex == "F", "Female", "Male"),
                   levels = c("Female", "Male"))
long$gene <- factor(long$gene, levels = genes_show)

sig_ann <- bind_rows(lapply(genes_show, function(g) {
  bind_rows(lapply(c("Female", "Male"), function(s) {
    result_df <- if (s == "Female") df_f else df_m
    vals <- long$count[long$gene == g & long$Sex == s]
    p <- result_df$padj[match(g, result_df$gene)]
    span <- max(diff(range(vals)), abs(max(vals)) * 0.15, 1e-9)
    data.frame(
      gene = g, Sex = s, padj = p,
      y = max(vals) + span * FIG_BRACKET_OFFSET,
      tick = span * FIG_BRACKET_TICK
    )
  }))
}))
sig_ann <- sig_ann %>% filter(is.finite(padj), padj < FIG_ALPHA)
sig_ann$label <- vapply(sig_ann$padj, function(p) {
  p_text <- if (p < 0.001) "p < 0.001" else if (p < 0.01) {
    sprintf("p = %.3f", p)
  } else {
    sprintf("p = %.2f", p)
  }
  paste0(fig_stars(p), "\n", p_text)
}, character(1))
sig_ann$gene <- factor(sig_ann$gene, levels = genes_show)
sig_ann$Sex <- factor(sig_ann$Sex, levels = c("Female", "Male"))

facet_brackets <- function(d, x1 = 1, x2 = 2) {
  if (nrow(d) == 0) return(NULL)
  d$x1 <- x1
  d$x2 <- x2
  d$xmid <- (x1 + x2) / 2
  list(
    geom_segment(data = d, aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_segment(data = d, aes(x = x1, xend = x1, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_segment(data = d, aes(x = x2, xend = x2, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH),
    geom_text(data = d, aes(x = xmid, y = y + tick * 0.6, label = label),
              inherit.aes = FALSE, vjust = 0, size = FIG_ANNOT_SIZE,
              family = FIG_FONT)
  )
}

set.seed(FIG_SEED)
p_genes <- ggplot(long, aes(Genotype, count, fill = Genotype, color = Genotype)) +
  fig_summary() + fig_points() + facet_brackets(sig_ann) +
  scale_fill_manual(
    values = fig_fills(c("dMIC60-WT", "dMIC60-CS")), guide = "none"
  ) +
  scale_color_manual(
    values = fig_colors(c("dMIC60-WT", "dMIC60-CS")), guide = "none"
  ) +
  scale_y_continuous(expand = expansion(
    mult = c(FIG_Y_EXPAND[1], FIG_Y_EXPAND[2] + 0.08)
  )) +
  facet_grid(
    gene ~ Sex, scales = "free_y", switch = "y",
    labeller = labeller(Sex = c(
      "Female" = "\u2640 Female",
      "Male" = "\u2642 Male"
    ))
  ) +
  labs(
    x = NULL, y = "Normalized counts", title = "Individual gene expression",
    caption = paste0(
      "Drosophila bulk RNA-seq: female WT n = 3, CS n = 2; ",
      "male WT n = 3, CS n = 3 biological libraries.\n",
      "CSF3 excluded before filtering, normalization, model fitting, and testing.\n",
      "DESeq2 adjusted p shown only where padj < ", FIG_ALPHA, "."
    )
  ) +
  theme_fig() +
  theme(
    panel.border = element_rect(fill = NA, color = FIG_FACET_BORDER_COLOR,
                                linewidth = FIG_LINE_WIDTH),
    panel.spacing.y = unit(FIG_FACET_GAP_PT, "pt"),
    strip.background = element_rect(fill = FIG_FACET_STRIP_FILL,
                                    color = FIG_FACET_BORDER_COLOR,
                                    linewidth = FIG_LINE_WIDTH),
    strip.text.y.left = element_text(angle = 0, face = "italic"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )
fig_save(p_genes, out("fig2b_genes"), width = 7.2, height = 6.7)

cat("Regenerated CSF3-excluded RNA-seq panels for Fig1, Fig2, and FigS1.\n")
cat("Run generate_string_network.R separately to rebuild STRING network views.\n")
cat(all_caption, "\n")
