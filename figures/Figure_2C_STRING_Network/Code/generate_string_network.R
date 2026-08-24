#!/usr/bin/env Rscript

# Female dMIC60-CS versus dMIC60-WT STRING views with CSF3 excluded before
# filtering, normalization, model fitting, differential expression, GSEA,
# module selection, and STRING construction.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Dm.eg.db)
  library(AnnotationDbi)
  library(rbioapi)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggforce)
  library(geomtextpath)
  library(cowplot)
  library(scales)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
package_dir <- normalizePath(file.path(figure_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))

count_path <- file.path(figure_dir, "Original_Data", "count_matrix_symbol.csv")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- function(name) file.path(output_dir, name)

excluded_samples <- "CSF3"
padj_cutoff <- 0.05
lfc_cutoff <- 1
string_species <- 7227
string_score <- 400
genes_per_module <- 15
drop_cg <- TRUE
cg_pattern <- "^CG[0-9]+$"

# -----------------------------------------------------------------------------
# Recalculate the analysis after excluding CSF3 from the raw matrix.
# -----------------------------------------------------------------------------
count_mat <- read.csv(count_path, row.names = 1, check.names = FALSE)
count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"
if (!all(excluded_samples %in% colnames(count_mat))) {
  stop("Excluded sample not found in the count matrix: ", excluded_samples)
}
count_mat <- count_mat[, !colnames(count_mat) %in% excluded_samples, drop = FALSE]

samples <- colnames(count_mat)
geno <- ifelse(grepl("^CS", samples), "CS", "WR")
sex <- ifelse(substr(samples, 3, 3) == "F", "F", "M")
coldata <- data.frame(
  row.names = samples,
  genotype = factor(geno, levels = c("WR", "CS")),
  sex = factor(sex, levels = c("F", "M")),
  group = factor(paste0(geno, sex))
)

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)
res_female <- results(dds, contrast = c("group", "CSF", "WRF"))
df_f <- as.data.frame(res_female)
df_f$gene <- rownames(df_f)
df_f <- df_f[order(df_f$padj), ]

write.csv(df_f, out("deseq2_female_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)
write.csv(
  data.frame(sample = rownames(coldata), group = coldata$group),
  out("samples_used_excluding_CSF3.csv"), row.names = FALSE
)

sig_up <- df_f %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff)
up_entrez <- unique(na.omit(mapIds(
  org.Dm.eg.db, keys = sig_up$gene, column = "ENTREZID",
  keytype = "SYMBOL", multiVals = "first"
)))
ora <- enrichGO(
  gene = up_entrez, OrgDb = org.Dm.eg.db, keyType = "ENTREZID",
  ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE
)
write.csv(as.data.frame(ora), out("GO_female_up_excluding_CSF3.csv"), row.names = FALSE)

# A hard up-DEG cutoff has no significant GO over-representation after CSF3 is
# removed, so modules are selected from GSEA across every tested gene.
rank_df <- df_f %>% filter(is.finite(stat)) %>% select(gene, stat)
rank_df$ENTREZID <- mapIds(
  org.Dm.eg.db, keys = rank_df$gene, column = "ENTREZID",
  keytype = "SYMBOL", multiVals = "first"
)
rank_df <- rank_df %>%
  filter(!is.na(ENTREZID)) %>%
  arrange(desc(abs(stat))) %>%
  distinct(ENTREZID, .keep_all = TRUE)
gene_list <- rank_df$stat
names(gene_list) <- rank_df$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

set.seed(FIG_SEED)
gsea <- gseGO(
  geneList = gene_list, OrgDb = org.Dm.eg.db, keyType = "ENTREZID",
  ont = "ALL", minGSSize = 10, maxGSSize = 500,
  pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE
)
gsea_df <- as.data.frame(gsea)
write.csv(gsea_df, out("GSEA_GO_female_CSvsWT_excluding_CSF3.csv"), row.names = FALSE)

# The most significant GO results are highly redundant. These representatives
# consolidate them into three distinct, interpretable programs.
module_defs <- tribble(
  ~description, ~go_id, ~direction,
  "mitochondrial translation/ribosome", "GO:0032543", "down in CS",
  "innate immune/antimicrobial defense", "GO:0042742", "down in CS",
  "oxidoreductase activity", "GO:0016491", "up in CS"
)

missing_terms <- setdiff(module_defs$go_id, gsea_df$ID)
if (length(missing_terms) > 0) {
  stop("Required representative GSEA terms are missing: ", paste(missing_terms, collapse = ", "))
}

module_summary <- module_defs %>%
  left_join(
    gsea_df %>% select(ID, GO_description = Description, setSize, NES, p.adjust, core_enrichment),
    by = c("go_id" = "ID")
  )
if (any(module_summary$p.adjust >= 0.05)) {
  stop("At least one selected module is not significant at BH-adjusted p < 0.05.")
}

# Map leading-edge Entrez IDs to symbols and retain the strongest 15 per module.
module_memberships <- module_summary %>%
  select(description, go_id, direction, NES, p.adjust, core_enrichment) %>%
  separate_rows(core_enrichment, sep = "/") %>%
  rename(ENTREZID = core_enrichment)
module_memberships$gene <- mapIds(
  org.Dm.eg.db, keys = module_memberships$ENTREZID, column = "SYMBOL",
  keytype = "ENTREZID", multiVals = "first"
)
module_memberships <- module_memberships %>%
  filter(!is.na(gene)) %>%
  left_join(df_f %>% select(gene, log2FoldChange, padj, stat), by = "gene") %>%
  filter(is.finite(stat)) %>%
  mutate(direction_score = sign(NES) * stat) %>%
  group_by(description) %>%
  arrange(desc(direction_score), .by_group = TRUE) %>%
  slice_head(n = genes_per_module) %>%
  ungroup()
if (drop_cg) {
  module_memberships <- module_memberships %>% filter(!grepl(cg_pattern, gene))
}

gene_modules <- module_memberships %>%
  group_by(gene) %>%
  summarise(
    module_memberships = paste(sort(unique(description)), collapse = "; "),
    membership_count = n_distinct(description),
    log2FoldChange = first(log2FoldChange),
    padj = first(padj),
    stat = first(stat),
    .groups = "drop"
  ) %>%
  mutate(description = case_when(
    grepl("mitochondrial translation", module_memberships, fixed = TRUE) ~
      "mitochondrial translation/ribosome",
    grepl("innate immune", module_memberships, fixed = TRUE) ~
      "innate immune/antimicrobial defense",
    TRUE ~ "oxidoreductase activity"
  ))

module_summary$n_leading_edge_selected <- vapply(
  module_summary$description,
  function(x) sum(module_memberships$description == x),
  integer(1)
)
module_summary$selection <- paste0(
  "Up to ", genes_per_module,
  " leading-edge genes ranked in the direction of the GSEA effect"
)
write.csv(
  module_summary %>% select(-core_enrichment),
  out("GSEA_top_biological_modules_excluding_CSF3.csv"), row.names = FALSE
)

# -----------------------------------------------------------------------------
# STRING mapping and interaction retrieval.
# -----------------------------------------------------------------------------
pm <- rba_string_map_ids(ids = gene_modules$gene, species = string_species)
if (!is.data.frame(pm) || nrow(pm) == 0) stop("STRING ID mapping returned no results.")
pm$gene <- if ("queryItem" %in% names(pm)) pm$queryItem else pm$preferredName
pm <- pm %>% select(stringId, gene) %>% distinct(gene, .keep_all = TRUE)

dc2 <- gene_modules %>% inner_join(pm, by = "gene") %>% distinct(stringId, .keep_all = TRUE)
if (nrow(dc2) < 3) stop("Too few selected genes mapped to STRING.")
int_net <- rba_string_interactions_network(
  ids = dc2$stringId, species = string_species, required_score = string_score
)
edges <- data.frame(from = int_net$stringId_A, to = int_net$stringId_B) %>%
  filter(from %in% dc2$stringId, to %in% dc2$stringId) %>%
  distinct()

module_levels <- module_defs$description
dc2$description <- factor(dc2$description, levels = module_levels)
dc2 <- dc2 %>% arrange(description, gene)
module_cols <- c(
  "mitochondrial translation/ribosome" = "#9B8AC4",
  "innate immune/antimicrobial defense" = "#E5A24A",
  "oxidoreductase activity" = "#70B99A"
)
down_col <- unname(FIG_COLORS["control"])
up_col <- unname(FIG_COLORS["mutant"])
fc_limit <- max(2, ceiling(max(abs(dc2$log2FoldChange), na.rm = TRUE)))

caption_text <- paste0(
  "Female Drosophila bulk RNA-seq: WT n = 3, CS n = 2 biological libraries; ",
  "CSF3 excluded before all analysis steps.\n",
  "Nodes: up to ", genes_per_module,
  " GSEA leading-edge genes per consolidated module; BH-adjusted GSEA p < 0.05."
)

# -----------------------------------------------------------------------------
# Ring view.
# -----------------------------------------------------------------------------
dc2$idx <- seq_len(nrow(dc2)) - 1
dc2$xnode <- seq(0, 2 * pi, length.out = nrow(dc2) + 1)[-1]
node_r <- 0.55
label_gap <- 0.06
label_r <- node_r + label_gap
gene_lab_size <- 2.5
per_char_reach <- 0.04
band_gap <- 0.001
band_width <- 0.06
band_label_gap <- 0.10
end_gap_physical <- 0.013
text_reach <- per_char_reach * max(nchar(dc2$gene))
group_band_inner <- label_r + text_reach + band_gap
group_band_outer <- group_band_inner + band_width
group_label_r <- group_band_outer + band_label_gap

dc2 <- dc2 %>% mutate(
  y = node_r * cos(xnode), x = node_r * sin(xnode),
  ang_deg = atan2(cos(xnode), sin(xnode)) * 180 / pi,
  flip = ang_deg > 90 | ang_deg < -90,
  lab_ang = ifelse(flip, ang_deg + 180, ang_deg),
  lhjust = ifelse(flip, 1, 0),
  label_x = label_r * sin(xnode), label_y = label_r * cos(xnode)
)
links <- edges %>%
  left_join(dc2 %>% select(stringId, x, y), by = c("from" = "stringId")) %>%
  rename(xfrom = x, yfrom = y) %>%
  left_join(dc2 %>% select(stringId, x, y), by = c("to" = "stringId")) %>%
  rename(xto = x, yto = y) %>%
  filter(if_all(c(xfrom, yfrom, xto, yto), is.finite))

node_step <- 2 * pi / nrow(dc2)
band_mid_r <- (group_band_inner + group_band_outer) / 2
arc_gap <- min(end_gap_physical / band_mid_r, node_step * 0.48)
rect0 <- dc2 %>%
  group_by(description) %>%
  summarise(
    lo = min(xnode) - node_step / 2 + arc_gap,
    hi = max(xnode) + node_step / 2 - arc_gap,
    xmid = (min(xnode) + max(xnode)) / 2,
    .groups = "drop"
  )
rect <- rect0 %>% rowwise() %>% do({
  r <- .
  if (r$hi > 2 * pi) {
    tibble(description = r$description,
           xmin = c(r$lo, 0), xmax = c(2 * pi, r$hi - 2 * pi), xmid = r$xmid)
  } else if (r$lo < 0) {
    tibble(description = r$description,
           xmin = c(0, r$lo + 2 * pi), xmax = c(r$hi, 2 * pi), xmid = r$xmid)
  } else {
    tibble(description = r$description, xmin = r$lo, xmax = r$hi, xmid = r$xmid)
  }
}) %>% ungroup()
rect_lab <- rect0 %>% select(description, xmid)

p1_limit <- label_r + text_reach + 0.10
p1 <- ggplot(dc2, aes(x, y)) +
  geom_segment(
    data = links,
    aes(x = xfrom, xend = xto, y = yfrom, yend = yto),
    inherit.aes = FALSE, linewidth = 0.4, color = "gray65", alpha = 0.28
  ) +
  geom_point(aes(fill = log2FoldChange), shape = 21, size = 5, color = "grey25") +
  geom_point(
    data = subset(dc2, membership_count > 1),
    shape = 21, fill = NA, color = "black", size = 6, stroke = 0.8
  ) +
  scale_fill_gradient2(
    low = down_col, mid = "white", high = up_col, midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish,
    name = expression(log[2] ~ FC ~ (female~CS/WT)),
    guide = guide_colorbar(
      barwidth = 7, barheight = 0.9, direction = "horizontal",
      title.position = "top", title.hjust = 0.5
    )
  ) +
  geom_text(
    aes(label = gene, x = label_x, y = label_y, angle = lab_ang, hjust = lhjust),
    size = gene_lab_size, fontface = "italic", family = FIG_FONT
  ) +
  coord_fixed() + theme_void() +
  theme(legend.position = c(0.5, 0.5)) +
  scale_y_continuous(limits = c(-p1_limit, p1_limit), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-p1_limit, p1_limit), expand = c(0, 0))

p2_limit <- group_label_r + 0.06
p2 <- ggplot() +
  geom_rect(
    data = rect,
    aes(xmin = xmin, xmax = xmax, ymin = group_band_inner,
        ymax = group_band_outer, fill = description)
  ) +
  geom_textpath(
    data = rect_lab,
    aes(x = xmid, y = group_label_r, label = description, color = description),
    linetype = 0, size = 4.0, fontface = "bold", upright = TRUE,
    family = FIG_FONT
  ) +
  scale_fill_manual(values = module_cols, guide = "none") +
  scale_color_manual(values = module_cols, guide = "none") +
  coord_polar() + theme_void() +
  scale_y_continuous(limits = c(-p2_limit, p2_limit), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0))

ring <- ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(p2, x = 0.04, y = 0.08, width = 0.92, height = 0.82) +
  draw_plot(p1, x = 0.04, y = 0.08, width = 0.92, height = 0.82) +
  draw_label(
    "Top female biological programs\nafter excluding CSF3",
    x = 0.04, y = 0.97, hjust = 0, vjust = 1,
    size = FIG_TITLE_SIZE, fontface = "bold", color = "black", lineheight = 0.95,
    fontfamily = FIG_FONT
  ) +
  draw_label(
    caption_text, x = 0.04, y = 0.015, hjust = 0, vjust = 0,
    size = FIG_CAPTION_SIZE, color = "grey35", fontfamily = FIG_FONT
  )
fig_save(
  ring, out("fig2_string_ring_female_excluding_CSF3"),
  width = 7.8, height = 7.8
)

# -----------------------------------------------------------------------------
# Clustered STRING network view.
# -----------------------------------------------------------------------------
nodes2 <- dc2 %>%
  transmute(
    name = stringId, gene, description, module_memberships,
    membership_count, log2FoldChange, padj, stat
  )
gt <- as_tbl_graph(graph_from_data_frame(edges, vertices = nodes2, directed = FALSE)) %>%
  mutate(description = factor(description, levels = module_levels))
n_cl <- length(module_levels)
cl_pos <- tibble(
  description = factor(module_levels, levels = module_levels),
  ca = seq(0, 2 * pi, length.out = n_cl + 1)[-1],
  cx = cos(ca) * 10, cy = sin(ca) * 10
)
lay <- create_layout(gt, layout = "circle") %>% left_join(cl_pos, by = "description")
for (grp in module_levels) {
  i <- which(lay$description == grp)
  ang <- seq(0, 2 * pi, length.out = length(i) + 1)[-1]
  lay$x[i] <- lay$cx[i] + cos(ang) * 2.7
  lay$y[i] <- lay$cy[i] + sin(ang) * 2.7
}
ell <- lay %>%
  group_by(description) %>%
  summarise(
    cx = mean(x), cy = mean(y),
    w = max(x) - min(x) + 2.6,
    h = max(y) - min(y) + 2.6,
    .groups = "drop"
  )
lab <- ell %>% mutate(
  lx = cx,
  above = cy >= 0,
  ly = ifelse(above, cy + h / 2 + 1.0, cy - h / 2 - 1.0),
  vjust = ifelse(above, 0, 1)
)

p_network <- ggraph(lay) +
  geom_ellipse(
    data = ell,
    aes(x0 = cx, y0 = cy, a = w / 2, b = h / 2, fill = description, angle = 0),
    alpha = 0.20, color = NA
  ) +
  scale_fill_manual(values = module_cols, guide = "none") +
  geom_edge_link(color = "gray65", alpha = 0.24, linewidth = 0.45) +
  geom_node_point(aes(color = log2FoldChange), size = 5) +
  geom_node_point(
    aes(filter = membership_count > 1), shape = 21,
    fill = NA, color = "black", size = 6, stroke = 0.8
  ) +
  scale_color_gradient2(
    low = down_col, mid = "white", high = up_col, midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish,
    name = expression(log[2] ~ FC ~ (female~CS/WT))
  ) +
  geom_node_text(
    aes(label = gene), repel = TRUE, max.overlaps = Inf,
    size = FIG_ANNOT_SIZE - 0.9, fontface = "italic",
    family = FIG_FONT, color = "black"
  ) +
  geom_text(
    data = lab,
    aes(x = lx, y = ly, label = description, vjust = vjust),
    fontface = "bold", size = FIG_ANNOT_SIZE,
    family = FIG_FONT, hjust = 0.5
  ) +
  coord_fixed() +
  expand_limits(x = c(-15.5, 15.5), y = c(-15.5, 15.5)) +
  labs(
    title = "Top female biological programs after excluding CSF3",
    subtitle = paste0(
      "STRING interactions among GSEA leading-edge genes; interaction score >= ",
      string_score, "\nNode color shows female CS-versus-WT log2 fold change"
    ),
    caption = caption_text
  ) +
  theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(
      size = FIG_TITLE_SIZE, face = "bold", hjust = 0, color = "black",
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE, hjust = 0, color = "grey35",
      margin = margin(b = 6)
    ),
    plot.caption = element_text(
      size = FIG_CAPTION_SIZE, hjust = 0, color = "grey35",
      margin = margin(t = 6)
    ),
    legend.title = element_text(size = FIG_LEGEND_TEXT_SIZE),
    legend.text = element_text(size = FIG_LEGEND_TEXT_SIZE),
    plot.margin = margin(
      t = FIG_MARGIN_PT[["t"]], r = FIG_MARGIN_PT[["r"]],
      b = FIG_MARGIN_PT[["b"]], l = FIG_MARGIN_PT[["l"]]
    )
  )
fig_save(
  p_network, out("fig2_string_network_female_excluding_CSF3"),
  width = 8.8, height = 8.8
)

audit <- dc2 %>%
  select(
    gene, stringId, description, module_memberships, membership_count,
    log2FoldChange, padj, stat
  ) %>%
  mutate(
    comparison = "female dMIC60-CS versus female dMIC60-WT; CSF3 excluded",
    selection = paste0(
      "Up to ", genes_per_module,
      " GSEA leading-edge genes per consolidated significant module"
    )
  )
write.csv(audit, out("fig2_string_network_female_excluding_CSF3_audit.csv"), row.names = FALSE)

cat("CSF3 excluded before all analysis steps.\n")
cat("Female biological libraries: WT n = 3; CS n = 2.\n")
cat("Female up-DEGs (padj < 0.05, log2FC > 1): ", nrow(sig_up), "\n", sep = "")
cat("Significant up-DEG GO over-representation terms: ", nrow(as.data.frame(ora)), "\n", sep = "")
cat("Selected significant GSEA modules:\n")
print(module_summary %>% select(description, GO_description, NES, p.adjust, n_leading_edge_selected))
cat("STRING nodes: ", nrow(dc2), "; STRING edges: ", nrow(edges), "\n", sep = "")
