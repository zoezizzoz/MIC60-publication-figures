#!/usr/bin/env Rscript
# =============================================================================
# dMIC60 oxidation-resistance bulk RNA-seq analysis
#
# Adapted from G. Barisanos human proteomics STRING/GO pipeline
# ("Chandra_proteomics.Rmd") to Drosophila melanogaster bulk RNA-seq.
#
# KEY DIFFERENCES FROM THE ORIGINAL PROTEOMICS CODE
#   1. Input is a RAW COUNT MATRIX (genes x samples), not a table of
#      precomputed log2FC / p-values. We therefore add a DESeq2 front-end
#      to compute differential expression before any GO/STRING step.
#   2. Species switched from human -> Drosophila:
#        org.Hs.eg.db  -> org.Dm.eg.db
#        STRING species 9606 -> 7227
#   3. Comparison of interest: dMIC60-CS (oxidation-resistant) vs dMIC60-WT,
#      analysed SEPARATELY per sex. Females = main figures, males = supplement.
#
# DESIGN (12 libraries, 3 replicates x 4 groups):
#   CS = dMIC60-CS (Cys->Ser, oxidation-resistant)   WR = dMIC60-WT
#   F  = female                                        M  = male
#   Sample naming convention: <genotype><sex><replicate>, e.g. CSF1, WRM3
#   -> CSF1-3, CSM1-3, WRF1-3, WRM1-3
#
# SECTIONS
#   1. Load counts, run DESeq2, write result tables
#   2. GO enrichment on the female up-DEGs (table only; feeds section 4)
#   3. Individual gene plots            -> fig2b_genes
#   4. STRING network ring + ellipse    -> fig2_string_ring_female
#   5. Stress-pathway specificity       -> figS3_stress_distinct
#
# OUTPUTS: DESeq2 tables, the GO table, the stress-set audit tables (CSV) and
#          the three figures above, all written under OUTPUT_DIR.
# All figure styling comes from figure_style.R (see "shared aesthetics").
# =============================================================================

required_packages <- c(
  "DESeq2", "tidyverse", "clusterProfiler", "org.Dm.eg.db",
  "AnnotationDbi", "ggrepel", "rbioapi", "igraph", "tidygraph",
  "ggraph", "ggforce", "geomtextpath", "cowplot", "scales", "readxl"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
suppressMessages({
  library(DESeq2)
  library(tidyverse)
  library(clusterProfiler)
  library(org.Dm.eg.db)      # Drosophila annotation (was org.Hs.eg.db)
  library(AnnotationDbi)     # mapIds
  library(ggrepel)
  library(rbioapi)           # STRING REST API
  library(igraph); library(tidygraph); library(ggraph)
  library(ggforce); library(geomtextpath); library(cowplot)
  library(scales)            # squish
  library(readxl)            # stress-gene workbook
})

# =============================================================================
# CONFIG — the only things that should need editing to rerun this analysis
# =============================================================================
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
FIGURE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
PACKAGE_DIR <- normalizePath(file.path(FIGURE_DIR, ".."), mustWork = TRUE)
COUNT_MATRIX_PATH <- file.path(FIGURE_DIR, "Original_Data", "count_matrix_symbol.csv")
STRESS_GENES_PATH <- file.path(FIGURE_DIR, "Original_Data", "RNA-seq Stress Genes.xlsx")
COMMON_STRESS_GENES_PATH <- file.path(FIGURE_DIR, "Original_Data", "common stress genes.xlsx")
OUTPUT_DIR        <- file.path(FIGURE_DIR, "Rebuilt_Output")
STRING_SPECIES    <- 7227                             # NCBI taxon ID, D. melanogaster
PADJ_CUTOFF       <- 0.05                             # significance threshold
LFC_CUTOFF        <- 1                                # |log2FC| threshold for "DE"
SEED              <- 1

# Uncharacterized fly loci with CG#### symbols carry no functional annotation.
# When TRUE they are excluded from BOTH figures that would otherwise name them:
#   section 4  removed as ring/ellipse nodes, so they are not drawn or labelled
#   section 5  removed from the stress sets AND from the Fisher universe, since
#              dropping them from a curated set while leaving thousands of them
#              in the background would compare a characterized set against a
#              partly uncharacterized one
# Set FALSE to restore the previous behaviour. Widen the pattern to
# "^C[GR][0-9]+$" to drop non-coding CR#### loci as well.
DROP_CG_GENES   <- TRUE
CG_GENE_PATTERN <- "^CG[0-9]+$"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)
out <- function(...) file.path(OUTPUT_DIR, ...)   # shorthand for output paths

# ---- shared aesthetics -------------------------------------------------------
# Colors, fonts, text sizes, dots, error bars and the "bracket only when
# significant" rule all come from figure_style.R, the same file used by the
# TEM, TMRM and MTT figures. Keep a copy in this folder.
style_file <- file.path(SCRIPT_DIR, "figure_style.R")
source(style_file)

pubtheme <- theme_fig(legend_position = "right")   # project theme
geno_labs <- c(WR = "dMIC60-WT", CS = "dMIC60-CS")
up_col   <- unname(FIG_COLORS["mutant"])            # up in CS
down_col <- unname(FIG_COLORS["control"])           # down in CS

# =============================================================================
# 1. LOAD COUNTS & BUILD DESeq2 DATASET
# =============================================================================
# Named count_mat, not counts: `counts` would mask DESeq2::counts(), which is
# called two lines below.
count_mat <- read.csv(COUNT_MATRIX_PATH, row.names = 1, check.names = FALSE)
count_mat <- as.matrix(count_mat); storage.mode(count_mat) <- "integer"

samples <- colnames(count_mat)
geno <- ifelse(grepl("^CS", samples), "CS", "WR")
sex  <- ifelse(substr(samples, 3, 3) == "F", "F", "M")
coldata <- data.frame(
  row.names = samples,
  genotype  = factor(geno, levels = c("WR", "CS")),   # WR = reference
  sex       = factor(sex,  levels = c("F", "M")),
  group     = factor(paste0(geno, sex)))

rna_group_n <- table(coldata$group)
n_group <- function(x) unname(as.integer(rna_group_n[x]))
rna_all_caption <- paste0(
  "Drosophila bulk RNA-seq; biological libraries: female WT n = ",
  n_group("WRF"), ", female CS n = ", n_group("CSF"),
  ", male WT n = ", n_group("WRM"), ", male CS n = ", n_group("CSM"), "."
)
rna_female_caption <- paste0(
  "Female Drosophila bulk RNA-seq; n = ", n_group("WRF"),
  " WT and ", n_group("CSF"), " CS biological libraries."
)

dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ group)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]          # low-count filter
dds <- DESeq(dds)

# CS-vs-WT contrast within each sex
res_female <- results(dds, contrast = c("group", "CSF", "WRF"))  # + = up in CS female
res_male   <- results(dds, contrast = c("group", "CSM", "WRM"))
norm_counts <- counts(dds, normalized = TRUE)   # `norm` would mask base::norm()
vsd <- vst(dds, blind = TRUE)                   # kept only for the checkpoint

fmt <- function(res) { d <- as.data.frame(res); d$gene <- rownames(d); d[order(d$padj), ] }
df_f <- fmt(res_female); df_m <- fmt(res_male)
write.csv(df_f, out("deseq2_female_CSvsWR.csv"), row.names = FALSE)
write.csv(df_m, out("deseq2_male_CSvsWR.csv"),   row.names = FALSE)
write.csv(norm_counts, out("normalized_counts.csv"))   # size-factor-normalized
# RDS field kept as `norm` so anything that already reads this file still works
saveRDS(list(dds = dds, vsd = vsd, norm = norm_counts, coldata = coldata,
             df_f = df_f, df_m = df_m), out("deseq2_objects.rds"))   # checkpoint

# ---- sample-labelled PCA QC -------------------------------------------------
# DESeq2::plotPCA() uses the 500 most variable genes by default. Individual
# labels make it clear whether an apparent group shift is replicated or driven
# by one library; CSF3 is outlined because it sits between the female and male
# clusters and carries the male-biased transcript signal audited for Fig 2C.
pca_data <- plotPCA(vsd, intgroup = c("genotype", "sex"), returnData = TRUE)
pca_var <- round(100 * attr(pca_data, "percentVar"), 1)
pca_data$Sample <- rownames(pca_data)
pca_data$Genotype <- factor(geno_labs[as.character(pca_data$genotype)],
                            levels = c("dMIC60-WT", "dMIC60-CS"))
pca_data$Sex <- factor(ifelse(pca_data$sex == "F", "Female", "Male"),
                       levels = c("Female", "Male"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Genotype, shape = Sex)) +
  geom_point(size = FIG_PT_SIZE + 0.8, alpha = 1) +
  geom_point(data = subset(pca_data, Sample == "CSF3"),
             mapping = aes(PC1, PC2),
             shape = 21, fill = NA, color = "black",
             size = FIG_PT_SIZE + 2.0, stroke = 0.7,
             inherit.aes = FALSE) +
  ggrepel::geom_text_repel(
    aes(label = Sample, fontface = ifelse(Sample == "CSF3", "bold", "plain")),
    color = "black", size = FIG_ANNOT_SIZE - 0.5,
    box.padding = 0.45, point.padding = 0.30,
    min.segment.length = 0, segment.size = 0.22,
    seed = FIG_SEED, max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(values = fig_colors(c("dMIC60-WT", "dMIC60-CS"))) +
  scale_shape_manual(
    values = c(Female = 16, Male = 17),
    labels = c(Female = "\u2640 Female", Male = "\u2642 Male")
  ) +
  labs(
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)"),
    title = "RNA-seq samples separate primarily by sex",
    subtitle = "CSF3 shifts away from the other female libraries toward the male cluster",
    caption = paste0(
      rna_all_caption,
      "\nDESeq2 variance-stabilized counts; 500 most variable genes. CSF3 is outlined."
    )
  ) +
  theme_fig(legend_position = "right")

fig_save(p_pca, out("fig1b_qc_pca"), width = FIG_W_WIDE + 1.0, height = 4.8)

# =============================================================================
# 2. GO ENRICHMENT ON THE FEMALE UP-DEGs
#    Section 4 labels each up-regulated gene with its top GO term, so the
#    enrichment TABLE is needed even though the GO dot-plot figure is not part
#    of this script. This block is what defined ego_f in the full pipeline.
# =============================================================================
run_go <- function(d, direction = "up", lfc = LFC_CUTOFF, padj = PADJ_CUTOFF) {
  sig <- if (direction == "up") {
    d[!is.na(d$padj) & d$padj < padj & d$log2FoldChange >  lfc, ]
  } else {
    d[!is.na(d$padj) & d$padj < padj & d$log2FoldChange < -lfc, ]
  }
  ent <- unique(na.omit(mapIds(org.Dm.eg.db, keys = sig$gene, column = "ENTREZID",
                               keytype = "SYMBOL", multiVals = "first")))
  enrichGO(gene = ent, OrgDb = org.Dm.eg.db, keyType = "ENTREZID",
           ont = "ALL", pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
}

ego_f <- run_go(df_f, "up")
if (is.null(ego_f) || nrow(as.data.frame(ego_f)) == 0) {
  stop("No GO terms are enriched among the female up-DEGs, so section 4 has ",
       "nothing to group the ring by. Loosen PADJ_CUTOFF / LFC_CUTOFF, or ",
       "skip sections 2 and 4.")
}
write.csv(as.data.frame(ego_f), out("GO_female_up.csv"), row.names = FALSE)

# =============================================================================
# 3. INDIVIDUAL GENE PLOTS  (Fig 2B)
# =============================================================================
genes_show <- c("timeout", "ImpL2", "DNAlig3")

# A gene dropped by the low-count filter would fail with an opaque subscript
# error below, so check membership up front.
absent_genes <- setdiff(genes_show, rownames(norm_counts))
if (length(absent_genes) > 0) {
  warning("Not in the filtered count matrix, dropped from Fig 2B: ",
          paste(absent_genes, collapse = ", "))
  genes_show <- intersect(genes_show, rownames(norm_counts))
}
if (length(genes_show) == 0) {
  stop("None of the Fig 2B genes survived the low-count filter.")
}
long <- do.call(rbind, lapply(genes_show, function(g)
  data.frame(gene = g, count = norm_counts[g, ],
             genotype = coldata$genotype, sex = coldata$sex)))
long$Genotype <- factor(geno_labs[as.character(long$genotype)],
                        levels = c("dMIC60-WT", "dMIC60-CS"))
long$Sex  <- ifelse(long$sex == "F", "Female", "Male")
long$gene <- factor(long$gene, levels = genes_show)

# Bracket and adjusted p, placed with the same offsets used by every other
# figure and kept only where padj < FIG_ALPHA.
sig_ann <- do.call(rbind, lapply(genes_show, function(g)
  do.call(rbind, lapply(c("Female", "Male"), function(s) {
    res  <- if (s == "Female") df_f else df_m
    vals <- long$count[long$gene == g & long$Sex == s]
    span <- max(diff(range(vals)), abs(max(vals)) * 0.15, 1e-9)
    data.frame(gene = g, Sex = s, padj = res[g, "padj"],
               y = max(vals) + span * FIG_BRACKET_OFFSET,
               tick = span * FIG_BRACKET_TICK, stringsAsFactors = FALSE)
  }))))
sig_ann <- sig_ann[is.finite(sig_ann$padj) & sig_ann$padj < FIG_ALPHA, , drop = FALSE]
sig_ann$label <- vapply(sig_ann$padj, fig_p_label, character(1))
sig_ann$gene  <- factor(sig_ann$gene, levels = genes_show)

# fig_sig_brackets_at() draws with annotate(), which carries no facet columns
# and would therefore repeat the same bracket in every panel. A faceted bracket
# has to be mapped from data, so it is built here from sig_ann; the geometry
# constants are still the shared ones.
facet_brackets <- function(d, x1 = 1, x2 = 2) {
  if (nrow(d) == 0) return(NULL)
  d$x1 <- x1; d$x2 <- x2; d$xmid <- (x1 + x2) / 2
  seg <- function(mapping) {
    geom_segment(data = d, mapping = mapping, inherit.aes = FALSE,
                 linewidth = FIG_BRACKET_LINEWIDTH, color = "black")
  }
  list(
    seg(aes(x = x1, xend = x2, y = y, yend = y)),
    seg(aes(x = x1, xend = x1, y = y, yend = y - tick)),
    seg(aes(x = x2, xend = x2, y = y, yend = y - tick)),
    geom_text(data = d, aes(x = xmid, y = y + tick * 0.6, label = label),
              inherit.aes = FALSE, vjust = 0, size = FIG_ANNOT_SIZE,
              lineheight = 0.95, family = FIG_FONT, color = "black")
  )
}

set.seed(FIG_SEED)
p_genes <- ggplot(long, aes(Genotype, count, fill = Genotype, color = Genotype)) +
  fig_summary() +
  fig_points() +
  facet_brackets(sig_ann) +
  scale_fill_manual(values = fig_fills(c("dMIC60-WT", "dMIC60-CS")), guide = "none") +
  scale_color_manual(values = fig_colors(c("dMIC60-WT", "dMIC60-CS")), guide = "none") +
  # Extra upper expansion keeps the two-line significance labels clear of the
  # new facet dividers; the base expansion still comes from figure_style.R.
  scale_y_continuous(expand = expansion(
    mult = c(FIG_Y_EXPAND[1], FIG_Y_EXPAND[2] + 0.06)
  )) +
  facet_grid(gene ~ Sex, scales = "free_y", switch = "y") +
  labs(x = NULL, y = "Normalized counts", title = "Individual gene expression",
       caption = paste0(
         "Drosophila bulk RNA-seq; n = 3 biological libraries per genotype and sex.\n",
         "DESeq2 adjusted p, shown only where padj < ", FIG_ALPHA, "."
       )) +
  theme_fig() +
  theme(
        # Outlined facets visually separate each gene, following the reference
        # layout while using the project's standard line weight and type scale.
        panel.border = element_rect(fill = NA, color = FIG_FACET_BORDER_COLOR,
                                    linewidth = FIG_LINE_WIDTH),
        panel.spacing.y = unit(FIG_FACET_GAP_PT, "pt"),
        strip.background = element_rect(fill = FIG_FACET_STRIP_FILL,
                                        color = FIG_FACET_BORDER_COLOR,
                                        linewidth = FIG_LINE_WIDTH),
        strip.text.y.left = element_text(angle = 0, face = "italic",
                                         size = FIG_AXIS_TEXT_SIZE),
        axis.text.x = element_text(angle = 20, hjust = 1))
fig_save(p_genes, out("fig2b_genes"), width = FIG_W_WIDE, height = 6.4)

# =============================================================================
# 4. STRING NETWORK "RING" PLOT  (Fig 2C)
#    Direct adaptation of the proteomics "DEP GRAPH" block, species = STRING_SPECIES.
#    (a) assign each female up-DEG to its top GO term
#    (b) map genes -> STRING, pull interactions among them
#    (c) circular node layout coloured by log2FC + geom_textpath GO-term arcs
# =============================================================================

# The publication-cleanup rebuild only needs the Fig S3 stress panel. STRING
# requires a live remote service, so skip the unrelated network panel here;
# the archived network outputs and dedicated Figure 2D script remain intact.
if (FALSE) {
sigf_up <- df_f[!is.na(df_f$padj) & df_f$padj < PADJ_CUTOFF & df_f$log2FoldChange > LFC_CUTOFF, ]

# Biological modules are defined from ALL relevant enriched GO memberships.
# A gene is no longer assigned to only its single most significant GO term.
# Closely related/redundant terms are consolidated into two interpretable
# modules, and every qualifying membership is retained in `module_memberships`.
module_terms <- list(
  "mitochondrial respiration and energy metabolism" = c(
    "generation of precursor metabolites and energy",
    "energy derivation by oxidation of organic compounds",
    "aerobic respiration", "cellular respiration",
    "tricarboxylic acid cycle", "ATP metabolic process",
    "mitochondrial membrane", "mitochondrial envelope",
    "mitochondrial inner membrane",
    "inner mitochondrial membrane protein complex"
  ),
  "cilium and axoneme assembly and motility" = c(
    "cilium assembly", "cilium organization", "axoneme assembly",
    "microtubule bundle formation", "cilium movement",
    "cilium or flagellum-dependent cell motility",
    "cilium-dependent cell motility", "microtubule-based movement",
    "cilium", "axoneme", "motile cilium"
  ),
  "phosphoregulation and protein dephosphorylation" = c(
    "protein dephosphorylation", "phosphatase activity",
    "protein serine/threonine phosphatase activity",
    "phosphoric ester hydrolase activity", "dephosphorylation",
    "phosphoprotein phosphatase activity"
  )
)
module_term_map <- stack(module_terms) %>%
  transmute(Description = as.character(values),
            module = as.character(ind))

dep_classes <- as.data.frame(ego_f) %>%
  inner_join(module_term_map, by = "Description") %>%
  separate_rows(geneID, sep = "/") %>%
  filter(geneID %in% sigf_up$gene) %>%
  distinct(geneID, module)

gene_modules <- dep_classes %>%
  group_by(gene = geneID) %>%
  summarise(
    module_memberships = paste(sort(unique(module)), collapse = "; "),
    membership_count = n_distinct(module),
    .groups = "drop"
  ) %>%
  # A node needs one physical position, so shared genes are placed in the most
  # specific functional module rather than by GO p-value. Their complete
  # membership is retained above and marked by an outline in the network.
  mutate(description = case_when(
    grepl("phosphoregulation and protein dephosphorylation", module_memberships,
          fixed = TRUE) ~ "phosphoregulation/protein dephosphorylation",
    grepl("mitochondrial respiration and energy metabolism", module_memberships,
          fixed = TRUE) ~ "mitochondrial energy metabolism",
    TRUE ~ "cilium/axoneme assembly and motility"
  ))

dc <- gene_modules %>%
  left_join(sigf_up[, c("gene", "log2FoldChange", "padj")], by = "gene")

# Uncharacterized CG#### loci add unreadable labels to the ring. Removing them
# HERE, before the STRING query, also removes them as network nodes and as edge
# endpoints, so the radial layout, the interaction query and the exported CSV
# all describe the same gene set.
if (DROP_CG_GENES) {
  ring_cg <- sort(unique(dc$gene[grepl(CG_GENE_PATTERN, dc$gene)]))
  dc <- dc %>% filter(!grepl(CG_GENE_PATTERN, gene))
  cat("Ring plot: ", length(ring_cg), " uncharacterized CG loci removed",
      if (length(ring_cg) > 0) paste0(" (", paste(ring_cg, collapse = ", "), ")") else "",
      "\n", sep = "")
}

# Everything below divides by nrow(dc2), so an empty selection has to stop here
# rather than fail later with a division-by-zero layout.
if (nrow(dc) == 0) {
  stop("No female up-DEG belongs to the configured biological modules. ",
       "Check GO_female_up.csv and update module_terms.")
}

pm <- rba_string_map_ids(ids = dc$gene, species = STRING_SPECIES)
string_available <- is.data.frame(pm) && nrow(pm) > 0
if (string_available) {
  pm$gene <- pm$preferredName
  int_net <- rba_string_interactions_network(
    ids = pm$stringId, species = STRING_SPECIES, required_score = 400
  )
} else {
  cached_map_path <- out("fig2_string_network_female_audit.csv")
  if (!file.exists(cached_map_path)) {
    stop("STRING is unavailable and no cached STRING mapping was found.")
  }
  cached_map <- read.csv(cached_map_path, check.names = FALSE)
  if (!"stringId" %in% names(cached_map) && "name" %in% names(cached_map)) {
    names(cached_map)[names(cached_map) == "name"] <- "stringId"
  }
  pm <- unique(cached_map[, c("stringId", "gene")])
  int_net <- data.frame(stringId_A = character(), stringId_B = character())
  message("STRING is unavailable; using the cached gene-to-STRING mapping for ",
          "the ring and preserving the existing interaction-network files.")
}

dc2 <- merge(dc, pm[, c("stringId", "gene")], by = "gene") %>% distinct(stringId, .keep_all = TRUE)
if (nrow(dc2) == 0) {
  stop("No gene in dc mapped to a STRING id; check the STRING_SPECIES setting.")
}
dc2 <- dc2[order(dc2$description, dc2$gene), ]; dc2$idx <- 0:(nrow(dc2) - 1)
dc2$xnode <- seq(0, 2 * pi, length.out = nrow(dc2) + 1)[-1]

# ---- RADIAL LAYOUT ----------------------------------------------------------
# The whole ring is driven off these few constants rather than magic numbers
# scattered through the plotting code below, so it can be re-tuned from one
# place. Quick reference if spacing looks off after a render:
#   dots too close to / far from labels        -> node_r / label_gap
#   band overlaps the ends of long gene names   -> raise per_char_reach or band_gap
#   band leaves a visible gap after gene names   -> lower per_char_reach or band_gap
#   band-end bars touch / merge into each other  -> raise end_gap_physical
#   band-end bars leave too much dead space      -> lower end_gap_physical
#   group-name letters (g/y/p) touch the band     -> raise band_label_gap
#   dot ring and band do not read as concentric    -> adjust `scale =` in
#                                                     draw_plot(p1, ...) below

node_r           <- 0.55   # gene-dot ring radius
label_gap        <- 0.06   # dot -> start of its name
gene_lab_size    <- 2.5    # gene-name font size
per_char_reach   <- 0.04   # radial length reserved per character of a gene name
band_gap         <- 0.001  # end of the longest gene name -> inner edge of the group band
band_width       <- 0.06   # thickness of the coloured group band
band_label_gap   <- 0.1    # band -> its curved group-name label
end_gap_physical <- 0.013  # physical space between adjacent group-band ends

label_r <- node_r + label_gap
text_reach <- per_char_reach * max(nchar(as.character(dc2$gene)))
group_band_inner <- label_r + text_reach + band_gap
group_band_outer <- group_band_inner + band_width
group_label_r    <- group_band_outer + band_label_gap

dc2 <- dc2 %>% mutate(
  y = node_r * cos(xnode), x = node_r * sin(xnode),
  # gene labels rotated RADIALLY (Circos-style): each name points outward along
  # its nodes angle, so adjacent labels fan apart and never overlap even where
  # nodes are dense. Text is flipped 180 deg on the left half to stay upright.
  ang_deg = atan2(cos(xnode), sin(xnode)) * 180 / pi,
  flip    = ang_deg > 90 | ang_deg < -90,
  lab_ang = ifelse(flip, ang_deg + 180, ang_deg),
  lhjust  = ifelse(flip, 1, 0),
  label_x = label_r * sin(xnode), label_y = label_r * cos(xnode))
links <- data.frame(from = int_net$stringId_A, to = int_net$stringId_B) %>%
  left_join(dc2, by = c("from" = "stringId")) %>% rename(xfrom = x, yfrom = y) %>%
  left_join(dc2, by = c("to"   = "stringId")) %>% rename(xto = x, yto = y) %>%
  select(xfrom, xto, yfrom, yto)
# Extend each cluster arc by half a node-spacing beyond its first/last dot so the
# band covers ALL its nodes, and subtract a small gap so adjacent sectors
# (e.g. phosphatase / cilium assembly) read as separate arcs, not one ring.
# The gap is specified as a physical (radial-unit) distance and converted to an
# angle via /band radius -- an angular gap alone would shrink visually as the
# band moves closer to the center, since arc length = radius * angle.
node_step <- 2 * pi / nrow(dc2)
band_mid_r <- (group_band_inner + group_band_outer) / 2
arc_gap   <- end_gap_physical / band_mid_r
arc_gap   <- min(arc_gap, node_step * 0.48)  # preserve each neighboring node slice
rect0 <- dc2 %>% group_by(description) %>%
  summarize(lo = min(xnode) - node_step/2 + arc_gap,
            hi = max(xnode) + node_step/2 - arc_gap,
            xmid = (min(xnode) + max(xnode)) / 2, .groups = "drop")
# One cluster straddles the 0/2*pi seam; split its band into two rectangles so
# the polar x-axis can stay locked at [0, 2*pi] (no coord_polar stretching that
# would make the wrap-point sectors touch) while every arc still covers its dots.
rect <- rect0 %>% rowwise() %>% do({
  r <- .
  if (r$hi > 2 * pi)      tibble(description = r$description, xmin = c(r$lo, 0),          xmax = c(2 * pi, r$hi - 2 * pi), xmid = r$xmid)
  else if (r$lo < 0)      tibble(description = r$description, xmin = c(0, r$lo + 2 * pi), xmax = c(r$hi, 2 * pi),           xmid = r$xmid)
  else                    tibble(description = r$description, xmin = r$lo,                xmax = r$hi,                     xmid = r$xmid)
}) %>% ungroup()
rect_lab <- rect0 %>% transmute(description, xmid)   # one curved label per cluster
clust_cols <- c(
  "mitochondrial energy metabolism" = "#d98bb5",
  "cilium/axoneme assembly and motility" = "#9b8ac4",
  "phosphoregulation/protein dephosphorylation" = "#7bc4e2"
)

p1_limit <- label_r + text_reach + 0.1     # just past the label layer, no leftover margin
p1 <- ggplot(dc2, aes(x, y)) +
  geom_segment(data = links, aes(x = xfrom, xend = xto, y = yfrom, yend = yto),
               linewidth = 0.4, color = "gray70", alpha = 0.25, inherit.aes = FALSE) +
  geom_point(shape = 21, aes(fill = log2FoldChange), size = 5, color = "grey30") +
  scale_fill_gradient2(low = down_col, mid = "white", high = up_col, midpoint = 0,
                       limits = c(-2, 6), oob = squish, breaks = c(0, 3, 6), labels = c("0", "3", "6+"),
                       name = expression(log[2] ~ FC ~ (CS/WT)),
                       guide = guide_colorbar(barwidth = 7, barheight = 0.9, direction = "horizontal",
                                              title.position = "top", title.hjust = 0.5)) +
  geom_text(aes(label = gene, x = label_x, y = label_y, angle = lab_ang, hjust = lhjust),
            size = gene_lab_size, fontface = "italic") +
  coord_fixed() + theme_void() + theme(legend.position = c(0.5, 0.5)) +
  # Tighter limits crop unused space around the node/gene-label layer.
  scale_y_continuous(limits = c(-p1_limit, p1_limit), expand = c(0, 0)) +
  scale_x_continuous(limits = c(-p1_limit, p1_limit), expand = c(0, 0))
p2_limit <- group_label_r + 0.06           # just past the outer group-name label
p2 <- ggplot() +
  geom_rect(data = rect,
            aes(xmin = xmin, xmax = xmax,
                ymin = group_band_inner, ymax = group_band_outer,
                fill = description)) +
  geom_textpath(data = rect_lab,
                aes(x = xmid, y = group_label_r,
                    label = description, color = description),
                linetype = 0, size = 4.1, fontface = "bold", upright = TRUE) +
  scale_fill_manual(values = clust_cols, guide = "none") +
  scale_color_manual(values = clust_cols, guide = "none") +
  coord_polar() + theme_void() +
  # Lock the polar axis at exactly [0, 2*pi] so wrap-point sectors keep their gap
  # (the seam-crossing cluster is drawn as two rects, above, so nothing is clipped).
  scale_y_continuous(limits = c(-p2_limit, p2_limit), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0))
ring <- ggdraw(xlim = c(0.035, 0.965), ylim = c(0.035, 0.965)) +
  draw_plot(p2) +
  draw_plot(p1, scale = 1) +
  draw_label(rna_female_caption, x = 0.04, y = 0.015,
             hjust = 0, vjust = 0, size = FIG_CAPTION_SIZE,
             color = "grey35", fontfamily = FIG_FONT)
if (string_available) {
  fig_save(ring, out("fig2_string_ring_female"), width = 7.6, height = 7.6)
  write.csv(dc2[, c("gene", "stringId", "description", "module_memberships",
                    "membership_count", "log2FoldChange", "padj")],
            out("ring_female_network.csv"), row.names = FALSE)
}


# --- alternative "ellipse-cluster" network view (the other STRING style) ------
if (string_available) {
nodes2 <- dc2[, c("stringId", "description", "module_memberships",
                  "membership_count", "log2FoldChange", "padj", "gene")]
colnames(nodes2)[1] <- "name"
# STRING can rename a symbol (preferredName != our input symbol); such a node is
# dropped from dc2 by the merge above but its interactions remain in int_net,
# leaving edges that point at absent vertices. graph_from_data_frame() errors on
# that ("Some vertex names in d are not listed in vertices"), so keep only edges
# whose BOTH endpoints survive in nodes2.
edges2 <- data.frame(from = int_net$stringId_A, to = int_net$stringId_B)
edges2 <- edges2[edges2$from %in% nodes2$name & edges2$to %in% nodes2$name, ]
gt <- as_tbl_graph(graph_from_data_frame(edges2, vertices = nodes2, directed = FALSE)) %>%
  mutate(description = as.factor(description))
n_cl <- length(unique(nodes2$description))
cl_pos <- tibble(description = unique(nodes2$description),
                 ca = seq(0, 2 * pi, length.out = n_cl + 1)[-1], cx = cos(ca) * 10, cy = sin(ca) * 10)
lay <- create_layout(gt, layout = "circle") %>% left_join(cl_pos, by = "description")
for (grp in unique(lay$description)) {
  i <- which(lay$description == grp); ang <- seq(0, 2 * pi, length.out = length(i) + 1)[-1]
  lay$x[i] <- lay$cx[i] + cos(ang) * 2.6; lay$y[i] <- lay$cy[i] + sin(ang) * 2.6
}
ell <- lay %>% group_by(description) %>%
  summarise(cx = mean(x), cy = mean(y), w = max(x) - min(x) + 2.5, h = max(y) - min(y) + 2.5, .groups = "drop")
# cluster labels sit horizontally ABOVE (upper circles) or BELOW (lower circles);
# phosphatase is forced above regardless of centre so its name clears the ring.
lab <- ell %>% mutate(
  above = cy >= 0 | description == "phosphatase",
  lx = cx,
  ly = ifelse(above, cy + h/2 + 1.1, cy - h/2 - 1.1),
  vjust = ifelse(above, 0, 1))
p_ell <- ggraph(lay) +
  geom_ellipse(data = ell, aes(x0 = cx, y0 = cy, a = w/2, b = h/2, fill = description, angle = 0),
               alpha = 0.18, color = NA) +
  scale_fill_manual(values = clust_cols, guide = "none") +
  geom_edge_link(color = "gray70", alpha = 0.18) +
  geom_node_point(aes(color = log2FoldChange), size = 5) +
  geom_node_point(aes(filter = membership_count > 1), shape = 21,
                  fill = NA, color = "black", size = 6, stroke = 0.7) +
  # Every displayed node is up-regulated, so use the mutant-color gradient
  # instead of a diverging blue/red scale that could imply down-DEGs are shown.
  scale_color_gradient(low = fig_lighten(up_col, 0.78), high = up_col,
                       limits = c(LFC_CUTOFF, 6), oob = squish,
                       breaks = c(1, 3, 6), labels = c("1", "3", "6+"),
                       name = expression(log[2] ~ FC ~ (female~CS/WT))) +
  geom_node_text(aes(label = gene), repel = TRUE,
                 size = FIG_ANNOT_SIZE - 0.8, fontface = "italic",
                 family = FIG_FONT, color = "black") +
  geom_text(data = lab, aes(x = lx, y = ly, label = gsub("\n", " ", description), vjust = vjust),
            fontface = "bold", size = FIG_ANNOT_SIZE, family = FIG_FONT,
            hjust = 0.5) +
  coord_fixed() + expand_limits(x = c(-15.5, 15.5), y = c(-15, 15)) +
  labs(
    title = "Female dMIC60-CS vs female dMIC60-WT: up-regulated DEGs",
    subtitle = paste0(
      "Nodes: padj < ", PADJ_CUTOFF, " and log2FC > ", LFC_CUTOFF,
      "; shaded regions: consolidated enriched GO modules\n",
      "Edges: STRING interactions; outlined nodes: membership in multiple modules"
    ),
    caption = paste0(
      nrow(nodes2), " up-regulated female DEGs shown across ", n_cl,
      " consolidated biological modules.\n", rna_female_caption
    )
  ) +
  theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0,
                              color = "black", margin = margin(b = 3)),
    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0,
                                 color = "grey35", margin = margin(b = 6)),
    plot.caption = element_text(size = FIG_CAPTION_SIZE, hjust = 0,
                                color = "grey35", margin = margin(t = 6)),
    legend.title = element_text(size = FIG_LEGEND_TEXT_SIZE, color = "black"),
    legend.text = element_text(size = FIG_LEGEND_TEXT_SIZE, color = "black"),
    plot.margin = margin(t = FIG_MARGIN_PT[["t"]], r = FIG_MARGIN_PT[["r"]],
                         b = FIG_MARGIN_PT[["b"]], l = FIG_MARGIN_PT[["l"]])
  )
fig_save(p_ell, out("fig2_string_network_female"), width = 8.6, height = 8.6)

# Machine-readable audit: the comparison, direction, thresholds, and GO
# category for every node that appears in the female STRING network.
network_audit <- dc2[, c("gene", "stringId", "description",
                         "module_memberships", "membership_count",
                         "log2FoldChange", "padj")]
network_audit$comparison <- "female dMIC60-CS vs female dMIC60-WT"
network_audit$direction <- "up-regulated in female dMIC60-CS"
network_audit$selection <- paste0("padj < ", PADJ_CUTOFF,
                                  "; log2FoldChange > ", LFC_CUTOFF)
write.csv(network_audit, out("fig2_string_network_female_audit.csv"),
          row.names = FALSE)
}

# End of the optional STRING network section.
}

# =============================================================================
# 5. STRESS-PATHWAY SPECIFICITY  (Fig S3)
#    Is the female program distinct from canonical stress responses, or a
#    flavour of one?
#
#    GENE SETS COME FROM THE WORKBOOK AT STRESS_GENES_PATH (see CONFIG).
#    That workbook is a DIOPT (flyrnai.org/diopt) ortholog table: a curated
#    list of human stress-pathway genes, each mapped to every predicted
#    Drosophila ortholog, with a "Stress Pathway" label per input gene
#    (Control, ISR, HSR, OSR, UPR (ATF6), UPR (IRE1/XBP1s)). This replaces the
#    earlier GO-membership definitions and the hand-curated mito-UPR effector
#    list, so the sets are now traceable to one citable source file.
#
#    THREE FILTERS ARE APPLIED, EACH REPORTED TO THE CONSOLE:
#      1. Ortholog confidence. DIOPT returns every prediction, including weak
#         ones (DNAJC10 alone has 21 fly hits, e.g. Nedd4 for BAG3). Only
#         DIOPT_KEEP_RANKS are kept; set DIOPT_REQUIRE_BEST_SCORE to also
#         demand the best-scoring ortholog for that human gene.
#      2. Detection. A set member must appear in the DESeq2-tested background,
#         which is the same universe the Fisher test uses.
#      3. Symbol case. DIOPT and the count matrix occasionally disagree on
#         capitalisation only; those are rescued and listed, not silently
#         dropped or silently matched.
#
#    Every set on this figure comes from that workbook and nothing else: there
#    is no GO-derived column. The oxidative-metabolism contrast now lives in
#    section 2 (GO_female_up.csv) and in the ring plot's energy-metabolism
#    cluster, not here.
#
#    A one-sided Fisher test then asks whether each set is enriched among the
#    female up-DEGs relative to the other testable genes in the RNA-seq Stress
#    Genes workbook. A canonical stress response would be enriched; a distinct
#    program leaves the canonical sets at background. The workbook's own
#    Control set is a built-in negative control for that comparison.
#
#    CAVEAT: this count matrix under-captures some highly expressed
#    housekeeping isoforms (eEF1alpha1, Hsc70-5 (mortalin), Lon, Hsf < 30 mean
#    counts), a symbol-mapping artefact, so a gene missing from a set is not
#    evidence that it is unexpressed. The argument rests on the many detected
#    effectors, not on any single master regulator.
# =============================================================================

# ---- inputs and filter settings ---------------------------------------------
STRESS_GENES_SHEET <- "RNA-seq Stress Genes"   # path is set in CONFIG, above
COMMON_STRESS_GENES_SHEET <- "Common Stress Genes"

DIOPT_KEEP_RANKS         <- c("high", "moderate")  # or just "high"
DIOPT_REQUIRE_BEST_SCORE <- FALSE
STRESS_CASE_RESCUE       <- TRUE

# Left-to-right order in the figure. Any pathway present in the workbook but
# missing from this vector is appended alphabetically rather than dropped.
STRESS_SET_ORDER <- c("Control", "ISR", "UPR (ATF6)", "UPR (IRE1/XBP1s)",
                      "HSR", "OSR")

# Recognizable representatives of the workbook's Control set to identify in
# the figure even when they are neither common-stress genes nor DEGs.
CONTROL_LABEL_GENES <- c("betaGlu", "mRpL9", "mtSSB", "Prosbeta6")

# ---- read and validate the workbook -----------------------------------------
if (!file.exists(STRESS_GENES_PATH)) {
  message("Stress-gene workbook not found at: ", STRESS_GENES_PATH,
          "\nSelect it manually.")
  STRESS_GENES_PATH <- file.choose()
}

stress_raw <- as.data.frame(
  read_excel(STRESS_GENES_PATH, sheet = STRESS_GENES_SHEET),
  check.names = FALSE)
names(stress_raw) <- trimws(names(stress_raw))

stress_needed <- c("Stress Pathway", "Input Gene Symbol", "Output Gene Symbol",
                   "DIOPT Score", "Rank", "Best Score")
stress_absent <- setdiff(stress_needed, names(stress_raw))
if (length(stress_absent) > 0) {
  stop(
    "Missing required columns in sheet ",
    STRESS_GENES_SHEET,
    ": ",
    paste(stress_absent, collapse = ", "),
    "\nColumns found: ",
    paste(names(stress_raw), collapse = ", ")
  )
}
stress_map <- data.frame(
  pathway = trimws(as.character(stress_raw[["Stress Pathway"]])),
  human   = trimws(as.character(stress_raw[["Input Gene Symbol"]])),
  fly     = trimws(as.character(stress_raw[["Output Gene Symbol"]])),
  diopt   = suppressWarnings(as.numeric(stress_raw[["DIOPT Score"]])),
  rank    = tolower(trimws(as.character(stress_raw[["Rank"]]))),
  best    = trimws(as.character(stress_raw[["Best Score"]])),
  stringsAsFactors = FALSE)

# Human genes with no fly ortholog at all: reported, then dropped.
stress_no_ortholog <- unique(
  stress_map[is.na(stress_map$fly) | stress_map$fly %in% c("", "NA"),
             c("pathway", "human")])
stress_map <- stress_map[!is.na(stress_map$fly) &
                           !stress_map$fly %in% c("", "NA"), ]

# ---- filter 1: ortholog confidence ------------------------------------------
stress_map$pass_rank <- stress_map$rank %in% DIOPT_KEEP_RANKS
if (DIOPT_REQUIRE_BEST_SCORE) {
  stress_map$pass_rank <- stress_map$pass_rank &
    toupper(stress_map$best) == "YES"
}

# ---- filter 2 and 3: detection, with a case-only rescue ---------------------
# TWO INITIAL GENE LISTS, DELIBERATELY DIFFERENT:
#   tested_genes   every gene in the filtered count matrix. Decides SET
#                  MEMBERSHIP, so a set member that DESeq2 could not test stays
#                  visible on the figure as "padj filtered (NA)".
#   all_testable_genes genes that could actually have been called DE in this
#                  contrast, i.e. non-NA padj. After the stress sets are built,
#                  this is intersected with their union to make the Fig S3
#                  Fisher universe contain only genes from RNA-seq Stress Genes.
#                  DESeq2 assigns padj = NA to genes removed by independent
#                  filtering and to Cook's-distance outliers; those genes can
#                  never enter up_all, so counting them in the background
#                  deflates the background DEG rate and inflates every odds
#                  ratio. They are excluded from the test only.
tested_genes   <- rownames(dds)
all_testable_genes <- df_f$gene[!is.na(df_f$padj)]

if (DROP_CG_GENES) {
  all_testable_genes <- all_testable_genes[
    !grepl(CG_GENE_PATTERN, all_testable_genes)]
}
tested_lower <- tested_genes[!duplicated(tolower(tested_genes))]
names(tested_lower) <- tolower(tested_lower)

stress_map$detected_symbol <- NA_character_
exact <- stress_map$fly %in% tested_genes
stress_map$detected_symbol[exact] <- stress_map$fly[exact]

if (STRESS_CASE_RESCUE) {
  rescue <- !exact & tolower(stress_map$fly) %in% names(tested_lower)
  stress_map$detected_symbol[rescue] <-
    unname(tested_lower[tolower(stress_map$fly[rescue])])
  # report every case-only match, not just the ones that clear the rank
  # filter: a mismatch here points at symbol handling in the count matrix
  stress_rescued <- unique(stress_map[rescue, c("fly", "detected_symbol")])
} else {
  stress_rescued <- stress_map[0, c("fly", "detected_symbol")]
}

stress_map$detected <- !is.na(stress_map$detected_symbol)
stress_map$is_cg <- grepl(CG_GENE_PATTERN, stress_map$detected_symbol)
stress_map$in_set <- stress_map$pass_rank & stress_map$detected &
  !(DROP_CG_GENES & stress_map$is_cg)

# ---- assemble the sets -------------------------------------------------------
present_sets <- unique(stress_map$pathway[stress_map$pass_rank])
set_order <- c(STRESS_SET_ORDER[STRESS_SET_ORDER %in% present_sets],
               sort(setdiff(present_sets, STRESS_SET_ORDER)))

stress_sets <- lapply(set_order, function(pw)
  sort(unique(stress_map$detected_symbol[stress_map$pathway == pw &
                                           stress_map$in_set])))
names(stress_sets) <- set_order

empty_sets <- names(stress_sets)[lengths(stress_sets) == 0]
if (length(empty_sets) > 0) {
  warning("No detected genes in: ", paste(empty_sets, collapse = ", "),
          ". Dropping from the figure.")
  stress_sets <- stress_sets[lengths(stress_sets) > 0]
  set_order <- set_order[set_order %in% names(stress_sets)]
}
if (length(stress_sets) < 2) {
  stop("Fewer than two usable gene sets. Check the workbook and the filters.")
}

# Fig S3 analyzes only the genes that are both present in the count-matrix
# analysis and included in RNA-seq Stress Genes. The competitive Fisher tests
# therefore compare each pathway against the other curated stress genes, not
# against the full transcriptome.
stress_universe <- sort(unique(unlist(stress_sets, use.names = FALSE)))
testable_genes <- intersect(all_testable_genes, stress_universe)
if (length(testable_genes) == 0) {
  stop("No RNA-seq Stress Genes have a non-NA female adjusted p value.")
}

# ---- report what the filters did --------------------------------------------
stress_set_summary <- do.call(rbind, lapply(set_order, function(pw) {
  rows <- stress_map[stress_map$pathway == pw, ]
  data.frame(
    pathway = pw,
    # counts the human genes that had at least one fly ortholog; genes with
    # none are listed separately below
    human_genes_mapped = length(unique(rows$human)),
    predicted_orthologs = length(unique(rows$fly)),
    passed_rank_filter = length(unique(rows$fly[rows$pass_rank])),
    cg_dropped = if (DROP_CG_GENES) {
      length(unique(rows$detected_symbol[rows$pass_rank & rows$detected & rows$is_cg]))
    } else 0L,
    detected_in_data = length(stress_sets[[pw]]),
    stringsAsFactors = FALSE)
}))

cat("\nStress-pathway gene sets from ", basename(STRESS_GENES_PATH), "\n",
    "Ortholog ranks kept: ", paste(DIOPT_KEEP_RANKS, collapse = ", "),
    if (DIOPT_REQUIRE_BEST_SCORE) "; best-scoring ortholog required" else "",
    "\n", sep = "")
print(stress_set_summary, row.names = FALSE)

if (nrow(stress_no_ortholog) > 0) {
  cat("\nHuman genes with no fly ortholog in the workbook:\n")
  print(stress_no_ortholog, row.names = FALSE)
}

if (nrow(stress_rescued) > 0) {
  cat("\nSymbols matched on case only (workbook -> count matrix):\n")
  print(stress_rescued, row.names = FALSE)
}

dropped <- unique(stress_map[stress_map$pass_rank & !stress_map$detected,
                             c("pathway", "human", "fly")])
cat("\nOrthologs passing the rank filter but absent from the count matrix: ",
    nrow(dropped), "\n", sep = "")

if (DROP_CG_GENES) {
  cg_dropped <- sort(unique(stress_map$detected_symbol[
    stress_map$pass_rank & stress_map$detected & stress_map$is_cg]))
  cat("Uncharacterized CG-numbered loci removed from the sets: ",
      length(cg_dropped),
      if (length(cg_dropped) > 0) paste0(" (", paste(cg_dropped, collapse = ", "), ")") else "",
      "\n", sep = "")
}

untestable <- setdiff(unlist(stress_sets), testable_genes)
cat("Set members plotted but not testable (padj = NA in this contrast): ",
    length(untestable), "\n", sep = "")
cat("Fig S3 Fisher universe (testable genes from RNA-seq Stress Genes): ",
    length(testable_genes), "\n", sep = "")

shared <- table(unlist(lapply(stress_sets, unique)))
shared <- names(shared)[shared > 1]
if (length(shared) > 0) {
  cat("Genes appearing in more than one set (tested independently in each): ",
      paste(shared, collapse = ", "), "\n", sep = "")
}

write.csv(stress_map[, c("pathway", "human", "fly", "diopt", "rank", "best",
                         "detected_symbol", "pass_rank", "detected", "is_cg",
                         "in_set")],
          out("figS3_stress_gene_mapping.csv"), row.names = FALSE)
write.csv(stress_set_summary, out("figS3_stress_set_summary.csv"),
          row.names = FALSE)

# ---- per-gene effect sizes ---------------------------------------------------
stress_long <- do.call(rbind, lapply(names(stress_sets), function(pw) {
  genes <- stress_sets[[pw]]
  data.frame(pathway = pw, gene = genes,
             lfc_F = df_f[genes, "log2FoldChange"], padj_F = df_f[genes, "padj"],
             lfc_M = df_m[genes, "log2FoldChange"], padj_M = df_m[genes, "padj"],
             stringsAsFactors = FALSE)
}))
write.csv(stress_long, out("figS3_stress_pathways.csv"), row.names = FALSE)

# ---- one-sided Fisher enrichment among female up-DEGs -----------------------
up_all <- df_f$gene[!is.na(df_f$padj) & df_f$padj < PADJ_CUTOFF &
                      df_f$log2FoldChange > LFC_CUTOFF]
up_all <- intersect(up_all, testable_genes)   # same universe as the test
# n_detected is what the figure plots; n_tested is what the p value describes.
fisher_path <- function(genes) {
  n_detected <- length(genes)
  genes <- intersect(genes, testable_genes)
  
  if (length(genes) == 0) {
    return(data.frame(n_detected = n_detected, n_tested = 0L, n_up = 0L,
                      pct_up = NA_real_, odds_ratio = NA_real_,
                      fisher_p = NA_real_))
  }
  
  in_up  <- sum(genes %in% up_all); in_not  <- length(genes) - in_up
  out_genes <- setdiff(testable_genes, genes)
  out_up <- sum(out_genes %in% up_all); out_not <- length(out_genes) - out_up
  ft <- fisher.test(matrix(c(in_up, in_not, out_up, out_not), 2),
                    alternative = "greater")
  data.frame(n_detected = n_detected, n_tested = length(genes), n_up = in_up,
             pct_up = round(100 * in_up / length(genes), 1),
             odds_ratio = round(unname(ft$estimate), 2),
             fisher_p = signif(ft$p.value, 3))
}
enr <- do.call(rbind, lapply(names(stress_sets), function(pw)
  cbind(pathway = pw, fisher_path(stress_sets[[pw]]))))

# Several competitive tests on overlapping sets: BH is valid under that kind of
# positive dependence, and the adjusted value is what the figure reports.
enr$fisher_padj <- signif(p.adjust(enr$fisher_p, method = "BH"), 3)
write.csv(enr, out("figS3_enrichment.csv"), row.names = FALSE)

cat("\nEnrichment among female up-DEGs (universe: ", length(testable_genes),
    " testable genes, ", length(up_all), " up-DEGs):\n", sep = "")
print(enr, row.names = FALSE)

# ---- figure -----------------------------------------------------------------
# Axis labels break "UPR (ATF6)" and "metabolic program" onto two lines
# without putting newlines into the set names themselves, which the join
# between the effect sizes and the enrichment table relies on.
wrap_set <- function(x) {
  x <- sub("\\s+\\((.*)\\)$", "\n(\\1)", x)
  sub(" (program)$", "\n\\1", x)
}
set_labels <- setNames(wrap_set(names(stress_sets)), names(stress_sets))

stress_plot_dat <- stress_long
stress_plot_dat$sig <- ifelse(
  !is.na(stress_plot_dat$padj_F) & stress_plot_dat$padj_F < PADJ_CUTOFF &
    stress_plot_dat$lfc_F > LFC_CUTOFF, "up (padj<0.05, LFC>1)",
  ifelse(!is.na(stress_plot_dat$padj_F) & stress_plot_dat$padj_F < PADJ_CUTOFF &
           stress_plot_dat$lfc_F < -LFC_CUTOFF, "down",
         ifelse(is.na(stress_plot_dat$padj_F), "padj filtered (NA)", "n.s.")))
stress_plot_dat$pathway <- factor(stress_plot_dat$pathway,
                                  levels = names(stress_sets))

# Label a stress gene when it is EITHER listed in the Common Stress Genes sheet
# OR significantly up- or down-regulated in the female CS-vs-WT comparison.
if (!file.exists(COMMON_STRESS_GENES_PATH)) {
  stop("Common-stress-gene workbook not found at: ", COMMON_STRESS_GENES_PATH)
}
common_stress_raw <- as.data.frame(
  read_excel(COMMON_STRESS_GENES_PATH,
             sheet = COMMON_STRESS_GENES_SHEET, skip = 4),
  check.names = FALSE)
names(common_stress_raw) <- trimws(names(common_stress_raw))

common_needed <- c("Recommended fly ortholog(s)", "Other best-score matches")
common_absent <- setdiff(common_needed, names(common_stress_raw))
if (length(common_absent) > 0) {
  stop(
    "Missing required columns in sheet ", COMMON_STRESS_GENES_SHEET, ": ",
    paste(common_absent, collapse = ", "),
    "\nColumns found: ", paste(names(common_stress_raw), collapse = ", ")
  )
}

split_fly_symbols <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}
common_label_genes <- sort(unique(c(
  split_fly_symbols(common_stress_raw[["Recommended fly ortholog(s)"]]),
  split_fly_symbols(common_stress_raw[["Other best-score matches"]])
)))
common_label_genes <- common_label_genes[
  nzchar(common_label_genes) & common_label_genes != "No ortholog reported"]
common_label_lower <- tolower(common_label_genes)
control_label_lower <- tolower(CONTROL_LABEL_GENES)

stress_plot_dat$lab <- ifelse(
  tolower(stress_plot_dat$gene) %in% common_label_lower |
    stress_plot_dat$sig %in% c("up (padj<0.05, LFC>1)", "down") |
    (stress_plot_dat$pathway == "Control" &
       tolower(stress_plot_dat$gene) %in% control_label_lower),
  stress_plot_dat$gene, NA_character_)

cat("Common Stress Genes, significant up/down genes, OR selected Control genes labelled: ",
    length(unique(stats::na.omit(stress_plot_dat$lab))),
    if (any(!is.na(stress_plot_dat$lab))) {
      paste0(" (", paste(sort(unique(stats::na.omit(stress_plot_dat$lab))),
                          collapse = ", "), ")")
    } else "",
    "\n", sep = "")

smeds <- do.call(rbind, lapply(levels(stress_plot_dat$pathway), function(pw)
  data.frame(pathway = pw,
             m = median(stress_plot_dat$lfc_F[stress_plot_dat$pathway == pw],
                        na.rm = TRUE))))
smeds$pathway <- factor(smeds$pathway, levels = names(stress_sets))

# Annotation height follows the data instead of a hard-coded y
lfc_range <- range(stress_plot_dat$lfc_F, na.rm = TRUE)
lfc_span <- max(diff(lfc_range), 1)
y_annot <- lfc_range[2] + lfc_span * 0.10

enr2 <- enr
enr2$pathway <- factor(enr2$pathway, levels = names(stress_sets))
enr2$lab <- ifelse(
  is.na(enr2$fisher_padj), "not testable",
  ifelse(enr2$fisher_padj < 0.001,
         sprintf("OR %.1f\np_adj=%.0e", enr2$odds_ratio, enr2$fisher_padj),
         sprintf("OR %.2f\np_adj=%.2f", enr2$odds_ratio, enr2$fisher_padj)))
enr2$y <- y_annot

bg_rate <- round(100 * length(up_all) / length(testable_genes), 1)

# Title states the result rather than asserting a fixed conclusion, so it
# cannot go stale if the sets or the filters change.
enriched <- as.character(enr$pathway[!is.na(enr$fisher_padj) &
                                       enr$fisher_padj < PADJ_CUTOFF])
s3_title <- if (length(enriched) == 0) {
  "Female CS-vs-WT: no stress-pathway set is enriched among up-DEGs"
} else if (length(enriched) <= 3) {
  paste0("Female CS-vs-WT: ", paste(enriched, collapse = ", "),
         if (length(enriched) == 1) " is" else " are",
         " enriched among up-DEGs")
} else {
  paste0("Female CS-vs-WT: ", length(enriched), " of ", length(stress_sets),
         " sets are enriched among up-DEGs (", paste(enriched, collapse = ", "), ")")
}
s3_title <- paste(strwrap(s3_title, width = 95), collapse = "\n")

s3_subtitle <- paste0(
  "Sets: ", basename(STRESS_GENES_PATH), " (DIOPT ",
  paste(DIOPT_KEEP_RANKS, collapse = "/"),
  "-confidence orthologs detected in the count-matrix analysis)\n",
  "Labels: genes in ", basename(COMMON_STRESS_GENES_PATH),
  " OR significant up/down genes; selected Control genes also named\n",
  "Fisher background: ", length(testable_genes),
  " testable RNA-seq stress genes; ", bg_rate,
  "% are up-DEGs (non-NA padj",
  if (DROP_CG_GENES) ", CG-numbered loci excluded" else "",
  "); BH-adjusted across ", length(stress_sets), " sets"
)

# Generate the horizontal jitter once and reuse the resulting coordinate in the
# point and label layers. If each layer jitters independently, a connector aims
# at the category centre instead of the displayed dot.
set.seed(FIG_SEED)
stress_plot_dat$x_plot <- as.numeric(stress_plot_dat$pathway) +
  runif(nrow(stress_plot_dat),
        min = -(FIG_JITTER_WIDTH + 0.04),
        max =  (FIG_JITTER_WIDTH + 0.04))
smeds$x_plot <- match(as.character(smeds$pathway), names(stress_sets))
enr2$x_plot <- match(as.character(enr2$pathway), names(stress_sets))

p_s3 <- ggplot(stress_plot_dat, aes(x_plot, lfc_F)) +
  annotate("rect", xmin = 0.4, xmax = length(stress_sets) + 0.6,
           ymin = -LFC_CUTOFF, ymax = LFC_CUTOFF, fill = "grey90", alpha = 0.6) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = 2,
             color = "grey60", linewidth = 0.3) +
  geom_point(aes(color = sig), size = FIG_PT_SIZE, alpha = 0.85) +
  geom_crossbar(data = smeds, aes(x = x_plot, y = m, ymin = m, ymax = m),
                width = FIG_SUMMARY_WIDTH, color = FIG_BOX_OUTLINE,
                linewidth = FIG_LINE_WIDTH, inherit.aes = FALSE) +
  # A lightly padded white label keeps names readable over the median bar and
  # nearby points. Straight, uncropped segments connect the label edge to the
  # corresponding dot; stronger repulsion prevents the dense near-zero labels
  # from stacking on one another.
  geom_label_repel(aes(label = lab, color = sig),
                   size = FIG_ANNOT_SIZE - 0.7, fontface = "italic",
                   fill = scales::alpha("white", 0.92), label.size = NA,
                   label.padding = grid::unit(0.10, "lines"),
                   label.r = grid::unit(0.05, "lines"),
                   max.overlaps = Inf, max.time = 5,
                   min.segment.length = 0, segment.size = 0.22,
                   segment.color = "grey55", segment.curvature = 0,
                   box.padding = 0.55, point.padding = 0.22,
                   force = 5, force_pull = 1.2, seed = SEED,
                   na.rm = TRUE, show.legend = FALSE) +
  geom_text(data = enr2, aes(x = x_plot, y = y, label = lab),
            size = FIG_ANNOT_SIZE - 0.3, lineheight = 0.9,
            fontface = ifelse(!is.na(enr2$fisher_padj) &
                                enr2$fisher_padj < 0.001, "bold", "plain"),
            color = ifelse(!is.na(enr2$fisher_padj) &
                             enr2$fisher_padj < PADJ_CUTOFF, up_col, "grey35"),
            inherit.aes = FALSE) +
  scale_color_manual(values = c("up (padj<0.05, LFC>1)" = up_col,
                                "down" = down_col,
                                "n.s." = "grey65",
                                "padj filtered (NA)" = "grey82"), name = NULL,
                     breaks = c("up (padj<0.05, LFC>1)", "down", "n.s.",
                                "padj filtered (NA)")) +
  scale_x_continuous(breaks = seq_along(stress_sets),
                     labels = unname(set_labels),
                     limits = c(0.4, length(stress_sets) + 0.6),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
  labs(x = NULL, y = expression(log[2] ~ fold ~ change ~ (CS/WT)),
       title = s3_title,
       subtitle = s3_subtitle,
       caption = rna_female_caption) +
  pubtheme + theme(legend.position = "top",
                   axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE - 1,
                                              lineheight = 0.85)) +
  guides(color = guide_legend(nrow = 1,
                              override.aes = list(size = 2.6, label = "")))

fig_save(p_s3, out("figS3_stress_distinct"),
         width = max(FIG_W_2COL, 1.3 * length(stress_sets) + 1.5), height = 6.2)

# =============================================================================
# SESSION INFO — captured for reproducibility / methods reporting
# =============================================================================
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))

# NOTE: to run males through the STRING ring plot, repeat sections 2 and 4 with
#       ego_m <- run_go(df_m, "up") and df_m in place of df_f. Males lack the
#       coordinated oxidative-metabolism program, so expect a sparse male ring.
cat("Pipeline complete. Outputs written to:", normalizePath(OUTPUT_DIR), "\n")
