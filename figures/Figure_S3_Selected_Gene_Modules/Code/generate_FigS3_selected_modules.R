# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
})
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[1])))
panel_dir <- dirname(script_dir)
checkpoint_file <- file.path(panel_dir, "Supporting_Data", "deseq2_objects.rds")
output_dir <- file.path(panel_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint <- readRDS(checkpoint_file)
df_f <- checkpoint$df_f
coldata <- checkpoint$coldata
required_samples <- c("CSF1", "CSF2", "CSF3", "WRF1", "WRF2", "WRF3")
stopifnot(all(required_samples %in% rownames(coldata)))
stopifnot(sum(coldata$group == "CSF") == 3, sum(coldata$group == "WRF") == 3)
padj_cutoff <- 0.05
lfc_cutoff <- 0.58
module_genes <- list(FOXO = c("ImpL2", "Lsp1gamma", "Gadd45", "Lsp1alpha", "GstE1", "Atg6", "Ilp6", "puc", "Rbf", "Cat",
    "dap", "Sesn", "bmm", "Atg8a"), `DNA replication` = c("timeout", "DNAlig3", "PolD2", "Top1", "DNAlig1", "PolE2", "PCNA2",
    "Cdc45", "Mcm6", "Cdc6", "Mcm2", "PolE1", "DNAlig4", "Prim2", "Mcm5", "RfC4", "Orc4", "PolD1", "RnrL", "Orc1", "RPA1",
    "Mcm7", "Orc5", "PolA1", "RnrS", "Orc6", "RfC3", "Mcm10", "Orc2", "Prim1", "PolA2"), Spargel = c("COX4", "ERR", "Scox",
    "mtTFB2", "ATPsynbeta", "mRpL4", "srl"), AMPK = c("Sesn", "FASN3", "Sirt1", "Atg1", "bmm", "Atg13", "Atg8a", "ACC", "SNF4Agamma",
    "alc", "FASN2", "S6k"), Chromatin = c("mor", "nej", "Su(var)3-9", "esc", "Caf1-105", "Gcn5", "Acf", "Snr1", "Caf1-180",
    "Pc", "HDAC1", "Su(z)12", "HP1c", "Chrac-14", "His2Av"), Checkpoint = c("mei-41", "Blm", "mus304", "RecQ4", "grp", "tefu",
    "mre11", "14-3-3epsilon", "mus101", "Ku80"))
plot_data <- bind_rows(lapply(names(module_genes), function(module_name) {
    tibble(module = module_name, gene = module_genes[[module_name]], display_order = seq_along(module_genes[[module_name]]))
})) %>% left_join(df_f %>% select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj), by = "gene") %>% mutate(status = case_when(!is.na(padj) &
    padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "Up in CS", !is.na(padj) & padj < padj_cutoff & log2FoldChange <=
    -lfc_cutoff ~ "Down in CS", TRUE ~ "Not significant"), evidence = pmin(50, -log10(pmax(padj, .Machine$double.xmin))),
    module = factor(module, levels = names(module_genes))) %>% group_by(module) %>% mutate(gene_display = factor(gene, levels = rev(unique(gene)))) %>%
    ungroup()
if (anyNA(plot_data$log2FoldChange)) {
    stop("One or more requested module genes are absent from the DESeq2 table.")
}
dna_padding_label <- "__DNA_TOP_PADDING__"
gene_levels_with_padding <- c(levels(plot_data$gene_display), dna_padding_label)
plot_data_chart <- bind_rows(plot_data %>% mutate(gene_display = factor(as.character(gene_display), levels = gene_levels_with_padding)),
    plot_data %>% filter(module == "DNA replication", gene == "timeout") %>% mutate(gene = dna_padding_label, gene_display = factor(dna_padding_label,
        levels = gene_levels_with_padding), log2FoldChange = NA_real_, evidence = NA_real_))
write_csv(plot_data %>% select(module, gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, status), .rebuild_file(file.path(output_dir,
    "FigS3_supporting_data.csv")), na = "")
status_colors <- c(`Up in CS` = "#A8322B", `Down in CS` = "#4C72B0", `Not significant` = "#B8B8B8")
observed_status <- c("Up in CS", "Down in CS", "Not significant")
observed_status <- observed_status[observed_status %in% unique(plot_data$status)]
p <- ggplot(plot_data_chart, aes(log2FoldChange, gene_display)) + geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
    linewidth = 0.45) + geom_point(aes(size = evidence, color = status), alpha = 0.92, na.rm = TRUE) + facet_wrap(~module,
    scales = "free_y", ncol = 2) + scale_x_continuous(limits = c(-2.1, 4.65), breaks = c(-2, 0, 2, 4), expand = expansion(mult = c(0.01,
    0.02))) + scale_y_discrete(breaks = function(x) x[x != dna_padding_label], expand = expansion(add = 0.6)) + scale_color_manual(values = status_colors,
    breaks = observed_status, name = "DE status") + scale_size_continuous(range = c(1.5, 5.2), limits = c(0, 50), breaks = c(10,
    20, 30, 40, 50), name = expression(-log[10] ~ "adjusted P")) + labs(title = "Expression across pre-specified gene modules",
    subtitle = "Female CS vs WT; n = 3 biological libraries per genotype", x = expression(log[2] ~ "fold change (CS/WT)"),
    y = NULL, caption = paste0("Significant: BH-adjusted P < 0.05 and |log2 fold change| >= 0.58. ", "Point size encodes -log10 adjusted P (capped at 50).")) +
    theme_classic(base_size = 10, base_family = "Helvetica") + theme(plot.title = element_text(size = 12, face = "bold",
    margin = margin(b = 3)), plot.subtitle = element_text(size = 9.2, color = "grey30", margin = margin(b = 6)), plot.caption = element_text(size = 7.6,
    color = "grey30", hjust = 0, margin = margin(t = 7)), axis.title.x = element_text(size = 10.5, margin = margin(t = 7)),
    axis.text.x = element_text(size = 8.5, color = "black"), axis.text.y = element_text(size = 7.4, color = "black", face = "italic"),
    axis.line = element_line(linewidth = 0.45, color = "black"), axis.ticks = element_line(linewidth = 0.45, color = "black"),
    strip.background = element_rect(fill = "grey94", color = "grey25", linewidth = 0.5), strip.text = element_text(size = 9.2,
        face = "bold", margin = margin(t = 3, b = 3)), panel.border = element_rect(fill = NA, color = "grey25", linewidth = 0.5),
    panel.spacing = grid::unit(4, "pt"), legend.position = "right", legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8.2), legend.key.height = grid::unit(11, "pt"), plot.margin = margin(10, 12, 9, 10)) +
    guides(color = guide_legend(order = 1, override.aes = list(size = 3.2)), size = guide_legend(order = 2))
pdf_file <- file.path(output_dir, "FigS3_publication.pdf")
png_file <- file.path(output_dir, "FigS3_publication.png")
tiff_file <- file.path(output_dir, "FigS3_publication_600dpi.tiff")
ggsave(pdf_file, p, width = 7.2, height = 9.4, device = "pdf", useDingbats = FALSE)
ggsave(png_file, p, width = 7.2, height = 9.4, dpi = 600, bg = "white")
ggsave(tiff_file, p, width = 7.2, height = 9.4, dpi = 600, units = "in", compression = "lzw", bg = "white")
legend_text <- c("Figure S3. Expression changes across pre-specified gene modules.", paste0("Female dMIC60-CS and dMIC60-WT bulk RNA-seq libraries were compared using ",
    "DESeq2 (n = 3 biological libraries per genotype). Points show ", "DESeq2 log2 fold-change estimates for the indicated pre-specified genes. Red ",
    "denotes significantly increased expression in CS, blue denotes significantly ", "decreased expression in CS, and gray denotes genes not meeting both thresholds ",
    "(Benjamini-Hochberg-adjusted P < 0.05 and absolute log2 fold change >= 0.58). Point ", "size represents -log10 adjusted P and is capped at 50 for display. These are ",
    "targeted gene-level summaries and not pathway-enrichment tests."))
writeLines(legend_text, .rebuild_file(file.path(output_dir, "FigS3_figure_legend.txt")))
message("Created publication package in: ", output_dir)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
