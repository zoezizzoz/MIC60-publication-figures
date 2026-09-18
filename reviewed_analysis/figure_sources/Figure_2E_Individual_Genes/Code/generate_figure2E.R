# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(DESeq2)
    library(patchwork)
    library(tidyverse)
})
args_all <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
set.seed(FIG_SEED)
input_dir <- file.path(figure_dir, "Local_Rebuild_Data")
support_dir <- file.path(figure_dir, "Supporting_Data")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
final_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
analysis_root <- normalizePath(file.path(figure_dir,"../.."))
count_mat <- as.matrix(read.csv(file.path(analysis_root,"data/rnaseq_counts.csv"),row.names=1,check.names=FALSE))
female_de <- read.csv(file.path(analysis_root,"data/rnaseq_female.csv"))
male_de <- read.csv(file.path(analysis_root,"data/rnaseq_male.csv"))
female_samples <- colnames(count_mat)[grepl("F[123]$",colnames(count_mat))]
male_samples <- colnames(count_mat)[grepl("M[123]$",colnames(count_mat))]
# Reuse the normalized counts already fitted once in 01_rnaseq.R.
female_fit <- list(normalized=as.matrix(read.csv(file.path(analysis_root,"results/rnaseq/female_normalized_counts.csv"),row.names=1,check.names=FALSE)))
male_fit <- list(normalized=as.matrix(read.csv(file.path(analysis_root,"results/rnaseq/male_normalized_counts.csv"),row.names=1,check.names=FALSE)))
genes_show <- c("timeout", "ImpL2", "DNAlig3")
stopifnot(all(genes_show %in% rownames(count_mat)))
make_long <- function(norm, samples, sex_label) {
    bind_rows(lapply(genes_show, function(g) data.frame(gene = g, Sample = samples, count = as.numeric(norm[g, samples]), 
        Genotype = ifelse(grepl("^CS", samples), "dMIC60-CS", "dMIC60-WT"), Sex = sex_label, stringsAsFactors = FALSE)))
}
long <- bind_rows(make_long(female_fit$normalized, female_samples, "Female"), make_long(male_fit$normalized, male_samples, 
    "Male"))
long$Genotype <- factor(long$Genotype, levels = c("dMIC60-WT", "dMIC60-CS"))
long$Sex <- factor(long$Sex, levels = c("Female", "Male"))
long$gene <- factor(long$gene, levels = genes_show)
long$PaletteGroup <- factor(paste(long$Sex, long$Genotype, sep = " | "), levels = c("Female | dMIC60-WT", "Female | dMIC60-CS", 
    "Male | dMIC60-WT", "Male | dMIC60-CS"))
lookup_padj <- function(g, sex_label) {
    d <- if (sex_label == "Female") 
        female_de
    else male_de
    d$padj[match(g, d$gene)]
}
sig_ann <- bind_rows(lapply(genes_show, function(g) bind_rows(lapply(c("Female", "Male"), function(s) {
    vals <- long$count[long$gene == g & long$Sex == s]
    p <- lookup_padj(g, s)
    span <- max(diff(range(vals)), abs(max(vals)) * 0.15, 1e-09)
    data.frame(gene = g, Sex = s, padj = p, y = max(vals) + span * FIG_BRACKET_OFFSET, tick = span * FIG_BRACKET_TICK)
})))) %>% filter(is.finite(padj), padj < FIG_ALPHA)
sig_ann$label <- vapply(sig_ann$padj, function(p) {
    if (p < 0.001) 
        "italic(p)[plain(BH)] < 0.001"
    else if (p < 0.01) 
        sprintf("italic(p)[plain(BH)] == %.3f", p)
    else sprintf("italic(p)[plain(BH)] == %.2f", p)
}, character(1))
sig_ann$gene <- factor(sig_ann$gene, levels = genes_show)
sig_ann$Sex <- factor(sig_ann$Sex, levels = c("Female", "Male"))
facet_brackets <- function(d, x1 = 1, x2 = 2) {
    if (nrow(d) == 0) 
        return(NULL)
    d$x1 <- x1
    d$x2 <- x2
    d$xmid <- (x1 + x2)/2
    list(geom_segment(data = d, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH), 
        geom_segment(data = d, aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH), 
        geom_segment(data = d, aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, linewidth = FIG_BRACKET_LINEWIDTH), 
        geom_text(data = d, aes(x = xmid, y = y + tick * 0.6, label = label), inherit.aes = FALSE, vjust = 0, size = FIG_ANNOT_SIZE, 
            family = FIG_FONT, parse = TRUE))
}
p <- ggplot(long, aes(Genotype, count, fill = PaletteGroup)) + fig_summary() + fig_points(color = FIG_PT_COLOR) + facet_brackets(sig_ann) + 
    scale_fill_manual(values = c(`Female | dMIC60-WT` = unname(FIG_GROUP_COLORS["WT_F"]), `Female | dMIC60-CS` = unname(FIG_GROUP_COLORS["CS_F"]), `Male | dMIC60-WT` = unname(FIG_GROUP_COLORS["WT_M"]), 
        `Male | dMIC60-CS` = unname(FIG_GROUP_COLORS["CS_M"])), guide = "none") + scale_y_continuous(expand = expansion(mult = c(FIG_Y_EXPAND[1], 
    FIG_Y_EXPAND[2] + 0.08))) + facet_grid(gene ~ Sex, scales = "free_y", switch = "y", labeller = labeller(Sex = c(Female = "Female", 
    Male = "Male"))) + labs(x = NULL, y = "DESeq2 normalized counts", title = "Individual gene expression", caption = paste0("Female and male comparisons were fitted separately with DESeq2; dMIC60-WT n = 3 and dMIC60-CS n = 3 per sex.\n", 
    "All libraries, including CSF3, were analyzed without manual gene prefiltering. ", "Wald-test P values were BH-adjusted genome-wide; adjusted P values are shown where padj < 0.05.")) + 
    theme_fig() + theme(panel.border = element_rect(fill = NA, color = FIG_FACET_BORDER_COLOR, linewidth = FIG_LINE_WIDTH), 
    panel.spacing.y = unit(FIG_FACET_GAP_PT, "pt"), strip.background = element_rect(fill = FIG_FACET_STRIP_FILL, color = FIG_FACET_BORDER_COLOR, 
        linewidth = FIG_LINE_WIDTH), strip.text.y.left = element_text(angle = 0, face = "bold.italic", size = FIG_AXIS_TEXT_SIZE + 
        2), strip.text.x = element_text(size = FIG_AXIS_TEXT_SIZE), axis.text.x = element_text(angle = 20, hjust = 1))
fig_save(p, file.path(output_dir, "Fig2D_individual_genes"), width = 7.2, height = 6.7)
horizontal_gene_order <- c("timeout", "DNAlig3", "ImpL2")
make_female_horizontal_panel <- function(gene_name, show_y_title = FALSE) {
    panel_data <- droplevels(subset(long, Sex == "Female" & gene == gene_name))
    panel_data$gene <- factor(as.character(panel_data$gene), levels = gene_name)
    panel_ann <- droplevels(subset(sig_ann, Sex == "Female" & gene == gene_name))
    panel_ann$gene <- factor(as.character(panel_ann$gene), levels = gene_name)
    ggplot(panel_data, aes(Genotype, count, fill = PaletteGroup)) + fig_summary() + fig_points(color = FIG_PT_COLOR) + facet_brackets(panel_ann) + 
        scale_fill_manual(values = c(`Female | dMIC60-WT` = unname(FIG_GROUP_COLORS["WT_F"]), `Female | dMIC60-CS` = unname(FIG_GROUP_COLORS["CS_F"])), guide = "none") + 
        scale_x_discrete(labels = c("dMIC60-WT", "dMIC60-CS"), expand = expansion(add = FIG_X_EXPAND)) + scale_y_continuous(position = "left", expand = expansion(mult = c(FIG_Y_EXPAND[1],
        FIG_Y_EXPAND[2] + 0.08))) + facet_wrap(~gene, scales = "free_y", nrow = 1) + labs(x = NULL, 
        y = NULL) + theme_fig() + theme(panel.border = element_rect(fill = NA, color = FIG_FACET_BORDER_COLOR, linewidth = FIG_LINE_WIDTH), 
        strip.background = element_rect(fill = FIG_FACET_STRIP_FILL, color = FIG_FACET_BORDER_COLOR, linewidth = FIG_LINE_WIDTH), 
        strip.text.x = element_text(angle = 0, face = "bold.italic", size = FIG_AXIS_TEXT_SIZE), axis.text.x = element_text(angle = 0, 
            hjust = 0.5, face = "plain", lineheight = 0.95), axis.title.x = element_text(face = "plain", size = FIG_AXIS_TEXT_SIZE, 
            margin = margin(t = 1)), plot.margin = margin(3, 3, 3, 3))
}
horizontal_panels <- lapply(seq_along(horizontal_gene_order), function(i) {
    make_female_horizontal_panel(horizontal_gene_order[[i]], show_y_title = i == 1L)
})
shared_y_title <- wrap_elements(full = grid::textGrob("DESeq2 normalized counts", x = grid::unit(0.70, "npc"), y = grid::unit(0.59, "npc"), rot = 90, gp = grid::gpar(fontsize = FIG_AXIS_TITLE_SIZE)), clip = FALSE)
p_horizontal <- shared_y_title | wrap_plots(horizontal_panels, nrow = 1)
p_horizontal <- p_horizontal + plot_layout(widths = c(0.045, 1))
fig_save(p_horizontal, file.path(output_dir, "Fig2D_individual_genes_female_horizontal"), width = 580/72, height = 140/72)
write.csv(long, .rebuild_file(file.path(support_dir, "Fig2D_normalized_counts_plotted.csv")), row.names = FALSE)
write.csv(sig_ann, .rebuild_file(file.path(support_dir, "Fig2D_adjusted_p_annotations.csv")), row.names = FALSE)
write.csv(female_de, .rebuild_file(file.path(support_dir, "DEG_CSF_vs_WRF.csv")), row.names = FALSE)
write.csv(male_de, .rebuild_file(file.path(support_dir, "DEG_CSM_vs_WRM.csv")), row.names = FALSE)
invisible(TRUE)
invisible(TRUE)
invisible(TRUE)
invisible(TRUE)
cat("Figure 2D rebuilt from the corrected sex-specific DESeq2 analyses.\n")
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
