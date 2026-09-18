# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~"," ",commandArgs(FALSE),fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
    library(readxl)
})
args_all <- gsub("~+~"," ",commandArgs(trailingOnly=FALSE),fixed=TRUE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))
# Fig. S4A is reduced when placed in the assembled Illustrator figure, so use
# a larger panel-specific text scale to keep every label legible at final size.
FIG_AXIS_TITLE_SIZE <- 11
FIG_AXIS_TEXT_SIZE <- 10
FIG_LEGEND_TEXT_SIZE <- 10
FIG_ANNOT_SIZE <- 9 / (72.27 / 25.4)
set.seed(FIG_SEED)
input_dir <- file.path(figure_dir, "Original_Data")
support_dir <- file.path(figure_dir, "Supporting_Data")
output_dir <- file.path(figure_dir, "Rebuilt_Output")
final_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 0.58
SEED <- 1
up_col <- unname(FIG_HEATMAP_SCALE["high"])
down_col <- unname(FIG_GROUP_COLORS["WT_F"])
pubtheme <- theme_fig(legend_position = "right")
female <- read.csv(file.path(input_dir, "DEG_CSF_vs_WRF.csv"), check.names = FALSE)
male <- read.csv(file.path(input_dir, "DEG_CSM_vs_WRM.csv"), check.names = FALSE)
female <- female[!duplicated(female$gene), ]
male <- male[!duplicated(male$gene), ]
rownames(female) <- female$gene
rownames(male) <- male$gene
# Recover identical curated memberships from the supplied plotted gene table.
curated <- read.csv(file.path(figure_dir,"../../data/stress_genes.csv"))
old_sets <- data.frame(pathway=curated$pathway,detected_symbol=curated$gene,in_set=TRUE)
set_order <- c("Control", "ISR", "UPR (ATF6)", "UPR (IRE1/XBP1s)", "HSR", "OSR")
stress_sets <- lapply(set_order, function(pw) {
    keep <- old_sets$pathway == pw & old_sets$in_set & !is.na(old_sets$detected_symbol) & nzchar(old_sets$detected_symbol)
    sort(unique(old_sets$detected_symbol[keep]))
})
names(stress_sets) <- set_order
stress_long <- bind_rows(lapply(set_order, function(pw) {
    genes <- stress_sets[[pw]]
    data.frame(pathway = pw, gene = genes, lfc_F = female[genes, "log2FoldChange"], padj_F = female[genes, "padj"], lfc_M = male[genes, 
        "log2FoldChange"], padj_M = male[genes, "padj"], stringsAsFactors = FALSE)
}))
write.csv(stress_long, .rebuild_file(file.path(support_dir, "figS4_stress_pathways.csv")), row.names = FALSE)
stress_universe <- sort(unique(unlist(stress_sets, use.names = FALSE)))
testable_genes <- intersect(stress_universe, female$gene[!is.na(female$padj)])
up_all <- intersect(testable_genes, female$gene[!is.na(female$padj) & female$padj < PADJ_CUTOFF & female$log2FoldChange >= 
    LFC_CUTOFF])
fisher_path <- function(genes) {
    n_detected <- length(genes)
    genes <- intersect(genes, testable_genes)
    if (length(genes) == 0) 
        return(data.frame(n_detected = n_detected, n_tested = 0, n_up = 0, pct_up = NA, odds_ratio = NA, fisher_p = NA))
    in_up <- sum(genes %in% up_all)
    out_genes <- setdiff(testable_genes, genes)
    ft <- fisher.test(matrix(c(in_up, length(genes) - in_up, sum(out_genes %in% up_all), sum(!out_genes %in% up_all)), 2), 
        alternative = "greater")
    data.frame(n_detected = n_detected, n_tested = length(genes), n_up = in_up, pct_up = round(100 * in_up/length(genes), 
        1), odds_ratio = round(unname(ft$estimate), 2), fisher_p = ft$p.value)
}
enr <- bind_rows(lapply(set_order, function(pw) cbind(pathway = pw, fisher_path(stress_sets[[pw]]))))
enr$fisher_padj <- p.adjust(enr$fisher_p, method = "BH")
write.csv(enr, .rebuild_file(file.path(support_dir, "figS3_enrichment.csv")), row.names = FALSE)
stress_plot_dat <- stress_long %>% mutate(sig = case_when(!is.na(padj_F) & padj_F < PADJ_CUTOFF & lfc_F >= LFC_CUTOFF ~ "Up in CS", 
    !is.na(padj_F) & padj_F < PADJ_CUTOFF & lfc_F <= -LFC_CUTOFF ~ "Down in CS", is.na(padj_F) ~ "Not tested", TRUE ~ "Not significant"), 
    pathway = factor(pathway, levels = set_order))
common_raw <- as.data.frame(read_excel(file.path(input_dir, "common stress genes.xlsx"), sheet = "Common Stress Genes", skip = 4), 
    check.names = FALSE)
names(common_raw) <- trimws(names(common_raw))
split_symbols <- function(x) {
    x <- trimws(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    trimws(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}
common_genes <- unique(c(split_symbols(common_raw[["Recommended fly ortholog(s)"]]), split_symbols(common_raw[["Other best-score matches"]])))
control_labels <- c("betaGlu", "mRpL9", "mtSSB", "Prosbeta6")
stress_plot_dat$lab <- ifelse(tolower(stress_plot_dat$gene) %in% tolower(common_genes) | stress_plot_dat$sig %in% c("Up in CS", 
    "Down in CS") | (stress_plot_dat$pathway == "Control" & stress_plot_dat$gene %in% control_labels), stress_plot_dat$gene, NA_character_)
smeds <- stress_plot_dat %>% group_by(pathway) %>% summarise(m = median(lfc_F, na.rm = TRUE), .groups = "drop")
lfc_range <- range(stress_plot_dat$lfc_F, na.rm = TRUE)
y_annot <- lfc_range[2] + max(diff(lfc_range), 1) * 0.1
enr2 <- enr %>% mutate(pathway=factor(pathway,levels=set_order),y=y_annot)
enr2$lab <- vapply(seq_len(nrow(enr2)),function(i) {
 or <- if(is.infinite(enr2$odds_ratio[i])) "infinity" else sprintf('"%.2f"',enr2$odds_ratio[i])
 paste0('atop(plain(OR)==',or,',italic(p)[plain(BH)]=="',sprintf('%.2f',enr2$fisher_padj[i]),'")')
},character(1))
labels <- c(Control = "Non-stress\ncontrol", ISR = "ISR", `UPR (ATF6)` = "ER-UPR\n(ATF6)", `UPR (IRE1/XBP1s)` = "ER-UPR\n(IRE1/XBP1s)", 
    HSR = "HSR", OSR = "OSR")
set.seed(FIG_SEED)
stress_plot_dat$x_plot <- as.numeric(stress_plot_dat$pathway) + runif(nrow(stress_plot_dat), -(FIG_JITTER_WIDTH + 0.04), 
    FIG_JITTER_WIDTH + 0.04)
smeds$x_plot <- match(as.character(smeds$pathway), set_order)
enr2$x_plot <- match(as.character(enr2$pathway), set_order)
p <- ggplot(stress_plot_dat, aes(x_plot, lfc_F)) + annotate("rect", xmin = 0.4, xmax = length(set_order) + 0.6, ymin = -LFC_CUTOFF, 
    ymax = LFC_CUTOFF, fill = "grey90", alpha = 0.6) + geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) + geom_hline(yintercept = c(-LFC_CUTOFF, 
    LFC_CUTOFF), linetype = 2, color = "grey60", linewidth = 0.3) + geom_point(aes(color = sig), size = FIG_PT_SIZE, alpha = 1) + 
    geom_crossbar(data = smeds, aes(x = x_plot, y = m, ymin = m, ymax = m), width = FIG_SUMMARY_WIDTH, color = FIG_BOX_OUTLINE, 
        linewidth = FIG_LINE_WIDTH, inherit.aes = FALSE) + geom_label_repel(aes(label = lab, color = sig), size = FIG_ANNOT_SIZE, fontface = "italic", fill = scales::alpha("white", 0.92), linewidth = 0, label.padding = grid::unit(0.1, "lines"), 
    label.r = grid::unit(0.05, "lines"), max.overlaps = Inf, max.time = 10, min.segment.length = 0, segment.size = 0.22, segment.color = "grey55", 
    box.padding = 0.55, point.padding = 0.22, force = 5, force_pull = 1.2, seed = SEED, na.rm = TRUE, show.legend = FALSE) + 
    geom_text(data = enr2, aes(x = x_plot, y = y, label = lab), size = FIG_ANNOT_SIZE, parse=TRUE, family=FIG_FONT, lineheight = 0.9, color = ifelse(!is.na(enr2$fisher_padj) & 
        enr2$fisher_padj < PADJ_CUTOFF, up_col, "grey35"), inherit.aes = FALSE) + scale_color_manual(values = c(`Up in CS` = up_col, 
    `Down in CS` = down_col, `Not significant` = "grey65", `Not tested` = "grey82"), name = NULL, breaks = c("Up in CS", 
    "Down in CS", "Not significant", "Not tested")) + scale_x_continuous(breaks = seq_along(set_order), labels = labels, limits = c(0.4, 
    length(set_order) + 0.6), expand = expansion(mult = 0)) + scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) + 
    labs(x = NULL, y = expression(log[2]~"fold change (" * italic("dMIC60") * "-CS/" * italic("dMIC60") * "-WT)"), title = NULL, 
        subtitle = paste0("Corrected female-only DESeq2 analysis; ", length(testable_genes), " testable curated genes; BH-adjusted across six sets"), 
        caption = "Female dMIC60-WT n = 3 and dMIC60-CS n = 3 biological libraries; all female libraries, including CSF3, were analyzed.") + 
    pubtheme + theme(plot.title=element_text(family=FIG_FONT,face="bold",size=11,hjust=0), legend.position = "top", legend.box.spacing=grid::unit(2,"pt"),
      legend.text=element_text(size=FIG_LEGEND_TEXT_SIZE), legend.margin=margin(0,0,0,0),
      axis.title.y=element_text(size=FIG_AXIS_TITLE_SIZE), axis.text=element_text(size=FIG_AXIS_TEXT_SIZE),
      axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE, lineheight = 0.85),
      plot.margin=margin(t=12,r=6,b=6,l=6)) + 
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.6, label = "")))
fig_save(p, file.path(output_dir, "figS4_stress_distinct"), width = 548/72, height = 365/72)
write.csv(stress_plot_dat, .rebuild_file(file.path(support_dir, "figS3_plotted_values.csv")), row.names = FALSE)
write.csv(data.frame(metric = c("testable_stress_genes", "female_up_DEGs_in_stress_universe"), value = c(length(testable_genes), 
    length(up_all))), .rebuild_file(file.path(support_dir, "figS3_analysis_summary.csv")), row.names = FALSE)
invisible(TRUE)
invisible(TRUE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(output_dir, "sessionInfo.txt")))
print(enr, row.names = FALSE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
