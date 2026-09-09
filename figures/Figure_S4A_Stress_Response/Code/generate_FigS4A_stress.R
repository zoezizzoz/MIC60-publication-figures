# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggrepel)
    library(readxl)
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
final_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
PADJ_CUTOFF <- 0.05
LFC_CUTOFF <- 0.58
SEED <- 1
up_col <- "#B40426"
down_col <- "#3B4CC0"
pubtheme <- theme_fig(legend_position = "right")
female <- read.csv(file.path(input_dir, "DEG_CSF_vs_WRF.csv"), check.names = FALSE)
male <- read.csv(file.path(input_dir, "DEG_CSM_vs_WRM.csv"), check.names = FALSE)
female <- female[!duplicated(female$gene), ]
male <- male[!duplicated(male$gene), ]
rownames(female) <- female$gene
rownames(male) <- male$gene
old_sets <- read.csv(file.path(support_dir, "figS4_stress_gene_mapping.csv"), check.names = FALSE)
if (!all(c("pathway", "detected_symbol", "in_set") %in% names(old_sets))) {
    stop("The curated stress mapping lacks pathway/detected_symbol/in_set columns.")
}
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
        1), odds_ratio = round(unname(ft$estimate), 2), fisher_p = signif(ft$p.value, 3))
}
enr <- bind_rows(lapply(set_order, function(pw) cbind(pathway = pw, fisher_path(stress_sets[[pw]]))))
enr$fisher_padj <- signif(p.adjust(enr$fisher_p, method = "BH"), 3)
write.csv(enr, .rebuild_file(file.path(support_dir, "figS4_enrichment.csv")), row.names = FALSE)
stress_plot_dat <- stress_long %>% mutate(sig = case_when(!is.na(padj_F) & padj_F < PADJ_CUTOFF & lfc_F >= LFC_CUTOFF ~ "up (padj<0.05, |LFC|≥0.58)",
    !is.na(padj_F) & padj_F < PADJ_CUTOFF & lfc_F <= -LFC_CUTOFF ~ "down", is.na(padj_F) ~ "padj filtered (NA)", TRUE ~ "n.s."),
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
stress_plot_dat$lab <- ifelse(tolower(stress_plot_dat$gene) %in% tolower(common_genes) | stress_plot_dat$sig %in% c("up (padj<0.05, |LFC|≥0.58)",
    "down") | (stress_plot_dat$pathway == "Control" & stress_plot_dat$gene %in% control_labels), stress_plot_dat$gene, NA_character_)
smeds <- stress_plot_dat %>% group_by(pathway) %>% summarise(m = median(lfc_F, na.rm = TRUE), .groups = "drop")
lfc_range <- range(stress_plot_dat$lfc_F, na.rm = TRUE)
y_annot <- lfc_range[2] + max(diff(lfc_range), 1) * 0.1
enr2 <- enr %>% mutate(pathway = factor(pathway, levels = set_order), lab = ifelse(is.na(fisher_padj), "not testable", ifelse(fisher_padj <
    0.001, ifelse(is.infinite(odds_ratio), sprintf("OR ∞\nadjusted P=%.0e", fisher_padj), sprintf("OR %.1f\nadjusted P=%.0e",
    odds_ratio, fisher_padj)), ifelse(is.infinite(odds_ratio), sprintf("OR ∞\nadjusted P=%.2f", fisher_padj), sprintf("OR %.2f\nadjusted P=%.2f",
    odds_ratio, fisher_padj)))), y = y_annot)
labels <- c(Control = "Non-stress\ncontrol", ISR = "ISR", `UPR (ATF6)` = "ER-UPR\n(ATF6)", `UPR (IRE1/XBP1s)` = "ER-UPR\n(IRE1/XBP1s)",
    HSR = "HSR", OSR = "OSR")
set.seed(FIG_SEED)
stress_plot_dat$x_plot <- as.numeric(stress_plot_dat$pathway) + runif(nrow(stress_plot_dat), -(FIG_JITTER_WIDTH + 0.04),
    FIG_JITTER_WIDTH + 0.04)
smeds$x_plot <- match(as.character(smeds$pathway), set_order)
enr2$x_plot <- match(as.character(enr2$pathway), set_order)
p <- ggplot(stress_plot_dat, aes(x_plot, lfc_F)) + annotate("rect", xmin = 0.4, xmax = length(set_order) + 0.6, ymin = -LFC_CUTOFF,
    ymax = LFC_CUTOFF, fill = "grey90", alpha = 0.6) + geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) + geom_hline(yintercept = c(-LFC_CUTOFF,
    LFC_CUTOFF), linetype = 2, color = "grey60", linewidth = 0.3) + geom_point(aes(color = sig), size = FIG_PT_SIZE, alpha = 0.85) +
    geom_crossbar(data = smeds, aes(x = x_plot, y = m, ymin = m, ymax = m), width = FIG_SUMMARY_WIDTH, color = FIG_BOX_OUTLINE,
        linewidth = FIG_LINE_WIDTH, inherit.aes = FALSE) + geom_label_repel(aes(label = lab, color = sig), size = FIG_ANNOT_SIZE -
    0.7, fontface = "italic", fill = scales::alpha("white", 0.92), label.size = NA, label.padding = grid::unit(0.1, "lines"),
    label.r = grid::unit(0.05, "lines"), max.overlaps = Inf, max.time = 5, min.segment.length = 0, segment.size = 0.22, segment.color = "grey55",
    box.padding = 0.55, point.padding = 0.22, force = 5, force_pull = 1.2, seed = SEED, na.rm = TRUE, show.legend = FALSE) +
    geom_text(data = enr2, aes(x = x_plot, y = y, label = lab), size = FIG_ANNOT_SIZE - 0.3, lineheight = 0.9, color = ifelse(!is.na(enr2$fisher_padj) &
        enr2$fisher_padj < PADJ_CUTOFF, up_col, "grey35"), inherit.aes = FALSE) + scale_color_manual(values = c(`up (padj<0.05, |LFC|≥0.58)` = up_col,
    down = down_col, n.s. = "grey65", `padj filtered (NA)` = "grey82"), name = NULL, breaks = c("up (padj<0.05, |LFC|≥0.58)",
    "down", "n.s.", "padj filtered (NA)")) + scale_x_continuous(breaks = seq_along(set_order), labels = labels, limits = c(0.4,
    length(set_order) + 0.6), expand = expansion(mult = 0)) + scale_y_continuous(expand = expansion(mult = c(0.03, 0.14))) +
    labs(x = NULL, y = expression(log[2] ~ fold ~ change ~ (("dMIC60-CS" ~ "♀")/("dMIC60-WT" ~ "♀"))), title = "Female dMIC60-CS-versus-WT targeted stress-response signatures",
        subtitle = paste0("Corrected female-only DESeq2 analysis; ", length(testable_genes), " testable curated genes; BH-adjusted across six sets"),
        caption = "Female dMIC60-WT n = 3 and dMIC60-CS n = 3 biological libraries; all female libraries, including CSF3, were analyzed.") +
    pubtheme + theme(legend.position = "top", axis.text.x = element_text(size = FIG_AXIS_TEXT_SIZE - 1, lineheight = 0.85)) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.6, label = "")))
fig_save(p, file.path(output_dir, "figS4_stress_distinct"), width = max(FIG_W_2COL, 1.3 * length(set_order) + 1.5), height = 6.2)
write.csv(stress_plot_dat, .rebuild_file(file.path(support_dir, "figS4_plotted_values.csv")), row.names = FALSE)
write.csv(data.frame(metric = c("testable_stress_genes", "female_up_DEGs_in_stress_universe"), value = c(length(testable_genes),
    length(up_all))), .rebuild_file(file.path(support_dir, "figS4_analysis_summary.csv")), row.names = FALSE)
invisible(TRUE)
invisible(TRUE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(output_dir, "sessionInfo.txt")))
print(enr, row.names = FALSE)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
