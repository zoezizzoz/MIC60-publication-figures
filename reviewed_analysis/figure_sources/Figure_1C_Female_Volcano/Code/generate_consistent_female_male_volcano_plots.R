# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
})
args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
    normalizePath("RNA-seq")
}
source(file.path(script_dir,"figure_style.R"))
analysis_root <- normalizePath(file.path(script_dir,"../../.."))
output_dir <- file.path(analysis_root,"results/full_figures")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
lfc_thresh <- 0.58
padj_thresh <- 0.05
read_de <- function(filename, sex) {
    data <- read.csv(file.path(analysis_root,"data",if(sex=="Female") "rnaseq_female.csv" else "rnaseq_male.csv"),row.names=1,check.names=FALSE)
    data %>% filter(!is.na(padj)) %>% mutate(sex = sex, neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)), direction = case_when(padj < 
        padj_thresh & log2FoldChange >= lfc_thresh ~ "Up in CS", padj < padj_thresh & log2FoldChange <= -lfc_thresh ~ "Down in CS", 
        TRUE ~ "Not significant"))
}
female <- read_de("DEG_female_CS_vs_WT.csv", "Female")
male <- read_de("DEG_male_CS_vs_WT.csv", "Male")
stopifnot(nrow(female) == 7456L, sum(female$direction == "Up in CS") == 98L, sum(female$direction == "Down in CS") == 108L, 
    nrow(male) == 9695L, sum(male$direction == "Up in CS") == 540L, sum(male$direction == "Down in CS") == 505L)
combined <- bind_rows(female, male)
x_limit <- ceiling(max(abs(combined$log2FoldChange), na.rm = TRUE))
y_limit <- ceiling(max(combined$neg_log10_padj, na.rm = TRUE) * 1.05)
colors <- c(`Up in CS` = "#D62728", `Down in CS` = "#4D4398", `Not significant` = "#B5B5B5")
make_volcano <- function(data, sex) {
    panel_title <- if (sex == "Female") expression(italic("dMIC60") * "-null female") else expression(italic("dMIC60") * "-null male")
    colors <- c(`Up in CS`=unname(FIG_GROUP_COLORS[if(sex=="Female") "CS_F" else "CS_M"]),
      `Down in CS`=unname(FIG_GROUP_COLORS[if(sex=="Female") "WT_F" else "WT_M"]),`Not significant`="#B5B5B5")
    up_n <- sum(data$direction == "Up in CS")
    down_n <- sum(data$direction == "Down in CS")
    testable_n <- nrow(data)
    labels <- bind_rows(data %>% filter(direction == "Up in CS") %>% arrange(padj) %>% slice_head(n = 6), data %>% filter(direction == 
        "Down in CS") %>% arrange(padj) %>% slice_head(n = 6))
    p <- ggplot(data, aes(log2FoldChange, neg_log10_padj)) + geom_point(aes(color = direction), size = 0.55, alpha = 0.8) + 
        geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dotdash", color = "#9E9E9E", linewidth = FIG_LINE_WIDTH) + 
        geom_hline(yintercept = -log10(padj_thresh), linetype = "dotdash", color = "#9E9E9E", linewidth = FIG_LINE_WIDTH) + geom_label_repel(data = labels, 
        aes(label = gene), color = "black", fill = "white", fontface = "italic", family=FIG_FONT, size = FIG_ANNOT_SIZE, box.padding = 0.45, point.padding = 0.2, 
        min.segment.length = 0, segment.color = "black", segment.size = 0.35, max.overlaps = Inf, max.time = 3, force = 1.1, 
        ylim = c(0, y_limit * 0.96), seed = 60, show.legend = FALSE) + scale_color_manual(values = colors, breaks = c("Up in CS", 
        "Down in CS", "Not significant")) + annotate("text", x = -x_limit * 0.95, y = y_limit*.96, label = sprintf("Down: %d (%.1f%%)", 
        down_n, 100 * down_n/testable_n), hjust = 0, color = colors[["Down in CS"]], size = FIG_ANNOT_SIZE) + annotate("text", x = x_limit * 
        0.95, y = y_limit*.96, label = sprintf("Up: %d (%.1f%%)", up_n, 100 * up_n/testable_n), hjust = 1, color = colors[["Up in CS"]], 
        size = FIG_ANNOT_SIZE) + coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(0, y_limit), clip = "on") + labs(title = panel_title, subtitle = expression(paste("|log"[2], " FC| cutoff: 0.58   adjusted ", italic(p), " cutoff: 0.05")), 
        x = expression(log[2]~"fold change (" * italic("dMIC60") * "-CS/" * italic("dMIC60") * "-WT)"), y = expression(-log[10](italic(p)[plain(BH)])), 
        caption = bquote(.(testable_n)~"genes with adjusted "*italic(p)*" values"), color = NULL) + theme_fig() + theme(plot.title=element_text(size=8,face="bold",hjust=0),plot.subtitle=element_text(size=7),plot.caption=element_text(size=7),plot.margin=margin(4,4,4,4))

    p
}
female_plot <- make_volcano(female, "Female")
male_plot <- make_volcano(male, "Male")
fig_save(female_plot,file.path(output_dir,"Fig1C_Female_Volcano"),width=260/72,height=255/72)
fig_save(male_plot,file.path(output_dir,"FigS1C_Male_Volcano"),width=260/72,height=255/72)
.rebuild_dir <- output_dir
summary <- data.frame(Figure = c("Figure 1C", "Figure S1C"), Sex = c("Female", "Male"), Contrast = c("Female CS vs WT", "Male CS vs WT"), 
    Testable_genes = c(nrow(female), nrow(male)), Up_in_CS = c(sum(female$direction == "Up in CS"), sum(male$direction == 
        "Up in CS")), Down_in_CS = c(sum(female$direction == "Down in CS"), sum(male$direction == "Down in CS")), LFC_threshold = lfc_thresh, 
    adjusted_P_threshold = padj_thresh)
write.csv(summary, .rebuild_file(file.path(output_dir, "volcano_DEG_count_summary.csv")), row.names = FALSE)
print(summary)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
