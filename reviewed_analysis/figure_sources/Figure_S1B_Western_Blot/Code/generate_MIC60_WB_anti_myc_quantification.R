# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(ggplot2)
    library(readxl)
})
command_args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) != 1L) {
    stop("Run this file with Rscript so its package-relative paths can be resolved.")
}
script_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
package_dir <- dirname(script_dir)
data_file <- file.path(package_dir, "Original_Data", "myc WB quantifications.xlsx")
output_dir <- file.path(package_dir, "Rebuilt_Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(script_dir, "figure_style.R"))
FIG_SHOW_NS <- TRUE
plot_data <- read.csv(file.path(package_dir,"../../../Analysis/data/western_blot.csv"))
plot_data$genotype <- factor(plot_data$genotype,levels=c("dMIC60WT","dMIC60CS"))
pairs <- merge(subset(plot_data,genotype=="dMIC60WT"),subset(plot_data,genotype=="dMIC60CS"),by="blot_date",suffixes=c("_WT","_CS"))
wt <- pairs$normalized_MIC60_Myc_over_ATP5B_WT
cs <- pairs$normalized_MIC60_Myc_over_ATP5B_CS
x_positions <- c(dMIC60WT = 1, dMIC60CS = 2)
plot_data$x_position <- unname(x_positions[as.character(plot_data$genotype)])
paired_test <- t.test(wt, cs, paired = TRUE, alternative = "two.sided")
p_label <- sprintf('atop(bold("ns"), italic(p) == "%.3f")', paired_test$p.value)
figure <- ggplot(plot_data, aes(x = x_position, y = normalized_MIC60_Myc_over_ATP5B, group = genotype, fill = genotype)) + 
    fig_summary(p_value = paired_test$p.value) + fig_points(color = "black") + fig_sig_bracket(p_value = paired_test$p.value, 
    values = plot_data$normalized_MIC60_Myc_over_ATP5B, x1 = x_positions[["dMIC60WT"]], x2 = x_positions[["dMIC60CS"]], label = p_label,parse_label=TRUE) + 
    fig_scale_fill(c("dMIC60WT", "dMIC60CS")) + scale_x_continuous(breaks = unname(x_positions), labels = c("dMIC60-WT", "dMIC60-CS"), 
    limits = c(0.65, 2.35), expand = expansion(mult = 0)) + fig_scale_y("MIC60-Myc / ATP5B\n(normalized intensity)", limits = c(0, 
    NA), breaks = seq(0, 2, 0.5)) + labs(x = NULL, title = "MIC60-Myc protein abundance", subtitle = "Three matched blots; 5 larvae pooled per replicate", 
    caption = paste0("n = 3 biological replicates per genotype; each point represents one blot.\n", "Each replicate contains lysate pooled from 5 larvae.\n", 
        "Boxes: median/IQR.\n", "Whiskers: 1.5 x IQR; two-sided paired t-test.")) + coord_cartesian(clip = "off") + theme_fig() + 
    theme(axis.title = element_text(size = FIG_AXIS_TITLE_SIZE), axis.text = element_text(size = FIG_AXIS_TEXT_SIZE), axis.text.x = element_text(size = 8), plot.margin = margin(t = 10, 
        r = 9, b = 9, l = 12))
output_base <- file.path(output_dir, "MIC60_WB_anti_myc_quantification")
fig_save(figure, output_base, width = 130/72, height = 140/72, formats = c("png", "pdf"))
write.csv(plot_data, .rebuild_file(file.path(output_dir, "MIC60_WB_graph_ready_data.csv")), row.names = FALSE)
summary_statistics <- data.frame(comparison = "dMIC60WT vs dMIC60CS", n_biological_replicates_per_genotype = length(wt), 
    larvae_pooled_per_replicate = 5L, dMIC60WT_mean = mean(wt), dMIC60WT_median = median(wt), dMIC60WT_sd = sd(wt), dMIC60CS_mean = mean(cs), 
    dMIC60CS_median = median(cs), dMIC60CS_sd = sd(cs), paired_t_statistic = unname(paired_test$statistic), paired_t_df = unname(paired_test$parameter), 
    paired_t_p_value = paired_test$p.value, stringsAsFactors = FALSE)
write.csv(summary_statistics, .rebuild_file(file.path(output_dir, "MIC60_WB_summary_statistics.csv")), row.names = FALSE)
message("Wrote publication outputs to: ", output_dir)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
