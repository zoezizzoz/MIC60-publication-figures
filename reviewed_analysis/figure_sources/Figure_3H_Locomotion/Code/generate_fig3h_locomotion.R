# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub('~+~',' ',commandArgs(FALSE),fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(ggplot2)
    library(readxl)
})
args <- gsub('~+~',' ',commandArgs(trailingOnly = FALSE),fixed=TRUE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "generate_fig3h_locomotion.R"
code_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
project_dir <- dirname(code_dir)
source(file.path(code_dir, "figure_style.R"))
FIG_ANNOT_SIZE <- 5.2
input_path <- file.path(project_dir, "Original_Data", "Fig3H-locomotion.xlsx")
final_dir <- file.path(project_dir, "Rebuilt_Output")
rebuilt_dir <- file.path(project_dir, "Rebuilt_Output")
support_dir <- file.path(project_dir, "Supporting_Data")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rebuilt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(support_dir, recursive = TRUE, showWarnings = FALSE)
raw <- read_excel(input_path)
if (ncol(raw) != 2L) stop("Expected two columns in Fig3H-locomotion.xlsx.")
groups <- c("WT", "CS")
dat <- do.call(rbind, lapply(seq_along(groups), function(i) {
    values <- suppressWarnings(as.numeric(raw[[i]]))
    values <- values[is.finite(values)]
    data.frame(group = factor(groups[[i]], levels = groups), performance_index = values, observation = seq_along(values))
}))
x_positions <- c(WT = 1, CS = 1.46)
dat$x_position <- unname(x_positions[as.character(dat$group)])
box_stats <- do.call(rbind, lapply(groups, function(group_name) {
    values <- dat$performance_index[dat$group == group_name]
    five <- fivenum(values)
    data.frame(group = factor(group_name, levels = groups), x_position = unname(x_positions[[group_name]]), ymin = five[[1]], 
        lower = five[[2]], middle = five[[3]], upper = five[[4]], ymax = five[[5]])
}))
wt <- dat$performance_index[dat$group == "WT"]
cs <- dat$performance_index[dat$group == "CS"]
wilcox_result <- wilcox.test(wt, cs, exact = FALSE, correct = TRUE)
all_values <- c(wt, cs)
assignments <- combn(seq_along(all_values), length(wt))
observed_difference <- abs(mean(cs) - mean(wt))
permuted_differences <- apply(assignments, 2, function(ix) abs(mean(all_values[ix]) - mean(all_values[-ix])))
p_exact <- mean(permuted_differences >= observed_difference - 1e-12)
stopifnot(ncol(assignments) == 70L, abs(p_exact - 2/70) < 1e-12)
p_reference <- p_exact
tests <- data.frame(comparison = "WT vs CS; Day 2", test = c("Reference annotation (supplied image)", "Wilcoxon rank-sum test with continuity correction", 
    "Exact two-sided permutation of mean difference", "Welch two-sample t-test", "Student two-sample t-test"), p_value = c(p_reference, wilcox_result$p.value, p_exact, t.test(wt, cs)$p.value,
    t.test(wt, cs, var.equal = TRUE)$p.value), stringsAsFactors = FALSE)
write.csv(dat, .rebuild_file(file.path(support_dir, "Fig3H_locomotion_plotted_values.csv")), row.names = FALSE)
write.csv(tests, .rebuild_file(file.path(support_dir, "Fig3H_locomotion_statistical_audit.csv")), row.names = FALSE)
write.csv(box_stats, .rebuild_file(file.path(support_dir, "Fig3H_locomotion_Tukey_hinges.csv")), row.names = FALSE)
p <- ggplot(dat, aes(x = x_position, y = performance_index, group = group, fill = group)) + geom_segment(data = box_stats, 
    aes(x = x_position, xend = x_position, y = ymin, yend = lower), inherit.aes = FALSE, linewidth = 0.8, color = "grey20") + 
    geom_segment(data = box_stats, aes(x = x_position, xend = x_position, y = upper, yend = ymax), inherit.aes = FALSE, linewidth = 0.8, 
        color = "grey20") + geom_rect(data = box_stats, aes(xmin = x_position - 0.15, xmax = x_position + 0.15, ymin = lower, 
    ymax = upper, fill = group), inherit.aes = FALSE, linewidth = 0.8, color = "grey20") + geom_segment(data = box_stats, 
    aes(x = x_position - 0.15, xend = x_position + 0.15, y = middle, yend = middle), inherit.aes = FALSE, linewidth = 0.8, 
    color = "grey20") + geom_segment(data = box_stats, aes(x = x_position - 0.08, xend = x_position + 0.08, y = ymin, yend = ymin), 
    inherit.aes = FALSE, linewidth = 0.8, color = "black") + geom_segment(data = box_stats, aes(x = x_position - 0.08, xend = x_position + 
    0.08, y = ymax, yend = ymax), inherit.aes = FALSE, linewidth = 0.8, color = "black") + geom_point(position = position_jitter(width = 0.035, 
    height = 0, seed = FIG_SEED), shape = FIG_PT_SHAPE, size = 2.6, alpha = 0.8, color = FIG_PT_COLOR, show.legend = FALSE) + 
    fig_sig_bracket(p_value = p_reference, values = dat$performance_index, x1 = x_positions[["WT"]], x2 = x_positions[["CS"]], 
        label = expression(italic(p) == 0.0286)) + fig_scale_fill(groups) + scale_x_continuous(breaks = unname(x_positions), 
    labels = c("dMIC60-WT", "dMIC60-CS"), limits = c(0.7, 1.76), expand = expansion(mult = 0)) + fig_scale_y(name = "Performance Index", 
    limits = c(0, 0.42), breaks = seq(0, 0.4, 0.1)) + labs(x = NULL) + guides(fill = "none") + theme_fig(base_size = 19) + 
    theme(axis.title.y = element_text(size = 22, margin = margin(r = 9)), axis.text.y = element_text(size = 18), axis.text.x = element_text(size = 16, 
        margin = margin(t = 7)), plot.margin = margin(t = 10, r = 9, b = 9, l = 12))
for (out_dir in c(final_dir, rebuilt_dir)) {
    fig_save(p, file.path(out_dir, "Fig3H_locomotion_day2"), width = 4, height = 4.5)
}
message("Created Fig3H Day 2 locomotion plot and statistical audit.")
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
