# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(stringr)
})
script_arg <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
figure_dir <- dirname(script_dir)
source(file.path(script_dir,"figure_style.R"))
plot_data <- read_csv(file.path(figure_dir, "Supporting_Data", "FigS2_plotted_pathways.csv"), show_col_types = FALSE) %>% 
    mutate(Direction = if_else(NES > 0, "Higher in CS", "Higher in WT"), facet_label = recode(comparison, CSF_vs_WTF = "'Female'",
        CSM_vs_WTM = "'Male'"), PaletteGroup = if_else(comparison == "CSF_vs_WTF", if_else(NES>0,"CS_F","WT_F"),if_else(NES>0,"CS_M","WT_M")), pathway_key = paste(comparison, str_wrap(Description, 32), sep = "___")) %>%
    arrange(facet_label, NES) %>% mutate(pathway_key = factor(pathway_key, levels = unique(pathway_key)))
# Display-only refinement: input pathways, NES, adjusted P and order are unchanged.
pathway_labels <- function(keys) {
    labels <- sub("^.*___", "", keys)
    full <- gsub("\n", " ", labels, fixed = TRUE)
    line_breaks <- c(
        "Metabolism of xenobiotics by cytochrome P450" = "Metabolism of xenobiotics\nby cytochrome P450",
        "Drug metabolism - cytochrome P450" = "Drug metabolism -\ncytochrome P450",
        "Drug metabolism - other enzymes" = "Drug metabolism -\nother enzymes",
        "Glycine, serine and threonine metabolism" = "Glycine, serine and\nthreonine metabolism",
        "Biosynthesis of unsaturated fatty acids" = "Biosynthesis of\nunsaturated fatty acids"
    )
    match_index <- match(full, names(line_breaks))
    labels[!is.na(match_index)] <- unname(line_breaks[match_index[!is.na(match_index)]])
    labels
}
plot <- ggplot(plot_data, aes(NES, pathway_key)) +
    geom_vline(xintercept = c(-2, -1, 1, 2), colour = "grey92", linewidth = 0.25) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.55) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_key, colour = PaletteGroup),
                 linewidth = 0.7, lineend = "round") +
    geom_point(aes(size = -log10(p.adjust), colour = PaletteGroup)) +
    facet_wrap(~facet_label, scales = "free_y", ncol = 2, labeller = label_parsed) +
    scale_y_discrete(labels = pathway_labels,
                     expand = expansion(add = 0.55)) +
    scale_x_continuous(breaks = seq(-2, 2, 1), expand = expansion(mult = c(0.06, 0.06))) +
    scale_colour_manual(values = FIG_GROUP_COLORS, guide = "none") +
    scale_size_continuous(name = expression(-log[10](italic(p)[plain(BH)])),
                          range = c(1, 3), breaks = c(2, 4, 6)) +
    labs(x = "Normalized enrichment score (NES)", y = NULL) +
    theme_fig(legend_position = "bottom") +
    theme(
        strip.text = element_text(size = 9, face = "bold", margin = margin(b = 7)),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 8, lineheight = 0.95, margin = margin(r = 5)),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = "grey55", linewidth = 0.3),
        axis.title.x = element_text(size = 8, lineheight = 1.35, margin = margin(t = 5)),
        legend.title = element_text(size = 7), legend.text = element_text(size = 7),
        legend.key.width = grid::unit(12, "pt"),
        legend.spacing.x = grid::unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(3, "pt"),
        panel.spacing.x = grid::unit(15, "pt"),
        plot.margin = margin(6, 5, 0, 5)
    ) +
    guides(size = guide_legend(nrow = 1, override.aes = list(colour = "grey30")))
# Add one direction key beneath each facet, aligned to the two sides of zero.
# These are layout grobs and cannot change the plotted data or scale limits.
measurement_pdf <- tempfile(fileext = ".pdf")
grDevices::quartz(type = "pdf", file = measurement_pdf, width = 580/72, height = 352/72, family = FIG_FONT)
built <- ggplot_build(plot)
plot_table <- ggplotGrob(plot)
panels <- plot_table$layout[grepl("^panel", plot_table$layout$name), ]
panels <- panels[order(panels$l), ]
stopifnot(nrow(panels) == 2L)
axis_row <- max(plot_table$layout$b[grepl("^axis-b", plot_table$layout$name)])
plot_table <- gtable::gtable_add_rows(plot_table, grid::unit(20, "pt"), pos = axis_row)
for (i in seq_len(nrow(panels))) {
    span <- built$layout$panel_params[[i]]$x.range
    zero <- (0 - span[1]) / diff(span)
    direction_key <- grid::grobTree(
        grid::textGrob("← Higher in WT", x = zero / 2, y = 0.5,
            gp = grid::gpar(fontfamily = FIG_FONT, fontsize = 9, fontface = "bold", col = "black")),
        grid::textGrob("Higher in CS →", x = (1 + zero) / 2, y = 0.5,
            gp = grid::gpar(fontfamily = FIG_FONT, fontsize = 9, fontface = "bold", col = "black"))
    )
    plot_table <- gtable::gtable_add_grob(plot_table, direction_key,
        t = axis_row + 1, l = panels$l[i], r = panels$r[i],
        clip = "off", name = paste0("direction-key-", i))
}
grDevices::dev.off()
unlink(measurement_pdf)
fig_save(plot_table, file.path(figure_dir, "Rebuilt_Output", "FigS2_KEGG_GSEA"),
         width = 580/72, height = 352/72)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
