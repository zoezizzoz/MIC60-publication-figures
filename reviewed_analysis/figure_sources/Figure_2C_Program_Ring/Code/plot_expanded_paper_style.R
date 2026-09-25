# Portable rebuild. Outputs are isolated from final artwork and input snapshots.
.rebuild_args <- grep("^--file=", gsub("~+~", " ", commandArgs(FALSE), fixed=TRUE), value = TRUE)
.rebuild_panel <- dirname(dirname(normalizePath(sub("^--file=", "", .rebuild_args[1]))))
.rebuild_dir <- file.path(.rebuild_panel, "Rebuilt_Output")
dir.create(.rebuild_dir, recursive = TRUE, showWarnings = FALSE)
.rebuild_file <- function(path) file.path(.rebuild_dir, basename(path))
suppressPackageStartupMessages({
    library(tidyverse)
    library(cowplot)
    library(geomtextpath)
    library(ggrepel)
    library(scales)
    library(grid)
})
args_all <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed=TRUE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) dirname(normalizePath(sub("^--file=", "", file_arg))) else getwd()
figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
support_dir <- file.path(figure_dir, "Supporting_Data")
v1_dir <- file.path(figure_dir, "Rebuilt_Output")
set.seed(20260810)
font_family <- ""
up_color <- "#881840"
down_color <- "#455A9F"
program_colors <- c(`ATP synthesis / proton transport` = "#FBB05B", `Complex I / NADH dehydrogenase` = "#B39DDB", `Mitochondrial membrane / transport` = "#ACD372", 
    `Redox / electron transfer` = "#7BC4E2", `Peroxisome / microbody` = "#D8A35D", `Immune / host defense` = "#E79AAA", `Heat / UV response` = "#BE9BD3", 
    `Mitochondrial ribosome` = "#79B9AA", `Protease regulation` = "#F0A667", `Hydrolase / lysozyme` = "#91AFD7", `Extracellular space` = "#B7B7B7")
up_members <- read.csv(file.path(support_dir, "V1_up_28_GO_gene_memberships.csv"), check.names = FALSE)
up_edges <- read.csv(file.path(support_dir, "V1_up_STRING_edges.csv"), check.names = FALSE)
up_terms <- read.csv(file.path(support_dir, "GO_up_28_BH_significant_terms.csv"), check.names = FALSE)
save_both <- function(plot, stem, width, height, dpi = 600) {
    ggsave(paste0(stem, ".png"), plot, width = width, height = height, units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
    ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, units = "in", device = grDevices::pdf, bg = "white", 
        limitsize = FALSE)
}
make_pathway_key <- function(terms, heading, text_size = 2.8) {
    x <- terms %>% arrange(factor(collapsed_program, levels = names(program_colors)), p.adjust) %>% mutate(row = rev(seq_len(n())), 
        display = paste0(Description, "  [", ID, "]"))
    ggplot(x, aes(y = row)) + geom_point(aes(x = 0, color = collapsed_program), size = 2) + geom_text(aes(x = 0.035, label = display, 
        color = collapsed_program), hjust = 0, size = text_size, family = font_family) + scale_color_manual(values = program_colors, 
        guide = "none") + scale_x_continuous(limits = c(-0.02, 1), expand = c(0, 0)) + scale_y_continuous(limits = c(0, nrow(x) + 
        2), expand = c(0, 0)) + labs(title = heading, subtitle = "Every listed term has BH-adjusted P < 0.05; colors match the grouped display fields") + 
        theme_void(base_family = font_family) + theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 9, 
        color = "grey35"), plot.margin = margin(25, 15, 25, 10))
}
summarise_nodes <- function(members, program_levels) {
    membership <- members %>% distinct(gene, collapsed_program)
    genes <- members %>% group_by(gene) %>% summarise(log2FoldChange = first(log2FoldChange), padj = first(padj), membership_count = n_distinct(collapsed_program), 
        memberships = paste(program_levels[program_levels %in% unique(collapsed_program)], collapse = "; "), .groups = "drop") %>% 
        arrange(desc(membership_count), gene)
    counts <- setNames(rep(0L, length(program_levels)), program_levels)
    assigned <- character(nrow(genes))
    for (i in seq_len(nrow(genes))) {
        eligible <- program_levels[program_levels %in% membership$collapsed_program[membership$gene == genes$gene[i]]]
        chosen <- eligible[which.min(counts[eligible])]
        assigned[i] <- chosen
        counts[chosen] <- counts[chosen] + 1L
    }
    genes %>% mutate(program = factor(assigned, levels = program_levels), significance = -log10(pmax(padj, .Machine$double.xmin)))
}
edge_table <- function(edges) edges %>% transmute(from = from_gene, to = to_gene, combined_score = combined_score) %>% distinct()
fc_limit_for <- function(nodes) max(2, ceiling(max(abs(nodes$log2FoldChange), na.rm = TRUE)))
make_ring <- function(nodes, edges, program_levels, program_palette, gene_size = 3, program_size = 5, node_r = 0.56, label_gap = 0.045, 
    bar_clearance = 0.045, legend_background = "transparent", legend_outside = FALSE, radial_compression = 1, link_color = "grey68", 
    link_alpha = 0.24, link_width = 0.38, legend_title_size = 9.5, legend_text_size = 9.5, fc_legend_title = "log2 FC (female dMIC60-CS/dMIC60-WT)", 
    legend_bar_width = 6.5, legend_bar_height = 0.7) {
    nodes <- nodes %>% arrange(program, desc(abs(log2FoldChange)), gene) %>% mutate(idx = row_number() - 1L, xnode = seq(0, 
        2 * pi, length.out = n() + 1)[-1] - pi/2)
    label_r <- node_r + label_gap
    text_reach <- 0.038 * max(nchar(nodes$gene)) * gene_size/2.7
    band_inner <- label_r + text_reach + bar_clearance
    band_outer <- band_inner + 0.065
    band_label_r <- band_outer + 0.105
    nodes <- nodes %>% mutate(y = node_r * cos(xnode), x = node_r * sin(xnode), angle = atan2(cos(xnode), sin(xnode)) * 180/pi, 
        flip = angle > 90 | angle < -90, label_angle = ifelse(flip, angle + 180, angle), label_hjust = ifelse(flip, 1, 0), 
        label_x = label_r * sin(xnode), label_y = label_r * cos(xnode))
    links <- edge_table(edges) %>% left_join(nodes %>% select(gene, x, y), by = c(from = "gene")) %>% rename(xfrom = x, yfrom = y) %>% 
        left_join(nodes %>% select(gene, x, y), by = c(to = "gene")) %>% rename(xto = x, yto = y) %>% filter(if_all(c(xfrom, 
        yfrom, xto, yto), is.finite))
    node_step <- 2 * pi/nrow(nodes)
    band_mid <- (band_inner + band_outer)/2
    arc_gap <- min(0.013/band_mid, node_step * 0.46)
    arc0 <- nodes %>% group_by(program) %>% summarise(lo = min(xnode) - node_step/2 + arc_gap, hi = max(xnode) + node_step/2 - 
        arc_gap, xmid = (min(xnode) + max(xnode))/2, .groups = "drop")
    arcs <- arc0 %>% rowwise() %>% do({
        z <- .
        if (z$hi > 2 * pi) {
            tibble(program = z$program, xmin = c(z$lo, 0), xmax = c(2 * pi, z$hi - 2 * pi))
        }
        else if (z$lo < 0) {
            tibble(program = z$program, xmin = c(0, z$lo + 2 * pi), xmax = c(z$hi, 2 * pi))
        }
        else tibble(program = z$program, xmin = z$lo, xmax = z$hi)
    }) %>% ungroup()
    labels <- arc0 %>% transmute(program, xmid = xmid%%(2 * pi), display_label = ring_program_label(program))
    fc_limit <- fc_limit_for(nodes)
    inner_limit <- (label_r + text_reach + 0.1) * radial_compression
    inner <- ggplot(nodes, aes(x, y)) + geom_segment(data = links, aes(x = xfrom, xend = xto, y = yfrom, yend = yto), inherit.aes = FALSE, 
        linewidth = link_width, color = link_color, alpha = link_alpha) + geom_point(aes(fill = log2FoldChange, size = significance), 
        shape = 21, color = "grey25", stroke = 0.55) + geom_point(data = subset(nodes, membership_count > 1), aes(size = significance), 
        shape = 21, fill = NA, color = "black", stroke = 0.8, show.legend = FALSE) + scale_fill_gradientn(colours = c(down_color, 
        "#FFFFFF", up_color), values = c(0, 0.5, 1), limits = c(-fc_limit, fc_limit), oob = squish, 
        name = fc_legend_title, guide = guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = 0.5, 
            barwidth = legend_bar_width, barheight = legend_bar_height)) + scale_size_continuous(range = c(3.2, 7.2), name = expression("DESeq2 adjusted "*italic(p)), 
        breaks = c(1.5, 2, 3), labels = c("0.03", "0.01", "0.001")) + geom_text(aes(label = gene, x = label_x, y = label_y, 
        angle = label_angle, hjust = label_hjust), size = gene_size, fontface = "italic", family = font_family) + coord_fixed() + 
        scale_x_continuous(limits = c(-inner_limit, inner_limit), expand = c(0, 0)) + scale_y_continuous(limits = c(-inner_limit, 
        inner_limit), expand = c(0, 0)) + theme_void(base_family = font_family) + theme(legend.position = c(0.5, 0.5), legend.justification = c(0.5, 
        0.5), legend.box = "vertical", legend.box.just = "center", legend.title = element_text(size = legend_title_size), 
        legend.text = element_text(size = legend_text_size), legend.background = element_rect(fill = legend_background, color = NA), 
        legend.box.background = element_rect(fill = legend_background, color = NA), legend.box.margin = margin(1, 2, 1, 2), 
        legend.spacing.y = unit(0.05, "cm"))
    outer_limit <- (band_label_r + 0.1) * radial_compression
    outer <- ggplot() + geom_rect(data = arcs, aes(xmin = xmin, xmax = xmax, ymin = band_inner, ymax = band_outer, fill = program)) + 
        geom_textpath(data = labels, aes(x = xmid, y = band_label_r, label = display_label), linetype = 0, size = program_size, 
            fontface = "bold", upright = TRUE, family = font_family, color = "black") + scale_fill_manual(values = program_palette, 
        guide = "none") + coord_polar() + scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0)) + scale_y_continuous(limits = c(-outer_limit, 
        outer_limit), expand = c(0, 0)) + theme_void(base_family = font_family)
    if (legend_outside) {
        key <- cowplot::get_legend(inner + theme(legend.position = "right", legend.justification = "center", legend.background = element_rect(fill = "white", 
            color = NA), legend.box.background = element_rect(fill = "white", color = NA)))
        ring_panel <- ggdraw() + draw_plot(outer, 0.015, 0.015, 0.97, 0.97) + draw_plot(inner + theme(legend.position = "none"), 
            0.015, 0.015, 0.97, 0.97)
        return(ggdraw() + draw_plot(ring_panel, x = 0, y = 0, width = 0.726, height = 1) + draw_grob(key, x = 0.615, y = 0.23, 
            width = 0.3, height = 0.54))
    }
    ggdraw() + draw_plot(outer, 0.015, 0.015, 0.97, 0.97) + draw_plot(inner, 0.015, 0.015, 0.97, 0.97)
}
program_centres <- function(program_levels, top_row_y = 6, second_top_y = top_row_y, top_x = c(-6, 6), bottom_x = c(-10, 
    0, 10)) {
    k <- length(program_levels)
    if (k == 5) {
        xy <- tibble(cx = c(top_x, bottom_x), cy = c(top_row_y, second_top_y, -6, -6, -6))
    }
    else if (k == 6) {
        xy <- tibble(cx = rep(c(-10, 0, 10), 2), cy = rep(c(6, -6), each = 3))
    }
    else {
        theta <- seq(pi/2, pi/2 + 2 * pi, length.out = k + 1)[-(k + 1)]
        xy <- tibble(cx = 9 * cos(theta), cy = 7 * sin(theta))
    }
    bind_cols(tibble(program = factor(program_levels, levels = program_levels)), xy)
}
wrap_program <- function(x) {
    dplyr::recode(as.character(x), `ATP synthesis / proton transport` = "ATP synthesis /\nproton transport", `Complex I / NADH dehydrogenase` = "Complex I /\nNADH dehydrogenase", 
        `Mitochondrial membrane / transport` = "Mitochondrial membrane /\ntransport", `Redox / electron transfer` = "Redox /\nelectron transfer", 
        `Peroxisome / microbody` = "Peroxisome /\nmicrobody", `Immune / host defense` = "Immune /\nhost defense", `Heat / UV response` = "Heat / UV\nresponse", 
        `Mitochondrial ribosome` = "Mitochondrial\nribosome", `Protease regulation` = "Protease\nregulation", `Hydrolase / lysozyme` = "Hydrolase /\nlysozyme", 
        `Extracellular space` = "Extracellular\nspace")
}
ring_program_label <- function(x) {
    dplyr::recode(as.character(x), `Energy metabolism` = "Energy metabolism", `Mitochondrial respiration & ATP production` = "Mitochondrial respiration & ATP production", 
        `Mitochondrial membrane & protein import` = "Mitochondrial membrane & protein import", `NAD(P)-linked redox metabolism` = "NAD(P)-linked redox metabolism", 
        `ATP synthesis / proton transport` = "ATP synthesis", `Complex I / NADH dehydrogenase` = "Complex I", `Mitochondrial membrane / transport` = "Mitochondrial membrane", 
        `Redox / electron transfer` = "Redox / electron transfer", `Peroxisome / microbody` = "Peroxisome / microbody", `Immune / host defense` = "Immune defense", 
        `Heat / UV response` = "Heat / UV", `Mitochondrial ribosome` = "Mito ribosome", `Protease regulation` = "Protease", 
        `Hydrolase / lysozyme` = "Hydrolase", `Extracellular space` = "Extracellular")
}
grouped_layout <- function(nodes, program_levels, top_row_y = 6, second_top_y = top_row_y, top_x = c(-6, 6), bottom_x = c(-10, 
    0, 10), center_singletons = FALSE) {
    centres <- program_centres(program_levels, top_row_y = top_row_y, second_top_y = second_top_y, top_x = top_x, bottom_x = bottom_x)
    out <- nodes %>% left_join(centres, by = c(program = "program")) %>% mutate(x = NA_real_, y = NA_real_)
    for (p in program_levels) {
        idx <- which(out$program == p)
        if (length(idx) == 1 && center_singletons) {
            out$x[idx] <- out$cx[idx]
            out$y[idx] <- out$cy[idx]
            next
        }
        theta <- pi/2 + 2 * pi * (seq_along(idx) - 1)/length(idx)
        radius <- case_when(length(idx) <= 3 ~ 2.3, length(idx) >= 12 ~ 4.35, length(idx) >= 8 ~ 4.1, TRUE ~ 3.35)
        out$x[idx] <- out$cx[idx] + radius * cos(theta)
        out$y[idx] <- out$cy[idx] + radius * sin(theta)
    }
    out %>% mutate(dx = x - cx, dy = y - cy, dist = pmax(sqrt(dx^2 + dy^2), 1e-06), label_x = x + 0.58 * dx/dist, label_y = y + 
        0.58 * dy/dist, label_hjust = ifelse(abs(dx) < 0.3, 0.5, ifelse(dx > 0, 0, 1)), label_vjust = ifelse(abs(dx) < 0.3, 
        ifelse(dy > 0, 0, 1), 0.5))
}
make_grouped_network <- function(nodes, edges, program_levels, program_palette, title = NULL, context = NULL, context_edges = NULL, 
    gene_size = 2.65, program_label_size = 4.8, legend_title_size = 9.5, legend_text_size = 9.5, field_alpha = 0.18, link_color = "grey66", 
    link_alpha_range = c(0.16, 0.58), node_outline = "black", compact_small_fields = FALSE, top_row_y = 6, second_top_y = top_row_y, 
    top_x = c(-6, 6), bottom_x = c(-10, 0, 10), center_singletons = FALSE, manual_label_offsets = NULL) {
    layout <- grouped_layout(nodes, program_levels, top_row_y = top_row_y, second_top_y = second_top_y, top_x = top_x, bottom_x = bottom_x, 
        center_singletons = center_singletons)
    program_counts <- layout %>% count(program, name = "node_count")
    ellipses <- program_centres(program_levels, top_row_y = top_row_y, second_top_y = second_top_y, top_x = top_x, bottom_x = bottom_x) %>% 
        left_join(program_counts, by = "program") %>% mutate(diameter = case_when(compact_small_fields & node_count == 1 ~ 
        5.8, compact_small_fields & node_count <= 3 ~ 7, node_count >= 12 ~ 11, node_count >= 8 ~ 10.6, TRUE ~ 9), w = diameter, 
        h = diameter)
    labels <- ellipses %>% mutate(label = wrap_program(program), lx = cx, ly = ifelse(cy > 0, cy + h/2 + 0.55, cy - h/2 - 
        0.55), vjust = ifelse(cy > 0, 0, 1))
    all_edges <- if (is.null(context_edges)) 
        edge_table(edges)
    else context_edges %>% transmute(from = from_gene, to = to_gene, combined_score = combined_score) %>% distinct()
    if (!is.null(context)) {
        seed_positions <- layout %>% select(gene, x, y)
        ctx <- context %>% filter(node_type == "First-shell context")
        ctx$context_program <- NA_character_
        for (i in seq_len(nrow(ctx))) {
            nb <- unique(c(all_edges$to[all_edges$from == ctx$gene[i]], all_edges$from[all_edges$to == ctx$gene[i]]))
            supported <- as.character(layout$program[layout$gene %in% nb])
            if (length(supported)) 
                ctx$context_program[i] <- names(sort(table(supported), decreasing = TRUE))[1]
        }
        ctx$context_program[is.na(ctx$context_program)] <- program_levels[1]
        ctx <- ctx %>% left_join(program_centres(program_levels, top_row_y = top_row_y, second_top_y = second_top_y, top_x = top_x, 
            bottom_x = bottom_x) %>% mutate(context_program = as.character(program)) %>% select(context_program, cx, cy), 
            by = "context_program") %>% mutate(x = NA_real_, y = NA_real_)
        for (pname in unique(ctx$context_program)) {
            idx <- which(ctx$context_program == pname)
            theta <- pi/2 + 2 * pi * (seq_along(idx) - 1)/length(idx)
            context_radius <- ifelse(length(idx) >= 5, 1.75, 1.35)
            ctx$x[idx] <- ctx$cx[idx] + context_radius * cos(theta)
            ctx$y[idx] <- ctx$cy[idx] + context_radius * sin(theta)
        }
        ctx <- ctx %>% mutate(ddx = x - cx, ddy = y - cy, dd = pmax(sqrt(ddx^2 + ddy^2), 1e-06), label_x = x + 0.38 * ddx/dd, 
            label_y = y + 0.38 * ddy/dd)
    }
    else ctx <- tibble(gene = character(), x = double(), y = double(), label_x = double(), label_y = double())
    positions <- bind_rows(layout %>% select(gene, x, y), ctx %>% select(gene, x, y))
    links <- all_edges %>% left_join(positions, by = c(from = "gene")) %>% rename(xfrom = x, yfrom = y) %>% left_join(positions, 
        by = c(to = "gene")) %>% rename(xto = x, yto = y) %>% filter(if_all(c(xfrom, yfrom, xto, yto), is.finite))
    fc_limit <- fc_limit_for(nodes)
    if (is.null(manual_label_offsets)) {
        manual_labels <- tibble(gene = character(), x = double(), y = double(), dx = double(), dy = double())
    }
    else {
        manual_label_offsets <- manual_label_offsets %>% transmute(gene, offset_x = dx, offset_y = dy)
        manual_labels <- layout %>% inner_join(manual_label_offsets, by = "gene") %>% mutate(label_x = x + offset_x, label_y = y + 
            offset_y)
    }
    gene_label_layers <- lapply(seq_len(nrow(ellipses)), function(i) {
        pname <- as.character(ellipses$program[i])
        b <- ellipses[i, ]
        geom_text_repel(data = layout %>% filter(as.character(program) == pname, !gene %in% manual_labels$gene), aes(x = x, 
            y = y, label = gene), size = gene_size, fontface = "italic", family = font_family, box.padding = 0.65, point.padding = 0.5, 
            min.segment.length = 0, segment.color = link_color, segment.size = 0.25, force = 10, force_pull = 0.1, max.overlaps = Inf, 
            max.time = 20, max.iter = 2e+05, seed = 20260810 + i, xlim = c(b$cx - b$w/2 + 0.75, b$cx + b$w/2 - 0.75), ylim = c(b$cy - 
                b$h/2 + 0.75, b$cy + b$h/2 - 0.75))
    })
    p <- ggplot() + ggforce::geom_ellipse(data = ellipses, aes(x0 = cx, y0 = cy, a = w/2, b = h/2, fill = program, angle = 0), 
        alpha = field_alpha, color = NA) + geom_segment(data = links, aes(x = xfrom, xend = xto, y = yfrom, yend = yto, alpha = combined_score), 
        color = link_color, linewidth = 0.45) + geom_text(data = labels, aes(x = lx, y = ly, label = label, vjust = vjust), 
        hjust = 0.5, fontface = "bold", size = program_label_size, lineheight = 0.9, family = font_family) + scale_fill_manual(values = program_palette, 
        guide = "none")
    p <- p + ggnewscale::new_scale_fill() + geom_point(data = layout, aes(x, y, fill = log2FoldChange, size = significance), 
        shape = 21, color = node_outline, stroke = 0.65) + geom_point(data = subset(layout, membership_count > 1), aes(x, 
        y, size = significance), shape = 21, fill = NA, color = "black", stroke = 0.9, show.legend = FALSE) + scale_fill_gradientn(colours = c(down_color, 
        "#FFFFFF", up_color), values = c(0, 0.5, 1), limits = c(-fc_limit, fc_limit), oob = squish, 
        name = "log2 FC (female dMIC60-CS/dMIC60-WT)") + scale_size_continuous(range = c(3.2, 7.2), name = expression("DESeq2 adjusted "*italic(p)), 
        breaks = c(1.5, 2, 3), labels = c("0.03", "0.01", "0.001")) + scale_alpha_continuous(range = link_alpha_range, guide = "none")
    p <- p + gene_label_layers
    if (nrow(manual_labels)) {
        p <- p + geom_text(data = manual_labels, aes(x = label_x, y = label_y, label = gene), size = gene_size, fontface = "italic", 
            family = font_family)
    }
    if (nrow(ctx)) {
        p <- p + geom_point(data = ctx, aes(x, y), shape = 21, size = 4.2, fill = "grey82", color = "grey54", stroke = 0.8) + 
            geom_text_repel(data = ctx, aes(x = x, y = y, label = gene), size = 2.7, color = "grey40", family = font_family, 
                box.padding = 0.65, point.padding = 0.45, min.segment.length = 0, segment.color = "grey68", segment.size = 0.22, 
                force = 8, force_pull = 0.08, max.overlaps = Inf, max.time = 20, max.iter = 180000, seed = 20260811, xlim = c(-15.2, 
                  12.5), ylim = c(-8.9, 8.9))
    }
    p + coord_fixed(clip = "off") + expand_limits(x = c(-16, 16), y = c(-13.5, 13.5)) + labs(title = title) + theme_void(base_family = font_family) + 
        theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5), legend.title = element_text(size = legend_title_size), 
            legend.text = element_text(size = legend_text_size), legend.position = "right", plot.margin = margin(8, 10, 8, 
                8))
}
up_levels <- c("ATP synthesis / proton transport", "Complex I / NADH dehydrogenase", "Mitochondrial membrane / transport", 
    "Redox / electron transfer", "Peroxisome / microbody")
down_levels <- c("Immune / host defense", "Heat / UV response", "Mitochondrial ribosome", "Protease regulation", "Hydrolase / lysozyme", 
    "Extracellular space")
up_nodes <- summarise_nodes(up_members, up_levels)
v1_ring <- make_ring(up_nodes, up_edges, up_levels, program_colors[up_levels], gene_size = 3.75, program_size = 5.6, node_r = 0.56, 
    label_gap = 0.07, bar_clearance = 0.04, legend_background = "transparent", legend_outside = FALSE, radial_compression = 1.13, 
    link_color = "grey50", link_alpha = 0.4, link_width = 0.42, legend_title_size = 8.2, legend_text_size = 7.8, fc_legend_title = "log2 FC (female CS/WT)", 
    legend_bar_width = 4.8, legend_bar_height = 0.55)
save_both(v1_ring, file.path(v1_dir, "Fig2B_program_ring"), 12.4, 9)
writeLines(capture.output(sessionInfo()), .rebuild_file(file.path(.rebuild_dir, "R_sessionInfo.txt")))
