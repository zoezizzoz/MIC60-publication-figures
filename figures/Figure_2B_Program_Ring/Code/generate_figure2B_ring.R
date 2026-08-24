#!/usr/bin/env Rscript

# Figure 2B paper-aligned ring and Figure 2C STRING views using the project's established
# ggplot/ggraph/geomtextpath/cowplot figure framework. Statistical values and
# category memberships are read from the current all-sample female analysis.

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggforce)
  library(geomtextpath)
  library(cowplot)
  library(scales)
  library(ggrepel)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}

figure_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "figure_style.R"))

publication_dir <- file.path(figure_dir, "Rebuilt_Output")
dir.create(publication_dir, recursive = TRUE, showWarnings = FALSE)

member_path <- file.path(
  figure_dir, "Supporting_Data",
  "Figure2_paper_aligned_enriched_DEGs_all_members.csv"
)
edge_path <- file.path(
  figure_dir, "Supporting_Data", "Figure2_STRING_edges.csv"
)

members <- read.csv(member_path, check.names = FALSE)
cached_edges <- read.csv(edge_path, check.names = FALSE)

theme_levels <- c(
  "Energy metabolism",
  "Mitochondrial respiration & ATP production",
  "Mitochondrial membrane & protein import",
  "NAD(P)-linked redox metabolism"
)
theme_colors <- c(
  # Match the project's GO reference palette: BP orange, MF blue, CC green.
  # Purple is retained for the mixed respiration/ATP program so all four
  # paper-aligned programs remain visually distinct.
  "Energy metabolism" = "#FBB05B",
  "Mitochondrial respiration & ATP production" = "#B39DDB",
  "NAD(P)-linked redox metabolism" = "#7BC4E2",
  "Mitochondrial membrane & protein import" = "#ACD372"
)

stopifnot(all(theme_levels %in% unique(members$theme)))

# A node has one physical location in a ring/network. The primary assignments
# below balance the four visual sectors; every GO membership remains in the
# membership column and shared nodes receive a black outline.
primary_groups <- list(
  "Energy metabolism" = c(
    "CG7335", "Gpdh2", "CG7755", "Pglym87", "mAcon2", "CG32026"
  ),
  "Mitochondrial respiration & ATP production" = c(
    "Vha68-3", "Hex-t2", "CG3483", "SdhCL", "ND-24L", "SdhAL"
  ),
  "NAD(P)-linked redox metabolism" = c(
    "CG7140", "CG10748", "ND-51L1", "ND-B14.5AL"
  ),
  "Mitochondrial membrane & protein import" = c(
    "CG4701", "CG3092", "Tengl1", "Tengl3", "Tengl4", "ttm3", "CG17991"
  )
)
primary_map <- stack(primary_groups) %>%
  transmute(gene = as.character(values), description = as.character(ind)) %>%
  filter(gene %in% members$gene)

# Retain the established assignments for the original display subset, then
# assign every additional significant member to one of its supported programs.
# Choosing the least-populated eligible program keeps the full ring legible
# without inventing memberships.
membership_lookup <- members %>% distinct(gene, theme)
assignment_counts <- setNames(
  vapply(theme_levels, function(x) sum(primary_map$description == x), integer(1)),
  theme_levels
)
unassigned_genes <- sort(setdiff(unique(members$gene), primary_map$gene))
for (gene_id in unassigned_genes) {
  eligible <- theme_levels[theme_levels %in% membership_lookup$theme[
    membership_lookup$gene == gene_id
  ]]
  if (length(eligible) == 0) stop("No supported paper-aligned group for ", gene_id)
  chosen <- eligible[which.min(assignment_counts[eligible])]
  primary_map <- bind_rows(
    primary_map,
    tibble(gene = gene_id, description = chosen)
  )
  assignment_counts[chosen] <- assignment_counts[chosen] + 1L
}

nodes <- members %>%
  group_by(gene) %>%
  summarise(
    log2FoldChange = first(log2FoldChange),
    padj = first(padj),
    module_memberships = paste(
      theme_levels[theme_levels %in% unique(theme)], collapse = "; "
    ),
    membership_count = n_distinct(theme),
    .groups = "drop"
  ) %>%
  left_join(primary_map, by = "gene") %>%
  mutate(
    description = factor(description, levels = theme_levels),
    significance = -log10(pmax(padj, .Machine$double.xmin))
  )

if (anyNA(nodes$description)) {
  stop("At least one plotted gene lacks a primary paper-aligned group.")
}
if (any(!vapply(seq_len(nrow(nodes)), function(i) {
  as.character(nodes$description[i]) %in% strsplit(
    nodes$module_memberships[i], "; ", fixed = TRUE
  )[[1]]
}, logical(1)))) {
  stop("A primary visual assignment is not supported by its GO membership.")
}

edges <- cached_edges %>%
  transmute(from = from_gene, to = to_gene) %>%
  filter(from %in% nodes$gene, to %in% nodes$gene) %>%
  distinct()

down_col <- unname(FIG_COLORS["control"])
up_col <- unname(FIG_COLORS["mutant"])
fc_limit <- max(2, ceiling(max(abs(nodes$log2FoldChange), na.rm = TRUE)))

fc_fill_scale <- function(horizontal = FALSE) {
  guide <- if (horizontal) {
    guide_colorbar(
      barwidth = 6.5, barheight = 0.75, direction = "horizontal",
      title.position = "top", title.hjust = 0.5, order = 1
    )
  } else {
    guide_colorbar(
      barheight = 4.2, title.position = "top", title.hjust = 0, order = 1
    )
  }
  scale_fill_gradient2(
    low = down_col, mid = "white", high = up_col, midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish,
    name = expression(log[2] ~ FC ~ ("\u2640"~CS/WT)), guide = guide
  )
}

fc_color_scale <- function() {
  scale_color_gradient2(
    low = down_col, mid = "white", high = up_col, midpoint = 0,
    limits = c(-fc_limit, fc_limit), oob = squish,
    name = expression(log[2] ~ FC ~ ("\u2640"~CS/WT)),
    guide = guide_colorbar(
      barheight = 4.2, title.position = "top", title.hjust = 0, order = 1
    )
  )
}

significance_scale <- function() {
  scale_size_continuous(
    range = c(3.2, 7.2),
    breaks = c(1.5, 2, 3),
    labels = c("0.03", "0.01", "0.001"),
    name = "DESeq2 adjusted P",
    guide = guide_legend(title.position = "top", title.hjust = 0.5, order = 2)
  )
}

caption_text <- paste0(
  "Female Drosophila bulk RNA-seq: WT n = 3 and CS n = 3 pooled biological libraries; ",
  "five 96-hour pupae per pool.\n",
  "All genes with DESeq2 padj < 0.05 and membership in a significantly enriched, ",
  "functionally grouped mitochondrial GO program are displayed."
)

# ggplot2 text sizes are specified in millimetres; 2.845276 points = 1 mm.
program_label_size <- FIG_ANNOT_SIZE + 0.65 + 5 / 2.845276

# -----------------------------------------------------------------------------
# Ring view — project-native geometry and typography.
# -----------------------------------------------------------------------------
ring_nodes <- nodes %>%
  arrange(description, desc(abs(log2FoldChange)), gene) %>%
  mutate(
    idx = row_number() - 1L,
    # Rotate the ring one quarter-turn counterclockwise so program positions
    # match the corresponding quadrants in the Figure 2D STRING network.
    xnode = seq(0, 2 * pi, length.out = n() + 1)[-1] - pi / 2
  )

node_r <- 0.56
label_gap <- 0.045
band_text_gap <- label_gap
label_r <- node_r + label_gap
gene_size <- FIG_ANNOT_SIZE - 0.5
text_reach <- 0.038 * max(nchar(ring_nodes$gene))
band_inner <- label_r + text_reach + band_text_gap
band_outer <- band_inner + 0.065
band_label_r <- band_outer + 0.105

ring_nodes <- ring_nodes %>% mutate(
  y = node_r * cos(xnode), x = node_r * sin(xnode),
  angle = atan2(cos(xnode), sin(xnode)) * 180 / pi,
  flip = angle > 90 | angle < -90,
  label_angle = ifelse(flip, angle + 180, angle),
  label_hjust = ifelse(flip, 1, 0),
  label_x = label_r * sin(xnode), label_y = label_r * cos(xnode),
  special_label_x = (label_r - 0.015) * sin(xnode),
  special_label_y = (label_r - 0.015) * cos(xnode)
)

ring_links <- edges %>%
  left_join(ring_nodes %>% select(gene, x, y), by = c("from" = "gene")) %>%
  rename(xfrom = x, yfrom = y) %>%
  left_join(ring_nodes %>% select(gene, x, y), by = c("to" = "gene")) %>%
  rename(xto = x, yto = y) %>%
  filter(if_all(c(xfrom, yfrom, xto, yto), is.finite))

node_step <- 2 * pi / nrow(ring_nodes)
band_mid <- (band_inner + band_outer) / 2
arc_gap <- min(0.013 / band_mid, node_step * 0.48)
ring_rect0 <- ring_nodes %>%
  group_by(description) %>%
  summarise(
    lo = min(xnode) - node_step / 2 + arc_gap,
    hi = max(xnode) + node_step / 2 - arc_gap,
    xmid = (min(xnode) + max(xnode)) / 2,
    .groups = "drop"
  )
ring_rect <- ring_rect0 %>% rowwise() %>% do({
  z <- .
  if (z$hi > 2 * pi) {
    tibble(description = z$description,
           xmin = c(z$lo, 0), xmax = c(2 * pi, z$hi - 2 * pi), xmid = z$xmid)
  } else if (z$lo < 0) {
    tibble(description = z$description,
           xmin = c(0, z$lo + 2 * pi), xmax = c(z$hi, 2 * pi), xmid = z$xmid)
  } else {
    tibble(description = z$description, xmin = z$lo, xmax = z$hi, xmid = z$xmid)
  }
}) %>% ungroup()
ring_labels <- ring_rect0 %>%
  transmute(description, xmid = xmid %% (2 * pi))

inner_limit <- label_r + text_reach + 0.10
ring_inner <- ggplot(ring_nodes, aes(x, y)) +
  geom_segment(
    data = ring_links,
    aes(x = xfrom, xend = xto, y = yfrom, yend = yto),
    inherit.aes = FALSE, linewidth = 0.38, color = "grey68", alpha = 0.24
  ) +
  geom_point(
    aes(fill = log2FoldChange, size = significance),
    shape = 21, color = "grey25", stroke = 0.55
  ) +
  geom_point(
    data = subset(ring_nodes, membership_count > 1),
    aes(size = significance), shape = 21, fill = NA,
    color = "black", stroke = 0.75, show.legend = FALSE
  ) +
  fc_fill_scale(horizontal = TRUE) +
  significance_scale() +
  geom_text(
    data = subset(ring_nodes, gene != "ND-B14.5AL"),
    aes(label = gene, x = label_x, y = label_y,
        angle = label_angle, hjust = label_hjust),
    size = gene_size, fontface = "italic", family = FIG_FONT
  ) +
  geom_text(
    data = subset(ring_nodes, gene == "ND-B14.5AL"),
    aes(label = gene, x = special_label_x, y = special_label_y,
        angle = label_angle, hjust = label_hjust),
    size = gene_size - 0.40, fontface = "italic", family = FIG_FONT
  ) +
  coord_fixed() +
  scale_x_continuous(limits = c(-inner_limit, inner_limit), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-inner_limit, inner_limit), expand = c(0, 0)) +
  theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT) +
  theme(
    legend.position = c(0.5, 0.5),
    legend.box = "vertical",
    legend.title = element_text(size = FIG_LEGEND_TEXT_SIZE),
    legend.text = element_text(size = FIG_LEGEND_TEXT_SIZE),
    legend.spacing.y = unit(2, "pt")
  )

outer_limit <- band_label_r + 0.08
ring_outer <- ggplot() +
  geom_rect(
    data = ring_rect,
    aes(xmin = xmin, xmax = xmax, ymin = band_inner, ymax = band_outer,
        fill = description)
  ) +
  geom_textpath(
    data = ring_labels,
    aes(x = xmid, y = band_label_r, label = description),
    linetype = 0, size = program_label_size, fontface = "bold",
    upright = TRUE, family = FIG_FONT, color = "black"
  ) +
  scale_fill_manual(values = theme_colors, guide = "none") +
  coord_polar() +
  scale_x_continuous(limits = c(0, 2 * pi), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-outer_limit, outer_limit), expand = c(0, 0)) +
  theme_void(base_family = FIG_FONT)

ring <- ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  draw_plot(ring_outer, x = 0.015, y = 0.015, width = 0.97, height = 0.97) +
  draw_plot(ring_inner, x = 0.015, y = 0.015, width = 0.97, height = 0.97)

# -----------------------------------------------------------------------------
# STRING view — project-native grouped ellipse layout.
# -----------------------------------------------------------------------------
graph <- as_tbl_graph(graph_from_data_frame(
  edges, vertices = nodes %>% rename(name = gene), directed = FALSE
)) %>% mutate(description = factor(description, levels = theme_levels))

module_centres <- tibble(
  description = factor(theme_levels, levels = theme_levels),
  cx = c(-7.8, 7.8, -7.8, 7.8),
  cy = c(7.0, 7.0, -7.0, -7.0)
)
layout <- create_layout(graph, layout = "circle") %>%
  left_join(module_centres, by = "description")

for (grp in theme_levels) {
  idx <- which(layout$description == grp)
  # Use the same centered node-ring radius in every program field.
  local_radius <- 3.60
  # Start every ring with one node at 12 o'clock. This prevents two labels
  # from straddling the narrow top-center space while keeping uniform angles.
  theta <- pi / 2 + 2 * pi * (seq_along(idx) - 1) / length(idx)
  layout$x[idx] <- layout$cx[idx] + local_radius * cos(theta)
  layout$y[idx] <- layout$cy[idx] + local_radius * sin(theta)
}

# Place every label radially outside its node. Using explicit positions avoids
# ggrepel occasionally allowing text to graze a large significance-scaled node.
layout <- layout %>% mutate(
  label_dx = x - cx,
  label_dy = y - cy,
  label_distance = pmax(sqrt(label_dx^2 + label_dy^2), 1e-6),
  horizontal_label = abs(label_dx) >= abs(label_dy),
  # Match ggplot2's area-based 3.2--7.2 mm size scale exactly. The panel uses
  # approximately 0.172 data units per mm at the publication export size.
  significance_unit = rescale(
    significance,
    to = c(0, 1),
    from = range(significance, na.rm = TRUE)
  ),
  node_size_mm = 3.2 + (7.2 - 3.2) * sqrt(significance_unit),
  node_radius_data = 0.5 * node_size_mm * 0.172,
  # Use a visually calibrated four-pixel clearance at the normal displayed
  # figure size. A literal four pixels in the 600-DPI export is effectively
  # lost when the publication figure is scaled for viewing.
  label_gap_data = 0.16,
  # Italic vertical/diagonal labels extend slightly inward from their ggplot
  # anchor; compensate for that font metric so the visible ink keeps the gap.
  font_anchor_correction = ifelse(horizontal_label, 0, 0.10),
  label_offset = node_radius_data + label_gap_data + font_anchor_correction,
  label_x = x + label_offset * label_dx / label_distance + case_when(
    name %in% c("CG12229", "CG17300", "CG10962", "CG1724") ~ -0.20,
    name %in% c("mAcon2", "Vha68-3", "UQCR-14L", "ttm3") ~ 0.20,
    TRUE ~ 0
  ),
  label_y = y + label_offset * label_dy / label_distance + case_when(
    name %in% c("ATPsynbetaL", "ATP8", "CG10748", "CG15458") ~ 0.16,
    name %in% c(
      "CG12229", "mAcon2", "CG17300", "Vha68-3",
      "CG10962", "UQCR-14L", "CG1724", "ttm3"
    ) ~ 0.02,
    TRUE ~ 0
  ),
  label_hjust = case_when(
    horizontal_label & label_dx >= 0 ~ 0,
    horizontal_label & label_dx < 0 ~ 1,
    TRUE ~ 0.5
  ),
  label_vjust = case_when(
    !horizontal_label & label_dy >= 0 ~ 0,
    !horizontal_label & label_dy < 0 ~ 1,
    TRUE ~ 0.5
  )
)

# Increase every program-field radius uniformly. At the intended assembled
# panel size this is approximately two display pixels and places the full
# outward-facing gene labels inside their colored program fields.
program_radius_padding_data <- 0.70
ellipses <- layout %>%
  group_by(description) %>%
  summarise(
    cx = mean(x), cy = mean(y),
    # Keep every program field identically sized and centered even though the
    # four programs contain slightly different numbers of genes.
    # Give long gene labels enough interior margin while keeping the four
    # identically sized program fields separated.
    w = 13.6 + 2 * program_radius_padding_data,
    h = 13.6 + 2 * program_radius_padding_data,
    .groups = "drop"
  )
module_labels <- ellipses %>% mutate(
  lx = cx,
  ly = ifelse(cy >= 0, cy + h / 2 + 0.9, cy - h / 2 - 0.9),
  vjust = ifelse(cy >= 0, 0, 1)
)

network <- ggraph(layout) +
  geom_ellipse(
    data = ellipses,
    aes(x0 = cx, y0 = cy, a = w / 2, b = h / 2,
        fill = description, angle = 0),
    alpha = 0.18, color = NA
  ) +
  scale_fill_manual(values = theme_colors, guide = "none") +
  geom_edge_link(color = "grey67", alpha = 0.28, linewidth = 0.42) +
  geom_node_point(
    aes(color = log2FoldChange, size = significance)
  ) +
  geom_node_point(
    aes(
      filter = membership_count > 1,
      size = significance,
      shape = "Multiple-program membership"
    ),
    fill = NA, color = "black", stroke = 0.75
  ) +
  fc_color_scale() +
  significance_scale() +
  scale_shape_manual(
    values = c("Multiple-program membership" = 21),
    name = NULL,
    guide = guide_legend(
      override.aes = list(size = 4.2, fill = NA, color = "black", stroke = 0.75),
      order = 3
    )
  ) +
  geom_text(
    data = layout,
    aes(x = label_x, y = label_y, label = name,
        hjust = label_hjust, vjust = label_vjust),
    inherit.aes = FALSE,
    size = FIG_ANNOT_SIZE - 0.55, fontface = "italic",
    family = FIG_FONT, color = "black"
  ) +
  geom_text(
    data = module_labels,
    aes(x = lx, y = ly, label = description, vjust = vjust),
    inherit.aes = FALSE, hjust = 0.5,
    fontface = "bold", family = FIG_FONT,
    size = program_label_size, color = "black"
  ) +
  coord_fixed(clip = "off") +
  expand_limits(x = c(-16.5, 16.5), y = c(-15.2, 15.2)) +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme_void(base_size = FIG_BASE_SIZE, base_family = FIG_FONT) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = FIG_LEGEND_TEXT_SIZE),
    legend.text = element_text(size = FIG_LEGEND_TEXT_SIZE),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(
      t = FIG_MARGIN_PT[["t"]], r = FIG_MARGIN_PT[["r"]],
      b = FIG_MARGIN_PT[["b"]], l = FIG_MARGIN_PT[["l"]]
    )
  )

# Preserve the handmade drafts once, then replace both the paper-aligned and
# canonical outputs with the project-native figures.
for (stem in c("Figure2_paper_aligned_STRING_ring", "Figure2_paper_aligned_STRING_network")) {
  for (extension in c("png", "pdf")) {
    current <- file.path(publication_dir, paste0(stem, ".", extension))
    backup <- file.path(publication_dir, paste0(stem, "_python_draft.", extension))
    if (file.exists(current) && !file.exists(backup)) file.copy(current, backup)
    cg_backup <- file.path(publication_dir, paste0(stem, "_with_CG_genes.", extension))
    if (file.exists(current) && !file.exists(cg_backup)) file.copy(current, cg_backup)
    no_cg_backup <- file.path(publication_dir, paste0(stem, "_without_CG_genes.", extension))
    if (file.exists(current) && !file.exists(no_cg_backup)) file.copy(current, no_cg_backup)
    top10_backup <- file.path(publication_dir, paste0(stem, "_top10_per_program.", extension))
    if (file.exists(current) && !file.exists(top10_backup)) file.copy(current, top10_backup)
  }
}

fig_save(
  ring, file.path(publication_dir, "Figure2_paper_aligned_STRING_ring"),
  width = 9.0, height = 9.0
)
fig_save(
  network, file.path(publication_dir, "Figure2_paper_aligned_STRING_network"),
  width = 9.2, height = 7.4
)

file.copy(
  file.path(publication_dir, "Figure2_paper_aligned_STRING_ring.png"),
  file.path(publication_dir, "Figure_2B_Program_Ring.png"), overwrite = TRUE
)
file.copy(
  file.path(publication_dir, "Figure2_paper_aligned_STRING_ring.pdf"),
  file.path(publication_dir, "Figure_2B_Program_Ring.pdf"), overwrite = TRUE
)
file.copy(
  file.path(publication_dir, "Figure2_paper_aligned_STRING_network.png"),
  file.path(publication_dir, "Figure_2C_STRING_Network.png"), overwrite = TRUE
)
file.copy(
  file.path(publication_dir, "Figure2_paper_aligned_STRING_network.pdf"),
  file.path(publication_dir, "Figure_2C_STRING_Network.pdf"), overwrite = TRUE
)

node_export <- nodes %>%
  mutate(
    comparison = "female dMIC60-CS versus female dMIC60-WT",
    sample_size = "WT n=3; CS n=3; five 96-hour pupae per pooled library",
    selection = paste0(
      "All genes with DESeq2 padj < 0.05 and functionally grouped mitochondrial GO membership; ",
      "no per-program display cap"
    )
  )
write.csv(
  node_export,
  file.path(publication_dir, "Figure2_paper_aligned_STRING_nodes.csv"),
  row.names = FALSE
)
write.csv(
  edges,
  file.path(publication_dir, "Figure2_paper_aligned_STRING_edges.csv"),
  row.names = FALSE
)

cat(
  "Saved project-native paper-aligned ring and STRING network with ",
  nrow(nodes), " unique genes and ", nrow(edges), " cached edges.\n",
  sep = ""
)
