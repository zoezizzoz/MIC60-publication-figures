library(damr)
library(sleepr)
library(ggetho)
library(data.table)
library(ggplot2)

# ggpubr, rstatix and dplyr are no longer needed: the within-phase
# Wilcoxon tests and the brackets are computed here and drawn with the
# shared helpers, so every bracket in the project looks the same.

# ============================================================
# 1. Load data
# ============================================================

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  normalizePath(getwd())
}
FIGURE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
PACKAGE_DIR <- normalizePath(file.path(FIGURE_DIR, ".."), mustWork = TRUE)
DATA_DIR <- file.path(FIGURE_DIR, "Original_Data")
FIG_DIR <- file.path(FIGURE_DIR, "Rebuilt_Output")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Metadata file that also contains the W1118 flies.
metadata <- fread(file.path(DATA_DIR, "sleep_metadata.csv"))
metadata <- link_dam_metadata(metadata, result_dir = DATA_DIR)

dt <- load_dam(metadata, FUN = sleepr::sleep_dam_annotation)
dt_curated <- curate_dead_animals(dt)

# ============================================================
# 2. Keep only female flies, only the three genotypes we plot
#    (WR = WT rescue control, CS = oxidation-resistant rescue,
#     W1118 = background strain). The separate literal "WT"
#     group in the file is intentionally excluded here.
# ============================================================

keep_geno <- c("WR", "CS", "W1118")
dt_female <- dt_curated[xmv(sex) == "Female" & xmv(genotype) %in% keep_geno]

# ============================================================
# 3. Shared style: labels, order, colours, theme
#    Colours, fonts, dot and error bar settings, chart sizes and the
#    "bracket only when significant" rule all live in figure_style.R,
#    the same file used by the TEM, TMRM, MTT and RNA-seq figures.
# ============================================================

source(file.path(SCRIPT_DIR, "figure_style.R"))

# Match the TMRM publication figure's vector export: Helvetica-family Quartz
# text with Unicode support (including the female symbol).
if (identical(Sys.info()[["sysname"]], "Darwin") &&
    isTRUE(capabilities("aqua"))) {
  FIG_PDF_DEVICE <- function(filename, width, height, ...) {
    grDevices::quartz(
      type = "pdf", file = filename, width = width, height = height, ...
    )
  }
}

# WR = dMIC60 wild-type rescue (control, shown as "WT");
# CS = dMIC60 oxidation-resistant mutant rescue (experimental);
# W1118 = background strain.
genotype_labels <- c(
  "WR" = "dMIC60WT \u2640",
  "CS" = "dMIC60CS \u2640"
)
geno_order      <- c("WR", "CS")   # WT (rescue control) always on the left

fly_n_caption <- function(data, phase_col = NULL) {
  d <- unique(as.data.table(data)[genotype %in% geno_order,
                                  c("id", "genotype", phase_col), with = FALSE])
  if (is.null(phase_col)) {
    counts <- d[, .(n = uniqueN(id)), by = genotype]
    n_for <- function(g) counts[as.character(genotype) == g, n]
    return(stringr::str_wrap(
      paste0("Adult female flies monitored individually; WT n = ",
             n_for("WR"), ", CS n = ", n_for("CS"), "."), width = 72
    ))
  }

  counts <- d[, .(n = uniqueN(id)), by = c(phase_col, "genotype")]
  phases <- unique(as.character(d[[phase_col]]))
  phase_text <- vapply(phases, function(ph) {
    n_for <- function(g) counts[get(phase_col) == ph &
                                 as.character(genotype) == g, n]
    paste0(ph, ": WT n = ", n_for("WR"), ", CS n = ", n_for("CS"))
  }, character(1))
  stringr::str_wrap(
    paste0("Adult female flies monitored individually; ",
           paste(phase_text, collapse = "; "), "."), width = 72
  )
}

# Light and dark phase shading for the ZT profiles. These are lighting
# conditions rather than genotypes, so they keep their own colours.
phase_fill <- c(on = "#FFF6C2", off = "#D9D9D9")
phase_bar  <- c(on = "#FFD700", off = "black")

# Helper: remove outliers by genotype (1.5 x IQR)
remove_outliers_iqr <- function(dt, value_col, group_cols) {
  dt[
    ,
    {
      x   <- get(value_col)
      q1  <- quantile(x, 0.25, na.rm = TRUE)
      q3  <- quantile(x, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      .SD[x >= (q1 - 1.5 * iqr) & x <= (q3 + 1.5 * iqr)]
    },
    by = group_cols
  ]
}

# Helper: a consistently formatted boxplot across all genotypes.
# A bracket, stars and the p value are shown only when p < FIG_ALPHA.
two_group_boxplot <- function(data, yvar, title, ylab, ylim = NULL) {
  plot_data <- as.data.table(data)[genotype %in% geno_order]
  plot_data[, genotype := factor(genotype, levels = geno_order)]
  
  group_1 <- plot_data[genotype == geno_order[1]][[yvar]]
  group_2 <- plot_data[genotype == geno_order[2]][[yvar]]
  comparison_p <- wilcox.test(group_1, group_2)$p.value
  
  set.seed(FIG_SEED)
  
  ggplot(plot_data,
         aes(x = genotype, y = .data[[yvar]], fill = genotype)) +
    fig_summary(p_value = comparison_p) +
    fig_points(color = FIG_PT_COLOR) +
    fig_sig_bracket(comparison_p, plot_data[[yvar]]) +
    fig_scale_x_group(genotype_labels) +
    fig_scale_fill(geno_order, labels = genotype_labels, guide = "none") +
    fig_scale_y(ylab, limits = NULL) +
    labs(x = NULL, title = stringr::str_wrap(title, width = 60),
         caption = fly_n_caption(plot_data)) +
    theme_fig() +
    coord_cartesian(ylim = ylim, clip = "off")
}

# Helper: within-phase WR vs CS Wilcoxon test, returned with the x
# positions of the two dodged boxes so a bracket spans the pair it
# actually compares.
phase_pvals <- function(data, yvar, dodge_width = FIG_DODGE_WIDTH) {
  d <- as.data.table(data)[genotype %in% geno_order]
  d[, genotype := factor(genotype, levels = geno_order)]
  
  phases <- if (is.factor(d$phase)) levels(d$phase) else sort(unique(d$phase))
  phases <- phases[phases %in% unique(as.character(d$phase))]
  
  offset <- dodge_width / 4   # two boxes per phase
  
  do.call(rbind, lapply(seq_along(phases), function(i) {
    g1 <- d[as.character(phase) == phases[i] & genotype == geno_order[1]][[yvar]]
    g2 <- d[as.character(phase) == phases[i] & genotype == geno_order[2]][[yvar]]
    
    p_value <- if (length(g1) > 1 && length(g2) > 1) {
      wilcox.test(g1, g2)$p.value
    } else {
      NA_real_
    }
    
    data.frame(
      phase = phases[i],
      x1 = i - offset,
      x2 = i + offset,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }))
}

# Helper: grouped (day vs night) dodged boxplot. The boxplot layer
# explicitly supplies the genotype key. A bracket, stars and the p value
# appear only for phase-specific comparisons with p < FIG_ALPHA.
phase_boxplot <- function(summary_dt, yvar, title, ylab, ylim = NULL) {
  plot_data <- as.data.table(summary_dt)[genotype %in% geno_order]
  plot_data[, genotype := factor(genotype, levels = geno_order)]
  
  comparisons <- phase_pvals(plot_data, yvar)
  
  set.seed(FIG_SEED)
  
  ggplot(
    plot_data,
    aes(x = phase, y = .data[[yvar]], fill = genotype)
  ) +
    fig_summary(
      width = 0.6,
      position = position_dodge(width = FIG_DODGE_WIDTH),
      show_legend = TRUE
    ) +
    geom_point(
      shape = FIG_PT_SHAPE,
      size = FIG_PT_SIZE,
      alpha = FIG_PT_ALPHA,
      color = FIG_PT_COLOR,
      show.legend = FALSE,
      position = position_jitterdodge(
        jitter.width = FIG_JITTER_WIDTH,
        jitter.height = 0,
        dodge.width = FIG_DODGE_WIDTH,
        seed = FIG_SEED
      )
    ) +
    fig_sig_brackets_at(comparisons, plot_data[[yvar]]) +
    fig_scale_fill(
      geno_order,
      labels = genotype_labels,
      guide = guide_legend(title = NULL)
    ) +
    fig_scale_y(ylab, limits = NULL) +
    scale_x_discrete(
      name = NULL,
      expand = expansion(add = FIG_X_EXPAND)
    ) +
    labs(
      title = stringr::str_wrap(title, width = 40),
      fill = NULL,
      caption = fly_n_caption(plot_data, "phase")
    ) +
    theme_fig(legend_position = "top") +
    theme(
      legend.direction = "horizontal",
      legend.justification = "center"
    ) +
    coord_cartesian(ylim = ylim, clip = "off")
}

# Helper: actual consecutive 48-hour ZT profile (sleep OR activity).
# The SEM ribbon is the profile equivalent of the error bars used in the
# categorical panels, so it takes the genotype colour.
build_profile <- function(perfly, ylab, title, legend_position = c(0.5, 0.12)) {
  prof <- perfly[, .(m = mean(val), sem = sd(val) / sqrt(.N)), by = .(genotype, zbin)]
  prof[, genotype := factor(genotype, levels = geno_order)]
  
  ymax   <- max(prof$m + prof$sem, na.rm = TRUE) * 1.05
  bar_lo <- ymax * 1.02
  bar_hi <- ymax * 1.08
  ld <- data.frame(xmin = c(0, 12, 24, 36), xmax = c(12, 24, 36, 48),
                   phase = c("on", "off", "on", "off"))
  
  ggplot() +
    # Light and dark shading, drawn without a fill mapping so the fill
    # scale stays available for the genotype ribbons.
    geom_rect(data = ld, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax),
              fill = phase_fill[ld$phase], alpha = 0.5) +
    geom_ribbon(data = prof,
                aes(x = zbin, ymin = m - sem, ymax = m + sem, fill = genotype),
                alpha = FIG_RIBBON_ALPHA) +
    geom_rect(data = ld, aes(xmin = xmin, xmax = xmax, ymin = bar_lo, ymax = bar_hi),
              fill = phase_bar[ld$phase], colour = "black", linewidth = 0.3) +
    geom_line(data = prof, aes(x = zbin, y = m, colour = genotype),
              linewidth = FIG_PROFILE_LINEWIDTH) +
    scale_fill_manual(values = fig_colors(geno_order), guide = "none") +
    fig_scale_color(geno_order, labels = genotype_labels) +
    scale_x_continuous(name = "ZT (h)", breaks = seq(0, 48, 6),
                       limits = c(0, 48), expand = c(0, 0)) +
    scale_y_continuous(name = ylab, expand = c(0, 0)) +
    coord_cartesian(ylim = c(0, bar_hi), clip = "off") +
    labs(title = stringr::str_wrap(title, width = 60),
         caption = fly_n_caption(perfly)) +
    theme_fig(legend_position = legend_position) +
    theme(legend.direction = "horizontal",
          legend.background = element_rect(fill = scales::alpha("white", 0.6),
                                           colour = NA))
}

# ============================================================
# 4. Build a plain data.table, trim endpoint, add ZT + phases
# ============================================================

dt_plot <- as.data.table(dt_female)
meta_dt <- as.data.table(meta(dt_female))
if (!"genotype" %in% names(dt_plot)) {
  dt_plot <- merge(dt_plot, meta_dt[, .(id, genotype)], by = "id", all.x = TRUE)
}
dt_plot[, genotype := factor(genotype, levels = geno_order)]

dt_plot[, hour := t / 3600]
dt_plot <- dt_plot[hour < 60]                 # drop jumpy endpoint

# Full 60-hour dataset used for all summary plots and their statistics.
# Experiment starts 22:00 -> lights ON at t = 12 h, so ZT0 = lights on.
dt_plot[, ZT   := ((hour - 12) %% 24)]
dt_plot[, zbin := floor(ZT * 2) / 2]          # 30-min bins

# lights-off (night) = ZT 12-24 ; lights-on (day) = ZT 0-12
dt_plot[, phase := ifelse(ZT >= 12, "Night (lights off)", "Day (lights on)")]

night_dt <- dt_plot[ZT >= 12]
day_dt   <- dt_plot[ZT <  12]

# Separate consecutive 48-hour dataset used only for the time-series profiles.
# Keeping this separate ensures that all bar/box summary plots still use 60 h.
dt_profile <- copy(dt_plot[hour >= 12 & hour < 60])
dt_profile[, zbin := floor((hour - 12) * 2) / 2]

# ============================================================
#  SLEEP FIGURES
# ============================================================

# ---- FIGURE 1 - Sleep profile (min / 30 min) across ZT ----------------------
perfly_sleep <- dt_profile[genotype %in% geno_order,
                           .(val = mean(asleep, na.rm = TRUE) * 30),   # min per 30-min bin
                           by = .(genotype, id, zbin)]

fig1 <- build_profile(perfly_sleep,
                      ylab  = "Sleep (min / 30 min)",
                      title = "CS females shift sleep toward the night")
fig_save(fig1, file.path(FIG_DIR, "Fig1_female_sleep_profile"),
         width = FIG_W_2COL, height = FIG_H_WIDE)

# ---- FIGURE 2 - Total sleep (percent asleep, full 0-60 h) -------------------
total_sleep <- dt_plot[genotype %in% geno_order,
                       .(percent_asleep = mean(asleep, na.rm = TRUE) * 100),
                       by = .(id, genotype)]
total_sleep <- remove_outliers_iqr(total_sleep, "percent_asleep", "genotype")

fig2 <- two_group_boxplot(total_sleep, "percent_asleep",
                          title = "Total daily sleep is unchanged",
                          ylab  = "Percent asleep", ylim = c(0, 100))
fig_save(fig2, file.path(FIG_DIR, "Fig2_female_total_sleep"),
         width = FIG_W_1COL, height = FIG_H_1COL)

pairwise.wilcox.test(total_sleep$percent_asleep, total_sleep$genotype)

# ---- FIGURE 3 - Daytime vs nighttime sleep ---------------------------------
day_sleep   <- day_dt[,   .(percent_asleep = mean(asleep, na.rm = TRUE) * 100),
                      by = .(id, genotype)][, phase := "Day (lights on)"]
night_sleep <- night_dt[, .(percent_asleep = mean(asleep, na.rm = TRUE) * 100),
                        by = .(id, genotype)][, phase := "Night (lights off)"]

phase_sleep <- rbind(day_sleep, night_sleep)
phase_sleep <- remove_outliers_iqr(phase_sleep, "percent_asleep", c("genotype", "phase"))
phase_sleep[, phase := factor(phase, levels = c("Day (lights on)", "Night (lights off)"))]

fig3 <- phase_boxplot(phase_sleep, "percent_asleep",
                      title = "Sleep is redistributed from day to night in CS females",
                      ylab  = "Percent asleep", ylim = c(0, 105))
fig_save(fig3, file.path(FIG_DIR, "Fig3_female_day_vs_night_sleep"),
         width = FIG_W_WIDE, height = FIG_H_WIDE)

# ============================================================
#  ACTIVITY FIGURES  (same set as sleep, mirrored)
# ============================================================

# ---- FIGURE 4 - Activity profile (beam crosses / min) across ZT ------------
perfly_act <- dt_profile[genotype %in% geno_order,
                         .(val = mean(activity, na.rm = TRUE)),   # mean beam crosses / min per bin
                         by = .(genotype, id, zbin)]

fig4 <- build_profile(perfly_act,
                      ylab  = "Activity (beam crosses / min)",
                      title = "CS females are less active during the night",
                      legend_position = c(0.5, 0.8))
fig_save(fig4, file.path(FIG_DIR, "Fig4_female_activity_profile"),
         width = FIG_W_2COL, height = FIG_H_WIDE)

# ---- FIGURE 5 - Total activity (mean beam crosses per fly, full 0-60 h) -----
total_act <- dt_plot[genotype %in% geno_order,
                     .(activity = mean(activity, na.rm = TRUE)),
                     by = .(id, genotype)]
total_act <- remove_outliers_iqr(total_act, "activity", "genotype")

fig5 <- two_group_boxplot(total_act, "activity",
                          title = "Total daily activity is unchanged",
                          ylab  = "Mean activity (beam crosses / min)",
                          ylim  = c(0, NA))
fig_save(fig5, file.path(FIG_DIR, "Fig5_female_total_activity"),
         width = FIG_W_1COL, height = FIG_H_1COL)

pairwise.wilcox.test(total_act$activity, total_act$genotype)

# ---- FIGURE 6 - Daytime vs nighttime activity ------------------------------
day_act   <- day_dt[,   .(activity = mean(activity, na.rm = TRUE)),
                    by = .(id, genotype)][, phase := "Day (lights on)"]
night_act <- night_dt[, .(activity = mean(activity, na.rm = TRUE)),
                      by = .(id, genotype)][, phase := "Night (lights off)"]

phase_act <- rbind(day_act, night_act)
phase_act <- remove_outliers_iqr(phase_act, "activity", c("genotype", "phase"))
phase_act[, phase := factor(phase, levels = c("Day (lights on)", "Night (lights off)"))]

fig6 <- phase_boxplot(phase_act, "activity",
                      title = "Nighttime activity is selectively reduced in CS females",
                      ylab  = "Mean activity (beam crosses / min)",
                      ylim  = c(0, NA))
fig_save(fig6, file.path(FIG_DIR, "Fig6_female_day_vs_night_activity"),
         width = FIG_W_WIDE, height = FIG_H_WIDE)
