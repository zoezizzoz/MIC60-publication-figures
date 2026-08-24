# ============================================================
# figure_style.R
# Shared figure style for the dMIC60 project
#
# Source this at the top of every analysis script:
#
#   source("figure_style.R")            # same folder
#   source(file.path(dirname(sys.frame(1)$ofile), "figure_style.R"))
#
# Everything that controls the look of a figure lives in the
# SETTINGS block below. Change a value once here and every
# figure in the project changes with it.
#
# What is standardized:
#   colors, fills, group labels
#   font family and every text size
#   error bars (color, linewidth, cap width)
#   sample dots (shape, size, alpha, jitter)
#   time-series profiles (line width, ribbon transparency)
#   panel spacing, axis lines, ticks, margins
#   saved chart size and resolution
#   significance annotation: a bracket and p value are drawn
#   ONLY when p < FIG_ALPHA
# ============================================================

library(ggplot2)


# ============================================================
# SETTINGS
# ============================================================

# ---- typography --------------------------------------------
# "" uses the device default sans font, which is the portable
# choice. Set to "Helvetica" or "Arial" for final figures.
FIG_FONT <- ""

FIG_BASE_SIZE        <- 11
FIG_TITLE_SIZE       <- 12
FIG_SUBTITLE_SIZE    <- 9
FIG_CAPTION_SIZE     <- 8
FIG_AXIS_TITLE_SIZE  <- 11
FIG_AXIS_TEXT_SIZE   <- 10
FIG_LEGEND_TEXT_SIZE <- 9.5

# geom_text sizes are in mm; 3.2 mm is about 9 pt
FIG_ANNOT_SIZE <- 3.2

# Publication panel policy: panel titles, interpretive subtitles, sample-size
# prose, and methods captions belong in the manuscript figure legend rather
# than inside the graph.  Analysis scripts may still call labs() so their
# intended wording remains recoverable, but theme_fig() suppresses those
# elements in exported panels.
FIG_SHOW_PLOT_TITLE    <- FALSE
FIG_SHOW_PLOT_SUBTITLE <- FALSE
FIG_SHOW_PLOT_CAPTION  <- FALSE


# ---- colors -------------------------------------------------
# One palette for the whole project: control blue, mutant red.
FIG_COLORS <- c(
  control = "#4C72B0",
  mutant  = "#A8322B"
)

# Every genotype name used anywhere in the project maps onto one
# of the two roles above. Add new names here, not in the scripts.
FIG_ROLE_OF <- c(
  "WR"        = "control",
  "WT"        = "control",
  "MIC60WT"   = "control",
  "dMIC60-WT" = "control",
  "dMIC60WT"  = "control",
  "dMIC60WR"  = "control",
  "CS"        = "mutant",
  "CSF"       = "mutant",
  "MIC60CS"   = "mutant",
  "dMIC60-CS" = "mutant",
  "dMIC60CS"  = "mutant"
)

# How far fills are lightened toward white (0 = full color, 1 = white)
FIG_FILL_LIGHTEN <- 0.45

FIG_BOX_ALPHA   <- 1
FIG_BOX_OUTLINE <- "grey25"

FIG_REF_LINE_COLOR     <- "grey60"
FIG_REF_LINE_TYPE      <- "dashed"
FIG_REF_LINE_LINEWIDTH <- 0.4


# ---- error bars ---------------------------------------------
FIG_EB_COLOR     <- "black"
FIG_EB_LINEWIDTH <- 0.5
FIG_EB_CAP       <- 0.20   # horizontal end-cap width in x-axis units
FIG_EB_TYPE      <- "sem"  # "sem", "sd", or "iqr"


# ---- sample dots --------------------------------------------
FIG_PT_SHAPE       <- 16
FIG_PT_COLOR       <- "black"
FIG_PT_SIZE        <- 1.9  # normal n
FIG_PT_SIZE_DENSE  <- 1.1  # hundreds of points in one column
FIG_PT_ALPHA       <- 0.75
FIG_PT_ALPHA_DENSE <- 0.45
FIG_JITTER_WIDTH   <- 0.12
FIG_SEED           <- 2026

# Marker for a nested summary value, e.g. one median per image
FIG_MARKER_SHAPE  <- 23
FIG_MARKER_SIZE   <- 2.6
FIG_MARKER_FILL   <- "white"
FIG_MARKER_COLOR  <- "black"
FIG_MARKER_STROKE <- 0.7


# ---- time-series profiles ----------------------------------
# Shared appearance for sleep/activity profiles and other
# continuous traces with uncertainty ribbons.
FIG_PROFILE_LINEWIDTH <- 0.8
FIG_RIBBON_ALPHA      <- 0.20


# ---- geometry and spacing -----------------------------------
FIG_LINE_WIDTH      <- 0.5   # axis lines and ticks
FIG_TICK_LENGTH_PT  <- 2.5
FIG_SUMMARY_WIDTH   <- 0.35  # box / crossbar width
FIG_DODGE_WIDTH     <- 0.75  # spacing for grouped boxes, points, and brackets
FIG_X_EXPAND        <- 0.40  # padding either side of the outer groups
FIG_Y_EXPAND        <- c(0.02, 0.16)
FIG_MARGIN_PT       <- c(t = 8, r = 10, b = 8, l = 8)

# Faceted gene panels can opt into these dividers without changing the
# borderless default used by the project's single-panel figures.
FIG_FACET_BORDER_COLOR <- "grey35"
FIG_FACET_STRIP_FILL   <- "grey92"
FIG_FACET_GAP_PT       <- 3


# ---- saved chart size ---------------------------------------
FIG_W_1COL <- 4.5   # single two-group panel
FIG_H_1COL <- 4.8
FIG_W_WIDE <- 5.4   # panel with a legend or a continuous x axis
FIG_H_WIDE <- 4.2
FIG_W_2COL <- 9.0   # two panels across
FIG_H_2COL <- 8.8
FIG_DPI    <- 600

# Built-in PDF output avoids platform-specific Cairo DLL failures.
# Trade-off: grDevices::pdf cannot encode non-ASCII glyphs, so the female
# symbol in the genotype labels becomes a dot in the PDF. PNG output is
# unaffected. Switch to grDevices::cairo_pdf on a machine where Cairo works
# if those glyphs are needed in the vector file.
FIG_PDF_DEVICE <- grDevices::pdf


# ---- statistics display -------------------------------------
FIG_ALPHA   <- 0.05   # bracket and p value appear only below this
FIG_SHOW_NS <- FALSE  # TRUE writes "ns" instead of drawing nothing

FIG_BRACKET_LINEWIDTH <- 0.5
FIG_BRACKET_OFFSET    <- 0.10  # bracket height above the data, as a
FIG_BRACKET_TICK      <- 0.035 # fraction of the plotted data span

# Summary geometry for every categorical panel:
#   "box"  median and interquartile range
#   "bars" mean with capped error bars
FIG_SUMMARY_STYLE <- "box"

# Literal reading of "error bars only where the difference is
# significant". Off by default because hiding the spread of a
# non-significant group is hard to defend in a figure legend.
# Set TRUE to drop error bars from non-significant comparisons.
FIG_ERRORBARS_ONLY_IF_SIG <- FALSE


# ============================================================
# COLOR HELPERS
# ============================================================

fig_lighten <- function(color, amount = FIG_FILL_LIGHTEN) {
  channels <- grDevices::col2rgb(color)
  mixed <- channels + (255 - channels) * amount
  
  grDevices::rgb(
    mixed[1, ],
    mixed[2, ],
    mixed[3, ],
    maxColorValue = 255
  )
}

fig_role <- function(group_names) {
  roles <- unname(FIG_ROLE_OF[as.character(group_names)])
  
  if (anyNA(roles)) {
    stop(
      "No color role defined for: ",
      paste(unique(group_names[is.na(roles)]), collapse = ", "),
      ". Add it to FIG_ROLE_OF in figure_style.R."
    )
  }
  
  roles
}

# Line / outline color for each group
fig_colors <- function(group_names) {
  out <- unname(FIG_COLORS[fig_role(group_names)])
  names(out) <- as.character(group_names)
  out
}

# Fill color for each group
fig_fills <- function(group_names, lighten = FIG_FILL_LIGHTEN) {
  out <- fig_lighten(fig_colors(group_names), lighten)
  names(out) <- as.character(group_names)
  out
}

fig_scale_fill <- function(group_names, labels = waiver(), ...) {
  scale_fill_manual(
    values = fig_fills(group_names),
    limits = as.character(group_names),
    labels = labels,
    ...
  )
}

fig_scale_color <- function(group_names, labels = waiver(), ...) {
  scale_color_manual(
    values = fig_colors(group_names),
    limits = as.character(group_names),
    labels = labels,
    ...
  )
}

fig_scale_x_group <- function(labels = waiver()) {
  scale_x_discrete(
    labels = labels,
    expand = expansion(add = FIG_X_EXPAND)
  )
}

fig_scale_y <- function(name = waiver(), limits = c(0, NA), ...) {
  scale_y_continuous(
    name = name,
    limits = limits,
    expand = expansion(mult = FIG_Y_EXPAND),
    ...
  )
}


# ============================================================
# THEME
# ============================================================

theme_fig <- function(
    base_size = FIG_BASE_SIZE,
    base_family = FIG_FONT,
    legend_position = "none"
) {
  
  theme_classic(
    base_size = base_size,
    base_family = base_family
  ) +
    
    theme(
      legend.position = legend_position,
      legend.justification = "center",
      legend.key.size = unit(9, "pt"),
      legend.title = element_blank(),
      legend.text = element_text(
        size = FIG_LEGEND_TEXT_SIZE,
        color = "black"
      ),
      legend.margin = margin(b = 2),
      
      axis.line = element_line(
        linewidth = FIG_LINE_WIDTH,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = FIG_LINE_WIDTH,
        color = "black"
      ),
      
      axis.ticks.length = unit(FIG_TICK_LENGTH_PT, "pt"),
      
      axis.text = element_text(
        size = FIG_AXIS_TEXT_SIZE,
        color = "black"
      ),
      
      axis.text.x = element_text(
        size = FIG_AXIS_TEXT_SIZE,
        color = "black",
        lineheight = 0.9,
        margin = margin(t = 4)
      ),
      
      axis.title = element_text(
        size = FIG_AXIS_TITLE_SIZE,
        color = "black"
      ),
      
      axis.title.y = element_text(margin = margin(r = 7)),
      axis.title.x = element_text(margin = margin(t = 6)),
      
      plot.title = if (FIG_SHOW_PLOT_TITLE) element_text(
        size = FIG_TITLE_SIZE,
        face = "bold",
        color = "black",
        hjust = 0,
        lineheight = 1.0,
        margin = margin(b = 3)
      ) else element_blank(),
      
      plot.subtitle = if (FIG_SHOW_PLOT_SUBTITLE) element_text(
        size = FIG_SUBTITLE_SIZE,
        color = "grey35",
        hjust = 0,
        lineheight = 1.05,
        margin = margin(b = 6)
      ) else element_blank(),
      
      plot.caption = if (FIG_SHOW_PLOT_CAPTION) element_text(
        size = FIG_CAPTION_SIZE,
        color = "grey35",
        hjust = 0,
        lineheight = 1.1,
        margin = margin(t = 6)
      ) else element_blank(),
      
      strip.background = element_blank(),
      strip.text = element_text(
        size = FIG_AXIS_TEXT_SIZE,
        color = "black"
      ),
      
      panel.spacing = unit(8, "pt"),
      
      plot.margin = margin(
        t = FIG_MARGIN_PT[["t"]],
        r = FIG_MARGIN_PT[["r"]],
        b = FIG_MARGIN_PT[["b"]],
        l = FIG_MARGIN_PT[["l"]]
      )
    )
}


# ============================================================
# LAYERS
# ============================================================

# Reference line, e.g. the control mean at 100 percent
fig_reference_line <- function(yintercept) {
  geom_hline(
    yintercept = yintercept,
    linetype = FIG_REF_LINE_TYPE,
    color = FIG_REF_LINE_COLOR,
    linewidth = FIG_REF_LINE_LINEWIDTH
  )
}


# One dot per observation
fig_points <- function(dense = FALSE, width = FIG_JITTER_WIDTH, ...) {
  geom_point(
    position = position_jitter(
      width = width,
      height = 0,
      seed = FIG_SEED
    ),
    shape = FIG_PT_SHAPE,
    size = if (dense) FIG_PT_SIZE_DENSE else FIG_PT_SIZE,
    alpha = if (dense) FIG_PT_ALPHA_DENSE else FIG_PT_ALPHA,
    show.legend = FALSE,
    ...
  )
}


# Open marker for a nested summary value, e.g. one median per image
fig_markers <- function(data, mapping, width = 0.045, ...) {
  geom_point(
    data = data,
    mapping = mapping,
    inherit.aes = FALSE,
    position = position_jitter(
      width = width,
      height = 0,
      seed = FIG_SEED
    ),
    shape = FIG_MARKER_SHAPE,
    size = FIG_MARKER_SIZE,
    fill = FIG_MARKER_FILL,
    color = FIG_MARKER_COLOR,
    stroke = FIG_MARKER_STROKE,
    show.legend = FALSE,
    ...
  )
}


# Error bar interval functions
fig_interval <- function(type = FIG_EB_TYPE) {
  switch(
    type,
    
    sem = function(x) {
      x <- x[is.finite(x)]
      m <- mean(x)
      s <- stats::sd(x) / sqrt(length(x))
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    
    sd = function(x) {
      x <- x[is.finite(x)]
      m <- mean(x)
      s <- stats::sd(x)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    
    iqr = function(x) {
      x <- x[is.finite(x)]
      q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
      data.frame(y = q[2], ymin = q[1], ymax = q[3])
    },
    
    stop("FIG_EB_TYPE must be 'sem', 'sd', or 'iqr'.")
  )
}


# Standalone capped error bar, for plots that supply their own
# summary table (the dose-response curve, for example)
fig_errorbar <- function(
    data = NULL,
    mapping = NULL,
    width = FIG_EB_CAP,
    ...
) {
  geom_errorbar(
    data = data,
    mapping = mapping,
    width = width,
    linewidth = FIG_EB_LINEWIDTH,
    color = FIG_EB_COLOR,
    show.legend = FALSE,
    ...
  )
}


# The summary geometry for a categorical panel.
# style = "box"  median and interquartile range
# style = "bars" mean with capped error bars
# p is only used when FIG_ERRORBARS_ONLY_IF_SIG is TRUE.
fig_summary <- function(
    style = FIG_SUMMARY_STYLE,
    p_value = NA_real_,
    width = FIG_SUMMARY_WIDTH,
    error = FIG_EB_TYPE,
    position = "identity",
    show_legend = FALSE
) {
  
  hide_bars <- FIG_ERRORBARS_ONLY_IF_SIG &&
    (!is.finite(p_value) || p_value >= FIG_ALPHA)
  
  if (identical(style, "box")) {

    # Draw whiskers as an explicit error-bar layer. This guarantees visible
    # horizontal end caps even with ggplot2 versions where geom_boxplot()
    # does not draw staples by default.
    whiskers <- stat_boxplot(
      geom = "errorbar",
      width = FIG_EB_CAP,
      linewidth = FIG_EB_LINEWIDTH,
      color = FIG_EB_COLOR,
      position = position,
      show.legend = FALSE
    )

    box <- geom_boxplot(
      width = width,
      outlier.shape = NA,
      alpha = FIG_BOX_ALPHA,
      linewidth = FIG_LINE_WIDTH,
      color = FIG_BOX_OUTLINE,
      position = position,
      show.legend = show_legend
    )
    
    # Whiskers are the error bars of a box plot, so honor the
    # same switch by flattening them when they must be hidden.
    if (hide_bars) {
      box <- geom_boxplot(
        width = width,
        outlier.shape = NA,
        alpha = FIG_BOX_ALPHA,
        linewidth = FIG_LINE_WIDTH,
        color = FIG_BOX_OUTLINE,
        coef = 0,
        position = position,
        show.legend = show_legend
      )

      return(list(box))
    }

    return(list(whiskers, box))
  }
  
  if (identical(style, "bars")) {
    
    crossbar <- stat_summary(
      fun = mean,
      geom = "crossbar",
      width = width,
      linewidth = FIG_LINE_WIDTH,
      color = FIG_BOX_OUTLINE,
      alpha = FIG_BOX_ALPHA,
      position = position,
      show.legend = show_legend
    )
    
    if (hide_bars) {
      return(list(crossbar))
    }
    
    bars <- stat_summary(
      fun.data = fig_interval(error),
      geom = "errorbar",
      width = FIG_EB_CAP,
      linewidth = FIG_EB_LINEWIDTH,
      color = FIG_EB_COLOR,
      position = position,
      show.legend = show_legend
    )
    
    return(list(crossbar, bars))
  }
  
  stop("FIG_SUMMARY_STYLE must be 'box' or 'bars'.")
}


# ============================================================
# SIGNIFICANCE ANNOTATION
# ============================================================

fig_stars <- function(p_value) {
  if (!is.finite(p_value)) {
    return("")
  }
  
  if (p_value < 0.0001) {
    "****"
  } else if (p_value < 0.001) {
    "***"
  } else if (p_value < 0.01) {
    "**"
  } else if (p_value < FIG_ALPHA) {
    "*"
  } else {
    "ns"
  }
}

fig_p_label <- function(p_value) {
  if (!is.finite(p_value)) {
    return("")
  }
  
  if (p_value >= FIG_ALPHA) {
    return("ns")
  }
  
  paste0(
    fig_stars(p_value),
    "\np = ",
    format.pval(p_value, digits = 2, eps = 0.001)
  )
}


# A bracket plus p value, drawn only when p < FIG_ALPHA.
# `values` is the vector of plotted y values, used to place the
# bracket the same distance above the data in every panel.
# Returns NULL when the comparison is not significant, and
# `plot + NULL` is a no-op, so callers need no if statement.
fig_sig_bracket <- function(
    p_value,
    values,
    x1 = 1,
    x2 = 2,
    offset = FIG_BRACKET_OFFSET,
    label = NULL
) {
  
  if (!is.finite(p_value)) {
    return(NULL)
  }
  
  if (p_value >= FIG_ALPHA && !FIG_SHOW_NS) {
    return(NULL)
  }
  
  values <- values[is.finite(values)]
  
  y_max <- max(values)
  y_min <- min(values)
  
  span <- max(
    y_max - y_min,
    abs(y_max) * 0.15,
    1e-9
  )
  
  y_bracket <- y_max + span * offset
  y_tick <- span * FIG_BRACKET_TICK
  
  if (is.null(label)) {
    label <- fig_p_label(p_value)
  }
  
  list(
    annotate(
      "segment",
      x = x1,
      xend = x2,
      y = y_bracket,
      yend = y_bracket,
      linewidth = FIG_BRACKET_LINEWIDTH,
      color = "black"
    ),
    
    annotate(
      "segment",
      x = c(x1, x2),
      xend = c(x1, x2),
      y = y_bracket,
      yend = y_bracket - y_tick,
      linewidth = FIG_BRACKET_LINEWIDTH,
      color = "black"
    ),
    
    annotate(
      "text",
      x = (x1 + x2) / 2,
      y = y_bracket + y_tick * 0.6,
      label = label,
      vjust = 0,
      size = FIG_ANNOT_SIZE,
      lineheight = 0.95,
      family = FIG_FONT,
      color = "black"
    )
  )
}


# Draw brackets at explicit x positions for grouped comparisons.
# `comparisons` must contain x1, x2, and p_value columns. Each row
# is passed through fig_sig_bracket(), so FIG_ALPHA and FIG_SHOW_NS
# are applied exactly as they are for a single comparison.
fig_sig_brackets_at <- function(
    comparisons,
    values,
    offset = FIG_BRACKET_OFFSET
) {
  
  required_columns <- c("x1", "x2", "p_value")
  missing_columns <- setdiff(required_columns, names(comparisons))
  
  if (length(missing_columns) > 0) {
    stop(
      "comparisons is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (nrow(comparisons) == 0) {
    return(NULL)
  }
  
  layers <- lapply(seq_len(nrow(comparisons)), function(i) {
    fig_sig_bracket(
      p_value = comparisons$p_value[[i]],
      values = values,
      x1 = comparisons$x1[[i]],
      x2 = comparisons$x2[[i]],
      offset = offset
    )
  })
  
  unlist(layers, recursive = FALSE)
}


# ============================================================
# SAVING
# ============================================================

fig_save <- function(
    plot,
    path_base,
    width = FIG_W_1COL,
    height = FIG_H_1COL,
    dpi = FIG_DPI,
    formats = c("png", "pdf")
) {
  
  path_base <- sub("\\.(png|pdf)$", "", path_base)
  
  if ("png" %in% formats) {
    ggsave(
      filename = paste0(path_base, ".png"),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  }
  
  if ("pdf" %in% formats) {
    # A label the PDF device cannot encode raises one warning per glyph per
    # draw, which buries anything else in "50 or more warnings". Collapse
    # them into a single note; see FIG_PDF_DEVICE above for the fix.
    glyph_note_shown <- FALSE
    
    withCallingHandlers(
      ggsave(
        filename = paste0(path_base, ".pdf"),
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = FIG_PDF_DEVICE,
        bg = "white"
      ),
      warning = function(w) {
        if (grepl("conversion failure|mbcsToSbcs", conditionMessage(w))) {
          if (!glyph_note_shown) {
            glyph_note_shown <<- TRUE
            message(
              basename(paste0(path_base, ".pdf")),
              ": a non-ASCII glyph was replaced with a dot by the PDF ",
              "device. The PNG is unaffected."
            )
          }
          invokeRestart("muffleWarning")
        }
      }
    )
  }
  
  invisible(path_base)
}
