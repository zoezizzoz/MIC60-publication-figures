library(survival)
library(survminer)
library(tidyverse)
library(readxl)   # for Excel files

# ---- Read data ----
# Point this at your file. Works for .csv, .tsv/.txt, or .xlsx/.xls.
# Needs a group column ("Type") and a survival-time column (e.g. "Days").
data_file <- "Fly Lifespan Streamlined.csv"   # <-- change to your filename

ext <- tolower(tools::file_ext(data_file))
df <- switch(ext,
  csv        = readr::read_csv(data_file, show_col_types = FALSE),
  tsv        = readr::read_tsv(data_file, show_col_types = FALSE),
  txt        = readr::read_tsv(data_file, show_col_types = FALSE),
  xlsx       = readxl::read_excel(data_file),
  xls        = readxl::read_excel(data_file),
  stop("Unsupported file type: ", ext)
)

# Find the survival-time column (accepts "Days", "Days Alive", etc.)
day_col <- names(df)[str_detect(tolower(names(df)), "day")][1]
df <- df %>% rename(days = all_of(day_col))
df$days <- suppressWarnings(as.numeric(df$days))  # blanks/text -> NA

# Blank Days cells -> NA. Here we DROP them (no death recorded).
# If a blank instead means "still alive / censored at end of experiment",
# replace the two lines below with:
#   df$event <- ifelse(is.na(df$days), 0, 1)          # 0 = censored
#   df$days  <- ifelse(is.na(df$days), max(df$days, na.rm=TRUE), df$days)
df <- df %>% filter(!is.na(days))
df$event <- 1                                          # every recorded obs = death

# ---- Remove outliers (Tukey's IQR fences, computed WITHIN each group) ----
# A point is an outlier if it falls outside [Q1 - k*IQR, Q3 + k*IQR].
# k = 1.5 is the standard "outlier" fence; use 3.0 for "extreme" only.
remove_outliers <- TRUE      # set FALSE to keep all points
iqr_k           <- 1.5       # fence multiplier

if (remove_outliers) {
  df <- df %>%
    group_by(Type) %>%
    mutate(
      .q1  = quantile(days, 0.25),
      .q3  = quantile(days, 0.75),
      .iqr = .q3 - .q1,
      .lo  = .q1 - iqr_k * .iqr,
      .hi  = .q3 + iqr_k * .iqr,
      .outlier = days < .lo | days > .hi
    ) %>%
    ungroup()

  cat(sprintf("Removed %d outlier(s):\n", sum(df$.outlier)))
  print(df %>% filter(.outlier) %>% select(Type, days))

  df <- df %>% filter(!.outlier) %>% select(-starts_with("."))
}

# Group order controls colour assignment (WT black, PR green, CS blue)
df$Type <- factor(df$Type, levels = c("WT", "PR", "CS"))
df$Type <- droplevels(df$Type)

# ---- Kaplan-Meier fit + log-rank test ----
fit  <- survfit(Surv(days, event) ~ Type, data = df)
lr   <- survdiff(Surv(days, event) ~ Type, data = df)
pval <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

# ---- Plot ----
# Colors consistent to rest of graphs
pal <- c(WT = "#00BFC4", PR = "#33cc33", CS = "#F8766D")[levels(df$Type)]

g <- ggsurvplot(
  fit, data = df, conf.int = FALSE, censor = FALSE,
  palette = unname(pal), size = 1.1,
  legend = c(0.72, 0.9), legend.title = "", legend.labs = levels(df$Type),
  xlab = "Days", ylab = "Percent survival",
  xlim = c(0, 30), break.x.by = 5, ylim = c(0, 100),
  surv.scale = "percent", ggtheme = theme_classic(base_size = 16)
)

med <- summary(fit)$table[, "median"]
pct <- round(100 * (med["Type=CS"] - med["Type=WT"]) / med["Type=WT"], 2)
sub <- sprintf("log-rank test (WT Vs CS), P=%.4f\nlife span %+.2f%%", pval, pct)

g$plot <- g$plot +
  scale_y_continuous(labels = function(x) x * 100, breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(caption = sub) +
  theme(
    axis.line   = element_line(linewidth = 0.9, colour = "black"),
    axis.ticks  = element_line(linewidth = 0.9, colour = "black"),
    axis.text   = element_text(colour = "black", face = "bold"),
    axis.title  = element_text(face = "bold"),
    legend.text = element_text(size = 13),
    plot.caption = element_text(hjust = 0, size = 13),
    plot.margin  = margin(10, 15, 10, 10)
  )

ggsave("survival_curve.png", plot = g$plot, width = 6.2, height = 5.2, dpi = 300, bg = "white")
