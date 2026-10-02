# ==============================================================================
# DomiDen | 04_run_and_plot.R   <<< RUN THIS SCRIPT <<<
# Single entry point: it loads 03_setup_config.R (packages, config, data, models),
# runs the dominance analysis for every dataset, and saves tables to results/
# and figures to figures/.
# Before running: edit the config block in scripts/03_setup_config.R.
#
# VERSION 1.0 | Last update 01.10.2026
# ==============================================================================
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)
source(here("scripts", "03_setup_config.R"))   # packages, config, data, models

# -----------------------------------------------------------------------------
# Helper: remove transformation prefixes to recover raw column names
# -----------------------------------------------------------------------------
strip_prefix <- function(x) sub("^(log_|arcsin_sqrt_|logit_)", "", x)

dir_results <- here("results"); dir.create(dir_results, showWarnings = FALSE)
dir_figures <- here("figures"); dir.create(dir_figures, showWarnings = FALSE)


# ==============================================================================
# 1. RUN THE ANALYSIS (class level and predictor level, all datasets)
# ==============================================================================
run_all <- function(fun) {
  res <- lapply(names(datasets), function(nm)
    fun(df = datasets[[nm]], model_list = Models_ensmbl,
        dataset_name = nm, cfg = config))
  bind_rows(lapply(res, `[[`, "dominance_summary"))
}

all_dominance <- run_all(Mdls_dominance_Analysis_OLS_AIC)                  # classes
absolute_dominance_predictorLvl <- run_all(Mdls_dominance_Analysis_OLS_AIC_predictor)  # predictors


# ==============================================================================
# 2. LABELS, ORDER AND COLORS (all derived from config)
# ==============================================================================
ds_names    <- names(config$datasets)
class_keys  <- names(config$pool)
class_label <- setNames(gsub("_", " ", class_keys), class_keys)

# Match a user-supplied named vector against the expected names.
# Missing entries get a default; unknown entries are dropped with a warning.
resolve_named <- function(user, expected, default, what) {
  out <- setNames(default, expected)
  if (!is.null(user)) {
    unknown <- setdiff(names(user), expected)
    if (length(unknown) > 0)
      warning(what, ": ignoring names not found: ",
              paste(unknown, collapse = ", "), call. = FALSE)
    keep <- intersect(names(user), expected)
    out[keep] <- user[keep]
  }
  out
}

auto_pal <- function(n) grDevices::hcl.colors(n, "Dark 3")

ds_colors <- resolve_named(config$dataset_colors, ds_names,
                           auto_pal(length(ds_names)), "dataset_colors")

cl_colors <- resolve_named(config$class_colors, class_keys,
                           auto_pal(length(class_keys)), "class_colors")
names(cl_colors) <- class_label[names(cl_colors)]

all_raw <- unlist(config$pool, use.names = FALSE)
pred_labels <- resolve_named(config$predictor_labels, all_raw,
                             all_raw, "predictor_labels")
pred_label <- function(raw) unname(pred_labels[raw])

# Shared y axis upper limit: a round number above the largest CI
y_max   <- ceiling(max(c(all_dominance$CI_upper,
                         absolute_dominance_predictorLvl$CI_upper), na.rm = TRUE) / 10) * 10
y_breaks <- seq(0, y_max, by = 10)

# Shared plot theme
theme_domiden <- function(aspect, ytext_size) {
  theme_bw() +
    theme(
      aspect.ratio = aspect,
      text = element_text(family = "sans"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.x = element_line(linewidth = 0.6, color = "grey90"),
      panel.grid.minor.x = element_line(linewidth = 0.1, color = "grey87"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      plot.title    = element_text(size = 12, colour = "grey30"),
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      axis.text.y   = element_text(size = ytext_size, colour = "grey40"),
      axis.text.x   = element_text(size = 8, colour = "grey40"),
      axis.ticks    = element_blank(),
      axis.title.x  = element_text(size = 9, colour = "grey40"),
      axis.title.y  = element_blank(),
      plot.margin   = unit(c(0.1, 0.2, 0.1, 0.1), "cm"),
      legend.position = "none"
    )
}


# ==============================================================================
# 3. CLASS-LEVEL TABLE AND FIGURE
# ==============================================================================
all_dominance$Class   <- factor(class_label[all_dominance$Class],
                                levels = rev(class_label))   # first class on top
all_dominance$Dataset <- factor(all_dominance$Dataset, levels = ds_names)
all_dominance <- all_dominance %>% group_by(Dataset) %>%
  arrange(desc(Class), .by_group = TRUE)

write.csv(all_dominance, file.path(dir_results, "DominanceResults_byClass.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

DA_plot <- ggplot(all_dominance,
                  aes(x = Class, y = Mean_Importance,
                      fill = Dataset, shape = Dataset)) +
  geom_segment(aes(xend = Class, y = 0, yend = Mean_Importance),   # stem
               linewidth = 1, color = "grey70") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),             # 95% CI
                width = 0.2, linewidth = 0.23, color = "grey10") +
  geom_point(size = 2.2, color = "grey20", alpha = 0.7, stroke = 0.3) +  # head
  scale_fill_manual(values = ds_colors) +
  scale_shape_manual(values = c(22, 21, 25)) +
  facet_wrap(~Dataset) +
  coord_flip() +
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = c(0, 0)) +
  labs(title = "Drivers of Denudation: Dominance Analysis",
       y = "General Dominance \n(% of Explained Variance)", x = "") +
  theme_domiden(aspect = 0.9, ytext_size = 10)

print(DA_plot)
ggsave(file.path(dir_figures, "Fig_DomiDen_byClass.pdf"), DA_plot,
       width = 13, height = 7.7, units = "cm", dpi = 300)


# ==============================================================================
# 4. PREDICTOR-LEVEL TABLE AND FIGURE
# ==============================================================================
absolute_dominance_predictorLvl$Class <-
  class_label[absolute_dominance_predictorLvl$Class]

# Axis order: classes in pool order, predictors in pool order inside each class,
# with an empty "gap" level between classes (first item ends up on top).
levels_pred <- character(0); gap_pos <- integer(0)
for (i in seq_along(class_keys)) {
  levels_pred <- c(levels_pred, rev(pred_label(config$pool[[class_keys[i]]])))
  if (i < length(class_keys)) {
    levels_pred <- c(levels_pred, strrep(" ", i))   # unique blank label
    gap_pos <- c(gap_pos, length(levels_pred))
  }
}
levels_pred <- rev(levels_pred); gap_pos <- length(levels_pred) + 1 - gap_pos

absolute_dominance_predictorLvl$Predictor <- factor(
  pred_label(strip_prefix(absolute_dominance_predictorLvl$Predictor)),
  levels = levels_pred)
absolute_dominance_predictorLvl$Dataset <-
  factor(absolute_dominance_predictorLvl$Dataset, levels = ds_names)
absolute_dominance_predictorLvl <- absolute_dominance_predictorLvl %>%
  group_by(Dataset) %>% arrange(desc(Predictor), .by_group = TRUE)

write.csv(absolute_dominance_predictorLvl,
          file.path(dir_results, "DominanceResults_byPredictor.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

DomiDen_predictor_plot <- ggplot(
  absolute_dominance_predictorLvl,
  aes(x = Predictor, y = Mean_Importance, fill = Dataset,
      shape = Dataset, color = Class)) +
  geom_segment(aes(xend = Predictor, y = 0, yend = Mean_Importance, group = Class),
               position = position_dodge(width = 0.2), linewidth = 1.9) +
  geom_vline(xintercept = gap_pos, linewidth = 0.4, color = "grey80") + # class separators
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.3, linewidth = 0.23, color = "grey10") +
  geom_point(size = 2.2, color = "grey20", alpha = 0.7, stroke = 0.3) +
  scale_color_manual(values = cl_colors) +
  scale_fill_manual(values = ds_colors) +
  scale_shape_manual(values = c(22, 21, 25)) +
  facet_wrap(~Dataset) +
  coord_flip() +
  scale_x_discrete(drop = FALSE) +   # keep the empty gap levels
  scale_y_continuous(breaks = y_breaks, limits = c(0, y_max), expand = c(0, 0)) +
  labs(title = "Predictor-Level Dominance Analysis",
       subtitle = "LMG Absolute Importance (Mean across top models)",
       y = "General Dominance\n(% of Total Variance Explained)", x = "") +
  theme_domiden(aspect = 1.85, ytext_size = 8)

print(DomiDen_predictor_plot)
fig_height <- max(8, 0.55 * length(levels_pred) + 3)
ggsave(file.path(dir_figures, "Fig_DomiDen_byPredictor.pdf"), DomiDen_predictor_plot,
       width = 12.5, height = fig_height, units = "cm", dpi = 300)


# ==============================================================================
# 5. PRINT SUMMARY TABLES
# ==============================================================================
print(knitr::kable(all_dominance %>%
                     select(Dataset, Class, Mean_Importance, SD_Importance, SE_Importance)))
print(knitr::kable(absolute_dominance_predictorLvl %>%
                     select(Dataset, Class, Predictor, Mean_Importance, SD_Importance, SE_Importance)))
message("Done. Tables in results/, figures in figures/.")

