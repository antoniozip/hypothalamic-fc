#!/usr/bin/env Rscript
# First-level statistics: neuron-as-observation mixed-effects models.
#
# Model: lmer(metric ~ condition * region + (1 | animal/neuron), data = df)
#
# Reads all metrics CSV files and fits separate models per metric x estimator.
# Includes convergence handling (bobyqa optimizer) with progressive fallbacks:
#   1. Full interaction: condition * region + (1 | animal/neuron)
#   2. Main effects only: condition + region + (1 | animal/neuron)
#   3. Random-intercept only: condition + (1 | animal)
# Output: results/stats_first_level_<estimator>.csv

library(here)
library(lmerTest)
library(lme4)

args <- commandArgs(trailingOnly = TRUE)
estimator <- "glmcc"
for (i in seq_along(args)) {
  if (args[i] == "--estimator") estimator <- args[i + 1]
}

animals <- c(
  "171019", "171207", "171208", "171213", "180110", "180111",
  "180131", "180221", "180228", "180302", "180419", "180420", "180423"
)
conditions <- c("ongoing", "lightON")
cond_aliases <- list(ongoing = c("ongoing", "ongoing_bis"), lightON = c("lightON"))

# --- Load metrics ---
metrics_dir <- here("results", estimator)
all_data <- data.frame()

for (animal in animals) {
  for (condition in conditions) {
    for (alias in cond_aliases[[condition]]) {
      fpath <- file.path(metrics_dir, paste0("metrics_", animal, "_", alias, ".csv"))
      if (file.exists(fpath)) break
    }
    if (!file.exists(fpath)) next
    df <- read.csv(fpath, stringsAsFactors = FALSE)
    df$animal <- as.character(animal)
    df$condition <- condition
    all_data <- rbind(all_data, df)
  }
}

if (nrow(all_data) == 0) {
  cat("No metrics files found for estimator:", estimator, "\n")
  cat("Run compute_graph_metrics.R first.\n")
  quit(save = "no")
}

all_data$condition <- factor(all_data$condition, levels = c("ongoing", "lightON"))
all_data$region <- factor(all_data$region)
all_data$animal <- factor(all_data$animal)
all_data$neuron <- with(all_data, interaction(animal, neuron_id, drop = TRUE))

metrics <- c("node_strength", "clustering_coefficient", "local_efficiency", "hub_score")
results_table <- data.frame()

# --- Model fitting with progressive fallbacks ---
fit_model <- function(metric_name, df, estimator_name) {
  forms <- list(
    interaction = as.formula(paste(metric_name, "~ condition * region + (1 | animal/neuron)")),
    main_effects = as.formula(paste(metric_name, "~ condition + region + (1 | animal/neuron)")),
    intercept_only = as.formula(paste(metric_name, "~ condition + (1 | animal)"))
  )
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

  for (level in names(forms)) {
    model <- tryCatch(
      lmer(forms[[level]], data = df, REML = TRUE, control = ctrl),
      error = function(e) NULL,
      warning = function(w) {
        # bobyqa can produce convergence warnings even when fit is fine;
        # try Nelder-Mead as fallback before giving up
        tryCatch(
          lmer(forms[[level]], data = df, REML = TRUE,
               control = lmerControl(optimizer = "Nelder_Mead",
                                     optCtrl = list(maxfun = 200000))),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(model)) {
      cat(sprintf("  %-20s | %-14s | %s (level=%s)\n",
                  metric_name, estimator_name, "converged", level))
      return(list(model = model, level = level))
    }
  }
  cat(sprintf("  %-20s | %-14s | FAILED all levels\n", metric_name, estimator_name))
  return(NULL)
}

cat("\n=== Fitting first-level models (estimator:", estimator, ") ===\n")
cat(sprintf("  %-20s | %-14s | %s\n", "metric", "estimator", "status"))
cat(strrep("-", 60), "\n")

for (metric in metrics) {
  if (!metric %in% colnames(all_data)) next

  result <- fit_model(metric, all_data, estimator)
  if (is.null(result)) {
    cat(sprintf("  Model failed for %s: not enough factor levels or data\n", metric))
    next
  }

  model <- result$model
  model_level <- result$level

  cat(sprintf("\n--- %s (%s) [level=%s] ---\n", metric, estimator, model_level))
  print(summary(model))
  anova_res <- tryCatch(
    anova(model, type = "III", ddf = "Kenward-Roger"),
    error = function(e) {
      cat("  ANOVA failed, using type II:\n")
      anova(model, type = "II")
    }
  )
  print(anova_res)

  # Fixed effects
  fixed <- summary(model)$coefficients
  fixed_df <- data.frame(
    metric = metric,
    estimator = estimator,
    model_level = model_level,
    term = rownames(fixed),
    estimate = fixed[, "Estimate"],
    std_error = fixed[, "Std. Error"],
    df = if ("df" %in% colnames(fixed)) fixed[, "df"] else NA,
    t_value = if ("t value" %in% colnames(fixed)) fixed[, "t value"] else NA,
    p_value = fixed[, "Pr(>|t|)"],
    row.names = NULL
  )
  results_table <- rbind(results_table, fixed_df)

  # ANOVA table
  anova_df <- data.frame(
    metric = metric,
    estimator = estimator,
    model_level = model_level,
    term = rownames(anova_res),
    num_df = anova_res[["NumDF"]],
    den_df = anova_res[["DenDF"]],
    f_value = anova_res[["F value"]],
    p_value = anova_res[["Pr(>F)"]],
    row.names = NULL
  )
  results_table <- rbind(results_table, data.frame(
    metric = metric, estimator = estimator, model_level = model_level,
    term = paste0("ANOVA_", anova_df$term),
    estimate = anova_df$f_value,
    std_error = NA, df = anova_df$num_df, t_value = NA, p_value = anova_df$p_value
  ))
}

out_path <- here("results", paste0("stats_first_level_", estimator, ".csv"))
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results_table, out_path, row.names = FALSE)
cat("\nResults written to", out_path, "\n")

# --- Model diagnostics summary ---
cat("\n=== Model diagnostics ===\n")
for (metric in intersect(metrics, colnames(all_data))) {
  result <- fit_model(metric, all_data, estimator)
  if (!is.null(result)) {
    model <- result$model
    res <- residuals(model)
    shp <- tryCatch(shapiro.test(sample(res, min(5000, length(res)))), error = function(e) NULL)
    cat(sprintf("  %-20s | level=%-14s | AIC=%.1f | BIC=%.1f | n_obs=%d | Shapiro_W=%.3f (p=%.3f)\n",
                metric, result$level,
                AIC(model), BIC(model), nobs(model),
                if (!is.null(shp)) shp$statistic else NA,
                if (!is.null(shp)) shp$p.value else NA))
  }
}
cat("\n")
