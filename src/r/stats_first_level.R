#!/usr/bin/env Rscript
# First-level statistics: neuron-as-observation mixed-effects models.
#
# lmer(metric ~ condition * region + (1 | animal/neuron), data = df)
#
# Reads all metrics CSV files and fits separate models per metric x estimator.
# Output: results/stats_first_level.csv

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

daynight <- c(
  "171019" = "D", "171207" = "N", "171208" = "N", "171213" = "N",
  "180110" = "N", "180111" = "N", "180131" = "D", "180221" = "N",
  "180228" = "N", "180302" = "D", "180419" = "D", "180420" = "D",
  "180423" = "D"
)

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
    df$day_night <- daynight[animal]
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

for (metric in metrics) {
  if (!metric %in% colnames(all_data)) next
  model1 <- tryCatch(
    lmer(as.formula(paste(metric, "~ condition + (1 | animal/neuron)")),
         data = all_data, REML = TRUE),
    error = function(e) {
      tryCatch(
        lmer(as.formula(paste(metric, "~ condition + (1 | animal)")),
             data = all_data, REML = TRUE),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(model1)) {
    cat("Model failed for", metric, ":", "not enough factor levels\n")
    next
  }

  cat("\n=== ", metric, " (", estimator, ") ===\n")
  print(summary(model1))
  anova_res <- anova(model1, type = "III", ddf = "Kenward-Roger")
  print(anova_res)

  fixed <- summary(model1)$coefficients
  fixed_df <- data.frame(
    metric = metric,
    estimator = estimator,
    term = rownames(fixed),
    estimate = fixed[, "Estimate"],
    std_error = fixed[, "Std. Error"],
    df = fixed[, "df"],
    t_value = fixed[, "t value"],
    p_value = fixed[, "Pr(>|t|)"],
    row.names = NULL
  )
  results_table <- rbind(results_table, fixed_df)

  anova_df <- data.frame(
    metric = metric,
    estimator = estimator,
    term = rownames(anova_res),
    num_df = anova_res[["NumDF"]],
    den_df = anova_res[["DenDF"]],
    f_value = anova_res[["F value"]],
    p_value = anova_res[["Pr(>F)"]],
    row.names = NULL
  )
  results_table <- rbind(results_table, data.frame(
    metric = metric, estimator = estimator,
    term = paste0("ANOVA_", anova_df$term),
    estimate = anova_df$f_value,
    std_error = NA, df = anova_df$num_df, t_value = NA, p_value = anova_df$p_value
  ))
}

out_path <- here("results", paste0("stats_first_level_", estimator, ".csv"))
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.csv(results_table, out_path, row.names = FALSE)
cat("\nResults written to", out_path, "\n")
