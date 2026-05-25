#!/usr/bin/env Rscript
# Second-level statistics: day x night x condition x region mixed-effects.
#
# lmer(density ~ condition * region * day_night + (1 | animal), data = df_density)
#
# Pre-specified contrasts:
# 1. LightON - Ongoing | Day
# 2. LightON - Ongoing | Night
# 3. (LightON - Ongoing | Night) - (LightON - Ongoing | Day)
# 4. Region-pair-specific tests

library(here)
library(lmerTest)
library(lme4)
library(emmeans)

args <- commandArgs(trailingOnly = TRUE)
estimator <- "glmcc"
for (i in seq_along(args)) {
  if (args[i] == "--estimator") estimator <- args[i + 1]
}

animals <- c(
  "day1", "night1", "night2", "night3", "night4", "night5",
  "day2", "night6", "night7", "day3", "day4", "day5", "day6"
)
conditions <- c("ongoing", "lightON")
cond_aliases <- list(ongoing = c("ongoing", "ongoing_bis"), lightON = c("lightON"))
daynight <- c(
  "day1" = "D", "night1" = "N", "night2" = "N", "night3" = "N",
  "night4" = "N", "night5" = "N", "day2" = "D", "night6" = "N",
  "night7" = "N", "day3" = "D", "day4" = "D", "day5" = "D",
  "day6" = "D"
)

density_dir <- here("results", estimator)
all_density <- data.frame()

for (animal in animals) {
  for (condition in conditions) {
    for (alias in cond_aliases[[condition]]) {
      fpath <- file.path(density_dir, paste0("density_long_", animal, "_", alias, ".csv"))
      if (file.exists(fpath)) break
    }
    if (!file.exists(fpath)) next
    df <- read.csv(fpath, stringsAsFactors = FALSE)
    df$condition <- condition  # standardise to "ongoing" / "lightON"
    df$day_night <- daynight[animal]
    all_density <- rbind(all_density, df)
  }
}

if (nrow(all_density) == 0) {
  cat("No density files found for estimator:", estimator, "\n")
  cat("Run aggregate_density.R first.\n")
  quit(save = "no")
}

all_density$condition <- factor(all_density$condition, levels = c("ongoing", "lightON"))
all_density$day_night <- factor(all_density$day_night, levels = c("D", "N"))
all_density$source <- factor(all_density$source)
all_density$target <- factor(all_density$target)
all_density$animal <- factor(all_density$animal)
all_density$region_pair <- with(all_density, interaction(source, target, sep = "-"))

formula <- density ~ condition * region_pair * day_night + (1 | animal)
model <- tryCatch(
  lmer(formula, data = all_density, REML = TRUE),
  error = function(e) {
    cat("Full model failed:", e$message, "\n")
    cat("Fitting reduced model: density ~ condition * day_night + (1 | animal)\n")
    tryCatch(
      lmer(density ~ condition * day_night + (1 | animal),
           data = all_density, REML = TRUE),
      error = function(e2) {
        cat("Reduced model also failed:", e2$message, "\n")
        return(NULL)
      }
    )
  }
)

if (is.null(model)) {
  quit(save = "no")
}

cat("\n=== Second-level model (", estimator, ") ===\n")
print(summary(model))
anova_res <- anova(model, type = "III", ddf = "Kenward-Roger")
print(anova_res)

emm <- emmeans(model, pairwise ~ condition | day_night, adjust = "fdr")
cat("\n=== Contrast: LightON - Ongoing | Day ===\n")
print(emm[[1]][1, ])

cat("\n=== Contrast: LightON - Ongoing | Night ===\n")
print(emm[[1]][2, ])

emm_interaction <- emmeans(model, ~ condition * day_night, adjust = "fdr")
cat("\n=== Contrasts with FDR correction ===\n")
print(pairs(emm_interaction))

out_path <- here("results", paste0("stats_second_level_", estimator, ".csv"))
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

fixed <- summary(model)$coefficients
summary_df <- data.frame(
  estimator = estimator,
  term = rownames(fixed),
  estimate = fixed[, "Estimate"],
  std_error = fixed[, "Std. Error"],
  t_value = fixed[, "t value"],
  row.names = NULL
)
write.csv(summary_df, out_path, row.names = FALSE)
cat("\nResults written to", out_path, "\n")
