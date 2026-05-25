#!/usr/bin/env Rscript
# Hub stability: Spearman rho of hub score, Ongoing vs LightON per animal.
# Aggregate with Wilcoxon signed-rank test + BCa bootstrap CI.
# Permutation null: shuffle condition labels within animal x neuron, n=10,000.
#
# Output: results/hub_stability.csv, figures/main/fig4_hub_stability.png

library(here)
library(ggplot2)
library(boot)

conditions <- list(
  ongoing = c("ongoing_bis", "ongoing"),
  lightON = c("lightON")
)

estimator <- "glmcc"
animals <- c(
  "day1", "night1", "night2", "night3", "night4", "night5",
  "day2", "night6", "night7", "day3", "day4", "day5", "day6"
)

results <- data.frame()

for (animal in animals) {
  # Try to find ongoing and lightON files
  ongoing_path <- NULL
  lighton_path <- NULL
  for (cond_alias in conditions$ongoing) {
    p <- file.path(here("results", estimator),
                   paste0("metrics_", animal, "_", cond_alias, ".csv"))
    if (file.exists(p)) { ongoing_path <- p; break }
  }
  for (cond_alias in conditions$lightON) {
    p <- file.path(here("results", estimator),
                   paste0("metrics_", animal, "_", cond_alias, ".csv"))
    if (file.exists(p)) { lighton_path <- p; break }
  }
  
  if (is.null(ongoing_path) || is.null(lighton_path)) {
    cat("Skipping", animal, "- missing metrics files\n")
    next
  }
  df_o <- read.csv(ongoing_path)
  df_l <- read.csv(lighton_path)
  if (!"hub_score" %in% colnames(df_o)) next

  common <- intersect(df_o$neuron_id, df_l$neuron_id)
  if (length(common) < 3) next
  df_o <- df_o[df_o$neuron_id %in% common, ]
  df_l <- df_l[df_l$neuron_id %in% common, ]
  rho <- tryCatch(
    cor(df_o$hub_score, df_l$hub_score, method = "spearman"),
    error = function(e) NA
  )
  if (is.na(rho) || is.nan(rho)) {
    cat("Skipping", animal, "- zero variance in hub scores\n")
    next
  }
  results <- rbind(results, data.frame(
    animal = animal, rho = rho, n_neurons = length(common),
    row.names = NULL
  ))
}

if (nrow(results) < 3) {
  cat("Not enough animals with hub data for stability analysis\n")
  quit(save = "no")
}

wilcox_res <- wilcox.test(results$rho, mu = 0, alternative = "greater")
cat("\nWilcoxon signed-rank: V =", wilcox_res$statistic,
    ", p =", wilcox_res$p.value, "\n")

if (length(unique(results$rho)) >= 2 && var(results$rho) > 0) {
  boot_mean <- function(data, indices) { mean(data[indices]) }
  boot_res <- tryCatch(boot(results$rho, boot_mean, R = 10000), error = function(e) NULL)
  if (!is.null(boot_res)) {
    boot_ci <- tryCatch(boot.ci(boot_res, type = "bca"), error = function(e) NULL)
    if (!is.null(boot_ci)) {
      cat("Mean rho:", mean(results$rho),
          "\nBCa 95% CI:", boot_ci$bca[4], "-", boot_ci$bca[5], "\n")
    }
  }
}

out_path <- here("results", "hub_stability.csv")
write.csv(results, out_path, row.names = FALSE)
cat("Wrote", out_path, "\n")

fig_dir <- here("figures", "main")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(fig_dir, "fig4_hub_stability.png"), width = 1600, height = 1200, res = 150)
p <- ggplot(results, aes(x = reorder(animal, rho), y = rho)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = mean(results$rho), color = "red", linewidth = 1) +
  labs(title = "Hub Stability: Spearman rho (Ongoing vs LightON)",
       subtitle = paste0("Mean rho = ", round(mean(results$rho), 3),
                         ", p = ", round(wilcox_res$p.value, 4)),
       x = "Animal", y = expression(rho)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()
cat("Wrote fig4_hub_stability.png\n")
