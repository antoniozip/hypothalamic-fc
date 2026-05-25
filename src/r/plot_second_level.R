#!/usr/bin/env Rscript
# Day/Night forest plot from second-level mixed-effects model.
#
# Reads stats_second_level_glmcc.csv and plots condition + region_pair
# fixed effects as a forest plot.
#
# Output: figures/main/figS7_second_level_forest.png

library(here)
library(ggplot2)
library(dplyr)

dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)

stats_file <- here("results", "stats_second_level_glmcc.csv")
if (!file.exists(stats_file)) {
  cat("No second-level stats file. Run stats_second_level.R first.\n")
  quit(save = "no")
}

df <- read.csv(stats_file, stringsAsFactors = FALSE)

# Filter to fixed effects with valid SE
fixed <- df[!grepl("^ANOVA_|nonEst", df$term) & !is.na(df$std_error) & df$std_error > 0, ]
fixed <- fixed[!is.na(fixed$estimate), ]

if (nrow(fixed) == 0) {
  cat("No valid fixed effects in second-level stats.\n")
  quit(save = "no")
}

# Classify terms
fixed$term_type <- "Region Pair"
fixed$term_type[grepl("condition", fixed$term)] <- "Condition"
fixed$term_type[grepl("Intercept", fixed$term)] <- "Intercept"

# Sort by estimate
fixed$term <- factor(fixed$term, levels = fixed$term[order(fixed$estimate)])

# Add significance
fixed$sig <- ""
fixed$sig[fixed$p_value < 0.05]  <- "*"
fixed$sig[fixed$p_value < 0.01]  <- "**"
fixed$sig[fixed$p_value < 0.001] <- "***"

# Direction
fixed$direction <- ifelse(fixed$estimate > 0, "positive", "negative")

cat(sprintf("Plotting %d fixed effects\n", nrow(fixed)))

p <- ggplot(fixed, aes(x = estimate, y = term, color = term_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(xmin = estimate - 1.96 * std_error,
                     xmax = estimate + 1.96 * std_error),
                width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = sig,
                x = estimate + sign(estimate) *
                  (1.96 * std_error + max(abs(estimate)) * 0.06)),
            size = 3, show.legend = FALSE, color = "#333333") +
  scale_color_manual(
    values = c(Condition = "#b2182b", `Region Pair` = "#444444", Intercept = "#999999"),
    name = "Term Type"
  ) +
  labs(
    x = "Fixed Effect Estimate ± 95% CI",
    y = "",
    title = "Fig S7: Second-level mixed-effects model (region-pair density)",
    subtitle = "lmer(density ~ condition + region_pair + (1 | animal)) — GLMCC, CCG-validated",
    caption = paste0(
      "NOTE: Model unstable due to incomplete ongoing data (5/13 animals missing). ",
      "Non-estimable contrasts excluded.\n",
      "* p<0.05  ** p<0.01  *** p<0.001"
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.caption = element_text(hjust = 0, size = 8, color = "#888888"),
    legend.position = "bottom"
  )

# Height based on number of terms
n_terms <- nrow(fixed)
plot_height <- max(6, n_terms * 0.22)

ggsave(here("figures", "main", "figS7_second_level_forest.png"),
       p, width = 11, height = plot_height, dpi = 150, limitsize = FALSE)

cat(sprintf("  figS7_second_level_forest.png (%d terms)\n", n_terms))
