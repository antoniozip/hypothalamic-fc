#!/usr/bin/env Rscript
# Tier 2 Analyses — Surrogate Null Models & Small-World Propensity
#
# 6. Surrogate null models: Are graph metrics higher than expected by chance?
#    Compare observed metrics against degree-preserving random networks.
# 7. Small-world: Compute sigma (σ) and omega (ω) indices.
#
# Output:
#   figures/main/figT2_null_models.png
#   figures/main/figT2_small_world.png
#   results/glmcc/null_models.csv
#   results/glmcc/small_world.csv

library(here)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)

dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)

animals <- c(
  "171019","171207","171208","171213","180110","180111",
  "180131","180221","180228","180302","180419","180420","180423"
)
daynight <- c(
  "171019"="D","171207"="N","171208"="N","171213"="N",
  "180110"="N","180111"="N","180131"="D","180221"="N",
  "180228"="N","180302"="D","180419"="D","180420"="D","180423"="D"
)

N_NULL <- 200  # null networks per animal

# ═══════════════════════════════════════════════════════════════════
# 6. SURROGATE NULL MODELS
# ═══════════════════════════════════════════════════════════════════
cat("=== 6. Surrogate Null Models ===\n")

null_results <- data.frame()

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0
  g <- graph_from_adjacency_matrix(adj_bin, mode = "directed")
  g_und <- as_undirected(g, mode = "collapse")
  n <- vcount(g)

  # --- Observed metrics ---
  obs_mod <- modularity(cluster_louvain(g_und))
  obs_cc  <- transitivity(g_und, type = "global")
  obs_pl  <- mean_distance(g_und, directed = FALSE)

  # --- Null distribution ---
  null_mod <- numeric(N_NULL)
  null_cc  <- numeric(N_NULL)
  null_pl  <- numeric(N_NULL)

  for (s in seq_len(N_NULL)) {
    # Degree-preserving rewiring
    g_rand <- rewire(g, keeping_degseq(niter = vcount(g) * 20))
    g_rand_und <- as_undirected(g_rand, mode = "collapse")

    null_mod[s] <- modularity(cluster_louvain(g_rand_und))
    null_cc[s]  <- transitivity(g_rand_und, type = "global")
    null_pl[s]  <- if (is_connected(g_rand_und)) {
      mean_distance(g_rand_und, directed = FALSE)
    } else {
      NA
    }
  }

  # --- Z-scores ---
  z_mod <- (obs_mod - mean(null_mod, na.rm = TRUE)) /
            sd(null_mod, na.rm = TRUE)
  z_cc  <- (obs_cc  - mean(null_cc, na.rm = TRUE)) /
            sd(null_cc, na.rm = TRUE)
  z_pl  <- (obs_pl  - mean(null_pl, na.rm = TRUE)) /
            sd(null_pl, na.rm = TRUE)

  # --- p-values ---
  p_mod <- sum(null_mod >= obs_mod, na.rm = TRUE) / sum(!is.na(null_mod))
  p_cc  <- sum(null_cc  >= obs_cc,  na.rm = TRUE) / sum(!is.na(null_cc))
  p_pl  <- sum(null_pl  >= obs_pl,  na.rm = TRUE) / sum(!is.na(null_pl))

  null_results <- rbind(null_results, data.frame(
    animal = animal, day_night = daynight[animal], n_neurons = n,
    obs_modularity = obs_mod, null_modularity_mean = mean(null_mod, na.rm = TRUE),
    z_modularity = z_mod, p_modularity = p_mod,
    obs_clustering = obs_cc, null_clustering_mean = mean(null_cc, na.rm = TRUE),
    z_clustering = z_cc, p_clustering = p_cc,
    obs_path_length = obs_pl, null_path_length_mean = mean(null_pl, na.rm = TRUE),
    z_path_length = z_pl, p_path_length = p_pl,
    row.names = NULL
  ))

  cat(sprintf("  %s: Z_mod=%.2f Z_cc=%.2f Z_pl=%.2f\n",
              animal, z_mod, z_cc, z_pl))
}

write.csv(null_results, here("results", "glmcc", "null_models.csv"),
          row.names = FALSE)

# Summary
cat(sprintf("\n  Mean Z_mod: %.2f ± %.2f (%.0f/13 p<0.05)\n",
            mean(null_results$z_modularity), sd(null_results$z_modularity),
            sum(null_results$p_modularity < 0.05)))
cat(sprintf("  Mean Z_cc:  %.2f ± %.2f (%.0f/13 p<0.05)\n",
            mean(null_results$z_clustering), sd(null_results$z_clustering),
            sum(null_results$p_clustering < 0.05)))
cat(sprintf("  Mean Z_pl:  %.2f ± %.2f (%.0f/13 p<0.05)\n",
            mean(null_results$z_path_length), sd(null_results$z_path_length),
            sum(null_results$p_path_length < 0.05)))

# --- Null models figure ---
# Panel A: Z-score distribution
null_long <- null_results %>%
  pivot_longer(c(z_modularity, z_clustering, z_path_length),
               names_to = "metric", values_to = "z_score") %>%
  mutate(metric = recode(metric,
    z_modularity = "Modularity",
    z_clustering = "Clustering",
    z_path_length = "Path Length"))

p1a <- ggplot(null_long, aes(x = metric, y = z_score, fill = metric)) +
  geom_boxplot(alpha = 0.5, outlier.size = 1.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888") +
  geom_hline(yintercept = 1.96, linetype = "dotted", color = "#b2182b", alpha = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "", y = "Z-score vs Degree-Preserving Null",
       title = "T6a: Graph Metric Z-Scores (vs 200 null networks)",
       subtitle = "Z > 1.96 = significantly higher than chance") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

# Panel B: Observed vs Null per metric
null_scatter <- null_results %>%
  pivot_longer(c(obs_modularity, obs_clustering, obs_path_length),
               names_to = "obs_metric", values_to = "obs_value") %>%
  mutate(
    null_value = case_when(
      obs_metric == "obs_modularity" ~ null_modularity_mean,
      obs_metric == "obs_clustering" ~ null_clustering_mean,
      obs_metric == "obs_path_length" ~ null_path_length_mean
    ),
    metric_label = case_when(
      obs_metric == "obs_modularity" ~ "Modularity",
      obs_metric == "obs_clustering" ~ "Clustering",
      obs_metric == "obs_path_length" ~ "Path Length"
    )
  )

p1b <- ggplot(null_scatter,
              aes(x = null_value, y = obs_value, color = day_night)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888") +
  facet_wrap(~ metric_label, scales = "free") +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = "Null Mean (degree-preserving)", y = "Observed Value",
       title = "T6b: Observed vs Null Graph Metrics",
       subtitle = "Points above line = higher than expected by degree alone",
       color = "Phase") +
  theme_minimal(base_size = 11)

fig6 <- ggpubr::ggarrange(p1a, p1b, ncol = 1, heights = c(1, 1.3))
ggsave(here("figures", "main", "figT2_null_models.png"),
       fig6, width = 10, height = 10, dpi = 150)


# ═══════════════════════════════════════════════════════════════════
# 7. SMALL-WORLD PROPENSITY
# ═══════════════════════════════════════════════════════════════════
cat("\n=== 7. Small-World Propensity ===\n")

sw_results <- data.frame()

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0
  g <- graph_from_adjacency_matrix(adj_bin, mode = "directed")
  g_und <- as_undirected(g, mode = "collapse")
  n <- vcount(g)

  if (!is_connected(g_und)) {
    cat(sprintf("  %s: disconnected, skipping\n", animal))
    next
  }

  # Observed metrics
  C_obs <- transitivity(g_und, type = "global")
  L_obs <- mean_distance(g_und, directed = FALSE)

  # Equivalent random network (Erdos-Renyi, same n and density)
  density_g <- edge_density(g_und)
  g_rand <- sample_gnp(n, density_g)
  # Ensure connected
  attempts <- 0
  while (!is_connected(g_rand) && attempts < 100) {
    g_rand <- sample_gnp(n, density_g)
    attempts <- attempts + 1
  }
  C_rand <- transitivity(g_rand, type = "global")
  L_rand <- mean_distance(g_rand, directed = FALSE)

  # Equivalent lattice network
  # Build a ring lattice with same number of edges
  m <- ecount(g_und)
  # Average degree
  k_avg <- round(2 * m / n)
  if (k_avg %% 2 == 1) k_avg <- k_avg + 1  # must be even for watts.strogatz
  g_latt <- sample_smallworld(1, n, k_avg / 2, 0)  # p=0 = lattice
  C_latt <- transitivity(g_latt, type = "global")
  L_latt <- mean_distance(g_latt, directed = FALSE)

  # Small-world indices (Telesford et al. 2011)
  # σ = (C/C_rand) / (L/L_rand)
  # ω = (L_rand/L) - (C/C_latt)
  sigma <- (C_obs / C_rand) / (L_obs / L_rand)
  omega <- (L_rand / L_obs) - (C_obs / C_latt)

  sw_results <- rbind(sw_results, data.frame(
    animal = animal, day_night = daynight[animal], n_neurons = n,
    C_obs = C_obs, C_rand = C_rand, C_latt = C_latt,
    L_obs = L_obs, L_rand = L_rand, L_latt = L_latt,
    sigma = sigma, omega = omega,
    row.names = NULL
  ))

  cat(sprintf("  %s: σ=%.3f ω=%.3f (C=%.3f L=%.2f)\n",
              animal, sigma, omega, C_obs, L_obs))
}

write.csv(sw_results, here("results", "glmcc", "small_world.csv"),
          row.names = FALSE)

cat(sprintf("\n  Mean σ: %.3f ± %.3f\n",
            mean(sw_results$sigma), sd(sw_results$sigma)))
cat(sprintf("  Mean ω: %.3f ± %.3f\n",
            mean(sw_results$omega), sd(sw_results$omega)))
cat(sprintf("  Small-world (σ>1): %d/%d animals\n",
            sum(sw_results$sigma > 1), nrow(sw_results)))
cat(sprintf("  Small-world (ω≈0): %d/%d animals (|ω|<0.5)\n",
            sum(abs(sw_results$omega) < 0.5), nrow(sw_results)))

# --- Small-world figure ---
# Panel A: σ and ω bar chart
sw_long <- sw_results %>%
  pivot_longer(c(sigma, omega), names_to = "index", values_to = "value")
sw_long$animal <- factor(sw_long$animal, levels = animals)

p2a <- ggplot(sw_long, aes(x = animal, y = value, fill = index)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(sigma = "#2166ac", omega = "#b2182b"),
                    labels = c(sigma = "σ (sigma)", omega = "ω (omega)")) +
  geom_hline(data = data.frame(index = "sigma", yint = 1),
             aes(yintercept = yint), linetype = "dashed", color = "#888888") +
  geom_hline(data = data.frame(index = "omega", yint = 0),
             aes(yintercept = yint), linetype = "dashed", color = "#888888") +
  labs(x = "Animal", y = "Index Value",
       title = "T7a: Small-World Indices",
       subtitle = paste0("sigma > 1 and omega ~ 0 = small-world. ",
                         "Mean sigma = ", round(mean(sw_results$sigma), 2),
                         ", Mean omega = ", round(mean(sw_results$omega), 3)),
       fill = "Index") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Clustering vs Path Length (normalized)
sw_norm <- sw_results %>%
  mutate(
    C_norm = C_obs / C_rand,
    L_norm = L_obs / L_rand,
    is_sw = sigma > 1 & abs(omega) < 0.5
  )

p2b <- ggplot(sw_norm, aes(x = C_norm, y = L_norm, color = day_night)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_text(aes(label = animal), size = 3, vjust = -1, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#888888") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#888888") +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = expression(C[obs] / C[rand] ~ "(normalized clustering)"),
       y = expression(L[obs] / L[rand] ~ "(normalized path length)"),
       title = "T7b: Normalized Clustering vs Path Length",
       subtitle = "Small-world: C_norm >> 1, L_norm ≈ 1 (top-left quadrant)",
       color = "Phase") +
  theme_minimal(base_size = 11)

# Panel C: C vs L ternary plot-like comparison
p2c <- ggplot(sw_results, aes(x = C_obs, y = L_obs, color = day_night)) +
  geom_point(size = 3.5, alpha = 0.8) +
  geom_point(aes(x = C_rand, y = L_rand), shape = 17, size = 3, alpha = 0.5,
             color = "#888888") +
  geom_segment(aes(xend = C_rand, yend = L_rand), alpha = 0.3,
               linetype = "dotted") +
  geom_text(aes(label = animal), size = 3, vjust = -1, show.legend = FALSE) +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = "Clustering Coefficient (C)", y = "Characteristic Path Length (L)",
       title = "T7c: C vs L (arrows point to equivalent random network)",
       subtitle = "Triangles = random; circles = observed. Long arrows = more small-world.",
       color = "Phase") +
  theme_minimal(base_size = 11)

fig7 <- ggpubr::ggarrange(p2a, ggpubr::ggarrange(p2b, p2c, ncol = 2),
                          ncol = 1, heights = c(1, 1.1))
ggsave(here("figures", "main", "figT2_small_world.png"),
       fig7, width = 12, height = 10, dpi = 150)

cat("\n=== Tier 2 complete ===\n")
