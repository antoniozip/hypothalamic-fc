#!/usr/bin/env Rscript
# Tier 1 Analyses — Rich Club, Edge Consistency, Region Decomposition, Centrality
#
# All use existing validated adjacency CSVs + neurons.csv.
#
# Output:
#   figures/main/figT1_rich_club.png
#   figures/main/figT1_edge_consistency.png
#   figures/main/figT1_region_decomposition.png
#   figures/main/figT1_centrality_convergence.png
#   results/glmcc/rich_club.csv
#   results/glmcc/edge_consistency.csv
#   results/glmcc/centrality.csv

library(here)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

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

neurons <- read.csv(here("data", "processed", "neurons.csv"), stringsAsFactors = FALSE)

# ═══════════════════════════════════════════════════════════════════
# 1. RICH-CLUB ANALYSIS
# ═══════════════════════════════════════════════════════════════════
cat("\n=== 1. Rich-Club Analysis ===\n")

N_NULL <- 100
rich_club_all <- data.frame()

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0
  g <- graph_from_adjacency_matrix(adj_bin, mode = "directed")
  n <- vcount(g)

  # Observed rich-club coefficient
  degrees <- degree(g, mode = "all")
  max_k <- max(degrees)
  k_levels <- seq(1, max_k)

  phi_obs <- numeric(length(k_levels))
  for (ki in seq_along(k_levels)) {
    k <- k_levels[ki]
    rich_nodes <- which(degrees > k)
    if (length(rich_nodes) < 2) {
      phi_obs[ki] <- NA
      next
    }
    sg <- induced_subgraph(g, rich_nodes)
    e_rich <- ecount(sg)
    e_max <- length(rich_nodes) * (length(rich_nodes) - 1)
    phi_obs[ki] <- if (e_max > 0) e_rich / e_max else NA
  }

  # Null: degree-preserving rewired networks
  phi_null <- matrix(NA, nrow = length(k_levels), ncol = N_NULL)
  for (s in seq_len(N_NULL)) {
    g_rand <- rewire(g, keeping_degseq(niter = vcount(g) * 10))
    deg_rand <- degree(g_rand, mode = "all")
    for (ki in seq_along(k_levels)) {
      k <- k_levels[ki]
      rich_nodes <- which(deg_rand > k)
      if (length(rich_nodes) < 2) next
      sg <- induced_subgraph(g_rand, rich_nodes)
      e_rich <- ecount(sg)
      e_max <- length(rich_nodes) * (length(rich_nodes) - 1)
      phi_null[ki, s] <- if (e_max > 0) e_rich / e_max else NA
    }
  }

  # Normalized rich-club: φ_norm = φ_obs / mean(φ_null)
  phi_null_mean <- rowMeans(phi_null, na.rm = TRUE)
  phi_norm <- phi_obs / phi_null_mean
  phi_norm[phi_null_mean == 0 | is.na(phi_obs)] <- NA

  # Significance: p = fraction of null >= observed
  p_values <- numeric(length(k_levels))
  for (ki in seq_along(k_levels)) {
    null_vals <- phi_null[ki, ]
    null_vals <- null_vals[!is.na(null_vals)]
    if (length(null_vals) > 0 && !is.na(phi_obs[ki])) {
      p_values[ki] <- sum(null_vals >= phi_obs[ki]) / length(null_vals)
    } else {
      p_values[ki] <- NA
    }
  }

  # Store
  for (ki in seq_along(k_levels)) {
    rich_club_all <- rbind(rich_club_all, data.frame(
      animal = animal, day_night = daynight[animal],
      k = k_levels[ki], n_nodes = n,
      phi_obs = phi_obs[ki], phi_null_mean = phi_null_mean[ki],
      phi_norm = phi_norm[ki], p_value = p_values[ki],
      row.names = NULL
    ))
  }

  cat(sprintf("  %s: max_k=%d, significant levels: %d\n",
              animal, max_k, sum(p_values < 0.05, na.rm = TRUE)))
}

write.csv(rich_club_all, here("results", "glmcc", "rich_club.csv"), row.names = FALSE)

# --- Rich-club figure ---
# Panel A: φ_norm curves per animal
rc_valid <- rich_club_all[!is.na(rich_club_all$phi_norm) &
                          is.finite(rich_club_all$phi_norm), ]

p1a <- ggplot(rc_valid, aes(x = k, y = phi_norm, color = animal, group = animal)) +
  geom_line(linewidth = 0.7, alpha = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#888888") +
  labs(x = "Degree (k)", y = expression(phi[norm](k)),
       title = "T1a: Normalized Rich-Club Coefficient",
       subtitle = expression(phi[norm] > 1 ~ "indicates rich-club organization")) +
  scale_color_viridis_d(option = "D") +
  theme_minimal(base_size = 11)

# Panel B: Mean φ_norm ± SE across animals with significance
rc_group <- rc_valid %>%
  group_by(k) %>%
  summarise(
    mean_phi = mean(phi_norm, na.rm = TRUE),
    se_phi = sd(phi_norm, na.rm = TRUE) / sqrt(n()),
    n_sig = sum(p_value < 0.05, na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  ) %>%
  filter(n_total >= 3)

p1b <- ggplot(rc_group, aes(x = k, y = mean_phi)) +
  geom_ribbon(aes(ymin = mean_phi - se_phi, ymax = mean_phi + se_phi),
              fill = "#6baed6", alpha = 0.3) +
  geom_line(color = "#08519c", linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#888888") +
  geom_point(aes(size = n_sig / n_total), color = "#b2182b", alpha = 0.7) +
  labs(x = "Degree (k)", y = expression(phi[norm](k) ~ "(mean ± SE)"),
       title = "T1b: Group Rich-Club (mean across animals)",
       subtitle = "Red points: fraction of animals with significant rich-club at degree k",
       size = "Sig. fraction") +
  theme_minimal(base_size = 11)

# Panel C: φ_obs vs φ_null for top k levels
rc_top <- rc_valid %>%
  group_by(animal) %>%
  slice_max(k, n = 3) %>%
  ungroup()

p1c <- ggplot(rc_top, aes(x = phi_null_mean, y = phi_obs, color = day_night)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888") +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = expression(phi[null]), y = expression(phi[obs]),
       title = "T1c: Observed vs Null Rich-Club (top-3 degree levels)",
       subtitle = "Points above line = richer than chance",
       color = "Phase") +
  theme_minimal(base_size = 11)

fig1 <- ggpubr::ggarrange(p1a, p1b, p1c, ncol = 1, heights = c(1, 1, 1.1))
ggsave(here("figures", "main", "figT1_rich_club.png"),
       fig1, width = 9, height = 14, dpi = 150)


# ═══════════════════════════════════════════════════════════════════
# 2. EDGE CONSISTENCY (CORE CONNECTOME)
# ═══════════════════════════════════════════════════════════════════
cat("\n=== 2. Edge Consistency ===\n")

# Build region-pair edge consistency matrix
# For each animal, extract binary adjacency
# For each region pair (r1, r2), compute: fraction of possible edges that exist
# Then average across animals

all_regions <- sort(unique(neurons$region6))
n_regions <- length(all_regions)
consistency_mat <- matrix(0, n_regions, n_regions,
                          dimnames = list(all_regions, all_regions))
count_mat <- matrix(0, n_regions, n_regions,
                    dimnames = list(all_regions, all_regions))

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0

  animal_neurons <- neurons[neurons$animal == animal, ]
  if (nrow(animal_neurons) == 0) next
  regions_animal <- animal_neurons$region6[
    match(seq_len(nrow(adj)), animal_neurons$neuron_id)]
  regions_animal[is.na(regions_animal)] <- "unknown"

  for (r1 in all_regions) {
    idx1 <- which(regions_animal == r1)
    if (length(idx1) == 0) next
    for (r2 in all_regions) {
      idx2 <- which(regions_animal == r2)
      if (length(idx2) == 0) next
      submat <- adj_bin[idx1, idx2, drop = FALSE]
      n_possible <- length(idx1) * length(idx2)
      if (r1 == r2) n_possible <- length(idx1) * (length(idx1) - 1)
      if (n_possible > 0) {
        consistency_mat[r1, r2] <- consistency_mat[r1, r2] +
          sum(submat) / n_possible
        count_mat[r1, r2] <- count_mat[r1, r2] + 1
      }
    }
  }
}

# Average across animals
consistency_avg <- consistency_mat / pmax(count_mat, 1)
consistency_avg[count_mat < 3] <- NA  # require at least 3 animals

cat("  Mean edge consistency:", round(mean(consistency_avg, na.rm = TRUE), 3), "\n")

# Save
consistency_df <- as.data.frame(as.table(consistency_avg))
names(consistency_df) <- c("source", "target", "consistency")
consistency_df <- consistency_df[!is.na(consistency_df$consistency), ]
write.csv(consistency_df, here("results", "glmcc", "edge_consistency.csv"),
          row.names = FALSE)

# Top consistent edges
consistency_df <- consistency_df[order(-consistency_df$consistency), ]
cat("  Top 5 most consistent region pairs:\n")
for (i in 1:min(5, nrow(consistency_df))) {
  cat(sprintf("    %s ↔ %s: %.3f\n",
              consistency_df$source[i], consistency_df$target[i],
              consistency_df$consistency[i]))
}

# --- Edge consistency figure ---
# Heatmap
cons_mat_plot <- consistency_avg
cons_mat_plot[is.na(cons_mat_plot)] <- 0

p2a <- ggplot(as.data.frame(as.table(cons_mat_plot)),
              aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c",
                      limits = c(0, max(cons_mat_plot, na.rm = TRUE)),
                      na.value = "grey90") +
  labs(x = "Source Region", y = "Target Region",
       title = "T2a: Edge Consistency Across Animals",
       subtitle = "Fraction of within-region-pair edges present (mean across animals)",
       fill = "Consistency") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bar chart of consistency by region
region_cons <- data.frame(region = all_regions,
                          mean_consistency = rowMeans(consistency_avg, na.rm = TRUE))
region_cons <- region_cons[order(-region_cons$mean_consistency), ]
region_cons$region <- factor(region_cons$region, levels = region_cons$region)

p2b <- ggplot(region_cons, aes(x = region, y = mean_consistency)) +
  geom_col(fill = "#6baed6", alpha = 0.85) +
  labs(x = "", y = "Mean Edge Consistency",
       title = "T2b: Mean Edge Consistency by Region") +
  coord_flip() +
  theme_minimal(base_size = 11)

fig2 <- ggpubr::ggarrange(p2a, p2b, ncol = 1, heights = c(1.3, 1))
ggsave(here("figures", "main", "figT1_edge_consistency.png"),
       fig2, width = 9, height = 11, dpi = 150)


# ═══════════════════════════════════════════════════════════════════
# 3. WITHIN vs BETWEEN REGION DECOMPOSITION
# ═══════════════════════════════════════════════════════════════════
cat("\n=== 3. Within vs Between Region Decomposition ===\n")

decomp_all <- data.frame()

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_bin <- (abs(adj) > 0) * 1
  diag(adj_bin) <- 0
  n <- nrow(adj)

  animal_neurons <- neurons[neurons$animal == animal, ]
  regions_animal <- animal_neurons$region6[
    match(seq_len(n), animal_neurons$neuron_id)]
  regions_animal[is.na(regions_animal)] <- "unknown"

  total_edges <- sum(adj_bin)
  within_edges <- 0
  between_edges <- 0

  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) next
      if (adj_bin[i, j] == 1) {
        if (regions_animal[i] == regions_animal[j]) {
          within_edges <- within_edges + 1
        } else {
          between_edges <- between_edges + 1
        }
      }
    }
  }

  decomp_all <- rbind(decomp_all, data.frame(
    animal = animal, day_night = daynight[animal],
    n_neurons = n,
    total_edges = total_edges,
    within_edges = within_edges,
    between_edges = between_edges,
    pct_within = within_edges / total_edges * 100,
    pct_between = between_edges / total_edges * 100,
    row.names = NULL
  ))

  cat(sprintf("  %s: %.1f%% within, %.1f%% between\n",
              animal, within_edges / total_edges * 100,
              between_edges / total_edges * 100))
}

cat(sprintf("\n  Group mean: %.1f%% within, %.1f%% between\n",
            mean(decomp_all$pct_within), mean(decomp_all$pct_between)))

# --- Decomposition figure ---
decomp_long <- decomp_all %>%
  pivot_longer(c(pct_within, pct_between),
               names_to = "type", values_to = "percentage") %>%
  mutate(type = recode(type, pct_within = "Within-Region",
                       pct_between = "Between-Region"))

decomp_long$animal <- factor(decomp_long$animal, levels = animals)

p3 <- ggplot(decomp_long, aes(x = animal, y = percentage, fill = type)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c(`Within-Region` = "#2166ac",
                                `Between-Region` = "#b2182b"),
                    name = "Edge Type") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#888888") +
  labs(x = "Animal", y = "Percentage of Total Edges",
       title = "T3: Within-Region vs Between-Region Connectivity",
       subtitle = paste0("Group mean: ",
                         round(mean(decomp_all$pct_within), 1),
                         "% within, ",
                         round(mean(decomp_all$pct_between), 1),
                         "% between")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("figures", "main", "figT1_region_decomposition.png"),
       p3, width = 9, height = 5, dpi = 150)


# ═══════════════════════════════════════════════════════════════════
# 4. CENTRALITY CONVERGENCE
# ═══════════════════════════════════════════════════════════════════
cat("\n=== 4. Centrality Convergence ===\n")

centrality_all <- data.frame()

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  adj_abs <- abs(adj)
  diag(adj_abs) <- 0
  g <- graph_from_adjacency_matrix(adj_abs, mode = "directed",
                                   weighted = TRUE, diag = FALSE)
  n <- vcount(g)

  # Compute 4 centrality measures
  s <- strength(g, mode = "all")
  b <- betweenness(g, directed = TRUE, weights = 1 / (E(g)$weight + 0.001))
  # Eigenvector centrality needs undirected
  g_und <- as_undirected(g, mode = "collapse")
  eig <- tryCatch(eigen_centrality(g_und)$vector,
                  error = function(e) rep(NA, n))
  pr <- page_rank(g, directed = TRUE)$vector

  # Store
  animal_neurons <- neurons[neurons$animal == animal, ]
  regions_animal <- animal_neurons$region6[
    match(seq_len(n), animal_neurons$neuron_id)]
  regions_animal[is.na(regions_animal)] <- "unknown"

  for (i in seq_len(n)) {
    centrality_all <- rbind(centrality_all, data.frame(
      animal = animal, day_night = daynight[animal],
      neuron_id = i, region = regions_animal[i],
      strength = s[i],
      betweenness = b[i],
      eigenvector = eig[i],
      pagerank = pr[i],
      row.names = NULL
    ))
  }

  # Within-animal Spearman correlations
  cor_mat <- cor(cbind(s, b, eig, pr), method = "spearman",
                 use = "pairwise.complete.obs")
  cat(sprintf("  %s: mean ρ=%.3f\n", animal,
              mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)))
}

write.csv(centrality_all, here("results", "glmcc", "centrality.csv"),
          row.names = FALSE)

# --- Centrality convergence figure ---
# Panel A: Correlation heatmap per animal
animals_with_data <- unique(centrality_all$animal)
cor_summary <- data.frame()
for (a in animals_with_data) {
  a_dat <- centrality_all[centrality_all$animal == a, ]
  cmat <- cor(a_dat[, c("strength", "betweenness", "eigenvector", "pagerank")],
              method = "spearman", use = "pairwise.complete.obs")
  measures <- colnames(cmat)
  for (i in 1:(length(measures)-1)) {
    for (j in (i+1):length(measures)) {
      cor_summary <- rbind(cor_summary, data.frame(
        animal = a,
        pair = paste(measures[i], measures[j], sep = " × "),
        rho = cmat[i, j],
        row.names = NULL
      ))
    }
  }
}

p4a <- ggplot(cor_summary, aes(x = pair, y = rho)) +
  geom_boxplot(fill = "#6baed6", alpha = 0.5, outlier.size = 1) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "#08519c") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#888888") +
  labs(x = "Centrality Measure Pair", y = "Spearman ρ",
       title = "T4a: Centrality Convergence Across Animals",
       subtitle = paste0("Mean ρ = ",
                         round(mean(cor_summary$rho, na.rm = TRUE), 3))) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Panel B: Strength vs Betweenness scatter (all neurons pooled)
set.seed(42)
sample_idx <- sample(nrow(centrality_all), min(5000, nrow(centrality_all)))
centrality_sample <- centrality_all[sample_idx, ]
centrality_sample <- centrality_sample[
  !is.na(centrality_sample$betweenness) &
  centrality_sample$betweenness > 0, ]

p4b <- ggplot(centrality_sample,
              aes(x = strength, y = betweenness, color = day_night)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = "Node Strength (log)", y = "Betweenness (log)",
       title = "T4b: Strength vs Betweenness (all neurons)",
       color = "Phase") +
  theme_minimal(base_size = 11)

fig4 <- ggpubr::ggarrange(p4a, p4b, ncol = 1, heights = c(1, 1.2))
ggsave(here("figures", "main", "figT1_centrality_convergence.png"),
       fig4, width = 10, height = 9, dpi = 150)

cat("\n=== Tier 1 analyses #1-4 complete ===\n")
