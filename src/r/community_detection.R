#!/usr/bin/env Rscript
# Community detection on CCG-validated adjacency matrices.
#
# Runs Louvain community detection on each animal's validated adjacency,
# generates per-animal community plots and cross-animal consistency analysis.
#
# Output:
#   figures/main/figS3_community_modularity.png  — modularity bar chart
#   figures/main/figS3_community_<animal>.png    — per-animal community plots
#   results/glmcc/communities.csv                 — community assignments

library(here)
library(igraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)

animals <- c(
  "day1","night1","night2","night3","night4","night5",
  "day2","night6","night7","day3","day4","day5","day6"
)
daynight <- c(
  "day1"="D","night1"="N","night2"="N","night3"="N",
  "night4"="N","night5"="N","day2"="D","night6"="N",
  "night7"="N","day3"="D","day4"="D","day5"="D","day6"="D"
)

# ─── Load neurons.csv for region labels ──────────────────────────
neurons <- read.csv(here("data", "processed", "neurons.csv"), stringsAsFactors = FALSE)

# ─── Community detection per animal ──────────────────────────────
all_communities <- data.frame()
modularity_scores <- data.frame()

cat("Running community detection...\n")

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  n_units <- nrow(adj)

  # Build undirected weighted graph from symmetric adjacency
  adj_sym <- (adj + t(adj)) / 2
  diag(adj_sym) <- 0
  g <- graph_from_adjacency_matrix(adj_sym, mode = "max",
                                   weighted = TRUE, diag = FALSE)

  # Louvain community detection
  comm <- cluster_louvain(g, weights = E(g)$weight)
  modularity_score <- modularity(comm)

  # Get region labels
  animal_neurons <- neurons[neurons$animal == animal, ]
  if (nrow(animal_neurons) > 0) {
    region_vec <- animal_neurons$region6[match(seq_len(n_units),
                                               animal_neurons$neuron_id)]
    region_vec[is.na(region_vec)] <- "unknown"
  } else {
    region_vec <- rep("unknown", n_units)
  }

  # Record
  n_comm <- length(unique(comm$membership))
  cat(sprintf("  %s (%s): %d communities, modularity=%.3f, n=%d\n",
              animal, daynight[animal], n_comm, modularity_score, n_units))

  modularity_scores <- rbind(modularity_scores, data.frame(
    animal = animal,
    day_night = daynight[animal],
    n_neurons = n_units,
    n_communities = n_comm,
    modularity = modularity_score,
    row.names = NULL
  ))

  # Save community assignments
  comm_df <- data.frame(
    animal = animal,
    day_night = daynight[animal],
    neuron_id = seq_len(n_units),
    region = region_vec,
    community = comm$membership,
    row.names = NULL
  )
  all_communities <- rbind(all_communities, comm_df)
}

# Save community assignments
write.csv(all_communities,
          here("results", "glmcc", "communities.csv"),
          row.names = FALSE)

# ─── Fig S3a: Modularity bar chart ────────────────────────────────
cat("Generating modularity summary...\n")

modularity_scores$animal <- factor(modularity_scores$animal, levels = animals)
modularity_scores$day_night <- factor(modularity_scores$day_night, levels = c("D", "N"))

p_mod <- ggplot(modularity_scores, aes(x = animal, y = modularity)) +
  geom_col(aes(fill = day_night), width = 0.7, alpha = 0.85) +
  geom_text(aes(label = paste0("Q=", round(modularity, 3)),
                y = modularity + 0.02), size = 3, angle = 90, hjust = 0) +
  scale_fill_manual(values = c(D = "#2166ac", N = "#b2182b"),
                    labels = c(D = "Day", N = "Night")) +
  labs(x = "Animal", y = "Modularity (Q)",
       title = "Fig S3a: Community structure modularity (Louvain, GLMCC lightON)",
       subtitle = paste0("Mean Q = ", round(mean(modularity_scores$modularity), 3),
                         ", Range = [", round(min(modularity_scores$modularity), 3),
                         ", ", round(max(modularity_scores$modularity), 3), "]"),
       fill = "Phase") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(here("figures", "main", "figS3_community_modularity.png"),
       p_mod, width = 10, height = 5, dpi = 150)

# ─── Fig S3b: Community count vs neuron count ────────────────────
p_counts <- ggplot(modularity_scores,
                   aes(x = n_neurons, y = n_communities, color = day_night)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, color = "grey40") +
  geom_text(aes(label = animal), size = 3, vjust = -1, show.legend = FALSE) +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = "Number of Neurons", y = "Number of Communities",
       title = "Fig S3b: Community count scales with neuron count",
       color = "Phase") +
  theme_minimal(base_size = 12)

ggsave(here("figures", "main", "figS3_community_vs_neurons.png"),
       p_counts, width = 8, height = 5, dpi = 150)

# ─── Fig S3c: Per-animal community-region heatmaps ───────────────
cat("Generating per-animal community plots...\n")

for (animal in animals) {
  animal_comm <- all_communities[all_communities$animal == animal, ]
  if (nrow(animal_comm) == 0) next

  # Build confusion matrix: region × community
  tbl <- table(animal_comm$region, animal_comm$community)
  # Normalize by region
  tbl_norm <- sweep(tbl, 1, rowSums(tbl), "/")
  tbl_norm[is.na(tbl_norm)] <- 0

  # Order regions by total neuron count
  region_order <- names(sort(rowSums(tbl), decreasing = TRUE))
  tbl_norm <- tbl_norm[region_order, , drop = FALSE]

  # Convert to long format
  tbl_long <- as.data.frame(as.table(tbl_norm))
  names(tbl_long) <- c("region", "community", "proportion")

  dn <- daynight[animal]
  n_comm <- ncol(tbl)
  n_neu <- sum(tbl)

  p <- ggplot(tbl_long, aes(x = factor(community), y = region, fill = proportion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1)) +
    labs(x = "Community", y = "Region",
         title = paste0(animal, " (", dn,
                        ") — Community × Region"),
         subtitle = paste0(n_neu, " neurons, ", n_comm,
                          " communities, Q = ",
                          round(modularity_scores$modularity[
                            modularity_scores$animal == animal
                          ], 3)),
         fill = "Proportion") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0, size = 7),
          panel.grid = element_blank())

  ggsave(here("figures", "main",
              paste0("figS3_community_", animal, ".png")),
         p, width = max(6, n_comm * 0.8), height = 4.5, dpi = 150)
}

cat("\nCommunity detection complete.\n")
