#!/usr/bin/env Rscript
# Compute per-animal graph metrics from validated adjacency matrices.
# Output: results/<estimator>/metrics_<animal>_<condition>.csv
#
# Usage: Rscript compute_graph_metrics.R --animal day1 --condition lightON --estimator glmcc

library(igraph)
library(here)

args <- commandArgs(trailingOnly = TRUE)
animal <- NULL
condition <- NULL
estimator <- "glmcc"

for (i in seq_along(args)) {
  if (args[i] == "--animal") animal <- args[i + 1]
  if (args[i] == "--condition") condition <- args[i + 1]
  if (args[i] == "--estimator") estimator <- args[i + 1]
}

if (is.null(animal) || is.null(condition)) {
  stop("Usage: Rscript compute_graph_metrics.R --animal <id> --condition <cond> [--estimator glmcc|te]")
}

adj_path <- here("results", estimator, paste0("validated_adj_", animal, "_", condition, ".csv"))
adj <- as.matrix(read.csv(adj_path, header = FALSE))
g <- graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE, diag = FALSE)

n_nodes <- vcount(g)
neuron_id <- seq_len(n_nodes)

strength <- strength(g, mode = "all")
clustering <- transitivity(g, type = "local")
loc_eff <- numeric(n_nodes)
for (v in seq_len(n_nodes)) {
  neighbors_v <- neighbors(g, v, mode = "all")
  if (length(neighbors_v) > 0) {
    sg <- induced_subgraph(g, neighbors_v)
    E(sg)$weight <- abs(E(sg)$weight)
    d <- distances(sg, weights = E(sg)$weight)
    d[d == 0] <- NA
    loc_eff[v] <- 1 / mean(1 / d, na.rm = TRUE)
    if (is.infinite(loc_eff[v]) || is.na(loc_eff[v])) loc_eff[v] <- 0
  }
}
clustering[is.na(clustering)] <- 0
loc_eff[is.na(loc_eff)] <- 0

hub_score <- (hits_scores(g))$hub
z_score <- (strength - mean(strength)) / sd(strength)
z_score[is.na(z_score)] <- 0

regions_path <- here("data", "processed", "neurons.csv")
if (file.exists(regions_path)) {
  neurons_df <- read.csv(regions_path, stringsAsFactors = FALSE)
  neurons_df <- neurons_df[neurons_df$animal == animal, ]
  region_vec <- rep("unknown", n_nodes)
  for (row_idx in seq_len(nrow(neurons_df))) {
    nid <- neurons_df$neuron_id[row_idx]
    if (nid >= 1 && nid <= n_nodes) {
      region_vec[nid] <- neurons_df$region6[row_idx]
    }
  }
} else {
  region_vec <- rep("unknown", n_nodes)
}

out_df <- data.frame(
  animal = animal,
  condition = condition,
  neuron_id = neuron_id,
  region = region_vec,
  node_strength = strength,
  clustering_coefficient = clustering,
  local_efficiency = loc_eff,
  hub_score = hub_score,
  z_score = z_score,
  row.names = NULL
)

out_dir <- here("results", estimator)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0("metrics_", animal, "_", condition, ".csv"))
write.csv(out_df, out_path, row.names = FALSE)
cat("Wrote", out_path, "\n")
