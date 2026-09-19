#!/usr/bin/env Rscript
# Compute per-animal graph metrics from validated adjacency matrices.
# Output: results/<estimator>/metrics_<animal>_<condition>.csv
#
# Usage: Rscript compute_graph_metrics.R --animal day1 --condition lightON --estimator glmcc

library(igraph)
library(here)
source(here::here("src", "r", "compat.R"))

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
# Local efficiency, Latora & Marchiori (2001): the mean inverse shortest-path
# distance among a node's neighbours. Two corrections over the previous version,
# which were self-cancelling for a uniform complete neighbour subgraph and so went
# unnoticed (see tests/test_local_efficiency.R):
#   1. Edge length is 1/|weight|. Coupling strength is not a distance -- a stronger
#      connection must be a SHORTER path, not a longer one.
#   2. Efficiency is mean(1/d), not 1/mean(1/d). The latter is the harmonic mean of
#      distances, which inverts the metric a second time.
# |weight| is used because GLMCC weights are signed; magnitude carries the coupling
# strength and the sign is an E/I label, not a path length.
loc_eff <- numeric(n_nodes)
for (v in seq_len(n_nodes)) {
  neighbors_v <- neighbors(g, v, mode = "all")
  if (length(neighbors_v) < 2) next          # a pair needs two neighbours
  sg <- induced_subgraph(g, neighbors_v)
  if (ecount(sg) == 0) next                  # no paths among them
  E(sg)$weight <- 1 / abs(E(sg)$weight)
  d <- distances(sg, weights = E(sg)$weight)
  diag(d) <- NA                              # exclude self-pairs
  # Disconnected pairs give d = Inf, hence 1/d = 0, which is the correct
  # contribution: unreachable neighbours add no efficiency.
  loc_eff[v] <- mean(1 / d, na.rm = TRUE)
  if (!is.finite(loc_eff[v])) loc_eff[v] <- 0
}
clustering[is.na(clustering)] <- 0
loc_eff[is.na(loc_eff)] <- 0

hub_score <- (hits_scores(g))$hub
z_score <- (strength - mean(strength)) / sd(strength)
z_score[is.na(z_score)] <- 0

# Region labels are mapped positionally: matrix row i is neuron_id i. That
# holds only while the adjacency covers every unit in neurons.csv. Loaders that
# drop spike-free units silently shift every later index, so a rank mismatch is
# a hard error rather than something to paper over with "unknown".
regions_path <- here("data", "processed", "neurons.csv")
if (file.exists(regions_path)) {
  neurons_df <- read.csv(regions_path, stringsAsFactors = FALSE)
  neurons_df <- neurons_df[neurons_df$animal == animal, ]

  if (nrow(neurons_df) != n_nodes) {
    stop(sprintf(
      paste0("Rank mismatch for %s/%s: adjacency has %d nodes but neurons.csv ",
             "lists %d neurons. Positional region mapping would mislabel every ",
             "neuron after the first missing unit. Re-run validation so the ",
             "matrix covers all units (see scripts/revalidate_all.py)."),
      animal, condition, n_nodes, nrow(neurons_df)
    ))
  }

  missing_ids <- setdiff(seq_len(n_nodes), neurons_df$neuron_id)
  if (length(missing_ids) > 0) {
    stop(sprintf("Animal %s: neurons.csv has no row for neuron_id %s",
                 animal, paste(missing_ids, collapse = ", ")))
  }

  region_vec <- rep("unknown", n_nodes)
  for (row_idx in seq_len(nrow(neurons_df))) {
    nid <- neurons_df$neuron_id[row_idx]
    region_vec[nid] <- neurons_df$region6[row_idx]
  }
} else {
  stop(sprintf("Missing %s - region labels are required", regions_path))
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
