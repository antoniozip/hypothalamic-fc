#!/usr/bin/env Rscript
# Region-level aggregation: compute connection density between region pairs,
# correcting for unequal neuron counts across regions.
#
# Density(r,s) = (# significant edges between r,s) / (n_r * n_s)
# Output: results/<estimator>/density_<animal>_<condition>.csv

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
  stop("Usage: Rscript aggregate_density.R --animal <id> --condition <cond> [--estimator glmcc|te]")
}

adj_path <- here("results", estimator, paste0("validated_adj_", animal, "_", condition, ".csv"))
adj <- as.matrix(read.csv(adj_path, header = FALSE))

metrics_path <- here("results", estimator, paste0("metrics_", animal, "_", condition, ".csv"))
if (file.exists(metrics_path)) {
  metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
} else {
  regions_path <- here("data", "processed", "regions.csv")
  reg <- read.csv(regions_path, stringsAsFactors = FALSE)
  reg <- reg[reg$animal == animal, ]
  metrics <- data.frame(neuron_id = seq_len(nrow(adj)), region = reg$region)
}

regions <- unique(metrics$region)
n_regions <- length(regions)
density_mat <- matrix(0, nrow = n_regions, ncol = n_regions)
rownames(density_mat) <- regions
colnames(density_mat) <- regions

for (ri in seq_len(n_regions)) {
  for (rj in seq_len(n_regions)) {
    if (ri == rj) next
    idx_i <- which(metrics$region == regions[ri])
    idx_j <- which(metrics$region == regions[rj])
    n_possible <- length(idx_i) * length(idx_j)
    if (n_possible == 0) next
    submat <- adj[idx_i, idx_j, drop = FALSE]
    n_significant <- sum(submat != 0)
    density_mat[ri, rj] <- n_significant / n_possible
  }
}

out_dir <- here("results", estimator)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, paste0("density_", animal, "_", condition, ".csv"))
write.csv(density_mat, out_path)
cat("Wrote", out_path, "\n")

long <- data.frame(
  animal = animal,
  condition = condition,
  source = rep(regions, each = n_regions),
  target = rep(regions, times = n_regions),
  density = as.vector(density_mat),
  row.names = NULL
)
long <- long[long$source != long$target, ]
long_path <- file.path(out_dir, paste0("density_long_", animal, "_", condition, ".csv"))
write.csv(long, long_path, row.names = FALSE)
