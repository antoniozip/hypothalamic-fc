#!/usr/bin/env Rscript
# Circos plots of region-pair connectivity density.
#
# Reads density CSV files and generates:
#   1. Per-animal circos plots (13 lightON)
#   2. Day/Night aggregated circos (2×2: Day/Night × LightON/Ongoing)
#
# Output: figures/main/figS2_circos_*.png

library(here)
library(circlize)
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

# ─── Load all density data ───────────────────────────────────────
all_density <- data.frame()
for (animal in animals) {
  fp <- here("results", "glmcc", paste0("density_long_", animal, "_lightON.csv"))
  if (!file.exists(fp)) next
  df <- read.csv(fp, stringsAsFactors = FALSE)
  df$animal <- animal
  df$day_night <- daynight[animal]
  all_density <- rbind(all_density, df)
}

if (nrow(all_density) == 0) {
  cat("No density data found. Run aggregate_density.R first.\n")
  quit(save = "no")
}

# ─── Helper: build adjacency matrix from long-format density ─────
build_matrix <- function(df, region_order = NULL) {
  regions <- sort(unique(c(df$source, df$target)))
  if (!is.null(region_order)) {
    regions <- intersect(region_order, regions)
  }
  n <- length(regions)
  mat <- matrix(0, n, n)
  rownames(mat) <- colnames(mat) <- regions
  for (i in seq_len(nrow(df))) {
    src <- df$source[i]
    tgt <- df$target[i]
    if (src %in% regions && tgt %in% regions) {
      mat[src, tgt] <- df$density[i]
    }
  }
  mat
}

# ─── Helper: draw single circos ──────────────────────────────────
draw_circos <- function(mat, title_text, out_path, color_map = NULL) {
  regions <- rownames(mat)
  if (is.null(color_map)) {
    color_map <- setNames(
      RColorBrewer::brewer.pal(min(length(regions), 9), "Set3")[seq_along(regions)],
      regions
    )
  }

  png(out_path, width = 1200, height = 1200, res = 150)
  par(mar = c(1, 1, 1, 1))

  circos.clear()
  circos.par(gap.after = c(rep(2, length(regions) - 1), 10),
             cell.padding = c(0.02, 0, 0.02, 0))

  # Make symmetric by averaging
  mat_sym <- (mat + t(mat)) / 2
  diag(mat_sym) <- 0

  chordDiagram(mat_sym,
               grid.col = color_map,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = c("name", "grid"),
               annotationTrackHeight = c(0.05, 0.05),
               preAllocateTracks = list(
                 track.height = 0.1
               ))

  title(title_text, cex = 0.9)
  circos.clear()
  dev.off()
  cat("  ", basename(out_path), "\n")
}

# ─── Fig S2a: Per-animal circos (lightON) ────────────────────────
cat("Generating per-animal circos plots...\n")

# Build a consistent region color map across all animals
all_regions <- sort(unique(c(all_density$source, all_density$target)))
n_colors <- length(all_regions)
region_colors <- setNames(
  if (n_colors <= 9) {
    RColorBrewer::brewer.pal(n_colors, "Set3")
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(n_colors)
  },
  all_regions
)

for (animal in animals) {
  animal_dens <- all_density[all_density$animal == animal, ]
  if (nrow(animal_dens) == 0) next

  mat <- build_matrix(animal_dens)
  if (nrow(mat) < 2) next

  dn <- daynight[animal]
  out_path <- here("figures", "main",
                   paste0("figS2_circos_", animal, "_lightON.png"))
  draw_circos(mat,
              paste0(animal, " (", dn, ") — LightON"),
              out_path, region_colors)
}

# ─── Fig S2b: Day/Night aggregated circos (2×2) ─────────────────
cat("Generating Day/Night aggregated circos...\n")

# Aggregate by day_night and condition
agg <- all_density %>%
  group_by(source, target, day_night) %>%
  summarise(density = mean(density, na.rm = TRUE), .groups = "drop")

for (dn in c("D", "N")) {
  agg_dn <- agg[agg$day_night == dn, ]
  if (nrow(agg_dn) == 0) next

  mat <- build_matrix(agg_dn)
  if (nrow(mat) < 2) next

  label <- if (dn == "D") "Day" else "Night"
  out_path <- here("figures", "main",
                   paste0("figS2_circos_aggregated_", label, "_lightON.png"))
  draw_circos(mat,
              paste0(label, " animals — LightON (mean density)"),
              out_path, region_colors)
}

cat("\nAll circos plots generated.\n")
