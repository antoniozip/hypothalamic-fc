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
  "day1","night1","night2","night3","night4","night5",
  "day2","night6","night7","day3","day4","day5","day6"
)
daynight <- c(
  "day1"="D","night1"="N","night2"="N","night3"="N",
  "night4"="N","night5"="N","day2"="D","night6"="N",
  "night7"="N","day3"="D","day4"="D","day5"="D","day6"="D"
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

# ─── Helper: draw single circos with legends ────────────────────
draw_circos <- function(mat, title_text, out_path, color_map = NULL) {
  regions <- rownames(mat)
  if (is.null(color_map)) {
    color_map <- setNames(
      RColorBrewer::brewer.pal(min(length(regions), 9), "Set3")[seq_along(regions)],
      regions
    )
  }

  # Symmetrize the matrix
  mat_sym <- (mat + t(mat)) / 2
  diag(mat_sym) <- 0

  # Density range for colorbar
  density_vals <- mat_sym[upper.tri(mat_sym)]
  density_vals <- density_vals[density_vals > 0]
  max_dens <- if (length(density_vals) > 0) max(density_vals) else 0.1

  # Open PNG with extra width for legend
  png(out_path, width = 1600, height = 1200, res = 150)
  
  # Layout: left 70% for circos, right 30% for legend
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  # --- Panel 1: Circos ---
  par(mar = c(1, 1, 3, 1))
  circos.clear()
  circos.par(
    gap.after = c(rep(2, length(regions) - 1), 10),
    cell.padding = c(0.02, 0, 0.02, 0)
  )

  chordDiagram(mat_sym,
               grid.col = color_map,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               annotationTrack = c("name", "grid"),
               annotationTrackHeight = c(0.05, 0.05),
               preAllocateTracks = list(track.height = 0.1))
  title(title_text, cex = 0.8, line = 1)
  circos.clear()

  # --- Panel 2: Legends ---
  par(mar = c(14, 2, 4, 2))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))

  # -- Region color legend --
  y_top <- 0.95
  y_step <- 0.06
  n_regions <- length(regions)
  
  text(0.1, y_top + 0.03, "Brain Regions", pos = 4, font = 2, cex = 1.0, col = "#333333")
  
  for (i in seq_along(regions)) {
    y_pos <- y_top - (i - 1) * y_step
    region_name <- regions[i]
    region_color <- color_map[region_name]
    
    # Color swatch
    rect(0.05, y_pos - 0.02, 0.25, y_pos + 0.02,
         col = region_color, border = "grey40", lwd = 0.8)
    # Region name
    text(0.30, y_pos, region_name, pos = 4, cex = 0.8, col = "#222222")
  }

  # -- Density colorbar --
  y_bar_bottom <- y_top - n_regions * y_step - 0.10
  y_bar_height <- 0.03
  
  text(0.1, y_bar_bottom + 0.07, "Connection Density",
       pos = 4, font = 2, cex = 0.9, col = "#333333")

  # Gradient bar
  n_steps <- 100
  bar_left <- 0.08
  bar_width <- 0.65
  bar_y <- y_bar_bottom
  
  density_colors <- colorRampPalette(c("#f7fbff", "#6baed6", "#08519c"))(n_steps)
  step_width <- bar_width / n_steps
  for (j in seq_len(n_steps)) {
    rect(bar_left + (j - 1) * step_width, bar_y - y_bar_height/2,
         bar_left + j * step_width, bar_y + y_bar_height/2,
         col = density_colors[j], border = NA)
  }
  rect(bar_left, bar_y - y_bar_height/2,
       bar_left + bar_width, bar_y + y_bar_height/2,
       col = NA, border = "grey40", lwd = 0.8)

  # Density labels
  text(bar_left, bar_y - 0.05, "0", cex = 0.6, col = "#555555")
  text(bar_left + bar_width, bar_y - 0.05,
       sprintf("%.2f", round(max_dens, 2)), cex = 0.6, col = "#555555")
  text(bar_left + bar_width/2, bar_y - 0.05,
       "density", cex = 0.7, col = "#555555")

  # -- Arrow direction legend --
  y_arrow <- y_bar_bottom - 0.15
  text(0.1, y_arrow + 0.03, "Arrow Direction",
       pos = 4, font = 2, cex = 0.8, col = "#333333")
  text(0.35, y_arrow, "→ outgoing from source region",
       cex = 0.65, col = "#666666")
  text(0.35, y_arrow - 0.05, "Link width ∝ connection density",
       cex = 0.65, col = "#666666")

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
