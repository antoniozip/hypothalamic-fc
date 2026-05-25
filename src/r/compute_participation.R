#!/usr/bin/env Rscript
# Participation coefficient and second-level hub analysis.
#
# Computes participation coefficient (Guimera & Amaral 2005) from
# validated adjacency matrices + community assignments.
#
# P_i = 1 - Σ_c (k_ic / k_i)²
#
# Output:
#   figures/main/figS4_participation_coefficient.png
#   results/glmcc/participation_coefficient.csv

library(here)
library(igraph)
library(ggplot2)
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

# ─── Load community assignments ──────────────────────────────────
comm_file <- here("results", "glmcc", "communities.csv")
if (!file.exists(comm_file)) {
  stop("Run community_detection.R first to generate communities.csv")
}
communities <- read.csv(comm_file, stringsAsFactors = FALSE)

# ─── Compute participation coefficient per animal ────────────────
all_pc <- data.frame()

cat("Computing participation coefficients...\n")

for (animal in animals) {
  adj_path <- here("results", "glmcc",
                   paste0("validated_adj_", animal, "_lightON.csv"))
  if (!file.exists(adj_path)) next

  adj <- as.matrix(read.csv(adj_path, header = FALSE))
  n_units <- nrow(adj)

  # Get community assignments for this animal
  animal_comm <- communities[communities$animal == animal, ]
  if (nrow(animal_comm) == 0) next

  comm_vec <- animal_comm$community
  n_comm <- length(unique(comm_vec))

  # Per-neuron degree within each community
  k_total <- rowSums(adj != 0)  # binary degree
  k_total[k_total == 0] <- 1    # avoid division by zero

  pc_vec <- numeric(n_units)
  for (i in seq_len(n_units)) {
    ci <- comm_vec[i]
    # Degree within own community
    k_within <- sum(adj[i, comm_vec == ci] != 0)
    # Participation: how distributed are connections across communities?
    pc_sum <- 0
    for (c in unique(comm_vec)) {
      k_ic <- sum(adj[i, comm_vec == c] != 0)
      pc_sum <- pc_sum + (k_ic / k_total[i])^2
    }
    pc_vec[i] <- 1 - pc_sum
  }

  # Get region labels
  animal_neurons <- animal_comm
  region_vec <- animal_neurons$region

  pc_df <- data.frame(
    animal = animal,
    day_night = daynight[animal],
    neuron_id = seq_len(n_units),
    region = region_vec,
    community = comm_vec,
    participation_coefficient = pc_vec,
    degree = k_total,
    row.names = NULL
  )
  all_pc <- rbind(all_pc, pc_df)
}

# Save
write.csv(all_pc, here("results", "glmcc", "participation_coefficient.csv"),
          row.names = FALSE)

cat(sprintf("  Computed PC for %d neurons across %d animals\n",
            nrow(all_pc), length(unique(all_pc$animal))))

# ─── Fig S4a: PC distribution by region ──────────────────────────
cat("Generating participation coefficient figures...\n")

all_pc$day_night <- factor(all_pc$day_night, levels = c("D", "N"))
all_pc$animal <- factor(all_pc$animal, levels = animals)

# Region-level PC (mean ± SE)
region_pc <- all_pc %>%
  group_by(region, day_night) %>%
  summarise(
    mean_pc = mean(participation_coefficient, na.rm = TRUE),
    se_pc = sd(participation_coefficient, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)  # at least 3 neurons per group

p1 <- ggplot(region_pc, aes(x = reorder(region, mean_pc), y = mean_pc,
                              fill = day_night)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_pc - se_pc, ymax = mean_pc + se_pc),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c(D = "#2166ac", N = "#b2182b"),
                    labels = c(D = "Day", N = "Night")) +
  labs(x = "Hypothalamic Region", y = "Mean Participation Coefficient (± SE)",
       title = "Fig S4a: Participation coefficient by region and phase",
       subtitle = "PC ≈ 0 = hub within community; PC ≈ 1 = connector across communities",
       fill = "Phase") +
  coord_flip() +
  theme_minimal(base_size = 12)

ggsave(here("figures", "main", "figS4_participation_by_region.png"),
       p1, width = 9, height = 5, dpi = 150)

# ─── Fig S4b: PC vs degree (hub classification) ──────────────────
# Classify nodes: PC > 0.3 = connector hub, PC < 0.3 = provincial hub
all_pc$hub_type <- ifelse(all_pc$participation_coefficient > 0.3,
                          "Connector Hub", "Provincial Hub")

# Normalize degree within animal (z-score)
all_pc <- all_pc %>%
  group_by(animal) %>%
  mutate(degree_z = (degree - mean(degree)) / sd(degree)) %>%
  ungroup()
all_pc$degree_z[is.na(all_pc$degree_z)] <- 0

p2 <- ggplot(all_pc, aes(x = degree_z, y = participation_coefficient,
                          color = day_night)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(D = "#2166ac", N = "#b2182b"),
                     labels = c(D = "Day", N = "Night")) +
  labs(x = "Degree (z-score within animal)", y = "Participation Coefficient",
       title = "Fig S4b: Hub classification (Guimera & Amaral 2005)",
       subtitle = paste0("Provincial hubs: ",
                         round(mean(all_pc$hub_type == "Provincial Hub") * 100, 1),
                         "%, Connector hubs: ",
                         round(mean(all_pc$hub_type == "Connector Hub") * 100, 1), "%"),
       color = "Phase") +
  theme_minimal(base_size = 12)

ggsave(here("figures", "main", "figS4_hub_classification.png"),
       p2, width = 8, height = 6, dpi = 150)

# ─── Print summary ───────────────────────────────────────────────
cat("\nParticipation coefficient summary:\n")
cat(sprintf("  Mean PC: %.3f ± %.3f\n",
            mean(all_pc$participation_coefficient),
            sd(all_pc$participation_coefficient)))
cat(sprintf("  Provincial hubs (PC < 0.3): %.1f%%\n",
            mean(all_pc$hub_type == "Provincial Hub") * 100))
cat(sprintf("  Connector hubs (PC > 0.3): %.1f%%\n",
            mean(all_pc$hub_type == "Connector Hub") * 100))

# Per-region summary
cat("\nPer-region mean PC:\n")
region_summary <- all_pc %>%
  group_by(region) %>%
  summarise(mean_pc = mean(participation_coefficient),
            n = n(), .groups = "drop") %>%
  arrange(desc(mean_pc))
for (i in seq_len(nrow(region_summary))) {
  cat(sprintf("  %-20s: PC=%.3f (n=%d)\n",
              region_summary$region[i],
              region_summary$mean_pc[i],
              region_summary$n[i]))
}

cat("\nDone.\n")
