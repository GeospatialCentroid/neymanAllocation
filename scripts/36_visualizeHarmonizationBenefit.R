# scripts/36_visualizeHarmonizationBenefit.R
# Purpose: Quantify and visualize the reduction in "phantom change" after normalization

pacman::p_load(dplyr, readr, ggplot2, tidyr)

# 1. Load Data
export_dir <- "~/trueNAS/work/neymanSampling/data/products/groundTruthSamples/scenarios_2020"
grid_stats_path <- file.path(export_dir, "grid_level_tof.csv")

if (!file.exists(grid_stats_path)) {
  stop("grid_level_tof.csv not found. Please run scripts/32_summarizeTOFArea.R first.")
}

grid_stats <- read_csv(grid_stats_path)

# 2. Process Data to calculate Delta from Baseline
# We need to find the 'base' year for each grid ID to use as a reference
processed_stats <- grid_stats %>%
  separate(grid_id, into = c("id", "year"), sep = "_", extra = "drop") %>%
  mutate(year = as.numeric(year))

# Calculate baseline area for each ID (using 'base' normalization type)
# Filter for baseline areas that have >10% tree cover
baselines <- processed_stats %>%
  filter(normalization == "base") %>%
  group_by(id) %>%
  summarize(
    baseline_area = mean(grid_tof_area, na.rm = TRUE),
    grid_total_area = mean(grid_total_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter((baseline_area / grid_total_area) > 0.10)

# Join back and calculate absolute error relative to baseline
comparison <- processed_stats %>%
  filter(normalization != "base") %>%
  inner_join(baselines, by = "id") %>%
  mutate(
    abs_error = abs(grid_tof_area - baseline_area),
    rel_error = (abs_error / baseline_area) * 100,
    absolute_change = grid_tof_area - baseline_area,
    rel_change = ((grid_tof_area - baseline_area) / baseline_area) * 100
  ) %>%
  filter(!is.infinite(rel_error), !is.na(rel_error))

# 3. Generate Comparison Plot (Visual 6)
p_error_comparison <- ggplot(comparison, aes(x = normalization, y = rel_error, fill = normalization)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  coord_cartesian(ylim = c(0, 100)) + # Focus on the bulk of the data
  scale_fill_manual(values = c("normalized" = "#2E6F40", "unnormalized" = "#A63603")) +
  labs(
    title = "Impact of Radiometric Harmonization on Classification Stability",
    subtitle = "Comparing relative error (%) against baseline year (>10% tree cover)",
    x = "Processing Method",
    y = "Relative Error (%) from Baseline",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# 4. Generate Absolute Change Plot (Visual 7)
p_abs_change <- ggplot(comparison, aes(x = normalization, y = absolute_change, fill = normalization)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  scale_fill_manual(values = c("normalized" = "#2E6F40", "unnormalized" = "#A63603")) +
  labs(
    title = "Directional Bias in Classification: Normalized vs. Unnormalized",
    subtitle = "Absolute change in TOF area (m²) from baseline (>10% tree cover)",
    x = "Processing Method",
    y = "Absolute Change from Baseline (m²)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# 5. Generate Directional Relative Change Plot (Visual 8)
p_rel_change <- ggplot(comparison, aes(x = normalization, y = rel_change, fill = normalization)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.5) +
  coord_cartesian(ylim = c(-100, 100)) + 
  scale_fill_manual(values = c("normalized" = "#2E6F40", "unnormalized" = "#A63603")) +
  labs(
    title = "Relative Directional Change: Normalized vs. Unnormalized",
    subtitle = "Relative change (%) from baseline (>10% tree cover)",
    x = "Processing Method",
    y = "Relative Change from Baseline (%)",
    fill = "Method"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# 6. Save and Print
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("outputs/figures/normalization_benefit_comparison.png", p_error_comparison, width = 8, height = 6)
ggsave("outputs/figures/normalization_absolute_change.png", p_abs_change, width = 8, height = 6)
ggsave("outputs/figures/normalization_relative_change.png", p_rel_change, width = 8, height = 6)
print("Plots saved to outputs/figures/")

# Summary Stats
summary_stats <- comparison %>%
  group_by(normalization) %>%
  summarize(
    median_rel_error = median(rel_error, na.rm = TRUE),
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    sd_rel_error = sd(rel_error, na.rm = TRUE),
    median_abs_change = median(absolute_change, na.rm = TRUE),
    mean_abs_change = mean(absolute_change, na.rm = TRUE),
    sd_abs_change = sd(absolute_change, na.rm = TRUE),
    median_rel_change = median(rel_change, na.rm = TRUE),
    mean_rel_change = mean(rel_change, na.rm = TRUE),
    sd_rel_change = sd(rel_change, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)
write_csv(summary_stats, file.path(export_dir, "normalization_summary_stats.csv"))
