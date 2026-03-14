# ==============================================================================
# 08_statewide_model_comparison.R
# Purpose: Generate state-level summaries comparing UNET vs Original Models.
#          Visualizes Area-Weighted TOF across years and top Neyman
#          Allocation predictors for slide decks.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------
OUTPUT_DIR <- file.path(DERIVED_DIR, "presentation_visuals")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# File Paths (Adjust if your summary names differ slightly)
orig_stats_path <- file.path(
  DERIVED_DIR,
  "summaries",
  "ALL_MLRA_summary_stats.csv"
)
unet_stats_path <- file.path(
  DERIVED_DIR,
  "summaries_unet",
  "ALL_MLRA_UNET_summary_stats.csv"
)

orig_neyman_path <- file.path(
  DERIVED_DIR,
  "projected_sampleEstimates",
  "Neyman_Class_summary.csv"
)
unet_neyman_path <- file.path(
  DERIVED_DIR,
  "projected_sampleEstimates_unet",
  "UNET_Neyman_Class_summary.csv"
)

# --- 2. PART 1: AREA-WEIGHTED TOF OVER TIME -----------------------------------
message("Processing Area-Weighted TOF...")

df_orig <- read_csv(orig_stats_path, show_col_types = FALSE) %>%
  mutate(Model = "RF")
df_unet <- read_csv(unet_stats_path, show_col_types = FALSE) %>%
  mutate(Model = "UNET")

combined_stats <- bind_rows(df_orig, df_unet)

state_tof <- combined_stats %>%
  filter(!is.na(TOF), !is.na(Total_Area_Ha)) %>%
  group_by(year, Model) %>%
  summarise(
    Total_State_Area_Ha = sum(Total_Area_Ha),
    State_Weighted_TOF = sum(TOF * Total_Area_Ha) / sum(Total_Area_Ha),
    .groups = "drop"
  )

# Export Summary Data
write_csv(state_tof, file.path(OUTPUT_DIR, "State_Weighted_TOF_Summary.csv"))
print(state_tof)

# Plot: TOF Variability over years
p_tof <- ggplot(
  state_tof,
  aes(x = factor(year), y = State_Weighted_TOF, fill = Model)
) +
  geom_col(position = "dodge", color = "black", alpha = 0.8) +
  geom_text(
    aes(label = sprintf("%.2f%%", State_Weighted_TOF)),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c("RF" = "#377EB8", "UNET" = "#4DAF4A")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "State-wide TOF Coverage (2010 - 2020)",
    subtitle = "Comparing Area-Weighted TOF (%) between RF and UNET Models with 10m Resolution",
    x = "Year",
    y = "Area-Weighted TOF (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  file.path(OUTPUT_DIR, "State_Weighted_TOF_Comparison.png"),
  p_tof,
  width = 8,
  height = 6,
  bg = "white"
)


# --- 3. PART 2: TOP LAND COVER PREDICTORS -------------------------------------
message("Processing Top Neyman Predictors...")

orig_neyman <- read_csv(orig_neyman_path, show_col_types = FALSE) %>%
  mutate(Model = "RF")
unet_neyman <- read_csv(unet_neyman_path, show_col_types = FALSE) %>%
  mutate(Model = "UNET")

combined_neyman <- bind_rows(orig_neyman, unet_neyman)

# Plot: Efficiency Gain by Predictor
p_neyman <- ggplot(
  combined_neyman,
  aes(
    x = reorder(Neyman_Class, Avg_Efficiency_Gain_Pct),
    y = Avg_Efficiency_Gain_Pct,
    fill = Model
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    color = "black",
    alpha = 0.8,
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(values = c("RF" = "#377EB8", "UNET" = "#4DAF4A")) +
  geom_text(
    aes(label = sprintf("%.1f%% (n=%d)", Avg_Efficiency_Gain_Pct, Count)),
    position = position_dodge(width = 0.8),
    hjust = -0.1,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Top Stratification Predictors by Model",
    subtitle = "Average Efficiency Gain (%) over Simple Random Sampling (n = MLRA Count)",
    x = "Land Cover Class (Dominant Variable)",
    y = "Average Efficiency Gain (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  file.path(OUTPUT_DIR, "Neyman_Predictors_Comparison.png"),
  p_neyman,
  width = 9,
  height = 6,
  bg = "white"
)

message("All presentation visual summaries saved to: ", OUTPUT_DIR)
