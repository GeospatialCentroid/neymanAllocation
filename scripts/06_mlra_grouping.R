# ==============================================================================
# 06_mlra_grouping.R
# Purpose: Group MLRAs by their dominant Neyman Allocation predictor
#          based on empirical NLCD thresholds. Outputs a table per year
#          defining which MLRAs get pooled for variance profiling.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(tidyr)
library(readr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------

# TOGGLE: Set to TRUE for UNET processing, FALSE for original processing
USE_UNET <- TRUE

if (USE_UNET) {
  SUMMARY_DIR <- file.path(DERIVED_DIR, "summaries_unet")
  INPUT_FILE <- file.path(SUMMARY_DIR, "ALL_MLRA_UNET_summary_stats.csv")

  # Create a new directory specifically for the variance profiling workflow
  OUTPUT_DIR <- file.path(DERIVED_DIR, "variance_profiling_unet")
  out_name <- "UNET_MLRA_Neyman_Groups.csv"
} else {
  SUMMARY_DIR <- file.path(DERIVED_DIR, "summaries")
  INPUT_FILE <- file.path(SUMMARY_DIR, "ALL_MLRA_summary_stats.csv")

  OUTPUT_DIR <- file.path(DERIVED_DIR, "variance_profiling")
  out_name <- "MLRA_Neyman_Groups.csv"
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Land cover columns to check for the "Otherwise" condition
NLCD_COLS <- c(
  "Water",
  "Developed",
  "Barren",
  "Forest",
  "Shrubland",
  "Herbaceous",
  "Cultivated",
  "Wetlands"
)

# --- 2. LOAD DATA -------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(paste("Summary file not found at:", INPUT_FILE))
}

df_summary <- readr::read_csv(INPUT_FILE, show_col_types = FALSE)

# --- 3. APPLY GROUPING LOGIC --------------------------------------------------
message(paste0(
  "Applying Neyman Allocation grouping logic (Mode: ",
  ifelse(USE_UNET, "UNET", "Original"),
  ")..."
))

# Step A: Identify the strictly dominant NLCD class for the "Otherwise" fallback
dominant_lc <- df_summary %>%
  dplyr::select(MLRA_ID, year, dplyr::all_of(NLCD_COLS)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(NLCD_COLS),
    names_to = "Dominant_Class",
    values_to = "Pct"
  ) %>%
  dplyr::group_by(MLRA_ID, year) %>%
  # Grab the single highest percentage class for each MLRA/Year
  dplyr::slice_max(order_by = Pct, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(MLRA_ID, year, Dominant_Class)

# Step B: Apply the hierarchical logic
grouped_mlras <- df_summary %>%
  dplyr::left_join(dominant_lc, by = c("MLRA_ID", "year")) %>%
  dplyr::mutate(
    Target_Neyman_Class = dplyr::case_when(
      Forest >= 1.0 ~ "Forest",
      Forest < 1.0 & Wetlands >= 0.5 ~ "Wetlands",
      Forest < 1.0 & Wetlands < 0.5 & Cultivated >= 50.0 ~ "Cultivated",
      # Otherwise fallback
      TRUE ~ Dominant_Class
    )
  ) %>%
  # Select relevant columns for the summary table
  dplyr::select(
    MLRA_ID,
    year,
    Forest_Pct = Forest,
    Wetlands_Pct = Wetlands,
    Cultivated_Pct = Cultivated,
    Dominant_Class,
    Target_Neyman_Class
  ) %>%
  dplyr::arrange(year, Target_Neyman_Class, MLRA_ID)

# --- 4. EXPORT & SUMMARIZE ----------------------------------------------------

out_path <- file.path(OUTPUT_DIR, out_name)
readr::write_csv(grouped_mlras, out_path)
message(paste("Grouping table saved to:", out_path))

# Print a clean summary to the console so you can see the groups instantly
summary_counts <- grouped_mlras %>%
  dplyr::group_by(year, Target_Neyman_Class) %>%
  dplyr::summarise(
    MLRA_Count = dplyr::n(),
    MLRAs_in_Group = paste(MLRA_ID, collapse = ", "),
    .groups = "drop"
  )

message("\n=== Summary of Aggregation Groups by Year ===")
print(summary_counts)
