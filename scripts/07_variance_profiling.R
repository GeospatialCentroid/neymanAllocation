# ==============================================================================
# 07_variance_profiling.R
# Purpose: Aggregates MLRAs by their Target Neyman Class and Year, applies
#          zero_kmeans_5 stratification, and calculates the Universal Variance
#          Profile (S_h and class boundaries) for predictive sampling.
# ==============================================================================

source("scripts/00_config.R")
source("src/neymanHelperFunctions.R")
library(dplyr)
library(readr)
library(tidyr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------

# TOGGLE: Set to TRUE for UNET processing, FALSE for original processing
USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling_unet")

  input_suffix <- "_UNET_master_dataset.csv"
  groups_file <- file.path(PROFILING_DIR, "UNET_MLRA_Neyman_Groups.csv")
  out_profile_name <- "UNET_Universal_Variance_Profiles.csv"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling")

  input_suffix <- "_master_dataset.csv"
  groups_file <- file.path(PROFILING_DIR, "MLRA_Neyman_Groups.csv")
  out_profile_name <- "Universal_Variance_Profiles.csv"
}

if (!file.exists(groups_file)) {
  stop(paste("Grouping file not found at:", groups_file))
}

# --- 2. LOAD GROUPS -----------------------------------------------------------

df_groups <- readr::read_csv(groups_file, show_col_types = FALSE)

# Identify all unique Year + Target Neyman Class combinations
grouping_combos <- df_groups %>%
  dplyr::distinct(year, Target_Neyman_Class) %>%
  dplyr::arrange(year, Target_Neyman_Class)

results_list <- list()

# --- 3. AGGREGATE & PROFILE LOOP ----------------------------------------------
message(paste0(
  "Generating Universal Variance Profiles (Mode: ",
  ifelse(USE_UNET, "UNET", "Original"),
  ")..."
))

for (i in 1:nrow(grouping_combos)) {
  yr <- grouping_combos$year[i]
  target_class <- grouping_combos$Target_Neyman_Class[i]

  # Identify which MLRAs belong to this specific pool
  pool_mlras <- df_groups %>%
    dplyr::filter(year == yr, Target_Neyman_Class == target_class) %>%
    dplyr::pull(MLRA_ID)

  message(sprintf(
    "   Pooling %d MLRAs for Year %d | Class: %s",
    length(pool_mlras),
    yr,
    target_class
  ))

  # Load and combine all master datasets for these MLRAs, filtering to the specific year
  pooled_data_list <- list()
  for (m_id in pool_mlras) {
    data_path <- file.path(
      INPUT_DYNAMIC_DIR,
      paste0("MLRA_", m_id, input_suffix)
    )
    if (file.exists(data_path)) {
      df_mlra <- readr::read_csv(data_path, show_col_types = FALSE) %>%
        dplyr::filter(year == yr, !is.na(TOF))

      pooled_data_list[[length(pooled_data_list) + 1]] <- df_mlra
    }
  }

  if (length(pooled_data_list) == 0) {
    warning(paste(
      "      No dynamic data found for",
      target_class,
      "in year",
      yr
    ))
    next
  }

  # Create the massive combined dataframe
  pooled_df <- dplyr::bind_rows(pooled_data_list)

  # Apply standard stratification (zero_kmeans, k=5) on the combined data
  stratified_df <- apply_stratification(
    pooled_df,
    target_class,
    "zero_kmeans",
    5
  )

  if (is.null(stratified_df)) {
    warning(paste(
      "      Stratification failed for",
      target_class,
      "in year",
      yr
    ))
    next
  }

  # Ensure no tiny strata break the math
  stratified_df <- merge_small_strata(
    stratified_df,
    target_class,
    min_size = 50
  )

  # Extract the Variance Profile
  profile <- stratified_df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      N_h = dplyr::n(),
      S_h = sd(TOF, na.rm = TRUE),
      Min_Boundary = min(!!rlang::sym(target_class), na.rm = TRUE),
      Max_Boundary = max(!!rlang::sym(target_class), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      year = yr,
      Neyman_Class = target_class,
      # Safety catch: if variance is exactly 0, assign a tiny number so Neyman math works
      S_h = ifelse(is.na(S_h) | S_h == 0, 1e-6, S_h)
    ) %>%
    dplyr::select(
      year,
      Neyman_Class,
      strata,
      N_h,
      S_h,
      Min_Boundary,
      Max_Boundary
    )

  results_list[[length(results_list) + 1]] <- profile
}

# --- 4. EXPORT RESULTS --------------------------------------------------------

if (length(results_list) > 0) {
  final_profiles <- dplyr::bind_rows(results_list) %>%
    dplyr::arrange(year, Neyman_Class, strata)

  out_path <- file.path(PROFILING_DIR, out_profile_name)
  readr::write_csv(final_profiles, out_path)

  message(paste(
    "\nVariance profiles successfully generated and saved to:",
    out_path
  ))
} else {
  warning("No profiles were successfully generated.")
}
