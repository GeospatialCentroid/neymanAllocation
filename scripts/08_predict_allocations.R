# ==============================================================================
# 08_predict_allocations.R
# Purpose: Applies the Universal Variance Profiles to distribute a target
#          sample budget across 1km grids in an MLRA using Neyman Allocation.
#          Outputs a stratum-level sampling guide per MLRA.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)
library(tidyr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------

# TOGGLE: Set to TRUE for UNET processing, FALSE for original processing
USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling_unet")
  SUMMARY_DIR <- file.path(DERIVED_DIR, "projected_sampleEstimates_unet")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "final_allocations_unet")

  input_suffix <- "_UNET_master_dataset.csv"
  profile_file <- file.path(
    PROFILING_DIR,
    "UNET_Universal_Variance_Profiles.csv"
  )
  groups_file <- file.path(PROFILING_DIR, "UNET_MLRA_Neyman_Groups.csv")
  budget_file <- file.path(SUMMARY_DIR, "UNET_Neyman_Class_summary.csv")
  out_suffix <- "_UNET_predicted_allocations.csv"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling")
  SUMMARY_DIR <- file.path(DERIVED_DIR, "projected_sampleEstimates")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "final_allocations")

  input_suffix <- "_master_dataset.csv"
  profile_file <- file.path(PROFILING_DIR, "Universal_Variance_Profiles.csv")
  groups_file <- file.path(PROFILING_DIR, "MLRA_Neyman_Groups.csv")
  budget_file <- file.path(SUMMARY_DIR, "Neyman_Class_summary.csv")
  out_suffix <- "_predicted_allocations.csv"
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# --- 2. LOAD DEPENDENCIES -----------------------------------------------------

message(paste0(
  "Loading predictive models (Mode: ",
  ifelse(USE_UNET, "UNET", "Original"),
  ")..."
))

if (!all(file.exists(profile_file, groups_file, budget_file))) {
  stop(
    "Missing prerequisite files. Ensure Scripts 05, 06, and 07 have been run."
  )
}

df_profiles <- readr::read_csv(profile_file, show_col_types = FALSE)
df_groups <- readr::read_csv(groups_file, show_col_types = FALSE)
df_budgets <- readr::read_csv(budget_file, show_col_types = FALSE)

YEARS <- unique(df_groups$year)
mlra_ids <- if (!is.null(TARGET_MLRA_IDS)) {
  TARGET_MLRA_IDS
} else {
  unique(df_groups$MLRA_ID)
}

# --- 3. ALLOCATION LOOP -------------------------------------------------------

for (m_id in mlra_ids) {
  message(paste("\n--- Generating Sampling Guide for MLRA:", m_id, "---"))

  data_path <- file.path(INPUT_DYNAMIC_DIR, paste0("MLRA_", m_id, input_suffix))
  if (!file.exists(data_path)) {
    warning(paste("   Dynamic data missing for MLRA", m_id))
    next
  }

  df_mlra <- readr::read_csv(data_path, show_col_types = FALSE)
  mlra_allocations_list <- list()

  for (yr in YEARS) {
    # 1. Identify Target Class & Budget for this specific Year
    group_info <- df_groups %>% dplyr::filter(MLRA_ID == m_id, year == yr)
    if (nrow(group_info) == 0) {
      next
    }

    target_class <- group_info$Target_Neyman_Class[1]

    budget_info <- df_budgets %>%
      dplyr::filter(year == yr, Neyman_Class == target_class)
    if (nrow(budget_info) == 0) {
      next
    }

    # We use the Conservative Projected Samples as our target effort ($n$)
    target_n <- ceiling(budget_info$conservative_Projected_Samples[1])

    # 2. Extract the Universal Profile for this Year & Class
    prof <- df_profiles %>%
      dplyr::filter(year == yr, Neyman_Class == target_class) %>%
      dplyr::arrange(strata)
    if (nrow(prof) == 0) {
      next
    }

    # 3. Create strict bin boundaries based on the profile
    # Ensures no gaps exist between the min/max of neighboring strata
    breaks <- c(-Inf, prof$Max_Boundary)
    breaks[length(breaks)] <- Inf # Top boundary is strictly open to catch extreme outliers

    # 4. Map the MLRA's 1km grids into the Universal Strata
    df_yr <- df_mlra %>% dplyr::filter(year == yr)

    if (!target_class %in% names(df_yr)) {
      warning(paste("   Column", target_class, "missing in dataset. Skipping."))
      next
    }

    df_yr$predicted_strata <- cut(
      df_yr[[target_class]],
      breaks = breaks,
      include.lowest = TRUE,
      labels = prof$strata
    )

    # 5. Calculate Neyman Allocations
    allocation_table <- df_yr %>%
      dplyr::group_by(predicted_strata) %>%
      dplyr::summarise(
        N_h_actual = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        predicted_strata = as.numeric(as.character(predicted_strata))
      ) %>%
      dplyr::left_join(
        prof %>% dplyr::select(strata, S_h),
        by = c("predicted_strata" = "strata")
      ) %>%
      # Safety catch: if a stratum has grids but no borrowed variance, assign a tiny variance
      dplyr::mutate(
        S_h = ifelse(is.na(S_h), 1e-6, S_h),
        Weight = N_h_actual * S_h
      )

    sum_weights <- sum(allocation_table$Weight)

    # Perform the Neyman distribution
    allocation_table <- allocation_table %>%
      dplyr::mutate(
        Target_Budget = target_n,
        # The core Neyman Formula: n_h = n * (N_h * S_h) / sum(N_i * S_i)
        n_h_target = round((Weight / sum_weights) * target_n),

        # Safety Rules:
        # 1. Never allocate more samples than there are actual grids (N_h_actual).
        # 2. Try to allocate at least 2 samples to prove variance, unless N_h_actual < 2.
        n_h_final = pmin(
          N_h_actual,
          pmax(n_h_target, ifelse(N_h_actual >= 2, 2, N_h_actual))
        ),

        year = yr,
        Neyman_Class = target_class
      ) %>%
      dplyr::select(
        year,
        Neyman_Class,
        Strata = predicted_strata,
        Available_Grids = N_h_actual,
        Borrowed_S_h = S_h,
        Target_Budget,
        Assigned_Samples = n_h_final
      )

    mlra_allocations_list[[
      length(mlra_allocations_list) + 1
    ]] <- allocation_table
    message(sprintf(
      "      Year %d | Class: %s | Budget: %d -> Allocated: %d",
      yr,
      target_class,
      target_n,
      sum(allocation_table$Assigned_Samples)
    ))
  }

  # Export the sampling guide for this MLRA
  if (length(mlra_allocations_list) > 0) {
    final_mlra_allocations <- dplyr::bind_rows(mlra_allocations_list)
    out_path <- file.path(OUTPUT_DIR, paste0("MLRA_", m_id, out_suffix))
    readr::write_csv(final_mlra_allocations, out_path)
  }
}

message(
  "\nPrediction complete. All sampling guides have been generated in: ",
  OUTPUT_DIR
)


# Define the 3 MLRAs you want to look at (replace with your actual IDs)
my_mlras <- c("78", "86", "150")
sim_results <- list.files(OUTPUT_DIR, full.names = TRUE) |>
  readr::read_csv() |>
  dplyr::bind_rows()
