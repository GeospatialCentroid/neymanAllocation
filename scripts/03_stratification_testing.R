# ==============================================================================
# 03_stratification_testing.R
# Purpose: Rank and export ALL tested stratification methods per MLRA per Year.
# ==============================================================================
setwd("~/trueNAS/work/neymanSampling")
source("scripts/00_config.R")
source("src/neymanHelperFunctions.R")

# --- 1. SETUP -----------------------------------------------------------------

# TOGGLE: Set to TRUE for UNET processing, FALSE for original processing
USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  OUTPUT_TEST_DIR <- file.path(DERIVED_DIR, "method_testing_unet")
  input_suffix <- "_UNET_master_dataset.csv"
  output_suffix <- "_all_methods_comparison_unet.csv"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  OUTPUT_TEST_DIR <- file.path(DERIVED_DIR, "method_testing")
  input_suffix <- "_master_dataset.csv"
  output_suffix <- "_all_methods_comparison.csv"
}

if (!dir.exists(OUTPUT_TEST_DIR)) {
  dir.create(OUTPUT_TEST_DIR, recursive = TRUE)
}

STRAT_VARS <- c(
  "riparian_pct",
  "Forest",
  "Cultivated",
  "Wetlands",
  "Water",
  "Developed",
  "Barren",
  "Shrubland",
  "Herbaceous"
)
METHODS <- c("kmeans", "zero_kmeans", "quantile", "zero_quantile")
STRATA_COUNTS <- c(3, 4, 5)
YEARS <- c(2010, 2016, 2020)

# --- 2. CORE METRIC FUNCTION --------------------------------------------------

calculate_theoretical_efficiency <- function(df) {
  stats <- df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      Nh = n(),
      Sh = sd(TOF, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Sh = ifelse(is.na(Sh) | Sh == 0, 1e-6, Sh))

  return(sum(stats$Nh * stats$Sh))
}

# --- 3. PROCESSING LOOP -------------------------------------------------------

mlra_ids <- if (!is.null(TARGET_MLRA_IDS)) TARGET_MLRA_IDS else ALL_MLRA_IDS

for (m_id in mlra_ids) {
  message(paste0("\n--- Analyzing All Strategies for MLRA: ", m_id, " ---"))
  message(paste0("    Mode: ", ifelse(USE_UNET, "UNET", "Original")))

  input_file <- file.path(
    INPUT_DYNAMIC_DIR,
    paste0("MLRA_", m_id, input_suffix)
  )

  if (!file.exists(input_file)) {
    warning(paste("Input file missing for MLRA", m_id, "at path:", input_file))
    next
  }

  df_full <- readr::read_csv(input_file, show_col_types = FALSE)
  results_list <- list()

  for (yr in YEARS) {
    df_test <- df_full %>% dplyr::filter(year == yr, !is.na(TOF))

    if (nrow(df_test) == 0) {
      next
    }

    baseline_score <- calculate_theoretical_efficiency(
      df_test %>% dplyr::mutate(strata = 1)
    )

    for (var_name in STRAT_VARS) {
      for (meth in METHODS) {
        for (k in STRATA_COUNTS) {
          df_strat <- apply_stratification(df_test, var_name, meth, k)
          if (is.null(df_strat)) {
            next
          }

          df_strat <- merge_small_strata(df_strat, var_name, min_size = 50)
          current_score <- calculate_theoretical_efficiency(df_strat)
          gain_pct <- (baseline_score - current_score) / baseline_score * 100

          results_list[[paste(var_name, meth, k, yr, sep = "_")]] <- data.frame(
            MLRA = m_id,
            year = yr,
            Variable = var_name,
            Method = meth,
            K = k,
            Theoretical_Score = current_score,
            Efficiency_Gain_Pct = round(gain_pct, 2)
          )
        }
      }
    }
  }

  if (length(results_list) > 0) {
    all_results_df <- dplyr::bind_rows(results_list) %>%
      dplyr::arrange(year, desc(Efficiency_Gain_Pct))

    out_path <- file.path(OUTPUT_TEST_DIR, paste0("MLRA_", m_id, output_suffix))
    readr::write_csv(all_results_df, out_path)
    message(paste(
      "Exported",
      nrow(all_results_df),
      "stratification results to:",
      out_path
    ))
  }
}
message("\nMethodology testing complete.")
