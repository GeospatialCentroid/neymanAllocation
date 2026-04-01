# ==============================================================================
# Neyman Allocation Sampling Estimates and projection (All Years)
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)

# --- 1. CONFIGURATION ---------------------------------------------------------

USE_UNET <- TRUE

if (USE_UNET) {
  SRS_DIR <- file.path(DERIVED_DIR, "simulation_results_unet")
  NEYMAN_DIR <- file.path(DERIVED_DIR, "method_testing_unet")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "projected_sampleEstimates_unet")
  SUMMARY_DIR <- file.path(DERIVED_DIR, "summaries_unet")

  srs_match_str <- "_baseline_weighted_results_unet.csv"
  neyman_match_str <- "_all_methods_comparison_unet.csv"
  summary_file_name <- "ALL_MLRA_UNET_summary_stats.csv"
  out_prefix <- "UNET_"
} else {
  SRS_DIR <- file.path(DERIVED_DIR, "simulation_results")
  NEYMAN_DIR <- file.path(DERIVED_DIR, "method_testing")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "projected_sampleEstimates")
  SUMMARY_DIR <- file.path(DERIVED_DIR, "summaries")

  srs_match_str <- "_baseline_weighted_results.csv"
  neyman_match_str <- "_all_methods_comparison.csv"
  summary_file_name <- "ALL_MLRA_summary_stats.csv"
  out_prefix <- ""
}

if (!dir.exists(OUTPUT_SIM_DIR)) {
  dir.create(OUTPUT_SIM_DIR, recursive = TRUE)
}

srsSample <- list.files(
  SRS_DIR,
  pattern = paste0("*", srs_match_str, "$"),
  full.names = TRUE
)
neymanAllocations <- list.files(
  NEYMAN_DIR,
  pattern = paste0("*", neyman_match_str, "$"),
  full.names = TRUE
)
summary_file_path <- file.path(SUMMARY_DIR, summary_file_name)

if (!file.exists(summary_file_path)) {
  stop(paste("Summary stats file missing at:", summary_file_path))
}

mlraData <- readr::read_csv(summary_file_path, show_col_types = FALSE) %>%
  dplyr::select(MLRA_ID, year, Total_Area_Ha, TOF)

YEARS <- c(2010, 2016, 2020)
results_list <- list()

for (i in TARGET_MLRA_IDS) {
  message(paste0(
    "Processing MLRA ",
    i,
    " - Mode: ",
    ifelse(USE_UNET, "UNET", "Original")
  ))

  baseline_path <- srsSample[grepl(
    paste0("MLRA_", i, srs_match_str),
    srsSample
  )]
  neyman_path <- neymanAllocations[grepl(
    paste0("MLRA_", i, neyman_match_str),
    neymanAllocations
  )]

  if (length(baseline_path) == 0 || length(neyman_path) == 0) {
    warning(paste0("   Missing required files for MLRA ", i))
    next
  }

  srsInput_full <- readr::read_csv(baseline_path, show_col_types = FALSE)
  neyman_full <- readr::read_csv(neyman_path, show_col_types = FALSE)

  srsInput_full <- readr::read_csv(baseline_path, show_col_types = FALSE)
  neyman_full <- readr::read_csv(neyman_path, show_col_types = FALSE)

  # --- NEW SAFETY CHECKS ---
  if (!"year" %in% names(srsInput_full)) {
    stop(paste(
      "\nERROR: 'year' column missing in:",
      baseline_path,
      "\n-> Please re-run Script 04 for this MLRA."
    ))
  }
  if (!"year" %in% names(neyman_full)) {
    stop(paste(
      "\nERROR: 'year' column missing in:",
      neyman_path,
      "\n-> Please re-run Script 03 for this MLRA."
    ))
  }
  # -------------------------

  for (yr in YEARS) {
    coreStats <- mlraData %>% dplyr::filter(MLRA_ID == i, year == yr)
    if (nrow(coreStats) == 0) {
      next
    }

    srs <- srsInput_full %>%
      dplyr::filter(year == yr, Passed_Both == TRUE) %>%
      dplyr::slice(1)

    if (nrow(srs) == 0) {
      next
    }

    neyman <- neyman_full %>%
      dplyr::filter(year == yr, Method == "zero_kmeans") %>%
      dplyr::arrange(desc(Efficiency_Gain_Pct)) %>%
      dplyr::slice(1)

    if (nrow(neyman) == 0) {
      next
    }

    projected_samples <- ceiling(
      srs$Sample_N * (1 - (neyman$Efficiency_Gain_Pct / 100))
    )

    results_list[[paste(i, yr, sep = "_")]] <- data.frame(
      MLRA_ID = i,
      year = yr,
      Known_TOF = coreStats$TOF,
      Total_Area_Ha = coreStats$Total_Area_Ha,
      est_1km_cells = nrow(srsInput_full %>% dplyr::filter(year == yr)) * 100,
      SRS_Samples = srs$Sample_N,
      Neyman_Class = neyman$Variable,
      Efficiency_Gain_Pct = neyman$Efficiency_Gain_Pct,
      Projected_Neyman_Samples = projected_samples,
      measureSampleDensity = srs$Sample_N / (coreStats$Total_Area_Ha / 100),
      projectedSampleDensity = projected_samples /
        (coreStats$Total_Area_Ha / 100)
    )
  }
}

if (length(results_list) > 0) {
  allResults_df <- dplyr::bind_rows(results_list)
  out_estimates <- file.path(
    OUTPUT_SIM_DIR,
    paste0(out_prefix, "MLRA_sample_estimates.csv")
  )
  readr::write_csv(allResults_df, out_estimates)

  summary_df <- allResults_df %>%
    dplyr::group_by(year, Neyman_Class) %>%
    dplyr::summarise(
      Avg_Efficiency_Gain_Pct = mean(Efficiency_Gain_Pct, na.rm = TRUE),
      Avg_Projected_Samples = mean(Projected_Neyman_Samples, na.rm = TRUE),
      std_Projected_Samples = sd(Projected_Neyman_Samples, na.rm = TRUE),
      conservative_Projected_Samples = Avg_Projected_Samples +
        ifelse(is.na(std_Projected_Samples), 0, std_Projected_Samples),
      Count = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(year, desc(Avg_Efficiency_Gain_Pct))

  out_summary <- file.path(
    OUTPUT_SIM_DIR,
    paste0(out_prefix, "Neyman_Class_summary.csv")
  )
  readr::write_csv(summary_df, out_summary)
  message(paste0("Saved multi-year estimates & summary."))
}
