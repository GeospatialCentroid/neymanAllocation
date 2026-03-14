# ==============================================================================
# 09_test_predicted_allocations.R
# Purpose: Simulates sampling 100 times using the predicted Neyman allocations.
#          Evaluates if the generated sample budgets and borrowed variances
#          successfully hit the +/- 10% accuracy target (80% of the time) AND
#          the 95% Confidence Interval coverage target (90% of the time).
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)
library(tidyr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------

USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling_unet")
  ALLOCATION_DIR <- file.path(DERIVED_DIR, "final_allocations_unet")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "validation_results_unet")

  input_suffix <- "_UNET_master_dataset.csv"
  profile_file <- file.path(
    PROFILING_DIR,
    "UNET_Universal_Variance_Profiles.csv"
  )
  allocation_suffix <- "_UNET_predicted_allocations.csv"
  out_file <- "UNET_Final_Validation_Report.csv"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling")
  ALLOCATION_DIR <- file.path(DERIVED_DIR, "final_allocations")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "validation_results")

  input_suffix <- "_master_dataset.csv"
  profile_file <- file.path(PROFILING_DIR, "Universal_Variance_Profiles.csv")
  allocation_suffix <- "_predicted_allocations.csv"
  out_file <- "Final_Validation_Report.csv"
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Simulation Parameters
SIM_REPS <- 100
ACCURACY_TOLERANCE <- 0.10
ACCURACY_PASS_RATE <- 0.80
COVERAGE_PASS_RATE <- 0.90
AREA_COL <- "grid_area"

# --- 2. HELPER FUNCTION -------------------------------------------------------

#' Calculate Stratified Area-Weighted Mean and CI (Combined Ratio Estimator)
analyze_stratified_weighted_sample <- function(
  sample_df,
  pop_strata_sizes,
  N_total,
  area_col = "grid_area"
) {
  strata_stats <- sample_df %>%
    dplyr::group_by(strata) %>%
    dplyr::summarise(
      n_h = dplyr::n(),
      sum_y_h = sum(TOF * !!rlang::sym(area_col)),
      sum_x_h = sum(!!rlang::sym(area_col)),
      s2_res_h = {
        # Prevent division by zero if sum_x_h is 0
        if (sum_x_h == 0) {
          0
        } else {
          r_h <- sum_y_h / sum_x_h
          y <- TOF * !!rlang::sym(area_col)
          x <- !!rlang::sym(area_col)
          res <- y - (r_h * x)
          ifelse(dplyr::n() > 1, sum(res^2) / (dplyr::n() - 1), 0)
        }
      },
      mean_x_h = mean(!!rlang::sym(area_col)),
      .groups = "drop"
    ) %>%
    dplyr::left_join(pop_strata_sizes, by = "strata")

  # Calculate point estimate
  numerator <- sum(
    (strata_stats$N_h / N_total) * (strata_stats$sum_y_h / strata_stats$n_h),
    na.rm = TRUE
  )
  X_bar_U <- sum(
    (strata_stats$N_h / N_total) * strata_stats$mean_x_h,
    na.rm = TRUE
  )

  if (X_bar_U == 0) {
    return(c(mean = NA, ci_l = NA, ci_u = NA))
  }

  R_hat_st <- numerator / X_bar_U

  # Calculate variance of the combined ratio estimator
  var_R_hat_st <- (1 / (X_bar_U^2)) *
    sum(
      ((strata_stats$N_h / N_total)^2) *
        (1 - (strata_stats$n_h / strata_stats$N_h)) *
        (strata_stats$s2_res_h / strata_stats$n_h),
      na.rm = TRUE
    )

  # Handle potential negative tiny variances due to floating point math
  var_R_hat_st <- max(0, var_R_hat_st)
  se <- sqrt(var_R_hat_st)

  # 95% Confidence Interval
  ci_lower <- R_hat_st - (1.96 * se)
  ci_upper <- R_hat_st + (1.96 * se)

  return(c(mean = R_hat_st, ci_l = ci_lower, ci_u = ci_upper))
}

# --- 3. LOAD DEPENDENCIES -----------------------------------------------------

message(paste0(
  "Starting Validation Simulation (Mode: ",
  ifelse(USE_UNET, "UNET", "Original"),
  ")..."
))

if (!file.exists(profile_file)) {
  stop("Variance profiles missing. Run Script 07 first.")
}
df_profiles <- readr::read_csv(profile_file, show_col_types = FALSE)

mlra_ids <- if (!is.null(TARGET_MLRA_IDS)) TARGET_MLRA_IDS else ALL_MLRA_IDS
results_list <- list()

# --- 4. SIMULATION LOOP -------------------------------------------------------

for (m_id in mlra_ids) {
  # Load the MLRA population data
  data_path <- file.path(INPUT_DYNAMIC_DIR, paste0("MLRA_", m_id, input_suffix))
  alloc_path <- file.path(
    ALLOCATION_DIR,
    paste0("MLRA_", m_id, allocation_suffix)
  )

  if (!file.exists(data_path) || !file.exists(alloc_path)) {
    warning(paste("   Missing data or allocation guide for MLRA", m_id))
    next
  }

  df_mlra <- readr::read_csv(data_path, show_col_types = FALSE)
  df_alloc <- readr::read_csv(alloc_path, show_col_types = FALSE)

  YEARS <- unique(df_alloc$year)

  for (yr in YEARS) {
    message(paste("   Simulating MLRA:", m_id, "| Year:", yr))

    # 1. Isolate the year and target class
    pop_df <- df_mlra %>% dplyr::filter(year == yr, !is.na(TOF))
    alloc_yr <- df_alloc %>% dplyr::filter(year == yr)

    if (nrow(pop_df) == 0 || nrow(alloc_yr) == 0) {
      next
    }

    target_class <- alloc_yr$Neyman_Class[1]
    target_budget <- alloc_yr$Target_Budget[1]

    # 2. Get True Mean
    TRUE_MEAN <- weighted.mean(pop_df$TOF, pop_df[[AREA_COL]], na.rm = TRUE)
    TOLERANCE_VAL <- TRUE_MEAN * ACCURACY_TOLERANCE
    N_total <- nrow(pop_df)

    # 3. Apply the Universal Bins to the Population
    prof <- df_profiles %>%
      dplyr::filter(year == yr, Neyman_Class == target_class) %>%
      dplyr::arrange(strata)
    breaks <- c(-Inf, prof$Max_Boundary)
    breaks[length(breaks)] <- Inf

    pop_df$strata <- as.numeric(as.character(cut(
      pop_df[[target_class]],
      breaks = breaks,
      include.lowest = TRUE,
      labels = prof$strata
    )))

    # Calculate population sizes for each stratum (N_h)
    pop_strata_sizes <- pop_df %>% dplyr::count(strata, name = "N_h")

    # 4. Run 100 Simulations
    accurate_count <- 0
    covered_count <- 0
    sim_estimates <- numeric(SIM_REPS)

    for (r in 1:SIM_REPS) {
      set.seed(r)

      # Draw samples according to Assigned_Samples in the guide
      drawn_sample <- pop_df %>%
        dplyr::left_join(
          alloc_yr %>% dplyr::select(Strata, Assigned_Samples),
          by = c("strata" = "Strata")
        ) %>%
        dplyr::group_by(strata) %>%
        # Use filter/sample to prevent slice_sample parameter errors
        dplyr::filter(
          dplyr::row_number() %in%
            sample(dplyr::n(), min(dplyr::first(Assigned_Samples), dplyr::n()))
        ) %>%
        dplyr::ungroup()

      # Calculate stratified weighted estimate and CI
      est <- analyze_stratified_weighted_sample(
        drawn_sample,
        pop_strata_sizes,
        N_total,
        AREA_COL
      )

      est_mean <- est["mean"]
      sim_estimates[r] <- est_mean

      # Check Accuracy and Coverage bounds
      if (!is.na(est_mean)) {
        if (abs(est_mean - TRUE_MEAN) <= TOLERANCE_VAL) {
          accurate_count <- accurate_count + 1
        }
        if (est["ci_l"] <= TRUE_MEAN && est["ci_u"] >= TRUE_MEAN) {
          covered_count <- covered_count + 1
        }
      }
    }

    # 5. Record Results
    accuracy_rate <- accurate_count / SIM_REPS
    coverage_rate <- covered_count / SIM_REPS

    results_list[[length(results_list) + 1]] <- data.frame(
      MLRA_ID = m_id,
      Year = yr,
      Neyman_Class = target_class,
      Total_Budget_Allocated = target_budget,
      Actual_Samples_Drawn = sum(alloc_yr$Assigned_Samples),
      True_Weighted_Mean = TRUE_MEAN,
      Avg_Simulated_Mean = mean(sim_estimates, na.rm = TRUE),
      Accuracy_Pass_Rate = accuracy_rate,
      Coverage_Pass_Rate = coverage_rate,
      Passed_Validation = (accuracy_rate >= ACCURACY_PASS_RATE) &&
        (coverage_rate >= COVERAGE_PASS_RATE)
    )
  }
}

# --- 5. EXPORT RESULTS --------------------------------------------------------

if (length(results_list) > 0) {
  final_report <- dplyr::bind_rows(results_list)
  out_path <- file.path(OUTPUT_DIR, out_file)
  readr::write_csv(final_report, out_path)

  message(paste("\nValidation complete. Report saved to:", out_path))

  # Print quick summary
  pass_count <- sum(final_report$Passed_Validation)
  total_count <- nrow(final_report)
  message(sprintf(
    "Overall Success: %d out of %d configurations (%.1f%%) passed BOTH the 80%% Accuracy AND 90%% Coverage thresholds.",
    pass_count,
    total_count,
    (pass_count / total_count) * 100
  ))
} else {
  warning("No validations were successfully completed.")
}
