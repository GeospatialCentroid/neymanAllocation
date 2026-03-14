# ==============================================================================
# 04_random_sampling_test.R (Weighted Version)
# Purpose: Baseline Stress-Test using Simple Random Sampling (SRS) for all years.
# ==============================================================================

source("scripts/00_config.R")

# --- 1. CONFIGURATION ---------------------------------------------------------

USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results_unet")
  input_suffix <- "_UNET_master_dataset.csv"
  output_suffix <- "_baseline_weighted_results_unet.csv"
  method_name <- "BASELINE_SRS_WEIGHTED_UNET"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results")
  input_suffix <- "_master_dataset.csv"
  output_suffix <- "_baseline_weighted_results.csv"
  method_name <- "BASELINE_SRS_WEIGHTED"
}

if (!dir.exists(OUTPUT_SIM_DIR)) {
  dir.create(OUTPUT_SIM_DIR, recursive = TRUE)
}

SIM_REPS <- 100
SAMPLE_STEP <- 100
ACCURACY_TOLERANCE <- 0.10
ACCURACY_PASS_RATE <- 0.80
COVERAGE_PASS_RATE <- 0.90
AREA_COL <- "grid_area"
YEARS <- c(2010, 2016, 2020)

# --- 2. HELPER FUNCTIONS ------------------------------------------------------

analyze_weighted_sample <- function(
  sample_df,
  N_total,
  area_col = "grid_area"
) {
  n <- nrow(sample_df)
  x <- sample_df[[area_col]]
  y <- sample_df$TOF * x

  R_hat <- sum(y) / sum(x)
  fpc <- 1 - (n / N_total)
  x_bar <- mean(x)
  residuals <- y - (R_hat * x)
  s2_res <- sum(residuals^2) / (n - 1)
  se <- sqrt(fpc / n) * (1 / x_bar) * sqrt(s2_res)

  ci_lower <- R_hat - (1.96 * se)
  ci_upper <- R_hat + (1.96 * se)

  return(c(mean = R_hat, ci_l = ci_lower, ci_u = ci_upper))
}

# --- 3. MAIN SIMULATION LOOP --------------------------------------------------

mlra_ids <- if (!is.null(TARGET_MLRA_IDS)) TARGET_MLRA_IDS else ALL_MLRA_IDS

for (m_id in mlra_ids) {
  message(paste0(
    "\n--- Starting Weighted Random Sampling Baseline: ",
    m_id,
    " ---"
  ))
  data_file <- file.path(INPUT_DYNAMIC_DIR, paste0("MLRA_", m_id, input_suffix))

  if (!file.exists(data_file)) {
    warning(paste("Data missing for MLRA", m_id))
    next
  }

  df_full <- readr::read_csv(data_file, show_col_types = FALSE)
  if (!AREA_COL %in% names(df_full)) {
    stop("Error: Area column not found.")
  }

  sim_results_list <- list()

  for (yr in YEARS) {
    pop_df <- df_full %>% dplyr::filter(year == yr, !is.na(TOF))
    if (nrow(pop_df) == 0) {
      next
    }

    TRUE_MEAN <- weighted.mean(pop_df$TOF, pop_df[[AREA_COL]], na.rm = TRUE)
    TOTAL_POP <- nrow(pop_df)
    TOLERANCE_VAL <- TRUE_MEAN * ACCURACY_TOLERANCE

    message(paste0(
      "   Year: ",
      yr,
      " | Pop: ",
      TOTAL_POP,
      " | Mean: ",
      round(TRUE_MEAN, 2),
      "%"
    ))

    max_n <- TOTAL_POP

    for (n_total in seq(100, max_n, SAMPLE_STEP)) {
      sim_outcomes <- data.frame(
        mean_est = numeric(SIM_REPS),
        covered = logical(SIM_REPS),
        accurate = logical(SIM_REPS)
      )

      for (r in 1:SIM_REPS) {
        set.seed(r)
        drawn_sample <- pop_df %>% dplyr::slice_sample(n = n_total)
        est <- analyze_weighted_sample(drawn_sample, TOTAL_POP, AREA_COL)

        sim_outcomes$mean_est[r] <- est["mean"]
        sim_outcomes$accurate[r] <- abs(est["mean"] - TRUE_MEAN) <=
          TOLERANCE_VAL
        sim_outcomes$covered[r] <- (est["ci_l"] <= TRUE_MEAN) &
          (est["ci_u"] >= TRUE_MEAN)
      }

      accuracy_rate <- mean(sim_outcomes$accurate)
      coverage_rate <- mean(sim_outcomes$covered)
      pass_both <- (accuracy_rate >= ACCURACY_PASS_RATE) &&
        (coverage_rate >= COVERAGE_PASS_RATE)

      sim_results_list[[length(sim_results_list) + 1]] <- data.frame(
        MLRA = m_id,
        year = yr,
        Method_Name = method_name,
        Sample_N = n_total,
        Avg_Estimate = mean(sim_outcomes$mean_est),
        Accuracy_Rate = accuracy_rate,
        Coverage_Rate = coverage_rate,
        Passed_Accuracy_Test = accuracy_rate >= ACCURACY_PASS_RATE,
        Passed_CI_Test = coverage_rate >= COVERAGE_PASS_RATE,
        Passed_Both = pass_both
      )
    }
  }

  if (length(sim_results_list) > 0) {
    final_sim_df <- dplyr::bind_rows(sim_results_list)
    out_file <- file.path(OUTPUT_SIM_DIR, paste0("MLRA_", m_id, output_suffix))
    readr::write_csv(final_sim_df, out_file)
    message(paste("   Saved results to:", out_file))
  }
}
