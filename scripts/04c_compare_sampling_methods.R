# ==============================================================================
# 06_compare_sampling_methods.R
# Purpose: Compare the required sample sizes between Random (SRS) and
#          Systematic sampling to evaluate efficiency gains.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# --- 1. CONFIGURATION ---------------------------------------------------------

USE_UNET <- TRUE

if (USE_UNET) {
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results_unet")
  srs_pattern <- "_baseline_weighted_results_unet\\.csv$"
  sys_pattern <- "_systematic_weighted_results_unet\\.csv$"
} else {
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results")
  srs_pattern <- "_baseline_weighted_results\\.csv$"
  sys_pattern <- "_systematic_weighted_results\\.csv$"
}

# --- 2. HELPER FUNCTIONS ------------------------------------------------------

#' Extract First Passing N for Milestones
#' @param df The simulation dataframe
#' @param method_label A clean label for the methodology
extract_milestones <- function(df, method_label) {
  # Milestone 1: 80% Accuracy
  m1 <- df %>%
    dplyr::filter(Passed_Accuracy_80 == TRUE) %>%
    dplyr::group_by(MLRA, year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "1. 80% Accuracy")

  # Milestone 2: 95% Accuracy
  m2 <- df %>%
    dplyr::filter(Passed_Accuracy_95 == TRUE) %>%
    dplyr::group_by(MLRA, year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "2. 95% Accuracy")

  # Milestone 3: 95% Acc + 95% CI
  m3 <- df %>%
    dplyr::filter(Passed_Both == TRUE) %>%
    dplyr::group_by(MLRA, year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "3. 95% Acc + 95% CI")

  dplyr::bind_rows(m1, m2, m3) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Method = method_label) %>%
    dplyr::select(MLRA, year, Method, Milestone, Required_N = Sample_N)
}

# --- 3. MAIN COMPARISON LOGIC -------------------------------------------------

compare_methods <- function(sim_dir, target_mlra) {
  # 1. Load the data
  srs_files <- list.files(sim_dir, pattern = srs_pattern, full.names = TRUE)
  sys_files <- list.files(sim_dir, pattern = sys_pattern, full.names = TRUE)

  if (length(srs_files) == 0 || length(sys_files) == 0) {
    stop(
      "Missing simulation results. Ensure both scripts 04 and 05 have been run."
    )
  }

  srs_df <- readr::read_csv(srs_files, show_col_types = FALSE) %>%
    dplyr::filter(MLRA == target_mlra)

  sys_df <- readr::read_csv(sys_files, show_col_types = FALSE) %>%
    dplyr::filter(MLRA == target_mlra)

  if (nrow(srs_df) == 0 || nrow(sys_df) == 0) {
    stop(paste(
      "Data for MLRA",
      target_mlra,
      "not found in one or both datasets."
    ))
  }

  # 2. Extract milestones for both
  srs_milestones <- extract_milestones(srs_df, "Random (SRS)")
  sys_milestones <- extract_milestones(sys_df, "Systematic")

  # 3. Combine and calculate efficiency
  combined_df <- dplyr::bind_rows(srs_milestones, sys_milestones)

  # Pivot wider to calculate the exact difference
  wide_comparison <- combined_df %>%
    tidyr::pivot_wider(
      names_from = Method,
      values_from = Required_N
    ) %>%
    dplyr::mutate(
      Samples_Saved = `Random (SRS)` - Systematic,
      Efficiency_Gain_Pct = round((Samples_Saved / `Random (SRS)`) * 100, 1)
    ) %>%
    dplyr::arrange(year, Milestone)

  # 4. Generate Comparative Plot
  p <- ggplot(
    combined_df,
    aes(x = factor(year), y = Required_N, fill = Method)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      alpha = 0.85
    ) +
    facet_wrap(~Milestone, scales = "free_y") +
    scale_fill_manual(
      values = c("Random (SRS)" = "#B0BEC5", "Systematic" = "#20B2AA"),
      name = "Sampling Method"
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Sampling Efficiency Comparison - MLRA", target_mlra),
      subtitle = "Lower bars indicate a more efficient methodology",
      x = "Simulation Year",
      y = "Required Sample Size (N)"
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12)
    )

  return(list(
    plot = p,
    summary_table = wide_comparison,
    raw_data = combined_df
  ))
}

# --- 4. EXECUTION & OUTPUT ----------------------------------------------------

# Pick ONE MLRA to compare
my_mlra <- "78"

results <- compare_methods(
  sim_dir = OUTPUT_SIM_DIR,
  target_mlra = my_mlra
)

# 1. Display the Plot (Visual Comparison)
print(results$plot)

# 2. View the Tabular Summary (Exact numbers & efficiency gains)
cat("\n--- Efficiency Summary Table ---\n")
print(as.data.frame(results$summary_table))
