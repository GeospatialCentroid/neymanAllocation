# ==============================================================================
# 04c_compare_sampling_methods.R
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
    dplyr::filter(as.character(MLRA) == as.character(target_mlra))
  
  sys_df <- readr::read_csv(sys_files, show_col_types = FALSE) %>%
    dplyr::filter(as.character(MLRA) == as.character(target_mlra))
  
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
    facet_wrap(~Milestone, scales = "fixed") +
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

summarize_all_mlras <- function(sim_dir) {
  # 1. Load ALL data
  srs_files <- list.files(sim_dir, pattern = srs_pattern, full.names = TRUE)
  sys_files <- list.files(sim_dir, pattern = sys_pattern, full.names = TRUE)
  
  if (length(srs_files) == 0 || length(sys_files) == 0) {
    stop("Missing simulation results for global summary.")
  }
  
  # Load and explicitly filter by ALL_MLRA_IDS
  srs_df <- readr::read_csv(srs_files, show_col_types = FALSE)
  sys_df <- readr::read_csv(sys_files, show_col_types = FALSE)
  
  if (exists("ALL_MLRA_IDS")) {
    srs_df <- srs_df %>% dplyr::filter(as.character(MLRA) %in% as.character(ALL_MLRA_IDS))
    sys_df <- sys_df %>% dplyr::filter(as.character(MLRA) %in% as.character(ALL_MLRA_IDS))
  } else {
    warning("ALL_MLRA_IDS not found. Proceeding without filtering MLRAs.")
  }
  
  # 2. Extract milestones
  srs_milestones <- extract_milestones(srs_df, "Random (SRS)")
  sys_milestones <- extract_milestones(sys_df, "Systematic")
  
  combined_df <- dplyr::bind_rows(srs_milestones, sys_milestones)
  
  # 3. Pivot wider to calculate the exact difference per MLRA/Year
  wide_comparison <- combined_df %>%
    tidyr::pivot_wider(
      names_from = Method,
      values_from = Required_N
    ) 
  
  # 4a. Aggregate across all MLRAs (Grouped by Year AND Milestone)
  yearly_summary <- wide_comparison %>%
    dplyr::group_by(year, Milestone) %>%
    dplyr::summarize(
      Mean_SRS_N = round(mean(`Random (SRS)`, na.rm = TRUE), 1),
      Mean_Sys_N = round(mean(Systematic, na.rm = TRUE), 1),
      SD_Systematic = round(sd(Systematic, na.rm = TRUE), 2),
      Min_Systematic  = min(Systematic, na.rm = TRUE),
      Max_Systematic  = max(Systematic, na.rm = TRUE),
      N_MLRAs    = dplyr::n_distinct(MLRA),
      .groups    = "drop"
    ) %>%
    dplyr::arrange(year, Milestone)
  
  # 4b. Aggregate across ALL Years AND MLRAs (Grouped by Milestone ONLY)
  overall_summary <- wide_comparison %>%
    dplyr::group_by(Milestone) %>%
    dplyr::summarize(
      Mean_SRS_N = round(mean(`Random (SRS)`, na.rm = TRUE), 1),
      Mean_Sys_N = round(mean(Systematic, na.rm = TRUE), 1),
      SD_Systematic = round(sd(Systematic, na.rm = TRUE), 2),
      Min_Systematic  = min(Systematic, na.rm = TRUE),
      Max_Systematic  = max(Systematic, na.rm = TRUE),
      N_MLRAs    = dplyr::n_distinct(MLRA),
      .groups    = "drop"
    ) %>%
    dplyr::arrange(Milestone)
  
  return(list(
    yearly = yearly_summary,
    overall = overall_summary
  ))
}

#' Plot Average Systematic Samples Required per MLRA
#' @param sim_dir Directory containing the simulation results
plot_systematic_mlra_averages <- function(sim_dir) {
  # 1. Load ONLY systematic data
  sys_files <- list.files(sim_dir, pattern = sys_pattern, full.names = TRUE)
  
  if (length(sys_files) == 0) {
    stop("Missing systematic simulation results for plot.")
  }
  
  sys_df <- readr::read_csv(sys_files, show_col_types = FALSE)
  
  # Explicitly filter by ALL_MLRA_IDS
  if (exists("ALL_MLRA_IDS")) {
    sys_df <- sys_df %>% dplyr::filter(as.character(MLRA) %in% as.character(ALL_MLRA_IDS))
  }
  
  # 2. Extract milestones specifically for systematic
  sys_milestones <- extract_milestones(sys_df, "Systematic")
  
  # 3. Filter for only the two requested milestones and calculate the mean per MLRA
  target_milestones <- c("3. 95% Acc + 95% CI") # "2. 95% Accuracy",
  
  plot_data <- sys_milestones %>%
    dplyr::filter(Milestone %in% target_milestones) %>%
    dplyr::group_by(MLRA, Milestone) %>%
    dplyr::summarize(Avg_Required_N = mean(Required_N, na.rm = TRUE), .groups = "drop")
  
  # Ensure MLRA is treated as a factor and ordered logically (numeric if possible)
  plot_data$MLRA <- factor(plot_data$MLRA, levels = sort(as.numeric(unique(as.character(plot_data$MLRA)))))
  
  # 4. Generate the Bar Chart
  p <- ggplot(plot_data, aes(x = MLRA, y = Avg_Required_N, fill = Milestone)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.85) +
    facet_wrap(~Milestone, scales = "fixed", ncol = 1) +
    scale_fill_viridis_d(
      option = "mako",
      begin = 0.3,
      end = 0.8
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Average Required Systematic Samples per Target MLRA",
      subtitle = "Values averaged across all simulated years",
      x = "MLRA ID",
      y = "Average Required Sample Size (N)"
    ) +
    theme(
      legend.position = "none", # Legend not needed since the facets label the milestones
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(p)
}

# --- 4. EXECUTION & OUTPUT ----------------------------------------------------

# Note: We assume ALL_MLRA_IDS has been defined by sourcing "scripts/00_config.R".
# You can test a single MLRA here (ideally one that exists in your config list)
my_mlra <- "87"

results <- compare_methods(
  sim_dir = OUTPUT_SIM_DIR,
  target_mlra = my_mlra
)

# 1. Display the Plot (Visual Comparison for Target MLRA)
print(results$plot)

# 2. View the Tabular Summary (Exact numbers & efficiency gains for Target MLRA)
cat(paste("\n--- Efficiency Summary Table: MLRA", my_mlra, "---\n"))
print(as.data.frame(results$summary_table))

# 3. Calculate and print the global summary across ALL Configured MLRAs
global_stats <- summarize_all_mlras(sim_dir = OUTPUT_SIM_DIR)

cat("\n====================================================================\n")
cat("GLOBAL EFFICIENCY SUMMARY (BY YEAR) - TARGET MLRAs ONLY\n")
cat("====================================================================\n")
print(as.data.frame(global_stats$yearly))
cat("====================================================================\n\n")

cat("\n====================================================================\n")
cat("GLOBAL EFFICIENCY SUMMARY (ALL YEARS COMBINED) - TARGET MLRAs ONLY\n")
cat("====================================================================\n")
print(as.data.frame(global_stats$overall))
cat("====================================================================\n\n")

# 4. Display the Bar Chart of Average Required Systematic Samples per Target MLRA
avg_plot <- plot_systematic_mlra_averages(sim_dir = OUTPUT_SIM_DIR)
print(avg_plot)