# ==============================================================================
# 09_test_predicted_allocations.R
# Purpose: Simulates sampling 100 times using the predicted Neyman allocations.
#          Evaluates if the generated sample budgets and borrowed variances
#          successfully hit the +/- 10% accuracy target (80% of the time) AND
#          the 95% Confidence Interval coverage target (90% of the time).
# ==============================================================================

source("scripts/00_config.R")
source("src/neymanHelperFunctions.R")
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
        # Expand the CI boundaries by a microscopic floating-point tolerance
        if (
          (est["ci_l"] - 1e-9) <= TRUE_MEAN && (est["ci_u"] + 1e-9) >= TRUE_MEAN
        ) {
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

library(ggplot2)
library(scales)

plot_validation_success <- function(report_path) {
  # Check if file exists
  if (!file.exists(report_path)) {
    stop("Validation report not found. Please check the file path.")
  }

  # Read the validation data
  df <- readr::read_csv(report_path, show_col_types = FALSE)

  # -------------------------------------------------------------------------
  # PLOT 1: Overall Pass Rate (% of MLRAs that Passed Validation)
  # -------------------------------------------------------------------------
  summary_df <- df %>%
    dplyr::group_by(Year, Neyman_Class) %>%
    dplyr::summarize(
      Total_MLRAs = dplyr::n(),
      Passed = sum(Passed_Validation == TRUE, na.rm = TRUE),
      Pass_Rate = Passed / Total_MLRAs,
      .groups = "drop"
    )

  p1 <- ggplot(
    summary_df,
    aes(x = factor(Year), y = Pass_Rate, fill = Neyman_Class)
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    geom_text(
      aes(label = scales::percent(Pass_Rate, accuracy = 1)),
      position = position_dodge(width = 0.8),
      vjust = -0.5,
      size = 3.5
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.1)) +
    labs(
      title = "Overall Validation Success Rate by Class and Year",
      subtitle = "Percentage of MLRAs meeting BOTH Accuracy and Coverage targets",
      x = "Year",
      y = "Success Rate (% of MLRAs Passed)",
      fill = "Neyman Class"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank())

  # -------------------------------------------------------------------------
  # PLOT 2: Distribution of Accuracy and Coverage Rates across all MLRAs
  # -------------------------------------------------------------------------
  # Reshape data to plot Accuracy and Coverage side-by-side
  df_long <- df %>%
    dplyr::select(
      MLRA_ID,
      Year,
      Neyman_Class,
      Accuracy_Pass_Rate,
      Coverage_Pass_Rate
    ) %>%
    tidyr::pivot_longer(
      cols = c(Accuracy_Pass_Rate, Coverage_Pass_Rate),
      names_to = "Metric",
      values_to = "Rate"
    ) %>%
    dplyr::mutate(
      Metric = dplyr::recode(
        Metric,
        "Accuracy_Pass_Rate" = "Accuracy Rate (Target: >= 80%)",
        "Coverage_Pass_Rate" = "Coverage Rate (Target: >= 90%)"
      )
    )

  # Create threshold dataframe for the dashed lines
  thresholds <- data.frame(
    Metric = c(
      "Accuracy Rate (Target: >= 80%)",
      "Coverage Rate (Target: >= 90%)"
    ),
    yint = c(0.80, 0.90)
  )

  p2 <- ggplot(df_long, aes(x = Neyman_Class, y = Rate, fill = factor(Year))) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15),
      alpha = 0.5,
      size = 1.5,
      color = "darkgray"
    ) +
    facet_wrap(~Metric) +
    geom_hline(
      data = thresholds,
      aes(yintercept = yint),
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    # Removed the limits from scale_y_continuous
    scale_y_continuous(labels = scales::percent_format()) +
    # Added coord_cartesian to visually zoom the y-axis to 75% - 100%
    coord_cartesian(ylim = c(0.75, 1.0)) +
    labs(
      title = "Distribution of Simulation Pass Rates Across MLRAs",
      subtitle = "Red dashed lines indicate the minimum passing threshold",
      x = "Neyman Class",
      y = "Simulation Pass Rate",
      fill = "Year"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")

  # Print the plots to the graphics device
  print(p1)
  print(p2)

  # Optional: Save plots automatically to the output directory
  # out_dir <- dirname(report_path)
  # ggsave(file.path(out_dir, "Validation_Overall_Success_Bar.png"), p1, width = 8, height = 6)
  # ggsave(file.path(out_dir, "Validation_Rate_Distributions.png"), p2, width = 10, height = 6)

  return(invisible(list(summary_plot = p1, distribution_plot = p2)))
}

# --- Example of how to call it at the bottom of script 09 ---
# (Assuming `out_path` is already defined in your environment from Step 5)
plot_validation_success(out_path)


# --- 6. VISUALIZE GENERALIZED SUCCESS SUMMARY ---------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

plot_generalized_success <- function(report_path) {
  # Check if file exists
  if (!file.exists(report_path)) {
    stop("Validation report not found. Please check the file path.")
  }

  df <- readr::read_csv(report_path, show_col_types = FALSE)

  # 1. Aggregate to the MLRA level
  # We average the simulation pass rates across all years/classes for each MLRA
  mlra_summary <- df %>%
    dplyr::group_by(MLRA_ID) %>%
    dplyr::summarize(
      Mean_Accuracy = mean(Accuracy_Pass_Rate, na.rm = TRUE),
      Mean_Coverage = mean(Coverage_Pass_Rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Pass_Accuracy = Mean_Accuracy >= 0.80,
      Pass_Coverage = Mean_Coverage >= 0.90
    )

  # 2. Count the successes dynamically
  total_mlras <- nrow(mlra_summary)
  acc_passed <- sum(mlra_summary$Pass_Accuracy)
  cov_passed <- sum(mlra_summary$Pass_Coverage)

  # 3. Create a summary dataframe for the plot
  plot_data <- data.frame(
    Test = c("Accuracy Requirement (≥ 80%)", "Coverage Requirement (≥ 90%)"),
    Passed = c(acc_passed, cov_passed),
    Failed = c(total_mlras - acc_passed, total_mlras - cov_passed),
    Total = c(total_mlras, total_mlras)
  ) %>%
    # Convert to long format for a stacked bar chart
    tidyr::pivot_longer(
      cols = c(Passed, Failed),
      names_to = "Status",
      values_to = "Count"
    ) %>%
    # Ensure "Passed" is plotted first (on the left)
    dplyr::mutate(Status = factor(Status, levels = c("Failed", "Passed")))

  # 4. Build a generalized, high-impact horizontal bar chart
  p <- ggplot(plot_data, aes(x = Count, y = Test, fill = Status)) +
    geom_bar(stat = "identity", width = 0.4) +

    # Add large text labels inside the "Passed" bars to drive the message home
    geom_text(
      data = plot_data %>% dplyr::filter(Status == "Passed"),
      aes(x = Count / 2, label = paste0(Count, " / ", Total, " MLRAs Passed")),
      color = "white",
      fontface = "bold",
      size = 6
    ) +

    # Clean, professional colors (Green for pass, light gray for fail to keep focus on success)
    scale_fill_manual(values = c("Passed" = "#2CA02C", "Failed" = "#E3E3E3")) +

    # Formatting and aesthetics
    scale_x_continuous(
      breaks = seq(0, total_mlras, by = 5),
      limits = c(0, total_mlras)
    ) +
    labs(
      title = "Universal Variance Profile: Overall MLRA Validation",
      subtitle = "Simulation results testing the Neyman allocation methodology across the landscape",
      x = "Number of MLRAs",
      y = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 15) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(color = "gray40", size = 12),
      panel.grid.major.y = element_blank(), # Remove horizontal grid lines
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(face = "bold", size = 14)
    )

  print(p)

  # Optional: Save it
  # ggsave(file.path(dirname(report_path), "Generalized_Validation_Summary.png"), p, width = 10, height = 5)

  return(invisible(p))
}

# --- Example of how to call it at the bottom of script 09 ---
plot_generalized_success(out_path)
