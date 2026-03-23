# ==============================================================================
# main.R
# Purpose: Master execution script for the Neyman Allocation TOF Pipeline.
#          Controls global configurations and sequential script execution.
# ==============================================================================

# --- 1. GLOBAL CONFIGURATION --------------------------------------------------

# Set working directory (Ensure this points to your project root)
setwd("~/trueNAS/work/neymanSampling")

# SET YOUR EXTENT ("NEBRASKA" or "GREAT_PLAINS")
PROCESSING_EXTENT <- "GREAT_PLAINS"
## Define which LRRs to include from the lower 48m, does not apply to Nebraska processing
TARGET_LRR <- c("F", "G", "H")

# Ensure the config file is loaded first
source("scripts/00_config.R")

# --- 2. EXECUTION TOGGLES -----------------------------------------------------
# Set these to TRUE to execute a step, or FALSE to skip it.
# This allows you to re-run later steps without waiting for earlier processing.
RUN_00a_GENERATE_GRIDS <- FALSE # Generate 1km grids from 100km parents (only needs to be done once)
RUN_01_STATIC_PROCESSING <- TRUE # Process base NLCD/Grid data
RUN_02_DYNAMIC_PROCESSING <- FALSE # Extract TOF from model tiles
RUN_03_STRATIFICATION_TESTING <- FALSE # Evaluate all stratification methods
RUN_04_BASELINE_SAMPLING <- FALSE # SRS baseline testing
# RUN_04_ITERATIVE_SAMPLING <- FALSE # (Optional) Detailed Neyman testing
RUN_05_GATHER_RESULTS <- FALSE # Calculate Sample Budgets
RUN_06_MLRA_GROUPING <- FALSE # Apply hierarchical logic for pooling
RUN_07_VARIANCE_PROFILING <- FALSE # Generate Universal Variance Profiles
RUN_08_PREDICT_ALLOCATIONS <- FALSE # Generate Nebraska sampling guides
RUN_09_TEST_PREDICTIONS <- FALSE # 100x Simulation on Nebraska MLRAs

# NEW: Validating on a completely novel, un-trained MLRA
RUN_10_NOVEL_VALIDATION <- TRUE

# --- 3. SCRIPT EXECUTION ------------------------------------------------------

message("=== Starting Neyman Allocation Pipeline ===")
start_time <- Sys.time()
if (RUN_00a_GENERATE_GRIDS) {
  message("\n>>> Running 00a_generate_1km_grids.R")
  source("scripts/00a_generate_1kmGrids.R")
}

if (RUN_01_STATIC_PROCESSING) {
  message("\n>>> Running 01_static_processing.R")
  source("scripts/01_static_processing.R")
}

if (RUN_02_DYNAMIC_PROCESSING) {
  message("\n>>> Running 02_dynamic_processing.R")
  source("scripts/02_dynamic_processing.R")
}

if (RUN_03_STRATIFICATION_TESTING) {
  message("\n>>> Running 03_stratification_testing.R")
  source("scripts/03_stratification_testing.R")
}

if (RUN_04_BASELINE_SAMPLING) {
  message("\n>>> Running 04_random_sampling_test.R")
  source("scripts/04_random_sampling_test.R")
}

if (RUN_04_ITERATIVE_SAMPLING) {
  message("\n>>> Running 04_iterative_sampling_test.R")
  source("scripts/04_iterative_sampling_test.R")
}

if (RUN_05_GATHER_RESULTS) {
  message("\n>>> Running 05_gatherResults.R")
  # Note: Generate Summary Stats usually runs before this step
  source("scripts/07_generate_summary_stats.R")
  source("scripts/05_gatherResults.R")
}

if (RUN_06_MLRA_GROUPING) {
  message("\n>>> Running 06_mlra_grouping.R")
  source("scripts/06_mlra_grouping.R")
}

if (RUN_07_VARIANCE_PROFILING) {
  message("\n>>> Running 07_variance_profiling.R")
  source("scripts/07_variance_profiling.R")
}

if (RUN_08_PREDICT_ALLOCATIONS) {
  message("\n>>> Running 08_predict_allocations.R")
  source("scripts/08_predict_allocations.R")
}

if (RUN_09_TEST_PREDICTIONS) {
  message("\n>>> Running 09_test_predicted_allocations.R")
  source("scripts/09_test_predicted_allocations.R")
}

if (RUN_10_NOVEL_VALIDATION) {
  message("\n>>> Running 10_novel_mlra_validation.R")
  source("scripts/10_novel_mlra_validation.R")
}

# --- 4. FINISH ----------------------------------------------------------------

end_time <- Sys.time()
run_time <- round(difftime(end_time, start_time, units = "mins"), 2)
message(sprintf(
  "\n=== Pipeline Complete. Total Time: %s minutes ===",
  run_time
))
