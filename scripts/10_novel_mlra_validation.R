# ==============================================================================
# 10_novel_mlra_validation.R
# Purpose: Validates predictive sampling methodology on a novel, fully-mapped MLRA.
# ==============================================================================

source("scripts/00_config.R")
source("src/neymanHelperFunctions.R")
library(dplyr)
library(terra)
library(sf)
library(exactextractr)
library(readr)

# --- 1. SETTINGS & PATHS ------------------------------------------------------

NOVEL_MLRA_ID <- "62"
TARGET_YEAR <- 2020

# Inputs
RAW_TILES_DIR <- file.path(
  INPUT_DIR,
  paste0("MLRA", NOVEL_MLRA_ID, "_classified_maps")
)
STATIC_ATTR_FILE <- file.path(
  DERIVED_DIR,
  "static_attributes",
  paste0("MLRA_", NOVEL_MLRA_ID, "_static_attributes.csv")
)

# Models
PROFILING_DIR <- file.path(DERIVED_DIR, "variance_profiling_unet")
SUMMARY_DIR <- file.path(DERIVED_DIR, "projected_sampleEstimates_unet")

profile_file <- file.path(PROFILING_DIR, "UNET_Universal_Variance_Profiles.csv")
budget_file <- file.path(SUMMARY_DIR, "UNET_Neyman_Class_summary.csv")

# Outputs
OUTPUT_DIR <- file.path(DERIVED_DIR, "novel_validation")
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

out_report_file <- file.path(
  OUTPUT_DIR,
  paste0("MLRA_", NOVEL_MLRA_ID, "_Novel_Validation_Report.csv")
)

# Simulation Parameters
SIM_REPS <- 100
ACCURACY_TOLERANCE <- 0.10
ACCURACY_PASS_RATE <- 0.80
COVERAGE_PASS_RATE <- 0.90
AREA_COL <- "grid_area"

# --- 2. EXTRACT TOF FROM RAW TILES & CACHE ------------------------------------
message(paste("\n--- Processing Novel MLRA:", NOVEL_MLRA_ID, "---"))

if (!file.exists(STATIC_ATTR_FILE)) {
  stop("Static attributes missing for MLRA ", NOVEL_MLRA_ID)
}
static_df <- readr::read_csv(STATIC_ATTR_FILE, show_col_types = FALSE)

# Fix NLCD Column Names (Strips 'y2020_' prefix)
static_df <- static_df %>%
  dplyr::rename_with(~ gsub(paste0("^y", TARGET_YEAR, "_"), "", .x))

# Define Cache File Path
CACHE_FILE <- file.path(
  OUTPUT_DIR,
  paste0("MLRA_", NOVEL_MLRA_ID, "_TOF_cache.csv")
)

if (file.exists(CACHE_FILE)) {
  message(
    "--> Cache found! Loading extracted TOF directly (Skipping VRT extraction)..."
  )
  df_mlra <- readr::read_csv(CACHE_FILE, show_col_types = FALSE)
} else {
  message("--> No cache found. Beginning heavy extraction...")

  # Find all TIF tiles
  tiles <- list.files(RAW_TILES_DIR, pattern = "\\.tif$", full.names = TRUE)
  if (length(tiles) == 0) {
    stop("No TIF tiles found in: ", RAW_TILES_DIR)
  }

  message(sprintf(
    "Found %d tiles. Building Virtual Raster (VRT)...",
    length(tiles)
  ))
  vrt_path <- file.path(RAW_TILES_DIR, "temp_mosaic.vrt")
  vrt <- terra::vrt(tiles, vrt_path, overwrite = TRUE)

  # Load 1km Grids for this MLRA
  grid_1km <- terra::vect(STATIC_INPUTS$grid_1km)
  mlra_grids_sf <- sf::st_as_sf(grid_1km[grid_1km$MLRA_ID == NOVEL_MLRA_ID, ])

  # Force polygons to match Raster CRS exactly before extracting
  raster_crs <- terra::crs(vrt)
  if (sf::st_crs(mlra_grids_sf) != sf::st_crs(raster_crs)) {
    message(
      "Aligning CRS of grids to match 1-meter rasters for maximum speed..."
    )
    mlra_grids_sf <- sf::st_transform(mlra_grids_sf, raster_crs)
  }

  message("Extracting TOF coverage to 1km grids (This may take a while)...")
  extracted_frac <- exactextractr::exact_extract(
    vrt,
    mlra_grids_sf,
    fun = "mean",
    progress = TRUE
  )

  # Clean up temporary VRT
  if (file.exists(vrt_path)) {
    file.remove(vrt_path)
  }

  # Build the temporary Master Dataset
  df_mlra <- static_df %>%
    dplyr::mutate(
      id = mlra_grids_sf$id,
      year = TARGET_YEAR,
      TOF = extracted_frac * 100
    ) %>%
    dplyr::filter(!is.na(TOF))

  # Save to cache so we never have to do this again
  message("Saving extraction results to cache for future runs...")
  readr::write_csv(df_mlra, CACHE_FILE)
}


# --- 3. DETERMINE TARGET NEYMAN CLASS -----------------------------------------
message("Determining Target Neyman Class...")

# If the static attributes contain raw NLCD codes instead of grouped names, group them here.
if (!"Forest" %in% names(df_mlra)) {
  df_mlra <- df_mlra %>%
    dplyr::mutate(
      Forest = rowSums(across(any_of(c("41", "42", "43"))), na.rm = TRUE),
      Wetlands = rowSums(across(any_of(c("90", "95"))), na.rm = TRUE),
      Cultivated = rowSums(across(any_of(c("81", "82"))), na.rm = TRUE),
      Water = rowSums(across(any_of(c("11"))), na.rm = TRUE),
      Developed = rowSums(
        across(any_of(c("21", "22", "23", "24"))),
        na.rm = TRUE
      ),
      Barren = rowSums(across(any_of(c("31"))), na.rm = TRUE),
      Shrubland = rowSums(across(any_of(c("52"))), na.rm = TRUE),
      Herbaceous = rowSums(across(any_of(c("71"))), na.rm = TRUE)
    )
}

# Safety check in case some NLCD classes don't exist in this MLRA
expected_classes <- c(
  "Forest",
  "Wetlands",
  "Cultivated",
  "Water",
  "Developed",
  "Barren",
  "Shrubland",
  "Herbaceous"
)
available_classes <- intersect(expected_classes, names(df_mlra))

if (length(available_classes) == 0) {
  stop(
    "CRITICAL ERROR: No NLCD land cover classes could be found or calculated. Please check the column names in your static attributes file."
  )
}

means <- df_mlra %>%
  dplyr::summarise(across(
    dplyr::all_of(available_classes),
    ~ weighted.mean(.x, w = grid_area, na.rm = TRUE)
  ))

dominant_class <- means %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = "Class",
    values_to = "Pct"
  ) %>%
  dplyr::slice_max(order_by = Pct, n = 1, with_ties = FALSE) %>%
  dplyr::pull(Class)

target_class <- dplyr::case_when(
  "Forest" %in% names(means) && means$Forest >= 1.0 ~ "Forest",
  "Wetlands" %in% names(means) && means$Wetlands >= 0.5 ~ "Wetlands",
  "Cultivated" %in% names(means) && means$Cultivated >= 50.0 ~ "Cultivated",
  TRUE ~ dominant_class
)

# --- 4. PREPARE PREDICTION MODELS ---------------------------------------------
df_profiles <- readr::read_csv(profile_file, show_col_types = FALSE) %>%
  dplyr::filter(year == TARGET_YEAR, Neyman_Class == target_class)
df_budgets <- readr::read_csv(budget_file, show_col_types = FALSE) %>%
  dplyr::filter(year == TARGET_YEAR, Neyman_Class == target_class)

target_budget <- ceiling(df_budgets$conservative_Projected_Samples[1])

breaks <- c(-Inf, df_profiles$Max_Boundary)
breaks[length(breaks)] <- Inf

df_mlra$strata <- as.numeric(as.character(cut(
  df_mlra[[target_class]],
  breaks = breaks,
  include.lowest = TRUE,
  labels = df_profiles$strata
)))

pop_strata_sizes <- df_mlra %>% dplyr::count(strata, name = "N_h")

# Calculate Allocations
alloc_table <- df_mlra %>%
  dplyr::group_by(strata) %>%
  dplyr::summarise(N_h_actual = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(
    df_profiles %>% dplyr::select(strata, S_h),
    by = "strata"
  ) %>%
  dplyr::mutate(
    S_h = ifelse(is.na(S_h), 1e-6, S_h),
    Weight = N_h_actual * S_h
  )

sum_weights <- sum(alloc_table$Weight)

alloc_table <- alloc_table %>%
  dplyr::mutate(
    n_h_target = round((Weight / sum_weights) * target_budget),
    n_h_final = pmin(
      N_h_actual,
      pmax(n_h_target, ifelse(N_h_actual >= 2, 2, N_h_actual))
    )
  )

# --- 5. VALIDATION SIMULATION -------------------------------------------------
message(sprintf(
  "Simulating 100 sample draws (Class: %s | Budget: %d)...",
  target_class,
  target_budget
))

TRUE_MEAN <- weighted.mean(df_mlra$TOF, df_mlra[[AREA_COL]], na.rm = TRUE)
TOLERANCE_VAL <- TRUE_MEAN * ACCURACY_TOLERANCE
N_total <- nrow(df_mlra)

accurate_count <- 0
covered_count <- 0
sim_estimates <- numeric(SIM_REPS)

for (r in 1:SIM_REPS) {
  set.seed(r)

  drawn_sample <- df_mlra %>%
    dplyr::left_join(
      alloc_table %>% dplyr::select(strata, n_h_final),
      by = "strata"
    ) %>%
    dplyr::group_by(strata) %>%
    dplyr::filter(
      dplyr::row_number() %in%
        sample(dplyr::n(), min(dplyr::first(n_h_final), dplyr::n()))
    ) %>%
    dplyr::ungroup()

  est <- analyze_stratified_weighted_sample(
    drawn_sample,
    pop_strata_sizes,
    N_total,
    AREA_COL
  )

  est_mean <- est["mean"]
  sim_estimates[r] <- est_mean

  if (!is.na(est_mean)) {
    if (abs(est_mean - TRUE_MEAN) <= TOLERANCE_VAL) {
      accurate_count <- accurate_count + 1
    }
    if (est["ci_l"] <= TRUE_MEAN && est["ci_u"] >= TRUE_MEAN) {
      covered_count <- covered_count + 1
    }
  }
}

# --- 6. EXPORT REPORT ---------------------------------------------------------
accuracy_rate <- accurate_count / SIM_REPS
coverage_rate <- covered_count / SIM_REPS

report <- data.frame(
  MLRA_ID = NOVEL_MLRA_ID,
  Year = TARGET_YEAR,
  Neyman_Class = target_class,
  Total_Budget_Allocated = target_budget,
  Actual_Samples_Drawn = sum(alloc_table$n_h_final),
  True_Weighted_Mean = TRUE_MEAN,
  Avg_Simulated_Mean = mean(sim_estimates, na.rm = TRUE),
  Accuracy_Pass_Rate = accuracy_rate,
  Coverage_Pass_Rate = coverage_rate,
  Passed_Validation = (accuracy_rate >= ACCURACY_PASS_RATE) &&
    (coverage_rate >= COVERAGE_PASS_RATE)
)

readr::write_csv(report, out_report_file)
message("\nValidation complete! Results saved to: ", out_report_file)
print(report)
