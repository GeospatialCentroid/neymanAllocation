# ==============================================================================
# Establish_stratifiedGrid.R
# Purpose: Develop the standardized stratified sample of LLR using a specified
#          number of features per MLRA. Supports multiple independent draw sizes.
# Output:  Comprehensive dataframes containing the 1km cell ID, MLRA ID,
#          and assigned LLR ID for all sampled grids, split by draw size.
# ==============================================================================

source("scripts/00_config.R")
source("src/sampleGridsFunctions.R")
library(tmap)
library(dplyr)
library(sf)
library(purrr)
library(readr)

# --- 1. SETUP & LOCAL PATHS ---------------------------------------------------
# Define sampling parameters
TARGET_LLR <- "G"
TARGET_SIZE <- 1400
message("Loading spatial inputs...")

lrr_id_path <- "data/derived/mlra/lower48MLRA.gpkg"
mlra_grid_path <- "data/derived/grids/GreatPlains_1km_mlra.gpkg"
g100 <- sf::st_read("data/raw/grid100km_aea.gpkg")
systematicSamples <- list.files("data/products/systematicSampleSelection/", full.names = TRUE)
nlcd_paths <- list.files("data/raw/nlcd/", full.names = TRUE)

# select out the draw that is specific files with site information 
siteIDs <- systematicSamples[grepl(pattern = paste0(TARGET_LLR,"_draw_", TARGET_SIZE), 
                                   x = systematicSamples)] |>
  read_csv()


# generate the AOI objects for the selected sites 
# --- 2. GENERATE AOIS & BOUNDING BOXES ----------------------------------------
message("Generating AOIs and extracting bounding boxes...")
siteExport <- paste0("data/derived/sampleBBOX/bbox_", TARGET_LLR, "_", TARGET_SIZE,".csv")
if(!file.exists(siteExport)){
  # Ensure unique IDs are extracted from the loaded siteIDs
  # (Note: adjust 'id' if your CSV uses a different column name for the grid cell ID)
  unique_sites <- unique(siteIDs$id) 
  
  # Set up parallel backend to speed up spatial operations
  workers <- max(1, future::availableCores() - 10)
  future::plan(future::multisession, workers = workers)
  
  # Utilize furrr to distribute the bounding box extraction
  bboxes <- furrr::future_map_dfr(unique_sites, function(current_id) {
    
    # 1. Generate AOI using your custom function
    # Note: If getAOI() strictly expects an argument named 'grid100', pass it as such, 
    # or update the function to be more generic (e.g., 'input_grid').
    aoi_native <- getAOI(grid100 = g100, id = current_id)
    
    # 2. Transform the AOI to the target CRS
    # If your upcoming NLCD/TCC rasters are in Albers (EPSG:5070), you might want to 
    # skip this transform or change 4326 to match your raster native CRS to avoid 
    # reprojection on the fly later.
    aoi_proj <- sf::st_transform(aoi_native, crs = 4326)
    
    # 3. Extract bounding box coordinates
    bbox <- sf::st_bbox(aoi_proj)
    
    # 4. Return as a clean, single-row data frame
    data.frame(
      id = current_id,
      xmin = as.numeric(bbox["xmin"]),
      ymin = as.numeric(bbox["ymin"]),
      xmax = as.numeric(bbox["xmax"]),
      ymax = as.numeric(bbox["ymax"]),
      stringsAsFactors = FALSE
    )
  }, .options = furrr::furrr_options(seed = TRUE))
  
  # Safely close the parallel backend
  future::plan(future::sequential)
  
  # Join the BBox coordinates back to the target sites dataframe
  sites_with_bbox <- siteIDs |>
    dplyr::left_join(bboxes, by = "id")
  
  message("Bounding box extraction complete.")
  write_csv(sites_with_bbox, siteExport)
}else{
  sites_with_bbox <- read_csv(siteExport)
}

# --- 3. SPATIAL METRICS EXTRACTION & EXPORT -----------------------------------
message("\n--- Starting NLCD & TCC Processing ---")

# Define the dynamic output path for the final metrics
metricsExport <- paste0("data/derived/metrics/nlcd_tcc_metrics_", TARGET_LLR, "_", TARGET_SIZE, ".csv")

if (!file.exists(metricsExport)) {
  message("No cached metrics found. Beginning parallel extraction...")
  
  # Ensure the output directory exists
  dir.create(dirname(metricsExport), recursive = TRUE, showWarnings = FALSE)
  
  # Separate NLCD and TCC file paths based on naming conventions
  nlcd_files <- list(
    "2010" = nlcd_paths[grepl("nlcd.*2010", nlcd_paths, ignore.case = TRUE)][1],
    "2016" = nlcd_paths[grepl("nlcd.*2016", nlcd_paths, ignore.case = TRUE)][1],
    "2020" = nlcd_paths[grepl("nlcd.*2020", nlcd_paths, ignore.case = TRUE)][1]
  )
  
  tcc_files <- list(
    "2010" = nlcd_paths[grepl("tcc.*2010", nlcd_paths, ignore.case = TRUE)][1],
    "2016" = nlcd_paths[grepl("tcc.*2016", nlcd_paths, ignore.case = TRUE)][1],
    "2020" = nlcd_paths[grepl("tcc.*2020", nlcd_paths, ignore.case = TRUE)][1]
  )
  
  # Initialize workers
  workers <- max(1, future::availableCores() - 10)
  future::plan(future::multisession, workers = workers)
  
  tictoc::tic("Total Metrics Extraction Time")
  
  # Run parallel extraction
  spatial_results <- furrr::future_pmap_dfr(
    list(
      current_id = sites_with_bbox$id,
      xmin = sites_with_bbox$xmin,
      ymin = sites_with_bbox$ymin,
      xmax = sites_with_bbox$xmax,
      ymax = sites_with_bbox$ymax
    ),
    .f = safe_metrics_worker, # Assuming this function is defined above in your script
    nlcd_paths = nlcd_files,
    tcc_paths = tcc_files,
    .options = furrr::furrr_options(seed = TRUE, packages = c("terra", "sf"))
  )
  
  future::plan(future::sequential)
  tictoc::toc()
  
  # Join metrics back to the main spatial dataframe
  final_sites <- sites_with_bbox |>
    dplyr::left_join(spatial_results, by = "id")
  
  # Export the final dataset
  readr::write_csv(final_sites, metricsExport)
  message(sprintf("Metrics successfully extracted and saved to: %s", metricsExport))
  
} else {
  message(sprintf("Cached metrics found. Loading from: %s", metricsExport))
  final_sites <- readr::read_csv(metricsExport, show_col_types = FALSE)
}




TARGET_TOTAL_SAMPLE <- 200

# --- 1. DATA PREPARATION ---
# Melt the dataset so every unique site/year combination is an independent row
long_sites <- final_sites |>
  tidyr::pivot_longer(
    cols = matches("p_forest_|mean_tcc_"),
    names_to = c(".value", "year"),
    names_pattern = "(.*)_(\\d{4})"
  ) |>
  dplyr::mutate(site_year_id = paste(id, year, sep = "_")) |>
  dplyr::filter(!is.na(p_forest) & !is.na(mean_tcc))

# --- 2. ASSIGN INDEPENDENT CLASSES ---
classified_pop <- long_sites |>
  dplyr::mutate(
    # TCC Classes (Used for Neyman TCC and Equal TCC tracks)
    tcc_class = dplyr::case_when(
      mean_tcc == 0 ~ 0,
      mean_tcc > 0 & mean_tcc <= 10 ~ 1,
      mean_tcc > 10 & mean_tcc <= 25 ~ 2,
      mean_tcc > 25 & mean_tcc <= 50 ~ 3,
      mean_tcc > 50 ~ 4,
      TRUE ~ NA_real_
    ),
    # Forest Classes (Used for Neyman Forest track)
    forest_class = dplyr::case_when(
      p_forest == 0 ~ 0,
      p_forest > 0 & p_forest <= 10 ~ 1,
      p_forest > 10 & p_forest <= 25 ~ 2,
      p_forest > 25 & p_forest <= 50 ~ 3,
      p_forest > 50 ~ 4,
      TRUE ~ NA_real_
    )
  ) |>
  dplyr::filter(!is.na(tcc_class) & !is.na(forest_class))


# --- 3. ALLOCATION CALCULATIONS (WITH BASELINE) ---

#' Helper Function: Calculate Neyman Weights and Targets with a Baseline Minimum
calc_neyman_with_baseline <- function(df, class_col, var_col, target_n, baseline_per_class = 5) {
  
  # Determine total baseline budget to reserve
  n_classes <- dplyr::n_distinct(df[[class_col]])
  total_baseline <- baseline_per_class * n_classes
  remaining_budget <- target_n - total_baseline
  
  if (remaining_budget < 0) {
    stop("Baseline budget exceeds the total target sample size.")
  }
  
  metrics <- df |>
    dplyr::group_by(!!rlang::sym(class_col)) |>
    dplyr::summarize(
      N_h = dplyr::n(),
      S_h = sd(!!rlang::sym(var_col), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      S_h = ifelse(is.na(S_h) | S_h == 0, 0.001, S_h), # Prevent zero-variance collapse
      weight = N_h * S_h,
      # Distribute the remaining budget using Neyman weights
      neyman_portion = round(remaining_budget * (weight / sum(weight))),
      # The final target is the baseline plus the Neyman portion
      target = baseline_per_class + neyman_portion
    )
  
  # Correct any rounding disparities to force exactly target_n
  diff <- target_n - sum(metrics$target)
  if (diff != 0) {
    # Apply the difference to the class that received the largest Neyman allocation
    max_idx <- which.max(metrics$neyman_portion)
    metrics$target[max_idx] <- metrics$target[max_idx] + diff
  }
  
  return(metrics)
}

# Track A: Neyman Allocation on NLCD Forest (with baseline 5)
alloc_neyman_forest <- calc_neyman_with_baseline(classified_pop, "forest_class", "p_forest", TARGET_TOTAL_SAMPLE, baseline_per_class = 5)

# Track B: Neyman Allocation on Tree Canopy Cover (with baseline 5)
alloc_neyman_tcc <- calc_neyman_with_baseline(classified_pop, "tcc_class", "mean_tcc", TARGET_TOTAL_SAMPLE, baseline_per_class = 5)

# Track C: Equal Allocation on Tree Canopy Cover (remains unchanged)
alloc_equal_tcc <- data.frame(
  tcc_class = sort(unique(classified_pop$tcc_class)),
  target = TARGET_TOTAL_SAMPLE / 5
)

# Print the NLCD Forest targets to verify the 0 class received its baseline
print("NLCD Forest Targets:")
print(alloc_neyman_forest |> dplyr::select(forest_class, neyman_portion, target))

# --- 4. MLRA-BALANCED EXTRACTION ---
library(dplyr)
library(tidyr)
library(rlang)

#' Dynamic MLRA-Balanced Card Dealer Sampling
#'
#' @param available_pool The full dataframe of available sites (must contain MLRA_ID, id, year, and the class_col).
#' @param allocation_df The dataframe containing the target allocations (must contain the class_col and a 'target' column).
#' @param class_col A string specifying the column containing the stratification class (e.g., "forest_class").
#' @return A dataframe of the final selected sites.
draw_dynamic_sample <- function(available_pool, allocation_df, class_col) {
  
  # 1. Prepare the Target Deck (Highest class first)
  # rlang::sym() allows us to pass column names as strings into dplyr functions
  groupOrder <- allocation_df |>
    dplyr::arrange(dplyr::desc(!!rlang::sym(class_col))) |>
    tidyr::uncount(weights = target) |>
    dplyr::pull(!!rlang::sym(class_col))
  
  # 2. Initialize Tracking Variables
  selected_sites <- list()
  
  # Create a named vector to track how many sites each MLRA has received
  mlra_names <- unique(available_pool$MLRA_ID)
  mlra_counts <- setNames(rep(0, length(mlra_names)), mlra_names)
  
  # 3. Dynamic Card Dealer Loop
  for (current_class in groupOrder) {
    
    # Filter available pool for the required class
    candidates <- available_pool |> 
      dplyr::filter(!!rlang::sym(class_col) == current_class)
    
    if (nrow(candidates) == 0) {
      warning(paste("Global shortage: No grids left for", class_col, "class", current_class))
      next 
    }
    
    # Identify which MLRAs have this class available
    eligible_mlras <- unique(candidates$MLRA_ID)
    
    # Find the eligible MLRA with the lowest current sample count
    current_eligible_counts <- mlra_counts[as.character(eligible_mlras)]
    
    # Break ties randomly
    min_count <- min(current_eligible_counts)
    tied_mlras <- names(current_eligible_counts[current_eligible_counts == min_count])
    target_mlra <- if (length(tied_mlras) > 1) sample(tied_mlras, 1) else tied_mlras
    
    # Randomly select exactly ONE site from the chosen MLRA
    draw <- candidates |> 
      dplyr::filter(MLRA_ID == target_mlra) |> 
      dplyr::slice_sample(n = 1)
    
    # Store the drawn site
    selected_sites[[length(selected_sites) + 1]] <- draw
    
    # Remove the selected site from the available pool using exact 'id' and 'year' match
    available_pool <- available_pool |> 
      dplyr::filter(!(id == draw$id & year == draw$year))
    
    # Increment the sample count for the winning MLRA
    mlra_counts[as.character(target_mlra)] <- mlra_counts[as.character(target_mlra)] + 1
  }
  
  # 4. Finalize Dataset
  final_sample <- dplyr::bind_rows(selected_sites)
  return(final_sample)
}
# --- EXECUTE THE THREE SAMPLING TRACKS ---

# Set a seed globally so the draws are reproducible, but do it outside the function 
# so each track gets a unique sequence of random numbers.
set.seed(2026)

# Track A: NLCD Forest Neyman Allocation
message("\nProcessing Neyman Forest Allocation...")
sample_neyman_forest <- draw_dynamic_sample(
  available_pool = classified_pop, 
  allocation_df = alloc_neyman_forest, 
  class_col = "forest_class"
)

# Track B: Tree Canopy Cover Neyman Allocation
message("Processing Neyman TCC Allocation...")
sample_neyman_tcc <- draw_dynamic_sample(
  available_pool = classified_pop, 
  allocation_df = alloc_neyman_tcc, 
  class_col = "tcc_class"
)

# Track C: Tree Canopy Cover Equal Allocation
message("Processing Equal TCC Allocation...")
sample_equal_tcc <- draw_dynamic_sample(
  available_pool = classified_pop, 
  allocation_df = alloc_equal_tcc, 
  class_col = "tcc_class"
)

# --- VERIFICATION ---
message("\n--- Final Output Summaries ---")
message(sprintf("Neyman Forest Total Sites: %d", nrow(sample_neyman_forest)))
message(sprintf("Neyman TCC Total Sites:    %d", nrow(sample_neyman_tcc)))
message(sprintf("Equal TCC Total Sites:     %d", nrow(sample_equal_tcc)))

# Check the MLRA distributions for one of the outputs
print("MLRA Distribution for Neyman Forest:")
print(table(sample_neyman_forest$MLRA_ID))



