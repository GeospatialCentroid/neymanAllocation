pacman::p_load(dplyr, terra, sf, readr, furrr, future)

# Load custom processing functions
source("src/processingModelData.R")
source("src/sampleSizeValidation.R")

# ------------------------------------------------------------------------------
# 2. Data Ingestion & Configuration
# ------------------------------------------------------------------------------
# Load MLRA boundaries and spatial data
mlras <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg")

# Directory with the imagery  
image_dir <- "/mnt/unraid_naip/LLR G/modelOutputs"
images <- list.files(image_dir, pattern = "\\.tif$", full.names = FALSE)

# Details on how to group the results 
summaryGroups <- read_csv("data/products/groundTruthSamples/scenarios_2020/all_scenarios_combined_200.csv", show_col_types = FALSE)

# Define export directory and file paths
export_dir <- "~/trueNAS/work/neymanSampling/data/products/groundTruthSamples/scenarios_2020"
grid_stats_path <- file.path(export_dir, "grid_level_tof.csv")
scenario_summary_path <- file.path(export_dir, "scenario_tof_summary.csv")

# ------------------------------------------------------------------------------
# 3 & 4. Data Processing and Export (Conditional Execution)
# ------------------------------------------------------------------------------

if (file.exists(grid_stats_path)){
  grid_stats <- read_csv(grid_stats_path)
}else{
  message("Outputs not found. Setting up parallel processing...")
  
  # Define the parallel plan. 
  # Leaving one core free prevents the system from locking up entirely during heavy loads.
  n_cores <- parallel::detectCores() - 16
  future::plan(future::multisession, workers = n_cores)
  
  message(sprintf("Running TOF grid calculation across %s workers...", n_cores))
  
  if(file.exists())
  
  # Generate the grid-level tibble using furrr
  grid_stats <- furrr::future_map_dfr(images, function(img_name) {
    r <- terra::rast(file.path(image_dir, img_name))
    
    # Calculate cell area (using mask = TRUE to avoid bounding box NA expansion)
    cell_areas <- terra::cellSize(r, unit = "ha", mask = TRUE)
    total_area_val <- terra::global(cell_areas, "sum", na.rm = TRUE)[[1]]
    
    # Calculate TOF area
    tof_area_raster <- r * cell_areas
    tof_area_val <- terra::global(tof_area_raster, "sum", na.rm = TRUE)[[1]]
    
    tibble(
      grid_id = tools::file_path_sans_ext(img_name),
      grid_total_area = total_area_val,
      grid_tof_area = tof_area_val
    )
  }, .options = furrr::furrr_options(seed = TRUE)) # Seed ensures reproducibility if custom functions use RNG
  
  # Shut down the parallel workers to free up memory
  future::plan(future::sequential)
  
  message("Writing grid-level TOF data...")
  write_csv(grid_stats, grid_stats_path)
  
  # ------------------------------------------------------------------------------
  # Data Parsing and Joining
  # ------------------------------------------------------------------------------
  message("Parsing grid IDs and joining to summary groups...")
  
  # 1. Rename the year column in summaryGroups to year_Neyman
  summaryGroups_clean <- summaryGroups %>%
    dplyr::rename(year_Neyman = year)
  
  # 2. Split the grid_id and map the actual year to the target year
  grid_stats_clean <- grid_stats %>%
    tidyr::separate(col = grid_id, 
                    into = c("id", "year", "suffix"), 
                    sep = "_", 
                    extra = "drop") %>%
    dplyr::select(-suffix) %>%
    # Ensure year is numeric for mapping
    dplyr::mutate(
      year = as.numeric(year),
      target_year = dplyr::case_when(
        year %in% c(2010, 2011, 2012, 2013) ~ 2012,
        year %in% c(2014, 2015, 2016, 2017) ~ 2016,
        year %in% c(2018, 2019, 2020, 2021) ~ 2020,
        TRUE ~ NA_real_ # Catch-all for any unexpected years
      )
    )
  
  # 3. Join the datasets by the 'id' key
  joined_data <- summaryGroups_clean %>%
    dplyr::left_join(grid_stats_clean, by = "id")
  
  
  # ------------------------------------------------------------------------------
  # Summarization: Per Individual Year (Actual Year)
  # ------------------------------------------------------------------------------
  message("Calculating scenario estimates per actual year...")
  
  scenario_year_estimates <- joined_data %>%
    dplyr::group_by(scenario, year) %>%
    dplyr::summarize(
      total_scenario_area = sum(grid_total_area, na.rm = TRUE),
      total_scenario_tof = sum(grid_tof_area, na.rm = TRUE),
      area_weighted_tof = total_scenario_tof / total_scenario_area,
      aoi_count = dplyr::n(),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------------------------
  # Summarization: Per Target Year Block (2012, 2016, 2020)
  # ------------------------------------------------------------------------------
  message("Calculating scenario estimates per target year block...")
  
  scenario_target_year_estimates <- joined_data %>%
    # Filter out any records that didn't map to a target year (just in case)
    dplyr::filter(!is.na(target_year)) %>% 
    dplyr::group_by(scenario, target_year) %>%
    dplyr::summarize(
      total_scenario_area = sum(grid_total_area, na.rm = TRUE),
      total_scenario_tof = sum(grid_tof_area, na.rm = TRUE),
      area_weighted_tof = total_scenario_tof / total_scenario_area,
      aoi_count = dplyr::n(),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------------------------
  # Summarization: All Years Combined
  # ------------------------------------------------------------------------------
  message("Calculating scenario estimates across all years...")
  
  scenario_combined_estimates <- joined_data %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarize(
      total_scenario_area = sum(grid_total_area, na.rm = TRUE),
      total_scenario_tof = sum(grid_tof_area, na.rm = TRUE),
      area_weighted_tof = total_scenario_tof / total_scenario_area,
      aoi_count = dplyr::n(),
      .groups = "drop"
    )
  
  # ------------------------------------------------------------------------------
  # Export
  # ------------------------------------------------------------------------------
  message("Writing scenario summary documentation...")
  write_csv(scenario_year_estimates, file.path(export_dir, "scenario_tof_summary_actual_year.csv"))
  write_csv(scenario_target_year_estimates, file.path(export_dir, "scenario_tof_summary_target_year.csv"))
  write_csv(scenario_combined_estimates, file.path(export_dir, "scenario_tof_summary_combined.csv"))
}
