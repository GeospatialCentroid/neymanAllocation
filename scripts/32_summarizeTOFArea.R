pacman::p_load(dplyr, terra, sf, readr, furrr, future, readr, ggplot2)

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
}
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

  
  
  # visualize the results 
  # Load necessary libraries
  library(readr)
  library(dplyr)
  library(ggplot2)
  
  # Function 1: Compare absolute area_weighted_tof across scenarios
  plot_absolute_tof <- function(file_path) {
    # Read the CSV data
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Generate a horizontal bar chart ordered by TOF value
    p <- ggplot(df, aes(x = reorder(scenario, area_weighted_tof), y = area_weighted_tof)) +
      geom_col(fill = "steelblue") +
      coord_flip() + 
      labs(
        title = "Area Weighted TOF by Scenario",
        x = "Scenario",
        y = "Area Weighted TOF"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8) # Adjust text size for readability
      )
    
    return(p)
  }
  
  # Function 2: Visualize the difference relative to a baseline scenario
  plot_tof_difference <- function(file_path, baseline_scenario_name) {
    # Read the CSV data
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Extract the baseline value
    baseline_value <- df %>% 
      filter(scenario == baseline_scenario_name) %>% 
      pull(area_weighted_tof)
    
    if (length(baseline_value) == 0) {
      stop("Baseline scenario not found in the dataset.")
    }
    
    # Calculate the difference and determine if it is positive or negative
    df_diff <- df %>%
      mutate(
        tof_diff = area_weighted_tof - baseline_value,
        diff_type = ifelse(tof_diff > 0, "Increase", "Decrease")
      ) %>%
      filter(scenario != baseline_scenario_name) # Optional: remove baseline from plot
    
    # Generate a diverging bar chart
    p <- ggplot(df_diff, aes(x = reorder(scenario, tof_diff), y = tof_diff, fill = diff_type)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("Increase" = "darkgreen", "Decrease" = "darkred")) +
      labs(
        title = paste("Difference in Area Weighted TOF vs", baseline_scenario_name),
        x = "Scenario",
        y = "Difference in Area Weighted TOF",
        fill = "Change"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 8)
      )
    
    return(p)
  }
  

  
  # Function 1: Line chart showing TOF trends over time for each scenario
  plot_tof_trend <- function(file_path) {
    # Read the CSV data
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Generate a line chart with points
    p <- ggplot(df, aes(x = year, y = area_weighted_tof, color = scenario, group = scenario)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      labs(
        title = "Area Weighted TOF Trends Over Time",
        x = "Year",
        y = "Area Weighted TOF",
        color = "Scenario"
      ) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    return(p)
  }
  
  # Function 2: Faceted bar chart comparing scenarios separated by year
  
  plot_tof_faceted_by_year <- function(file_path) {
    # Read the CSV data
    df <- read_csv(file_path, show_col_types = FALSE)
    
    # Ensure year is treated as a factor for clear faceting
    df$year <- as.factor(df$target_year)
    
    # Define the exact desired order from top to bottom
    scenario_order <- c(
      "base_neyman_forest",
      "forest_max_100_strat_forest",
      "forest_max_50_strat_forest",
      "forest_max_20_strat_forest",
      "forest_max_0_strat_forest",   # Included based on earlier dataset
      "base_neyman_tcc",
      "tcc_max_100_strat_tcc",
      "tcc_max_50_strat_tcc",
      "tcc_max_20_strat_tcc",
      "tcc_max_0_strat_tcc"          # Included for consistency
    )
    
    # Identify scenarios that are actually in the dataset to avoid empty factor levels
    actual_scenarios <- unique(df$scenario)
    
    # Get the ordered ones that exist, and append any unexpected ones to the end
    ordered_levels <- intersect(scenario_order, actual_scenarios)
    other_levels <- setdiff(actual_scenarios, scenario_order)
    
    # Reverse the levels because coord_flip() plots the first level at the bottom
    final_levels <- rev(c(ordered_levels, other_levels))
    
    # Apply the ordered factor to the scenario column
    df$scenario <- factor(df$scenario, levels = final_levels)
    
    # Calculate the global maximum TOF to explicitly set the uniform axis limit
    max_tof <- max(df$area_weighted_tof, na.rm = TRUE)
    
    # Generate a bar chart faceted by year
    # Note: 'reorder(scenario, area_weighted_tof)' is replaced with just 'scenario'
    p <- ggplot(df, aes(x = scenario, y = area_weighted_tof, fill = scenario)) +
      geom_col() +
      coord_flip() +
      # Lock the visual x-axis (mapped to y) to be exactly the same across all facets
      scale_y_continuous(limits = c(0, max_tof * 1.05)) + 
      facet_wrap(~ year, scales = "fixed") + 
      labs(
        title = "Area Weighted TOF by Scenario Across Years",
        x = "Scenario",
        y = "Area Weighted TOF"
      ) +
      theme_minimal() +
      theme(
        legend.position = "none", 
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold")
      )
    
    return(p)
  }
  
  
  # set paths to files 
  tof_combined <- file.path(export_dir, "scenario_tof_summary_combined.csv")
  tof_year <- file.path(export_dir, "scenario_tof_summary_target_year.csv")
  
  # 1. Plot the absolute values
  plot_1 <- plot_absolute_tof(tof_combined)
  print(plot_1)
  
  # 2. Plot the difference relative to the 'base_neyman_forest' scenario
  plot_1b <- plot_tof_difference(tof_combined, "base_neyman_forest")
  print(plot_1b)
  
  # 2. Plot the faceted bar charts broken out by year
  faceted_plot <- plot_tof_faceted_by_year(tof_year)
  print(faceted_plot)
  
  
  