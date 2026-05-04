# Load required libraries
library(sf)
library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(stringr)

# 1. Source your custom spatial generation functions
source("src/sampleGridsFunctions.R")

# =====================================================================
# DATA PREPARATION FUNCTION
# =====================================================================
#' Convert CSV sample IDs into spatial sf objects using generateAOI.R
#' @param sample_csv Path to your CSV
#' @param grid100_obj Your base 100km grid sf object (required by getAOI)
prepare_sample_grids <- function(sample_csv, grid100_obj) {
  # Read the CSV
  samples <- read_csv(sample_csv, show_col_types = FALSE)

  # Iterate over every ID in the CSV and generate the 1km grid geometry
  # using your getAOI function
  grid_list <- purrr::map(samples$id, function(current_id) {
    # Call your getAOI function
    # (Assuming it returns an sf polygon for the 1km grid)
    aoi_grid <- getAOI(grid100 = grid100_obj, id = current_id)

    # Assign the ID back into the spatial object for joining
    aoi_grid$id <- current_id
    return(aoi_grid)
  })

  # Bind the list of grids into a single spatial sf object
  spatial_grids <- dplyr::bind_rows(grid_list)

  # Join the MLRA_ID and LLR_ID from the original CSV onto the spatial grids
  spatial_grids <- dplyr::left_join(spatial_grids, samples, by = "id")

  return(spatial_grids)
}


# =====================================================================
# FUNCTION 1: MAP A SPECIFIC MLRA
# =====================================================================
#' Produces a map for a specific MLRA ID zooming into the 1km sample boundaries
#' @param target_mlra_id The specific MLRA ID you want to plot (e.g., 56)
#' @param mlra_boundaries The sf object containing all MLRAs
#' @param sample_grids The sf object containing the 1km spatial sample grids
plot_specific_mlra <- function(target_mlra_id, mlra_boundaries, sample_grids) {
  # Filter boundaries for the specific MLRA
  # (Note: Replace "MLRARSYM" with your exact column name in the gpkg if different)
  mlra_subset <- mlra_boundaries %>%
    filter(MLRARSYM == target_mlra_id | MLRA_ID == target_mlra_id)

  # Filter the spatial grids for the specific MLRA
  grid_subset <- sample_grids %>%
    filter(MLRA_ID == target_mlra_id)

  # Generate Map
  ggplot() +
    # Plot the specific MLRA region
    geom_sf(
      data = mlra_subset,
      fill = "lightblue",
      color = "darkblue",
      alpha = 0.4
    ) +
    # Plot the 1km sample grids inside it
    geom_sf(
      data = grid_subset,
      fill = "red",
      color = "black",
      alpha = 0.8,
      size = 0.5
    ) +
    theme_minimal() +
    labs(
      title = paste("Sampling Strategy Map - MLRA:", target_mlra_id),
      subtitle = paste(
        "Showing",
        nrow(grid_subset),
        "selected 1km sample grid(s)"
      ),
      x = "Longitude",
      y = "Latitude"
    )
}


# =====================================================================
# FUNCTION 2: MAP ALL LOCATIONS (AS POINTS) ACROSS ALL MLRAS
# =====================================================================
#' Produces a large-scale map displaying 1km grids as points over MLRA boundaries
#' @param mlra_boundaries The sf object containing all MLRAs
#' @param sample_grids The sf object containing the 1km spatial sample grids
plot_all_mlras <- function(mlra_boundaries, sample_grids) {
  # Convert the 1km grid polygons into points (centroids)
  # so they render properly on a wide-scale map
  grid_points <- st_centroid(sample_grids)

  # Generate Map
  ggplot() +
    # Plot all MLRA boundaries
    geom_sf(
      data = mlra_boundaries,
      fill = "whitesmoke",
      color = "darkgray",
      size = 0.2
    ) +
    # Plot the sample grids as points
    geom_sf(data = grid_points, color = "red", size = 1.5, alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Complete Sampling Strategy Map",
      subtitle = "All 1km Sample Locations shown as points across MLRA boundaries",
      x = "Longitude",
      y = "Latitude"
    )
}

# =====================================================================
# EXECUTION WORKFLOW
# =====================================================================
# To run this script, follow these steps:

# 1. Ensure you have yo"ur base 100km grid loaded in your environment!
my_grid100 <- st_read("data/raw/grid100km_aea.gpkg")

# 2. Load MLRA Geopackage
mlra_sf <- st_read("data/derived/mlra/lower48MLRA.gpkg")

# 3. Generate the spatial objects from the CSV using your functions
sample_spatial_grids <- prepare_sample_grids(
  sample_csv = "Ldata/products/systematicSampleSelection/selectedSample.csv",
  grid100_obj = my_grid100
)

# 4. Create the specific map (For example, MLRA ID 56)
# map_mlra_56 <- plot_specific_mlra(target_mlra_id = 56,
#                                   mlra_boundaries = mlra_sf,
#                                   sample_grids = sample_spatial_grids)
# print(map_mlra_56)

# 5. Create the comprehensive map (Grids converted to points)
# map_all <- plot_all_mlras(mlra_boundaries = mlra_sf,
#                           sample_grids = sample_spatial_grids)
# print(map_all)
