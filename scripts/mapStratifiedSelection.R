# Load required libraries
library(sf)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(leaflet)

# 1. Source your custom spatial generation functions
source("src/sampleGridsFunctions.R")

# =====================================================================
# DATA PREPARATION FUNCTION
# =====================================================================
#' Convert CSV sample IDs into spatial sf objects using your custom getAOI
#' @param sample_csv Path to your CSV
#' @param grid100_obj Your base 100km grid sf object (required by getAOI)
prepare_sample_grids <- function(sample_csv, grid100_obj) {
  # Read the CSV
  samples <- read_csv(sample_csv, show_col_types = FALSE)

  # Iterate over every ID in the CSV and generate the 1km grid geometry
  grid_list <- purrr::map(samples$id, function(current_id) {
    # Call your getAOI function
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
# FUNCTION 1: MAP A SPECIFIC MLRA (LEAFLET VERSION)
# =====================================================================
#' Produces an interactive leaflet map for a specific MLRA ID
#' @param target_mlra_id The specific MLRA ID you want to plot (e.g., 56)
#' @param mlra_boundaries The sf object containing all MLRAs
#' @param sample_grids The sf object containing the 1km spatial sample grids
plot_specific_mlra <- function(target_mlra_id, mlra_boundaries, sample_grids) {
  # Filter boundaries and transform to WGS84 for Leaflet
  mlra_subset <- mlra_boundaries %>%
    filter(MLRARSYM == target_mlra_id | MLRA_ID == target_mlra_id) %>%
    st_transform(4326)

  # Filter grids and transform to WGS84 for Leaflet
  grid_subset <- sample_grids %>%
    filter(MLRA_ID == target_mlra_id) %>%
    st_transform(4326)

  # Generate Leaflet Map
  leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    # Add the MLRA boundary
    addPolygons(
      data = mlra_subset,
      fillColor = "#ADD8E6",
      color = "#00008B",
      weight = 2,
      fillOpacity = 0.3,
      popup = ~ paste("<b>MLRA ID:</b>", target_mlra_id)
    ) %>%
    # Add the 1km Sample Grids
    addPolygons(
      data = grid_subset,
      fillColor = "red",
      color = "black",
      weight = 1,
      fillOpacity = 0.8,
      popup = ~ paste("<b>Sample ID:</b>", id, "<br>", "<b>LLR ID:</b>", LLR_ID)
    ) %>%
    addControl(
      html = paste("<b>MLRA:", target_mlra_id, "Sampling Strategy</b>"),
      position = "topright"
    )
}


# =====================================================================
# FUNCTION 2: MAP ALL LOCATIONS (LEAFLET VERSION)
# =====================================================================
#' Produces a large-scale interactive map displaying 1km grids as circle markers
#' @param mlra_boundaries The sf object containing all MLRAs
#' @param sample_grids The sf object containing the 1km spatial sample grids
plot_all_mlras <- function(mlra_boundaries, sample_grids) {
  # Transform all boundaries to WGS84 for Leaflet
  mlra_wgs84 <- mlra_boundaries %>%
    st_transform(4326)

  # Convert the 1km grid polygons into points (centroids), then to WGS84
  grid_points <- sample_grids %>%
    st_centroid() %>%
    st_transform(4326)

  # Generate Leaflet Map
  leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    # Add all MLRA boundaries
    addPolygons(
      data = mlra_wgs84,
      fillColor = "#f5f5f5",
      color = "#a9a9a9",
      weight = 1,
      fillOpacity = 0.2,
      popup = ~ paste("<b>MLRA ID:</b>", MLRARSYM)
    ) %>%
    # Add the sample grids as circle markers
    addCircleMarkers(
      data = grid_points,
      color = "red",
      radius = 4,
      stroke = FALSE,
      fillOpacity = 0.7,
      popup = ~ paste(
        "<b>Sample ID:</b>",
        id,
        "<br>",
        "<b>MLRA ID:</b>",
        MLRA_ID,
        "<br>",
        "<b>LLR ID:</b>",
        LLR_ID
      )
    ) %>%
    addControl(
      html = "<b>Complete Sampling Strategy</b>",
      position = "topright"
    )
}

# =====================================================================
# EXECUTION WORKFLOW
# =====================================================================
# To run this script, follow these steps:

# 1. Load your base 100km grid from the raw data folder
my_grid100 <- st_read("data/raw/grid100km_aea.gpkg")

# 2. Load MLRA Geopackage from the derived data folder
mlra_sf <- st_read("data/derived/mlra/lower48MLRA.gpkg")

# 3. Generate the spatial objects from the CSV in the products folder
sample_spatial_grids <- prepare_sample_grids(
  sample_csv = "data/products/systematicSampleSelection/selectedSample.csv",
  grid100_obj = my_grid100
)

# 4. Create the specific map (For example, MLRA ID 56)
map_mlra_56 <- plot_specific_mlra(
  target_mlra_id = 56,
  mlra_boundaries = mlra_sf,
  sample_grids = sample_spatial_grids
)
print(map_mlra_56) # Uncomment to view in RStudio viewer

# 5. Create the comprehensive map
map_all <- plot_all_mlras(
  mlra_boundaries = mlra_sf,
  sample_grids = sample_spatial_grids
)
print(map_all) # Uncomment to view in RStudio viewer
