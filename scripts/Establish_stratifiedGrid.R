# ==============================================================================
# Establish_stratifiedGrid.R
# Purpose: Develop the standardized stratified sample of LLR using a specified
#          number of features per MLRA.
# Output:  A comprehensive dataframe containing the 1km cell ID, MLRA ID,
#          and assigned LLR ID for all sampled grids.
# ==============================================================================

source("scripts/00_config.R")

library(tmap)
library(dplyr)
library(sf)

# --- 1. SETUP & LOCAL PATHS ---------------------------------------------------
# Define sampling parameters
DRAW_SIZE <- 1400
TARGET_LLR <- "F"

message("Loading spatial inputs...")

lrr_id_path <- "data/derived/mlra/lower48MLRA.gpkg"
mlra_grid_path <- "data/derived/grids/GreatPlains_1km_mlra.gpkg"

lrrID <- sf::st_read(lrr_id_path, quiet = TRUE)
mlras <- sf::st_read(mlra_grid_path, quiet = TRUE)
grid100 <- sf::st_read("data/raw/grid100km_aea.gpkg", quiet = TRUE)

message(paste("Subsetting MLRAs for LLR:", TARGET_LLR))

llr_target_ids <- lrrID %>%
  dplyr::filter(LRRSYM == TARGET_LLR) %>%
  dplyr::pull(MLRA_ID)

mlras_target <- mlras %>%
  dplyr::filter(MLRA_ID %in% llr_target_ids)


# --- 2. HELPER FUNCTIONS ------------------------------------------------------
#' Draw Spatially Balanced Sample per MLRA and format as Dataframe
#' Now uses dynamic 1km grid ID extraction via the getAOI workflow
get_mlra_sample_df <- function(
  spatial_grid,
  target_mlra,
  n_desired,
  llr_val,
  grid100_obj,
  seed_val = 1234
) {
  # Subset spatial data to the current MLRA (keeps sf geometry active)
  mlra_subset <- spatial_grid[spatial_grid$MLRA_ID == target_mlra, ]

  if (nrow(mlra_subset) == 0) {
    message(paste("   No grids found for MLRA", target_mlra, "- Skipping."))
    return(data.frame(
      id = character(0),
      MLRA_ID = character(0),
      LLR_ID = character(0)
    ))
  }

  message(paste("   Drawing spatial sample for MLRA", target_mlra, "..."))

  set.seed(seed_val)

  # 1. Generate spatially balanced points across the MLRA area
  sampled_points_sfc <- sf::st_sample(
    x = sf::st_union(mlra_subset),
    size = n_desired,
    type = "regular"
  )

  # 2. Convert the geometry list into a proper sf data frame
  sampled_points_sf <- sf::st_sf(geometry = sampled_points_sfc)

  message("   Extracting 1km grid IDs for sampled points...")

  # 3. Use the adapted getAOI function to extract full 1km grid IDs
  extracted_grids_df <- extract_grid_ids_for_points(
    sampled_points = sampled_points_sf,
    grid100 = grid100_obj
  )

  # 4. Clean out any NA results using dplyr to avoid the drop=TRUE vector bug
  extracted_grids_df <- dplyr::filter(extracted_grids_df, !is.na(grid_id))

  # 5. Construct and return the final dataframe
  result_df <- data.frame(
    id = extracted_grids_df$grid_id,
    MLRA_ID = target_mlra,
    LLR_ID = llr_val,
    stringsAsFactors = FALSE
  )

  # 6. Ensure no duplicate 1km grids were selected
  result_df <- unique(result_df)

  return(result_df)
}
extract_grid_ids_for_points <- function(sampled_points, grid100) {
  # 1. Ensure the points are in the equal-area projection used by getAOI
  if (
    is.na(sf::st_crs(sampled_points)) || sf::st_crs(sampled_points)$epsg != 5070
  ) {
    message("Transforming sampled points to EPSG:5070...")
    sampled_points <- sf::st_transform(sampled_points, crs = "EPSG:5070")
  }

  # Add a temporary unique identifier to track points during the loop
  sampled_points$temp_pt_id <- 1:nrow(sampled_points)

  message(paste("Processing", nrow(sampled_points), "points..."))

  # 2. Iterate over every point in the sf object
  # Using purrr::map_dfr to automatically bind the results into a dataframe
  grid_results <- purrr::map_dfr(1:nrow(sampled_points), function(i) {
    # Isolate the single point feature for this iteration
    single_point <- sampled_points[i, ]

    # Intersect with the 100km grid
    g1 <- grid100[single_point, ]

    # Handle edge cases where a point might fall outside the 100km grid
    if (nrow(g1) == 0) {
      return(data.frame(
        temp_pt_id = single_point$temp_pt_id,
        grid_id = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    # Run the nested grid generation from getAOI
    # filter and generate to new area 50k
    t1 <- buildSubGrids(grids = g1, cell_size = 50000, aoi = g1)[single_point, ]
    # filter and generate to new area 10k
    t2 <- buildSubGrids(grids = t1, cell_size = 10000, aoi = t1)[single_point, ]
    # filter and generate to new area 2k
    t3 <- buildSubGrids(grids = t2, cell_size = 2000, aoi = t2)[single_point, ]
    # generate 1km grids
    t4 <- buildSubGrids(grids = t3, cell_size = 1000, aoi = t3)[single_point, ]

    # Extract the full grid ID string from the final 1km feature
    final_id <- as.data.frame(t4)$id

    # Return a single-row dataframe for this point
    return(data.frame(
      temp_pt_id = single_point$temp_pt_id,
      grid_id = final_id,
      stringsAsFactors = FALSE
    ))
  })

  # 3. Clean up and format the final output
  final_df <- sampled_points |>
    # Drop the spatial geometry to return a pure dataframe
    sf::st_drop_geometry() |>
    # Join the newly generated IDs back to the original point attributes
    dplyr::left_join(grid_results, by = "temp_pt_id") |>
    # Remove the temporary tracking ID
    dplyr::select(-temp_pt_id)

  message("Grid ID extraction complete.")
  return(final_df)
}

# --- 3. MAIN PROCESSING LOOP --------------------------------------------------

if (length(llr_target_ids) == 0) {
  stop("No valid MLRA IDs found for the target LLR.")
}

message(paste(
  "\nFound",
  length(llr_target_ids),
  "MLRAs to process for sampling."
))

sample_list <- lapply(llr_target_ids, function(m_id) {
  get_mlra_sample_df(
    spatial_grid = mlras_target,
    target_mlra = m_id,
    n_desired = DRAW_SIZE,
    llr_val = TARGET_LLR,
    grid100_obj = grid100
  )
})


final_sample_df <- dplyr::bind_rows(sample_list)

message("\nSampling process complete.")
message(paste("Total grids sampled:", nrow(final_sample_df)))
# export
readr::write_csv(
  final_sample_df,
  "data/products/systematicSampleSelection/selectedSample_lrr_F_05_2026.csv"
)
# --- 4. VISUALIZATION QA/QC ---------------------------------------------------

test_mlra_id <- 56

if (test_mlra_id %in% llr_target_ids) {
  message(paste("\nGenerating QA map for MLRA", test_mlra_id, "..."))

  test_mlra_spatial <- mlras_target[mlras_target$MLRA_ID == test_mlra_id, ]
  test_mlra_samples <- final_sample_df[
    final_sample_df$MLRA_ID == test_mlra_id,
  ]

  test_mlra_spatial$inSample <- test_mlra_spatial$id %in% test_mlra_samples$id

  tmap::tmap_mode("view")

  qa_map <- tmap::tm_shape(test_mlra_spatial) +
    tmap::tm_polygons(
      col = "inSample",
      palette = c("FALSE" = "#cccccc", "TRUE" = "#e31a1c"),
      alpha = 0.7,
      border.col = "white",
      border.alpha = 0.3,
      title = paste("Sampled 1km Grids - MLRA", test_mlra_id)
    )

  print(qa_map)
}

# view results 
r1 <- readr::read_csv("data/products/systematicSampleSelection/selectedSample_lrr_F_05_2026.csv")

