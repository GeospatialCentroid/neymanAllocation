# ==============================================================================
# Establish_stratifiedGrid.R
# Purpose: Develop the standardized stratified sample of LLR using a specified
#          number of features per MLRA. Supports multiple independent draw sizes.
# Output:  Comprehensive dataframes containing the 1km cell ID, MLRA ID,
#          and assigned LLR ID for all sampled grids, split by draw size.
# ==============================================================================

source("scripts/00_config.R")
pacman::p_load("tmap","dplyr","sf","purrr","readr")

# --- 1. SETUP & LOCAL PATHS ---------------------------------------------------
# Define sampling parameters
DRAW_SIZES <- c(500, 650, 1400) # Pass multiple draw sizes here
TARGET_LLR <- "G"

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
get_mlra_sample_df <- function(
    spatial_grid,
    target_mlra,
    n_desired,
    llr_val,
    grid100_obj,
    seed_val = 1234
) {
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
  
  sampled_points_sfc <- sf::st_sample(
    x = sf::st_union(mlra_subset),
    size = n_desired,
    type = "regular"
  )
  
  sampled_points_sf <- sf::st_sf(geometry = sampled_points_sfc)
  
  message("   Extracting 1km grid IDs for sampled points...")
  
  extracted_grids_df <- extract_grid_ids_for_points(
    sampled_points = sampled_points_sf,
    grid100 = grid100_obj
  )
  
  extracted_grids_df <- dplyr::filter(extracted_grids_df, !is.na(grid_id))
  
  result_df <- data.frame(
    id = extracted_grids_df$grid_id,
    MLRA_ID = target_mlra,
    LLR_ID = llr_val,
    stringsAsFactors = FALSE
  )
  
  result_df <- unique(result_df)
  
  return(result_df)
}

extract_grid_ids_for_points <- function(sampled_points, grid100) {
  if (
    is.na(sf::st_crs(sampled_points)) || sf::st_crs(sampled_points)$epsg != 5070
  ) {
    message("Transforming sampled points to EPSG:5070...")
    sampled_points <- sf::st_transform(sampled_points, crs = "EPSG:5070")
  }
  
  sampled_points$temp_pt_id <- 1:nrow(sampled_points)
  
  message(paste("Processing", nrow(sampled_points), "points..."))
  
  grid_results <- purrr::map_dfr(1:nrow(sampled_points), function(i) {
    single_point <- sampled_points[i, ]
    g1 <- grid100[single_point, ]
    
    if (nrow(g1) == 0) {
      return(data.frame(
        temp_pt_id = single_point$temp_pt_id,
        grid_id = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    
    t1 <- buildSubGrids(grids = g1, cell_size = 50000, aoi = g1)[single_point, ]
    t2 <- buildSubGrids(grids = t1, cell_size = 10000, aoi = t1)[single_point, ]
    t3 <- buildSubGrids(grids = t2, cell_size = 2000, aoi = t2)[single_point, ]
    t4 <- buildSubGrids(grids = t3, cell_size = 1000, aoi = t3)[single_point, ]
    
    final_id <- as.data.frame(t4)$id
    
    return(data.frame(
      temp_pt_id = single_point$temp_pt_id,
      grid_id = final_id,
      stringsAsFactors = FALSE
    ))
  })
  
  final_df <- sampled_points |>
    sf::st_drop_geometry() |>
    dplyr::left_join(grid_results, by = "temp_pt_id") |>
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

# Environment tracking to hold the last processed dataset for QA step
last_sample_df <- NULL
last_draw_size <- NULL

# Iterate through each independent draw size
for (current_draw_size in DRAW_SIZES) {
  
  message(paste("\n--- STARTING WORKFLOW FOR DRAW SIZE:", current_draw_size, "---"))
  
  sample_list <- lapply(llr_target_ids, function(m_id) {
    get_mlra_sample_df(
      spatial_grid = mlras_target,
      target_mlra = m_id,
      n_desired = current_draw_size,
      llr_val = TARGET_LLR,
      grid100_obj = grid100
    )
  })
  
  final_sample_df <- dplyr::bind_rows(sample_list)
  
  message(paste("Sampling complete for Draw Size", current_draw_size, "."))
  message(paste("Total grids sampled:", nrow(final_sample_df)))
  
  # Programmatically generate export path based on LLR and Draw Size variables
  export_filename <- sprintf(
    "selectedSample_lrr_%s_draw_%s_05_2026.csv", 
    TARGET_LLR, 
    current_draw_size
  )
  export_path <- file.path("data/products/systematicSampleSelection", export_filename)
  
  readr::write_csv(final_sample_df, export_path)
  message(paste("Exported dataset to:", export_path))
  
  # Cache references for the QA/QC block below
  last_sample_df <- final_sample_df
  last_draw_size <- current_draw_size
}

# --- 4. VISUALIZATION QA/QC ---------------------------------------------------
# Runs QA utilizing the final processed iteration from the loop above

test_mlra_id <- llr_target_ids[1]

if (test_mlra_id %in% llr_target_ids && !is.null(last_sample_df)) {
  message(sprintf("\nGenerating QA map for MLRA %s using Draw Size %s...", test_mlra_id, last_draw_size))
  
  test_mlra_spatial <- mlras_target[mlras_target$MLRA_ID == test_mlra_id, ]
  test_mlra_samples <- last_sample_df[last_sample_df$MLRA_ID == test_mlra_id, ]
  
  test_mlra_spatial$inSample <- test_mlra_spatial$id %in% test_mlra_samples$id
  
  tmap::tmap_mode("view")
  
  qa_map <- tmap::tm_shape(test_mlra_spatial) +
    tmap::tm_polygons(
      col = "inSample",
      palette = c("FALSE" = "#cccccc", "TRUE" = "#e31a1c"),
      alpha = 0.7,
      border.col = "white",
      border.alpha = 0.3,
      title = paste("Sampled 1km Grids - MLRA", test_mlra_id, "(Draw:", last_draw_size, ")")
    )
  
  print(qa_map)
}