# ==============================================================================
# 00a_generate_1km_grids.R
# Purpose: Generates standardized 1km grids from a 100km parent grid.
#          Iteratively processes by MLRA and appends to a master GeoPackage
#          to prevent RAM overflow when scaling.
# ==============================================================================

source("scripts/00_config.R")
source("src/sampleGridsFunctions.R") # Ensure this points to where you saved the functions

# --- 1. SETTINGS & INPUTS -----------------------------------------------------

# Paths for the parent structures
GRID_100KM_PATH <- file.path(INPUT_DIR, "grid100km_aea.gpkg")
LOWER48_MLRA_PATH <- file.path(DERIVED_DIR, "mlra/lower48MLRA.gpkg")

# Dynamically pull the output path from the config file's active state
OUTPUT_GRID_PATH <- STATIC_INPUTS$grid_1km

# Ensure output directory exists
if (!dir.exists(dirname(OUTPUT_GRID_PATH))) {
  dir.create(dirname(OUTPUT_GRID_PATH), recursive = TRUE)
}

# --- 2. LOAD DATA -------------------------------------------------------------

message("Loading 100km Parent Grids and MLRA Boundaries...")
parent_100km <- sf::st_read(GRID_100KM_PATH, quiet = TRUE)
mlra_poly <- sf::st_read(LOWER48_MLRA_PATH, quiet = TRUE)

# The buildGrids() function hardcodes a transformation to EPSG:5070 to ensure
# perfectly square metric grids. We must force our inputs to 5070 to match.
TARGET_CRS <- sf::st_crs(5070)

if (sf::st_crs(parent_100km) != TARGET_CRS) {
  parent_100km <- sf::st_transform(parent_100km, TARGET_CRS)
}

if (sf::st_crs(mlra_poly) != TARGET_CRS) {
  mlra_poly <- sf::st_transform(mlra_poly, TARGET_CRS)
}

# Ensure the parent grid has an 'id' column for the buildSubGrids function
if (!"id" %in% names(parent_100km)) {
  parent_100km <- parent_100km %>% dplyr::mutate(id = dplyr::row_number())
}

# Determine which MLRAs to process.
target_mlras <- if (!is.null(ALL_MLRA_IDS)) {
  ALL_MLRA_IDS
} else {
  unique(mlra_poly$MLRA_ID)
}

# --- 3. GENERATE GRIDS ITERATIVELY --------------------------------------------

message(paste("Generating 1km grids for", length(target_mlras), "MLRAs..."))

# If the output file already exists, we want to know so we don't duplicate data
file_exists_flag <- file.exists(OUTPUT_GRID_PATH)

for (m_id in target_mlras) {
  message(sprintf("   Processing MLRA: %s...", m_id))

  # 1. Isolate the specific MLRA polygon
  current_mlra <- mlra_poly %>% dplyr::filter(MLRA_ID == m_id)

  # does the file exist 
  if(!file.exists(paste0("data/derived/grids/mlra_",m_id,"_mlra.gpkg"))){
    if (nrow(current_mlra) == 0) {
      warning(paste("MLRA", m_id, "not found in polygon data. Skipping."))
      next
    }
    
    # 2. Select only the 100km parent grids that intersect this MLRA
    intersecting_100km <- select100km(parent_100km, current_mlra)
    
    if (nrow(intersecting_100km) == 0) {
      next
    }
    
    # 3. Generate 1km SubGrids
    # Uses your function: buildSubGrids(grids, cell_size, aoi)
    # cell_size is 1000 meters (1km)
    
    # sub_grids_1km <- buildSubGrids(
    #   grids = intersecting_100km,
    #   cell_size = 1000,
    #   aoi = current_mlra
    # )
    
    t1 <- buildSubGrids(grids = intersecting_100km, cell_size = 50000, aoi = intersecting_100km)
    # filter and generate to new area 10k
    t2 <- buildSubGrids(grids = t1, cell_size = 10000, aoi = t1)
    # filter and generate to new area 2k
    t3 <- buildSubGrids(grids = t2, cell_size = 2000, aoi = t2)
    # generate 1km grids
    t4 <- buildSubGrids(grids = t3, cell_size = 1000, aoi = t3)
    
    
    # 4. Attach MLRA_ID and format for export
    sub_grids_1km <- t4 %>%
      dplyr::mutate(MLRA_ID = m_id) %>%
      # crop the feature to the MLRA area 
      sf::st_intersection(current_mlra) |> 
      dplyr::select(id, MLRA_ID)
    
    
    # 5. Append to the Master GeoPackage
    # By appending, we save memory and build the dataset safely piece by piece
    sf::st_write(
      obj = sub_grids_1km,
      dsn = paste0("data/derived/grids/mlra_",m_id,"_mlra.gpkg"),
      layer = "grid_1km",
      append = file_exists_flag, # Append if true, overwrite if false
      quiet = TRUE
    )
    
    # After the first successful write, ensure all subsequent writes are appended
    file_exists_flag <- TRUE
  }
  
 
}

message("\n=== Grid Generation Complete! ===")
message("Master 1km Grid saved to: ", OUTPUT_GRID_PATH)
