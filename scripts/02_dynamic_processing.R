# ==============================================================================
# 02_dynamic_processing.R
# Purpose: Attribute 1km grids with dynamic Change Over Time (COT) data.
#          It links 1km grids to 12-mile COT rasters, extracts pixel counts,
#          summarizes them into %TOF per year, and merges with static attributes.
# Output:  A master dataset (Long Format) for Neyman Allocation.
# ==============================================================================

# Source configuration from the scripts folder
source("scripts/00_config.R")

# --- 1. SETUP & PATHS ---------------------------------------------------------

# Define where the Step 1 (Static) outputs are located
STATIC_DIR <- file.path(DERIVED_DIR, "static_attributes")

# Define where the Master datasets will be saved
OUTPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")

if (!dir.exists(OUTPUT_DYNAMIC_DIR)) {
  dir.create(OUTPUT_DYNAMIC_DIR, recursive = TRUE)
}

# Define the directory containing COT raster tiles
# Note: This expects the COT data to be in data/raw/cot_10meter
COT_DIR <- file.path(INPUT_DIR, "cot_10meter")
# UNET COT
COT_DIR_UNET <- file.path(INPUT_DIR, "cot_unet")


# --- 2. LOAD SHARED INPUTS ----------------------------------------------------
message("Loading shared spatial inputs...")

# Load Vectors using paths defined in 00_config.R
mlra_poly <- terra::vect(STATIC_INPUTS$mlra)
grid_1km <- terra::vect(STATIC_INPUTS$grid_1km)
grid_12mile <- terra::vect(STATIC_INPUTS$grid_12mile)


# --- 3. HELPER FUNCTIONS ------------------------------------------------------

#' Identify All Intersecting Parent Grids
#'
#' Finds every 12-mile grid that touches a 1km grid.
#' Returns a "Long" dataframe (One 1km grid may appear multiple times).
#'
#' @param child_vect Terra vector of 1km grids (The "Child")
#' @param parent_vect Terra vector of 12-mile grids (The "Parent")
#' @return Dataframe with child_id and parent_id columns
get_all_intersections <- function(
  child_vect,
  parent_vect,
  child_id_col = "id",
  parent_id_col = "uid"
) {
  message("   ...mapping intersections (1km <-> 12mile)...")

  # 1. Coordinate System Check
  if (terra::crs(child_vect) != terra::crs(parent_vect)) {
    warning("CRS mismatch! Reprojecting parent grid to match child...")
    parent_vect <- terra::project(parent_vect, terra::crs(child_vect))
  }

  # 2. Fast Relation Check
  relation_matrix <- terra::relate(
    child_vect,
    parent_vect,
    relation = "intersects",
    pairs = TRUE
  )

  # 3. Map Indices to Actual IDs
  # Check if parent column exists, if not assume it's the first column or specific logic
  if (!parent_id_col %in% names(parent_vect)) {
    stop(paste0("Column '", parent_id_col, "' not found in parent vector."))
  }

  df_links <- data.frame(
    child_id = child_vect[[child_id_col]][relation_matrix[, 1], ],
    parent_id = parent_vect[[parent_id_col]][relation_matrix[, 2], ]
  )

  names(df_links) <- c(child_id_col, parent_id_col)
  return(df_links)
}

#' Count COT Pixels (0-9)
#'
#' Custom summary function for exactextractr.
#' Counts coverage of specific pixel values (0-9).
#' @param df Dataframe from exact_extract with 'value' and 'coverage_fraction'
count_cot_pixels <- function(df) {
  # 1. Filter valid range 0-9
  valid_mask <- !is.na(df$value) & df$value >= 0 & df$value <= 9
  vals <- df$value[valid_mask]
  covs <- df$coverage_fraction[valid_mask]

  # 2. Sum coverage (Effective Pixel Count)
  counts <- tapply(covs, vals, sum)

  # 3. Format Output: Ensure we return 0 for missing classes
  out_vec <- setNames(numeric(10), paste0("count_", 0:9)) # count_0 ... count_9

  if (length(counts) > 0) {
    # Match names found (e.g., "1") to output names ("count_1")
    out_vec[paste0("count_", names(counts))] <- counts
  }
  return(out_vec)
}


# --- 4. PRE-PROCESSING LOOP: SPATIAL RELATIONS --------------------------------

available_ids <- unique(mlra_poly$MLRA_ID)

if (!is.null(TARGET_MLRA_IDS)) {
  mlra_ids <- intersect(available_ids, TARGET_MLRA_IDS)
  missing_ids <- setdiff(TARGET_MLRA_IDS, available_ids)
  if (length(missing_ids) > 0) {
    warning(paste(
      "Requested IDs missing:",
      paste(missing_ids, collapse = ", ")
    ))
  }
} else {
  mlra_ids <- available_ids
}

if (length(mlra_ids) == 0) {
  stop("No valid MLRA IDs found.")
}

message(paste(
  "Processing",
  length(mlra_ids),
  "MLRAs:",
  paste(mlra_ids, collapse = ", ")
))

# Loop A: Generate Spatial Relations
for (m_id in mlra_ids) {
  out_file <- file.path(
    OUTPUT_DYNAMIC_DIR,
    paste0("MLRA_", m_id, "_12m_spatial_relations.csv")
  )

  if (file.exists(out_file)) {
    next
  }

  mlra_grids <- grid_1km[grid_1km$MLRA_ID == m_id, ]

  vals <- get_all_intersections(
    child_vect = mlra_grids,
    parent_vect = grid_12mile,
    child_id_col = "id",
    parent_id_col = "Unique_ID"
  )

  readr::write_csv(vals, out_file)
  message(paste0("   Created Spatial Relations: ", out_file))
}


# --- 5. MAIN PROCESSING LOOP: EXTRACTION & MERGE ------------------------------

for (m_id in mlra_ids) {
  message(paste0("\nGenerating Master Dataset for MLRA: ", m_id))

  # --- A. Setup Paths ---
  spatial_ref_file <- file.path(
    OUTPUT_DYNAMIC_DIR,
    paste0("MLRA_", m_id, "_12m_spatial_relations.csv")
  )
  static_attr_file <- file.path(
    STATIC_DIR,
    paste0("MLRA_", m_id, "_static_attributes.csv")
  )
  final_out_file <- file.path(
    OUTPUT_DYNAMIC_DIR,
    paste0("MLRA_", m_id, "_master_dataset.csv")
  )

  if (file.exists(final_out_file)) {
    message("   Master dataset exists. Skipping.")
    next
  }

  # --- B. Build VRT for the Whole MLRA ---
  ref_df <- readr::read_csv(spatial_ref_file, show_col_types = FALSE)
  needed_ids <- unique(ref_df$Unique_ID)

  # Construct File Paths (Pattern: COT_10m_[ID].tif)
  # Uses the COT_DIR relative path defined in Setup
  expected_files <- file.path(COT_DIR, paste0("COT_10m_", needed_ids, ".tif"))
  valid_files <- expected_files[file.exists(expected_files)]

  if (length(valid_files) == 0) {
    warning(paste0(
      "   No matching COT files found for MLRA ",
      m_id,
      " in ",
      COT_DIR
    ))
    next
  }

  message(paste0("   Building VRT from ", length(valid_files), " tiles..."))
  r_vrt <- terra::vrt(valid_files)

  # --- C. Extract Data (Batch Process) ---
  mlra_grids_sf <- sf::st_as_sf(grid_1km[grid_1km$MLRA_ID == m_id, ])

  message("   Extracting pixel counts (0-9)...")

  cot_results_matrix <- exactextractr::exact_extract(
    r_vrt,
    mlra_grids_sf,
    fun = count_cot_pixels,
    summarize_df = TRUE,
    progress = TRUE,
    max_cells_in_memory = 3e8
  )

  # Transpose and format
  cot_matrix_t <- t(cot_results_matrix)
  cot_df <- as.data.frame(cot_matrix_t) %>%
    dplyr::mutate(id = mlra_grids_sf$id) %>%
    dplyr::select(id, everything())

  # --- D. Summarize TOF Percentages ---
  message("   Summarizing TOF percentages for 2010, 2016, 2020...")

  cot_summary <- cot_df %>%
    dplyr::mutate(
      total_pixels = rowSums(dplyr::pick(starts_with("count_")), na.rm = TRUE),
      sum_2010 = count_1 + count_4 + count_6 + count_9,
      sum_2016 = count_3 + count_4 + count_8 + count_9,
      sum_2020 = count_5 + count_6 + count_8 + count_9
    ) %>%
    dplyr::mutate(
      tof_pct_2010 = dplyr::if_else(
        total_pixels > 0,
        (sum_2010 / total_pixels) * 100,
        0
      ),
      tof_pct_2016 = dplyr::if_else(
        total_pixels > 0,
        (sum_2016 / total_pixels) * 100,
        0
      ),
      tof_pct_2020 = dplyr::if_else(
        total_pixels > 0,
        (sum_2020 / total_pixels) * 100,
        0
      )
    ) %>%
    dplyr::select(id, tof_pct_2010, tof_pct_2016, tof_pct_2020)

  # --- E. Merge, Clean, and Pivot to Long Format ---
  message("   Merging and formatting Master Dataset...")

  # Ensure static attributes exist before attempting read
  if (!file.exists(static_attr_file)) {
    warning(paste("Static attributes file not found:", static_attr_file))
    next
  }

  static_df <- readr::read_csv(static_attr_file, show_col_types = FALSE)

  master_df_long <- static_df %>%
    dplyr::left_join(cot_summary, by = "id") %>%
    dplyr::select(
      id,
      MLRA_ID,
      grid_area,
      riparian_pct = riparian_class_1,
      matches(
        "nlcd_\\d{4}_(Water|Developed|Barren|Forest|Shrubland|Herbaceous|Cultivated|Wetlands)"
      ),
      starts_with("tof_pct")
    ) %>%
    dplyr::rename_with(
      .fn = ~ gsub("tof_pct_(\\d{4})", "nlcd_\\1_TOF", .x),
      .cols = starts_with("tof_pct")
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("nlcd_"),
      names_pattern = "nlcd_(\\d{4})_(.*)",
      names_to = c("year", ".value")
    ) %>%
    dplyr::mutate(year = as.numeric(year))

  readr::write_csv(master_df_long, final_out_file)
  message(paste0("   Saved Long-Format Master: ", final_out_file))
}

message("\nDynamic processing complete.")
