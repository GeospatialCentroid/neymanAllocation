# copy of the static processing script to quickly summarize NLCD for areas outside of nebraska 
# not intended to be a standard part of the workflow. 


# ==============================================================================
# 01_static_processing.R
# Purpose: Attribute 1km vector grids with NLCD and Riparian raster data.
#          It handles calculating raw class percentages and aggregating them
#          into primary land cover groups (Forest, Water, etc.).
# Output:  A combined CSV of static attributes for each 1km grid cell per MLRA.
# ==============================================================================

# Update: Source from the correct 'scripts' directory relative to project root
source("scripts/00_config.R")

# --- 1. SETUP & LOCAL PATHS ---------------------------------------------------
# Define local output directory for static attributes
# Note: DERIVED_DIR is defined in 00_config.R as "data/derived"
OUTPUT_STATIC_DIR <- file.path(DERIVED_DIR, "mlra_nlcd_summaryArea")

# Ensure this specific output folder exists (it was not created by the migration script)
if (!dir.exists(OUTPUT_STATIC_DIR)) {
    dir.create(OUTPUT_STATIC_DIR, recursive = TRUE)
    message(paste("Created output directory:", OUTPUT_STATIC_DIR))
}

# --- 2. LOAD INPUTS -----------------------------------------------------------
message("Loading spatial inputs...")

# Load Vector Data (Paths defined in 00_config.R)
mlra_polyAll <- terra::vect("data/derived/mlra/lower48MLRA.gpkg")
# select the mlra of interest 
mlra_poly <- mlra_polyAll[mlra_polyAll$LRRSYM %in% c("F", "G","H"), ]

# Load Raster Data
r_riparian <- terra::rast(STATIC_INPUTS$riparian)

# Load NLCD rasters into a named list for iteration
r_nlcd_list <- list(
    "2010" = terra::rast(STATIC_INPUTS$nlcd$y2010),
    "2016" = terra::rast(STATIC_INPUTS$nlcd$y2016),
    "2020" = terra::rast(STATIC_INPUTS$nlcd$y2020)
)

# --- 3. HELPER FUNCTIONS ------------------------------------------------------

#' Aggregate NLCD Columns into Primary Groups
aggregate_nlcd_classes <- function(df, prefix) {
    # Helper to sum columns safely (returns 0 if no columns match pattern)
    sum_cols <- function(pattern) {
        cols <- df %>%
            dplyr::select(matches(paste0("^", prefix, "_class_", pattern, "$")))

        if (ncol(cols) == 0) {
            return(0)
        }
        rowSums(cols, na.rm = TRUE)
    }

    # create new aggregated columns
    df %>%
        dplyr::mutate(
            !!paste0(prefix, "_Water") := sum_cols("11|12"),
            !!paste0(prefix, "_Developed") := sum_cols("21|22|23|24"),
            !!paste0(prefix, "_Barren") := sum_cols("31"),
            !!paste0(prefix, "_Forest") := sum_cols("41|42|43"),
            !!paste0(prefix, "_Shrubland") := sum_cols("51|52"),
            !!paste0(prefix, "_Herbaceous") := sum_cols("71|72|73|74"),
            !!paste0(prefix, "_Cultivated") := sum_cols("81|82"),
            !!paste0(prefix, "_Wetlands") := sum_cols("90|95")
        )
}

#' Calculate Class Percentages for Vector Polygons
calculate_raster_stats <- function(vect_layer, raster_layer, prefix) {
    # Ensure CRS alignment
    if (terra::crs(raster_layer) != terra::crs(vect_layer)) {
        vect_layer <- terra::project(vect_layer, terra::crs(raster_layer))
    }

    message(paste0("   ...extracting ", prefix))
    # get the total area of the mlra 
    
    # Extract values (returns ID=row_index and Value)
    extract_df <- terra::extract(raster_layer, vect_layer, fun = NULL)

    val_col <- names(extract_df)[2]

    # Summarize to percentages (Wide Format)
    stats_df <- extract_df %>%
        dplyr::group_by(ID, !!rlang::sym(val_col)) %>%
        dplyr::summarise(count = n(), .groups = "drop_last") %>%
        dplyr::mutate(total = sum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(percent = (count / total) * 100) %>%
        dplyr::select(ID, class = !!rlang::sym(val_col), percent) %>%
        tidyr::pivot_wider(
            names_from = class,
            values_from = percent,
            values_fill = 0,
            names_prefix = paste0(prefix, "_class_")
        )

    # Attach back to original data for context reference 
    result <- vect_layer %>%
      as.data.frame() |>
      dplyr::bind_cols(stats_df) %>%
      dplyr::select(-ID)

    return(result)
}


# --- 4. MAIN PROCESSING LOOP --------------------------------------------------

available_ids <- unique(mlra_poly$MLRA_ID)

for (m_id in available_ids) {
    message(paste0("\nProcessing MLRA: ", m_id))

    out_file <- file.path(
        OUTPUT_STATIC_DIR,
        paste0("MLRA_", m_id, "_static_attributes.csv")
    )

    if (file.exists(out_file)) {
        message("   Output exists. Skipping.")
        next
    }

    mlra_subset <- mlra_poly[mlra_poly$MLRA_ID == m_id, ]
    if (nrow(mlra_subset) == 0) {
        message("   No grids found in this MLRA. Skipping.")
        next
    }

    message(paste0("   Found ", nrow(mlra_subset), " grid cells."))

    # B. Riparian Extraction - current file is limited to nebraska 
    # riparian_stats <- calculate_raster_stats(
    #     vect_layer = mlra_subset,
    #     raster_layer = r_riparian,
    #     prefix = "riparian"
    # )

    # C. NLCD Extraction
    nlcd_results <- list()
    for (year in names(r_nlcd_list)) {
        message(paste0("   Processing NLCD ", year))
        raw_stats <- calculate_raster_stats(
            vect_layer = mlra_subset,
            raster_layer = r_nlcd_list[[year]],
            prefix = paste0("nlcd_", year)
        )
        agg_stats <- aggregate_nlcd_classes(raw_stats, paste0("nlcd_", year))
        nlcd_results[[year]] <- agg_stats
    }
    # get total area 
    mlra_subset$grid_area <- terra::expanse(mlra_subset,unit = "km")
      
    # D. Merge All Attributes
    combined_df <- as.data.frame(mlra_subset) %>%
        dplyr::select(MLRA_ID, grid_area)
    # 
    for (year in names(nlcd_results)) {
        sel <- nlcd_results[[year]] |>
          dplyr::select(-c("MLRARSYM","MLRA_NAME","LRRSYM","LRR_NAME"))
        #bind the data 
        combined_df <- dplyr::left_join(
            combined_df,
            sel,
            by = "MLRA_ID"
        )
    }
    # E. Save Output
    readr::write_csv(combined_df, out_file)
    message(paste0("   Saved: ", out_file))
}

message("\nStatic processing complete.")


# read in the files 
# summarize so that we can assign a neyman allocation class using the following logic 
# if forest >= 1% use forest 
# if forest is < 1% and wetlands >= 0.5% use wetlands
# if forest is < 1% and wetlands is < 0.5% use class with the highest land cover percentage

# R
dir <- "data/derived/mlra_nlcd_summaryArea"
summarize_neyman_from_static <- function(dir) {
  files <- list.files(dir, pattern = "MLRA_.*_static_attributes\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    return(tibble::tibble(
      MLRA_ID = integer(),
      forest_percent = numeric(),
      wetlands_percent = numeric(),
      most_abundant_class = character(),
      most_abundant_percent = numeric(),
      neyman_allocation = character()
    ))
  }

  df <- purrr::map_dfr(files, readr::read_csv, show_col_types = FALSE)

  # Detect available NLCD years and primary group columns
  pattern <- "nlcd_(\\d{4})_([A-Za-z]+)$"
  matches <- stringr::str_match(names(df), pattern)
  matches <- as.data.frame(matches, stringsAsFactors = FALSE)
  col_info <- tibble::tibble(
    col = names(df),
    year = matches[,2],
    group = matches[,3]
  ) %>%
    filter(!is.na(year) & !is.na(group))

  # If no nlcd columns found, return empty summary
  if (nrow(col_info) == 0) {
    stop("No nlcd_*_GROUP columns found in files.")
  }

  # Choose most recent year available
  chosen_year <- max(as.integer(unique(col_info$year)))
  groups <- c("Water","Developed","Barren","Forest","Shrubland","Herbaceous","Cultivated","Wetlands")
  chosen_cols <- paste0("nlcd_", chosen_year, "_", groups)

  # Ensure grid_area exists
  if (!"grid_area" %in% names(df)) {
    stop("grid_area column not found in data (required for area-weighting).")
  }

  # Replace missing chosen columns with 0 if not present
  for (col in chosen_cols) {
    if (!col %in% names(df)) {
      df[[col]] <- 0
    }
  }

  # Area-weighted MLRA aggregation (weighted mean of per-grid percent)
  agg <- df %>%
    mutate(MLRA_ID = as.integer(MLRA_ID)) %>%
    group_by(MLRA_ID) %>%
    summarize(
      across(all_of(chosen_cols),
             ~ if (sum(grid_area, na.rm = TRUE) == 0) 0 else sum(.x * grid_area, na.rm = TRUE) / sum(grid_area, na.rm = TRUE)),
      .groups = "drop"
    )

  # Determine most abundant class and apply Neyman logic
  result <- agg %>%
    rowwise() %>%
    mutate(
      forest_percent = .data[[paste0("nlcd_", chosen_year, "_Forest")]],
      wetlands_percent = .data[[paste0("nlcd_", chosen_year, "_Wetlands")]],
      # find the most abundant primary class and its percent
      most_abundant = {
        vals <- c_across(all_of(chosen_cols))
        names(vals) <- groups
        idx <- which.max(vals)
        tibble::tibble(class = names(vals)[idx], pct = vals[idx])
      } %>% list(),
      most_abundant_class = most_abundant$class,
      most_abundant_percent = most_abundant$pct,
      neyman_allocation = dplyr::case_when(
        forest_percent >= 1 ~ "Forest",
        forest_percent < 1 & wetlands_percent >= 0.5 ~ "Wetlands",
        TRUE ~ most_abundant_class
      )
    ) %>%
    ungroup() %>%
    select(
      MLRA_ID,
      forest_percent,
      wetlands_percent,
      most_abundant_class,
      most_abundant_percent,
      neyman_allocation
    )

  # Ensure numeric columns are numeric
  result <- result %>%
    mutate(
      forest_percent = as.numeric(forest_percent),
      wetlands_percent = as.numeric(wetlands_percent),
      most_abundant_percent = as.numeric(most_abundant_percent)
    )

  return(result)
}


dir <- "data/derived/mlra_nlcd_summaryArea"
neymanSummary <- summarize_neyman_from_static(dir)
neymanSummary
# export 
write_csv(neymanSummary, "data/derived/mlra_nlcd_summaryArea/neyman_allocation_summary.csv")


# using the mlra_poly area
mp <- mlra_poly
mp$grid_area <- terra::expanse(mp,unit = "km")
mp$train_current <- mp$grid_area / 500
mp$test_density <- mp$grid_area / 1500

nTest <- sum(mp$test_density, na.rm = TRUE)
nTrain <- sum(mp$train_current, na.rm = TRUE)
totalOrignal <- round(nTest + nTrain)

fProp <- 200 / totalOrignal
mp$test_200 <- mp$test_density * fProp
mp$train_200 <- mp$train_current * fProp
#average test and train per mlra
sum(mp$test_200)/nrow(mp)
sum(mp$train_200)/nrow(mp)


# process for all MLRAs, in the great plains 
mlra_polyGP <- mlra_polyAll[
  mlra_polyAll$LRR_NAME %in% c(
    "Northern Great Plains Spring Wheat Region",
    "Western Great Plains Range and Irrigated Region",
    "Central Great Plains Winter Wheat and Range Region"
  ),
]

# remove area already used in nebraska 
new_mlra <- mlra_polyGP[!mlra_polyGP$MLRA_ID %in% ALL_MLRA_IDS, ]

# calculate area,
new_mlra$grid_area <- terra::expanse(new_mlra,unit = "km")
new_mlra$train_density <- new_mlra$grid_area / 500
new_mlra$test_density <- new_mlra$grid_area / 1500

nTest_new <- sum(new_mlra$test_density, na.rm = TRUE)
nTrain_new <- sum(new_mlra$train_density, na.rm = TRUE)
totalNew <- round(nTest_new + nTrain_new)

# estimate the proportional elemnt 
propTest <- 200 / totalNew
new_mlra$test_200 <- new_mlra$test_density * propTest
new_mlra$train_200 <- new_mlra$train_density * propTest

# average test and train per mlra 
sum(new_mlra$test_200)/nrow(new_mlra)
sum(new_mlra$train_200)/nrow(new_mlra)
