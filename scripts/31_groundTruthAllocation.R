# ==============================================================================
# 31_groundTruthAllocation.R
# Purpose: Develop a standardized stratified sample of a Land Resource Region (LRR)
#          using a specified number of features per Major Land Resource Area (MLRA).
#          Supports dynamic stratification thresholds and multiple allocation tracks.
# Outputs: Comprehensive dataframes containing grid IDs, MLRA IDs, and assigned
#          classes for sampled grids, exported to disk. Also generates an interactive
#          Leaflet map for spatial review.
# ==============================================================================

# --- 0. LIBRARIES & SETUP -----------------------------------------------------
source("scripts/00_config.R")
source("src/sampleGridsFunctions.R")

pacman::p_load(
  tmap, dplyr, tidyr, rlang, sf, purrr, readr, 
  future, furrr, tictoc, terra, leaflet) 

# Global Parameters
TARGET_LLR <- "G"
TARGET_SIZE <- 1400
TARGET_TOTAL_SAMPLE <- 200

# File Paths
lrr_id_path <- "data/derived/mlra/lower48MLRA.gpkg"
mlra_grid_path <- "data/derived/grids/GreatPlains_1km_mlra.gpkg"
g100 <- sf::st_read("data/raw/grid100km_aea.gpkg", quiet = TRUE)
nlcd_paths <- list.files("data/raw/nlcd/", full.names = TRUE)

siteExport <- paste0("data/derived/sampleBBOX/bbox_", TARGET_LLR, "_", TARGET_SIZE,".csv")
metricsExport <- paste0("data/derived/metrics/nlcd_tcc_metrics_", TARGET_LLR, "_", TARGET_SIZE, ".csv")

# --- 1. CORE FUNCTIONS --------------------------------------------------------

#' Calculate Neyman Weights and Targets with a Baseline Minimum
#' 
#' @description Allocates a total sample size across strata based on the variance 
#' within each stratum (Neyman allocation). Ensures every stratum receives a 
#' guaranteed baseline number of samples before variance-based allocation occurs.
#'
#' @param df A dataframe containing the population data.
#' @param class_col A string representing the column name containing the stratification classes.
#' @param var_col A string representing the column name of the continuous variable used to calculate variance.
#' @param target_n Numeric. The total number of samples to allocate across all classes.
#' @param baseline_per_class Numeric. The minimum number of samples guaranteed to each class. Default is 10.
#'
#' @return A summary dataframe containing the sample size targets (`target`) for each class.
calc_neyman_with_baseline <- function(df, class_col, var_col, target_n, baseline_per_class = 10) {
  
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
      neyman_portion = round(remaining_budget * (weight / sum(weight))),
      target = baseline_per_class + neyman_portion
    )
  
  # Correct rounding disparities to enforce exact target_n
  diff <- target_n - sum(metrics$target)
  if (diff != 0) {
    max_idx <- which.max(metrics$neyman_portion)
    metrics$target[max_idx] <- metrics$target[max_idx] + diff
  }
  
  return(metrics)
}

#' Balanced MLRA & Spatial K-Means Sampling (with Fallback)
#'
#' @description Ensures total sample counts remain balanced across MLRAs (e.g., ~10-11 sites) 
#' while utilizing spatial k-means clustering within each MLRA. Includes a fallback 
#' for overlapping grid geometries.
#' 
#' @param available_pool A dataframe of available sites containing bounding box coordinates.
#' @param allocation_df A dataframe containing the target allocations.
#' @param class_col A string specifying the column containing the stratification class.
#'
#' @return A dataframe of the final selected, spatially distributed sites.
draw_balanced_spatial_sample <- function(available_pool, allocation_df, class_col) {
  
  # 1. Calculate centroids for spatial clustering
  pool_sf <- available_pool |>
    dplyr::mutate(
      x_cent = (xmin + xmax) / 2,
      y_cent = (ymin + ymax) / 2
    )
  
  selected_site_ids <- c() 
  selected_sites <- list()
  
  mlra_names <- as.character(unique(pool_sf$MLRA_ID))
  global_mlra_counts <- setNames(rep(0, length(mlra_names)), mlra_names)
  
  allocation_df <- allocation_df |> dplyr::arrange(dplyr::desc(!!rlang::sym(class_col)))
  
  for (i in seq_len(nrow(allocation_df))) {
    current_class <- allocation_df[[class_col]][i]
    target_n <- allocation_df$target[i]
    
    candidates <- pool_sf |> 
      dplyr::filter(!!rlang::sym(class_col) == current_class) |>
      dplyr::filter(!id %in% selected_site_ids)
    
    if (nrow(candidates) == 0) {
      warning(sprintf("Global shortage: No candidate grids left for class %s", current_class))
      next
    }
    
    # 2. Card Dealer Quota Assignment
    eligible_mlras <- as.character(unique(candidates$MLRA_ID))
    class_mlra_quota <- setNames(rep(0, length(eligible_mlras)), eligible_mlras)
    
    temp_mlra_counts <- global_mlra_counts[eligible_mlras]
    candidates_avail <- candidates |> dplyr::count(MLRA_ID) |> tibble::deframe()
    
    for (draw_idx in seq_len(target_n)) {
      valid_mlras <- names(class_mlra_quota)[class_mlra_quota < candidates_avail[names(class_mlra_quota)]]
      
      if (length(valid_mlras) == 0) {
        warning(sprintf("Ran out of unique grids for class %s before hitting target.", current_class))
        break
      }
      
      current_eligible_counts <- temp_mlra_counts[valid_mlras]
      min_count <- min(current_eligible_counts)
      tied_mlras <- names(current_eligible_counts[current_eligible_counts == min_count])
      
      target_mlra <- if (length(tied_mlras) > 1) sample(tied_mlras, 1) else tied_mlras
      
      class_mlra_quota[target_mlra] <- class_mlra_quota[target_mlra] + 1
      temp_mlra_counts[target_mlra] <- temp_mlra_counts[target_mlra] + 1
    }
    
    # 3. Spatial Sub-sampling per MLRA via K-Means based on Quota
    for (t_mlra in names(class_mlra_quota)) {
      k <- class_mlra_quota[t_mlra]
      
      if (k <= 0) next
      
      mlra_candidates <- candidates |> dplyr::filter(MLRA_ID == t_mlra)
      
      if (k >= nrow(mlra_candidates)) {
        draw <- mlra_candidates
      } else {
        coords <- mlra_candidates |> dplyr::select(x_cent, y_cent)
        n_unique_coords <- nrow(unique(coords))
        
        # Check for overlapping/duplicate geometries to prevent kmeans failure
        if (k >= n_unique_coords) {
          draw <- mlra_candidates |> dplyr::slice_sample(n = k)
        } else {
          km <- kmeans(coords, centers = k, nstart = 10)
          mlra_candidates$cluster <- km$cluster
          
          draw <- mlra_candidates |>
            dplyr::group_by(cluster) |>
            dplyr::slice_sample(n = 1) |>
            dplyr::ungroup() |>
            dplyr::select(-cluster)
        }
      }
      
      selected_sites[[length(selected_sites) + 1]] <- draw
      selected_site_ids <- c(selected_site_ids, draw$id)
      global_mlra_counts[t_mlra] <- global_mlra_counts[t_mlra] + nrow(draw)
    }
  }
  
  final_draw <- dplyr::bind_rows(selected_sites) |> 
    dplyr::select(-x_cent, -y_cent) 
  
  return(final_draw)
}


# --- 2. BBOX GENERATION & METRICS EXTRACTION ----------------------------------
message("\n--- Preparing Spatial Data & Metrics ---")

# Load Systematic Sample Sites
systematicSamples <- list.files("data/products/systematicSampleSelection/", full.names = TRUE)
siteIDs <- systematicSamples[grepl(pattern = paste0(TARGET_LLR,"_draw_", TARGET_SIZE), x = systematicSamples)] |>
  readr::read_csv(show_col_types = FALSE)

# Generate Bounding Boxes
if (!file.exists(siteExport)) {
  unique_sites <- unique(siteIDs$id) 
  future::plan(future::multisession, workers = max(1, future::availableCores() - 10))
  
  bboxes <- furrr::future_map_dfr(unique_sites, function(current_id) {
    aoi_proj <- sf::st_transform(getAOI(grid100 = g100, id = current_id), crs = 4326)
    bbox <- sf::st_bbox(aoi_proj)
    data.frame(
      id = current_id, xmin = as.numeric(bbox["xmin"]), ymin = as.numeric(bbox["ymin"]),
      xmax = as.numeric(bbox["xmax"]), ymax = as.numeric(bbox["ymax"]), stringsAsFactors = FALSE
    )
  }, .options = furrr::furrr_options(seed = TRUE))
  
  future::plan(future::sequential)
  
  sites_with_bbox <- siteIDs |> dplyr::left_join(bboxes, by = "id")
  readr::write_csv(sites_with_bbox, siteExport)
} else {
  sites_with_bbox <- readr::read_csv(siteExport, show_col_types = FALSE)
}

# Extract NLCD & TCC Metrics
if (!file.exists(metricsExport)) {
  dir.create(dirname(metricsExport), recursive = TRUE, showWarnings = FALSE)
  
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
  
  future::plan(future::multisession, workers = max(1, future::availableCores() - 10))
  tictoc::tic("Metrics Extraction Time")
  
  spatial_results <- furrr::future_pmap_dfr(
    list(current_id = sites_with_bbox$id, xmin = sites_with_bbox$xmin, ymin = sites_with_bbox$ymin, xmax = sites_with_bbox$xmax, ymax = sites_with_bbox$ymax),
    .f = safe_metrics_worker, 
    nlcd_paths = nlcd_files, tcc_paths = tcc_files,
    .options = furrr::furrr_options(seed = TRUE, packages = c("terra", "sf"))
  )
  
  future::plan(future::sequential)
  tictoc::toc()
  
  final_sites <- sites_with_bbox |> dplyr::left_join(spatial_results, by = "id")
  readr::write_csv(final_sites, metricsExport)
} else {
  final_sites <- readr::read_csv(metricsExport, show_col_types = FALSE)
}


# --- 3. DYNAMIC STRATIFICATION ------------------------------------------------
message("\n--- Calculating Dynamic Thresholds ---")

# Define your target year for the sample
target_year <- "2020"

# Pivot dataset and isolate the single target year
long_sites <- final_sites |>
  tidyr::pivot_longer(
    cols = matches("p_forest_|mean_tcc_"),
    names_to = c(".value", "year"),
    names_pattern = "(.*)_(\\d{4})"
  ) |>
  dplyr::filter(year == target_year) |>  # <-- Apply the temporal filter here
  dplyr::mutate(site_year_id = paste(id, year, sep = "_")) |>
  dplyr::filter(!is.na(p_forest) & !is.na(mean_tcc))

# Calculate optimal k-means breaks for forest cover > 0

# Calculate optimal k-means breaks for forest cover > 0
non_zero_forest <- long_sites |> dplyr::filter(p_forest > 0) |> dplyr::pull(p_forest)

set.seed(2026) 
km_res <- kmeans(non_zero_forest, centers = 4, nstart = 25)
forest_breaks <- tapply(non_zero_forest, km_res$cluster, max) |> sort()

t1 <- forest_breaks[1]
t2 <- forest_breaks[2]
t3 <- forest_breaks[3]

message(sprintf("Calculated Forest Thresholds: 0, %.2f, %.2f, %.2f", t1, t2, t3))

# Assign final classes
classified_pop <- long_sites |>
  dplyr::mutate(
    tcc_class = dplyr::case_when(
      mean_tcc == 0 ~ 0, mean_tcc > 0 & mean_tcc <= 10 ~ 1,
      mean_tcc > 10 & mean_tcc <= 25 ~ 2, mean_tcc > 25 & mean_tcc <= 50 ~ 3,
      mean_tcc > 50 ~ 4, TRUE ~ NA_real_
    ),
    forest_class = dplyr::case_when(
      p_forest == 0 ~ 0, p_forest > 0 & p_forest <= t1 ~ 1,
      p_forest > t1 & p_forest <= t2 ~ 2, p_forest > t2 & p_forest <= t3 ~ 3,
      p_forest > t3 ~ 4, TRUE ~ NA_real_
    )
  ) |>
  dplyr::filter(!is.na(tcc_class) & !is.na(forest_class))


# --- 4. ALLOCATION TARGETS ----------------------------------------------------
# Define targets across the three sampling tracks
alloc_neyman_forest <- calc_neyman_with_baseline(classified_pop, "forest_class", "p_forest", TARGET_TOTAL_SAMPLE, 10)
alloc_neyman_tcc <- calc_neyman_with_baseline(classified_pop, "tcc_class", "mean_tcc", TARGET_TOTAL_SAMPLE, 10)
alloc_equal_tcc <- data.frame(tcc_class = sort(unique(classified_pop$tcc_class)), target = TARGET_TOTAL_SAMPLE / 5)


# --- 5. SAMPLE EXTRACTION -----------------------------------------------------
message("\n--- Executing Sample Draws ---")
set.seed(2026)

sample_neyman_forest <- draw_balanced_spatial_sample(classified_pop, alloc_neyman_forest, "forest_class")
sample_neyman_tcc <- draw_balanced_spatial_sample(classified_pop, alloc_neyman_tcc, "tcc_class")
sample_equal_tcc <- draw_balanced_spatial_sample(classified_pop, alloc_equal_tcc, "tcc_class")


# --- 6. EXPORT FINAL ALLOCATIONS ----------------------------------------------
message("\n--- Exporting Final Datasets ---")

export_dir <- "data/products/groundTruthSamples"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(sample_neyman_forest, file.path(export_dir, "sample_neyman_forest_200.csv"))
readr::write_csv(sample_neyman_tcc, file.path(export_dir, "sample_neyman_tcc_200.csv"))
readr::write_csv(sample_equal_tcc, file.path(export_dir, "sample_equal_tcc_200.csv"))


# --- 7. LEAFLET MAP GENERATION ------------------------------------------------
message("\n--- Generating Leaflet Map ---")

# Convert bounding box coordinates to sf point geometries
sites_sf <- sample_neyman_forest |>
  dplyr::mutate(lon = (xmin + xmax) / 2, lat = (ymin + ymax) / 2) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Load and prepare MLRA boundaries
mlra_sf <- sf::st_read(lrr_id_path, quiet = TRUE) |> 
  dplyr::filter(LRRSYM == "G") |>
  sf::st_transform(4326)

# Prep Custom HTML Info Panel
class_descriptions <- c(
  "0" = "0%", "1" = sprintf(">0 to %.1f%%", t1), "2" = sprintf(">%.1f to %.1f%%", t1, t2),
  "3" = sprintf(">%.1f to %.1f%%", t2, t3), "4" = sprintf(">%.1f%%", t3)
)

summary_df <- data.frame(forest_class = as.numeric(names(class_descriptions)), desc = class_descriptions) |>
  dplyr::left_join(sites_sf |> sf::st_drop_geometry() |> dplyr::count(forest_class), by = "forest_class") |>
  dplyr::mutate(n = ifelse(is.na(n), 0, n))

legend_rows <- paste0("<div style='margin-bottom: 6px;'><b>Class ", summary_df$forest_class, "</b> (", summary_df$desc, "): ", summary_df$n, " sites</div>")
custom_html_panel <- paste0(
  "<div style='background: white; padding: 16px; border-radius: 5px; box-shadow: 0 0 15px rgba(0,0,0,0.2); font-size: 14px;'>",
  "<div style='font-size: 16px; font-weight: bold; margin-bottom: 4px;'>Neyman Allocation Data</div>",
  "<div style='font-style: italic; color: #555; margin-bottom: 12px;'>Forest Cover Thresholds</div>",
  "<div style='border-top: 1px solid #ccc; margin-bottom: 12px;'></div>",
  paste(legend_rows, collapse = ""),
  "</div>"
)

# Color Palettes
mlra_pal <- leaflet::colorFactor(palette = "Set3", domain = mlra_sf$MLRA_NAME)
discrete_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
pal <- leaflet::colorFactor(palette = discrete_colors, domain = sites_sf$forest_class)

# Generate Map Base
map <- leaflet::leaflet() |>
  leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron, group = "Light Map (Default)") |>
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Imagery (Satellite)") |>
  leaflet::addPolygons(
    data = mlra_sf, color = "#444444", weight = 1, fillColor = ~mlra_pal(MLRA_NAME),
    fillOpacity = 0.4, label = ~MLRA_NAME,
    highlightOptions = leaflet::highlightOptions(weight = 3, color = "#000000", fillOpacity = 0.6, bringToFront = FALSE)
  )

# Add Layered Markers for Toggling
classes <- sort(unique(sites_sf$forest_class))
for (cls in classes) {
  class_data <- sites_sf |> dplyr::filter(forest_class == cls)
  map <- map |>
    leaflet::addCircleMarkers(
      data = class_data, group = paste("Class", cls), radius = 5,
      color = ~pal(forest_class), stroke = TRUE, weight = 1, fillOpacity = 0.9,
      popup = ~paste0("<b>Site ID:</b> ", id, "<br><b>MLRA ID:</b> ", MLRA_ID, "<br><b>Forest Class:</b> ", forest_class, "<br><b>Year:</b> ", year),
      label = ~paste("Class:", forest_class)
    )
}

# Add Controls and Legend
map <- map |>
  leaflet::addControl(html = custom_html_panel, position = "topright") |>
  leaflet::addLegend(position = "bottomright", pal = pal, values = sites_sf$forest_class, title = "Forest Class", opacity = 1) |>
  leaflet::addLayersControl(
    baseGroups = c("Light Map (Default)", "Imagery (Satellite)"),
    overlayGroups = paste("Class", classes),
    options = leaflet::layersControlOptions(collapsed = FALSE),
    position = "topleft" 
  )

map