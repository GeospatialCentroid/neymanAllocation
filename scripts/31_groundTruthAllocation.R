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
TARGET_TOTAL_SAMPLE <- 100

# File Paths
lrr_id_path <- "data/derived/mlra/lower48MLRA.gpkg"
mlra_grid_path <- "data/derived/grids/GreatPlains_1km_mlra.gpkg"
g100 <- sf::st_read("data/raw/grid100km_aea.gpkg", quiet = TRUE)
nlcd_paths <- list.files("data/raw/nlcd/", full.names = TRUE)

siteExport <- paste0("data/derived/sampleBBOX/bbox_", TARGET_LLR, "_", TARGET_SIZE,".csv")
metricsExport <- paste0("data/derived/metrics/nlcd_tcc_metrics_", TARGET_LLR, "_", TARGET_SIZE, ".csv")

# --- 1. CORE FUNCTIONS --------------------------------------------------------

#' Calculate Neyman Weights with a Fixed Class-0 Minimum
#' 
#' @description Allocates a total sample size across strata based on the variance 
#' within each stratum (Neyman allocation). If class 0 is present, it receives a
#' fixed minimum first and the remaining sample budget is distributed across the
#' remaining classes using Neyman weights.
#'
#' @param df A dataframe containing the population data.
#' @param class_col A string representing the column name containing the stratification classes.
#' @param var_col A string representing the column name of the continuous variable used to calculate variance.
#' @param target_n Numeric. The total number of samples to allocate across all classes.
#' @param class_0_baseline Numeric. The minimum number of samples guaranteed to class 0. Default is 10.
#'
#' @return A summary dataframe containing the sample size targets (`target`) for each class.
calc_neyman_with_baseline <- function(df, class_col, var_col, target_n, class_0_baseline = 10) {
  
  metrics <- df |>
    dplyr::group_by(!!rlang::sym(class_col)) |>
    dplyr::summarize(
      N_h = dplyr::n(),
      S_h = sd(!!rlang::sym(var_col), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      S_h = ifelse(is.na(S_h) | S_h == 0, 0.001, S_h),
      weight = N_h * S_h,
      neyman_portion = 0,
      target = 0
    )
  
  has_zero_class <- any(metrics[[class_col]] == 0)
  n_non_zero_classes <- sum(metrics[[class_col]] != 0)
  
  if (has_zero_class && n_non_zero_classes == 0) {
    metrics$target[metrics[[class_col]] == 0] <- target_n
    return(metrics)
  }
  
  fixed_baseline <- if (has_zero_class) class_0_baseline else 0
  remaining_budget <- target_n - fixed_baseline
  
  if (remaining_budget < 0) {
    stop("Class 0 baseline exceeds the total target sample size.")
  }
  
  if (has_zero_class) {
    metrics$target[metrics[[class_col]] == 0] <- class_0_baseline
    neyman_idx <- metrics[[class_col]] != 0
  } else {
    neyman_idx <- rep(TRUE, nrow(metrics))
  }
  
  neyman_weights <- metrics$weight[neyman_idx]
  metrics$neyman_portion[neyman_idx] <- round(remaining_budget * (neyman_weights / sum(neyman_weights)))
  metrics$target[neyman_idx] <- metrics$target[neyman_idx] + metrics$neyman_portion[neyman_idx]
  
  # Correct rounding disparities to enforce exact target_n
  diff <- target_n - sum(metrics$target)
  if (diff != 0) {
    max_idx <- which(metrics[[class_col]] != 0)[which.max(metrics$neyman_portion[metrics[[class_col]] != 0])]
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
      # kmeans appoarch 
      # if (k >= nrow(mlra_candidates)) {
      #   draw <- mlra_candidates
      # } else {
      #   coords <- mlra_candidates |> dplyr::select(x_cent, y_cent)
      #   n_unique_coords <- nrow(unique(coords))
      #   
      #   # Check for overlapping/duplicate geometries to prevent kmeans failure
      #   if (k >= n_unique_coords) {
      #     draw <- mlra_candidates |> dplyr::slice_sample(n = k)
      #   } else {
      #     km <- kmeans(coords, centers = k, nstart = 25)
      #     mlra_candidates$cluster <- km$cluster
      #     
      #     draw <- mlra_candidates |>
      #       dplyr::group_by(cluster) |>
      #       dplyr::slice_sample(n = 1) |>
      #       dplyr::ungroup() |>
      #       dplyr::select(-cluster)
      #   }
      # }
      ## distributional along the lat long values 
      if (k >= nrow(mlra_candidates)) {
        draw <- mlra_candidates
      } else {
        # 1. Determine how many slices we need on each axis to get at least 'k' bins
        grid_dim <- ceiling(sqrt(k))
        
        # 2. Use ntile to divide the X and Y distributions into equal percentiles
        mlra_candidates <- mlra_candidates |>
          dplyr::mutate(
            x_bin = dplyr::ntile(x_cent, n = grid_dim),
            y_bin = dplyr::ntile(y_cent, n = grid_dim),
            spatial_bin = paste(x_bin, y_bin, sep = "_") # Create a unique ID for each intersection
          )
        
        # 3. Randomly select 'k' unique bins from those thatactually contain candidate sites
        populated_bins <- unique(mlra_candidates$spatial_bin)
        
        # If we have more populated bins than we need, sample exactly k bins
        if (length(populated_bins) > k) {
          selected_bins <- sample(populated_bins, k)
        } else {
          # If the MLRA shape creates fewer populated bins than k, use them all
          # and we will over-sample some bins to hit 'k'
          selected_bins <- populated_bins
        }
        
        # 4. Draw sites based on the selected bins
        draw <- mlra_candidates |>
          dplyr::filter(spatial_bin %in% selected_bins) |>
          dplyr::group_by(spatial_bin) |>
          # If we have fewer bins than k, slice_sample needs to pull more than 1 from some bins
          # We use a weight or simple slice to hit the exact k target
          dplyr::slice_sample(n = ceiling(k / length(selected_bins))) |> 
          dplyr::ungroup() |>
          dplyr::slice_sample(n = k) |> # Final strict crop to ensure we return exactly k sites
          dplyr::select(-x_bin, -y_bin, -spatial_bin)
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

#' Tag Validation Splits (Stratified)
#'
#' @description Assigns a training or validation tag to rows in a dataframe, 
#' ensuring that the split fraction is applied proportionally across each class.
#'
#' @param df A dataframe containing the sampled sites.
#' @param class_col A string specifying the column with the stratification class.
#' @param val_fraction Numeric value between 0 and 1 indicating the proportion of validation sites.
#'
#' @return The original dataframe with an appended `split_tag` column.
tag_validation_splits <- function(df, class_col, val_fraction = 0.50) {
  
  df |>
    dplyr::group_by(!!rlang::sym(class_col)) |>
    dplyr::mutate(
      # Determine the exact number of validation sites needed for this specific class
      n_val = round(dplyr::n() * val_fraction),
      
      # Create a vector of exact length containing the required number of Validation/Training tags,
      # then shuffle it using sample() so the tags are randomly assigned to the rows.
      split_tag = sample(
        c(rep("Validation", n_val[1]), 
          rep("Training", dplyr::n() - n_val[1]))
      )
    ) |>
    dplyr::ungroup()
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

# --- 3. SPATIAL METRICS EXTRACTION (NLCD & TCC) -------------------------------
message("\n--- Starting Parallel NLCD & TCC Extraction ---")

extract_metrics_worker <- function(current_id, xmin, ymin, xmax, ymax, nlcd_paths, tcc_paths) {
  require(terra, quietly = TRUE)
  require(sf, quietly = TRUE)
  
  bbox_4326 <- sf::st_bbox(
    c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    crs = 4326
  )
  poly_4326 <- sf::st_as_sfc(bbox_4326)
  
  results <- data.frame(id = current_id, stringsAsFactors = FALSE)
  target_years <- unique(c(names(nlcd_paths), names(tcc_paths)))
  
  for (year in target_years) {
    n_path <- nlcd_paths[[year]]
    t_path <- tcc_paths[[year]]
    
    if (!is.null(n_path) && !is.na(n_path) && file.exists(n_path)) {
      r_nlcd <- terra::rast(n_path)
      poly_native_n <- sf::st_transform(poly_4326, crs = terra::crs(r_nlcd))
      vect_native_n <- terra::vect(poly_native_n)
      
      nlcd_cropped <- terra::crop(r_nlcd, vect_native_n, mask = TRUE)
      freq_table <- terra::freq(nlcd_cropped)
      
      if (nrow(freq_table) == 0) {
        results[[paste0("p_forest_", year)]] <- 0
      } else {
        forest_pixels <- sum(freq_table$count[freq_table$value %in% c(41, 42, 43)], na.rm = TRUE)
        total_pixels <- sum(freq_table$count, na.rm = TRUE)
        results[[paste0("p_forest_", year)]] <- (forest_pixels / total_pixels) * 100
      }
    } else {
      results[[paste0("p_forest_", year)]] <- NA
    }
    
    if (!is.null(t_path) && !is.na(t_path) && file.exists(t_path)) {
      r_tcc <- terra::rast(t_path)
      poly_native_t <- sf::st_transform(poly_4326, crs = terra::crs(r_tcc))
      vect_native_t <- terra::vect(poly_native_t)
      
      tcc_cropped <- terra::crop(r_tcc, vect_native_t, mask = TRUE)
      mean_tcc <- terra::global(tcc_cropped, "mean", na.rm = TRUE)[[1]]
      results[[paste0("mean_tcc_", year)]] <- if (is.nan(mean_tcc)) 0 else mean_tcc
    } else {
      results[[paste0("mean_tcc_", year)]] <- NA
    }
  }
  
  return(results)
}

safe_metrics_worker <- purrr::possibly(
  .f = extract_metrics_worker,
  otherwise = function(current_id, ...) {
    data.frame(
      id = current_id,
      p_forest_2010 = NA, p_forest_2016 = NA, p_forest_2020 = NA,
      mean_tcc_2010 = NA, mean_tcc_2016 = NA, mean_tcc_2020 = NA,
      stringsAsFactors = FALSE
    )
  },
  quiet = FALSE
)

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
  
  workers <- max(1, future::availableCores() - 10)
  future::plan(future::multisession, workers = workers)
  tictoc::tic("Total Metrics Extraction Time")
  
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
  message("NLCD and TCC metrics extraction complete.")
} else {
  final_sites <- readr::read_csv(metricsExport, show_col_types = FALSE)
}


# --- 3 & 4 & 5. FLEXIBLE SCENARIO GENERATION ----------------------------------
message("\n--- Running Allocation Scenarios ---")

target_year <- "2020"

# Pivot dataset and isolate the target year
long_sites <- final_sites |>
  tidyr::pivot_longer(
    cols = matches("p_forest_|mean_tcc_"),
    names_to = c(".value", "year"),
    names_pattern = "(.*)_(\\d{4})"
  ) |>
  dplyr::filter(year == target_year) |>
  dplyr::mutate(site_year_id = paste(id, year, sep = "_")) |>
  dplyr::filter(!is.na(p_forest) & !is.na(mean_tcc))

#' Run Allocation Scenario
#' @param data The baseline long_sites dataframe.
#' @param filter_col The column to apply exclusion thresholds against.
#' @param threshold The maximum allowable value for the filter_col (e.g., 0, 20, 100). NA for no filter.
#' @param stratify_col The column used to calculate Neyman variance and k-means breaks.
run_allocation_scenario <- function(data, filter_col, threshold, stratify_col, target_n = TARGET_TOTAL_SAMPLE) {
  
  # 1. Apply Threshold Filters
  if (!is.na(threshold)) {
    pop_data <- data |> dplyr::filter(!!rlang::sym(filter_col) <= threshold)
  } else {
    pop_data <- data
  }
  
  if (nrow(pop_data) == 0) {
    warning("Filter resulted in 0 candidate sites.")
    return(NULL)
  }
  
  # 2. Check Variance for Stratification
  var_val <- var(pop_data[[stratify_col]], na.rm = TRUE)
  
  if (is.na(var_val) || var_val == 0) {
    message(sprintf("Zero variance in %s after filtering. Defaulting to single-class balanced draw.", stratify_col))
    
    classified_pop <- pop_data |> 
      dplyr::mutate(scenario_class = ifelse(!!rlang::sym(stratify_col) == 0, 0, 1))
    
    alloc_targets <- data.frame(
      scenario_class = unique(classified_pop$scenario_class),
      target = target_n
    )
    
  } else {
    # 3. Dynamic K-Means Breaks
    non_zero_vals <- pop_data |> 
      dplyr::filter(!!rlang::sym(stratify_col) > 0) |> 
      dplyr::pull(!!rlang::sym(stratify_col))
    
    # Protect against too few unique values for k=2
    k_centers <- min(2, length(unique(non_zero_vals))) 
    
    if (k_centers > 1) {
      km_res <- kmeans(non_zero_vals, centers = k_centers, nstart = 25)
      breaks <- tapply(non_zero_vals, km_res$cluster, max) |> sort()
      
      classified_pop <- pop_data |>
        dplyr::mutate(
          scenario_class = dplyr::case_when(
            !!rlang::sym(stratify_col) == 0 ~ 0,
            !!rlang::sym(stratify_col) > 0 & !!rlang::sym(stratify_col) <= breaks[1] ~ 1,
            length(breaks) >= 2 & !!rlang::sym(stratify_col) > breaks[1] ~ 2,
            TRUE ~ NA_real_
          )
        )
    } else {
      classified_pop <- pop_data |> 
        dplyr::mutate(scenario_class = ifelse(!!rlang::sym(stratify_col) == 0, 0, 1))
    }
    
    # 4. Neyman Allocation targets
    alloc_targets <- tryCatch({
      calc_neyman_with_baseline(classified_pop, "scenario_class", stratify_col, target_n, 10)
    }, error = function(e) {
      warning("Neyman baseline budget exceeded remaining filtered sites. Distributing equally.")
      classes <- unique(classified_pop$scenario_class)
      data.frame(scenario_class = classes, target = floor(target_n / length(classes)))
    })
  }
  
  # 5. Extract Sample using your existing spatial function
  draw_balanced_spatial_sample(classified_pop, alloc_targets, "scenario_class")
}

# --- DEFINE & EXECUTE SCENARIOS ---
# Define a grid of all the iterations you want to run
scenarios <- tibble::tribble(
  ~scenario_name,                ~filter_col,   ~threshold, ~stratify_col,
  "base_neyman_forest",          NA,            NA,         "p_forest",
  "base_neyman_tcc",             NA,            NA,         "mean_tcc",
  "forest_max_0_strat_forest",   "p_forest",    0,          "p_forest",
  "forest_max_20_strat_forest",  "p_forest",    20,         "p_forest",
  "forest_max_50_strat_forest",  "p_forest",    50,         "p_forest",
  "forest_max_100_strat_forest", "p_forest",    100,        "p_forest",
  "tcc_max_0_strat_tcc",         "mean_tcc",    0,          "mean_tcc",
  "tcc_max_20_strat_tcc",        "mean_tcc",    20,         "mean_tcc",
  "tcc_max_50_strat_tcc",        "mean_tcc",    50,         "mean_tcc",
  "tcc_max_100_strat_tcc",       "mean_tcc",    100,        "mean_tcc"
)
set.seed(2026)
results_list <- purrr::pmap(scenarios, function(scenario_name, filter_col, threshold, stratify_col) {
  message(paste("Running:", scenario_name))
  result <- run_allocation_scenario(long_sites, filter_col, threshold, stratify_col, TARGET_TOTAL_SAMPLE)
  
  if (!is.null(result)) {
    # Tag the resulting dataframe with the scenario name for easy tracking
    result$scenario <- scenario_name 
  }
  return(result)
})

# Apply Stratified Split (60/40)
results_list <- purrr::map(results_list, function(df) {
  if (!is.null(df)) {
    tag_validation_splits(df, class_col = "scenario_class", val_fraction = 0.60)
  }
})

# Name the list elements using the scenario names from your tribble
# This ensures purrr::iwalk uses the scenario name instead of the index number
names(results_list) <- scenarios$scenario_name

# --- 4. EXPORT DATA -----------------------------------------------------------
message("\n--- Exporting Final Datasets ---")

export_dir <- "data/products/groundTruthSamples/scenarios_2020"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# 4a. Export individual scenario files
number_groups <- 3
purrr::iwalk(results_list, function(df, name) {
  if (!is.null(df)) {
    # Added _3 to denote the 3-class stratification"
    readr::write_csv(df, file.path(export_dir, paste0(name, "_",TARGET_TOTAL_SAMPLE ,"_",number_groups,".csv")))
  }
})

# 4b. Combine all scenarios into a single data frame and export
message("Combining all scenarios into a single dataset...")
all_scenarios_df <- dplyr::bind_rows(results_list)

# Added _3 to the combined export path
combined_export_path <- paste0(export_dir, "/all_scenarios_combined_200_",number_groups,".csv")
readr::write_csv(all_scenarios_df, combined_export_path)

message(sprintf("Combined dataset exported to: %s", combined_export_path))


# 



# --- 5. LEAFLET MAP GENERATION ------------------------------------------------
message("\n--- Generating Leaflet Map ---")

all_sites <- dplyr::bind_rows(results_list) |>
  dplyr::mutate(lon = (xmin + xmax) / 2, lat = (ymin + ymax) / 2) |>
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326)

mlra_sf <- sf::st_read(lrr_id_path, quiet = TRUE) |> 
  dplyr::filter(LRRSYM == TARGET_LLR) |>
  sf::st_transform(4326)

discrete_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
pal <- leaflet::colorFactor(palette = discrete_colors, domain = 0:4)

map <- leaflet::leaflet() |>
  leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron, group = "Light Map") |>
  leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery, group = "Imagery") |>
  leaflet::addPolygons(
    data = mlra_sf, color = "#444444", weight = 1, fillOpacity = 0.1, label = ~MLRA_NAME,
    highlightOptions = leaflet::highlightOptions(weight = 3, color = "#000000", bringToFront = FALSE)
  )

scenario_names <- unique(all_sites$scenario)

for (scen in scenario_names) {
  scen_data <- all_sites |> dplyr::filter(scenario == scen)
  val_data <- scen_data |> dplyr::filter(split_tag == "Validation")
  
  if (nrow(val_data) > 0) {
    map <- map |>
      leaflet::addCircleMarkers(
        data = val_data, group = scen, radius = 8, color = "#222222", 
        stroke = FALSE, fillOpacity = 0.85, options = leaflet::pathOptions(clickable = FALSE) 
      )
  }
  
  map <- map |>
    leaflet::addCircleMarkers(
      data = scen_data, group = scen, radius = 5, color = ~pal(scenario_class),
      stroke = TRUE, weight = 1, fillOpacity = 0.9,
      popup = ~paste0(
        "<b>Site ID:</b> ", id, "<br><b>Split Tag:</b> ", split_tag, 
        "<br><b>Class:</b> ", scenario_class,
        "<br><b>Forest %:</b> ", round(p_forest, 2),
        "<br><b>Mean TCC %:</b> ", round(mean_tcc, 2)
      ),
      label = ~paste0(split_tag, " - Class: ", scenario_class)
    )
}

shape_legend_html <- paste0(
  "<div style='background: white; padding: 12px; border-radius: 4px; box-shadow: 0 0 10px rgba(0,0,0,0.1); font-size: 13px; font-family: sans-serif;'>",
  "<b>Dataset Splits</b><br>",
  "<div style='margin-top: 8px; display: flex; align-items: center;'>",
  "<span style='display:inline-block; width:10px; height:10px; border-radius:50%; background:#888; border:1px solid #333; margin-right:8px;'></span> Training Site (50%)</div>",
  "<div style='margin-top: 8px; display: flex; align-items: center;'>",
  "<span style='display:inline-block; width:16px; height:16px; border-radius:50%; background:#222; display:flex; justify-content:center; align-items:center; margin-right:5px;'>",
  "<span style='display:inline-block; width:10px; height:10px; border-radius:50%; background:#888; border:1px solid #333;'></span></span> Validation Site (50%)</div>",
  "</div>"
)

map <- map |>
  leaflet::addLegend(position = "bottomright", pal = pal, values = 0:4, title = "Allocation Class") |>
  leaflet::addControl(html = shape_legend_html, position = "bottomleft") |>
  leaflet::addLayersControl(
    baseGroups = c("Light Map", "Imagery"),
    overlayGroups = scenario_names,
    options = leaflet::layersControlOptions(collapsed = FALSE),
    position = "topleft" 
  ) |>
  leaflet::hideGroup(scenario_names[-1]) 

map
# --- 6. QUANTITATIVE SUMMARY & PLOTTING ---------------------------------------
message("\n--- Generating Comparison Metrics and Plots ---")

comparison_table <- all_sites |>
  sf::st_drop_geometry() |>
  dplyr::group_by(scenario) |>
  dplyr::summarize(
    total_sites_sampled = dplyr::n(),
    unique_mlras_represented = dplyr::n_distinct(MLRA_ID),
    mean_forest_sampled = mean(p_forest, na.rm = TRUE),
    sd_forest_sampled = sd(p_forest, na.rm = TRUE),
    mean_tcc_sampled = mean(mean_tcc, na.rm = TRUE),
    sd_tcc_sampled = sd(mean_tcc, na.rm = TRUE),
    class_0_count = sum(scenario_class == 0, na.rm = TRUE),
    class_1_count = sum(scenario_class == 1, na.rm = TRUE),
    class_2_count = sum(scenario_class == 2, na.rm = TRUE),
    class_3_count = sum(scenario_class == 3, na.rm = TRUE),
    class_4_count = sum(scenario_class == 4, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::relocate(scenario, total_sites_sampled, unique_mlras_represented)

readr::write_csv(comparison_table, file.path(export_dir, "scenario_comparison_summary.csv"))

scenario_order <- c(
  "base_neyman_forest",
  "forest_max_100_strat_forest",
  "forest_max_50_strat_forest",
  "forest_max_20_strat_forest",
  "forest_max_0_strat_forest",
  "base_neyman_tcc",
  "tcc_max_100_strat_tcc",
  "tcc_max_50_strat_tcc",
  "tcc_max_20_strat_tcc",
  "tcc_max_0_strat_tcc"
)

plot1_data <- comparison_table |>
  dplyr::mutate(
    target_metric = ifelse(grepl("forest", scenario, ignore.case = TRUE), "Forest Cover (%)", "Tree Canopy Cover (%)"),
    mean_val = ifelse(grepl("forest", scenario, ignore.case = TRUE), mean_forest_sampled, mean_tcc_sampled),
    sd_val = ifelse(grepl("forest", scenario, ignore.case = TRUE), sd_forest_sampled, sd_tcc_sampled),
    scenario = factor(scenario, levels = rev(scenario_order)) 
  )

p1 <- ggplot2::ggplot(plot1_data, ggplot2::aes(x = scenario, y = mean_val, fill = target_metric)) +
  ggplot2::geom_col(color = "black", alpha = 0.8) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(0, mean_val - sd_val), ymax = mean_val + sd_val), width = 0.2, alpha = 0.7) +
  ggplot2::coord_flip() +
  ggplot2::facet_wrap(~target_metric, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("Forest Cover (%)" = "#2ca25f", "Tree Canopy Cover (%)" = "#2b8cbe")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Mean and Standard Deviation of Targeted Cover Metric",
    x = "Scenario", y = "Mean Cover % (± SD)"
  ) +
  ggplot2::theme(legend.position = "none")

plot2_data <- comparison_table |>
  dplyr::select(scenario, dplyr::starts_with("class_")) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("class_"),
    names_to = "allocation_class",
    values_to = "count"
  ) |>
  dplyr::mutate(
    allocation_class = gsub("class_|_count", "", allocation_class),
    scenario = factor(scenario, levels = rev(scenario_order))
  )

p2 <- ggplot2::ggplot(plot2_data, ggplot2::aes(x = scenario, y = count, fill = allocation_class)) +
  ggplot2::geom_col(position = "stack", color = "black", linewidth = 0.2) +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c("0" = "#E41A1C", "1" = "#377EB8", "2" = "#4DAF4A", "3" = "#984EA3", "4" = "#FF7F00")) +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Distribution of Sample Sites Across Allocation Classes",
    x = "Scenario", y = "Number of Sites Drawn", fill = "Class"
  )

ggplot2::ggsave(file.path(export_dir, "plot_1_mean_sd_cover.png"), plot = p1, width = 10, height = 6, bg = "white")
ggplot2::ggsave(file.path(export_dir, "plot_2_class_distribution.png"), plot = p2, width = 10, height = 6, bg = "white")

p1
p2



