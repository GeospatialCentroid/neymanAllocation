# ------------------------------------------------------------------------------
# 1. Package Loading & Configuration
# ------------------------------------------------------------------------------
pacman::p_load(dplyr, terra, sf, readr, stringr, tidyr, stats)

# Define core evaluation parameters
EVAL_BANDS <- c("Blue", "NIR") #[cite: 3]
SAMPLE_SIZE <- 50000           #[cite: 3]
N_QUANTILES <- 1000            #[cite: 3]

# ------------------------------------------------------------------------------
# 2. Harmonization Engines: Sampled vs. Full Population
# ------------------------------------------------------------------------------

# Method A: Memory-Optimized (Sampled Population)[cite: 3]
match_histograms_sampled <- function(target_rast, ref_rast, n_quantiles = N_QUANTILES, sample_size = SAMPLE_SIZE) {
  probs <- seq(0, 1, length.out = n_quantiles)
  matched_layers <- list()
  
  for (i in seq_len(nlyr(target_rast))) {
    tgt_vals <- terra::spatSample(target_rast[[i]], size = sample_size, method = "regular", na.rm = TRUE)[[1]] #[cite: 3]
    ref_vals <- terra::spatSample(ref_rast[[i]], size = sample_size, method = "regular", na.rm = TRUE)[[1]] #[cite: 3]
    
    tgt_q <- stats::quantile(tgt_vals, probs = probs, na.rm = TRUE) #[cite: 3]
    ref_q <- stats::quantile(ref_vals, probs = probs, na.rm = TRUE) #[cite: 3]
    
    unique_idx <- !duplicated(tgt_q) #[cite: 3]
    tgt_q <- tgt_q[unique_idx] #[cite: 3]
    ref_q <- ref_q[unique_idx] #[cite: 3]
    
    rm(tgt_vals, ref_vals) #[cite: 3]
    
    matched_layer <- terra::app(target_rast[[i]], fun = function(x) { #[cite: 3]
      stats::approx(x = tgt_q, y = ref_q, xout = x, rule = 2)$y #[cite: 3]
    })
    
    names(matched_layer) <- names(target_rast)[i] #[cite: 3]
    matched_layers[[i]] <- matched_layer #[cite: 3]
  }
  
  out_rast <- terra::rast(matched_layers) #[cite: 3]
  rm(matched_layers) #[cite: 3]
  return(out_rast) #[cite: 3]
}

# Method B: Full Population (Unsampled)
match_histograms_full <- function(target_rast, ref_rast, n_quantiles = N_QUANTILES) {
  probs <- seq(0, 1, length.out = n_quantiles)
  matched_layers <- list()
  
  for (i in seq_len(nlyr(target_rast))) {
    # Extract all non-NA grid cells across the full raster array
    tgt_vals <- terra::values(target_rast[[i]], mat = FALSE, na.rm = TRUE)
    ref_vals <- terra::values(ref_rast[[i]], mat = FALSE, na.rm = TRUE)
    
    tgt_q <- stats::quantile(tgt_vals, probs = probs, na.rm = TRUE)
    ref_q <- stats::quantile(ref_vals, probs = probs, na.rm = TRUE)
    
    unique_idx <- !duplicated(tgt_q)
    tgt_q <- tgt_q[unique_idx]
    ref_q <- ref_q[unique_idx]
    
    rm(tgt_vals, ref_vals)
    
    matched_layer <- terra::app(target_rast[[i]], fun = function(x) {
      stats::approx(x = tgt_q, y = ref_q, xout = x, rule = 2)$y
    })
    
    names(matched_layer) <- names(target_rast)[i]
    matched_layers[[i]] <- matched_layer
  }
  
  out_rast <- terra::rast(matched_layers)
  rm(matched_layers)
  return(out_rast)
}

# ------------------------------------------------------------------------------
# 3. Time Test & Quantitative Distribution Evaluation
# ------------------------------------------------------------------------------
evaluate_method_differences <- function(target_rast, ref_rast, grid_id, target_year, ref_year) {
  
  # 1. Execute and Time Sampled Approach
  time_sampled <- system.time({
    out_sampled <- match_histograms_sampled(target_rast, ref_rast)
  })["elapsed"]
  
  # 2. Execute and Time Full Population Approach
  time_full <- system.time({
    out_full <- match_histograms_full(target_rast, ref_rast)
  })["elapsed"]
  
  # 3. Compute Quantitative Metrics Across Evaluation Bands
  band_metrics <- lapply(EVAL_BANDS, function(band) {
    # Extract full values for robust statistical comparison
    v_target  <- terra::values(target_rast[[band]], mat = FALSE, na.rm = TRUE)
    v_ref     <- terra::values(ref_rast[[band]], mat = FALSE, na.rm = TRUE)
    v_sampled <- terra::values(out_sampled[[band]], mat = FALSE, na.rm = TRUE)
    v_full    <- terra::values(out_full[[band]], mat = FALSE, na.rm = TRUE)
    
    # A. Distribution Divergence: Kolmogorov-Smirnov D-statistic between outputs
    ks_stat <- stats::ks.test(v_sampled, v_full)$statistic
    
    # B. Quantile RMSE between Sampled and Full Harmonized Outputs
    q_probs <- seq(0, 1, by = 0.01)
    q_samp  <- stats::quantile(v_sampled, probs = q_probs, na.rm = TRUE)
    q_full  <- stats::quantile(v_full, probs = q_probs, na.rm = TRUE)
    q_rmse  <- sqrt(mean((q_full - q_samp)^2))
    
    # C. Variability Analysis: Variance differentials relative to the Reference image
    var_ref     <- stats::var(v_ref)
    var_target  <- stats::var(v_target)
    var_sampled <- stats::var(v_sampled)
    var_full    <- stats::var(v_full)
    
    data.frame(
      grid_id = grid_id,
      target_year = target_year,
      ref_year = ref_year,
      band = band,
      time_sampled_sec = as.numeric(time_sampled),
      time_full_sec = as.numeric(time_full),
      time_diff_sec = as.numeric(time_full - time_sampled),
      ks_d_stat_outputs = as.numeric(ks_stat),
      quantile_rmse = as.numeric(q_rmse),
      var_initial_diff = as.numeric(var_target - var_ref),
      var_sampled_residual_diff = as.numeric(var_sampled - var_ref),
      var_full_residual_diff = as.numeric(var_full - var_ref),
      stringsAsFactors = FALSE
    )
  })
  
  rm(out_sampled, out_full)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE)
  
  return(dplyr::bind_rows(band_metrics))
}

# ------------------------------------------------------------------------------
# 4. Harmonization Scaling Evaluation (1 Image vs. 2 Images)
# ------------------------------------------------------------------------------
evaluate_scaling <- function(target_rasts, ref_rast, grid_id) {
  # Validate that at least two targets are available for the scaling test
  if (length(target_rasts) < 2) {
    warning(sprintf("[%s] Insufficient target images for 2-image scaling test.", grid_id))
    return(NULL)
  }
  
  # Test 1: Harmonize a single image to the reference location
  time_one_img <- system.time({
    out_1 <- match_histograms_sampled(target_rasts[[1]], ref_rast)
  })["elapsed"]
  rm(out_1)
  gc(verbose = FALSE)
  
  # Test 2: Harmonize two images sequentially to the single reference location
  time_two_img <- system.time({
    out_a <- match_histograms_sampled(target_rasts[[1]], ref_rast)
    out_b <- match_histograms_sampled(target_rasts[[2]], ref_rast)
  })["elapsed"]
  rm(out_a, out_b)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE)
  
  data.frame(
    grid_id = grid_id,
    time_one_image_sec = as.numeric(time_one_img),
    time_two_images_sec = as.numeric(time_two_img),
    scaling_factor = as.numeric(time_two_img / time_one_img),
    overhead_sec = as.numeric(time_two_img - (2 * time_one_img)),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Grid Processing Wrapper (Single AOI to List of AOIs)
# ------------------------------------------------------------------------------
process_grid_evaluation <- function(grid_id, files) {
  cat(sprintf("\n=== Evaluating Performance for Grid: %s ===\n", grid_id))
  
  # Load available years for the current grid
  raw_list <- list()
  for (f in files) {
    yr <- stringr::str_extract(basename(f), "\\d{4}(?=\\.tif$)") #[cite: 3]
    if (is.na(yr)) next #[cite: 3]
    
    r <- terra::rast(f) #[cite: 3]
    names(r) <- c("Red", "Green", "Blue", "NIR") #[cite: 2, 3]
    raw_list[[yr]] <- r #[cite: 3]
  }
  
  years <- sort(names(raw_list)) #[cite: 3]
  if (length(years) < 2) {
    warning(sprintf("[%s] Skipping: Requires at least 2 years for evaluation.", grid_id))
    return(list(methods = NULL, scaling = NULL))
  }
  
  # Assign the most recent year as the stable reference target
  ref_yr <- max(years)
  target_yrs <- setdiff(years, ref_yr)
  ref_rast <- raw_list[[ref_yr]]
  
  # A. Evaluate Method Differences (Sampled vs Full) on the earliest year
  eval_target_yr <- min(target_yrs)
  cat(sprintf(" -> Testing Sampled vs Full matching on target %s (Ref: %s)...\n", eval_target_yr, ref_yr))
  method_results <- evaluate_method_differences(
    target_rast = raw_list[[eval_target_yr]], 
    ref_rast = ref_rast, 
    grid_id = grid_id, 
    target_year = eval_target_yr, 
    ref_year = ref_yr
  )
  
  # B. Evaluate Scaling (1 Image vs 2 Images) if enough layers exist
  scaling_results <- NULL
  if (length(target_yrs) >= 2) {
    cat(" -> Testing execution scaling (1 vs 2 target images)...\n")
    scaling_targets <- list(raw_list[[target_yrs[1]]], raw_list[[target_yrs[2]]])
    scaling_results <- evaluate_scaling(scaling_targets, ref_rast, grid_id)
  } else {
    cat(" -> Notice: Only 1 target year available; skipping 2-image scaling test.\n")
  }
  
  rm(raw_list, ref_rast)
  gc(verbose = FALSE)
  
  return(list(methods = method_results, scaling = scaling_results))
}

# ------------------------------------------------------------------------------
# 6. Optimized Execution Across a List of AOIs / Grids (Directory-First)
# ------------------------------------------------------------------------------
naip_dir <- "/mnt/unraid_naip/LLR G/" #

# Step 1: Grab top-level directories only (avoids statting millions of files over network)
top_dirs <- list.dirs(naip_dir, recursive = FALSE, full.names = TRUE) #[cite: 2]

# Immediately exclude non-imagery directories to narrow the search scope
search_dirs <- top_dirs[!grepl("modelOutputs|metadata|ARD", top_dirs)] #

# Step 2: Define the targeted file retrieval function (matching Script 33)
get_grid_files <- function(grid_id, directories) {
  # Search ONLY for files matching the specific grid_id pattern within valid directories
  target_files <- list.files(
    directories, 
    pattern = paste0(grid_id, ".*\\.tif$"), #[cite: 2]
    full.names = TRUE, 
    recursive = TRUE
  )
  
  # Ensure residual non-imagery paths are filtered out
  return(target_files[!grepl("modelOutputs|metadata|/ARD/", target_files)]) #[cite: 2, 3]
}

# --- DISCOVERY OPTION A: Using a Pre-defined List of Target AOIs (Script 33 Approach) ---
# If you are pulling target grids from an summary CSV (e.g., target_grid_max, target_grid_mid)
# target_grids <- c(target_grid_max, target_grid_mid, target_grid_min)

# --- DISCOVERY OPTION B: Fast Directory-Only Discovery ---
# If you need to discover all grid_ids dynamically, scan directory nodes instead of file nodes.
# Scanning folders is orders of magnitude faster over NFS/SMB shares than listing .tif files.
cat("Scanning directory tree for grid folders...\n")
all_subdirs <- list.dirs(search_dirs, recursive = TRUE, full.names = TRUE)

# Extract unique grid IDs directly from folder names
discovered_grids <- unique(na.omit(
  stringr::str_extract(basename(all_subdirs), "\\d{4}-\\d-\\w+-\\w+-\\d")
))

# Select targets to process (e.g., taking the first 3 discovered grids for evaluation)
target_grids <- sort(discovered_grids)[1:min(3, length(discovered_grids))]

cat(sprintf("Discovered %d distinct grid folders. Starting evaluation on %d selected targets...\n", 
            length(discovered_grids), length(target_grids)))

# Step 3: Execute evaluation workflow using targeted file queries
evaluation_outputs <- lapply(target_grids, function(gid) {
  # Query files dynamically ONLY for the current grid_id
  grid_files <- get_grid_files(grid_id = gid, directories = search_dirs)
  
  if (length(grid_files) == 0) {
    warning(paste("No NAIP files found for grid:", gid)) #[cite: 2]
    return(NULL)
  }
  
  # Pass the targeted file list directly into our evaluation wrapper
  process_grid_evaluation(grid_id = gid, files = grid_files)
})

# Compile tabular outputs into consolidated dataframes
method_comparison_summary <- dplyr::bind_rows(lapply(evaluation_outputs, `[[`, "methods"))
scaling_comparison_summary <- dplyr::bind_rows(lapply(evaluation_outputs, `[[`, "scaling"))

# Display summary results to console
cat("\n=== Method Comparison Summary (Sampled vs Full) ===\n")
print(method_comparison_summary)

cat("\n=== Execution Scaling Summary (1 vs 2 Images) ===\n")
print(scaling_comparison_summary)