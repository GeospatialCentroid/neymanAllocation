# ------------------------------------------------------------------------------
# 1. Package Loading & Configuration
# ------------------------------------------------------------------------------
pacman::p_load(dplyr, terra, sf, readr, stringr, tidyr, stats)

KS_THRESHOLD <- 0.15 
EVAL_BANDS <- c("Blue", "NIR")

NAIP_DIR <- "/mnt/unraid_naip/LLR G"
ARD_DIR  <- file.path(NAIP_DIR, "ARD")

# ------------------------------------------------------------------------------
# 2. Histogram Specification Function
# ------------------------------------------------------------------------------
match_histograms <- function(target_rast, ref_rast, n_quantiles = 1000) {
  probs <- seq(0, 1, length.out = n_quantiles)
  matched_layers <- list()
  
  for (i in seq_len(nlyr(target_rast))) {
    tgt_vals <- terra::spatSample(target_rast[[i]], size = 50000, method = "regular", na.rm = TRUE)[[1]]
    ref_vals <- terra::spatSample(ref_rast[[i]], size = 50000, method = "regular", na.rm = TRUE)[[1]]
    
    tgt_q <- stats::quantile(tgt_vals, probs = probs, na.rm = TRUE)
    ref_q <- stats::quantile(ref_vals, probs = probs, na.rm = TRUE)
    
    unique_idx <- !duplicated(tgt_q)
    tgt_q <- tgt_q[unique_idx]
    ref_q <- ref_q[unique_idx]
    
    matched_layer <- terra::app(target_rast[[i]], fun = function(x) {
      stats::approx(x = tgt_q, y = ref_q, xout = x, rule = 2)$y
    })
    
    names(matched_layer) <- names(target_rast)[i]
    matched_layers[[i]] <- matched_layer
  }
  return(terra::rast(matched_layers))
}

# ------------------------------------------------------------------------------
# 3. Dynamic Grid Processing Function
# ------------------------------------------------------------------------------
process_grid_files <- function(files, grid_id, ard_dir = ARD_DIR) {
  grid_ard_dir <- file.path(ard_dir, grid_id)
  
  # Check if all 3 expected files/symlinks already exist in the ARD folder
  # If they do, we skip processing entirely to save time
  if (dir.exists(grid_ard_dir)) {
    existing_outputs <- list.files(grid_ard_dir, pattern = "\\.tif$")
    if (length(existing_outputs) >= 3) {
      cat(sprintf("[%s] -> Already fully processed in ARD. Skipping.\n", grid_id))
      return(invisible(TRUE))
    }
  }
  
  # Create directory if it's missing or partially filled
  dir.create(grid_ard_dir, recursive = TRUE, showWarnings = FALSE)
  
  raw_list <- list()
  sample_list <- list()
  file_map <- list()
  
  for (f in files) {
    # Extract the 4 digits immediately before .tif at the end of the filename
    yr <- stringr::str_extract(basename(f), "\\d{4}(?=\\.tif$)")
    if (is.na(yr)) next
    
    r <- terra::rast(f)
    names(r) <- c("Red", "Green", "Blue", "NIR")
    raw_list[[yr]] <- r
    file_map[[yr]] <- f
    sample_list[[yr]] <- terra::spatSample(r[[EVAL_BANDS]], size = 10000, method = "regular", na.rm = TRUE)
  }
  
  years <- sort(names(raw_list))
  if (length(years) < 3) {
    warning(sprintf("[%s] Skipping: Insufficient unique year layers found.", grid_id))
    return(invisible(NULL))
  }
  
  # Pairwise KS evaluation
  pairs <- combn(years, 2, simplify = FALSE)
  dist_matrix <- matrix(0, nrow = length(years), ncol = length(years), dimnames = list(years, years))
  
  for (p in pairs) {
    yr1 <- p[1]; yr2 <- p[2]
    max_d <- 0
    for (b in EVAL_BANDS) {
      d <- stats::ks.test(sample_list[[yr1]][[b]], sample_list[[yr2]][[b]])$statistic
      if (d > max_d) max_d <- d
    }
    dist_matrix[yr1, yr2] <- max_d
    dist_matrix[yr2, yr1] <- max_d
  }
  
  scores <- rowSums(dist_matrix)
  outlier_yr <- names(which.max(scores))
  consensus_yrs <- setdiff(years, outlier_yr)
  
  consensus_dist <- dist_matrix[consensus_yrs[1], consensus_yrs[2]]
  outlier_dist <- mean(c(dist_matrix[outlier_yr, consensus_yrs[1]], dist_matrix[outlier_yr, consensus_yrs[2]]))
  
  cat(sprintf("\n=== Processing Grid: %s ===\n", grid_id))
  cat(sprintf("Consensus Pair (%s vs %s) D-stat: %.4f\n", consensus_yrs[1], consensus_yrs[2], consensus_dist))
  cat(sprintf("Candidate Outlier (%s) Avg D-stat: %.4f\n", outlier_yr, outlier_dist))
  
  ref_yr <- max(consensus_yrs)
  
  for (yr in years) {
    target_filename <- paste0(grid_id, "_", yr, ".tif")
    target_path <- file.path(grid_ard_dir, target_filename)
    
    if (file.exists(target_path) || nzchar(Sys.readlink(target_path))) {
      file.remove(target_path)
    }
    
    if (yr == outlier_yr && outlier_dist > KS_THRESHOLD && outlier_dist > (2 * consensus_dist)) {
      cat(sprintf(" -> [NORMALIZE] Generating matched raster for %s (Reference: %s)...\n", yr, ref_yr))
      norm_rast <- match_histograms(raw_list[[yr]], raw_list[[ref_yr]])
      terra::writeRaster(norm_rast, target_path, overwrite = TRUE)
    } else {
      cat(sprintf(" -> [SYMLINK] Linking raw raster for stable year %s...\n", yr))
      file.symlink(from = file_map[[yr]], to = target_path)
    }
  }
  return(invisible(TRUE))
}

# ------------------------------------------------------------------------------
# 4. Recursive Discovery & Dynamic Grouping
# ------------------------------------------------------------------------------
cat("Scanning for NAIP imagery collections inside LLR G...\n")

# Recursively locate all TIF files, tracking down both flat grid folders and batch folders
all_tif_files <- list.files(NAIP_DIR, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
all_tif_files <- all_tif_files[!grepl("modelOutputs|metadata|/ARD/", all_tif_files)]

# Dynamically parse out the true grid ID from the file name instead of using folder name
# Matches the specific format: 4 digits, dash, 1 digit, dash, up to two letters/numbers per group
grid_mapping <- data.frame(file_path = all_tif_files) %>%
  mutate(
    file_name = basename(file_path),
    grid_id = str_extract(file_name, "\\d{4}-\\d-\\w+-\\w+-\\d")
  ) %>%
  filter(!is.na(grid_id))

# Group files together by their parsed true Grid ID
grouped_grids <- split(grid_mapping$file_path, grid_mapping$grid_id)

cat(sprintf("Discovered %d distinct grids across various batch and site folders.\n", length(grouped_grids)))

# Execute over sorted true Grid IDs
for (grid_id in sort(names(grouped_grids))) {
  files_to_process <- grouped_grids[[grid_id]]
  
  tryCatch({
    process_grid_files(files_to_process, grid_id)
  }, error = function(e) {
    cat(sprintf(" [ERROR] Failed processing grid %s: %s\n", grid_id, e$message))
  })
}

cat("\nARD Build Run Complete.\n")