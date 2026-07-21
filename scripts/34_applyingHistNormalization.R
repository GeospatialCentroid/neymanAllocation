# ------------------------------------------------------------------------------
# 1. Package Loading & Configuration
# ------------------------------------------------------------------------------
pacman::p_load(dplyr, terra, sf, readr, stringr, tidyr, stats, foreach, doParallel)

KS_THRESHOLD <- 0.15 
EVAL_BANDS <- c("Blue", "NIR")

NAIP_DIR <- "/mnt/unraid_naip/LLR G"
ARD_DIR  <- file.path(NAIP_DIR, "ARD")
LOG_DIR  <- file.path(ARD_DIR, "logs")

# Core Limitation for Memory Evaluation
NUM_CORES <- 5
cl <- parallel::makeCluster(NUM_CORES)
doParallel::registerDoParallel(cl)

# ------------------------------------------------------------------------------
# 2. Histogram Specification Function (Memory Optimized)
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
# 3. Dynamic Grid Processing Function (With Pre-Checks & Logging)
# ------------------------------------------------------------------------------
process_grid_files <- function(files, grid_id, ard_dir = ARD_DIR, log_dir = LOG_DIR) {
  grid_ard_dir <- file.path(ard_dir, grid_id)
  dir.create(grid_ard_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Step 1: Pre-Analysis Check for Existing Target Files
  # Parse out years from input file strings before processing to know what to look for
  expected_years <- sapply(files, function(f) stringr::str_extract(basename(f), "\\d{4}(?=\\.tif$)"))
  expected_years <- sort(unname(expected_years[!is.na(expected_years)]))
  
  if (length(expected_years) >= 3) {
    expected_outputs <- file.path(grid_ard_dir, paste0(grid_id, "_", expected_years, ".tif"))
    # Check if all files exist as standard files OR valid symlinks
    exists_check <- sapply(expected_outputs, function(p) file.exists(p) || nzchar(Sys.readlink(p)))
    
    if (all(exists_check)) {
      cat(sprintf("[%s] -> All expected outputs exist in ARD. Skipping analysis entirely.\n", grid_id))
      return(invisible(TRUE))
    }
  }
  
  # Step 2: Load Data Only After Verification Failure
  raw_list <- list()
  sample_list <- list()
  file_map <- list()
  
  for (f in files) {
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
    rm(raw_list, sample_list, file_map)
    gc(verbose = FALSE)
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
  
  rm(sample_list)
  gc(verbose = FALSE)
  
  scores <- rowSums(dist_matrix)
  outlier_yr <- names(which.max(scores))
  consensus_yrs <- setdiff(years, outlier_yr)
  
  consensus_dist <- dist_matrix[consensus_yrs[1], consensus_yrs[2]]
  outlier_dist <- mean(c(dist_matrix[outlier_yr, consensus_yrs[1]], dist_matrix[outlier_yr, consensus_yrs[2]]))
  
  cat(sprintf("\n=== Processing Grid: %s ===\n", grid_id))
  
  ref_yr <- max(consensus_yrs)
  log_entries <- data.frame()
  
  for (yr in years) {
    target_filename <- paste0(grid_id, "_", yr, ".tif")
    target_path <- file.path(grid_ard_dir, target_filename)
    
    if (file.exists(target_path) || nzchar(Sys.readlink(target_path))) {
      file.remove(target_path)
    }
    
    is_outlier <- (yr == outlier_yr && outlier_dist > KS_THRESHOLD && outlier_dist > (2 * consensus_dist))
    
    if (is_outlier) {
      cat(sprintf(" -> [NORMALIZE] Generating matched raster for %s (Reference: %s)...\n", yr, ref_yr))
      norm_rast <- match_histograms(raw_list[[yr]], raw_list[[ref_yr]])
      terra::writeRaster(norm_rast, target_path, overwrite = TRUE)
      rm(norm_rast)
      action_taken <- "NORMALIZED"
    } else {
      cat(sprintf(" -> [SYMLINK] Linking raw raster for stable year %s...\n", yr))
      file.symlink(from = file_map[[yr]], to = target_path)
      action_taken <- "SYMLINK"
    }
    
    # Track actions taken for the audit trail
    log_entries <- rbind(log_entries, data.frame(
      grid_id = grid_id,
      year = yr,
      action = action_taken,
      is_candidate_outlier = (yr == outlier_yr),
      d_stat_vs_ref = dist_matrix[yr, ref_yr],
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Write isolated log file to prevent parallel write conflicts
  if (nrow(log_entries) > 0) {
    readr::write_csv(log_entries, file.path(log_dir, paste0("log_", grid_id, ".csv")))
  }
  
  rm(raw_list, file_map, log_entries)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE)
  
  return(invisible(TRUE))
}

# ------------------------------------------------------------------------------
# 4. Recursive Discovery & Parallel Execution
# ------------------------------------------------------------------------------
cat("Scanning for NAIP imagery collections inside LLR G...\n")

all_tif_files <- list.files(NAIP_DIR, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
all_tif_files <- all_tif_files[!grepl("modelOutputs|metadata|/ARD/", all_tif_files)]

grid_mapping <- data.frame(file_path = all_tif_files) %>%
  mutate(
    file_name = basename(file_path),
    grid_id = str_extract(file_name, "\\d{4}-\\d-\\w+-\\w+-\\d")
  ) %>%
  filter(!is.na(grid_id))

grouped_grids <- split(grid_mapping$file_path, grid_mapping$grid_id)
grid_names <- sort(names(grouped_grids))

cat(sprintf("Discovered %d distinct grids across various batch and site folders.\n", length(grouped_grids)))
cat(sprintf("Starting parallel processing limited to %d cores...\n", NUM_CORES))

foreach(
  grid_id = grid_names,
  .packages = c("dplyr", "terra", "sf", "readr", "stringr", "tidyr", "stats"),
  .errorhandling = "pass"
) %dopar% {
  files_to_process <- grouped_grids[[grid_id]]
  tryCatch({
    process_grid_files(files_to_process, grid_id)
  }, error = function(e) {
    cat(sprintf(" [ERROR] Failed processing grid %s: %s\n", grid_id, e$message))
  })
}

parallel::stopCluster(cl)

# ------------------------------------------------------------------------------
# 5. Compile Normalization Log Records
# ------------------------------------------------------------------------------
cat("\nCompiling normalization log audit trail...\n")
log_files <- list.files(LOG_DIR, pattern = "^log_.*\\.csv$", full.names = TRUE)

if (length(log_files) > 0) {
  compiled_log <- lapply(log_files, readr::read_csv, show_col_types = FALSE) %>% 
    dplyr::bind_rows() %>%
    dplyr::arrange(grid_id, year)
  
  readr::write_csv(compiled_log, file.path(ARD_DIR, "normalization_log.csv"))
  
  # Clean up individual worker log files to keep the directory clean
  file.remove(log_files)
  
  total_normalized <- sum(compiled_log$action == "NORMALIZED")
  cat(sprintf("Audit complete. Saved record of %d normalized layers to %s.\n", 
              total_normalized, file.path(ARD_DIR, "normalization_log.csv")))
} else {
  cat("No new grids were processed during this run; no new log updates required.\n")
}

cat("\nARD Build Run Complete.\n")