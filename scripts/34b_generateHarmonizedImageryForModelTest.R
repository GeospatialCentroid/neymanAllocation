# ------------------------------------------------------------------------------
# 1. Package Loading & Configuration
# ------------------------------------------------------------------------------
pacman::p_load(dplyr, terra, sf, readr, stringr, tidyr, stats, foreach, doParallel)

# Define Storage Endpoints
NAIP_NETWORK_DIR <- "/mnt/unraid_naip/LLR G"
ARD_NETWORK_DIR  <- file.path(NAIP_NETWORK_DIR, "ARD")
LOCAL_EXPORT_DIR <- "data/derived/harmonization_model_ready"

# --- EXECUTION MODE TOGGLE ---
TEST_MODE  <- FALSE   # Set to TRUE for sequential testing; FALSE for full parallel run
TEST_SITES <- 3      # Number of sites to process during testing
MAX_SITES  <- 1000    # Number of sites to process during production run
NUM_CORES  <- 6      # Number of parallel workers for production run

# Processing Parameters for Double Harmonization Generation
SAMPLE_SIZE <- 50000
N_QUANTILES <- 1000

# Create local storage root
dir.create(LOCAL_EXPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Optimized Harmonization Engine with Quantile Caching (For Double Harmonization)
# ------------------------------------------------------------------------------
cache_reference_quantiles <- function(ref_rast, n_quantiles = N_QUANTILES, sample_size = SAMPLE_SIZE) {
  probs <- seq(0, 1, length.out = n_quantiles)
  cached_quantiles <- list()
  
  for (i in seq_len(nlyr(ref_rast))) {
    ref_vals <- terra::spatSample(ref_rast[[i]], size = sample_size, method = "regular", na.rm = TRUE)[[1]]
    ref_q <- stats::quantile(ref_vals, probs = probs, na.rm = TRUE)
    
    unique_idx <- !duplicated(ref_q)
    cached_quantiles[[names(ref_rast)[i]]] <- list(
      probs = probs[unique_idx],
      quantiles = ref_q[unique_idx]
    )
  }
  return(cached_quantiles)
}

apply_cached_harmonization <- function(target_rast, cached_ref, sample_size = SAMPLE_SIZE) {
  matched_layers <- list()
  
  for (i in seq_len(nlyr(target_rast))) {
    band_name <- names(target_rast)[i]
    ref_q_obj <- cached_ref[[band_name]]
    
    tgt_vals <- terra::spatSample(target_rast[[i]], size = sample_size, method = "regular", na.rm = TRUE)[[1]]
    tgt_q <- stats::quantile(tgt_vals, probs = ref_q_obj$probs, na.rm = TRUE)
    
    unique_idx <- !duplicated(tgt_q)
    tgt_q <- tgt_q[unique_idx]
    ref_q <- ref_q_obj$quantiles[unique_idx]
    
    matched_layer <- terra::app(target_rast[[i]], fun = function(x) {
      stats::approx(x = tgt_q, y = ref_q, xout = x, rule = 2)$y
    })
    
    names(matched_layer) <- band_name
    matched_layers[[i]] <- matched_layer
  }
  
  out_rast <- terra::rast(matched_layers)
  rm(matched_layers)
  return(out_rast)
}

# ------------------------------------------------------------------------------
# 3. Fast Network Discovery & In-Memory Pre-Indexing
# ------------------------------------------------------------------------------
cat("Scanning directory tree for physical ARD products...\n")

active_site_limit <- if (TEST_MODE) TEST_SITES else MAX_SITES

# Step 1: Grab top-level AOI folder paths inside the ARD directory
ard_dir_paths <- list.dirs(ARD_NETWORK_DIR, recursive = FALSE, full.names = TRUE)
ard_dir_paths <- ard_dir_paths[stringr::str_detect(basename(ard_dir_paths), "^\\d{4}-\\d-\\w+")]

if (length(ard_dir_paths) == 0) stop("No valid grid folders found in the network ARD directory.")

# Step 2: Check folder contents for physical files present (ignoring symlinks)
valid_grids <- c()
for (dir_path in sort(ard_dir_paths)) {
  tif_files <- list.files(dir_path, pattern = "\\.tif$", recursive = FALSE, full.names = TRUE)
  if (any(!nzchar(Sys.readlink(tif_files)))) {
    valid_grids <- c(valid_grids, basename(dir_path))
  }
  if (length(valid_grids) == active_site_limit) break
}

if (length(valid_grids) == 0) stop("No locations containing physical (non-symlinked) TIFF files found in ARD.")
target_grids <- valid_grids
cat(sprintf("Verified %d qualified ARD sites containing physical harmonized data products.\n", length(target_grids)))

# Step 3: Pre-index top-level directories into memory to prevent repetitive NAS scanning lag
network_top_dirs <- list.dirs(NAIP_NETWORK_DIR, recursive = FALSE, full.names = TRUE)
search_dirs <- network_top_dirs[!grepl("modelOutputs|metadata|ARD", network_top_dirs)]

cat("Pre-indexing network file tree into memory (runs once to eliminate network I/O lag)...\n")
all_network_files <- list.files(search_dirs, recursive = TRUE, full.names = TRUE)
all_network_files <- all_network_files[!grepl("modelOutputs|metadata|/ARD/|/old/", all_network_files, ignore.case = TRUE)]
cat(sprintf("Indexed %d total network files ready for fast in-memory querying.\n", length(all_network_files)))

# ------------------------------------------------------------------------------
# 4. Local Export & Harmonization Generation Engine
# ------------------------------------------------------------------------------
process_local_aoi_export <- function(grid_id, local_root = LOCAL_EXPORT_DIR) {
  
  # FAST IN-MEMORY QUERY: Replace slow recursive disk scans with string matching against our index
  grid_files <- grep(paste0("/", grid_id, ".*\\.tif$"), all_network_files, value = TRUE, ignore.case = TRUE)
  
  file_map <- list()
  for (f in grid_files) {
    yr <- stringr::str_extract(basename(f), "\\d{4}(?=\\.tif$)")
    if (!is.na(yr)) file_map[[yr]] <- f
  }
  
  years <- sort(names(file_map))
  if (length(years) < 3) {
    warning(sprintf("[%s] Skipping: Found only %d imagery years; 3 required.", grid_id, length(years)))
    return(invisible(FALSE))
  }
  
  aoi_local_dir <- file.path(local_root, grid_id)
  dir.create(aoi_local_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat(sprintf("\n=== Processing AOI: %s ===\n", grid_id))
  
  # Lazy-load imagery headers via terra
  raster_list <- list()
  for (yr in years) {
    r <- terra::rast(file_map[[yr]])
    names(r) <- c("Red", "Green", "Blue", "NIR")
    raster_list[[yr]] <- r
  }
  
  # ----------------------------------------------------------------------------
  # 1. EXPORT RAW AS-DOWNLOADED IMAGERY (Direct File Copy Speed)
  # ----------------------------------------------------------------------------
  cat(" -> [1/4] Checking and transferring raw downloaded imagery...\n")
  for (yr in years) {
    raw_out <- file.path(aoi_local_dir, paste0("raw_", yr, ".tif"))
    if (!file.exists(raw_out)) {
      # Direct OS-level byte copy bypasses GDAL memory decoding entirely
      file.copy(from = file_map[[yr]], to = raw_out, overwrite = TRUE)
      cat(sprintf("    [COPY] Transferred raw imagery for %s.\n", yr))
    } else {
      cat(sprintf("    [SKIP] Raw imagery for %s already exists locally.\n", yr))
    }
  }
  
  # ----------------------------------------------------------------------------
  # 2. EXPORT SINGLE HARMONIZED OUTLIER (Direct Copy from Network ARD)
  # ----------------------------------------------------------------------------
  cat(" -> [2/4] Checking network ARD for physical single harmonized raster...\n")
  
  ard_grid_dir <- file.path(ARD_NETWORK_DIR, grid_id)
  ard_files <- list.files(ard_grid_dir, pattern = "\\.tif$", full.names = TRUE)
  
  files_processed <- 0
  for (af in ard_files) {
    # Verify file is physically present on disk (not a symlink pointing to raw folders)
    if (!nzchar(Sys.readlink(af))) {
      outlier_yr <- stringr::str_extract(basename(af), "\\d{4}(?=\\.tif$)")
      
      if (!is.na(outlier_yr) && outlier_yr %in% years) {
        consensus_yrs <- setdiff(years, outlier_yr)
        ref_year_single <- max(consensus_yrs)
        
        single_out <- file.path(aoi_local_dir, paste0("harmonized_single_", outlier_yr, "_to_", ref_year_single, ".tif"))
        
        if (!file.exists(single_out)) {
          cat(sprintf("    [COPY] Confirmed physical ARD file for %s. Copying directly to local disk...\n", outlier_yr))
          file.copy(from = af, to = single_out, overwrite = TRUE)
        } else {
          cat(sprintf("    [SKIP] Single harmonized file for %s already exists locally.\n", outlier_yr))
        }
        files_processed <- files_processed + 1
      }
    }
  }
  
  if (files_processed == 0) cat("    [STABLE] No physical harmonized files found in ARD folder for this site.\n")
  
  # ----------------------------------------------------------------------------
  # 3. FIXED ANCHOR EVALUATION (Double Harmonization Generation)
  # ----------------------------------------------------------------------------
  cat(" -> [3/4] Checking and generating double-image harmonized dataset...\n")
  
  ref_year_double <- if ("2020" %in% years) "2020" else min(years)
  target_years_double <- setdiff(years, ref_year_double)
  
  if (length(target_years_double) >= 2) {
    all_double_exist <- all(sapply(target_years_double, function(tgt_yr) {
      file.exists(file.path(aoi_local_dir, paste0("harmonized_double_", tgt_yr, "_to_", ref_year_double, ".tif")))
    }))
    
    if (all_double_exist) {
      cat("    [SKIP] All double harmonized rasters already exist locally. Skipping calculation.\n")
    } else {
      ref_rast_double <- raster_list[[ref_year_double]]
      cached_double_quantiles <- cache_reference_quantiles(ref_rast_double)
      
      for (i in seq_along(target_years_double)) {
        tgt_yr <- target_years_double[i]
        double_out <- file.path(aoi_local_dir, paste0("harmonized_double_", tgt_yr, "_to_", ref_year_double, ".tif"))
        
        if (!file.exists(double_out)) {
          cat(sprintf("    [CALC] Generating double harmonization for %s -> %s...\n", tgt_yr, ref_year_double))
          harm_double <- apply_cached_harmonization(raster_list[[tgt_yr]], cached_double_quantiles)
          terra::writeRaster(harm_double, double_out, overwrite = TRUE)
          rm(harm_double)
        } else {
          cat(sprintf("    [SKIP] Double harmonized file for %s already exists locally.\n", tgt_yr))
        }
      }
      rm(cached_double_quantiles)
    }
  } else {
    cat(" -> Note: Fewer than 2 target years available; skipping double harmonization export.\n")
  }
  
  # ----------------------------------------------------------------------------
  # 4. EXPORT AOI BOUNDARY VECTOR FILES (Fast In-Memory Search)
  # ----------------------------------------------------------------------------
  cat(" -> [4/4] Searching for and copying AOI spatial boundary files...\n")
  
  grid_folders <- unique(dirname(c(unlist(file_map), ard_grid_dir)))
  
  # Immediate parent check
  aoi_files_local <- list.files(
    grid_folders,
    pattern = "\\.(shp|shx|dbf|prj|cpg|sbn|sbx|gpkg|geojson)$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  # Fast in-memory grep against pre-indexed network tree
  aoi_files_named <- grep(
    paste0("/", grid_id, ".*\\.(shp|shx|dbf|prj|cpg|sbn|sbx|gpkg|geojson)$"), 
    all_network_files, 
    value = TRUE, 
    ignore.case = TRUE
  )
  
  aoi_files <- unique(c(aoi_files_local, aoi_files_named))
  
  if (length(aoi_files) > 0) {
    for (af in aoi_files) {
      target_aoi <- file.path(aoi_local_dir, basename(af))
      if (!file.exists(target_aoi)) {
        file.copy(from = af, to = target_aoi, overwrite = TRUE)
        cat(sprintf("    [COPY] Copied boundary file: %s\n", basename(af)))
      } else {
        cat(sprintf("    [SKIP] Boundary file %s already exists locally.\n", basename(af)))
      }
    }
    cat(sprintf("    [AOI] Ensured %d spatial vector file(s) are present in local storage.\n", length(aoi_files)))
  } else {
    cat("    [AOI] Warning: No spatial vector boundary files found for this grid ID.\n")
  }
  
  rm(raster_list)
  terra::tmpFiles(remove = TRUE)
  gc(verbose = FALSE)
  
  return(invisible(TRUE))
}

# ------------------------------------------------------------------------------
# 5. Execution: Sequential Testing vs. Parallel Production
# ------------------------------------------------------------------------------
if (TEST_MODE) {
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf("[TEST MODE ENABLED] Running sequentially on %d test sites...\n", length(target_grids)))
  cat(sprintf("=======================================================\n"))
  
  for (grid_id in target_grids) {
    tryCatch({
      process_local_aoi_export(grid_id = grid_id)
    }, error = function(e) {
      cat(sprintf(" [ERROR] Failed processing AOI %s: %s\n", grid_id, e$message))
    })
  }
  
} else {
  cat(sprintf("\n=======================================================\n"))
  cat(sprintf("[PRODUCTION MODE] Starting parallel run across %d workers for %d sites...\n", NUM_CORES, length(target_grids)))
  cat(sprintf("=======================================================\n"))
  
  cl <- parallel::makeCluster(NUM_CORES)
  doParallel::registerDoParallel(cl)
  
  # Export the pre-indexed network tree to parallel workers
  parallel::clusterExport(cl, varlist = c("all_network_files", "ARD_NETWORK_DIR", "LOCAL_EXPORT_DIR", "cache_reference_quantiles", "apply_cached_harmonization", "SAMPLE_SIZE", "N_QUANTILES"))
  
  foreach(
    grid_id = target_grids,
    .packages = c("dplyr", "terra", "sf", "readr", "stringr", "tidyr", "stats"),
    .errorhandling = "pass"
  ) %dopar% {
    tryCatch({
      process_local_aoi_export(grid_id = grid_id)
    }, error = function(e) {
      cat(sprintf(" [ERROR] Failed processing AOI %s: %s\n", grid_id, e$message))
    })
  }
  
  parallel::stopCluster(cl)
}

cat(sprintf("\nRun complete. All model-ready files exported to: %s\n", LOCAL_EXPORT_DIR))