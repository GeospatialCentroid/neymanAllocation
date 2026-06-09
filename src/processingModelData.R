
# this hit a memory error with 20 cores 
process_model_data_fast <- function(aoi_paths, model_paths, local_output_dir, num_cores = 10) {
  
  # Ensure the local export directory exists
  if (!dir.exists(local_output_dir)) {
    dir.create(local_output_dir, recursive = TRUE)
  }
  
  message("Step 1: Building a high-speed lookup table for AOI paths...")
  
  # Build a clean dataframe index using a vectorized regex extraction to pull the ID
  id_regex <- "[0-9]+-[0-9]+-[a-zA-Z0-9]+-[0-9]+-[0-9]+"
  
  aoi_matches <- regexpr(id_regex, aoi_paths)
  valid_aois <- aoi_matches != -1
  
  aoi_lookup <- data.frame(
    aoi_id = regmatches(aoi_paths, aoi_matches),
    aoi_path = aoi_paths[valid_aois],
    stringsAsFactors = FALSE
  ) %>% 
    # Force a unique row-per-ID index to avoid duplication bugs on the vector side
    distinct(aoi_id, .keep_all = TRUE)
  
  message(sprintf("Lookup table complete. Found %d unique AOI features.", nrow(aoi_lookup)))
  
  # Set up the parallel processing backend
  message(sprintf("Step 2: Initializing multisession parallel workers on %d cores...", num_cores))
  plan(multisession, workers = num_cores)
  
  # 3. Process the model paths in parallel
  message("Step 3: Beginning batch processing...")
  
  results <- future_lapply(model_paths, function(model_path) {
    # We must load terra explicitly inside the isolated parallel worker environments
    library(terra)
    
    # Extract the matching ID from the current image path
    img_id_match <- regexpr("[0-9]+-[0-9]+-[a-zA-Z0-9]+-[0-9]+-[0-9]+", model_path)
    if (img_id_match == -1) {
      return(paste0("Error: Could not parse ID from ", basename(model_path)))
    }
    
    img_id <- regmatches(model_path, img_id_match)
    
    # Safe lookup selection forcing a single string pull [1], handling multi-year references smoothly
    match_path <- aoi_lookup$aoi_path[aoi_lookup$aoi_id == img_id][1]
    
    if (is.na(match_path) || length(match_path) == 0) {
      return(paste0("Skipped: No AOI matching ID ", img_id))
    }
    
    # Define local export name (Preserves the unique year string natively from the model_path)
    img_name <- basename(model_path)
    new_img_name <- sub("_[0-9\\.]+km_", "_1km_", img_name)
    local_export_path <- file.path(local_output_dir, new_img_name)
    
    # Skip if file already exists (makes the script execution resume-on-fail safe)
    if (file.exists(local_export_path)) {
      return(paste0("Skipped: Already exists ", new_img_name))
    }
    
    # Spatial Processing Block
    tryCatch({
      # Optimize terra configuration settings inside the worker thread
      terraOptions(gdal_cachemax = 512, tempdir = "/tmp") 
      
      v <- terra::vect(match_path)
      r <- terra::rast(model_path)
      
      # Ensure Coordinate Reference Systems (CRS) match
      if (terra::crs(v) != terra::crs(r)) {
        v <- terra::project(v, terra::crs(r))
      }
      
      # Spatial Crop and Mask
      r_cropped <- terra::crop(r, v)
      r_masked <- terra::mask(r_cropped, v)
      
      # Write out compressed locally to save disk I/O time
      terra::writeRaster(r_masked, local_export_path, 
                         overwrite = TRUE, 
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
      
      return(paste0("Success: ", new_img_name))
      
    }, error = function(e) {
      return(paste0("Error processing ID ", img_id, " (", img_name, "): ", e$message))
    })
  }, future.seed = TRUE)
  
  # Return back to a standard sequential execution plan
  plan(sequential)
  
  # Print execution diagnostics
  errors <- results[grepl("^Error", results)]
  successes <- results[grepl("^Success", results)]
  skipped <- results[grepl("^Skipped", results)]
  
  message("--------------------------------------------------")
  message(sprintf("Processing complete summary:\n Total Files: %d\n Successful Exports: %d\n Skipped Files: %d\n Failed Files: %d", 
                  length(model_paths), length(successes), length(skipped), length(errors)))
  message("--------------------------------------------------")
  
  if(length(errors) > 0) {
    message("First 20 recorded errors:")
    print(head(errors, 20))
  }
}
