
# this hit a memory error with 20 cores 
process_model_data_fast <- function(aoi_paths, model_paths, local_output_dir, num_cores = 10) {
  
  # Ensure the local export directory exists
  if (!dir.exists(local_output_dir)) {
    dir.create(local_output_dir, recursive = TRUE)
  }
  
  message("Step 1: Building a high-speed lookup table for AOI paths...")
  
  # The flexible regex pattern
  id_regex <- "[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+"
  
  aoi_matches <- regexpr(id_regex, aoi_paths)
  valid_aois <- aoi_matches != -1
  
  aoi_lookup <- data.frame(
    aoi_id = regmatches(aoi_paths, aoi_matches),
    aoi_path = aoi_paths[valid_aois],
    stringsAsFactors = FALSE
  ) %>% 
    distinct(aoi_id, .keep_all = TRUE)
  
  message(sprintf("Lookup table complete. Found %d unique AOI features.", nrow(aoi_lookup)))
  
  # Set up the parallel processing backend
  message(sprintf("Step 2: Initializing multisession parallel workers on %d cores...", num_cores))
  plan(multisession, workers = num_cores)
  
  # 3. Process the model paths in parallel
  message("Step 3: Beginning batch processing...")
  
  results <- future_lapply(model_paths, function(model_path) {
    library(terra)
    
    # --- FIXED HERE: Updated to use the identical flexible regex pattern ---
    img_id_match <- regexpr("[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+", model_path)
    if (img_id_match == -1) {
      return(paste0("Error: Could not parse ID from ", basename(model_path)))
    }
    
    img_id <- regmatches(model_path, img_id_match)
    
    # Safe lookup selection
    match_path <- aoi_lookup$aoi_path[aoi_lookup$aoi_id == img_id][1]
    
    if (is.na(match_path) || length(match_path) == 0) {
      return(paste0("Skipped: No AOI matching ID ", img_id))
    }
    
    # Define local export name (Preserves the unique year string)
    img_name <- basename(model_path)
    new_img_name <- sub("_[0-9\\.]+km_", "_1km_", img_name)
    local_export_path <- file.path(local_output_dir, new_img_name)
    
    # Skip if file already exists
    if (file.exists(local_export_path)) {
      return(paste0("Skipped: Already exists ", new_img_name))
    }
    
    # Spatial Processing Block
    tryCatch({
      # Direct temp files away from RAM-disk
      terraOptions(gdal_cachemax = 512, tempdir = "data/temp_terra") 
      
      v <- terra::vect(match_path)
      r <- terra::rast(model_path)
      
      if (terra::crs(v) != terra::crs(r)) {
        v <- terra::project(v, terra::crs(r))
      }
      
      r_cropped <- terra::crop(r, v)
      r_masked <- terra::mask(r_cropped, v)
      
      terra::writeRaster(r_masked, local_export_path, 
                         overwrite = TRUE, 
                         gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
      
      # Worker Memory Cleanup
      rm(v, r, r_cropped, r_masked)
      gc()
      terra::tmpFiles(remove = TRUE)
      
      return(paste0("Success: ", new_img_name))
      
    }, error = function(e) {
      gc()
      terra::tmpFiles(remove = TRUE)
      return(paste0("Error processing ID ", img_id, " (", img_name, "): ", e$message))
    })
  }, future.seed = TRUE)
  
  plan(sequential)
  
  # Diagnostics
  errors <- results[grepl("^Error", results)]
  successes <- results[grepl("^Success", results)]
  skipped <- results[grepl("^Skipped", results)]
  
  message("--------------------------------------------------")
  message(sprintf("Processing complete summary:\n Total Files: %d\n Successful Exports: %d\n Skipped Files: %d\n Failed Files: %d", 
                  length(model_paths), length(successes), length(skipped), length(errors)))
  message("--------------------------------------------------")
}

summarize_mlra_data <- function(target_mlra_id, index_df, mlra_sf, local_product_dir) {
  
  message(sprintf("Step 1: Filtering AOI IDs for MLRA: %s", target_mlra_id))
  
  # 1. Get subset of the index table for this specific MLRA
  mlra_index <- index_df %>% filter(MLRA_ID == target_mlra_id)
  
  if (nrow(mlra_index) == 0) {
    stop(sprintf("No AOI IDs found in the index table for MLRA_ID: %s", target_mlra_id))
  }
  
  # 2. Extract the specific MLRA spatial feature
  mlra_geom <- mlra_sf %>% filter(MLRA_ID == target_mlra_id)
  if (nrow(mlra_geom) == 0) stop(sprintf("MLRA_ID %s not found in the spatial MLRA dataset.", target_mlra_id))
  v_mlra <- terra::vect(mlra_geom)
  
  all_local_rasters <- list.files(local_product_dir, pattern = "\\.tif$", full.names = TRUE)
  
  message("Step 2: Processing matching raster data products via Vector intersection...")
  summary_list <- list()
  
  # Iterate through the dataframe rows to access BOTH the site ID and its .gpkg path
  for (i in 1:nrow(mlra_index)) {
    
    site_id <- mlra_index$id[i]
    
    # We pull the path to the original unmasked vector polygon
    aoi_path <- mlra_index$aoi_path[i] 
    
    matching_rasters <- all_local_rasters[grepl(site_id, all_local_rasters, fixed = TRUE)]
    if (length(matching_rasters) == 0) next
    
    tryCatch({
      # 3. Load the pure Vector polygon
      v_aoi <- terra::vect(aoi_path)
      
      if (terra::crs(v_aoi) != terra::crs(v_mlra)) {
        v_aoi <- terra::project(v_aoi, terra::crs(v_mlra))
      }
      
      # Intersect the 1km AOI with the MLRA boundary to trim any edges hanging outside the MLRA
      v_aoi_cropped <- terra::intersect(v_aoi, v_mlra)
      
      # 4. Calculate the true, mathematically perfect geometric area of the polygon
      total_geom_area_m2 <- sum(terra::expanse(v_aoi_cropped, unit="m"))
      
      # Process each raster year for this site
      for (raster_path in matching_rasters) {
        filename <- basename(raster_path)
        
        # --- FIXED REGEX (Reverted to _pred) ---
        year_match_idx <- regexpr("[0-9]{4}_pred", filename)
        year_val <- ifelse(year_match_idx != -1, 
                           sub("_pred", "", regmatches(filename, year_match_idx)), 
                           "Unknown")
        
        r <- terra::rast(raster_path)
        
        if (terra::crs(v_aoi_cropped) != terra::crs(r)) {
          v_aoi_proj <- terra::project(v_aoi_cropped, terra::crs(r))
        } else {
          v_aoi_proj <- v_aoi_cropped
        }
        
        # 5. Crop AND Mask the raster using the exact intersected vector shape
        r_cropped <- terra::crop(r, v_aoi_proj)
        r_masked <- terra::mask(r_cropped, v_aoi_proj)
        
        # Get physical size of a single pixel in this raster
        cell_size_m2 <- prod(terra::res(r_masked))
        if (terra::is.lonlat(r_masked)) {
          cell_size_m2 <- prod(terra::res(r_masked) * 111000) 
        }
        
        # 6. Tabulate ONLY the TOF features (Value 1)
        freq_table <- as.data.frame(terra::freq(r_masked))
        count_1 <- sum(freq_table$count[freq_table$value == 1], na.rm = TRUE)
        area_1_m2 <- count_1 * cell_size_m2
        
        # 7. Deduce Area 0 (Non-TOF / Masked Forest / Background)
        area_0_m2 <- total_geom_area_m2 - area_1_m2
        
        summary_list[[length(summary_list) + 1]] <- data.frame(
          MLRA_ID = target_mlra_id,
          AOI_Site_ID = site_id,
          Year = year_val,
          Total_Area_M2 = total_geom_area_m2,
          Area_Value_1_M2 = area_1_m2,
          Area_Value_0_M2 = area_0_m2,
          stringsAsFactors = FALSE
        )
        
        # --- MEMORY CLEANUP (Inner Loop) ---
        rm(r, r_cropped, r_masked)
      }
      
      # --- MEMORY CLEANUP (Outer Loop) ---
      rm(v_aoi, v_aoi_cropped)
      gc()
      terra::tmpFiles(remove = TRUE)
      
    }, error = function(e) {
      warning(sprintf("Error processing site %s: %s", site_id, e$message))
      
      # Ensure memory is released even if the site fails
      gc()
      terra::tmpFiles(remove = TRUE)
    })
  }
  
  if (length(summary_list) == 0) {
    return(data.frame())
  }
  
  return(do.call(rbind, summary_list))
}



batch_summarize_mlras <- function(mlra_ids, index_df, mlra_sf, local_product_dir) {
  
  message(sprintf("Starting batch processing for %d MLRA zones...", length(mlra_ids)))
  
  # map_dfr automatically binds the resulting data frames row-wise into a single master data frame
  master_summary_df <- purrr::map_dfr(mlra_ids, function(current_id) {
    
    message(sprintf("\n=== Processing MLRA_ID: %s ===", current_id))
    
    # Call your existing function for the individual MLRA
    # tryCatch ensures that if one MLRA fails (e.g., no rasters found), the batch continues
    tryCatch({
      mlra_res <- summarize_mlra_data(
        target_mlra_id = current_id,
        index_df = index_df,
        mlra_sf = mlra_sf,
        local_product_dir = local_product_dir
      )
      return(mlra_res)
      
    }, error = function(e) {
      warning(sprintf("Skipping MLRA_ID %s due to execution error: %s", current_id, e$message))
      return(data.frame()) # Return empty data frame so map_dfr skips it cleanly
    })
  })
  
  message("\nBatch processing complete.")
  return(master_summary_df)
}



batch_apply_nlcd_mask <- function(input_dir, nlcd_dir, output_dir, num_cores = 10) {
  
  # Ensure the local export directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ensure the temp directory exists to avoid RAM-disk crashes
  temp_dir_path <- "data/temp_terra"
  if (!dir.exists(temp_dir_path)) {
    dir.create(temp_dir_path, recursive = TRUE)
  }
  
  # 1. Dynamically list all target files in the input directory
  model_paths <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
  
  if (length(model_paths) == 0) {
    stop("No .tif files found in the specified input directory.")
  }
  
  message(sprintf("Found %d files to process. Initializing %d cores...", length(model_paths), num_cores))
  
  # Define NLCD paths based on your previous screenshot
  nlcd_2010_path <- file.path(nlcd_dir, "Annual_NLCD_LndCov_2010_CU_C1V1.tif")
  nlcd_2016_path <- file.path(nlcd_dir, "Annual_NLCD_LndCov_2016_CU_C1V1.tif")
  nlcd_2020_path <- file.path(nlcd_dir, "Annual_NLCD_LndCov_2020_CU_C1V1.tif")
  
  # 2. Set up parallel backend
  plan(multisession, workers = num_cores)
  
  # 3. Process files
  results <- future_lapply(model_paths, function(model_path) {
    # Must explicitly load terra in the isolated worker environment
    library(terra)
    
    filename <- basename(model_path)
    new_filename <- sub("\\.tif$", "_masked.tif", filename)
    export_path <- file.path(output_dir, new_filename)
    
    # Idempotency (Resume-on-fail safe)
    if (file.exists(export_path)) {
      return(paste0("Skipped (Already exists): ", new_filename))
    }
    
    # Temporal mapping logic to select the correct NLCD year
    year_match_idx <- regexpr("[0-9]{4}_pred", filename)
    if (year_match_idx == -1) return(paste0("Error: Missing year tag in ", filename))
    
    model_year <- sub("_pred", "", regmatches(filename, year_match_idx))
    
    target_nlcd_path <- switch(
      model_year,
      "2010" = nlcd_2010_path, "2011" = nlcd_2010_path, "2012" = nlcd_2010_path, "2013" = nlcd_2010_path,
      "2014" = nlcd_2016_path, "2015" = nlcd_2016_path, "2016" = nlcd_2016_path, "2017" = nlcd_2016_path,
      "2018" = nlcd_2020_path, "2019" = nlcd_2020_path, "2020" = nlcd_2020_path, "2021" = nlcd_2020_path,
      NULL
    )
    
    if (is.null(target_nlcd_path)) return(paste0("Skipped: No NLCD mapping for year ", model_year))
    
    tryCatch({
      # Reroute temp files to physical disk to prevent Out-Of-Memory host crashes
      terraOptions(gdal_cachemax = 512, tempdir = temp_dir_path) 
      
      # Load the pointers
      r_model <- rast(model_path)
      target_nlcd <- rast(target_nlcd_path)
      
      # Extent mapping and fast crop
      model_extent_poly <- as.polygons(ext(r_model), crs = crs(r_model))
      extent_proj <- project(model_extent_poly, crs(target_nlcd))
      nlcd_snippet <- crop(target_nlcd, extent_proj)
      
      # Evaluate Presence
      snippet_freq <- freq(nlcd_snippet)
      has_forest <- any(snippet_freq$value %in% c(41, 42, 43))
      
      if (!has_forest) {
        # Export a raster of zeros so downstream data architecture remains intact
        r_zero <- classify(r_model, cbind(-Inf, Inf, 0))
        writeRaster(r_zero, export_path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
        return(paste0("Processed (No Forest - Zeroed): ", new_filename))
      }
      
      # Execute Masking Math
      nlcd_aligned <- project(nlcd_snippet, r_model, method = "near")
      
      m <- c(0, 40, NA, 
             40, 43, 1, 
             43, 255, NA)
      rcl_matrix <- matrix(m, ncol = 3, byrow = TRUE)
      forest_mask <- classify(nlcd_aligned, rcl_matrix)
      
      r_masked <- mask(r_model, forest_mask)
      writeRaster(r_masked, export_path, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
      
      return(paste0("Processed (Masked): ", new_filename))
      
    }, error = function(e) {
      return(paste0("Error processing ", filename, ": ", e$message))
    })
  }, future.seed = TRUE)
  
  plan(sequential)
  
  # Diagnostics
  errors <- results[grepl("^Error", results)]
  successes <- results[grepl("^Processed", results)]
  skipped <- results[grepl("^Skipped", results)]
  
  message("--------------------------------------------------")
  message(sprintf("Masking Batch Complete:\n Total Files: %d\n Successfully Exported: %d\n Skipped: %d\n Errors: %d", 
                  length(model_paths), length(successes), length(skipped), length(errors)))
  message("--------------------------------------------------")
  
  if(length(errors) > 0) {
    message("First 20 recorded errors:")
    print(head(errors, 20))
  }
}
