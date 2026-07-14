# ==============================================================================
# 35_evaluatingHarmonizationOnModelOutputs.R
# Purpose: Compare raw vs. histogram-normalized (ARD) model outputs
# ==============================================================================

# 1. Package Loading & Environment Setup
# ==============================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(terra, sf, dplyr, readr, stringr, ggplot2, tidyterra, patchwork)

# ==============================================================================
# 2. Configuration & Paths
# ==============================================================================
naip_dir  <- "/mnt/unraid_naip/LLR G/"
model_dir <- "/mnt/unraid_naip/LLR G/modelOutputs"
ard_dir   <- "/mnt/unraid_naip/LLR G/ARD"

# Target features for normalization comparison
target_features <- tibble::tibble(
  aoi_id = c("1610-3-a-14-3", "1676-1-6-a-1", "1477-4-b-e-2", "1610-4-1-8-4", "1740-3-2-d-1"),
  baseline_year  = c("2016", "2016", "2015", "2016", "2011"),
  corrected_year = c("2012", "2012", "2012", "2012", "2015")
)

# ==============================================================================
# 3. Multi-Plot Generation Function
# ==============================================================================
generate_dual_comparison_plots <- function(grid_id, base_yr, corr_yr, naip_directory, model_directory, ard_directory) {
  
  # ----------------------------------------------------------------------------
  # A. Locate Files
  # ----------------------------------------------------------------------------
  
  # 1. Raw NAIP for corrected year
  naip_unnorm_file <- list.files(
    path = file.path(naip_directory, grid_id),
    pattern = paste0("naip_.*", grid_id, "_", corr_yr, "\\.tif$"),
    full.names = TRUE,
    ignore.case = TRUE
  )[1]

  # 2. ARD Normalized NAIP for corrected year
  naip_norm_file <- list.files(
    path = file.path(ard_directory, grid_id),
    pattern = paste0(grid_id, "_", corr_yr, "\\.tif$"),
    full.names = TRUE,
    ignore.case = TRUE
  )[1]
  
  # 3. Model Output: Baseline Year
  mod_base_file <- list.files(
    path = model_directory, 
    pattern = paste0(grid_id, "_", base_yr, "_pred\\.tif$"), 
    full.names = TRUE, 
    ignore.case = TRUE
  )[1]
  
  # 4. Model Output: Corrected Year - Unnormalized vs. Normalized
  mod_corr_unnorm_file <- list.files(
    path = model_directory, 
    pattern = paste0(grid_id, "_", corr_yr, "_pred_unnormalized\\.tif$"), 
    full.names = TRUE, 
    ignore.case = TRUE
  )[1]
  
  mod_corr_norm_file <- list.files(
    path = model_directory, 
    pattern = paste0(grid_id, "_", corr_yr, "_pred_normalized\\.tif$"), 
    full.names = TRUE, 
    ignore.case = TRUE
  )[1]
  
  # --- SAFETY CHECK ---
  required_files <- list(
    "NAIP Unnorm (2012)"  = naip_unnorm_file,
    "NAIP Norm (ARD)"     = naip_norm_file,
    "Model Base (2016)"   = mod_base_file,
    "Model Unnorm (2012)" = mod_corr_unnorm_file,
    "Model Norm (2012)"   = mod_corr_norm_file
  )
  
  for (name in names(required_files)) {
    if (is.na(required_files[[name]]) || length(required_files[[name]]) == 0) {
      stop(paste("CRITICAL ERROR: File missing for [", name, "] in AOI:", grid_id))
    }
  }
  
  # ----------------------------------------------------------------------------
  # B. Load Rasters
  # ----------------------------------------------------------------------------
  r_naip_unnorm <- terra::rast(naip_unnorm_file)
  r_naip_norm   <- terra::rast(naip_norm_file)
  
  if (terra::nlyr(r_naip_unnorm) >= 3) names(r_naip_unnorm)[1:3] <- c("Red", "Green", "Blue")
  if (terra::nlyr(r_naip_norm) >= 3)   names(r_naip_norm)[1:3]   <- c("Red", "Green", "Blue")
  
  r_mod_base        <- terra::rast(mod_base_file)
  r_mod_corr_unnorm <- terra::rast(mod_corr_unnorm_file)
  r_mod_corr_norm   <- terra::rast(mod_corr_norm_file)
  
  # ----------------------------------------------------------------------------
  # C. Spatial Alignment
  # ----------------------------------------------------------------------------
  # Load vector for specific AOI
  active_aoi_path <- file.path(naip_directory, grid_id, paste0("aoi-", grid_id, ".gpkg"))
  
  if (!file.exists(active_aoi_path)) {
    stop(paste("CRITICAL ERROR: Vector file not found for AOI ID:", grid_id, "at", active_aoi_path))
  }
  
  active_aoi <- terra::vect(active_aoi_path)
  
  if (nrow(active_aoi) == 0) {
    stop(paste("CRITICAL ERROR: No boundary found in vector file for AOI ID:", grid_id))
  }
  
  # 1. Align Model Layers
  aoi_vect_model <- terra::project(active_aoi, terra::crs(r_mod_base))
  
  r_mod_base        <- terra::mask(terra::crop(r_mod_base, aoi_vect_model), aoi_vect_model)
  r_mod_corr_unnorm <- terra::mask(terra::crop(r_mod_corr_unnorm, aoi_vect_model), aoi_vect_model)
  r_mod_corr_norm   <- terra::mask(terra::crop(r_mod_corr_norm, aoi_vect_model), aoi_vect_model)
  
  # 2. Align NAIP Layers (Corrected Year)
  aoi_vect_naip_unnorm <- terra::project(active_aoi, terra::crs(r_naip_unnorm))
  r_naip_unnorm_masked <- terra::mask(terra::crop(r_naip_unnorm, aoi_vect_naip_unnorm), aoi_vect_naip_unnorm)
  r_naip_unnorm_aligned <- terra::project(r_naip_unnorm_masked, r_mod_base, method = "bilinear")
  
  aoi_vect_naip_norm <- terra::project(active_aoi, terra::crs(r_naip_norm))
  r_naip_norm_masked <- terra::mask(terra::crop(r_naip_norm, aoi_vect_naip_norm), aoi_vect_naip_norm)
  r_naip_norm_aligned <- terra::project(r_naip_norm_masked, r_mod_base, method = "bilinear")
  
  # ----------------------------------------------------------------------------
  # D. Categorical Styling Setup
  # ----------------------------------------------------------------------------
  binary_palette <- c("Non-Canopy" = "#EAEAEA", "Canopy" = "#2E6F40")
  
  # Prepare model layers for plotting
  model_layers <- list(
    base    = r_mod_base,
    unnorm  = r_mod_corr_unnorm,
    norm    = r_mod_corr_norm
  )
  
  for (name in names(model_layers)) {
    levels(model_layers[[name]]) <- data.frame(ID = c(0, 1), Class = c("Non-Canopy", "Canopy"))
  }
  
  # ----------------------------------------------------------------------------
  # E. Figure 1: NAIP RGB Comparison
  # ----------------------------------------------------------------------------
  p_naip_raw <- ggplot() +
    geom_spatraster_rgb(data = r_naip_unnorm_aligned) +
    theme_minimal() +
    labs(title = paste0("Unnormalized NAIP (", corr_yr, ")")) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  
  p_naip_ard <- ggplot() +
    geom_spatraster_rgb(data = r_naip_norm_aligned) +
    theme_minimal() +
    labs(title = paste0("ARD Normalized NAIP (", corr_yr, ")")) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  
  fig_naip <- (p_naip_raw | p_naip_ard) +
    patchwork::plot_annotation(
      title = paste0("Radiometric Harmonization: NAIP RGB Imagery | AOI: ", grid_id),
      subtitle = paste0("Direct spatial comparison of raw versus histogram-normalized reflectance for ", corr_yr),
      theme = theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 11, color = "grey30"))
    )
  
  # ----------------------------------------------------------------------------
  # F. Figure 2: Classified Map Comparison
  # ----------------------------------------------------------------------------
  plot_model_layer <- function(raster_obj, plot_title) {
    ggplot() +
      geom_spatraster(data = raster_obj, aes(fill = Class)) +
      scale_fill_manual(values = binary_palette, na.value = "transparent", name = "Classification") +
      theme_minimal() +
      labs(title = plot_title) +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  }
  
  p_m1 <- plot_model_layer(model_layers$base,   paste0("Baseline (", base_yr, ")"))
  p_m2 <- plot_model_layer(model_layers$unnorm, paste0("Unnormalized (", corr_yr, ")"))
  p_m3 <- plot_model_layer(model_layers$norm,   paste0("Normalized (", corr_yr, ")"))
  
  # Assemble 
  fig_models <- (p_m1 | p_m2 | p_m3) +
    patchwork::plot_layout(nrow = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0("Multi-Temporal Canopy Classification Comparison | AOI: ", grid_id),
      subtitle = "Evaluating structural mask consistency across unnormalized and ARD normalized workflows",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "grey30"),
        legend.position = "right"
      )
    )
  
  # Render both figures
  print(fig_naip)
  print(fig_models)
  
  return(list(naip_comparison = fig_naip, model_comparison = fig_models))
}

# ==============================================================================
# 4. Execution Loop
# ==============================================================================
plot_outputs <- list()

for (i in 1:nrow(target_features)) {
  current_id   <- target_features$aoi_id[i]
  base_year    <- target_features$baseline_year[i]
  corr_year    <- target_features$corrected_year[i]
  
  cat("------------------------------------------------------------------------\n")
  cat("Generating comparison figures for AOI:", current_id, "\n")
  cat("Baseline Year:", base_year, "| Corrected Year:", corr_year, "\n")
  cat("------------------------------------------------------------------------\n")
  
  plot_outputs[[current_id]] <- generate_dual_comparison_plots(
    grid_id          = current_id,
    base_yr          = base_year,
    corr_yr          = corr_year,
    naip_directory   = naip_dir,
    model_directory  = model_dir,
    ard_directory    = ard_dir
  )
}
