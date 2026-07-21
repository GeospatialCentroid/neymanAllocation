# ------------------------------------------------------------------------------
# 1. Package Loading
# ------------------------------------------------------------------------------
pacman::p_load(dplyr, terra, sf, readr, furrr, future, stringr, tidyr, ggplot2, tidyterra, patchwork)

# ------------------------------------------------------------------------------
# 2. Data Ingestion & Configuration
# ------------------------------------------------------------------------------
# Load MLRA boundaries and spatial data
mlras <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg")

# Directory with the imagery  

naip_dir <- "/mnt/unraid_naip/LLR G/"
naips <- list.dirs(naip_dir, recursive = FALSE, full.names = TRUE)

# Directory with the predicted models
model_dir <- "/mnt/unraid_naip/LLR G/modelOutputs"

# Details on how to group the results 
summaryGroups <- read_csv("data/products/groundTruthSamples/scenarios_2020/all_scenarios_combined_200.csv", show_col_types = FALSE)

# Define export directory and file paths
export_dir <- "~/trueNAS/work/neymanSampling/data/products/groundTruthSamples/scenarios_2020"
individualSummary <- read_csv(file.path(export_dir, "grid_level_tof.csv"))
targetYearTOF <- read_csv(file.path(export_dir, "scenario_tof_summary_target_year.csv"))

# ------------------------------------------------------------------------------
# 3. Identify Grids with Highest, Middle, and Lowest Variance
# ------------------------------------------------------------------------------
grid_differences <- individualSummary %>%
  mutate(
    base_grid = str_extract(grid_id, "^[^_]+"),
    year = str_extract(grid_id, "\\d{4}")
  ) %>%
  group_by(base_grid) %>%
  summarize(
    min_area = min(grid_tof_area, na.rm = TRUE),
    max_area = max(grid_tof_area, na.rm = TRUE),
    area_diff = max_area - min_area,
    .groups = "drop"
  ) %>%
  arrange(desc(area_diff))

# quick visual of the distribution 
# 2. Define the threshold point
threshold_val <- 20

# 3. Generate the histogram plot
ggplot(grid_differences, aes(x = area_diff)) +
  # Create histogram (adjust binwidth based on your full dataset size)
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white", alpha = 0.8) +
  
  # Add vertical line for the 20% difference point
  geom_vline(xintercept = threshold_val, color = "firebrick", linetype = "dashed", linewidth = 1.2) +
  
  # Add a text annotation pointing out the threshold line
  annotate("text", x = threshold_val + 2, y = Inf, label = "Over 20% Difference", 
           color = "firebrick", vjust = 2, hjust = 0, fontface = "bold") +
  
  # Clean up labels and theme
  labs(
    title = "Distribution of Area Difference (area_diff)",
    subtitle = "Dashed line indicates the 20% difference threshold",
    x = "Area Difference",
    y = "Record Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

# Extract the three representative indexes
max_idx <- 10
min_idx <- nrow(grid_differences)
mid_idx <- round(nrow(grid_differences) / 2)

target_grid_max <- grid_differences$base_grid[max_idx]
target_grid_mid <- grid_differences$base_grid[mid_idx]
target_grid_min <- grid_differences$base_grid[min_idx]

cat("Selected Targets:\n")
cat("1. Highest Variance Grid [Row", max_idx, "]:", target_grid_max, "\n")
cat("2. Middle Variance Grid  [Row", mid_idx, "]:", target_grid_mid, "\n")
cat("3. Lowest Variance Grid  [Row", min_idx, "]:", target_grid_min, "\n\n")

# ------------------------------------------------------------------------------
# 4. Processing & Plot Generation Function Definition
# ------------------------------------------------------------------------------
generate_comparison_plots <- function(grid_id, naip_directory, model_directory, class_label) {
  
  # Find all NAIP files matching the target grid
  target_files <- list.files(naip_directory, pattern = paste0(grid_id, ".*\\.tif$"), full.names = TRUE, recursive = TRUE)
  naip_files <- target_files[!grepl(pattern = "modelOutputs", target_files)]
  
  if (length(naip_files) == 0) {
    warning(paste("No NAIP files found for grid:", grid_id))
    return(NULL)
  }
  
  df_list <- list()
  raster_list <- list()
  
  for (file in naip_files) {
    file_year <- str_extract(basename(file), "\\d{4}\\.tif") |>
      tools::file_path_sans_ext()
    
    r <- terra::rast(file)
    names(r) <- c("Red", "Green", "Blue", "NIR")
    raster_list[[file_year]] <- r
    
    ndvi <- (r[["NIR"]] - r[["Red"]]) / (r[["NIR"]] + r[["Red"]])
    names(ndvi) <- "NDVI"
    
    blue <- r[["Blue"]]
    names(blue) <- "Blue"
    
    pixel_data <- terra::spatSample(c(ndvi, blue), size = 10000, method = "regular", as.df = TRUE, na.rm = TRUE)
    pixel_data$Year <- as.factor(file_year) # Converted to factor for discrete color mapping in ggplot
    
    df_list[[file_year]] <- pixel_data
  }
  
  combined_pixel_data <- bind_rows(df_list)
  years_present <- sort(unique(as.character(combined_pixel_data$Year)))
  
  # ---------------------------------------------------------
  # Quantitative Differences (KS Test Output to Console)
  # ---------------------------------------------------------
  if (length(years_present) >= 2) {
    year1 <- combined_pixel_data %>% filter(Year == years_present[1])
    year2 <- combined_pixel_data %>% filter(Year == years_present[length(years_present)])
    
    ndvi_ks <- ks.test(year1$NDVI, year2$NDVI)
    blue_ks <- ks.test(year1$Blue, year2$Blue)
    
    cat("--- Statistics for", class_label, "Grid (", grid_id, ") ---\n")
    cat("NDVI KS Test D-statistic (", years_present[1], "vs", years_present[length(years_present)], "):", round(ndvi_ks$statistic, 4), "\n")
    cat("Blue KS Test D-statistic (", years_present[1], "vs", years_present[length(years_present)], "):", round(blue_ks$statistic, 4), "\n\n")
  }
  
  # ---------------------------------------------------------
  # Map Visualizations
  # ---------------------------------------------------------
  # 1. NAIP RGB Plots
  rgb_plots <- lapply(names(raster_list), function(yr) {
    ggplot() +
      geom_spatraster_rgb(data = raster_list[[yr]], r = 1, g = 2, b = 3) +
      theme_minimal() +
      labs(title = paste("NAIP", yr)) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  })
  
  # 2. Predicted Model Plots
  model_plots <- lapply(names(raster_list), function(yr) {
    model_file <- list.files(model_directory, pattern = paste0(grid_id, "_", yr, ".*\\.tif$"), full.names = TRUE)
    
    if (length(model_file) > 0) {
      mod_r <- terra::rast(model_file[1])
      ggplot() +
        geom_spatraster(data = mod_r) +
        scale_fill_viridis_c(na.value = "transparent", name = "Predicted\nArea") + 
        theme_minimal() +
        labs(title = paste("Model", yr)) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              legend.position = "bottom", legend.key.width = unit(1, "cm"))
    } else {
      ggplot() + 
        theme_void() + 
        labs(title = paste("Model", yr, "- Missing")) +
        theme(plot.title = element_text(hjust = 0.5))
    }
  })
  
  # ---------------------------------------------------------
  # Histogram Visualizations
  # ---------------------------------------------------------
  # 3. NDVI Histogram
  p_ndvi_hist <- ggplot(combined_pixel_data, aes(x = NDVI, fill = Year)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 50, color = "white", size = 0.2) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "NDVI Distribution", x = "NDVI", y = "Count") +
    theme(legend.position = "bottom")
  
  # 4. Blue Band Histogram
  p_blue_hist <- ggplot(combined_pixel_data, aes(x = Blue, fill = Year)) +
    geom_histogram(alpha = 0.5, position = "identity", bins = 50, color = "white", size = 0.2) +
    scale_fill_viridis_d() +
    theme_minimal() +
    labs(title = "Blue Band Distribution", x = "Blue Reflectance", y = "Count") +
    theme(legend.position = "bottom")
  
  # Combine histograms using patchwork
  p_histograms <- (p_ndvi_hist | p_blue_hist) +
    patchwork::plot_annotation(
      title = paste0(class_label, " Variance Site | Spectral Histograms | Grid: ", grid_id),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Combine RGB plots into a 1-row layout
  p_rgb <- patchwork::wrap_plots(rgb_plots, nrow = 1) + 
    patchwork::plot_annotation(
      title = paste0(class_label, " Variance Site | NAIP RGB | Grid: ", grid_id),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Combine Model plots into a separate 1-row layout
  p_model <- patchwork::wrap_plots(model_plots, nrow = 1) + 
    patchwork::plot_annotation(
      title = paste0(class_label, " Variance Site | Model Output | Grid: ", grid_id),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Print explicitly to the graphics device
  print(p_rgb)
  print(p_model)
  print(p_histograms)
  
  # Return as a named list so they can be accessed individually
  return(list(rgb = p_rgb, model = p_model, histograms = p_histograms))
}

# ------------------------------------------------------------------------------
# 5. Execution & Plot Management
# ------------------------------------------------------------------------------
print("Generating High Variance Plots...")
plots_max <- generate_comparison_plots(target_grid_max, naips, model_dir, "Highest")
plots_max
print("Generating Middle Variance Plots...")
# plots_mid <- generate_comparison_plots(target_grid_mid, naips, model_dir, "Middle")
 
print("Generating Low Variance Plots...")
# plots_min <- generate_comparison_plots(target_grid_min, naips, model_dir, "Lowest")
