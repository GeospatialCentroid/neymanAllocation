# ==============================================================================
# 02_dynamic_processing_unet.R
# Purpose: Attribute 1km grids with dynamic Tree Outside Forest (TOF) data
#          using yearly UNET mosaic rasters.
# Output:  A master dataset (Long Format) for Neyman Allocation.
# ==============================================================================

# Source configuration from the scripts folder
source("scripts/00_config.R")

# --- 1. SETUP & PATHS ---------------------------------------------------------

# Define where the Step 1 (Static) outputs are located
STATIC_DIR <- file.path(DERIVED_DIR, "static_attributes")

# Define where the UNET Master datasets will be saved
OUTPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")

if (!dir.exists(OUTPUT_DYNAMIC_DIR)) {
  dir.create(OUTPUT_DYNAMIC_DIR, recursive = TRUE)
}

# Define the directory containing UNET raster mosaics
UNET_DIR <- file.path(INPUT_DIR, "cot_unet")

# Define the years you want to process (matches your target NLCD years)
YEARS <- c(2010, 2016, 2020)

# --- 2. LOAD SHARED INPUTS ----------------------------------------------------
message("Loading shared spatial inputs...")

# Load Vectors using paths defined in 00_config.R
mlra_poly <- terra::vect(STATIC_INPUTS$mlra)
grid_1km <- terra::vect(STATIC_INPUTS$grid_1km)

# Note: grid_12mile is no longer needed since we are using seamless yearly mosaics.

# --- 3. MAIN PROCESSING LOOP: EXTRACTION & MERGE ------------------------------

available_ids <- unique(mlra_poly$MLRA_ID)

if (!is.null(TARGET_MLRA_IDS)) {
  mlra_ids <- intersect(available_ids, TARGET_MLRA_IDS)
  missing_ids <- setdiff(TARGET_MLRA_IDS, available_ids)
  if (length(missing_ids) > 0) {
    warning(paste(
      "Requested IDs missing:",
      paste(missing_ids, collapse = ", ")
    ))
  }
} else {
  mlra_ids <- available_ids
}

if (length(mlra_ids) == 0) {
  stop("No valid MLRA IDs found.")
}

for (m_id in mlra_ids) {
  message(paste0("\nGenerating UNET Master Dataset for MLRA: ", m_id))

  # --- A. Setup Paths ---
  static_attr_file <- file.path(
    STATIC_DIR,
    paste0("MLRA_", m_id, "_static_attributes.csv")
  )
  final_out_file <- file.path(
    OUTPUT_DYNAMIC_DIR,
    paste0("MLRA_", m_id, "_UNET_master_dataset.csv")
  )

  if (file.exists(final_out_file)) {
    message("   Master dataset exists. Skipping.")
    next
  }

  if (!file.exists(static_attr_file)) {
    warning(paste("   Static attributes file not found:", static_attr_file))
    next
  }

  # --- B. Setup MLRA Grids ---
  mlra_grids_sf <- sf::st_as_sf(grid_1km[grid_1km$MLRA_ID == m_id, ])

  # Initialize a summary dataframe with just the grid IDs
  cot_summary <- data.frame(id = mlra_grids_sf$id)

  # --- C. Extract Data Year-by-Year ---
  for (yr in YEARS) {
    raster_file <- file.path(UNET_DIR, paste0("mosaic_10_10_", yr, ".tif"))

    if (!file.exists(raster_file)) {
      warning(paste0("   Raster missing for year ", yr, ": ", raster_file))
      cot_summary[[paste0("tof_pct_", yr)]] <- NA
      next
    }

    message(paste0("   Extracting UNET TOF coverage for ", yr, "..."))
    r_unet <- terra::rast(raster_file)

    # ASSUMPTION: Raster is binary where 1 = TOF and 0 = Non-TOF.
    # Using fun = "mean" returns the area-weighted fraction of pixels equal to 1.
    # We multiply by 100 to convert the fraction to a percentage.
    extracted_frac <- exactextractr::exact_extract(
      r_unet,
      mlra_grids_sf,
      fun = "mean",
      progress = FALSE
    )

    cot_summary[[paste0("tof_pct_", yr)]] <- extracted_frac * 100
  }

  # --- D. Merge, Clean, and Pivot to Long Format ---
  message("   Merging and formatting UNET Master Dataset...")

  static_df <- readr::read_csv(static_attr_file, show_col_types = FALSE)

  master_df_long <- static_df %>%
    dplyr::left_join(cot_summary, by = "id") %>%
    dplyr::select(
      id,
      MLRA_ID,
      grid_area,
      riparian_pct = riparian_class_1,
      matches(
        "nlcd_\\d{4}_(Water|Developed|Barren|Forest|Shrubland|Herbaceous|Cultivated|Wetlands)"
      ),
      starts_with("tof_pct")
    ) %>%
    dplyr::rename_with(
      .fn = ~ gsub("tof_pct_(\\d{4})", "nlcd_\\1_TOF", .x),
      .cols = starts_with("tof_pct")
    ) %>%
    tidyr::pivot_longer(
      cols = starts_with("nlcd_"),
      names_pattern = "nlcd_(\\d{4})_(.*)",
      names_to = c("year", ".value")
    ) %>%
    dplyr::mutate(year = as.numeric(year))

  readr::write_csv(master_df_long, final_out_file)
  message(paste0("   Saved Long-Format Master: ", final_out_file))
}

message("\nDynamic UNET processing complete.")


# ==============================================================================
# --- 4. SUMMARY & VISUALIZATION -----------------------------------------------
# ==============================================================================

library(dplyr)
library(ggplot2)
library(readr)

#' Calculate, Plot, and Export TOF Summaries
#'
#' @param data_dir Directory containing the UNET master datasets.
#' @param plot_mlras Character vector of MLRA IDs to include in the plot.
#' @param save_outputs Logical; whether to save the results to disk (default: TRUE).
#' @return A list containing the summary table and the ggplot object.
summarize_and_plot_tof <- function(
  data_dir = OUTPUT_DYNAMIC_DIR,
  plot_mlras = c("78", "86", "150"),
  save_outputs = TRUE
) {
  message("\n--- Generating TOF Summary Table and Plot ---")

  # 1. Locate all dynamic datasets
  files <- list.files(
    data_dir,
    pattern = "_UNET_master_dataset\\.csv$",
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No dataset files found in the specified directory.")
  }

  # 2. Read and combine all files
  df_all <- lapply(files, readr::read_csv, show_col_types = FALSE) %>%
    dplyr::bind_rows()

  if (!"TOF" %in% names(df_all)) {
    stop("Column 'TOF' not found. Check the pivot logic in Step 3.")
  }

  # 3. Calculate Weighted Mean and Variance
  summary_table <- df_all %>%
    dplyr::filter(!is.na(TOF), !is.na(grid_area)) %>%
    dplyr::group_by(MLRA_ID, year) %>%
    dplyr::summarise(
      Grid_Count = dplyr::n(),
      Weighted_Mean_TOF = weighted.mean(TOF, w = grid_area, na.rm = TRUE),
      # Calculate weighted population variance
      Weighted_Var_TOF = sd(TOF),
      .groups = "drop"
    )

  # 4. Filter data specifically for the plot
  plot_data <- summary_table %>%
    dplyr::filter(MLRA_ID %in% plot_mlras)

  if (nrow(plot_data) == 0) {
    warning("None of the specified MLRAs for plotting were found in the data.")
    p <- NULL
  } else {
    # 5. Generate the Visualization
    p <- ggplot(
      plot_data,
      aes(x = factor(year), y = Weighted_Mean_TOF, fill = factor(MLRA_ID))
    ) +
      geom_col(
        position = position_dodge(width = 0.8),
        width = 0.7,
        color = "black",
        alpha = 0.85
      ) +
      geom_errorbar(
        aes(
          ymin = pmax(0, Weighted_Mean_TOF - sqrt(Weighted_Var_TOF)),
          ymax = Weighted_Mean_TOF + sqrt(Weighted_Var_TOF)
        ),
        position = position_dodge(width = 0.8),
        width = 0.25,
        alpha = 0.7
      ) +
      scale_fill_viridis_d(
        name = "MLRA",
        option = "mako",
        begin = 0.2,
        end = 0.8
      ) +
      theme_minimal(base_size = 14) +
      labs(
        title = "Weighted Mean TOF by Year and MLRA",
        subtitle = "Error bars represent ±1 Standard Deviation",
        x = "Year",
        y = "Weighted Mean TOF (%)"
      ) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold")
      )
  }

  # 6. Export Results to Disk
  if (save_outputs) {
    out_csv <- file.path(data_dir, "UNET_TOF_Summary_Table.csv")
    readr::write_csv(summary_table, out_csv)
    message(paste("   Saved summary table to:", out_csv))

    if (!is.null(p)) {
      out_png <- file.path(data_dir, "UNET_TOF_Summary_Plot.png")
      ggplot2::ggsave(
        filename = out_png,
        plot = p,
        width = 10, # Standard wide format
        height = 6,
        dpi = 300, # High resolution for presentations/reports
        bg = "white" # Prevents transparent backgrounds
      )
      message(paste("   Saved summary plot to:", out_png))
    }
  }

  return(list(table = summary_table, plot = p))
}

# --- Execute the Summary Function ---
# This will run automatically and save outputs to your DERIVED_DIR/dynamic_attributes_unet
tof_results <- summarize_and_plot_tof(
  data_dir = OUTPUT_DYNAMIC_DIR,
  plot_mlras = c("78", "86", "150"),
  save_outputs = TRUE
)

# Optional: Print to console during execution
print(tof_results$plot)
