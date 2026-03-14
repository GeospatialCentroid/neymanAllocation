# ==============================================================================
# 07_generate_summary_stats.R
# Purpose: Generate area-weighted summaries for Land Cover, TOF, and Riparian.
#          Includes Total Area calculations.
#          Creates visualizations comparing TOF against Top 3 Land Cover classes.
# ==============================================================================

source("scripts/00_config.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)

# --- 1. SETTINGS --------------------------------------------------------------

# TOGGLE: Set to TRUE for UNET processing, FALSE for original processing
USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "summaries_unet")
  file_pattern <- "MLRA_.*_UNET_master_dataset\\.csv"
  out_suffix <- "_UNET_summary_stats.csv"
  master_out <- "ALL_MLRA_UNET_summary_stats.csv"
  plot_prefix <- "UNET_"
} else {
  INPUT_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  OUTPUT_DIR <- file.path(DERIVED_DIR, "summaries")
  file_pattern <- "MLRA_.*_master_dataset\\.csv"
  out_suffix <- "_summary_stats.csv"
  master_out <- "ALL_MLRA_summary_stats.csv"
  plot_prefix <- ""
}

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# NLCD Columns to consider for "Top 3" ranking
NLCD_COLS <- c(
  "Water",
  "Developed",
  "Barren",
  "Forest",
  "Shrubland",
  "Herbaceous",
  "Cultivated",
  "Wetlands"
)

# --- 2. HELPER FUNCTIONS ------------------------------------------------------

#' Calculate Area-Weighted Means and Total Area
summarize_weighted <- function(df, mlra_id) {
  target_cols <- c("riparian_pct", "TOF", NLCD_COLS)
  target_cols <- intersect(target_cols, names(df))

  summary_df <- df %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(
      # 1. Calculate Total Area (Sum of grid_area)
      # Assuming grid_area is in square meters, convert to Hectares (Ha)
      Total_Area_Ha = sum(grid_area, na.rm = TRUE) * 0.0001,

      # 2. Calculate Weighted Means for percentages
      dplyr::across(
        dplyr::all_of(target_cols),
        ~ weighted.mean(.x, w = grid_area, na.rm = TRUE)
      ),
      .groups = "drop"
    ) %>%
    dplyr::mutate(MLRA_ID = mlra_id) %>%
    # Reorder: ID, Year, Area, then the percentages
    dplyr::select(MLRA_ID, year, Total_Area_Ha, dplyr::everything())

  return(summary_df)
}

#' Extract Top 3 Land Cover Classes + TOF
#' Reshapes summary data to focus on dominant landscape features.
get_top_features_long <- function(summary_df) {
  long_df <- list()

  for (i in 1:nrow(summary_df)) {
    row <- summary_df[i, ]

    # 1. Extract Top 3 NLCD
    nlcd_vals <- row %>%
      dplyr::select(dplyr::any_of(NLCD_COLS)) %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "Class",
        values_to = "Pct"
      ) %>%
      dplyr::arrange(desc(Pct)) %>%
      head(3)

    # 2. Extract TOF and Riparian
    special_vals <- row %>%
      dplyr::select(MLRA_ID, year, TOF, riparian_pct) %>%
      tidyr::pivot_longer(
        cols = c(TOF, riparian_pct),
        names_to = "Class",
        values_to = "Pct"
      )

    # 3. Combine
    combined <- dplyr::bind_rows(
      nlcd_vals %>% dplyr::mutate(Type = "Land Cover (Top 3)"),
      special_vals %>% dplyr::mutate(Type = "Target Feature")
    ) %>%
      dplyr::mutate(MLRA_ID = row$MLRA_ID, year = row$year)

    long_df[[i]] <- combined
  }
  return(dplyr::bind_rows(long_df))
}

#' Plot TOF vs Top 3 Land Cover
#' @param full_summary_df The combined summary dataframe of all MLRAs
#' @param target_ids Vector of MLRA IDs to display.
#' @param target_years Vector of Years to display (e.g., c(2020)). Defaults to NULL (All Years).
#' @param use_unet Boolean to update plot title.
plot_mlra_comparison <- function(
  full_summary_df,
  target_ids,
  target_years = NULL,
  use_unet = FALSE
) {
  # 1. Filter by ID
  plot_data <- full_summary_df %>%
    dplyr::filter(MLRA_ID %in% target_ids)

  # 2. Filter by Year (if argument provided)
  if (!is.null(target_years)) {
    plot_data <- plot_data %>% dplyr::filter(year %in% target_years)
  }

  if (nrow(plot_data) == 0) {
    message("No data found for the requested parameters.")
    return(NULL)
  }

  # 3. Prepare Top 3 Format
  top_features <- get_top_features_long(plot_data)

  plot_title <- ifelse(
    use_unet,
    "MLRA Characterization (UNET): TOF vs. Dominant Land Cover",
    "MLRA Characterization: TOF vs. Dominant Land Cover"
  )

  # 4. Create Plot
  p <- ggplot(top_features, aes(x = factor(year), y = Pct, fill = Class)) +
    geom_col(position = "dodge", color = "black", alpha = 0.9) +

    # SCALING: 'scales = "fixed"' ensures all Y-axes match exactly for comparison
    facet_wrap(~MLRA_ID, scales = "fixed", labeller = label_both) +

    scale_fill_viridis_d(option = "turbo", name = "Feature") +

    labs(
      title = plot_title,
      subtitle = "Showing weighted average for TOF, Riparian, and the Top 3 NLCD classes per area.",
      y = "Area-Weighted Average (%)",
      x = "Year"
    ) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )

  return(p)
}


# --- 3. MAIN EXECUTION --------------------------------------------------------

files <- list.files(INPUT_DIR, pattern = file_pattern, full.names = TRUE)

if (length(files) > 0) {
  message(paste(
    "Found",
    length(files),
    "MLRA datasets in",
    INPUT_DIR,
    "\nSummarizing..."
  ))

  all_summaries <- list()

  for (f in files) {
    f_name <- basename(f)
    m_id <- stringr::str_extract(f_name, "(?<=MLRA_)\\d+")

    # Read and Summarize
    df <- readr::read_csv(f, show_col_types = FALSE)
    summ <- summarize_weighted(df, m_id)

    # Save Individual
    out_name <- file.path(OUTPUT_DIR, paste0("MLRA_", m_id, out_suffix))
    readr::write_csv(summ, out_name)

    all_summaries[[length(all_summaries) + 1]] <- summ
  }

  # Save Master Table
  final_table <- dplyr::bind_rows(all_summaries)
  out_master_path <- file.path(OUTPUT_DIR, master_out)
  readr::write_csv(final_table, out_master_path)
  message(paste(
    "Summaries generated with Total Area (Ha). Saved to:",
    out_master_path
  ))

  # --- 4. PLOTTING EXAMPLES ---------------------------------------------------

  # Define subset from config, or use all
  target_plot_ids <- if (!is.null(TARGET_MLRA_IDS)) {
    TARGET_MLRA_IDS
  } else {
    unique(final_table$MLRA_ID)
  }

  # Example A: Plot All Years
  message("Generating plot for All Years...")
  p_all <- plot_mlra_comparison(
    final_table,
    target_plot_ids,
    target_years = NULL,
    use_unet = USE_UNET
  )
  if (!is.null(p_all)) {
    ggsave(
      file.path(OUTPUT_DIR, paste0(plot_prefix, "MLRA_Features_AllYears.png")),
      p_all,
      width = 12,
      height = 8
    )
  }

  # Example B: Plot Only 2020
  message("Generating plot for 2020 only...")
  p_2020 <- plot_mlra_comparison(
    final_table,
    target_plot_ids,
    target_years = c(2020),
    use_unet = USE_UNET
  )
  if (!is.null(p_2020)) {
    ggsave(
      file.path(OUTPUT_DIR, paste0(plot_prefix, "MLRA_Features_2020.png")),
      p_2020,
      width = 12,
      height = 8
    )
  }
} else {
  warning(paste("No master dataset files found in:", INPUT_DIR))
}
