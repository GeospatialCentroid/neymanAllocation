# ==============================================================================
# 04_random_sampling_test.R (Weighted Version)
# Purpose: Baseline Stress-Test using Simple Random Sampling (SRS) for all years.
# ==============================================================================

source("scripts/00_config.R")

# --- 1. CONFIGURATION ---------------------------------------------------------

USE_UNET <- TRUE

if (USE_UNET) {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes_unet")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results_unet")
  input_suffix <- "_UNET_master_dataset.csv"
  output_suffix <- "_baseline_weighted_results_unet.csv"
  method_name <- "BASELINE_SRS_WEIGHTED_UNET"
} else {
  INPUT_DYNAMIC_DIR <- file.path(DERIVED_DIR, "dynamic_attributes")
  OUTPUT_SIM_DIR <- file.path(DERIVED_DIR, "simulation_results")
  input_suffix <- "_master_dataset.csv"
  output_suffix <- "_baseline_weighted_results.csv"
  method_name <- "BASELINE_SRS_WEIGHTED"
}

if (!dir.exists(OUTPUT_SIM_DIR)) {
  dir.create(OUTPUT_SIM_DIR, recursive = TRUE)
}

SIM_REPS <- 100
SAMPLE_STEP <- 50
ACCURACY_TOLERANCE <- 0.10
ACCURACY_PASS_RATE_80 <- 0.80
ACCURACY_PASS_RATE_95 <- 0.95
COVERAGE_PASS_RATE_95 <- 0.95 # Added for 95% CI Coverage
AREA_COL <- "grid_area"
YEARS <- c(2010, 2016, 2020)

# --- 2. HELPER FUNCTIONS ------------------------------------------------------

analyze_weighted_sample <- function(
  sample_df,
  N_total,
  area_col = "grid_area"
) {
  n <- nrow(sample_df)
  x <- sample_df[[area_col]]
  y <- sample_df$TOF * x

  R_hat <- sum(y) / sum(x)
  fpc <- 1 - (n / N_total)
  x_bar <- mean(x)
  residuals <- y - (R_hat * x)
  s2_res <- sum(residuals^2) / (n - 1)
  se <- sqrt(fpc / n) * (1 / x_bar) * sqrt(s2_res)

  # 95% Confidence Interval
  ci_lower <- R_hat - (1.96 * se)
  ci_upper <- R_hat + (1.96 * se)

  return(c(mean = R_hat, ci_l = ci_lower, ci_u = ci_upper))
}

# --- 3. MAIN SIMULATION LOOP --------------------------------------------------

mlra_ids <- if (!is.null(TARGET_MLRA_IDS)) TARGET_MLRA_IDS else ALL_MLRA_IDS

for (m_id in mlra_ids) {
  message(paste0(
    "\n--- Starting Weighted Random Sampling Baseline: ",
    m_id,
    " ---"
  ))
  data_file <- file.path(INPUT_DYNAMIC_DIR, paste0("MLRA_", m_id, input_suffix))

  if (!file.exists(data_file)) {
    warning(paste("Data missing for MLRA", m_id))
    next
  }

  df_full <- readr::read_csv(data_file, show_col_types = FALSE)
  if (!AREA_COL %in% names(df_full)) {
    stop("Error: Area column not found.")
  }

  sim_results_list <- list()

  for (yr in YEARS) {
    pop_df <- df_full %>% dplyr::filter(year == yr, !is.na(TOF))
    if (nrow(pop_df) == 0) {
      next
    }

    TRUE_MEAN <- weighted.mean(pop_df$TOF, pop_df[[AREA_COL]], na.rm = TRUE)
    TOTAL_POP <- nrow(pop_df)
    TOLERANCE_VAL <- TRUE_MEAN * ACCURACY_TOLERANCE

    message(paste0(
      "   Year: ",
      yr,
      " | Pop: ",
      TOTAL_POP,
      " | Mean: ",
      round(TRUE_MEAN, 2),
      "%"
    ))

    max_n <- TOTAL_POP

    for (n_total in seq(100, max_n, SAMPLE_STEP)) {
      sim_outcomes <- data.frame(
        mean_est = numeric(SIM_REPS),
        covered = logical(SIM_REPS),
        accurate = logical(SIM_REPS)
      )

      for (r in 1:SIM_REPS) {
        set.seed(r)
        drawn_sample <- pop_df %>% dplyr::slice_sample(n = n_total)
        est <- analyze_weighted_sample(drawn_sample, TOTAL_POP, AREA_COL)

        sim_outcomes$mean_est[r] <- est["mean"]
        sim_outcomes$accurate[r] <- abs(est["mean"] - TRUE_MEAN) <=
          TOLERANCE_VAL
        sim_outcomes$covered[r] <- (est["ci_l"] <= TRUE_MEAN) &
          (est["ci_u"] >= TRUE_MEAN)
      }

      accuracy_rate <- mean(sim_outcomes$accurate)
      coverage_rate <- mean(sim_outcomes$covered)

      pass_80 <- accuracy_rate >= ACCURACY_PASS_RATE_80
      pass_95 <- accuracy_rate >= ACCURACY_PASS_RATE_95
      pass_cov <- coverage_rate >= COVERAGE_PASS_RATE_95
      pass_both <- pass_95 & pass_cov

      sim_results_list[[length(sim_results_list) + 1]] <- data.frame(
        MLRA = m_id,
        year = yr,
        Method_Name = method_name,
        Sample_N = n_total,
        Avg_Estimate = mean(sim_outcomes$mean_est),
        Accuracy_Rate = accuracy_rate,
        Coverage_Rate = coverage_rate,
        Passed_Accuracy_80 = pass_80,
        Passed_Accuracy_95 = pass_95,
        Passed_Coverage_95 = pass_cov,
        Passed_Both = pass_both
      )

      # Early exit: Stop sampling only once BOTH 95% accuracy and 95% CI coverage are achieved
      if (pass_both) {
        message(sprintf(
          "      -> Reached 95%% Acc & 95%% CI at N = %d. Moving to next feature.",
          n_total
        ))
        break
      }
    }
  }

  if (length(sim_results_list) > 0) {
    final_sim_df <- dplyr::bind_rows(sim_results_list)
    out_file <- file.path(OUTPUT_SIM_DIR, paste0("MLRA_", m_id, output_suffix))
    readr::write_csv(final_sim_df, out_file)
    message(paste("   Saved results to:", out_file))
  }
}

library(dplyr)
library(ggplot2)

#' Plot Milestone Sample Sizes for a Single MLRA
#'
#' @param sim_data The final output dataframe from your simulation script.
#' @param target_mlra A single MLRA ID (character or numeric) to visualize.
#' @param buffer An integer value to add to the passing sample sizes (default: 0).
#' @return A list containing the ggplot object and the summarized data.
plot_mlra_milestones <- function(
  sim_data,
  target_mlra,
  buffer = 0
) {
  # 1. Filter to the specific MLRA
  filtered_data <- sim_data %>%
    dplyr::filter(MLRA == target_mlra)

  if (nrow(filtered_data) == 0) {
    stop(paste("MLRA", target_mlra, "was not found in the dataset."))
  }

  # 2. Extract the first passing instance for each of the 3 milestones

  # Milestone 1: 80% Accuracy
  first_80 <- filtered_data %>%
    dplyr::filter(Passed_Accuracy_80 == TRUE) %>%
    dplyr::group_by(year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "1. 80% Accuracy Threshold")

  # Milestone 2: 95% Accuracy
  first_95 <- filtered_data %>%
    dplyr::filter(Passed_Accuracy_95 == TRUE) %>%
    dplyr::group_by(year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "2. 95% Accuracy Threshold")

  # Milestone 3: 95% Accuracy AND 95% Confidence Interval Coverage
  first_both <- filtered_data %>%
    dplyr::filter(Passed_Both == TRUE) %>%
    dplyr::group_by(year) %>%
    dplyr::arrange(Sample_N, .by_group = TRUE) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::mutate(Milestone = "3. 95% Acc + 95% CI Coverage")

  # Combine and format the data
  milestone_df <- dplyr::bind_rows(first_80, first_95, first_both) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Buffered_Sample_N = Sample_N + buffer) %>%
    dplyr::select(year, Milestone, Original_N = Sample_N, Buffered_Sample_N)

  # 3. Generate the Visualization
  p <- ggplot(
    milestone_df,
    aes(x = factor(year), y = Buffered_Sample_N, fill = Milestone)
  ) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      alpha = 0.85
    ) +
    # Using viridis for a clean, accessible, modern look
    scale_fill_viridis_d(
      name = "Milestone Reached",
      option = "mako",
      begin = 0.3,
      end = 0.9
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Required Sample Size Milestones - MLRA", target_mlra),
      x = "Simulation Year",
      y = if (buffer > 0) {
        paste("Buffered Sample Size (N +", buffer, ")")
      } else {
        "Sample Size (N)"
      }
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold")
    )

  return(list(plot = p, data = milestone_df))
}

# Example execution: Pick ONE MLRA to visualize
my_mlra <- "78"

# Note: Handled binding assuming files exist.
sim_results <- list.files(OUTPUT_SIM_DIR, full.names = TRUE) |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::bind_rows()

# Run the function (added a buffer of 0 by default, change to 500 if desired)
results <- plot_mlra_milestones(
  sim_data = sim_results,
  target_mlra = my_mlra,
  buffer = 0
)

# Display the plot
print(results$plot)

# View the exact stability numbers per milestone
print(results$data)




# ==============================================================================
# 4. PRESENTATION MAPPING FUNCTION (Simple Random Sampling)
# ==============================================================================
library(sf)
library(tigris) 

options(tigris_use_cache = TRUE) 

plot_random_map <- function(
    df,
    grid_sf,
    mlra_id,
    sample_size = 400,
    id_col = "id",
    gpkg_path = NULL
) {
  
  # 1. Filter CSV by Year
  if ("year" %in% names(df)) {
    target_year <- max(df$year, na.rm = TRUE)
    plot_df <- df %>% dplyr::filter(year == target_year)
  } else {
    plot_df <- df
  }
  
  # 2. Filter Spatial Grid to ONLY the cells that actually exist in the CSV
  mlra_grid <- grid_sf %>%
    dplyr::filter(.data[[id_col]] %in% plot_df[[id_col]])
  
  # 3. Determine Random Indices (THE MAIN DIFFERENCE)
  N_total <- nrow(mlra_grid)
  set.seed(42) # Keep seed fixed so the random map doesn't jitter between presentations
  
  if (sample_size >= N_total) {
    sampled_ids <- mlra_grid[[id_col]]
  } else {
    # Perform a simple random sample without replacement
    sampled_ids <- sample(mlra_grid[[id_col]], size = sample_size, replace = FALSE)
  }
  
  # 4. Assign plotting status
  mlra_grid <- mlra_grid %>%
    dplyr::mutate(Sample_Status = dplyr::case_when(
      .data[[id_col]] %in% sampled_ids ~ "Randomly Selected",
      TRUE ~ "Present in CSV Data"
    ))
  
  mlra_grid$Sample_Status <- factor(
    mlra_grid$Sample_Status,
    levels = c("Present in CSV Data", "Randomly Selected")
  )
  
  # 5. Generate the Map
  p <- ggplot()
  
  # Draw MLRA background and State Line
  if (!is.null(gpkg_path) && file.exists(gpkg_path)) {
    mlra_bound <- sf::st_read(gpkg_path, quiet = TRUE) %>%
      dplyr::filter(MLRA_ID == mlra_id)
    
    p <- p + geom_sf(data = mlra_bound, fill = "white", color = "black", linewidth = 0.8)
    
    # Fetch Nebraska State Line
    ne_state <- tigris::states(cb = TRUE, class = "sf", progress_bar = FALSE) %>%
      dplyr::filter(STUSPS == "NE") %>%
      sf::st_transform(sf::st_crs(mlra_bound)) 
    
    # Convert POLYGON to MULTILINESTRING before cropping
    ne_state_lines <- sf::st_cast(ne_state, "MULTILINESTRING")
    
    # Physically crop the lines to the MLRA's bounding box
    ne_state_cropped <- suppressWarnings(sf::st_crop(ne_state_lines, sf::st_bbox(mlra_bound)))
    
    # Add the cropped state line on top
    p <- p + geom_sf(
      data = ne_state_cropped, 
      fill = NA, 
      color = "#1f78b4", 
      linewidth = 1.2, 
      linetype = "dashed",
      alpha = 0.8
    )
  }
  
  # Draw the Grid Cells
  p <- p +
    geom_sf(
      data = mlra_grid,
      aes(fill = Sample_Status, color = Sample_Status),
      linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c("Present in CSV Data" = "gray90", "Randomly Selected" = "red")
    ) +
    scale_color_manual(
      values = c("Present in CSV Data" = "gray80", "Randomly Selected" = "darkred")
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Simple Random Sampling - MLRA", mlra_id),
      subtitle = paste("Sample Size: N =", sample_size, "| Dashed blue line indicates NE boundary."),
      fill = NULL, color = NULL
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}


# ==============================================================================
# WORKED EXAMPLE: Rendering the Random Sample Map for MLRA 86
# ==============================================================================

# Ensure required spatial libraries are loaded
library(sf)
library(tigris)
library(dplyr)
library(ggplot2)

# 1. Define target parameters
target_mlra <- "86"
n_samples <- 400

# 2. Load the Attribute Data (CSV)
# Assuming INPUT_DYNAMIC_DIR and input_suffix are set from your 00_config.R
csv_file_path <- file.path(INPUT_DYNAMIC_DIR, paste0("MLRA_", target_mlra, input_suffix))
mlra_data <- readr::read_csv(csv_file_path, show_col_types = FALSE)

# 3. Load the Spatial Grid Data (sf object)
# IMPORTANT: Replace this path with the actual file containing your spatial grid cells.
# This spatial file MUST contain a column that matches the ID column in your CSV.
grid_file_path <- "data/derived/grids/Nebraska_1km_mlra.gpkg"
mlra_grid_spatial <- sf::st_read(grid_file_path, quiet = TRUE)

# 4. Define the path to your MLRA boundaries (Optional, but needed to draw the Nebraska line)
# IMPORTANT: Replace this path with your actual MLRA boundary shapefile or geopackage.
mlra_boundary_path <- "path/to/your/mlra_boundaries.gpkg"

# 5. Execute the Plotting Function
# NOTE: Check your data to ensure "id" is the correct linking column. 
# If your column is named "grid_id" or "UID", change the id_col argument below.
random_sample_map <- plot_random_map(
  df = mlra_data,
  grid_sf = mlra_grid_spatial,
  mlra_id = target_mlra,
  sample_size = n_samples,
  id_col = "id", 
  gpkg_path = mlra_boundary_path
)

# 6. Display the Map
print(random_sample_map)

# 7. (Optional) Save the Map to your output directory
# ggsave(
#   filename = file.path(OUTPUT_SIM_DIR, paste0("MLRA_", target_mlra, "_Random_Map.png")),
#   plot = random_sample_map,
#   width = 10, 
#   height = 8, 
#   dpi = 300,
#   bg = "white"
# )

