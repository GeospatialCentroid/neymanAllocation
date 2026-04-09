# ==============================================================================
# Establish_stratifiedGrid.R
# Purpose: Develop the standardized stratified sample of LLR using a specified
#          number of features per MLRA.
# Output:  A comprehensive dataframe containing the 1km cell ID, MLRA ID,
#          and assigned LLR ID for all sampled grids.
# ==============================================================================

source("scripts/00_config.R")
source("src/sampleGridsFunctions.R")
source("src/systematicSampleFunctions.R")


# --- 1. SETUP & LOCAL PATHS ---------------------------------------------------
# Define sampling parameters
DRAW_SIZE <- 1400
TARGET_LLR <- "F"

message("Loading spatial inputs...")

# Load Vector Data (Assuming paths are built off a base dir like in 01_static)
# If these aren't in 00_config.R, you can define them relative to the project root here.
lrr_id_path <- "data/derived/mlra/lower48MLRA.gpkg"
mlra_grid_path <- "data/derived/grids/GreatPlains_1km_mlra.gpkg"

# quiet = TRUE suppresses the noisy sf load text to keep the console clean
lrrID <- sf::st_read(lrr_id_path, quiet = TRUE)
mlras <- sf::st_read(mlra_grid_path, quiet = TRUE)

message(paste("Subsetting MLRAs for LLR:", TARGET_LLR))

# Subset MLRA IDs based on target LLR
llr_target_ids <- lrrID %>%
  dplyr::filter(LRRSYM == TARGET_LLR) %>%
  dplyr::pull(MLRA_ID)

mlras_target <- mlras %>%
  dplyr::filter(MLRA_ID %in% llr_target_ids)


# --- 2. HELPER FUNCTIONS ------------------------------------------------------

#' Draw Systematic Sample per MLRA and format as Dataframe
get_mlra_sample_df <- function(
  spatial_data,
  target_mlra,
  n_desired,
  llr_val,
  seed_val = 1234
) {
  # Subset spatial data to the current MLRA
  mlra_subset <- spatial_data[spatial_data$MLRA_ID == target_mlra, ]

  # Skip if no data exists for this MLRA ID
  if (nrow(mlra_subset) == 0) {
    message(paste("   No grids found for MLRA", target_mlra, "- Skipping."))
    return(data.frame(
      id = character(0),
      MLRA_ID = character(0),
      LLR_ID = character(0)
    ))
  }

  message(paste(
    "   Drawing",
    n_desired,
    "samples for MLRA",
    target_mlra,
    "..."
  ))

  # Set seed for reproducible sampling
  set.seed(seed_val)

  # Draw the sample
  sampled_data <- draw_systematic_sample(
    df = as.data.frame(mlra_subset),
    n_desired = n_desired
  )

  # Construct and return the dataframe
  result_df <- data.frame(
    id = sampled_data$id,
    MLRA_ID = target_mlra,
    LLR_ID = llr_val,
    stringsAsFactors = FALSE
  )

  return(result_df)
}


# --- 3. MAIN PROCESSING LOOP --------------------------------------------------

if (length(llr_target_ids) == 0) {
  stop("No valid MLRA IDs found for the target LLR.")
}

message(paste(
  "\nFound",
  length(llr_target_ids),
  "MLRAs to process for sampling."
))

# Generate the samples across all target MLRAs using lapply
# (This mirrors the functional approach while keeping console output flowing)
sample_list <- lapply(llr_target_ids, function(m_id) {
  get_mlra_sample_df(
    spatial_data = mlras_target,
    target_mlra = m_id,
    n_desired = DRAW_SIZE,
    llr_val = TARGET_LLR
  )
})

# Bind the list of individual MLRA dataframes into one master dataframe
final_sample_df <- dplyr::bind_rows(sample_list)

message("\nSampling process complete.")
message(paste("Total grids sampled:", nrow(final_sample_df)))

# export results
readr::write_csv(final_sample_df, "data/products/systematicSampleSelection/")


# --- 4. VISUALIZATION QA/QC ---------------------------------------------------

# Set an MLRA ID here to test and visualize the output
test_mlra_id <- 62

if (test_mlra_id %in% llr_target_ids) {
  message(paste("\nGenerating QA map for MLRA", test_mlra_id, "..."))

  # Isolate the spatial data and the tabular sample data for the test MLRA
  test_mlra_spatial <- mlras_target[mlras_target$MLRA_ID == test_mlra_id, ]
  test_mlra_samples <- final_sample_df[
    final_sample_df$MLRA_ID == test_mlra_id,
  ]

  # Tag the spatial data with TRUE/FALSE if it was selected
  test_mlra_spatial$inSample <- test_mlra_spatial$id %in% test_mlra_samples$id

  tmap::tmap_mode("view")

  qa_map <- tmap::tm_shape(test_mlra_spatial) +
    tmap::tm_polygons(
      col = "inSample",
      palette = c("FALSE" = "#cccccc", "TRUE" = "#e31a1c"),
      fill_alpha = 0.7,
      border.col = "white",
      col_alpha = 0.3,
      title = paste("Sampled 1km Grids - MLRA", test_mlra_id)
    )

  # Print the map object
  print(qa_map)
}
