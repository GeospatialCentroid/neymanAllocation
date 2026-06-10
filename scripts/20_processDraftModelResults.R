# ==============================================================================
# Script: 20_processDraftModelResults.R
# Purpose: Aggregate model outputs to calculate the percent Trees Outside Forest (TOF).
#          Generates area-weighted averages for each MLRA per year.
#          Applies NLCD forest class masking for subsequent analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Initialization
# ------------------------------------------------------------------------------
pacman::p_load(terra, readr, future, future.apply, dplyr, purrr, stringr)

# Load custom processing functions
source("src/processingModelData.R")

# ------------------------------------------------------------------------------
# 2. Data Ingestion
# ------------------------------------------------------------------------------
# Load MLRA boundaries and spatial data
mlras <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg")

# Load LCC summary tables per MLRA
mlraSummaries <- list.files(
  path = "data/derived/mlra_nlcd_summaryArea/", 
  pattern = "MLRA_",
  full.names = TRUE
)

# Load systematic sampling locations
llr_f_sites <- read_csv(file = "data/products/systematicSampleSelection/selectedSample_lrr_F_05_2026.csv")

# Define network paths (Requires active Unraid share mount)
model_dir <- "/mnt/unraid_fileShare/NAIP/modelOutputs/may24_runs/212022/"
model_paths <- list.files(model_dir, full.names = TRUE)

aoi_dir <- "/mnt/unraid_fileShare/NAIP/"
aoi_paths <- list.files(path = aoi_dir, pattern = "\\.gpkg$", recursive = TRUE, full.names = TRUE)


# ------------------------------------------------------------------------------
# 3. Phase I: Local Data Generation
# ------------------------------------------------------------------------------
# Transfer and crop model data from the network share to the local drive.
# Kept at 8 cores to balance processing speed with memory constraints.
process_model_data_fast(
  aoi_paths = aoi_paths,
  model_paths = model_paths, 
  local_output_dir = "data/products/may2026ModelResults/",
  num_cores = 8 
)


# ------------------------------------------------------------------------------
# 4. Phase II: Index Preparation
# ------------------------------------------------------------------------------
# Bind the physical .gpkg file paths to the main sampling index for geometric area calculations.

# Build a lookup table using a flexible regex to handle alphanumeric ID structures
id_regex <- "[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+-[a-zA-Z0-9]+"
aoi_matches <- regexpr(id_regex, aoi_paths)
valid_aois <- aoi_matches != -1

aoi_lookup <- data.frame(
  id = regmatches(aoi_paths, aoi_matches), 
  aoi_path = aoi_paths[valid_aois],
  stringsAsFactors = FALSE
)|> 
  distinct(id, .keep_all = TRUE)

# Join the physical paths back to the sampling dataframe
llr_f_sites <- llr_f_sites|>
  left_join(aoi_lookup, by = "id")


# ------------------------------------------------------------------------------
# 5. Phase III: Spatial Summarization (Unmasked)
# ------------------------------------------------------------------------------
# Define the target MLRAs to process (Defaults to all unique MLRAs in the index)
target_mlras <- unique(llr_f_sites$MLRA_ID)

# Execute the batch summarization using vector-anchored geometry
final_master_df <- batch_summarize_mlras(
  mlra_ids = target_mlras,
  index_df = llr_f_sites,   
  mlra_sf = mlras,    
  local_product_dir = "data/products/may2026ModelResults/" 
)


# ------------------------------------------------------------------------------
# 6. Phase IV: Statistical Aggregation
# ------------------------------------------------------------------------------
# Calculate individual site percentages and assign representative temporal epochs
detailed_tof_df <- final_master_df|>
  mutate(
    Rep_Year = case_when(
      Year %in% c("2010", "2011", "2012", "2013") ~ 2012,
      Year %in% c("2014", "2015", "2016", "2017") ~ 2016,
      Year %in% c("2018", "2019", "2020", "2021") ~ 2020,
      TRUE ~ NA_real_ 
    ),
    Percent_TOF = (Area_Value_1_M2 / Total_Area_M2) * 100
  )

# Calculate the final area-weighted averages for each MLRA per epoch
summary_tof_df <- detailed_tof_df|>
  filter(!is.na(Rep_Year))|> 
  group_by(MLRA_ID, Rep_Year)|>
  summarize(
    Weighted_Avg_TOF_Percent = (sum(Area_Value_1_M2, na.rm = TRUE) / sum(Total_Area_M2, na.rm = TRUE)) * 100,
    Record_Count = n(),
    .groups = "drop" 
  )

# Optional: Print summaries to console for review
# print(head(detailed_tof_df))
# print(head(summary_tof_df))


# ------------------------------------------------------------------------------
# 7. Phase V: NLCD Forest Masking
# ------------------------------------------------------------------------------
# Generate a secondary batch of model results explicitly masked to NLCD forest classifications
batch_apply_nlcd_mask(
  input_dir = "data/products/may2026ModelResults/",
  nlcd_dir = "data/raw/nlcd/",
  output_dir = "data/products/may2026ModelsMasked/",
  num_cores = 4
)
# Execute the batch summarization using vector-anchored geometry
final_master_df_mask <- batch_summarize_mlras(
  mlra_ids = target_mlras,
  index_df = llr_f_sites,   
  mlra_sf = mlras,    
  local_product_dir = "data/products/may2026ModelsMasked/" 
)


# Calculate individual site percentages and assign representative temporal epochs
detailed_tof_df_mask <- final_master_df_mask|>
  mutate(
    Rep_Year = case_when(
      Year %in% c("2010", "2011", "2012", "2013") ~ 2012,
      Year %in% c("2014", "2015", "2016", "2017") ~ 2016,
      Year %in% c("2018", "2019", "2020", "2021") ~ 2020,
      TRUE ~ NA_real_ 
    ),
    Percent_TOF = (Area_Value_1_M2 / Total_Area_M2) * 100
  )

# Calculate the final area-weighted averages for each MLRA per epoch
summary_tof_df_mask <- detailed_tof_df_mask|>
  filter(!is.na(Rep_Year))|> 
  group_by(MLRA_ID, Rep_Year)|>
  summarize(
    Weighted_Avg_TOF_Percent = (sum(Area_Value_1_M2, na.rm = TRUE) / sum(Total_Area_M2, na.rm = TRUE)) * 100,
    Record_Count = n(),
    .groups = "drop" 
  )


### 

