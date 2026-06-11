# ==============================================================================
# Script: 20_processDraftModelResults.R
# Purpose: Aggregate model outputs to calculate the percent Trees Outside Forest (TOF).
#          Generates area-weighted averages for each MLRA per year.
#          Applies NLCD forest class masking for subsequent analysis.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup and Initialization
# ------------------------------------------------------------------------------
pacman::p_load(terra, readr, sf, future, future.apply, dplyr, purrr, stringr, ggplot2, tigris, gt)

# Load custom processing functions
source("src/processingModelData.R")
source("src/sampleSizeValidation.R")

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

# Define the target MLRAs to process
target_mlras <- unique(llr_f_sites$MLRA_ID)


# ------------------------------------------------------------------------------
# 5. Phase III & IV: Unmasked Spatial Summarization & Aggregation
# ------------------------------------------------------------------------------
checkpoint_unmasked_file <- "data/products/detailed_tof_unmasked.csv"

if (file.exists(checkpoint_unmasked_file)) {
  message("Checkpoint found. Loading pre-calculated UNMASKED summary data...")
  detailed_tof_df <- read_csv(checkpoint_unmasked_file, show_col_types = FALSE)
  
} else {
  message("No checkpoint found. Running heavy spatial summarization for UNMASKED data...")
  
  # Execute the batch summarization using vector-anchored geometry
  final_master_df <- batch_summarize_mlras(
    mlra_ids = target_mlras,
    index_df = llr_f_sites,   
    mlra_sf = mlras,    
    local_product_dir = "data/products/may2026ModelResults/" 
  )
  
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
  
  # Save checkpoint
  write_csv(detailed_tof_df, checkpoint_unmasked_file)
  message(sprintf("Unmasked data checkpoint saved to %s", checkpoint_unmasked_file))
}

# Calculate the final area-weighted averages for each MLRA per epoch
summary_tof_df <- detailed_tof_df|>
  filter(!is.na(Rep_Year))|> 
  group_by(MLRA_ID, Rep_Year)|>
  summarize(
    Weighted_Avg_TOF_Percent = (sum(Area_Value_1_M2, na.rm = TRUE) / sum(Total_Area_M2, na.rm = TRUE)) * 100,
    Record_Count = n(),
    .groups = "drop" 
  )


# ------------------------------------------------------------------------------
# 6. Phase V & VI: NLCD Forest Masking & Aggregation
# ------------------------------------------------------------------------------
# Generate secondary batch of model results explicitly masked to NLCD forest classifications
# (This function safely skips existing files internally)
batch_apply_nlcd_mask(
  input_dir = "data/products/may2026ModelResults/",
  nlcd_dir = "data/raw/nlcd/",
  output_dir = "data/products/may2026ModelsMasked/",
  num_cores = 4
)

checkpoint_mask_file <- "data/products/detailed_tof_mask.csv"

if (file.exists(checkpoint_mask_file)) {
  message("Checkpoint found. Loading pre-calculated MASKED summary data...")
  detailed_tof_df_mask <- read_csv(checkpoint_mask_file, show_col_types = FALSE)
  
} else {
  message("No checkpoint found. Running heavy spatial summarization for MASKED data...")
  
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
  
  # Save checkpoint
  write_csv(detailed_tof_df_mask, checkpoint_mask_file)
  message(sprintf("Masked data checkpoint saved to %s", checkpoint_mask_file))
}

# Calculate the final area-weighted averages for each MLRA per epoch
summary_tof_df_mask <- detailed_tof_df_mask|>
  filter(!is.na(Rep_Year))|> 
  group_by(MLRA_ID, Rep_Year)|>
  summarize(
    Weighted_Avg_TOF_Percent = (sum(Area_Value_1_M2, na.rm = TRUE) / sum(Total_Area_M2, na.rm = TRUE)) * 100,
    Record_Count = n(),
    .groups = "drop" 
  )

# ------------------------------------------------------------------------------
# 7. Sample Size Validation Simulation
# ------------------------------------------------------------------------------

# 1. Run the Spatially Balanced Systematic Simulation
validation_df <- simulate_sample_sizes(
  detailed_df = detailed_tof_df_mask, 
  min_sample = 100, 
  step_size = 100, 
  iterations = 50,
  tolerance_pct = 0.10 # 10% accuracy tolerance
)

# 2. Calculate Absolute Error directly (No joining needed anymore!)
evaluation_df <- validation_df %>%
  mutate(Absolute_Error = abs(Estimated_TOF_Percent - True_TOF_Percent))

# 3. Aggregate the milestone thresholds (Pass/Fail Rates)
milestone_summary_df <- validation_df %>%
  group_by(MLRA_ID, Rep_Year, Sample_Size) %>%
  summarize(
    Avg_Estimate = mean(Estimated_TOF_Percent, na.rm = TRUE),
    Accuracy_Rate = mean(Is_Accurate, na.rm = TRUE),
    Coverage_Rate = mean(Is_Covered, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Passed_Accuracy_80 = Accuracy_Rate >= 0.80,
    Passed_Accuracy_95 = Accuracy_Rate >= 0.95,
    Passed_Coverage_95 = Coverage_Rate >= 0.95,
    Passed_Both_95 = Passed_Accuracy_95 & Passed_Coverage_95,
    Passed_Both_80 = Passed_Accuracy_80 & Passed_Coverage_95
  ) %>%
  arrange(MLRA_ID, Rep_Year, desc(Sample_Size))

print(head(milestone_summary_df, 10))


# ------------------------------------------------------------------------------
# 8. Presentation Mapping & Tables
# ------------------------------------------------------------------------------

# --- 8A. Generate TOF Maps ---
map_full <- generate_tof_map(
  spatial_sf = mlras, 
  tof_df = summary_tof_df_mask, 
  target_year = 2020, 
  target_size = "Full"
)

map_600 <- generate_tof_map(
  spatial_sf = mlras, 
  tof_df = validation_df, 
  target_year = 2020, 
  target_size = 600
)

map_500 <- generate_tof_map(
  spatial_sf = mlras, 
  tof_df = validation_df, 
  target_year = 2020, 
  target_size = 500
)

ggsave("map_tof_full_2020.png", plot = map_full, width = 10, height = 6, dpi = 300)
ggsave("map_tof_600_2020.png", plot = map_600, width = 10, height = 6, dpi = 300)
ggsave("map_tof_500_2020.png", plot = map_500, width = 10, height = 6, dpi = 300)

# --- 8B. Process Total Forest Cover (NLCD) per Year ---
nlcd_forest_df <- purrr::map_dfr(mlraSummaries, function(file) {
  df <- read_csv(file, show_col_types = FALSE)
  
  # Calculate forest cover explicitly for each NLCD epoch
  df <- df %>%
    mutate(
      Forest_2010 = rowSums(across(any_of(c("nlcd_2010_class_41", "nlcd_2010_class_42", "nlcd_2010_class_43"))), na.rm = TRUE),
      Forest_2016 = rowSums(across(any_of(c("nlcd_2016_class_41", "nlcd_2016_class_42", "nlcd_2016_class_43"))), na.rm = TRUE),
      Forest_2020 = rowSums(across(any_of(c("nlcd_2020_class_41", "nlcd_2020_class_42", "nlcd_2020_class_43"))), na.rm = TRUE)
    ) %>%
    select(MLRA_ID, Forest_2010, Forest_2016, Forest_2020) %>%
    # Pivot to a long format so we can join it cleanly to the validation data
    tidyr::pivot_longer(
      cols = starts_with("Forest_"),
      names_to = "NLCD_Year",
      values_to = "Percent_Forest"
    ) %>%
    mutate(
      # Map the NLCD years to your modeling Rep_Years
      Rep_Year = case_when(
        NLCD_Year == "Forest_2010" ~ 2012,
        NLCD_Year == "Forest_2016" ~ 2016,
        NLCD_Year == "Forest_2020" ~ 2020
      )
    ) %>%
    select(MLRA_ID, Rep_Year, Percent_Forest)
})

# Calculate overall MLRA averages for finding extremes and mapping
mlra_forest_avg <- nlcd_forest_df %>%
  group_by(MLRA_ID) %>%
  summarize(Avg_Percent_Forest = mean(Percent_Forest, na.rm = TRUE), .groups = "drop")

# Identify Extremes using the multi-year average
max_forest_mlra <- mlra_forest_avg %>%
  filter(MLRA_ID %in% target_mlras) %>%
  filter(Avg_Percent_Forest == max(Avg_Percent_Forest, na.rm = TRUE)) %>% pull(MLRA_ID)

min_forest_mlra <- mlra_forest_avg %>%
  filter(MLRA_ID %in% target_mlras) %>%
  filter(Avg_Percent_Forest == min(Avg_Percent_Forest, na.rm = TRUE)) %>% pull(MLRA_ID)

target_table_mlras <- c(max_forest_mlra, min_forest_mlra)


# --- 8C. Presentation Table (With Year-Specific NLCD Cover) ---
table_truth <- summary_tof_df_mask %>%
  filter(MLRA_ID %in% target_table_mlras) %>%
  select(MLRA_ID, Rep_Year, TOF_Percent = Weighted_Avg_TOF_Percent) %>%
  mutate(Sample_Size = "Full (~1400)", Accuracy_Rate = NA_real_, Coverage_Rate = NA_real_)

table_sample <- milestone_summary_df %>%
  filter(MLRA_ID %in% target_table_mlras, Sample_Size == 500) %>%
  select(MLRA_ID, Rep_Year, TOF_Percent = Avg_Estimate, Accuracy_Rate, Coverage_Rate) %>%
  mutate(Sample_Size = "Sample (500)")

presentation_table <- bind_rows(table_truth, table_sample) %>%
  # Join the year-specific forest cover directly to the table data
  left_join(nlcd_forest_df, by = c("MLRA_ID", "Rep_Year")) %>%
  arrange(MLRA_ID, Rep_Year, desc(Sample_Size)) %>%
  mutate(MLRA_Label = case_when(
    MLRA_ID == max_forest_mlra ~ paste("MLRA", MLRA_ID, "(Highest Forest Cover)"),
    MLRA_ID == min_forest_mlra ~ paste("MLRA", MLRA_ID, "(Lowest Forest Cover)")
  )) %>%
  select(MLRA_Label, Rep_Year, Sample_Size, Percent_Forest, TOF_Percent, Accuracy_Rate, Coverage_Rate) %>%
  
  gt(groupname_col = "MLRA_Label") %>%
  tab_header(
    title = md("**Methodological Validation: Effect of Sample Size**"),
    subtitle = "Comparing Full Population vs. Spatially Balanced Sub-samples (n=500)"
  ) %>%
  cols_label(
    Percent_Forest = "NLCD Forest Cover",
    TOF_Percent = "Estimated TOF",
    Accuracy_Rate = "Accuracy Rate (±10%)",
    Coverage_Rate = "CI Coverage Rate"
  ) %>%
  fmt_percent(
    columns = c(Percent_Forest, TOF_Percent, Accuracy_Rate, Coverage_Rate),
    scale_values = FALSE, 
    decimals = 2
  ) %>%
  sub_missing(
    columns = c(Accuracy_Rate, Coverage_Rate),
    missing_text = "-" 
  ) %>%
  tab_options(
    heading.background.color = "#f4f4f4",
    column_labels.font.weight = "bold",
    row_group.background.color = "#e0e0e0",
    row_group.font.weight = "bold"
  )

print(presentation_table)


## --- 8D. Generate Total Forest Cover Map (Target MLRAs & Custom Bins with Labels) ---
# Ensure the map uses the new multi-year average dataframe
forest_map_data <- mlras %>%
  mutate(MLRA_ID = as.character(MLRA_ID)) %>%
  inner_join(
    mlra_forest_avg %>% 
      filter(MLRA_ID %in% target_mlras) %>% 
      mutate(MLRA_ID = as.character(MLRA_ID)), 
    by = "MLRA_ID"
  ) %>%
  mutate(
    Forest_Class = cut(
      Avg_Percent_Forest, 
      breaks = c(-Inf, 0.2, 0.5, 1.0, 2.0, Inf), 
      labels = c("< 0.2%", "0.2 - 0.5%", "0.5 - 1.0%", "1.0 - 2.0%", "> 2.0%")
    )
  )

options(tigris_use_cache = TRUE) 
usa_states <- tigris::states(cb = TRUE, class = "sf") %>% st_transform(st_crs(forest_map_data)) 
states_cropped <- st_crop(usa_states, st_bbox(forest_map_data))

map_forest <- ggplot() +
  geom_sf(data = states_cropped, fill = "grey95", color = "grey60", linewidth = 0.5) +
  geom_sf(data = forest_map_data, aes(fill = Forest_Class), color = "white", linewidth = 0.3) +
  
  # New Layer: Automatically calculates centroids and maps the MLRA_ID as text
  geom_sf_text(data = forest_map_data, aes(label = MLRA_ID), color = "black", size = 3.5, fontface = "bold") +
  
  scale_fill_viridis_d(option = "cividis", name = "NLCD Forest\nCover") + 
  theme_void() + 
  labs(
    title = "Total NLCD Forest Cover", 
    subtitle = "Average forested area across target epochs"
  ) +
  theme(
    legend.position = "right", 
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5), 
    plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0.5),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

ggsave("map_total_forest_cover_targets.png", plot = map_forest, width = 10, height = 6, dpi = 300)



# 1. Join the year-specific forest cover to your existing table data
lowest_passing_samples_with_forest <- lowest_passing_samples %>%
  left_join(nlcd_forest_df, by = c("MLRA_ID", "Rep_Year"))

# 2. Format as a presentation-ready table
passing_samples_table <- lowest_passing_samples_with_forest %>%
  gt(groupname_col = "MLRA_ID") %>%
  tab_header(
    title = md("**Minimum Required Sample Sizes**"),
    subtitle = "Lowest tested sample size to achieve 95% CI Coverage at target accuracies"
  ) %>%
  cols_label(
    Rep_Year = "Year",
    Percent_Forest = "NLCD Forest Cover",  # <-- Added forest cover label
    Min_Sample_80 = "N required for 80% Acc.",
    Min_Sample_95 = "N required for 95% Acc."
  ) %>%
  fmt_percent(                             # <-- Formats the decimal as a clean percentage (e.g., 1.50%)
    columns = Percent_Forest,
    scale_values = FALSE, 
    decimals = 2
  ) %>%
  sub_missing(
    columns = c(Min_Sample_80, Min_Sample_95),
    missing_text = "Did Not Pass"
  ) %>%
  tab_options(
    heading.background.color = "#f4f4f4",
    column_labels.font.weight = "bold"
  )

print(passing_samples_table)
