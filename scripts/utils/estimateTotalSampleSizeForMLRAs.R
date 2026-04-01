# Load necessary libraries (Make sure to install 'sf' if you haven't already: install.packages("sf"))
library(dplyr)
library(readr)
library(purrr)
library(fs)
library(sf)

# -------------------------------------------------------------------------
# 1. READ AND COMBINE CSVs
# -------------------------------------------------------------------------
# Update this directory to wherever your MLRA CSVs are located
csv_dir <- "data/derived/mlra_nlcd_summaryArea"

# Update glob to explicitly look for *_attributes.csv
csv_files <- dir_ls(csv_dir, glob = "*_attributes.csv")

all_mlra_data <- csv_files %>%
  map_dfr(read_csv, .id = "source_file")

# -------------------------------------------------------------------------
# 2. COMBINE NLCD CLASSES BASED ON MRLC LEGEND (NORMALIZED TO 100%)
# -------------------------------------------------------------------------
all_mlra_data <- all_mlra_data %>%
  mutate(
    # Step A: Calculate the raw sums.
    # Because 'matches()' captures columns for 2010, 2016, and 2020 simultaneously,
    # these raw sums might add up to ~300% per row.
    Forest_raw = rowSums(select(., matches("class_(41|42|43)")), na.rm = TRUE),
    Herbaceous_raw = rowSums(
      select(., matches("class_(71|72|73|74)")),
      na.rm = TRUE
    ),
    Developed_raw = rowSums(
      select(., matches("class_(21|22|23|24)")),
      na.rm = TRUE
    ),
    Cultivated_raw = rowSums(select(., matches("class_(81|82)")), na.rm = TRUE),
    Wetlands_raw = rowSums(select(., matches("class_(90|95)")), na.rm = TRUE),
    Water_raw = rowSums(select(., matches("class_(11|12)")), na.rm = TRUE),
    Barren_raw = rowSums(select(., matches("class_31")), na.rm = TRUE),
    Shrubland_raw = rowSums(select(., matches("class_(51|52)")), na.rm = TRUE),

    # Step B: Calculate the total percentage captured in this row
    Total_Coverage = Forest_raw +
      Herbaceous_raw +
      Developed_raw +
      Cultivated_raw +
      Wetlands_raw +
      Water_raw +
      Barren_raw +
      Shrubland_raw,

    # Step C: Normalize back to a true 0-100% scale
    Forest = (Forest_raw / Total_Coverage) * 100,
    Herbaceous = (Herbaceous_raw / Total_Coverage) * 100,
    Developed = (Developed_raw / Total_Coverage) * 100,
    Cultivated = (Cultivated_raw / Total_Coverage) * 100,
    Wetlands = (Wetlands_raw / Total_Coverage) * 100,
    Water = (Water_raw / Total_Coverage) * 100,
    Barren = (Barren_raw / Total_Coverage) * 100,
    Shrubland = (Shrubland_raw / Total_Coverage) * 100
  ) %>%
  # Keep only the cleaned, normalized percentage columns
  select(
    source_file,
    MLRA_ID,
    Forest,
    Herbaceous,
    Developed,
    Cultivated,
    Wetlands,
    Water,
    Barren,
    Shrubland
  )


# -------------------------------------------------------------------------
# 3. AGGREGATE DATA ACROSS YEARS (Mean per MLRA)
# -------------------------------------------------------------------------
mlra_aggregated <- all_mlra_data %>%
  group_by(MLRA_ID) %>%
  summarise(across(
    # Now, the only numeric columns left to average are Forest, Herbaceous, etc.
    .cols = where(is.numeric) & !matches("year|source_file"),
    .fns = ~ mean(.x, na.rm = TRUE)
  )) %>%
  ungroup()

# -------------------------------------------------------------------------
# 4. LOAD SPATIAL DATA (GPKG) & JOIN
# -------------------------------------------------------------------------
# Extract the attributes from the GPKG without keeping the heavy geometries
mlra_sf <- st_read("data/derived/mlra/lower48MLRA.gpkg", quiet = TRUE) %>%
  st_drop_geometry() %>%
  select(MLRA_ID, LRRSYM, LRR_NAME)

# Join spatial attributes to our aggregated stats
mlra_joined <- mlra_aggregated %>%
  left_join(mlra_sf, by = "MLRA_ID")

# -------------------------------------------------------------------------
# 5. APPLY SAMPLING STRATEGY LOGIC
# -------------------------------------------------------------------------
# Updated parameters based on your requirements
forest_threshold <- 0.5
herbaceous_threshold <- 85.0 # Assumed high threshold based on earlier MLRA 79 logic
# 900 when include all MLRA, or 600 when exclude those that would be sampled at 1400
mean_plus_sd_size <- 600

high_density_size <- 1400
num_years <- 3

mlra_scored <- mlra_joined %>%
  mutate(
    # Assign base sample size per year per MLRA
    base_sample_size = case_when(
      Herbaceous >= herbaceous_threshold |
        Forest <= forest_threshold ~ high_density_size,
      TRUE ~ mean_plus_sd_size
    ),

    # Calculate the total sample size across the 3 years for that individual MLRA
    mlra_total_samples = base_sample_size * num_years
  )

# -------------------------------------------------------------------------
# 6. GROUP BY LRR AND CALCULATE FINAL TOTAL SAMPLES
# -------------------------------------------------------------------------
# This creates the final summary grouped by the LRR Synonym
lrr_summary <- mlra_scored %>%
  group_by(LRRSYM, LRR_NAME) %>%
  summarise(
    Total_MLRAs_in_LRR = n(),
    Total_Samples_Required = sum(mlra_total_samples, na.rm = TRUE),
    Avg_Forest_Pct = mean(Forest, na.rm = TRUE),
    Avg_Herbaceous_Pct = mean(Herbaceous, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------------------------------------------------------
# 7. EXPORT THE RESULTS
# -------------------------------------------------------------------------
print(head(lrr_summary))

# Save the grouped summary out
write_csv(lrr_summary, "LRR_Sample_Requirements_Summary.csv")

# Save the individual MLRA table as a backup to inspect the individual MLRA math
write_csv(mlra_scored, "MLRA_Individual_Sample_Requirements.csv")
