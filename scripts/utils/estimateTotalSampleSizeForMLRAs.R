# =========================================================================
# estimateTotalSampleSizeForMLRAs.R (Updated)
# =========================================================================

library(dplyr)
library(readr)
library(purrr)
library(fs)
library(sf)
library(ggplot2)

# -------------------------------------------------------------------------
# 1. READ AND COMBINE CSVs
# -------------------------------------------------------------------------
csv_dir <- "data/derived/mlra_nlcd_summaryArea"
csv_files <- dir_ls(csv_dir, glob = "*_attributes.csv")

all_mlra_data <- csv_files %>%
  map_dfr(read_csv, .id = "source_file") %>%
  mutate(
    # Normalize NLCD classes
    Forest_raw = rowSums(select(., matches("class_(41|42|43)")), na.rm = TRUE),
    Herbaceous_raw = rowSums(select(., matches("class_(71|72|73|74)")), na.rm = TRUE),
    Developed_raw = rowSums(select(., matches("class_(21|22|23|24)")), na.rm = TRUE),
    Cultivated_raw = rowSums(select(., matches("class_(81|82)")), na.rm = TRUE),
    Wetlands_raw = rowSums(select(., matches("class_(90|95)")), na.rm = TRUE),
    Water_raw = rowSums(select(., matches("class_(11|12)")), na.rm = TRUE),
    Barren_raw = rowSums(select(., matches("class_31")), na.rm = TRUE),
    Shrubland_raw = rowSums(select(., matches("class_(51|52)")), na.rm = TRUE),
    
    Total_Coverage = Forest_raw + Herbaceous_raw + Developed_raw + Cultivated_raw + 
      Wetlands_raw + Water_raw + Barren_raw + Shrubland_raw,
    
    Forest = (Forest_raw / Total_Coverage) * 100
  ) %>%
  # Keep only what we need to save memory
  select(source_file, MLRA_ID, Forest)

# -------------------------------------------------------------------------
# 2. AGGREGATE DATA ACROSS YEARS 
# -------------------------------------------------------------------------
mlra_aggregated <- all_mlra_data %>%
  group_by(MLRA_ID) %>%
  summarise(Forest = mean(Forest, na.rm = TRUE), .groups = "drop")

# -------------------------------------------------------------------------
# 3. LOAD SPATIAL DATA (GPKG) & FILTER TO LRRs G & H
# -------------------------------------------------------------------------
# Note: We are keeping the geometry this time so we can map it later
mlra_sf <- st_read("data/derived/mlra/lower48MLRA.gpkg", quiet = TRUE) 

mlra_joined_sf <- mlra_sf %>%
  left_join(mlra_aggregated, by = "MLRA_ID") %>%
  filter(LRRSYM %in% c("G", "H")) # Isolate target regions

# -------------------------------------------------------------------------
# 4. APPLY THE 3 SAMPLING STRATEGIES 
# -------------------------------------------------------------------------
forest_threshold <- 0.2
high_density_size <- 1400
low_density_mean_sd <- 550
low_density_max <- 650
num_years <- 3

mlra_scored_sf <- mlra_joined_sf %>%
  mutate(
    Is_High_Target = Forest <= forest_threshold,
    
    # Base Budgets per Year
    Strat1_Base = high_density_size,
    Strat2_Base = ifelse(Is_High_Target, high_density_size, low_density_mean_sd),
    Strat3_Base = ifelse(Is_High_Target, high_density_size, low_density_max),
    
    # 3-Year Total Budgets
    Strat1_Total = Strat1_Base * num_years,
    Strat2_Total = Strat2_Base * num_years,
    Strat3_Total = Strat3_Base * num_years
  )

# -------------------------------------------------------------------------
# 5. TABULAR SUMMARY FOR LRRs G & H
# -------------------------------------------------------------------------
tabular_summary <- mlra_scored_sf %>%
  st_drop_geometry() %>% # Drop geometry to cleanly aggregate tabular data
  group_by(LRRSYM) %>%
  summarise(
    Total_MLRAs = n(),
    High_Density_Count = sum(Is_High_Target, na.rm = TRUE),
    Low_Density_Count = Total_MLRAs - High_Density_Count,
    
    # Total samples across the 3 years for the entire LRR
    Strat1_Cost = sum(Strat1_Total, na.rm = TRUE),
    Strat2_Cost = sum(Strat2_Total, na.rm = TRUE),
    Strat3_Cost = sum(Strat3_Total, na.rm = TRUE),
    .groups = "drop"
  )

# Add a 'Total' row to sum G and H together
total_row <- tabular_summary %>%
  summarise(
    LRRSYM = "TOTAL (G + H)",
    Total_MLRAs = sum(Total_MLRAs),
    High_Density_Count = sum(High_Density_Count),
    Low_Density_Count = sum(Low_Density_Count),
    Strat1_Cost = sum(Strat1_Cost),
    Strat2_Cost = sum(Strat2_Cost),
    Strat3_Cost = sum(Strat3_Cost)
  )

final_summary <- bind_rows(tabular_summary, total_row)

print("--- SAMPLING COST SUMMARY (3-YEAR ESTIMATES) ---")
print(final_summary)
write_csv(final_summary, "LRR_G_H_Cost_Comparison.csv")

# -------------------------------------------------------------------------
# 6. SPATIAL SUMMARY (MAP)
# -------------------------------------------------------------------------
# Mapping Strategy 2 (1400 vs 550). 
# If you want to map Strategy 3 instead, change 'fill = factor(Strat2_Base)' to 'Strat3_Base'

map_plot <- ggplot(data = mlra_scored_sf) +
  geom_sf(aes(fill = factor(Strat2_Base)), color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("550" = "#41B6C4", "650" = "#2C7FB8", "1400" = "#E31A1C"),
    name = "Sample Budget\n(per year)"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Sampling Density Assignment (LRRs G & H)",
    subtitle = "Rule: 1400 samples if Forest <= 0.2%",
    caption = "Assigned values represent annual baseline samples."
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text = element_blank() # Hide lat/lon text for a cleaner presentation map
  )

# Save the map to disk
ggsave("LRR_G_H_Sample_Density_Map.png", plot = map_plot, width = 8, height = 6, bg = "white")
cat("\nProcess Complete! Map exported to 'LRR_G_H_Sample_Density_Map.png'\n")