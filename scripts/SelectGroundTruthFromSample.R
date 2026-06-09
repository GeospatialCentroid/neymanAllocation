# 
pacman::p_load(dplyr, readr)

# for LLR F 
## read in existing sample locations 
## allocate a set number of sample per MLRA based on number of remaining sites 
## draw random sample from the each mlra grouping 
lrrF <- read_csv("data/products/systematicSampleSelection/selectedSample_lrr_F_05_2026.csv")
# drawing and additional 60 sites. 
sites_per_mlra <- ceiling(60/length(unique(lrrF$MLRA_ID)))
# random draw from each MLRA 
set.seed(1234)
sampled_data <- lrrF %>%
  group_by(MLRA_ID) %>%
  slice_sample(n = 6) %>%
  ungroup()
#export 
write_csv(sampled_data, "data/products/groundTruthSiteSelection/lrr_F_final60_GT_sites.csv")

# for LLR G and H 
## read in existing sample locations 
## attribute NLCD data to these locations 
## distribute samples among the MLRA's (proportionally?), then apply the neyman allocation to select from the groups 

