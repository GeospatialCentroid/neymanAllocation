# So we have to do a few important things at this script.
# aggregate the results of the model data into percent trees outside of forest. include the area in a spreadsheet 
# need to generate the estimated trees outside of forest for each M.L.R.A. for each year 
# report this as a map.
#  This includes developing the rather robust method or function for calculating the weighted average and estimated trees outside of forest that I can apply throughout the rest of the workflow.
# 
# With the known values in hand, I then want to conduct a sampling test where I randomly select 100 features at a time and calculate back out the trees outside of forest weighted average for each of the MLRA's.
# I need to append this with information about the land cover class for the MLRAs, particularly the forest cover class.


# final product 
## We've assumed that we can 



## procssing the input files. 

pacman::p_load(terra, readr, future, future.apply, dplyr)

# source functions 
source("src/processingModelData.R")


# read in required inputs 
## MLRA boundaries 
mlras <- sf::st_read("data/derived/mlra/lower48MLRA.gpkg")
## summary tables of LCC per MLRA 
mlraSummaries <- list.files(path = "data/derived/mlra_nlcd_summaryArea/", pattern = "MLRA_",
                            full.names = TRUE)
## samapling locations 
llr_f_sites <- read_csv(file = "data/products/systematicSampleSelection/selectedSample_lrr_F_05_2026.csv")


# systematic sample aoiP
# Example: shorter set of output to start testing the method 
## might need to remount the share sudo mount -t nfs 192.168.20.101:/mnt/user/fileShare /mnt/unraid_fileShare

model_dir <- "/mnt/unraid_fileShare/NAIP/modelOutputs/212022/"
model_dir <- "/mnt/unraid_fileShare/NAIP/modelOutputs/may24_runs/212022/"
model_paths <- list.files(model_dir, full.names = TRUE)
aoi_dir <- "/mnt/unraid_fileShare/NAIP/"
aoi_paths <- list.files(path = aoi_dir, pattern = ".gpkg", recursive = TRUE, full.names = TRUE)

process_model_data_fast(aoi_paths = aoi_paths, model_paths = model_paths, local_output_dir = "data/products/may2026ModelResults/")



