

library(terra)
source("scripts/00_config.R")

# Define the IDs to highlight
target_ids <- c(78, 86, 150)

s1 <- terra::vect(STATIC_INPUTS$mlra)

# Create a subset of the vector for just those IDs
# Note: We convert target_ids to numeric because the MLRA_ID column is <num>
s1_highlight <- s1[s1$MLRA_ID %in% as.numeric(target_ids), ]

# 1. Plot the Base Map (All MLRAs)
# col = fill color, border = outline color
plot(s1, col = "gray90", border = "gray60", 
     main = "Selected MLRAs (78, 86, 150)", axes = FALSE)

# 2. Add the Highlighted Layer (add = TRUE)
plot(s1_highlight, col = "tomato", border = "black", lwd = 2, add = TRUE)

# 3. Add Labels
# terra::text() automatically calculates centroids for polygons
text(s1_highlight, labels = s1_highlight$MLRA_ID, 
     col = "black", halo = TRUE, cex = 1.2, font = 2)



