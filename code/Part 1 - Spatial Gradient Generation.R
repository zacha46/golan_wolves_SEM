
########################################################################
## ANALYSIS CODE FOR: [PAPER TITLE]
## Authors: [AUTHORS]
## Last Updated: May 13, 2025
########################################################################
## This script performs spatial gradient analysis for culling pressure and land use
## The code creates spatially smoothed layers using circular kernels with exponential decay
## and extracts values at camera trap locations
########################################################################

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  raster,      # For raster operations
  rgdal,       # For spatial data operations
  tidyverse,   # For data manipulation and visualization
  sp,          # For spatial data structures
  rgeos,       # For geometric operations
  ggplot2      # For visualization
)

### Assign where you have stored your data
personal_path = "~/Dropbox/UQ 2024/Collaboration_Shlomo/Shlomo_R_code_for_git_via_zach/"


# Load raster (example: boar culling, nature reserves etc.) where all values are either 0 or 1
r <- raster(paste(personal_path,"boar_culling_raster.tif", sep = "")) ## NOT PRESENT! CANNOT TEST! 

# Replace NA with 0 and save corrected raster
r[is.na(r[])] <- 0
writeRaster(r, "D:/Golan Wolves GIS/0_1_Rasters/Shootings_10YR_TimeDecay_0_1.tif")

# Parameters for kernel creation
kernel_radius_m <- 7000       # Radius in meters
decay_factor <- 0.01          # Controls how quickly the effect diminishes with distance

# Define a circular kernel (7km radius) with exponential decay
d1 <- focalWeight(r, kernel_radius_m, "circle")
d1[d1 != 0] <- 1  # Normalize to binary circle

# Create exponential decay weights based on 25m resolution
d1_new_matrix <- outer(seq(-280, 280, 1), seq(-280, 280, 1), function(x, y) exp(-decay_factor * sqrt(x^2 + y^2)))
my_circle <- d1 * d1_new_matrix

# Visualize kernel shape
my_circle %>%
  as_tibble() %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  mutate(name = parse_number(name), rowname = as.numeric(rowname)) %>%
  ggplot(aes(x = name, y = rowname, fill = value)) +
  geom_tile()

# Apply spatial filter using custom kernel
the_focus <- focal(r, my_circle, fun = sum, na.rm = FALSE, pad = TRUE, padValue = 0)
plot(the_focus)
writeRaster(the_focus, "Boar Culling Focal.tif")


#Once you have all your rasters, stack them together and extract the values at the coordinates of each site

#Stack rasters
Rastack10 <- stack(r3,r7,r9,r10,r14)

#Load camera trap coordinates
pointCoordinates <- read.csv(r"(Camera Coordinates CSV.csv)")

#Convert to spatial points
coordinates(pointCoordinates) <- ~ Longitude+Latitude

#Extract raster values at camera locations
rasvalue <- raster::extract(Rastack, pointCoordinates)

# Combine point coordinates with extracted values
combinePointValue <- cbind(pointCoordinates,rasvalue)

# Save the combined data and voila! You have your site-level covariates file, for use in unmarked or other models
write.table(combinePointValue, file = "SiteCovs.csv", append=FALSE, sep= ",", row.names = FALSE, col.names = TRUE)




