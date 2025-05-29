
########################################################################
## ANALYSIS CODE FOR: [PAPER TITLE]
## Authors: [AUTHORS]
## Last Updated: May 27, 2025
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
#personal_path = "~/Dropbox/UQ 2024/Collaboration_Shlomo/Shlomo_R_code_for_git_via_zach/"
personal_path = "C:/Users/preissbloom/OneDrive - Tel-Aviv University/Golan Wolf Study/Paper - Trophic Thunder/Golan Codes and Dataframes for Github/"


# Load raster (example: boar culling, nature reserves etc.) where all values are either 0 or 1
r1 <- raster(paste(personal_path,"boar_culling_raster.tif", sep = "")) 
r2 <- raster(paste(personal_path,"nature_reserves_raster.tif", sep = "")) 

# Replace NA with 0 just in case some values are NA 
r1[is.na(r1[])] <- 0
r2[is.na(r2[])] <- 0


# Parameters for kernel creation
kernel_radius_m <- 7000       # Radius in meters
decay_factor <- 0.01          # Controls how quickly the effect diminishes with distance

# Define a circular kernel (7km radius) with exponential decay
d1 <- focalWeight(r1, kernel_radius_m, "circle")
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
the_focus_1 <- focal(r1, my_circle, fun = sum, na.rm = FALSE, pad = TRUE, padValue = 0)
plot(the_focus_1)
writeRaster(the_focus_1, "Boar Culling Focal.tif")

# Apply spatial filter using custom kernel
the_focus_2 <- focal(r2, my_circle, fun = sum, na.rm = FALSE, pad = TRUE, padValue = 0)
plot(the_focus_2)
writeRaster(the_focus_2, "Nature Reserves Culling Focal.tif")

#Once you have all your rasters, stack them together and extract the values at the coordinates of each site

#Stack rasters
Rastack <- stack(r1,r2)

#Load camera trap coordinates
pointCoordinates <- read.csv(r"(Camera Coordinates CSV.csv)")

#Convert to spatial points
coordinates(pointCoordinates) <- ~ Longitude + Latitude

#Extract raster values at camera locations
rasvalue <- raster::extract(Rastack, pointCoordinates)

# Combine point coordinates with extracted values
combinePointValue <- cbind(pointCoordinates,rasvalue)

# Save the combined data and voila! You have your site-level covariates file, for use in unmarked or other models
write.table(combinePointValue, file = "SiteCovs.csv", append=FALSE, sep= ",", row.names = FALSE, col.names = TRUE)




