#####PART 2: N-mixture models (Detection and State) for species abundance

# Install required packages if not already installed
required_packages <- c("unmarked", "data.table", "tidyverse", "raster", "ggplot2", "AICcmodavg")

installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load required libraries
library(unmarked)
library(data.table)
library(tidyverse)
library(raster)
library(ggplot2)
library(AICcmodavg)


## Make sure nothing is in your global environment
rm(list=ls(all=TRUE))

# Load camera trap detection and covariate data
# setwd("D:/ShlomoPB/OneDrive - Tel-Aviv University/Golan Wolf Study/Wolf Reports and Excels")   # for Shlomo
setwd("~/Dropbox/UQ 2024/Collaboration_Shlomo/Shlomo_R_code_for_git_via_zach/")                # for Zach 

#Import detection histories for the four species
# y_w <- read.csv("Wolf Detection History for Manuscript.csv") # not present
y_w <- read.csv("Wolf Detection History Unmarked Format.csv")
y_j <- read.csv("Canis_aureus_Abundance.csv")
y_b <- read.csv("Sus_scrofa_Abundance.csv")
y_g <- read.csv("Gazella_gazella_Abundance.csv")

#Load site and observation covariates
siteCovs <- read.csv("SiteCovs_060324.csv")
# obsCovs <- read.csv("Detection Covariates_For_MS.csv") # not present
obsCovs <- read.csv("Detection Covariates.csv")


#Make sure all siteCovs and obsCovs characters are factors for unmarked
siteCovs[sapply(siteCovs, is.character)] <- lapply(siteCovs[sapply(siteCovs, is.character)], as.factor)
obsCovs[sapply(obsCovs, is.character)] <- lapply(obsCovs[sapply(obsCovs, is.character)], as.factor)


# Construct unmarked frames
umf_pcount_w <- unmarkedFramePCount(y = y_w, siteCovs = siteCovs, obsCovs = obsCovs) # error --> cant proceed w/ testing
umf_pcount_j <- unmarkedFramePCount(y = y_j, siteCovs = siteCovs, obsCovs = obsCovs) # error --> cant proceed w/ testing
umf_pcount_b <- unmarkedFramePCount(y = y_b, siteCovs = siteCovs, obsCovs = obsCovs)
umf_pcount_g <- unmarkedFramePCount(y = y_g, siteCovs = siteCovs, obsCovs = obsCovs)

# Test WOLF null model
wm1 <- pcount(~1 ~1, umf_pcount_w, mixture = "ZIP")

# Test WOLF detection functions with individual and multiple covariates
models_detection <- list(
  Agri = pcount(~Agricultural.Presence ~1, umf_pcount_w, mixture = "ZIP"),
  Cattle = pcount(~Cattle.Sightings ~1, umf_pcount_w, mixture = "ZIP"),
  Model = pcount(~Camera.Model ~1, umf_pcount_w, mixture = "ZIP"),
  Delay = pcount(~Capture.Delay ~1, umf_pcount_w, mixture = "ZIP"),
  Veg = pcount(~Vegetation.Structure ~1, umf_pcount_w, mixture = "ZIP"),
  Path = pcount(~Path.Type ~1, umf_pcount_w, mixture = "ZIP"),
  BehindCam = pcount(~Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP"),
  Bottleneck = pcount(~Bottleneck. ~1, umf_pcount_w, mixture = "ZIP"),
  wm2 = pcount(~Path.Type + Camera.Model + Bottleneck. ~ 1, umf_pcount_w, mixture = "ZIP"),
  wm3 = pcount(~Path.Type + Camera.Model ~1, umf_pcount_w, mixture = "ZIP"),
  wm4 = pcount(~Path.Type + Camera.Model + Cattle.Sightings ~1, umf_pcount_w, mixture = "ZIP"),
  wm5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP"),
  wm6 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP"),
  wm7 = pcount(~Cattle.Sightings + Path.Type + Camera.Model ~1, umf_pcount_w, mixture = "ZIP"),
  wm8 = pcount(~Cattle.Sightings + Camera.Model + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP"),
  wm85 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP"),
  wm87 = pcount(~Cattle.Sightings + Camera.Model + Agricultural.Presence + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP"),
  wm9 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP"),
  wm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_w, mixture = "ZIP")
)

# Compare models
fitList(c(Null = list(wm1), models_detection)) %>%
  modSel(nullmod = "Null") %>%
  as("data.frame") %>%
  write.csv("Wolf_Detection_Functions_All.csv", row.names = FALSE)

# Jackal detection function

# Test JACKAL null model
jm1 <- pcount(~1 ~1, umf_pcount_j)

jackal_models_detection <- list(
  Cattle = pcount(~Cattle.Sightings ~1, umf_pcount_j),
  Model = pcount(~Camera.Model ~1, umf_pcount_j),
  Sensitivity = pcount(~Camera.Sensitivity ~1, umf_pcount_j),
  Delay = pcount(~Capture.Delay ~1, umf_pcount_j),
  Veg = pcount(~Vegetation.Structure ~1, umf_pcount_j),
  Path = pcount(~Path.Type ~1, umf_pcount_j),
  BehindCam = pcount(~Area.Behind.Camera ~1, umf_pcount_j),
  Bottleneck = pcount(~Bottleneck. ~1, umf_pcount_j),
  jm2 = pcount(~Path.Type + Camera.Model + Bottleneck. ~ 1, umf_pcount_j),
  jm3 = pcount(~Path.Type + Camera.Model ~1, umf_pcount_j),
  jm4 = pcount(~Path.Type + Camera.Model + Cattle.Sightings ~1, umf_pcount_j),
  jm5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_j),
  jm6 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_j),
  jm7 = pcount(~Cattle.Sightings + Path.Type + Camera.Model ~1, umf_pcount_j),
  jm8 = pcount(~Cattle.Sightings + Camera.Model + Bottleneck. ~1, umf_pcount_j),
  jm85 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_j),
  jm87 = pcount(~Cattle.Sightings + Camera.Model + Agricultural.Presence + Bottleneck. ~1, umf_pcount_j),
  jm9 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_j),
  jm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_j)
)

# Compare all jackal detection models
fitList(c(Null = list(jm1), jackal_models_detection)) %>%
  modSel(nullmod = "Null") %>%
  as("data.frame") %>%
  write.csv("Jackal_Detection_Functions_All.csv", row.names = FALSE)

# BOAR detection function

# Boar detection function

# Null model
bm1 <- pcount(~1 ~1, umf_pcount_b)

# Boar detection models
boar_models_detection <- list(
  Cattle = pcount(~Cattle.Sightings ~1, umf_pcount_b),
  Model = pcount(~Camera.Model ~1, umf_pcount_b),
  Sensitivity = pcount(~Camera.Sensitivity ~1, umf_pcount_b),
  Delay = pcount(~Capture.Delay ~1, umf_pcount_b),
  Veg = pcount(~Vegetation.Structure ~1, umf_pcount_b),
  Path = pcount(~Path.Type ~1, umf_pcount_b),
  BehindCam = pcount(~Area.Behind.Camera ~1, umf_pcount_b),
  Bottleneck = pcount(~Bottleneck. ~1, umf_pcount_b),
  Livestock = pcount(~Livestock.Present. ~1, umf_pcount_b),
  Agri = pcount(~Agricultural.Presence ~1, umf_pcount_b),
  
  bm2 = pcount(~Path.Type + Camera.Model ~1, umf_pcount_b),
  bm3 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_b),
  bm4 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_b),
  bm5 = pcount(~Path.Type + Camera.Model + Vegetation.Structure ~1, umf_pcount_b),
  bm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~1, umf_pcount_b)
)

# Compare all boar detection models
fitList(c(Null = list(bm1), boar_models_detection)) %>% ## error --> cant save
  modSel(nullmod = "Null") %>%
  as("data.frame") %>%
  write.csv("Boar_Detection_Models_All.csv", row.names = FALSE)



# GAZELLE detection function


# Test GAZELLE null model
g20 <- pcount(~1 ~1, umf_pcount_g)

gazelle_models_detection <- list(
  g21 = pcount(~Cattle.Sightings ~1, umf_pcount_g),
  g22 = pcount(~Camera.Model ~1, umf_pcount_g),
  g23 = pcount(~Camera.Sensitivity ~1, umf_pcount_g),
  g24 = pcount(~Livestock.Present. ~1, umf_pcount_g),
  g25 = pcount(~Capture.Delay ~1, umf_pcount_g),
  g26 = pcount(~Vegetation.Structure ~1, umf_pcount_g),
  g27 = pcount(~Path.Type ~1, umf_pcount_g),
  g28 = pcount(~Area.Behind.Camera ~1, umf_pcount_g),
  g29 = pcount(~Bottleneck. ~1, umf_pcount_g),
  g30 = pcount(~Agricultural.Presence ~1, umf_pcount_g),
  g3 = pcount(~Path.Type + Camera.Model ~1, umf_pcount_g),
  g5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_g),
  g8 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_g),
  g10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~1, umf_pcount_g),
  g11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure ~1, umf_pcount_g)
)

# Compare all gazelle detection models
fitList(c(Null = list(g20), gazelle_models_detection)) %>%
  modSel(nullmod = "Null") %>%
  as("data.frame") %>%
  write.csv("Gazelle_Detection_Models.csv", row.names = FALSE)




######## Generate site-level abundance estimates for each species at each site #####

# Define the wrapper function that takes the detection function model and species name as arguments

estimate_site_activity <- function(model, species_name, nsims = 500) {
  # Run empirical Bayes estimation
  re <- ranef(model)
  
  # Draw posterior samples
  ppd <- posteriorSamples(re, nsims = 500)
  samples <- ppd@samples
  samples_df <- as.data.frame(samples)
  
  # Function to calculate mean, SE, and 95% credible interval
  calculate_stats <- function(row_values) {
    mean_val <- mean(row_values)
    sd_val <- sd(row_values)
    ci_low <- quantile(row_values, 0.025)
    ci_high <- quantile(row_values, 0.975)
    return(c(mean = mean_val, SE = sd_val, ci_low = ci_low, ci_high = ci_high))
  }
  
  # Apply stats function to each site (row)
  stats_matrix <- t(apply(samples_df, 1, calculate_stats))
  
  # Convert to data frame
  activity_df <- as.data.frame(stats_matrix)
  
  # Add normalized activity column
  activity_df$activity_norm <- activity_df$mean / max(activity_df$mean)
  
  # Write to CSV
  out_filename <- paste0(species_name, "_null_abundance_with_SE.csv")
  write.csv(activity_df, out_filename, row.names = FALSE)
  
  # Return the data frame (optional)
  return(activity_df)
}

# Estimate site activity for WOLF, JACKAL, BOAR and GAZELLE

wolf_activity <- estimate_site_activity(Wolf_models_state$Null, "Wolf")
jackal_activity <- estimate_site_activity(jackal_models_state$Null, "Jackal")
boar_activity <- estimate_site_activity(boar_models_state$Null, "Boar")
gazelle_activity <- estimate_site_activity(gazelle_models_state$Null, "Gazelle")


#These estimates and their SEs get added to the SiteCovs data frame

###### Integrated N-mixture models for species abundance ######

# Build and compare WOLF state functions using the best detection function as the null model
Wolf_models_state <- list(
  
  Null   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera ~ 1, umf_pcount_w, mixture = "ZIP"),
  
  # Culling alone
  wm11   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year, umf_pcount_w, mixture = "ZIP"),
  
  # Culling x Minefields
  wm30   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Minefields, umf_pcount_w, mixture = "ZIP"),
  wm31   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Minefields, umf_pcount_w, mixture = "ZIP"),
  
  # Culling x Nature Reserves
  wm32   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm33   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm12   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Nature.Reserves + Open.Spaces, umf_pcount_w, mixture = "ZIP"),
  wm13   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Nature.Reserves + Built.and.Roads, umf_pcount_w, mixture = "ZIP"),
  
  # Culling + Minefields + Nature Reserves (+ Built & Roads in later models)
  wm34   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Minefields + Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm35   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Minefields * Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm36   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Minefields + Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm37   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Nature.Reserves + Minefields, umf_pcount_w, mixture = "ZIP"),
  wm47   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Minefields + Built.and.Roads, umf_pcount_w, mixture = "ZIP"),
  
  # Open Spaces interactions
  wm38   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Open.Spaces, umf_pcount_w, mixture = "ZIP"),
  wm39   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Open.Spaces, umf_pcount_w, mixture = "ZIP"),
  wm391  = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Open.Spaces + Nature.Reserves, umf_pcount_w, mixture = "ZIP"),
  wm392  = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Open.Spaces + Minefields, umf_pcount_w, mixture = "ZIP"),
  wm393  = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Open.Spaces + Built.and.Roads, umf_pcount_w, mixture = "ZIP"),
  
  # Other land uses
  wm41   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Proximity.to.Water, umf_pcount_w, mixture = "ZIP"),
  wm42   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Orchards, umf_pcount_w, mixture = "ZIP"),
  wm421  = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Orchards + Wolf.Culling.Last.1.Year, umf_pcount_w, mixture = "ZIP"),
  wm422  = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Orchards * Wolf.Culling.Last.1.Year, umf_pcount_w, mixture = "ZIP"),
  wm43   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Cropland, umf_pcount_w, mixture = "ZIP"),
  wm44   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year * Built.and.Roads, umf_pcount_w, mixture = "ZIP"),
  wm45   = pcount(~Vegetation.Structure + Path.Type + Camera.Model + Area.Behind.Camera 
                  ~ Wolf.Culling.Last.1.Year + Built.and.Roads, umf_pcount_w, mixture = "ZIP")
  
)

# Build fitList
Wolf_models_state_fit <- fitList(fits = Wolf_models_state)

# Model selection and export for WOLF state models
Wolf_models_state_selection <- modSel(Wolf_models_state_fit)

# View the top-ranked models
print(Wolf_models_state_selection)

# Export model selection table to CSV
wolf_state_model_selection_df <- as(Wolf_models_state_selection, "data.frame")
write.csv(wolf_state_model_selection_df, "Wolf_state_model_selection.csv", row.names = FALSE)


# Fit jackal models with consistent detection structure
Jackal_models_state <- list(
  
  jNull = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_j),
  
  j1 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling, umf_pcount_j),
  
  jm1 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance, umf_pcount_j),
  
  jm3 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling, umf_pcount_j),
  
  jm31 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling * Null.Wolf.Abundance, umf_pcount_j),
  
  jm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling + Null.Wolf.Abundance, umf_pcount_j),
  
  jm9 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Open.Spaces, umf_pcount_j),
  
  jm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves, umf_pcount_j),
  
  jm101 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling * Nature.Reserves + Null.Wolf.Abundance, umf_pcount_j),
  
  jm102 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling + Nature.Reserves * Null.Wolf.Abundance, umf_pcount_j),
  
  jm11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads, umf_pcount_j),
  
  jm12 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Minefields, umf_pcount_j),
  
  jm13 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads * Jackal.Culling, umf_pcount_j),
  
  jm15 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves * Null.Wolf.Abundance, umf_pcount_j),
  
  jm16 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads + Jackal.Culling, umf_pcount_j),
  
  jm17 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads + Null.Wolf.Abundance, umf_pcount_j),
  
  jm18 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads * Jackal.Culling + Null.Wolf.Abundance, umf_pcount_j),
  
  jm19 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Animal.Structures, umf_pcount_j),
  
  jm20 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Animal.Structures + Nature.Reserves + Open.Spaces + Jackal.Culling * Null.Wolf.Abundance, umf_pcount_j),
  
  jm21 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads + Jackal.Culling + Null.Wolf.Abundance, umf_pcount_j),
  
  jm22 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling + Nature.Reserves * Null.Wolf.Abundance, umf_pcount_j),
  
  jm23 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Jackal.Culling * Open.Spaces + Nature.Reserves * Null.Wolf.Abundance, umf_pcount_j)
)

# Fit and compare models
Jackal_models_state_fit <- fitList(fits = Jackal_models_state)
Jackal_models_state_selection <- modSel(Jackal_models_state_fit)

# Print and export
print(Jackal_models_state_selection)
write.csv(as(Jackal_models_state_selection, "data.frame"), "Jackal_state_model_selection.csv", row.names = FALSE)



# Fit boar abundance models with consistent detection structure
Boar_models_state <- list(
  
  bmNull = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_b),
  
  bm2 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Boar.Culling.Focal.Corrected, umf_pcount_b),
  
  bm3 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves, umf_pcount_b),
  
  bm4 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves * Null.Wolf.Abundance, umf_pcount_b),
  
  bm5 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Minefields, umf_pcount_b),
  
  bm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.jackal.abundance, umf_pcount_b),
  
  bm7 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance, umf_pcount_b),
  
  bm8 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + Boar.Culling.Focal.Corrected * Nature.Reserves + null.jackal.abundance, umf_pcount_b),
  
  bm9 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + Boar.Culling.Focal.Corrected + Nature.Reserves, umf_pcount_b),
  
  bm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Orchards, umf_pcount_b),
  
  bm11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cropland, umf_pcount_b),
  
  bm12 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Open.Spaces, umf_pcount_b),
  
  bm13 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads, umf_pcount_b),
  
  bm14 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cattle.Abundance, umf_pcount_b)
)

# Fit and compare models
Boar_models_state_fit <- fitList(fits = Boar_models_state)
Boar_models_state_selection <- modSel(Boar_models_state_fit)

# Print and export
print(Boar_models_state_selection)
write.csv(as(Boar_models_state_selection, "data.frame"), "Boar_state_model_selection.csv", row.names = FALSE)

#Gazelle abundance models

# Fit gazelle models with consistent detection structure
Gazelle_models_state <- list(
  
  gm0 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_g),
  
  gm1 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Predicted.Wolf.Abundance, umf_pcount_g),
  
  gm2 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance, umf_pcount_g),
  
  gm3 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Built.and.Roads, umf_pcount_g),
  
  gm4 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Open.Spaces, umf_pcount_g),
  
  gm5 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.jackal.abundance, umf_pcount_g),
  
  gm6 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Nature.Reserves, umf_pcount_g),
  
  gm7 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Firing.Zones, umf_pcount_g),
  
  gm8 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Orchards, umf_pcount_g),
  
  gm9 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cropland, umf_pcount_g),
  
  gm10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Cattle.Abundance, umf_pcount_g),
  
  gm11 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance * Nature.Reserves + null.jackal.abundance, umf_pcount_g),
  
  gm12 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + Nature.Reserves * null.jackal.abundance, umf_pcount_g),
  
  gm13 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ Null.Wolf.Abundance + null.boar.abundance + null.jackal.abundance, umf_pcount_g),
  
  gm14 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ null.boar.abundance * Nature.Reserves, umf_pcount_g)
)

# Fit and compare models
Gazelle_models_state_fit <- fitList(fits = Gazelle_models_state)
Gazelle_models_state_selection <- modSel(Gazelle_models_state_fit, nullmod = "gm0")

# Print and export
print(Gazelle_models_state_selection)
write.csv(as(Gazelle_models_state_selection, "data.frame"), "Gazelle_Abundance_Models.csv", row.names = FALSE)

#Now that we've ranked all the abundance models, we can check their goodness of fit

##compute observed chi-square, assess significance, and estimate c-hat
obs.boot <- Nmix.gof.test(gm11, nsim = 100)
obs.boot
print(obs.boot, digits.vals = 4, digits.chisq = 4)
