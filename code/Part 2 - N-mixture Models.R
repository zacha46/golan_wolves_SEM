#####PART 2: N-mixture models (Detection and State) for species abundance
## Last Updated: 03 June 2025

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
 setwd("C:/Users/preissbloom/OneDrive - Tel-Aviv University/Golan Wolf Study/Paper - Trophic Thunder/Golan Codes and Dataframes for Github")   # for Shlomo
#setwd("~/Dropbox/UQ 2024/Collaboration_Shlomo/Shlomo_R_code_for_git_via_zach/")                # for Zach 

#Import detection histories for the four species
y_w <- read.csv("Canis_Lupus_Abundance.csv")
y_j <- read.csv("Canis_aureus_Abundance.csv")
y_b <- read.csv("Sus_scrofa_Abundance.csv")
y_g <- read.csv("Gazella_gazella_Abundance.csv")

# Remove the first column (camera name) from y_w and y_j
y_w <- y_w[, -1]
y_j <- y_j[, -1]

#Load site and observation covariates
siteCovs <- read.csv("SiteCovs_060324.csv")
# obsCovs <- read.csv("Detection Covariates_For_MS.csv") # not present
obsCovs <- read.csv("ObsCovs_050524.csv")


#Make sure all siteCovs and obsCovs characters are factors for unmarked
siteCovs[sapply(siteCovs, is.character)] <- lapply(siteCovs[sapply(siteCovs, is.character)], as.factor)
obsCovs[sapply(obsCovs, is.character)] <- lapply(obsCovs[sapply(obsCovs, is.character)], as.factor)


# Construct unmarked frames
umf_pcount_w <- unmarkedFramePCount(y = y_w, siteCovs = siteCovs, obsCovs = obsCovs) 
umf_pcount_j <- unmarkedFramePCount(y = y_j, siteCovs = siteCovs, obsCovs = obsCovs) 
umf_pcount_b <- unmarkedFramePCount(y = y_b, siteCovs = siteCovs, obsCovs = obsCovs)
umf_pcount_g <- unmarkedFramePCount(y = y_g, siteCovs = siteCovs, obsCovs = obsCovs)

# Test WOLF null model
Wolf_null <- pcount(~1 ~1, umf_pcount_w, mixture = "ZIP")

# Test WOLF detection functions with individual and multiple covariates

  Agri = pcount(~Agricultural.Presence ~1, umf_pcount_w, mixture = "ZIP")
  Cattle = pcount(~Cattle.Sightings ~1, umf_pcount_w, mixture = "ZIP")
  Model = pcount(~Camera.Model ~1, umf_pcount_w, mixture = "ZIP")
  Delay = pcount(~Capture.Delay ~1, umf_pcount_w, mixture = "ZIP")
  Veg = pcount(~Vegetation.Structure ~1, umf_pcount_w, mixture = "ZIP")
  Path = pcount(~Path.Type ~1, umf_pcount_w, mixture = "ZIP")
  BehindCam = pcount(~Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP")
  Bottleneck = pcount(~Bottleneck. ~1, umf_pcount_w, mixture = "ZIP")
  wd2 = pcount(~Path.Type + Camera.Model + Bottleneck. ~ 1, umf_pcount_w, mixture = "ZIP")
  wd3 = pcount(~Path.Type + Camera.Model ~1, umf_pcount_w, mixture = "ZIP")
  wd4 = pcount(~Path.Type + Camera.Model + Cattle.Sightings ~1, umf_pcount_w, mixture = "ZIP")
  wd5 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP")
  wd6 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_w, mixture = "ZIP")
  wd7 = pcount(~Cattle.Sightings + Path.Type + Camera.Model ~1, umf_pcount_w, mixture = "ZIP")
  wd8 = pcount(~Cattle.Sightings + Camera.Model + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP")
  wd85 = pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP")
  wd87 = pcount(~Cattle.Sightings + Camera.Model + Agricultural.Presence + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP")
  wd9 = pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_w, mixture = "ZIP")
  wd10 = pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~ 1, umf_pcount_w, mixture = "ZIP")

  
  # Combine into a named list for model selection
  wolf_detection_models <- list(Null = Wolf_null, Agri = Agri, Cattle = Cattle, 
    Model = Model, Delay = Delay, Veg = Veg, Path = Path, BehindCam = BehindCam,
    Bottleneck = Bottleneck, wd2 = wd2, wd3 = wd3, wd4 = wd4, wd5 = wd5, wd6 = wd6,
    wd7 = wd7, wd8 = wd8, wd85 = wd85, wd87 = wd87, wd9 = wd9, wd10 = wd10 )
  
  # Rank models by AIC
  wolf_detection_models_ranked <- fitList(fits = wolf_detection_models) %>%
    modSel(nullmod = "Null") %>%
    as("data.frame")
  
  # Optional: write ranked models to CSV
  # write.csv(wolf_detection_models_ranked, "wolf_detection_models_ranked.csv", row.names = FALSE)  

# Jackal detection function
 
   # Test JACKAL null model
  Jackal_null <- pcount(~1 ~1, umf_pcount_j, mixture = "ZIP")
  
  # Test JACKAL detection functions with individual and multiple covariates
  j_Agri        <- pcount(~Agricultural.Presence ~1, umf_pcount_j, mixture = "ZIP")
  j_Cattle      <- pcount(~Cattle.Sightings ~1, umf_pcount_j, mixture = "ZIP")
  j_Model       <- pcount(~Camera.Model ~1, umf_pcount_j, mixture = "ZIP")
  j_Delay       <- pcount(~Capture.Delay ~1, umf_pcount_j, mixture = "ZIP")
  j_Veg         <- pcount(~Vegetation.Structure ~1, umf_pcount_j, mixture = "ZIP")
  j_Path        <- pcount(~Path.Type ~1, umf_pcount_j, mixture = "ZIP")
  j_BehindCam   <- pcount(~Area.Behind.Camera ~1, umf_pcount_j, mixture = "ZIP")
  j_Bottleneck  <- pcount(~Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  
  # Combined detection models
  jm2   <- pcount(~Path.Type + Camera.Model + Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  jd3   <- pcount(~Path.Type + Camera.Model ~1, umf_pcount_j, mixture = "ZIP")
  jd4   <- pcount(~Path.Type + Camera.Model + Cattle.Sightings ~1, umf_pcount_j, mixture = "ZIP")
  jd5   <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_j, mixture = "ZIP")
  jd6   <- pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_j, mixture = "ZIP")
  jd7   <- pcount(~Cattle.Sightings + Path.Type + Camera.Model ~1, umf_pcount_j, mixture = "ZIP")
  jd8   <- pcount(~Cattle.Sightings + Camera.Model + Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  jd85  <- pcount(~Cattle.Sightings + Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  jd87  <- pcount(~Cattle.Sightings + Camera.Model + Agricultural.Presence + Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  jd9   <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_j, mixture = "ZIP")
  jd10  <- pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~1, umf_pcount_j, mixture = "ZIP")
  
  
  jackal_detection_models <- list(
    Null = Jackal_null, j_Agri = j_Agri, j_Cattle = j_Cattle, j_Model = j_Model, 
    j_Delay = j_Delay, j_Veg = j_Veg, j_Path = j_Path, j_BehindCam = j_BehindCam, 
    j_Bottleneck = j_Bottleneck, jd2 = jd2, jd3 = jd3, jd4 = jd4, jd5 = jd5, 
    jd6 = jd6, jd7 = jd7, jd8 = jd8, jd85 = jd85, jd87 = jd87, jd9 = jd9, jd10 = jd10
  )
  
  jackal_detection_models_ranked <- fitList(fits = jackal_detection_models) %>%
    modSel(nullmod = "Null") %>%
    as("data.frame")
  

  # Optional: write ranked models to CSV
  # write.csv(jackal_detection_models_ranked, "jackal_detection_models_ranked.csv", row.names = FALSE)  
  
# BOAR detection function

  # Test BOAR null model
  Boar_null <- pcount(~1 ~1, umf_pcount_b, mixture = "ZIP")
  
  # Test BOAR detection functions
  boar_Cattle     <- pcount(~Cattle.Sightings ~1, umf_pcount_b, mixture = "ZIP")
  boar_Model      <- pcount(~Camera.Model ~1, umf_pcount_b, mixture = "ZIP")
  boar_Sens       <- pcount(~Camera.Sensitivity ~1, umf_pcount_b, mixture = "ZIP")
  boar_Delay      <- pcount(~Capture.Delay ~1, umf_pcount_b, mixture = "ZIP")
  boar_Veg        <- pcount(~Vegetation.Structure ~1, umf_pcount_b, mixture = "ZIP")
  boar_Path       <- pcount(~Path.Type ~1, umf_pcount_b, mixture = "ZIP")
  boar_BehindCam  <- pcount(~Area.Behind.Camera ~1, umf_pcount_b, mixture = "ZIP")
  boar_Bottleneck <- pcount(~Bottleneck. ~1, umf_pcount_b, mixture = "ZIP")
  boar_Livestock  <- pcount(~Livestock.Present. ~1, umf_pcount_b, mixture = "ZIP")
  boar_Agri       <- pcount(~Agricultural.Presence ~1, umf_pcount_b, mixture = "ZIP")
  bd2  <- pcount(~Path.Type + Camera.Model ~1, umf_pcount_b, mixture = "ZIP")
  bd3  <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_b, mixture = "ZIP")
  bd4  <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_b, mixture = "ZIP")
  bd5  <- pcount(~Path.Type + Camera.Model + Vegetation.Structure ~1, umf_pcount_b, mixture = "ZIP")
  bd6  <- pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~1, umf_pcount_b, mixture = "ZIP")
  
  # Combine into a named list for model selection
  boar_detection_models <- list(
    Null = Boar_null, Cattle = boar_Cattle, Model = boar_Model, Sens = boar_Sens,
    Delay = boar_Delay, Veg = boar_Veg, Path = boar_Path, BehindCam = boar_BehindCam,
    Bottleneck = boar_Bottleneck, Livestock = boar_Livestock, Agri = boar_Agri,
    bd2 = bd2, bd3 = bd3, bd4 = bd4, bd5 = bd5, bd6 = bd6
  )
  
  # Rank models by AIC
  boar_detection_models_ranked <- fitList(fits = boar_detection_models) %>%
    modSel(nullmod = "Null") %>%
    as("data.frame")
  
  # Optional: write ranked models to CSV
  # write.csv(boar_detection_models_ranked, "boar_detection_models_ranked.csv", row.names = FALSE)  

# GAZELLE detection function

  # Test GAZELLE null model
  Gazelle_null <- pcount(~1 ~1, umf_pcount_g, mixture = "ZIP")
  
  # Test GAZELLE detection functions
  gaz_Cattle     <- pcount(~Cattle.Sightings ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Model      <- pcount(~Camera.Model ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Sens       <- pcount(~Camera.Sensitivity ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Livestock  <- pcount(~Livestock.Present. ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Delay      <- pcount(~Capture.Delay ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Veg        <- pcount(~Vegetation.Structure ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Path       <- pcount(~Path.Type ~1, umf_pcount_g, mixture = "ZIP")
  gaz_BehindCam  <- pcount(~Area.Behind.Camera ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Bottleneck <- pcount(~Bottleneck. ~1, umf_pcount_g, mixture = "ZIP")
  gaz_Agri       <- pcount(~Agricultural.Presence ~1, umf_pcount_g, mixture = "ZIP")
  gd3  <- pcount(~Path.Type + Camera.Model ~1, umf_pcount_g, mixture = "ZIP")
  gd5  <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera ~1, umf_pcount_g, mixture = "ZIP")
  gd8  <- pcount(~Path.Type + Camera.Model + Area.Behind.Camera + Bottleneck. ~1, umf_pcount_g, mixture = "ZIP")
  gd10 <- pcount(~Path.Type + Camera.Model + Vegetation.Structure + Area.Behind.Camera ~1, umf_pcount_g, mixture = "ZIP")
  gd11 <- pcount(~Path.Type + Camera.Model + Vegetation.Structure ~1, umf_pcount_g, mixture = "ZIP")
  
  # Combine into a named list for model selection
  gazelle_detection_models <- list(
    Null = Gazelle_null, Cattle = gaz_Cattle, Model = gaz_Model, Sens = gaz_Sens,
    Livestock = gaz_Livestock, Delay = gaz_Delay, Veg = gaz_Veg, Path = gaz_Path,
    BehindCam = gaz_BehindCam, Bottleneck = gaz_Bottleneck, Agri = gaz_Agri,
    gd3 = gd3, gd5 = gd5, gd8 = gd8, gd10 = gd10, gd11 = gd11
  )
  
  # Rank models by AIC
  gazelle_detection_models_ranked <- fitList(fits = gazelle_detection_models) %>%
    modSel(nullmod = "Null") %>%
    as("data.frame")
  
  # Optional: write ranked models to CSV
  # write.csv(gazelle_detection_models_ranked, "gazelle_detection_models_ranked.csv", row.names = FALSE)  

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

# Estimate site activity for WOLF, JACKAL, BOAR and GAZELLE using the best detection function for each species

wolf_activity <- estimate_site_activity(wd10, "Wolf")
jackal_activity <- estimate_site_activity(jd10, "Jackal")
boar_activity <- estimate_site_activity(bd6, "Boar")
gazelle_activity <- estimate_site_activity(gd10, "Gazelle")

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
obs.boot <- Nmix.gof.test(wm393, nsim = 100)
obs.boot
print(obs.boot, digits.vals = 4, digits.chisq = 4)
