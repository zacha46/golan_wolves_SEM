## 3. Piecewise structural equation modeling (SEM) of trophic relationships

# Install required packages if not already installed
required_packages <- c("piecewiseSEM", "lme4", "ggplot2", "tidyverse")

installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load required libraries
library(piecewiseSEM)
library(lme4)
library(ggplot2)
library(tidyverse)
library(dplyr)

## Make sure nothing is in your global environment
rm(list=ls(all=TRUE))


## Set working directory
# setwd("C:/Users/preissbloom/OneDrive - Tel-Aviv University/Golan Wolf Study/Wolf Reports and Excels") # for Shlomo
setwd("~/Dropbox/UQ 2024/Collaboration_Shlomo/Shlomo_R_code_for_git_via_zach/")                       # for Zach 


##Import site-level abundance data, which is already in our sitecovs dataframe from Part 2
siteCovs <- read.csv("SiteCovs_060324.csv")


#define variables to use

siteCovs2 <- select(siteCovs,Jackal.Culling,Boar.Culling.Focal.Corrected,Wolf.Culling.Last.1.Year,Null.Wolf.Abundance, null.jackal.abundance,null.boar.abundance,Nature.Reserves,gazelleabundance,Wolf_10.SE,Jackal_10.SE,Gazelle_10.SE,Boar_10.SE)
new_colnames <- c("JackalCulling", "BoarCulling", "WolfCulling", "wolfabundance", "jackalabundance", "boarabundance", "NatureReserves", "gazelleabundance", "WolfSE", "JackalSE", "GazelleSE", "BoarSE")
colnames(siteCovs2) <- new_colnames

#Construct GLMs with culling and predator pressures, without interactions
wolf <- glm(wolfabundance ~ WolfCulling + NatureReserves, family = gaussian, weights = WolfSE, siteCovs2)
jackal <- glm(jackalabundance ~ wolfabundance + NatureReserves + JackalCulling , family = gaussian, weights = JackalSE, siteCovs2)
boar <- glm(boarabundance ~ wolfabundance + NatureReserves + BoarCulling + jackalabundance , family = gaussian, weights = BoarSE, siteCovs2)
gazelle <- glm(gazelleabundance ~ wolfabundance + NatureReserves + jackalabundance , family = gaussian, weights = GazelleSE, siteCovs2)

#Construct GLMs with culling and predator pressures, including interactions between land use and culling/predator pressures
wolf_1 <- glm(wolfabundance ~ WolfCulling * NatureReserves, family = gaussian, weights = WolfSE, siteCovs2)
jackal_1 <- glm(jackalabundance ~ wolfabundance + NatureReserves * JackalCulling , family = gaussian, weights = JackalSE, siteCovs2)
jackal_2 <- glm(jackalabundance ~ wolfabundance * NatureReserves + JackalCulling , family = gaussian, weights = JackalSE, siteCovs2)
boar_1 <- glm(boarabundance ~ wolfabundance + NatureReserves * BoarCulling + jackalabundance , family = gaussian, weights = BoarSE, siteCovs2)
boar_2 <- glm(boarabundance ~ wolfabundance * NatureReserves + BoarCulling + jackalabundance , family = gaussian, weights = BoarSE, siteCovs2)
gazelle_1 <- glm(gazelleabundance ~ wolfabundance * NatureReserves + jackalabundance , family = gaussian, weights = GazelleSE, siteCovs2)
gazelle_2 <- glm(gazelleabundance ~ wolfabundance + NatureReserves * jackalabundance , family = gaussian, weights = GazelleSE, siteCovs2)

#view Results of individual models

summary(wolf_1)
summary(jackal)
summary(boar)
summary(gazelle)

# Compare jackal models
AIC(jackal_1, jackal_2)

# Compare boar models
AIC(boar_1, boar_2)

# Compare gazelle models
AIC(gazelle_1, gazelle_2)

#Continue to constructing SEMs using the top GLMs for each species

#Construct SEMs with wolf model
model_1 <- psem(wolf, jackal, boar, gazelle)
model_2 <- psem(wolf_1, jackal_1, boar_1, gazelle_1)

#Construct SEMs without wolf model
model_3 <- psem(jackal, boar, gazelle)
model_4 <- psem(jackal_1, boar_1, gazelle_1)

#With culling x land use interaction but none for gazelles
model_5 <- psem(wolf_1, jackal_1, boar_1, gazelle)
model_6 <- psem(jackal_1, boar_1, gazelle)

# Summarize model outputs including path coefficients, significance, and fit statistics
# Fisher's C statistic tests for missing paths (independence claims)summary(model_1)
summary(model_2)
summary(model_3)
summary(model_4)
summary(model_5)
summary(model_6)

AIC(model_1, model_2, model_3, model_4, model_5, model_6)



#To interpret interactive effects, we will plot them using high and low levels of land protection.
# Plot the predicted relationship between jackal culling and jackal abundance
# at two levels of land use (low and high natural cover)
# Based on interaction model 'jackal_1'


# Define 1st and 75th Quantiles for NatureReserves (25th and 75th quantiles)
quantiles <- quantile(siteCovs$Nature.Reserves, probs = c(0.25, 0.75))

# Create a dataset for prediction across a range of JackalCulling values
prediction_data <- expand.grid(
  JackalCulling = seq(min(siteCovs2$JackalCulling), max(siteCovs2$JackalCulling), length.out = 100),
  NatureReserves = quantiles,
  wolfabundance = mean(siteCovs2$wolfabundance)  # Keeping wolfabundance constant (can vary if needed)
)

# Predict jackal abundance using the existing GLM model, including standard errors
predictions <- predict(jackal_1, newdata = prediction_data, type = "response", se.fit = TRUE)

# Add the predicted abundance and standard errors to the data
prediction_data$predicted_abundance <- predictions$fit
prediction_data$se.fit <- predictions$se.fit

# Calculate the 95% confidence intervals
prediction_data$lower_CI <- prediction_data$predicted_abundance - 1.96 * prediction_data$se.fit
prediction_data$upper_CI <- prediction_data$predicted_abundance + 1.96 * prediction_data$se.fit

# Convert NatureReserves to a factor for coloring in the plot
prediction_data$NatureReserves <- factor(prediction_data$NatureReserves, 
                                         labels = c("Less Natural (25th Quantile)", "More Natural (75th Quantile)"))

# Step 1: Subset predictions for the two quantiles (25th and 75th quantile)
pred_25th <- prediction_data[prediction_data$NatureReserves == "Less Natural (25th Quantile)", ]
pred_75th <- prediction_data[prediction_data$NatureReserves == "More Natural (75th Quantile)", ]

# Step 2: Find the intersection point (where the difference between the predicted abundances is zero)
# Calculate the absolute difference between predicted abundances for the two quantiles
difference <- abs(pred_25th$predicted_abundance - pred_75th$predicted_abundance)

# Find the index of the minimum difference (closest to zero)
intersection_index <- which.min(difference)

# Get the corresponding x-value (JackalCulling) and y-value (predicted_abundance) at the intersection
intersection_x <- pred_25th$JackalCulling[intersection_index]  # x-value where the intersection occurs
intersection_y <- pred_25th$predicted_abundance[intersection_index]  # y-value at the intersection

# Create the plot with confidence intervals, placing the legend inside the plot area
jackals <- ggplot(prediction_data, aes(x = JackalCulling, y = predicted_abundance, color = NatureReserves, fill = NatureReserves)) +
  geom_line(size = 1) +  # Lines for each level of NatureReserves
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, linetype = 0) +  # Confidence intervals
  scale_color_manual(values = c("Less Natural (25th Quantile)" = "red", "More Natural (75th Quantile)" = "aquamarine3")) +  # Aquamarine3 for green
  scale_fill_manual(values = c("Less Natural (25th Quantile)" = "red", "More Natural (75th Quantile)" = "aquamarine3")) +   # Aquamarine3 for fill
  labs(x = "Jackal Culling Intensity",
       y = "Predicted Jackal Abundance",
       fill = "Land Use:\n(p-value: 0.0153, β = -0.629)", 
       color = "Land Use:\n(p-value: 0.0153, β = -0.629)") +  # Combine fill and color labels
  ylim(0, 1) +  # Set y-axis limits between 0 and 1
  theme_minimal() +
  theme(legend.position = c(0.3, 0.8), 
        legend.background = element_rect(fill = "white", color = "black"),  # White background for the legend
        legend.title = element_text(size = 12),  # Adjust legend title size if needed
        axis.title = element_text(size = 14)) +  # Increase axis label size
  guides(fill = guide_legend(override.aes = list(alpha = 0.2))) +  # Combine legends
  coord_fixed(ratio = 1) +  # Make the plot area square
  geom_hline(yintercept = intersection_y, linetype = "dotted", color = "black", size = 1) +
  geom_text(aes(x = max(prediction_data$JackalCulling), y = intersection_y, label = "*Threshold (65th Quantile)"),
            color = "black", hjust = 1, vjust = 1.5, size = 4)


# Save the plot in high resolution (300 DPI)
ggsave("jackal_abundance_plot.png", plot = jackals, width = 6, height = 6, dpi = 300)

#Find threshold value of nature reserves at which the slope of culling is zero

# Extract coefficients from the jackal_1 model
coef_jackal <- coef(jackal_1)

# Coefficients for JackalCulling and its interaction with NatureReserves
beta_jackalCulling <- coef_jackal["JackalCulling"]
beta_interaction <- coef_jackal["NatureReserves:JackalCulling"]

# Solve for the NatureReserves value where the marginal effect of JackalCulling is zero
threshold_nature_reserves <- -beta_jackalCulling / beta_interaction

# Print the result
cat("Threshold value of NatureReserves where the slope of JackalCulling is zero:", threshold_nature_reserves, "\n")

#Find the quantile of nature reserves at which the slope of culling is zero
q_threshold <- ecdf(siteCovs$Nature.Reserves)(threshold_nature_reserves)


#Let's plot interactions: BOAR----

# Step 1: Calculate the threshold where the slope of BoarCulling is zero

# Extract coefficients from the boar_1 model
coef_boar <- coef(boar_1)

# Coefficients for BoarCulling and its interaction with NatureReserves
beta_boarCulling <- coef_boar["BoarCulling"]
beta_interaction_boar <- coef_boar["NatureReserves:BoarCulling"]

# Solve for the BoarCulling value where the marginal effect of BoarCulling is zero
threshold_boarCulling <- -beta_boarCulling / beta_interaction_boar

#Find the quantile of nature reserves at which the slope of culling is zero
q_threshold_boar <- ecdf(siteCovs$Nature.Reserves)(threshold_boarCulling)

# Define 1st and 75th Quantiles for NatureReserves (25th and 75th quantiles)

# Create a grid for predictions
prediction_data_boar <- expand.grid(
  BoarCulling = seq(min(siteCovs2$BoarCulling), max(siteCovs2$BoarCulling), length.out = 100),
  wolfabundance = mean(siteCovs2$wolfabundance),
  jackalabundance = mean(siteCovs2$jackalabundance),
  NatureReserves = quantile(siteCovs2$NatureReserves, probs = c(0.25, 0.75)) # 25th and 75th quantiles
)


# Predict boar abundance using the existing GLM model, including standard errors
predictions_boar <- predict(boar_1, newdata = prediction_data_boar, type = "response", se.fit = TRUE)

# Add the predicted abundance and standard errors to the data
prediction_data_boar$predicted_abundance <- predictions_boar$fit
prediction_data_boar$se.fit <- predictions_boar$se.fit

# Calculate the 95% confidence intervals
prediction_data_boar$lower_CI <- prediction_data_boar$predicted_abundance - 1.96 * prediction_data_boar$se.fit
prediction_data_boar$upper_CI <- prediction_data_boar$predicted_abundance + 1.96 * prediction_data_boar$se.fit

# Convert NatureReserves to a factor for coloring in the plot
prediction_data_boar$NatureReserves <- factor(prediction_data_boar$NatureReserves, 
                                              labels = c("Less Natural (25th Quantile)", "More Natural (75th Quantile)"))

# Step 1: Subset predictions for the two quantiles (25th and 75th quantile)
pred_25th_boar <- prediction_data_boar[prediction_data_boar$NatureReserves == "Less Natural (25th Quantile)", ]
pred_75th_boar <- prediction_data_boar[prediction_data_boar$NatureReserves == "More Natural (75th Quantile)", ]

# Step 2: Find the intersection point (where the difference between the predicted abundances is zero)
difference_boar <- abs(pred_25th_boar$predicted_abundance - pred_75th_boar$predicted_abundance)

# Find the index of the minimum difference (closest to zero)
intersection_index_boar <- which.min(difference_boar)

# Get the corresponding x-value (BoarCulling) and y-value (predicted_abundance) at the intersection
intersection_x_boar <- pred_25th_boar$BoarCulling[intersection_index_boar]  # x-value where the intersection occurs
intersection_y_boar <- pred_25th_boar$predicted_abundance[intersection_index_boar]  # y-value at the intersection

boars <- ggplot(prediction_data_boar, aes(x = BoarCulling, y = predicted_abundance, color = NatureReserves, fill = NatureReserves)) +
  geom_line(size = 1) +  # Lines for each level of NatureReserves
  geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), alpha = 0.2, linetype = 0) +  # Confidence intervals
  scale_color_manual(values = c("Less Natural (25th Quantile)" = "red", "More Natural (75th Quantile)" = "aquamarine3")) + 
  scale_fill_manual(values = c("Less Natural (25th Quantile)" = "red", "More Natural (75th Quantile)" = "aquamarine3")) +
  labs(x = "Boar Culling Intensity",
       y = "Predicted Boar Abundance",
       fill = "Land Use:\n(p-value: 0.0117, β = 0.4839)", 
       color = "Land Use:\n(p-value: 0.0117, β = 0.4839)") + 
  theme_minimal() +
  theme(legend.position = c(0.3, 0.8), 
        legend.background = element_rect(fill = "white", color = "black"),  # White background for the legend
        legend.title = element_text(size = 12), 
        axis.title = element_text(size = 14)) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.2))) +
  coord_fixed(ratio = 1) +  # Make the plot area square
  geom_hline(yintercept = intersection_y_boar, linetype = "dotted", color = "black", size = 1) +  # Horizontal line at intersection
  geom_text(aes(x = max(prediction_data_boar$BoarCulling), y = intersection_y_boar, label = "*Threshold (62nd Quantile)"),
            color = "black", hjust = 1, vjust = -0.5, size = 4) +
  coord_cartesian(ylim = c(0, 1))  # Ensures display limit without cutting off elements


# Save the plot in high resolution (300 DPI)
ggsave("boar_abundance_plot.png", plot = boars, width = 6, height = 6, dpi = 300)


#Plot for gazelles----


# Generate a sequence for NatureReserves and quantiles for wolfabundance
nature_reserves_seq <- seq(min(siteCovs2$NatureReserves), max(siteCovs2$NatureReserves), length.out = 100)
wolf_activity_quantiles <- quantile(siteCovs2$wolfabundance, c(0.25, 0.75))  # 25th and 75th quantiles

# Create a data frame with combinations of NatureReserves and wolfabundance quantiles
predictions_df <- expand.grid(NatureReserves = nature_reserves_seq, 
                              wolfabundance = wolf_activity_quantiles,
                              jackalabundance = mean(siteCovs2$jackalabundance))

# Generate predictions and confidence intervals
predictions_df$fit <- predict(gazelle_1, newdata = predictions_df, type = "response", se.fit = TRUE)$fit
predictions_df$ci.lb <- predictions_df$fit - 
  1.96 * predict(gazelle_1, newdata = predictions_df, type = "response", se.fit = TRUE)$se.fit
predictions_df$ci.ub <- predictions_df$fit + 
  1.96 * predict(gazelle_1, newdata = predictions_df, type = "response", se.fit = TRUE)$se.fit

# Convert wolfabundance to a factor with updated labels
predictions_df$wolfabundance <- factor(predictions_df$wolfabundance, 
                                       labels = c("Low Wolf Activity (25th Quantile)", 
                                                  "High Wolf Activity (75th Quantile)"))

# Set y-axis limits to 0-1
y_limits <- c(0, 1)
x_limits <- c(0,1)


# Step 1: Subset predictions for the two quantiles (25th and 75th quantile)
pred_25th <- predictions_df[predictions_df$wolfabundance == "Low Wolf Activity (25th Quantile)", ]
pred_75th <- predictions_df[predictions_df$wolfabundance == "High Wolf Activity (75th Quantile)", ]

# Step 2: Find the intersection point (where the difference between the predicted abundances is zero)
# Calculate the absolute difference between predicted abundances for the two quantiles
difference <- abs(pred_25th$fit - pred_75th$fit)

# Find the index of the minimum difference (closest to zero)
intersection_index <- which.min(difference)

# Get the corresponding x-value (NatureReserves) at the intersection
intersection_x_gazelle <- pred_25th$NatureReserves[intersection_index]

# Plot the interaction between Nature Reserves and wolf activity on gazelle abundance
gazelles <- ggplot(predictions_df, aes(x = NatureReserves, y = fit, color = wolfabundance, fill = wolfabundance)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), alpha = 0.2, colour=NA) +  # Confidence intervals
  scale_color_manual(values = c("darkgoldenrod2", "darkblue")) +  # Colors for the two wolf activity quantiles
  scale_fill_manual(values = c("darkgoldenrod2", "darkblue")) +   # Same colors for the fill
  labs(x = "Nature Reserves (Size and Proximity)", y = "Predicted Gazelle Abundance") +
  theme_minimal(base_size = 15) +
  theme(
    legend.position = c(0.3, 0.8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    panel.border = element_rect(color = "white", fill = NA, size = 1),  # Frame for the plot
    legend.background = element_rect(fill = "white", color = "black"),  # Black frame around legend
    aspect.ratio = 1  # Make the plot area square
  ) +
  ylim(y_limits) +  # Set the y-axis limits from 0 to 1
  guides(fill = guide_legend(title = "Wolf Activity:\n(p = 0.0067, β = 1.13)", 
                             labels = c("Low Wolf Activity (25th Quantile)", 
                                        "High Wolf Activity (75th Quantile)"))) +
  guides(color = guide_legend(title = "Wolf Activity:\n(p = 0.0067, β = 1.13)", 
                              labels = c("Low Wolf Activity (25th Quantile)", 
                                         "High Wolf Activity (75th Quantile)"))) +
  geom_vline(xintercept = intersection_x_gazelle, linetype = "dotted", color = "black", size = 1) +  # Add vertical line at intersection
  geom_text(aes(x = intersection_x_gazelle +0.025, y = 0.25, label = "*Threshold (0.114)"),
            color = "black", vjust = 1.5, hjust = 0, size = 4)  # Caption to the right of the vertical line
ggsave("gazelle_abundance_plot.png", plot = gazelles, width = 6, height = 6, dpi = 300)
