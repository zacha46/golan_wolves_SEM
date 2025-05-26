# Spatial Ecology Modeling Toolkit: Wolves, Culling, and Land Use

This repository contains R scripts and data for modeling species abundance, lethal management effects, and land-use gradients in a fragmented landscape. It focuses on wolves, jackals, wild boars, and gazelles in the Golan Heights.

We implement:
- Spatial smoothing of culling/protection rasters
- N-mixture models for abundance estimation
- Structural Equation Models (SEM) to understand trophic dynamics

---

## 🧭 Overview

This toolkit includes three integrated components:

### 1. Spatial Gradient Analysis
Uses binary raster layers and exponential-decay kernels to estimate spatial intensity of culling or protection. Outputs are used as covariates in downstream models.

### 2. N-mixture Models (via `unmarked`)
Estimates species abundance and detection probabilities from camera trap data. Includes model selection and empirical Bayes estimation.

### 3. Piecewise SEM
Tests direct and indirect interactions among species and the effects of human management using generalized linear models in a SEM framework.

---

## ✅ Quickstart Checklist

### Before You Begin

1. **Do you have repeated wildlife counts across multiple sites with covariates?**
   - No → Collect appropriate data.
   - Yes → Continue.

2. **Do you know how to use N-mixture models with the `unmarked` R package?**
   - No → Learn about detection matrices and covariate formatting.
   - Yes → Continue.

3. **Are you familiar with SEM in R?**
   - No → See Kéry & Royle (2016), Chapter 5.
   - Yes → You're ready to proceed.

---

## 🗂 Repository Structure

TO BE COMPLETED! REVISIT WHEN CODE AND DATA ARE READY! 

--

## 📦 Dependencies

All scripts auto-install required packages:

- `raster`, `sp`, `rgdal`, `rgeos`
- `tidyverse`, `ggplot2`, `data.table`
- `unmarked`, `AICcmodavg`
- `piecewiseSEM`, `lme4`, `dplyr`

---

## 📁 Required Inputs

### Raster Layers (binary, `0`/`1` values):
- `boar_culling_raster.tif`, etc.
- Smoothed layers (e.g., `Boar Culling Focal.tif`) are generated from these

### CSV Files:
- **Camera Coordinates**: `Camera Coordinates CSV.csv`
- **Detection History**:
  - `Wolf Detection History for Manuscript.csv`
  - `Canis_aureus_Abundance.csv`
  - `Sus_scrofa_Abundance.csv`
  - `Gazella_gazella_Abundance.csv`
- **Covariates**:
  - `SiteCovs_060324.csv` (site-level)
  - `Detection Covariates_For_MS.csv` (observation-level)

---

## 📊 Outputs

### From Gradient Analysis:
- Smoothed raster maps (`.tif`)
- Site-level covariates (`SiteCovs.csv`)

### From N-mixture Models:
- Model comparison tables (`*_Detection_Functions_All.csv`)
- Site-level abundance estimates with SE (`*_abundance_with_SE.csv`)
- Wolf state model selection (`Wolf_state_model_selection.csv`)

### From SEM:
- Model summaries (printed)
- SEM interaction plot (`jackal_abundance_plot.png`)

---

## ⚙️ Parameters to Adjust

- `kernel_radius_m`: Radius for spatial smoothing kernel (default: `7000`)
- `decay_factor`: Exponential decay factor for distance weighting (default: `0.01`)
- `setwd()`: Update file paths in each script to match your local setup

---

## 📖 Notes

- Detection histories should have counts per site per occasion.
- Covariates must match detection matrix dimensions.
- Categorical variables must be properly factorized.
- Empirical Bayes estimates provide posterior distributions of abundance.
- Interaction terms in SEM reveal conditional effects of management (e.g., culling × land protection).

---

## 📚 Citation

If you use this repository in your work, please cite the corresponding manuscript:

> Preiss-Bloom et al. (in prep). **Human disturbance thresholds determine the ecological role of an apex predator**.

---

## 🆘 Contact

For help adapting this code or questions about the methodology, contact the lead author.