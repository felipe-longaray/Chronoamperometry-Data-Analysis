# Chronoamperometry Data Analysis & Diffusion Fitting ⚡🧪

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Data processing: PalmSens](https://img.shields.io/badge/Instrument-PalmSens_PSTrace-lightgrey.svg)]()

This repository contains an automated R pipeline designed to process, clean, and analyze raw chronoamperometry data exported directly from PalmSens PSTrace software. 

Rather than relying on manual spreadsheet manipulation, this script securely handles data provenance, extracts specific transient data points, performs statistical group analysis (ANOVA/t-tests), and fits a power-law decay model to evaluate diffusion kinetics.

## Physics & Materials Context
In standard chronoamperometry, the current response for a diffusion-controlled process is described by the **Cottrell equation**, where current decays proportionally to the inverse square root of time (I ∝ t^(-0.5)). However, in complex materials or nanostructured systems, transport often deviates from classical Fickian behavior, resulting in **anomalous diffusion**. This script fits the transient current decay to a generalized power-law model:

**I = a * t^b**

By applying a logarithmic transformation (`ln(I) = ln(a) + b*ln(t)`), the script extracts:
* **a (Pre-exponential factor):** Related to the apparent rate constant and active surface area.
* **b (Decay exponent):** The kinetic signature. A value of b = -0.5 indicates ideal Cottrellian diffusion, while deviations indicate anomalous transport phenomena.

## Key Features
* **Automated Parsing:** Dynamically standardizes alternating Time/Current columns from `.xlsx` exports.
* **Log-Linear Regression:** Computes the `a` and `b` diffusion parameters and their associated R² values for every individual measurement.
* **Clustering Visualization:** Utilizes `ggforce` to plot and cluster the kinetic parameters (a vs. b), allowing for rapid visual differentiation of electrochemical behavior across different experimental conditions (e.g., Positive vs. Negative potentials).
* **Consolidated Reporting:** Automatically exports a multi-sheet Excel workbook containing cleaned data, statistical summaries, and model parameters.

## Dependencies
This script requires the following R packages:
* `readxl`
* `openxlsx`
* `dplyr`
* `ggplot2`
* `ggpubr`
* `ggforce`

## Usage
1. Place your raw PalmSens export `.xlsx` file in the `data/` directory.
2. Update the `data_path` variable in the script if necessary.
3. Run the complete script to generate the console output, the visualization plots, and the `Consolidated_Results.xlsx` export.

## Author
**Felipe Longaray Kadel**
*Undergraduate Materials Engineering Student*
