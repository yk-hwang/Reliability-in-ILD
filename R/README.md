# R

A single R script (datagen_combine_analysis.R) was developed to manipulate and verify the simulation conditions related to within-person reliability levels.

## Context & Citation

This code was used for the simulation study in the following manuscript:

[Manuscript in preparation]
Hwang, Y. K., Lee, Y. S. & Suk, H. W. (2025). Examining the impact of within-person reliability in dynamic structural equation modeling.

## Core Script

All functionalities are contained within the single R script: datagen_combine_analysis.R.

This script manages the entire workflow, from data generation to model fitting and parameter verification.

## Requirements

R (tested on version 4.5.0)

Mplus (tested on version 8.7)

R Packages:

MplusAutomation: To interface between R and Mplus.

(other packages, e.g., dplyr, stringr, purrr, readr, openxlsx)

## Workflow

The datagen_combine_analysis.R script executes the following sequential steps:

1. Data Generation: Calls Mplus to generate Intensive Longitudinal Data (ILD) with varying, pre-defined levels of within-person reliability.
2. Replication Setup: Collects and organizes the generated data to create 100 replication datasets.
3. Model Fitting: Uses MplusAutomation to automatically run a MEAR(1) model (Multilevel Error Autoregressive Model) on each of the 100 replication datasets in Mplus.
4. Parameter Extraction: Reads the Mplus output files (.out) back into R to extract the estimated person-specific parameters for each replication. This includes: Innovation variance, Measurement error variance, Autoregressive coefficient
5. Condition Verification: Calculates the mean and variance of the estimated within-person reliability parameters across replications. This final step verifies that the generated data successfully reflects the researcher-intended distribution of within-person reliability.

## Author

Y. K. Hwang

Contact: hwangyk96@gmail.com
