# RRSARMMI: Ridge Regularization for Spatial Auto-Regressive Models with Multicollinearity Issues

This repository contains the code and data necessary to reproduce the results presented in the paper "Ridge Regularization for Spatial Auto-regressive Models with Multicollinearity Issues" submitted to Advances in Statistical Analysis (AStA). 

Authors: Cristina Chavez-Chong, Cécile Hardouin, Ana Karina Fermin

## Table of Contents

- [Introduction](#introduction)
- [Usage](#usage)
- [Folder Structure](#folder-structure)
- [Data](#data)
- [Scripts](#scripts)
- [Reports](#reports)
- [License](#license)
- [Contact](#contact)

## Introduction

In this repository, we provide the resources necessary to reproduce the analysis and results presented in our paper "Ridge Regularization for Spatial Auto-regressive Models with Multicollinearity Issues" submitted to Advances in Statistical Analysis (AStA). The project provides two main resources:

- Two examples of our numerical experiments (RRSAR and RRSEM algorithms).
- The reproduction of the analysis of a dataset about COVID-19 in metropolitan France.

## Usage

The code in this repository was built in R version 3.6.2.  To use the provided scripts and functions, follow these steps:

1. Clone this GitHub repository to your local machine using the following command:
   ```bash
   git clone https://github.com/c0ra/RRSAR.git
   ```
2. Install Required Packages:
Navigate to the project directory and run the `requirements.R` script to install the required R packages. Use the following command in your R console or terminal:
   ```R
   source("requirements.R")
   ```
   This will install all the necessary packages automatically.
3. Explore the Scripts and Data:
Explore the contents of the repository, including R Markdown files, scripts, and data files, to reproduce the analyses or adapt them to your needs.
4. Run the R Markdown Reports:
Open them in RStudio or any R Markdown compatible environment and follow the instructions within each report.

## Folder Structure

The repository is organized as follows:

- `/data`: Contains the real dataset COVID-19.
- `/scripts`: Includes R scripts with functions to perform SFRR, RRSAR and RRSEM estimation methods.
- `/reports`: Contains R Markdown files that generate analysis reports and the .pdf files with full results of our numerical experiments.

## Data

The  `covid_data.shp` is a simple features file which contains 12 variables and a geometry column type Multi Polygon. The observations are collected from the 96 French *départements*, which are administrative units. The full meaning of the variables is in the appendix of our paper.

The variables are:

* `Department`: name of the French department.

* `Region`: name of the French region.

* `LnHosp`: logarithm of hospitalization rate.

* `LnPop` : logarithm of population density.

* `SLivR` : standard of living gap.

* `C3` : percentage of inhabitants in the three largest cities of the department.

* `Work` : rate of workers.

* `Inac` : inactive rate.

* `A65pls` : share of people aged 65 and over.

* `Emer` : rate of emergency services per 1000 habitants.

* `FDoc` : rate of doctors per 1000 habitants.

* `Urban`: indicator of an urban department (vs rural).

## Scripts

The repository includes three scripts:

- `requirements.R`: Includes the preliminary code needed to install the necessary libraries.
- `functions_numerical_experiments.R`: Includes the code for SFRR, RRSAR and RRSEM estimation procedures for the simulation study.
- `functions_covid.R`: Includes the code for SFRR, RRSAR and RRSEM estimation procedures in the context of the application.

The `functions_*` scripts are called in the main .Rmd files to run the SFRR, RRSAR and RRSEM estimation algorithms.

## Reports

The reports folder contains .html files for generated in Rmarkdown and .pdf files with the full results that are summarised in our paper. 

- `covid_19_application.Rmd`: Instructions to replicate the exploratory analysis performed on the COVID-19 dataset, as well as to apply RRSAR and RRSEM methodology to estimate all the parameters and to assess the importance of the covariates on the COVID-19 intensity.
- `example_simulation_estimation_sar.Rmd`: Instructions to generate one simulation following the SAR model structure and to estimate the model coefficients using the OLS, SAR, Ridge regression, SFRR and RRSAR procedures.
- `example_simulation_estimation_sar.html`: Report with instructions to generate one simulation following the SAR model structure and to estimate the model coefficients using the OLS, SAR, Ridge regression, SFRR and RRSAR procedures.
- - `example_simulation_estimation_sem.Rmd`: Instructions to generate one simulation following the SEM model structure and to estimate the model coefficients using the OLS, SEM, Ridge regression, SFRR and RRSAR procedures.
- `example_simulation_estimation_sem.html`: Report with instructions to generate one simulation following the SEM model structure and to estimate the model coefficients using the OLS, SEM, Ridge regression, SFRR and RRSAR procedures.
- `complete_results_sar.pdf`: Full results of our simulation study with data generated following the SAR model.
- `complete_results_sem.pdf`: Full results of our simulation study with data generated following the SEM model.

To reproduce specific results, run the corresponding `.Rmd` file. Make sure to adjust file paths and settings as needed.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

If you have any questions or need further information, please contact [Cristina Chavez](mailto:cristi0929@gmail.com).
