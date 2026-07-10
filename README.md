# Return values of significant wave heights and near-bed orbital velocity in strong currents : an application to the Alderney Race using the Resourcecode hindcast (1994-2024)

This repository contains the R scripts used to analyse sea states at the Alderney Race site. The scripts cover exploratory data description, stationary extreme-value analysis, and non-stationary EVGAM models with DIFORM contours in two and three dimensions. Figures and intermediate results are saved in `figures/` and `outputs/`.

The main dataset is `data_alderney_rscd_2024.csv`. Paths are managed through `paths.R`; an alternative data location can be provided with the `DATA_ALDERNEY_RSCD_2024_PATH` environment variable. The main pipelines can be run from the repository root with `source("nonstationary_EVA/04c_run_pipeline_DIFORM_3D.R")` and `source("nonstationary_EVA/05b_run_pipeline_DIFORM_2D.R")`.
