# Statistical modeling of significant wave height and near-bed orbital speed in strong currents: an application to the Alderney Race using the ResourceCode hindcast (1994–2026)

R code for the analyses and figures of the paper of the same title. The scripts cover the description of the site, a stationary peaks-over-threshold analysis, a non-stationary GAM–GPD model fitted with `evgam`, and three-dimensional Direct-IFORM (D-IFORM) environmental contours.

## Data

All analyses use `data_rscd_alderney_2026.csv`, hourly sea states at ResourceCode node 158348. The file contains, among others, significant wave height `hs`, peak period `tp`, peak direction `dp`, current components `ucur`/`vcur`, current speed `Cspd`, wave–current angle `angleDiff`, near-bed orbital speed `Os` (at 10 m above the seabed) and the tidal phase `fje`: slack (`étale`) when the current speed is below 1 m/s, otherwise ebb (`jusant`) when the current flows towards 180–360° and flood (`flot`) when it flows towards 0–180°. Only the common window from 1994-04-01 00:00 UTC to 2026-03-31 23:00 UTC is used; the CSV itself is not modified.

`importation/01_data_importation_example.R` shows how the hindcast variables and the orbital speed are extracted with the `resourcecode` package. The ResourceCode database is publicly available at https://resourcecode.ifremer.fr/.

Paths are handled by `paths.R`. The data file can be placed elsewhere and given through the `DATA_RSCD_2026_PATH` environment variable.

## Repository structure

| Folder | Scripts | Content |
|---|---|---|
| `importation/` | `01_data_importation_example.R`, `02_read_rscd_2026.R` | Data extraction example; reading and preparation of the CSV |
| `site_description/` | `00_location_map.R` to `04_time_series.R` | Location map, distributions, directions and tides, density plots, severe storms, time series |
| `stationary_EVA/` | `00_stationary_helpers.R` to `05_run_stationary_EVA.R` | Threshold selection, declustering sensitivity, stationary GPD fit and return levels |
| `nonstationary_EVA/` | `01_select_evgam.R`, `02_run_evgam_analysis.R`, `03_evgam_plots.R` | Covariate-dependent ALD threshold, GAM–GPD model, block cross-validation, partial effects and QQ plots |
| `contours/` | `00_diform_explanation.R` to `04_run_diform_3d.R`, `pot_decluster.cpp` | D-IFORM illustration in (Hs, Tp), three-dimensional contours and directional QQ diagnostics |

Figures are written to `figures/`, fitted objects to `outputs/`.

## Requirements

The code was run with R 4.5. Required packages: `data.table`, `dplyr`, `ggplot2`, `patchwork`, `cowplot`, `gridExtra`, `scales`, `viridis`, `viridisLite`, `ragg`, `sf`, `rnaturalearth`, `rnaturalearthdata`, `resourcecode`, `resourcecodedata`, `evd`, `extRemes`, `POT`, `evgam` and `Rcpp`. The declustering of the D-IFORM projections is written in C++ (`contours/pot_decluster.cpp`) and compiled by `Rcpp`, which requires a C++ compiler (Rtools on Windows).

## Reproducing the results

Run the scripts from the repository root.

```r
# Site description (Figs. 1, 3-5, 7-9)
source("site_description/00_location_map.R")
source("site_description/01_histograms_roses.R")
source("site_description/02_scatterplots.R")
source("site_description/03_storms.R")
source("site_description/04_time_series.R")

# Stationary analysis (Figs. 10-13): set response_variable to "hs" or "Os" at the top of the script
source("stationary_EVA/05_run_stationary_EVA.R")

# Non-stationary GAM-GPD model (Figs. 14-19)
source("nonstationary_EVA/02_run_evgam_analysis.R")
source("nonstationary_EVA/03_evgam_plots.R")

# D-IFORM contours (Figs. 20-24): set response_variable to "hs" or "Os" at the top of 04_run_diform_3d.R
source("contours/00_diform_explanation.R")
source("contours/04_run_diform_3d.R")
```

`02_run_evgam_analysis.R` runs the block cross-validation for both responses and saves the results to `outputs/evgam_selection_rscd_2026_hs.rds` and `outputs/evgam_selection_rscd_2026_os.rds` (about 15 min on 6 cores). Setting `RUN_MODEL_SELECTION <- FALSE` reloads these files instead of refitting.
