# AOI malaria forecasting

This folder contains the script used to reproduce the final EO-based malaria forecasting analysis for the selected Area of Interest (AOI) in Sudan and South Sudan.

## Script

### `02_aoi_malaria_forecasting.R`

The script implements the final M2 model described in the CHANGE project analysis.

M2 predicts the relative spatial distribution of malaria incidence in year `t+1` from environmental conditions observed in years `t` and `t-1`:

`spatial(t+1) = f[EO(t), EO(t-1)]`

The EO predictors are:

- monthly Soil Moisture
- monthly Precipitation

for both the current and previous year, resulting in 48 predictor variables.

The AOI-wide temporal component is forecast using persistence:

`mu(t+1) = mu(t)`

The predicted temporal and spatial components are then combined to reconstruct malaria incidence.

## Validation

Internal model performance is assessed using the original five-fold cross-validation structure, with complete years assigned to individual folds.

Three predictor years are excluded from model fitting and used for independent external testing:

- 2008 → 2009
- 2015 → 2016
- 2023 → 2024

The script also compares M2 with:

- spatial persistence
- AOI-wide temporal persistence without spatial information

## Requirements

The analysis uses the following R packages:

- `dplyr`
- `caret`
- `ranger`
- `doParallel`

The script should be run from the root directory of the repository.

## Inputs

- `data/forecasting/malaria_M2_training.rds`
- `data/forecasting/malaria_M2_external_test.rds`

## Outputs

Results are written to:

`results/forecasting/`