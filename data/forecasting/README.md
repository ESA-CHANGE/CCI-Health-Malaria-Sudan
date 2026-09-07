# Forecasting data

This folder contains the prepared datasets used to reproduce the final EO-based malaria forecasting model (M2) for the selected Area of Interest (AOI) in Sudan and South Sudan.

The datasets represent the final modelling support used in the analysis and do not include the original EO preprocessing workflow.

## Files

### `malaria_M2_training.rds`

Training dataset used for model fitting and internal cross-validation.

- 384,645 pixel-year observations
- Predictor years excluding 2008, 2015 and 2023
- 48 EO predictors:
  - 12 monthly Soil Moisture variables for year `t`
  - 12 monthly Precipitation variables for year `t`
  - 12 monthly Soil Moisture variables for year `t-1`
  - 12 monthly Precipitation variables for year `t-1`
- `target_spatial`: spatial component of malaria incidence in year `t+1`
- `spatial_t`: observed spatial component in year `t`, used as a persistence benchmark
- `mu_t` and `target_mu`: AOI-wide temporal components for years `t` and `t+1`
- `CV_fold`: original five-fold year-grouped cross-validation assignment

### `malaria_M2_external_test.rds`

Independent external test dataset.

The predictor years were excluded from model selection and training:

- 2008 → target year 2009
- 2015 → target year 2016
- 2023 → target year 2024

The dataset contains the same 48 EO predictors used by M2, together with the temporal and spatial components required to reconstruct total malaria incidence.

## Notes

These are analysis-ready datasets derived from the EO and epidemiological products described in the CHANGE project deliverables.

The repository starts from these prepared datasets and does not reproduce the complete upstream workflow used to download, harmonise, resample and preprocess the original Earth Observation products.