# MAP outcome uncertainty analysis

This folder contains the script used to propagate uncertainty in MAP malaria incidence estimates through the external evaluation of the final M2 model.

## Script

### `03_map_outcome_uncertainty.R`

The script evaluates how uncertainty in the MAP incidence target affects the apparent predictive performance of M2.

M2 predictions are kept fixed throughout the analysis. Alternative target surfaces are generated from the MAP uncertainty bounds.

For each pixel, a triangular distribution is defined on the `log1p` incidence scale using:

- MAP lower confidence limit as the minimum
- MAP central estimate as the mode
- MAP upper confidence limit as the maximum

To preserve spatial coherence, one random quantile is generated for each target year and Monte Carlo simulation and applied to all pixel-specific triangular distributions within that year.

Each simulated incidence surface is then re-centred around its AOI-wide mean to obtain the corresponding spatial anomaly.

Model performance is recalculated for each simulation using:

- RMSE
- MAE
- predictive R²
- bias

A total of 10,000 Monte Carlo simulations are performed for each of the three target years:

- 2009
- 2016
- 2024

## Requirements

The analysis uses the following R packages:

- `dplyr`
- `foreach`
- `doParallel`
- `ggplot2`

The script should be run from the root directory of the repository.

## Input

- `data/uncertainty/malaria_M2_MAP_uncertainty.rds`

## Outputs

Results are written to:

`results/uncertainty/`