# MAP outcome uncertainty data

This folder contains the prepared dataset used to evaluate the sensitivity of the final M2 malaria forecasting model to uncertainty in the Malaria Atlas Project (MAP) incidence estimates.

## File

### `malaria_M2_MAP_uncertainty.rds`

Analysis-ready dataset containing the external test observations for the three target years used in the final uncertainty analysis:

- 2009
- 2016
- 2024

The dataset contains 68,059 pixel-year observations and the following variables:

- `cell`: raster cell identifier
- `TargetYear`: malaria incidence target year
- `pred_spatial_M2`: spatial component predicted by the final M2 model
- `incidence_LCI`: MAP lower confidence limit
- `incidence_central`: MAP central incidence estimate
- `incidence_UCI`: MAP upper confidence limit

The MAP lower, central and upper incidence estimates are used to define pixel-specific triangular uncertainty distributions.

The dataset contains only the information required for the uncertainty analysis. The upstream processing of the original MAP raster products and the fitting of the M2 model are not reproduced here.