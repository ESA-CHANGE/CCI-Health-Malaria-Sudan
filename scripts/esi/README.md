# ESI analysis

This directory contains the **Environmental Suitability Index (ESI)** analysis workflow for the CHANGE UC5 Malaria case study.

The analysis compares annual ESI values with *Plasmodium falciparum* incidence estimates from the **Malaria Atlas Project (MAP)** for Sudan and South Sudan.

## Description

The workflow is implemented in:

```text
scripts/esi/01_esi_analysis.R
```

It reads annual country-level MAP–ESI summaries, defines a historical reference range, classifies subsequent annual observations as within or outside that range, produces the MAP–ESI comparison figure, and calculates agreement between the MAP and ESI classifications.

The analysis uses:

- **Reference period:** 2000–2008
- **Analysis period:** 2009–2024
- **Full plotting period:** 2000–2024

For each country, the minimum and maximum annual Q50 values observed during 2000–2008 are used to define separate reference ranges for MAP incidence and ESI.

For 2009–2024, each annual Q50 value is classified as:

- **Within range**
- **Outside range**

The classification is performed on the original Q50 values.

Standardised values are used only to place MAP and ESI on a common scale in the comparison figure. The country-specific centring and scaling parameters used for this graphical standardisation are calculated from the full annual series available in the input dataset.

## Dependencies

The script requires the following R packages:

```r
library(dplyr)
library(readr)
library(ggplot2)
library(irr)
```

Missing packages can be installed from CRAN, for example:

```r
install.packages(c("dplyr", "readr", "ggplot2", "irr"))
```

## Input data

The input dataset is expected at:

```text
data/esi/esi_map.rds
```

The script should be executed from the **repository root**, because the input path is defined relative to that directory.

The current workflow uses annual country-level variables including:

- `country`
- `year`
- `map_q05`
- `map_q50`
- `map_q95`
- `esi_q05`
- `esi_q50`
- `esi_q95`
- `map_q05_z`
- `map_q50_z`
- `map_q95_z`
- `esi_q05_z`
- `esi_q50_z`
- `esi_q95_z`

## Analysis workflow

The script performs the following steps.

### 1. Read the annual MAP–ESI dataset

The prepared annual dataset is loaded from:

```text
data/esi/esi_map.rds
```

### 2. Define analysis parameters

The temporal parameters are:

```r
year_min <- 2000L
year_max <- 2024L
reference_start <- 2000L
reference_end <- 2008L
analysis_start <- 2009L
analysis_end <- 2024L
```

### 3. Calculate graphical standardisation parameters

For each country, the script calculates the mean and standard deviation of the annual MAP Q50 and ESI Q50 series.

These parameters are used to express the historical reference-range limits on the same standardised scale used by the figure.

Standardisation is used for **graphical comparison only** and does not determine the within-range / outside-range classification.

### 4. Define the historical reference ranges

For each country and variable, the script calculates the minimum and maximum annual Q50 values observed during 2000–2008:

- MAP Q50 minimum and maximum
- ESI Q50 minimum and maximum

The reference limits are retained on the original scale for classification and are also transformed to the standardised plotting scale.

### 5. Classify annual observations

For each year from 2009 to 2024, MAP and ESI are evaluated independently against their corresponding country-specific historical ranges.

This produces two logical indicators:

```text
map_within_range
esi_within_range
```

Annual observations are subsequently labelled as `Within range` or `Outside range`.

### 6. Produce the MAP–ESI comparison figure

The figure compares annual MAP incidence and ESI for Sudan and South Sudan.

It includes:

- annual MAP Q50 trajectories;
- annual ESI Q50 trajectories;
- MAP Q05–Q95 spatial distributions;
- ESI Q05–Q95 spatial distributions;
- variable-specific 2000–2008 minimum–maximum reference limits;
- annual Q50 points from 2009 onward;
- green points for observations within the corresponding historical range;
- red points for observations outside the corresponding historical range.

The y-axis shows **standardised annual values** to allow MAP and ESI to be displayed on a common graphical scale.

The current script displays the figure on screen. The corresponding PNG output is intended to be stored under:

```text
results/esi/
```

Figure captions and descriptions of generated result files are documented in the results directory rather than in the analysis script.

### 7. Calculate concordance and Cohen's kappa

For each country, the 2009–2024 MAP and ESI classifications are converted to two categories:

```text
IN
OUT
```

The script then calculates:

- number of paired annual observations (`N`);
- number of concordant annual classifications (`agreement_n`);
- observed agreement;
- observed agreement percentage;
- unweighted Cohen's kappa;
- p-value associated with Cohen's kappa.

Cohen's kappa is calculated with:

```r
irr::kappa2()
```

The calculation is wrapped in `tryCatch()` so that a kappa value that cannot be estimated is returned as `NA` rather than stopping the workflow.

## Running the analysis

From the repository root, run:

```bash
Rscript scripts/esi/01_esi_analysis.R
```

At the current stage, the script:

1. displays the MAP–ESI comparison figure;
2. calculates the country-level concordance statistics.

## Outputs

Results from this workflow belong under:

```text
results/esi/
```

The results directory is intended to contain the exported comparison figure and any tabular outputs produced by the ESI analysis.

The full figure caption and result-specific documentation should be maintained in:

```text
results/esi/README.md
```

## Repository context

The ESI workflow is one of the three main analytical components of the UC5 Malaria repository:

```text
scripts/esi/
scripts/modelling/
scripts/uncertainty/
```

The root `README.md` provides the general overview of the complete use case.
