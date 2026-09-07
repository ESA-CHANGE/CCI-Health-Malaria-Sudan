# ============================================================
# CHANGE - Malaria Use Case
# MAP outcome uncertainty propagation for the final M2 model
#
# The analysis evaluates how uncertainty in MAP malaria
# incidence estimates affects the apparent predictive
# performance of the final EO-based M2 spatial model.
#
# M2 predictions are kept fixed.
#
# For each target year (2009, 2016, 2024), alternative MAP
# incidence surfaces are generated from pixel-specific
# triangular distributions defined by:
#
#   minimum = MAP lower confidence limit (LCI)
#   mode    = MAP central estimate
#   maximum = MAP upper confidence limit (UCI)
#
# A single random quantile is shared by all pixels within
# each target year and simulation, preserving spatial
# coherence while allowing pixel-specific uncertainty ranges.
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES AND SETTINGS
# ------------------------------------------------------------

library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)

set.seed(123)

data.dir <- file.path(
  "data",
  "uncertainty"
)

results.dir <- file.path(
  "results",
  "uncertainty"
)

dir.create(
  results.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Number of Monte Carlo simulations
B <- 10000L


# ------------------------------------------------------------
# 2. LOAD PREPARED DATA
# ------------------------------------------------------------

uncertainty.data <- readRDS(
  file.path(
    data.dir,
    "malaria_M2_MAP_uncertainty.rds"
  )
)

required.variables <- c(
  "cell",
  "TargetYear",
  "pred_spatial_M2",
  "incidence_LCI",
  "incidence_central",
  "incidence_UCI"
)

stopifnot(
  all(
    required.variables %in%
      names(uncertainty.data)
  )
)

stopifnot(
  nrow(uncertainty.data) == 68059
)

stopifnot(
  identical(
    sort(
      unique(
        uncertainty.data$TargetYear
      )
    ),
    c(
      2009L,
      2016L,
      2024L
    )
  )
)


# ------------------------------------------------------------
# 3. CHECK MAP UNCERTAINTY BOUNDS
# ------------------------------------------------------------

bounds.check <- uncertainty.data |>
  summarise(
    N = n(),
    
    central_below_LCI = sum(
      incidence_central < incidence_LCI
    ),
    
    central_above_UCI = sum(
      incidence_central > incidence_UCI
    ),
    
    UCI_below_LCI = sum(
      incidence_UCI < incidence_LCI
    )
  )

print(
  bounds.check
)

stopifnot(
  bounds.check$central_below_LCI == 0,
  bounds.check$central_above_UCI == 0,
  bounds.check$UCI_below_LCI == 0
)


# ------------------------------------------------------------
# 4. TRANSFORM MAP INCIDENCE TO LOG SCALE
# ------------------------------------------------------------

uncertainty.data <- uncertainty.data |>
  mutate(
    
    log_LCI =
      log1p(
        incidence_LCI
      ),
    
    log_central =
      log1p(
        incidence_central
      ),
    
    log_UCI =
      log1p(
        incidence_UCI
      )
    
  ) |>
  group_by(
    TargetYear
  ) |>
  mutate(
    
    # Relative spatial component of the MAP central estimate
    spatial_central =
      log_central -
      mean(
        log_central,
        na.rm = TRUE
      )
    
  ) |>
  ungroup()


# ------------------------------------------------------------
# 5. PERFORMANCE METRIC FUNCTION
# ------------------------------------------------------------

calculate_metrics <- function(
    observed,
    predicted
) {
  
  valid <- complete.cases(
    observed,
    predicted
  )
  
  observed <- observed[valid]
  predicted <- predicted[valid]
  
  error <- predicted - observed
  
  SSE <- sum(
    error^2
  )
  
  SST <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  data.frame(
    
    N = length(observed),
    
    RMSE = sqrt(
      mean(
        error^2
      )
    ),
    
    MAE = mean(
      abs(error)
    ),
    
    R2_predictive = if (
      SST > 0
    ) {
      
      1 - SSE / SST
      
    } else {
      
      NA_real_
    },
    
    R2_correlation = if (
      sd(observed) > 0 &&
      sd(predicted) > 0
    ) {
      
      cor(
        observed,
        predicted
      )^2
      
    } else {
      
      NA_real_
    },
    
    Bias = mean(
      error
    )
  )
}


# ------------------------------------------------------------
# 6. PERFORMANCE AGAINST MAP CENTRAL ESTIMATE
# ------------------------------------------------------------

Central.performance.by.year <-
  uncertainty.data |>
  group_by(
    TargetYear
  ) |>
  group_modify(
    ~ calculate_metrics(
      observed =
        .x$spatial_central,
      
      predicted =
        .x$pred_spatial_M2
    )
  ) |>
  ungroup()


print(
  Central.performance.by.year
)


# ------------------------------------------------------------
# 7. TRIANGULAR QUANTILE FUNCTION
# ------------------------------------------------------------

qtriangular_vector <- function(
    p,
    minimum,
    mode,
    maximum
) {
  
  stopifnot(
    length(minimum) == length(mode),
    length(mode) == length(maximum)
  )
  
  n <- length(minimum)
  
  # A single probability can be applied to all pixels.
  if (length(p) == 1L) {
    p <- rep(
      p,
      n
    )
  }
  
  stopifnot(
    length(p) == n
  )
  
  p <- pmin(
    pmax(
      p,
      0
    ),
    1
  )
  
  # Numerical protection
  mode <- pmin(
    pmax(
      mode,
      minimum
    ),
    maximum
  )
  
  output <- numeric(
    n
  )
  
  interval <-
    maximum -
    minimum
  
  # Degenerate distributions
  degenerate <-
    !is.finite(interval) |
    interval <= 0
  
  output[degenerate] <-
    minimum[degenerate]
  
  valid <- !degenerate
  
  if (any(valid)) {
    
    a <- minimum[valid]
    c <- mode[valid]
    b <- maximum[valid]
    u <- p[valid]
    
    F.mode <-
      (c - a) /
      (b - a)
    
    left <-
      u < F.mode
    
    sampled <- numeric(
      length(a)
    )
    
    if (any(left)) {
      
      sampled[left] <-
        a[left] +
        sqrt(
          u[left] *
            (b[left] - a[left]) *
            (c[left] - a[left])
        )
    }
    
    if (any(!left)) {
      
      sampled[!left] <-
        b[!left] -
        sqrt(
          (1 - u[!left]) *
            (b[!left] - a[!left]) *
            (b[!left] - c[!left])
        )
    }
    
    output[valid] <-
      sampled
  }
  
  output
}


# ------------------------------------------------------------
# 8. PREPARE YEAR-SPECIFIC INPUTS
# ------------------------------------------------------------

target.years <- sort(
  unique(
    uncertainty.data$TargetYear
  )
)

year.data <- lapply(
  target.years,
  function(yr) {
    
    d <- uncertainty.data |>
      filter(
        TargetYear == yr
      )
    
    list(
      
      TargetYear = yr,
      
      log_LCI =
        d$log_LCI,
      
      log_central =
        d$log_central,
      
      log_UCI =
        d$log_UCI,
      
      prediction =
        d$pred_spatial_M2
    )
  }
)

names(year.data) <-
  as.character(
    target.years
  )


# ------------------------------------------------------------
# 9. PRE-GENERATE MONTE CARLO QUANTILES
# ------------------------------------------------------------

# One independent quantile per target year and simulation.
# Within a year, the same quantile is applied to every pixel.

set.seed(123)

quantile.matrix <- matrix(
  runif(
    B *
      length(target.years)
  ),
  nrow = B,
  ncol = length(target.years)
)

colnames(
  quantile.matrix
) <- target.years


# ------------------------------------------------------------
# 10. FAST METRIC FUNCTION FOR MONTE CARLO
# ------------------------------------------------------------

calculate_metrics_fast <- function(
    observed,
    predicted
) {
  
  error <-
    predicted -
    observed
  
  SSE <- sum(
    error^2
  )
  
  SST <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  c(
    
    N =
      length(observed),
    
    RMSE =
      sqrt(
        mean(
          error^2
        )
      ),
    
    MAE =
      mean(
        abs(error)
      ),
    
    R2_predictive =
      if (
        SST > 0
      ) {
        
        1 - SSE / SST
        
      } else {
        
        NA_real_
      },
    
    Bias =
      mean(error)
  )
}


# ------------------------------------------------------------
# 11. ONE MONTE CARLO SIMULATION
# ------------------------------------------------------------

run_one_simulation <- function(
    simulation.id
) {
  
  yearly.results <- vector(
    "list",
    length(year.data)
  )
  
  for (
    j in seq_along(year.data)
  ) {
    
    d <- year.data[[j]]
    
    # Spatially coherent annual quantile
    u.year <-
      quantile.matrix[
        simulation.id,
        j
      ]
    
    simulated.log.incidence <-
      qtriangular_vector(
        
        p =
          u.year,
        
        minimum =
          d$log_LCI,
        
        mode =
          d$log_central,
        
        maximum =
          d$log_UCI
      )
    
    # Convert simulated incidence surface to spatial anomaly
    simulated.spatial <-
      simulated.log.incidence -
      mean(
        simulated.log.incidence
      )
    
    metrics <-
      calculate_metrics_fast(
        
        observed =
          simulated.spatial,
        
        predicted =
          d$prediction
      )
    
    yearly.results[[j]] <-
      data.frame(
        
        Simulation =
          simulation.id,
        
        TargetYear =
          d$TargetYear,
        
        Quantile =
          u.year,
        
        Target_spatial_SD =
          sd(
            simulated.spatial
          ),
        
        N =
          unname(
            metrics["N"]
          ),
        
        RMSE =
          unname(
            metrics["RMSE"]
          ),
        
        MAE =
          unname(
            metrics["MAE"]
          ),
        
        R2_predictive =
          unname(
            metrics["R2_predictive"]
          ),
        
        Bias =
          unname(
            metrics["Bias"]
          )
      )
  }
  
  do.call(
    rbind,
    yearly.results
  )
}


# ------------------------------------------------------------
# 12. RUN MONTE CARLO ANALYSIS
# ------------------------------------------------------------

available.cores <-
  parallel::detectCores()

if (
  is.na(available.cores)
) {
  available.cores <- 2L
}

n.cores <- max(
  1L,
  min(
    12L,
    available.cores - 1L
  )
)

cl <- parallel::makePSOCKcluster(
  n.cores
)

registerDoParallel(
  cl
)

simulation.results <- tryCatch(
  
  foreach(
    b = seq_len(B),
    .packages = character(0),
    .inorder = TRUE
  ) %dopar% {
    
    run_one_simulation(b)
  },
  
  finally = {
    
    parallel::stopCluster(cl)
    registerDoSEQ()
  }
)


MC.performance <-
  bind_rows(
    simulation.results
  )


# ------------------------------------------------------------
# 13. SUMMARISE MONTE CARLO RESULTS
# ------------------------------------------------------------

Final.summary <-
  MC.performance |>
  group_by(
    TargetYear
  ) |>
  summarise(
    
    RMSE_median =
      median(
        RMSE,
        na.rm = TRUE
      ),
    
    RMSE_low =
      quantile(
        RMSE,
        0.025,
        na.rm = TRUE
      ),
    
    RMSE_high =
      quantile(
        RMSE,
        0.975,
        na.rm = TRUE
      ),
    
    MAE_median =
      median(
        MAE,
        na.rm = TRUE
      ),
    
    MAE_low =
      quantile(
        MAE,
        0.025,
        na.rm = TRUE
      ),
    
    MAE_high =
      quantile(
        MAE,
        0.975,
        na.rm = TRUE
      ),
    
    R2_median =
      median(
        R2_predictive,
        na.rm = TRUE
      ),
    
    R2_low =
      quantile(
        R2_predictive,
        0.025,
        na.rm = TRUE
      ),
    
    R2_high =
      quantile(
        R2_predictive,
        0.975,
        na.rm = TRUE
      ),
    
    Negative_R2_percent =
      100 *
      mean(
        R2_predictive < 0,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  left_join(
    
    Central.performance.by.year |>
      transmute(
        
        TargetYear,
        
        RMSE_central =
          RMSE,
        
        MAE_central =
          MAE,
        
        R2_central =
          R2_predictive,
        
        Bias_central =
          Bias
      ),
    
    by = "TargetYear"
  )


print(
  Final.summary
)


# ------------------------------------------------------------
# 14. WRITE RESULTS
# ------------------------------------------------------------

write.csv(
  Final.summary,
  file.path(
    results.dir,
    "map_uncertainty_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  MC.performance,
  file.path(
    results.dir,
    "map_uncertainty_simulations.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. FINAL FIGURE
# ------------------------------------------------------------

Plot.data <- bind_rows(
  
  Final.summary |>
    transmute(
      TargetYear,
      Metric = "RMSE",
      Central = RMSE_central,
      Median = RMSE_median,
      Lower = RMSE_low,
      Upper = RMSE_high
    ),
  
  Final.summary |>
    transmute(
      TargetYear,
      Metric = "MAE",
      Central = MAE_central,
      Median = MAE_median,
      Lower = MAE_low,
      Upper = MAE_high
    ),
  
  Final.summary |>
    transmute(
      TargetYear,
      Metric = "Predictive R2",
      Central = R2_central,
      Median = R2_median,
      Lower = R2_low,
      Upper = R2_high
    )
)

Plot.data$Metric <- factor(
  Plot.data$Metric,
  levels = c(
    "RMSE",
    "MAE",
    "Predictive R2"
  )
)


Final.plot <- ggplot(
  Plot.data,
  aes(
    x = factor(TargetYear)
  )
) +
  geom_errorbar(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    width = 0.08,
    linewidth = 0.8,
    position =
      position_nudge(
        x = -0.07
      )
  ) +
  geom_point(
    aes(
      y = Median,
      shape = "Monte Carlo median"
    ),
    size = 3,
    position =
      position_nudge(
        x = -0.07
      )
  ) +
  geom_point(
    aes(
      y = Central,
      shape = "MAP central estimate"
    ),
    size = 3.3,
    position =
      position_nudge(
        x = 0.07
      )
  ) +
  facet_wrap(
    ~Metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_shape_manual(
    values = c(
      "Monte Carlo median" = 16,
      "MAP central estimate" = 18
    ),
    name = NULL
  ) +
  labs(
    title =
      "Sensitivity of M2 performance to MAP incidence uncertainty",
    
    subtitle =
      paste0(
        format(
          B,
          big.mark = ","
        ),
        " spatially coherent Monte Carlo simulations"
      ),
    
    x =
      "Target year",
    
    y =
      NULL,
    
    caption =
      paste(
        "Bars represent the 2.5th-97.5th percentile range.",
        "A common random quantile was applied to the",
        "pixel-specific triangular distributions within each year."
      )
  ) +
  theme_bw() +
  theme(
    legend.position =
      "bottom",
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    panel.grid.minor =
      element_blank()
  )


# ggsave(
#   filename = file.path(
#     results.dir,
#     "map_uncertainty_sensitivity.png"
#   ),
#   plot = Final.plot,
#   width = 10,
#   height = 4.5,
#   dpi = 300
# )


# ------------------------------------------------------------
# 16. FINAL SUMMARY
# ------------------------------------------------------------

cat(
  "\n=============================================\n",
  "MAP outcome uncertainty analysis complete\n",
  "=============================================\n\n",
  "Monte Carlo simulations: ",
  B,
  "\n\n",
  "Results written to:\n",
  results.dir,
  "\n",
  sep = ""
)
