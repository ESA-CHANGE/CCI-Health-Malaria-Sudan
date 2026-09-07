# ============================================================
# CHANGE - Malaria Use Case
# AOI malaria forecasting - Final M2 model
#
# M2 predicts the spatial component of malaria incidence at
# year t+1 using current and one-year-lagged Earth Observation
# variables:
#
#   s[p,t+1] = f(EO[p,t], EO[p,t-1])
#
# EO predictors:
#   - Soil Moisture (12 monthly variables)
#   - Precipitation (12 monthly variables)
#
# The AOI-wide temporal component is predicted by persistence:
#
#   mu[t+1] = mu[t]
#
# External predictor years:
#   2008, 2015, 2023
#
# Corresponding target years:
#   2009, 2016, 2024
# ============================================================


# ------------------------------------------------------------
# 1. PACKAGES AND SETTINGS
# ------------------------------------------------------------

library(dplyr)
library(caret)
library(ranger)
library(doParallel)

set.seed(123)

data.dir <- file.path(
  "data",
  "forecasting"
)

results.dir <- file.path(
  "results",
  "forecasting"
)

dir.create(
  results.dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. LOAD PREPARED DATA
# ------------------------------------------------------------

train.data <- readRDS(
  file.path(
    data.dir,
    "malaria_M2_training.rds"
  )
)

test.data <- readRDS(
  file.path(
    data.dir,
    "malaria_M2_external_test.rds"
  )
)


cat(
  "Training observations:",
  nrow(train.data),
  "\n"
)

cat(
  "External test observations:",
  nrow(test.data),
  "\n"
)


# Expected dataset sizes
stopifnot(
  nrow(train.data) == 384645,
  nrow(test.data) == 68059
)


# ------------------------------------------------------------
# 3. DEFINE M2 PREDICTORS
# ------------------------------------------------------------

SM.current <- grep(
  "^SM_m[0-9]{2}$",
  names(train.data),
  value = TRUE
)

PRE.current <- grep(
  "^PRE_m[0-9]{2}$",
  names(train.data),
  value = TRUE
)

SM.lag1 <- grep(
  "^SM_m[0-9]{2}_lag1$",
  names(train.data),
  value = TRUE
)

PRE.lag1 <- grep(
  "^PRE_m[0-9]{2}_lag1$",
  names(train.data),
  value = TRUE
)


M2.variables <- c(
  SM.current,
  PRE.current,
  SM.lag1,
  PRE.lag1
)

target.spatial <- "target_spatial"


# M2 must contain 48 EO predictors:
# 12 SM(t) + 12 PRE(t) + 12 SM(t-1) + 12 PRE(t-1)

stopifnot(
  length(SM.current) == 12,
  length(PRE.current) == 12,
  length(SM.lag1) == 12,
  length(PRE.lag1) == 12,
  length(M2.variables) == 48
)

stopifnot(
  all(
    c(
      M2.variables,
      target.spatial,
      "CV_fold"
    ) %in% names(train.data)
  )
)


# ------------------------------------------------------------
# 4. RECONSTRUCT THE ORIGINAL CROSS-VALIDATION FOLDS
# ------------------------------------------------------------

# CV_fold contains the original validation-fold assignment used
# during the final model comparison. This avoids regenerating
# the folds and ensures that the same year-grouped CV structure
# is reproduced.

fold.ids <- sort(
  unique(train.data$CV_fold)
)

stopifnot(
  identical(
    fold.ids,
    1:5
  )
)


cv.index.out <- lapply(
  fold.ids,
  function(f) {
    which(
      train.data$CV_fold == f
    )
  }
)

cv.index <- lapply(
  cv.index.out,
  function(validation.rows) {
    setdiff(
      seq_len(nrow(train.data)),
      validation.rows
    )
  }
)

names(cv.index) <- paste0(
  "Fold",
  fold.ids
)

names(cv.index.out) <- names(cv.index)


# Verify that whole years are kept together within each fold.

for (i in seq_along(cv.index)) {
  
  training.years <- unique(
    train.data$Year[
      cv.index[[i]]
    ]
  )
  
  validation.years <- unique(
    train.data$Year[
      cv.index.out[[i]]
    ]
  )
  
  stopifnot(
    length(
      intersect(
        training.years,
        validation.years
      )
    ) == 0
  )
  
  cat(
    names(cv.index)[i],
    "\n  Training years:",
    paste(sort(training.years), collapse = ", "),
    "\n  Validation years:",
    paste(sort(validation.years), collapse = ", "),
    "\n\n"
  )
}


# ------------------------------------------------------------
# 5. RANDOM FOREST TUNING GRID
# ------------------------------------------------------------

# Same tuning strategy used for the final model comparison.

make_tune_grid <- function(p) {
  
  mtry.values <- unique(
    pmax(
      1L,
      pmin(
        p,
        as.integer(
          round(
            c(
              p / 8,
              sqrt(p),
              p / 4,
              p / 2
            )
          )
        )
      )
    )
  )
  
  expand.grid(
    mtry = mtry.values,
    splitrule = c(
      "variance",
      "extratrees"
    ),
    min.node.size = 20L,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}


tune.grid <- make_tune_grid(
  length(M2.variables)
)

print(tune.grid)


# ------------------------------------------------------------
# 6. CROSS-VALIDATION CONTROL
# ------------------------------------------------------------

# Seed 102 corresponds to M2 in the original comparative
# modelling workflow.

set.seed(102)

cv.seeds <- vector(
  "list",
  length(cv.index) + 1L
)

for (i in seq_along(cv.index)) {
  
  cv.seeds[[i]] <- sample.int(
    1000000,
    nrow(tune.grid)
  )
}

cv.seeds[[length(cv.index) + 1L]] <-
  sample.int(
    1000000,
    1
  )


train.control <- trainControl(
  method = "cv",
  index = cv.index,
  indexOut = cv.index.out,
  savePredictions = "final",
  returnResamp = "final",
  returnData = FALSE,
  allowParallel = TRUE,
  verboseIter = TRUE,
  seeds = cv.seeds
)


# ------------------------------------------------------------
# 7. FIT M2
# ------------------------------------------------------------

available.cores <- parallel::detectCores()

if (is.na(available.cores)) {
  available.cores <- 2L
}

n.cores <- max(
  1L,
  min(
    12L,
    available.cores - 1L
  )
)

cl <- makePSOCKcluster(
  n.cores
)

registerDoParallel(cl)


set.seed(102)

rf.M2 <- tryCatch(
  
  train(
    x = train.data[
      ,
      M2.variables,
      drop = FALSE
    ],
    
    y = train.data[,
      target.spatial,
      drop = TRUE
    ],
    
    method = "ranger",
    
    trControl = train.control,
    
    tuneGrid = tune.grid,
    
    metric = "RMSE",
    
    num.trees = 300,
    
    importance = "permutation",
    
    num.threads = 2
  ),
  
  finally = {
    
    stopCluster(cl)
    registerDoSEQ()
    
  }
)


cat(
  "\nBest tuning parameters:\n"
)

print(
  rf.M2$bestTune
)

 
# ------------------------------------------------------------
# 8. EXTRACT OUT-OF-FOLD PREDICTIONS
# ------------------------------------------------------------

extract_oof_predictions <- function(
    fit,
    n
) {
  
  predictions <- fit$pred
  
  # Retain predictions corresponding to the selected
  # hyperparameter combination.
  
  for (v in names(fit$bestTune)) {
    
    predictions <- predictions[
      predictions[[v]] ==
        fit$bestTune[[v]],
      ,
      drop = FALSE
    ]
  }
  
  stopifnot(
    !anyDuplicated(
      predictions$rowIndex
    )
  )
  
  output <- rep(
    NA_real_,
    n
  )
  
  output[
    predictions$rowIndex
  ] <- predictions$pred
  
  output
}


train.data$pred_spatial_M2_oof <-
  extract_oof_predictions(
    rf.M2,
    nrow(train.data)
  )


# ------------------------------------------------------------
# 9. PERFORMANCE METRICS
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
  
  error <- observed - predicted
  
  sse <- sum(
    error^2
  )
  
  sst <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  r2.predictive <- if (
    sst > 0
  ) {
    
    1 - sse / sst
    
  } else {
    
    NA_real_
  }
  
  
  r2.correlation <- if (
    sd(observed) > 0 &&
    sd(predicted) > 0
  ) {
    
    cor(
      observed,
      predicted
    )^2
    
  } else {
    
    NA_real_
  }
  
  
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
    
    R2_predictive = r2.predictive,
    
    R2_correlation = r2.correlation,
    
    Bias = mean(
      predicted - observed
    )
  )
}


# ------------------------------------------------------------
# 10. CROSS-VALIDATED M2 PERFORMANCE
# ------------------------------------------------------------

CV.performance <- bind_rows(
  
  cbind(
    ModelID = "AOI_mean_only",
    Model = "AOI mean only",
    
    calculate_metrics(
      observed =
        train.data$target_spatial,
      
      predicted =
        rep(
          0,
          nrow(train.data)
        )
    )
  ),
  
  cbind(
    ModelID = "Spatial_persistence",
    Model = "Spatial persistence",
    
    calculate_metrics(
      observed =
        train.data$target_spatial,
      
      predicted =
        train.data$spatial_t
    )
  ),
  
  cbind(
    ModelID = "M2",
    Model = "M2: EO_t + EO_t-1",
    
    calculate_metrics(
      observed =
        train.data$target_spatial,
      
      predicted =
        train.data$pred_spatial_M2_oof
    )
  )
) |>
  arrange(RMSE)


print(
  CV.performance
)


write.csv(
  CV.performance,
  file.path(
    results.dir,
    "m2_cross_validation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. EXTERNAL TEST
# ------------------------------------------------------------

# Predictor years 2008, 2015 and 2023 were kept completely
# outside model selection and correspond to target years
# 2009, 2016 and 2024.

stopifnot(
  identical(
    sort(
      unique(
        test.data$Year
      )
    ),
    c(
      2008L,
      2015L,
      2023L
    )
  )
)

stopifnot(
  identical(
    sort(
      unique(
        test.data$TargetYear
      )
    ),
    c(
      2009L,
      2016L,
      2024L
    )
  )
)


# Spatial prediction from M2

test.data$pred_spatial_M2_raw <- predict(
  rf.M2,
  newdata = test.data[
    ,
    M2.variables,
    drop = FALSE
  ]
)


# Spatial persistence benchmark

test.data$pred_spatial_persistence_raw <-
  test.data$spatial_t


# The spatial component is defined as a deviation from the
# AOI-wide mean. Predictions are therefore centred within
# each target year.

test.data <- test.data |>
  group_by(
    TargetYear
  ) |>
  mutate(
    
    pred_spatial_M2 =
      pred_spatial_M2_raw -
      mean(
        pred_spatial_M2_raw,
        na.rm = TRUE
      ),
    
    pred_spatial_persistence =
      pred_spatial_persistence_raw -
      mean(
        pred_spatial_persistence_raw,
        na.rm = TRUE
      )
    
  ) |>
  ungroup()


# ------------------------------------------------------------
# 12. RECONSTRUCT TOTAL MALARIA INCIDENCE
# ------------------------------------------------------------

# Temporal prediction by persistence:
# predicted AOI-wide mean at t+1 = observed AOI-wide mean at t.

test.data$pred_mu <-
  test.data$mu_t


# Observed target log-incidence

test.data$target_log_incidence <-
  test.data$target_mu +
  test.data$target_spatial


# M2 reconstruction

test.data$pred_log_incidence_M2 <-
  test.data$pred_mu +
  test.data$pred_spatial_M2


# Temporal + spatial persistence benchmark

test.data$pred_log_incidence_persistence <-
  test.data$pred_mu +
  test.data$pred_spatial_persistence


# Temporal persistence only

test.data$pred_log_incidence_AOI_mean <-
  test.data$pred_mu


# Back-transform from log1p scale

test.data$target_incidence <-
  expm1(
    test.data$target_log_incidence
  )

test.data$pred_incidence_M2 <-
  pmax(
    0,
    expm1(
      test.data$pred_log_incidence_M2
    )
  )

test.data$pred_incidence_persistence <-
  pmax(
    0,
    expm1(
      test.data$pred_log_incidence_persistence
    )
  )

test.data$pred_incidence_AOI_mean <-
  pmax(
    0,
    expm1(
      test.data$pred_log_incidence_AOI_mean
    )
  )


# ------------------------------------------------------------
# 13. EXTERNAL TEST - OVERALL PERFORMANCE
# ------------------------------------------------------------

External.performance <- bind_rows(
  
  # Spatial component
  
  cbind(
    Scale = "Spatial component",
    ModelID = "M2",
    Model = "M2: EO_t + EO_t-1",
    
    calculate_metrics(
      test.data$target_spatial,
      test.data$pred_spatial_M2
    )
  ),
  
  cbind(
    Scale = "Spatial component",
    ModelID = "Spatial_persistence",
    Model = "Spatial persistence",
    
    calculate_metrics(
      test.data$target_spatial,
      test.data$pred_spatial_persistence
    )
  ),
  
  # Log-incidence
  
  cbind(
    Scale = "Log incidence",
    ModelID = "M2",
    Model = "Temporal persistence + M2",
    
    calculate_metrics(
      test.data$target_log_incidence,
      test.data$pred_log_incidence_M2
    )
  ),
  
  cbind(
    Scale = "Log incidence",
    ModelID = "Spatial_persistence",
    Model = "Temporal + spatial persistence",
    
    calculate_metrics(
      test.data$target_log_incidence,
      test.data$pred_log_incidence_persistence
    )
  ),
  
  cbind(
    Scale = "Log incidence",
    ModelID = "AOI_mean_only",
    Model = "Temporal persistence only",
    
    calculate_metrics(
      test.data$target_log_incidence,
      test.data$pred_log_incidence_AOI_mean
    )
  ),
  
  # Original incidence scale
  
  cbind(
    Scale = "Incidence",
    ModelID = "M2",
    Model = "Temporal persistence + M2",
    
    calculate_metrics(
      test.data$target_incidence,
      test.data$pred_incidence_M2
    )
  ),
  
  cbind(
    Scale = "Incidence",
    ModelID = "Spatial_persistence",
    Model = "Temporal + spatial persistence",
    
    calculate_metrics(
      test.data$target_incidence,
      test.data$pred_incidence_persistence
    )
  ),
  
  cbind(
    Scale = "Incidence",
    ModelID = "AOI_mean_only",
    Model = "Temporal persistence only",
    
    calculate_metrics(
      test.data$target_incidence,
      test.data$pred_incidence_AOI_mean
    )
  )
)


print(
  External.performance
)


write.csv(
  External.performance,
  file.path(
    results.dir,
    "m2_external_test.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. EXTERNAL TEST PERFORMANCE BY TARGET YEAR
# ------------------------------------------------------------

External.performance.by.year <-
  test.data |>
  group_by(
    Year,
    TargetYear
  ) |>
  group_modify(
    function(.x, .y) {
      
      calculate_metrics(
        observed =
          .x$target_log_incidence,
        
        predicted =
          .x$pred_log_incidence_M2
      )
    }
  ) |>
  ungroup()


print(
  External.performance.by.year
)


write.csv(
  External.performance.by.year,
  file.path(
    results.dir,
    "m2_external_test_by_year.csv"
  ),
  row.names = FALSE
)


write.csv(
  rf.M2$bestTune,
  file.path(
    results.dir,
    "m2_best_tuning.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. FINAL SUMMARY
# ------------------------------------------------------------

cat(
  "\n========================================\n",
  "M2 malaria forecasting analysis complete\n",
  "========================================\n\n",
  "Best tuning:\n",
  sep = ""
)

print(
  rf.M2$bestTune
)

cat(
  "\nResults written to:\n",
  results.dir,
  "\n"
)


