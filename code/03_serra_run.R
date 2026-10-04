# ============================================================
# SERRA analysis workflow
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#
# Main tasks:
# - load model inputs;
# - compute thermal time;
# - calibrate Gaussian and Levy-type formulations;
# - reuse cached calibration results when available;
# - generate final Monte Carlo predictions;
# - build observed-versus-predicted tables;
# - assemble the final calibration-results table.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "dplyr",
  "purrr",
  "readxl",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}


# ============================================================
# Source SERRA functions
# ============================================================
source(file.path("code", "01_serra_simulator.R"))
source(file.path("code", "02_serra_core.R"))

# ============================================================
# Project paths
# ============================================================

DATA_DIR <- "data"
CACHE_DIR <- "cache"

if (!dir.exists(CACHE_DIR)) {
  dir.create(
    CACHE_DIR,
    recursive = TRUE
  )
}


# ============================================================
# Input files
# ============================================================

OBSERVATION_FILE <- file.path(
  DATA_DIR,
  "by_variety_stage_means_2023.rds"
)

BASE_PARAMETER_FILE <- file.path(
  DATA_DIR,
  "base_params_used.rds"
)

CLIMATE_FILE <- file.path(
  DATA_DIR,
  "climat_fanaye1.xlsx"
)


required_input_files <- c(
  OBSERVATION_FILE,
  BASE_PARAMETER_FILE,
  CLIMATE_FILE
)

missing_input_files <- required_input_files[
  !file.exists(required_input_files)
]

if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "Missing input files:\n",
      paste(
        paste0(" - ", missing_input_files),
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# Load model inputs
# ============================================================

observations <- readRDS(
  OBSERVATION_FILE
)

base_params <- readRDS(
  BASE_PARAMETER_FILE
)

climate <- readxl::read_excel(
  CLIMATE_FILE
)


# ============================================================
# Validate climate inputs
# ============================================================

required_climate_columns <- c(
  "DAS",
  "TMIN",
  "TMAX",
  "PAR"
)

missing_climate_columns <- setdiff(
  required_climate_columns,
  names(climate)
)

if (length(missing_climate_columns) > 0) {
  stop(
    paste0(
      "Missing climate columns: ",
      paste(
        missing_climate_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


climate <- climate |>
  dplyr::mutate(
    DAS = as.integer(DAS)
  ) |>
  dplyr::arrange(DAS)


# ============================================================
# Climate forcing
# ============================================================

TMIN <- climate$TMIN
TMAX <- climate$TMAX
PAR <- climate$PAR
days <- climate$DAS


# ============================================================
# Thermal time
# ============================================================

thermal_time_table <- build_thermal_time_table(
  days = days,
  TMIN = TMIN,
  TMAX = TMAX,
  Tbase = base_params$P4
)


varieties <- sort(
  unique(
    observations$VAR
  )
)


# ============================================================
# Calibration helper
# ============================================================

run_calibration <- function(
    noise_kind,
    noise_args,
    m_ens_lai = 50,
    m_ens_biomass = 50,
    buffer_fraction = 0.10,
    k_step = 0.005,
    k_senB = 0.01,
    TT0_start_grid = c(
      200,
      400,
      600
    )
) {
  
  results <- purrr::map_dfr(
    varieties,
    function(variety) {
      
      fit_serra_by_variety(
        variety = variety,
        base_params = base_params,
        thermal_time_table = thermal_time_table,
        observations = observations,
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        days = days,
        noise_kind = noise_kind,
        noise_args = noise_args,
        m_ens_lai = m_ens_lai,
        m_ens_biomass = m_ens_biomass,
        buffer_fraction = buffer_fraction,
        k_step = k_step,
        k_senB = k_senB,
        TT0_start_grid = TT0_start_grid
      )
    }
  )
  
  results |>
    dplyr::mutate(
      RMSE_mix =
        0.7 * RMSE_B +
        0.3 * RMSE_L
    ) |>
    dplyr::select(
      VAR,
      Noise,
      P1,
      P2,
      P3,
      P5,
      TT0_LAI,
      P6,
      P7,
      P9,
      RMSE_B,
      MAE_B,
      R2_B,
      RMSE_L,
      MAE_L,
      R2_L,
      RMSE_mix,
      Obj_L,
      Obj_B,
      Conv_L,
      Conv_B
    )
}


# ============================================================
# Cached calibration helper
# ============================================================

load_or_run_calibration <- function(
    cache_file,
    noise_kind,
    noise_args,
    force_recalibration = FALSE
) {
  
  if (
    file.exists(cache_file) &&
    !force_recalibration
  ) {
    
    message(
      "Loading cached calibration: ",
      cache_file
    )
    
    return(
      readRDS(
        cache_file
      )
    )
  }
  
  
  message(
    "Running calibration for ",
    noise_kind,
    " formulation."
  )
  
  
  results <- run_calibration(
    noise_kind = noise_kind,
    noise_args = noise_args,
    m_ens_lai = 50,
    m_ens_biomass = 50
  )
  
  
  saveRDS(
    results,
    cache_file
  )
  
  
  results
}


# ============================================================
# Calibration cache files
# ============================================================

GAUSSIAN_CALIBRATION_FILE <- file.path(
  CACHE_DIR,
  "calibration_gaussian.rds"
)

LEVY_CALIBRATION_FILE <- file.path(
  CACHE_DIR,
  "calibration_levy.rds"
)


# ============================================================
# Recalibration control
# ============================================================

FORCE_RECALIBRATION <- FALSE


# ============================================================
# Gaussian calibration
# ============================================================

calibration_gaussian <- load_or_run_calibration(
  cache_file =
    GAUSSIAN_CALIBRATION_FILE,
  
  noise_kind =
    "gaussian",
  
  noise_args = list(
    sigma1 = 0.10,
    sigma2 = 0.05
  ),
  
  force_recalibration =
    FORCE_RECALIBRATION
)


# ============================================================
# Levy-type calibration
# ============================================================

calibration_levy <- load_or_run_calibration(
  cache_file =
    LEVY_CALIBRATION_FILE,
  
  noise_kind =
    "levy",
  
  noise_args = list(
    alpha = 1.5,
    beta = 0,
    scale = 0.20
  ),
  
  force_recalibration =
    FORCE_RECALIBRATION
)


# ============================================================
# Combined calibration results
# ============================================================

calibration_results <- dplyr::bind_rows(
  calibration_gaussian,
  calibration_levy
) |>
  dplyr::arrange(
    VAR,
    Noise
  )


# ============================================================
# Final Gaussian predictions
# ============================================================

predictions_gaussian <- purrr::map_dfr(
  varieties,
  function(variety) {
    
    predict_serra_variety(
      variety = variety,
      calibration_results =
        calibration_gaussian,
      base_params = base_params,
      TMIN = TMIN,
      TMAX = TMAX,
      PAR = PAR,
      days = days,
      noise_kind = "gaussian",
      noise_args = list(
        sigma1 = 0.10,
        sigma2 = 0.05
      ),
      m = 500,
      seed = 999,
      k_step = 0.005,
      k_senB = 0.01
    )
  }
)


# ============================================================
# Final Levy-type predictions
# ============================================================

predictions_levy <- purrr::map_dfr(
  varieties,
  function(variety) {
    
    predict_serra_variety(
      variety = variety,
      calibration_results =
        calibration_levy,
      base_params = base_params,
      TMIN = TMIN,
      TMAX = TMAX,
      PAR = PAR,
      days = days,
      noise_kind = "levy",
      noise_args = list(
        alpha = 1.5,
        beta = 0,
        scale = 0.20
      ),
      m = 500,
      seed = 999,
      k_step = 0.005,
      k_senB = 0.01
    )
  }
)


# ============================================================
# Combined prediction table
# ============================================================

prediction_results <- dplyr::bind_rows(
  predictions_gaussian,
  predictions_levy
)


# ============================================================
# Observed-versus-predicted tables
# ============================================================

comparison_gaussian <- purrr::map_dfr(
  varieties,
  function(variety) {
    
    prediction_variety <-
      predictions_gaussian |>
      dplyr::filter(
        VAR == variety
      )
    
    build_comparison_table(
      variety = variety,
      prediction_table =
        prediction_variety,
      observations =
        observations
    )
  }
)


comparison_levy <- purrr::map_dfr(
  varieties,
  function(variety) {
    
    prediction_variety <-
      predictions_levy |>
      dplyr::filter(
        VAR == variety
      )
    
    build_comparison_table(
      variety = variety,
      prediction_table =
        prediction_variety,
      observations =
        observations
    )
  }
)


comparison_results <- dplyr::bind_rows(
  comparison_gaussian,
  comparison_levy
)


# ============================================================
# Save reusable outputs
# ============================================================

saveRDS(
  thermal_time_table,
  file.path(
    CACHE_DIR,
    "thermal_time_table.rds"
  )
)

saveRDS(
  calibration_results,
  file.path(
    CACHE_DIR,
    "calibration_results.rds"
  )
)

saveRDS(
  prediction_results,
  file.path(
    CACHE_DIR,
    "prediction_results.rds"
  )
)

saveRDS(
  comparison_results,
  file.path(
    CACHE_DIR,
    "comparison_results.rds"
  )
)


# ============================================================
# Console summary
# ============================================================

print(
  calibration_results |>
    dplyr::select(
      VAR,
      Noise,
      P1,
      P2,
      P3,
      P5,
      TT0_LAI,
      P6,
      P7,
      P9,
      RMSE_B,
      MAE_B,
      R2_B,
      RMSE_L,
      MAE_L,
      R2_L,
      RMSE_mix
    )
)