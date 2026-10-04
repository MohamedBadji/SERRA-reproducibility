# ============================================================
# SERRA diagnostic analyses
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#   03_serra_run.R
#
# This file performs:
# - AIC/BIC comparison based on standardized residuals;
# - numerical terminal-stability diagnostic;
# - reuse of existing diagnostic results when available;
# - assembly of final diagnostic tables.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "dplyr",
  "purrr",
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
# Source analysis workflow
# ============================================================

source(file.path("code", "03_serra_run.R"))

# ============================================================
# Output directory
# ============================================================

RESULTS_DIR <- "results"

if (!dir.exists(RESULTS_DIR)) {
  dir.create(
    RESULTS_DIR,
    recursive = TRUE
  )
}


# ============================================================
# Diagnostic settings
# ============================================================

DIAGNOSTIC_MC_SIZE <- 500
DIAGNOSTIC_SEED <- 123
DIAGNOSTIC_EPS <- 1e-6

FORCE_DIAGNOSTIC_RECALCULATION <- FALSE


# ============================================================
# Final calibration table
# ============================================================

diagnostic_calibration_results <- calibration_results |>
  dplyr::filter(
    Noise %in% c(
      "gaussian",
      "levy"
    )
  )


# ============================================================
# Parameter reconstruction
# ============================================================

build_diagnostic_params <- function(
    calibration_row
) {
  
  params <- base_params
  
  params$P1 <-
    calibration_row$P1[[1]]
  
  params$P2 <-
    calibration_row$P2[[1]]
  
  params$P3 <-
    calibration_row$P3[[1]]
  
  params$P5 <-
    calibration_row$P5[[1]]
  
  params$TT0_LAI <-
    calibration_row$TT0_LAI[[1]]
  
  params$P6 <-
    calibration_row$P6[[1]]
  
  params$P7 <-
    calibration_row$P7[[1]]
  
  params$P9 <-
    calibration_row$P9[[1]]
  
  
  params
}


# ============================================================
# Diagnostic noise specification
# ============================================================

build_diagnostic_noise <- function(
    noise_kind
) {
  
  if (
    noise_kind == "gaussian"
  ) {
    
    return(
      build_noise_spec(
        noise_kind = "gaussian",
        noise_args = list(
          sigma1 = 0.10,
          sigma2 = 0.05
        )
      )
    )
  }
  
  
  if (
    noise_kind == "levy"
  ) {
    
    return(
      build_noise_spec(
        noise_kind = "levy",
        noise_args = list(
          alpha = 1.5,
          beta = 0,
          scale = 0.20
        )
      )
    )
  }
  
  
  stop(
    paste0(
      "Unsupported stochastic formulation: ",
      noise_kind
    ),
    call. = FALSE
  )
}


# ============================================================
# AIC / BIC diagnostic
# ============================================================

compute_aic_bic <- function(
    calibration_row
) {
  
  variety <-
    calibration_row$VAR[[1]]
  
  noise_kind <-
    calibration_row$Noise[[1]]
  
  
  obs_variety <- observations |>
    dplyr::filter(
      VAR == variety
    ) |>
    dplyr::arrange(
      DAS
    )
  
  
  prediction <- predict_serra_variety(
    variety = variety,
    calibration_results = calibration_row,
    base_params = base_params,
    TMIN = TMIN,
    TMAX = TMAX,
    PAR = PAR,
    days = days,
    noise_kind = noise_kind,
    noise_args =
      if (
        noise_kind == "gaussian"
      ) {
        list(
          sigma1 = 0.10,
          sigma2 = 0.05
        )
      } else {
        list(
          alpha = 1.5,
          beta = 0,
          scale = 0.20
        )
      },
    m = 500,
    seed = 999,
    k_step = 0.005,
    k_senB = 0.01
  )
  
  
  comparison <- obs_variety |>
    dplyr::left_join(
      prediction |>
        dplyr::select(
          DAS,
          B_mean,
          LAI_mean
        ),
      by = "DAS"
    ) |>
    dplyr::mutate(
      residual_biomass =
        BiomassDW_g - B_mean,
      residual_lai =
        LAI - LAI_mean
    )
  
  
  sd_biomass <- stats::sd(
    obs_variety$BiomassDW_g,
    na.rm = TRUE
  )
  
  sd_lai <- stats::sd(
    obs_variety$LAI,
    na.rm = TRUE
  )
  
  
  standardized_residuals <- c(
    comparison$residual_biomass /
      sd_biomass,
    comparison$residual_lai /
      sd_lai
  )
  
  
  standardized_residuals <-
    standardized_residuals[
      is.finite(
        standardized_residuals
      )
    ]
  
  
  n <- length(
    standardized_residuals
  )
  
  
  rss <- sum(
    standardized_residuals^2
  )
  
  
  sigma2_hat <- rss / n
  
  
  log_likelihood <-
    -0.5 *
    n *
    (
      log(
        2 * pi
      ) +
        log(
          sigma2_hat
        ) +
        1
    )
  
  
  k <- 8
  
  
  aic <-
    2 * k -
    2 * log_likelihood
  
  
  bic <-
    log(n) * k -
    2 * log_likelihood
  
  
  tibble::tibble(
    VAR = variety,
    Noise = noise_kind,
    n = n,
    k = k,
    logLik = log_likelihood,
    AIC = aic,
    BIC = bic
  )
}


# ============================================================
# Terminal-stability diagnostic
# ============================================================

compute_terminal_stability <- function(
    calibration_row,
    m = DIAGNOSTIC_MC_SIZE,
    seed = DIAGNOSTIC_SEED,
    eps = DIAGNOSTIC_EPS
) {
  
  variety <-
    calibration_row$VAR[[1]]
  
  noise_kind <-
    calibration_row$Noise[[1]]
  
  
  params <- build_diagnostic_params(
    calibration_row
  )
  
  
  noise <- build_diagnostic_noise(
    noise_kind
  )
  
  
  set.seed(
    seed
  )
  
  
  lambda_values <- numeric(
    m
  )
  
  
  for (
    i in seq_len(
      m
    )
  ) {
    
    simulation <- simulate_once_v3_senB(
      TMIN = TMIN,
      TMAX = TMAX,
      PAR = PAR,
      params = params,
      noise = noise,
      project = TRUE,
      seed = NULL,
      k_step = 0.005,
      k_senB = 0.01
    )
    
    
    biomass <- simulation$B
    
    
    terminal_index <- seq(
      floor(
        0.70 *
          length(
            biomass
          )
      ),
      length(
        biomass
      ) - 1
    )
    
    
    log_increment <- diff(
      log(
        biomass +
          eps
      )
    )
    
    
    lambda_values[i] <-
      mean(
        log_increment[
          terminal_index
        ],
        na.rm = TRUE
      )
  }
  
  
  tibble::tibble(
    VAR = variety,
    Noise = noise_kind,
    Lambda_mean =
      mean(
        lambda_values,
        na.rm = TRUE
      ),
    Lambda_sd =
      stats::sd(
        lambda_values,
        na.rm = TRUE
      ),
    S2_satisfied =
      mean(
        lambda_values,
        na.rm = TRUE
      ) < 0
  )
}


# ============================================================
# Result files
# ============================================================

AIC_BIC_FILE <- file.path(
  "results",
  "table6_aic_bic.rds"
)

STABILITY_FILE <- file.path(
  "results",
  "table9_s2.rds"
)

# ============================================================
# AIC / BIC: load or compute
# ============================================================

if (
  file.exists(
    AIC_BIC_FILE
  ) &&
  !FORCE_DIAGNOSTIC_RECALCULATION
) {
  
  message(
    "Loading existing AIC/BIC results: ",
    AIC_BIC_FILE
  )
  
  table6_aic_bic <- readRDS(
    AIC_BIC_FILE
  )
  
} else {
  
  message(
    "Computing AIC/BIC diagnostic."
  )
  
  
  table6_aic_bic <-
    diagnostic_calibration_results |>
    split(
      seq_len(
        nrow(
          diagnostic_calibration_results
        )
      )
    ) |>
    purrr::map_dfr(
      compute_aic_bic
    )
  
  
  saveRDS(
    table6_aic_bic,
    AIC_BIC_FILE
  )
}


# ============================================================
# Stability diagnostic: load or compute
# ============================================================

if (
  file.exists(
    STABILITY_FILE
  ) &&
  !FORCE_DIAGNOSTIC_RECALCULATION
) {
  
  message(
    "Loading existing stability diagnostic: ",
    STABILITY_FILE
  )
  
  table9_s2 <- readRDS(
    STABILITY_FILE
  )
  
} else {
  
  message(
    "Computing terminal-stability diagnostic."
  )
  
  
  table9_s2 <-
    diagnostic_calibration_results |>
    split(
      seq_len(
        nrow(
          diagnostic_calibration_results
        )
      )
    ) |>
    purrr::map_dfr(
      ~ compute_terminal_stability(
        calibration_row = .x,
        m = DIAGNOSTIC_MC_SIZE,
        seed = DIAGNOSTIC_SEED,
        eps = DIAGNOSTIC_EPS
      )
    )
  
  
  saveRDS(
    table9_s2,
    STABILITY_FILE
  )
}


# ============================================================
# Console output
# ============================================================

print(
  table6_aic_bic
)

print(
  table9_s2
)