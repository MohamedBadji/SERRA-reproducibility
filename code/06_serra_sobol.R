# ============================================================
# SERRA global sensitivity analysis
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#   03_serra_run.R
#
# This file performs:
# - construction of Sobol parameter domains;
# - Monte Carlo evaluation of RMSE_mix;
# - Jansen first-order and total-order sensitivity analysis;
# - reuse of cached Sobol results when available;
# - assembly of the final Sobol results table.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "dplyr",
  "sensitivity",
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
# Cache directory
# ============================================================

CACHE_DIR <- "cache"

if (!dir.exists(CACHE_DIR)) {
  dir.create(
    CACHE_DIR,
    recursive = TRUE
  )
}


# ============================================================
# Sobol settings
# ============================================================

SOBOL_N <- 1000
SOBOL_MARGIN <- 0.20
SOBOL_MC_SIZE <- 150
SOBOL_NBOOT <- 100
SOBOL_SEED <- 123

FORCE_SOBOL_RECALCULATION <- FALSE


# ============================================================
# Parameter names
# ============================================================

SOBOL_PARAMETERS <- c(
  "P1",
  "P2",
  "P3",
  "P5",
  "TT0_LAI",
  "P6",
  "P7",
  "P9"
)


# ============================================================
# Parameter-domain construction
# ============================================================

build_sobol_bounds <- function(
    calibration_row,
    margin = SOBOL_MARGIN
) {
  
  parameter_center <- calibration_row |>
    dplyr::select(
      dplyr::all_of(
        SOBOL_PARAMETERS
      )
    ) |>
    unlist(
      use.names = FALSE
    ) |>
    as.numeric()
  
  
  bounds <- tibble::tibble(
    param = SOBOL_PARAMETERS,
    center = parameter_center
  ) |>
    dplyr::mutate(
      
      lower = dplyr::case_when(
        
        param == "TT0_LAI" ~ 0,
        
        TRUE ~ pmax(
          center * (1 - margin),
          1e-8
        )
      ),
      
      upper = dplyr::case_when(
        
        param == "TT0_LAI" ~ 10,
        
        TRUE ~
          center * (1 + margin)
      )
    )
  
  
  if (
    any(!is.finite(bounds$lower)) ||
    any(!is.finite(bounds$upper))
  ) {
    stop(
      "Non-finite Sobol parameter bounds detected.",
      call. = FALSE
    )
  }
  
  
  if (
    any(
      bounds$upper <=
      bounds$lower
    )
  ) {
    stop(
      "Invalid Sobol parameter bounds detected.",
      call. = FALSE
    )
  }
  
  
  bounds
}


# ============================================================
# Uniform parameter sampling
# ============================================================

sample_parameter_matrix <- function(
    bounds,
    n
) {
  
  sampled_values <- matrix(
    NA_real_,
    nrow = n,
    ncol = nrow(bounds)
  )
  
  colnames(sampled_values) <-
    bounds$param
  
  
  for (
    j in seq_len(
      nrow(bounds)
    )
  ) {
    
    sampled_values[, j] <-
      stats::runif(
        n = n,
        min = bounds$lower[j],
        max = bounds$upper[j]
      )
  }
  
  
  as.data.frame(
    sampled_values
  )
}


# ============================================================
# Noise specification
# ============================================================

get_sobol_noise <- function(
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
# Sobol model evaluation
# ============================================================

evaluate_sobol_rmse <- function(
    parameter_matrix,
    calibration_row,
    m_ens = SOBOL_MC_SIZE
) {
  
  obs_variety <- observations |>
    dplyr::filter(
      VAR ==
        calibration_row$VAR[[1]]
    ) |>
    dplyr::arrange(
      DAS
    )
  
  
  observation_index <- match(
    obs_variety$DAS,
    days
  )
  
  
  if (
    anyNA(
      observation_index
    )
  ) {
    stop(
      paste0(
        "At least one observed DAS value is absent from the ",
        "simulation time axis for variety ",
        calibration_row$VAR[[1]],
        "."
      ),
      call. = FALSE
    )
  }
  
  
  noise <- get_sobol_noise(
    calibration_row$Noise[[1]]
  )
  
  
  output <- numeric(
    nrow(
      parameter_matrix
    )
  )
  
  
  for (
    i in seq_len(
      nrow(parameter_matrix)
    )
  ) {
    
    params <- base_params
    
    params$P1 <-
      parameter_matrix$P1[i]
    
    params$P2 <-
      parameter_matrix$P2[i]
    
    params$P3 <-
      parameter_matrix$P3[i]
    
    params$P5 <-
      parameter_matrix$P5[i]
    
    params$TT0_LAI <-
      parameter_matrix$TT0_LAI[i]
    
    params$P6 <-
      parameter_matrix$P6[i]
    
    params$P7 <-
      parameter_matrix$P7[i]
    
    params$P9 <-
      parameter_matrix$P9[i]
    
    
    ensemble <-
      simulate_ensemble_v3_senB(
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        params = params,
        noise = noise,
        m = m_ens,
        seed = 100 + i,
        k_step = 0.005,
        k_senB = 0.01
      )
    
    
    predicted_biomass <-
      ensemble$B$mean[
        observation_index
      ]
    
    predicted_lai <-
      ensemble$LAI$mean[
        observation_index
      ]
    
    
    rmse_biomass <- sqrt(
      mean(
        (
          obs_variety$BiomassDW_g -
            predicted_biomass
        )^2,
        na.rm = TRUE
      )
    )
    
    
    rmse_lai <- sqrt(
      mean(
        (
          obs_variety$LAI -
            predicted_lai
        )^2,
        na.rm = TRUE
      )
    )
    
    
    output[i] <-
      0.7 * rmse_biomass +
      0.3 * rmse_lai
  }
  
  
  output
}


# ============================================================
# Single Sobol analysis
# ============================================================

run_sobol_analysis <- function(
    variety,
    noise_kind,
    n = SOBOL_N,
    margin = SOBOL_MARGIN,
    m_ens = SOBOL_MC_SIZE,
    nboot = SOBOL_NBOOT,
    seed = SOBOL_SEED
) {
  
  calibration_row <- calibration_results |>
    dplyr::filter(
      VAR == variety,
      Noise == noise_kind
    ) |>
    dplyr::slice(1)
  
  
  if (
    nrow(
      calibration_row
    ) == 0
  ) {
    stop(
      paste0(
        "No calibration result found for ",
        variety,
        " / ",
        noise_kind,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  bounds <- build_sobol_bounds(
    calibration_row = calibration_row,
    margin = margin
  )
  
  
  set.seed(
    seed
  )
  
  
  X1 <- sample_parameter_matrix(
    bounds = bounds,
    n = n
  )
  
  
  X2 <- sample_parameter_matrix(
    bounds = bounds,
    n = n
  )
  
  
  model_function <- function(X) {
    
    result <- evaluate_sobol_rmse(
      parameter_matrix = X,
      calibration_row = calibration_row,
      m_ens = m_ens
    )
    
    stopifnot(
      is.numeric(result)
    )
    
    result
  }
  
  
  sobol_object <-
    sensitivity::soboljansen(
      model = model_function,
      X1 = X1,
      X2 = X2,
      nboot = nboot
    )
  
  
  first_order <-
    sobol_object$S
  
  total_order <-
    sobol_object$T
  
  
  sobol_table <- tibble::tibble(
    
    VAR =
      rep(
        variety,
        nrow(first_order)
      ),
    
    Noise =
      rep(
        noise_kind,
        nrow(first_order)
      ),
    
    param =
      rownames(
        first_order
      ),
    
    Si =
      first_order[
        ,
        "original"
      ],
    
    Si_low =
      first_order[
        ,
        "min. c.i."
      ],
    
    Si_high =
      first_order[
        ,
        "max. c.i."
      ],
    
    STi =
      total_order[
        ,
        "original"
      ],
    
    STi_low =
      total_order[
        ,
        "min. c.i."
      ],
    
    STi_high =
      total_order[
        ,
        "max. c.i."
      ]
  )
  
  
  list(
    bounds = bounds,
    sobol = sobol_table,
    object = sobol_object,
    X1 = X1,
    X2 = X2
  )
}


# ============================================================
# Cache helper
# ============================================================

load_or_run_sobol <- function(
    cache_file,
    variety,
    noise_kind,
    force_recalculation =
      FORCE_SOBOL_RECALCULATION
) {
  
  if (
    file.exists(cache_file) &&
    !force_recalculation
  ) {
    
    message(
      "Loading cached Sobol analysis: ",
      cache_file
    )
    
    return(
      readRDS(
        cache_file
      )
    )
  }
  
  
  message(
    "Running Sobol analysis for ",
    variety,
    " / ",
    noise_kind,
    "."
  )
  
  
  result <- run_sobol_analysis(
    variety = variety,
    noise_kind = noise_kind,
    n = SOBOL_N,
    margin = SOBOL_MARGIN,
    m_ens = SOBOL_MC_SIZE,
    nboot = SOBOL_NBOOT,
    seed = SOBOL_SEED
  )
  
  
  saveRDS(
    result,
    cache_file
  )
  
  
  result
}


# ============================================================
# Cache files
# ============================================================

SOBOL_V1_G_FILE <- file.path(
  CACHE_DIR,
  "sob_V1_G_1000.rds"
)

SOBOL_V1_L_FILE <- file.path(
  CACHE_DIR,
  "sob_V1_L_1000.rds"
)

SOBOL_V2_G_FILE <- file.path(
  CACHE_DIR,
  "sob_V2_G_1000.rds"
)

SOBOL_V2_L_FILE <- file.path(
  CACHE_DIR,
  "sob_V2_L_1000.rds"
)


# ============================================================
# Run or load the four final analyses
# ============================================================

sobol_V1_gaussian <- load_or_run_sobol(
  cache_file =
    SOBOL_V1_G_FILE,
  variety =
    "V1",
  noise_kind =
    "gaussian"
)


sobol_V1_levy <- load_or_run_sobol(
  cache_file =
    SOBOL_V1_L_FILE,
  variety =
    "V1",
  noise_kind =
    "levy"
)


sobol_V2_gaussian <- load_or_run_sobol(
  cache_file =
    SOBOL_V2_G_FILE,
  variety =
    "V2",
  noise_kind =
    "gaussian"
)


sobol_V2_levy <- load_or_run_sobol(
  cache_file =
    SOBOL_V2_L_FILE,
  variety =
    "V2",
  noise_kind =
    "levy"
)


# ============================================================
# Final Sobol tables
# ============================================================

sobol_table_V1_gaussian <-
  sobol_V1_gaussian$sobol

sobol_table_V1_levy <-
  sobol_V1_levy$sobol

sobol_table_V2_gaussian <-
  sobol_V2_gaussian$sobol

sobol_table_V2_levy <-
  sobol_V2_levy$sobol


sobol_results <- dplyr::bind_rows(
  sobol_table_V1_gaussian,
  sobol_table_V1_levy,
  sobol_table_V2_gaussian,
  sobol_table_V2_levy
) |>
  dplyr::arrange(
    VAR,
    Noise,
    param
  )


# ============================================================
# Save combined Sobol results
# ============================================================

saveRDS(
  sobol_results,
  file.path(
    CACHE_DIR,
    "sobol_results.rds"
  )
)


# ============================================================
# Save parameter domains
# ============================================================

sobol_bounds <- dplyr::bind_rows(
  
  sobol_V1_gaussian$bounds |>
    dplyr::mutate(
      VAR = "V1",
      Noise = "gaussian"
    ),
  
  sobol_V1_levy$bounds |>
    dplyr::mutate(
      VAR = "V1",
      Noise = "levy"
    ),
  
  sobol_V2_gaussian$bounds |>
    dplyr::mutate(
      VAR = "V2",
      Noise = "gaussian"
    ),
  
  sobol_V2_levy$bounds |>
    dplyr::mutate(
      VAR = "V2",
      Noise = "levy"
    )
) |>
  dplyr::select(
    VAR,
    Noise,
    param,
    center,
    lower,
    upper
  )


saveRDS(
  sobol_bounds,
  file.path(
    CACHE_DIR,
    "sobol_bounds.rds"
  )
)


# ============================================================
# Console output
# ============================================================

print(
  sobol_results
)

print(
  sobol_bounds
)