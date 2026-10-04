# ============================================================
# SERRA core functions
#
# Requires:
#   01_serra_simulator.R
#
# This file defines:
# - thermal-time construction;
# - stochastic-model specification;
# - calibration objective functions;
# - goodness-of-fit metrics;
# - observation-derived parameters;
# - sequential calibration;
# - final predictive simulation;
# - observed-versus-predicted tables.
# ============================================================


# ============================================================
# Dependency checks
# ============================================================

check_serra_dependencies <- function() {
  
  required_functions <- c(
    "make_noise",
    "simulate_ensemble_v3_senB"
  )
  
  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      mode = "function"
    )
  ]
  
  if (length(missing_functions) > 0) {
    stop(
      paste0(
        "Missing required functions: ",
        paste(missing_functions, collapse = ", "),
        ". Run 01_serra_simulator.R first."
      ),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ============================================================
# Thermal-time calculation
# ============================================================

build_thermal_time_table <- function(
    days,
    TMIN,
    TMAX,
    Tbase
) {
  
  stopifnot(
    length(days) == length(TMIN),
    length(days) == length(TMAX)
  )
  
  dTT <- pmax(
    (TMIN + TMAX) / 2 - Tbase,
    0
  )
  
  TT <- cumsum(
    c(
      0,
      dTT[-length(dTT)]
    )
  )
  
  tibble::tibble(
    DAS = days,
    TT = TT
  )
}


# ============================================================
# Stochastic-model specification
# ============================================================

build_noise_spec <- function(
    noise_kind = c("gaussian", "levy"),
    noise_args = list()
) {
  
  noise_kind <- match.arg(noise_kind)
  
  default_args <- list(
    kind = noise_kind,
    sigma1 = 0.10,
    sigma2 = 0.05,
    alpha = 1.5,
    beta = 0,
    scale = 0.20
  )
  
  noise_parameters <- utils::modifyList(
    default_args,
    noise_args
  )
  
  do.call(
    make_noise,
    noise_parameters
  )
}


# ============================================================
# LAI calibration objective
# ============================================================

objective_lai <- function(
    par,
    obs_tbl,
    base_params,
    P5_fixed,
    TMIN,
    TMAX,
    PAR,
    days,
    noise_kind = "gaussian",
    noise_args = list(),
    m_ens = 50,
    seed = 123,
    k_step = 0.005,
    k_senB = 0.01,
    var_min = 0.02,
    penalty_variance = 50,
    penalty_saturation = 20,
    saturation_fraction = 0.98,
    saturation_allowed = 0.25
) {
  
  params <- base_params
  
  params$P3 <- par[1]
  params$P5 <- P5_fixed
  params$TT0_LAI <- par[2]
  
  noise <- build_noise_spec(
    noise_kind = noise_kind,
    noise_args = noise_args
  )
  
  ensemble <- simulate_ensemble_v3_senB(
    TMIN = TMIN,
    TMAX = TMAX,
    PAR = PAR,
    params = params,
    noise = noise,
    m = m_ens,
    seed = seed,
    k_step = k_step,
    k_senB = k_senB
  )
  
  obs_index <- match(
    obs_tbl$DAS,
    days
  )
  
  if (anyNA(obs_index)) {
    stop(
      "At least one observed DAS value is absent from the simulation time axis.",
      call. = FALSE
    )
  }
  
  predicted_lai <- ensemble$LAI$mean[
    obs_index
  ]
  
  if (
    anyNA(predicted_lai) ||
    any(!is.finite(predicted_lai))
  ) {
    return(1e12)
  }
  
  rmse_lai <- sqrt(
    mean(
      (obs_tbl$LAI - predicted_lai)^2,
      na.rm = TRUE
    )
  )
  
  
  # Low-variability penalty
  
  predicted_variance <- stats::var(
    predicted_lai,
    na.rm = TRUE
  )
  
  variance_penalty <- if (
    !is.finite(predicted_variance)
  ) {
    
    penalty_variance * var_min
    
  } else if (
    predicted_variance < var_min
  ) {
    
    penalty_variance *
      (var_min - predicted_variance)
    
  } else {
    
    0
  }
  
  
  # Premature-saturation penalty
  
  saturation_rate <- mean(
    predicted_lai >=
      saturation_fraction *
      P5_fixed,
    na.rm = TRUE
  )
  
  saturation_penalty <- if (
    is.finite(saturation_rate) &&
    saturation_rate > saturation_allowed
  ) {
    
    penalty_saturation *
      (
        saturation_rate -
          saturation_allowed
      )
    
  } else {
    
    0
  }
  
  
  rmse_lai +
    variance_penalty +
    saturation_penalty
}


# ============================================================
# Biomass calibration objective
# ============================================================

objective_biomass <- function(
    par,
    obs_tbl,
    base_params,
    P3_fixed,
    P5_fixed,
    TT0_fixed,
    TMIN,
    TMAX,
    PAR,
    days,
    noise_kind = "gaussian",
    noise_args = list(),
    m_ens = 50,
    seed = 123,
    k_step = 0.005,
    k_senB = 0.01
) {
  
  params <- base_params
  
  params$P1 <- par[1]
  params$P2 <- par[2]
  params$P9 <- par[3]
  
  params$P3 <- P3_fixed
  params$P5 <- P5_fixed
  params$TT0_LAI <- TT0_fixed
  
  noise <- build_noise_spec(
    noise_kind = noise_kind,
    noise_args = noise_args
  )
  
  ensemble <- simulate_ensemble_v3_senB(
    TMIN = TMIN,
    TMAX = TMAX,
    PAR = PAR,
    params = params,
    noise = noise,
    m = m_ens,
    seed = seed,
    k_step = k_step,
    k_senB = k_senB
  )
  
  obs_index <- match(
    obs_tbl$DAS,
    days
  )
  
  if (anyNA(obs_index)) {
    stop(
      "At least one observed DAS value is absent from the simulation time axis.",
      call. = FALSE
    )
  }
  
  predicted_biomass <- ensemble$B$mean[
    obs_index
  ]
  
  if (
    anyNA(predicted_biomass) ||
    any(!is.finite(predicted_biomass))
  ) {
    return(1e12)
  }
  
  sqrt(
    mean(
      (
        obs_tbl$BiomassDW_g -
          predicted_biomass
      )^2,
      na.rm = TRUE
    )
  )
}


# ============================================================
# Goodness-of-fit metrics
# ============================================================

calculate_metrics <- function(
    observed,
    predicted
) {
  
  valid <- is.finite(observed) &
    is.finite(predicted)
  
  observed <- observed[valid]
  predicted <- predicted[valid]
  
  if (length(observed) < 2) {
    
    return(
      c(
        RMSE = NA_real_,
        MAE = NA_real_,
        R2 = NA_real_
      )
    )
  }
  
  rmse <- sqrt(
    mean(
      (observed - predicted)^2
    )
  )
  
  mae <- mean(
    abs(
      observed - predicted
    )
  )
  
  ss_res <- sum(
    (observed - predicted)^2
  )
  
  ss_tot <- sum(
    (
      observed -
        mean(observed)
    )^2
  )
  
  r2 <- if (ss_tot <= 0) {
    
    NA_real_
    
  } else {
    
    1 -
      ss_res /
      ss_tot
  }
  
  c(
    RMSE = rmse,
    MAE = mae,
    R2 = r2
  )
}


# ============================================================
# Observation-derived parameters
# ============================================================

get_fixed_p5 <- function(
    observations,
    multiplier = 1.00,
    p5_min = 1.5,
    p5_max = 6.0
) {
  
  lai_max <- suppressWarnings(
    max(
      observations$LAI,
      na.rm = TRUE
    )
  )
  
  if (!is.finite(lai_max)) {
    lai_max <- 2.5
  }
  
  pmin(
    pmax(
      lai_max * multiplier,
      p5_min
    ),
    p5_max
  )
}


derive_p6_p7 <- function(
    observations,
    thermal_time_table,
    buffer_fraction = 0.10
) {
  
  observed_thermal_time <- thermal_time_table[
    thermal_time_table$DAS %in%
      observations$DAS,
    ,
    drop = FALSE
  ]
  
  observed_thermal_time <- observed_thermal_time[
    order(
      observed_thermal_time$DAS
    ),
    ,
    drop = FALSE
  ]
  
  if (
    nrow(
      observed_thermal_time
    ) == 0
  ) {
    stop(
      paste(
        "P6 and P7 cannot be derived because",
        "no observed DAS values match the thermal-time table."
      ),
      call. = FALSE
    )
  }
  
  max_observed_TT <- max(
    observed_thermal_time$TT,
    na.rm = TRUE
  )
  
  if (!is.finite(max_observed_TT)) {
    stop(
      "Maximum thermal time at observation dates is not finite.",
      call. = FALSE
    )
  }
  
  P6 <- max_observed_TT
  
  P7 <- max_observed_TT *
    (
      1 -
        buffer_fraction
    )
  
  list(
    P6 = P6,
    P7 = P7
  )
}


# ============================================================
# Sequential calibration by variety
# ============================================================

fit_serra_by_variety <- function(
    variety,
    base_params,
    thermal_time_table,
    observations,
    TMIN,
    TMAX,
    PAR,
    days,
    noise_kind = c("gaussian", "levy"),
    noise_args = list(),
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
  
  check_serra_dependencies()
  
  noise_kind <- match.arg(
    noise_kind
  )
  
  obs_variety <- observations[
    observations$VAR == variety,
    ,
    drop = FALSE
  ]
  
  obs_variety <- obs_variety[
    order(
      obs_variety$DAS
    ),
    ,
    drop = FALSE
  ]
  
  if (nrow(obs_variety) == 0) {
    stop(
      paste0(
        "No observations found for variety ",
        variety,
        "."
      ),
      call. = FALSE
    )
  }
  
  
  # Thermal thresholds
  
  thresholds <- derive_p6_p7(
    observations = obs_variety,
    thermal_time_table = thermal_time_table,
    buffer_fraction = buffer_fraction
  )
  
  TT_end <- max(
    thermal_time_table$TT,
    na.rm = TRUE
  )
  
  P6 <- min(
    max(
      50,
      thresholds$P6
    ),
    TT_end
  )
  
  P7 <- min(
    max(
      50,
      thresholds$P7
    ),
    P6
  )
  
  
  local_params <- base_params
  
  local_params$P6 <- P6
  local_params$P7 <- P7
  
  
  # Maximum observed LAI
  
  P5_fixed <- get_fixed_p5(
    obs_variety
  )
  
  
  # ----------------------------------------------------------
  # Stage 1: LAI calibration
  # ----------------------------------------------------------
  
  lower_lai <- c(
    P3 = 1e-6,
    TT0 = 0
  )
  
  upper_lai <- c(
    P3 = 2.00,
    TT0 = TT_end
  )
  
  best_lai_fit <- NULL
  
  
  for (
    TT0_start in TT0_start_grid
  ) {
    
    TT0_start <- max(
      0,
      min(
        TT0_start,
        TT_end
      )
    )
    
    start_lai <- c(
      P3 = 0.05,
      TT0 = TT0_start
    )
    
    
    fit_lai <- stats::nlminb(
      start = start_lai,
      
      objective = function(par) {
        
        objective_lai(
          par = par,
          obs_tbl = obs_variety,
          base_params = local_params,
          P5_fixed = P5_fixed,
          TMIN = TMIN,
          TMAX = TMAX,
          PAR = PAR,
          days = days,
          noise_kind = noise_kind,
          noise_args = noise_args,
          m_ens = m_ens_lai,
          seed = 123,
          k_step = k_step,
          k_senB = k_senB
        )
      },
      
      lower = lower_lai,
      upper = upper_lai,
      
      control = list(
        eval.max = 1200,
        iter.max = 1200
      )
    )
    
    
    if (
      is.null(best_lai_fit) ||
      fit_lai$objective <
      best_lai_fit$objective
    ) {
      
      best_lai_fit <- fit_lai
    }
  }
  
  
  P3_hat <- as.numeric(
    best_lai_fit$par["P3"]
  )
  
  TT0_hat <- as.numeric(
    best_lai_fit$par["TT0"]
  )
  
  
  # ----------------------------------------------------------
  # Stage 2: biomass calibration
  # ----------------------------------------------------------
  
  lower_biomass <- c(
    P1 = 0.0005,
    P2 = 0.1,
    P9 = 1e-6
  )
  
  upper_biomass <- c(
    P1 = 2.00,
    P2 = 200,
    P9 = 0.020
  )
  
  start_biomass <- c(
    P1 = 0.02,
    P2 = 10.0,
    P9 = 0.002
  )
  
  
  fit_biomass <- stats::nlminb(
    start = start_biomass,
    
    objective = function(par) {
      
      objective_biomass(
        par = par,
        obs_tbl = obs_variety,
        base_params = local_params,
        P3_fixed = P3_hat,
        P5_fixed = P5_fixed,
        TT0_fixed = TT0_hat,
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        days = days,
        noise_kind = noise_kind,
        noise_args = noise_args,
        m_ens = m_ens_biomass,
        seed = 123,
        k_step = k_step,
        k_senB = k_senB
      )
    },
    
    lower = lower_biomass,
    upper = upper_biomass,
    
    control = list(
      eval.max = 2000,
      iter.max = 2000
    )
  )
  
  
  # ----------------------------------------------------------
  # Final parameter set
  # ----------------------------------------------------------
  
  params <- local_params
  
  params$P1 <- as.numeric(
    fit_biomass$par["P1"]
  )
  
  params$P2 <- as.numeric(
    fit_biomass$par["P2"]
  )
  
  params$P9 <- as.numeric(
    fit_biomass$par["P9"]
  )
  
  params$P3 <- P3_hat
  params$P5 <- P5_fixed
  params$TT0_LAI <- TT0_hat
  
  
  # ----------------------------------------------------------
  # Final goodness-of-fit evaluation
  # ----------------------------------------------------------
  
  noise <- build_noise_spec(
    noise_kind = noise_kind,
    noise_args = noise_args
  )
  
  ensemble <- simulate_ensemble_v3_senB(
    TMIN = TMIN,
    TMAX = TMAX,
    PAR = PAR,
    params = params,
    noise = noise,
    m = 500,
    seed = 999,
    k_step = k_step,
    k_senB = k_senB
  )
  
  obs_index <- match(
    obs_variety$DAS,
    days
  )
  
  if (anyNA(obs_index)) {
    stop(
      "At least one observed DAS value is absent from the simulation time axis.",
      call. = FALSE
    )
  }
  
  biomass_metrics <- calculate_metrics(
    observed =
      obs_variety$BiomassDW_g,
    predicted =
      ensemble$B$mean[
        obs_index
      ]
  )
  
  lai_metrics <- calculate_metrics(
    observed =
      obs_variety$LAI,
    predicted =
      ensemble$LAI$mean[
        obs_index
      ]
  )
  
  
  tibble::tibble(
    VAR = variety,
    Noise = noise_kind,
    
    P1 = params$P1,
    P2 = params$P2,
    P3 = params$P3,
    P5 = params$P5,
    TT0_LAI = params$TT0_LAI,
    P6 = params$P6,
    P7 = params$P7,
    P9 = params$P9,
    
    RMSE_B =
      biomass_metrics["RMSE"],
    
    MAE_B =
      biomass_metrics["MAE"],
    
    R2_B =
      biomass_metrics["R2"],
    
    RMSE_L =
      lai_metrics["RMSE"],
    
    MAE_L =
      lai_metrics["MAE"],
    
    R2_L =
      lai_metrics["R2"],
    
    Obj_L =
      best_lai_fit$objective,
    
    Obj_B =
      fit_biomass$objective,
    
    Conv_L =
      best_lai_fit$convergence,
    
    Conv_B =
      fit_biomass$convergence
  )
}


# ============================================================
# Predictive simulation for one variety
# ============================================================

predict_serra_variety <- function(
    variety,
    calibration_results,
    base_params,
    TMIN,
    TMAX,
    PAR,
    days,
    noise_kind = NULL,
    noise_args = list(),
    m = 500,
    seed = 999,
    k_step = 0.005,
    k_senB = 0.01
) {
  
  check_serra_dependencies()
  
  result_row <- calibration_results[
    calibration_results$VAR ==
      variety,
    ,
    drop = FALSE
  ]
  
  if (
    nrow(
      result_row
    ) == 0
  ) {
    stop(
      paste0(
        "No calibration result found for variety ",
        variety,
        "."
      ),
      call. = FALSE
    )
  }
  
  result_row <- result_row[1, ]
  
  
  params <- base_params
  
  params$P1 <- result_row$P1
  params$P2 <- result_row$P2
  params$P3 <- result_row$P3
  params$P5 <- result_row$P5
  params$TT0_LAI <- result_row$TT0_LAI
  params$P6 <- result_row$P6
  params$P7 <- result_row$P7
  params$P9 <- result_row$P9
  
  
  if (is.null(noise_kind)) {
    
    if (
      "Noise" %in%
      names(result_row)
    ) {
      
      noise_kind <- as.character(
        result_row$Noise[[1]]
      )
      
    } else {
      
      stop(
        paste(
          "'noise_kind' was not supplied and",
          "'Noise' is absent from the calibration results."
        ),
        call. = FALSE
      )
    }
  }
  
  
  noise <- build_noise_spec(
    noise_kind = noise_kind,
    noise_args = noise_args
  )
  
  
  ensemble <- simulate_ensemble_v3_senB(
    TMIN = TMIN,
    TMAX = TMAX,
    PAR = PAR,
    params = params,
    noise = noise,
    m = m,
    seed = seed,
    k_step = k_step,
    k_senB = k_senB
  )
  
  
  tibble::tibble(
    VAR = variety,
    Noise = noise_kind,
    DAS = days,
    
    B_mean =
      ensemble$B$mean,
    
    B_lo =
      ensemble$B$lo,
    
    B_hi =
      ensemble$B$hi,
    
    LAI_mean =
      ensemble$LAI$mean,
    
    LAI_lo =
      ensemble$LAI$lo,
    
    LAI_hi =
      ensemble$LAI$hi
  )
}


# ============================================================
# Observed-versus-predicted table
# ============================================================

build_comparison_table <- function(
    variety,
    prediction_table,
    observations
) {
  
  observed_variety <- observations[
    observations$VAR ==
      variety,
    ,
    drop = FALSE
  ]
  
  observed_variety <- observed_variety[
    ,
    c(
      "VAR",
      "Stage",
      "DAS",
      "LAI",
      "BiomassDW_g"
    ),
    drop = FALSE
  ]
  
  names(observed_variety)[
    names(observed_variety) ==
      "LAI"
  ] <- "LAI_obs"
  
  names(observed_variety)[
    names(observed_variety) ==
      "BiomassDW_g"
  ] <- "Biomass_obs"
  
  
  merged <- merge(
    prediction_table,
    observed_variety,
    by = c(
      "VAR",
      "DAS"
    ),
    all.x = TRUE,
    sort = FALSE
  )
  
  merged[
    order(
      merged$DAS
    ),
    ,
    drop = FALSE
  ]
}