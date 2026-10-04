# ============================================================
# SERRA Leave-One-Out cross-validation
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#   03_serra_run.R
#
# This file performs:
# - Leave-One-Out recalibration;
# - held-out prediction;
# - LOO performance summaries;
# - optional reuse of cached LOO results.
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
# Cache location
# ============================================================

CACHE_DIR <- "cache"

if (!dir.exists(CACHE_DIR)) {
  dir.create(
    CACHE_DIR,
    recursive = TRUE
  )
}


LOO_CACHE_FILE <- file.path(
  CACHE_DIR,
  "loo_results_full.rds"
)


FORCE_LOO_RECALCULATION <- FALSE


# ============================================================
# Single-variety LOO validation
# ============================================================

run_loo_variety <- function(
    variety,
    noise_kind,
    noise_args
) {
  
  obs_variety <- observations |>
    dplyr::filter(
      VAR == variety
    ) |>
    dplyr::arrange(
      DAS
    )
  
  
  if (nrow(obs_variety) < 2) {
    stop(
      paste0(
        "At least two observations are required for LOO validation: ",
        variety
      ),
      call. = FALSE
    )
  }
  
  
  purrr::map_dfr(
    seq_len(
      nrow(obs_variety)
    ),
    
    function(i) {
      
      training_data <- obs_variety[
        -i,
        ,
        drop = FALSE
      ]
      
      test_data <- obs_variety[
        i,
        ,
        drop = FALSE
      ]
      
      
      # Recalibration excluding the held-out observation
      
      loo_fit <- fit_serra_by_variety(
        variety = variety,
        base_params = base_params,
        thermal_time_table = thermal_time_table,
        observations = training_data,
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        days = days,
        noise_kind = noise_kind,
        noise_args = noise_args,
        m_ens_lai = 200,
        m_ens_biomass = 200,
        buffer_fraction = 0.10,
        k_step = 0.005,
        k_senB = 0.01
      )
      
      
      # Prediction of the held-out observation
      
      prediction <- predict_serra_variety(
        variety = variety,
        calibration_results = loo_fit,
        base_params = base_params,
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        days = days,
        noise_kind = noise_kind,
        noise_args = noise_args,
        m = 300,
        seed = 999,
        k_step = 0.005,
        k_senB = 0.01
      )
      
      
      held_out_prediction <- prediction |>
        dplyr::filter(
          DAS == test_data$DAS
        )
      
      
      tibble::tibble(
        VAR = variety,
        Noise = noise_kind,
        Stage = test_data$Stage,
        DAS = test_data$DAS,
        
        obs_B =
          test_data$BiomassDW_g,
        
        pred_B =
          held_out_prediction$B_mean,
        
        obs_LAI =
          test_data$LAI,
        
        pred_LAI =
          held_out_prediction$LAI_mean
      )
    }
  )
}


# ============================================================
# Complete LOO analysis
# ============================================================

run_full_loo <- function() {
  
  dplyr::bind_rows(
    
    run_loo_variety(
      variety = "V1",
      noise_kind = "gaussian",
      noise_args = list(
        sigma1 = 0.10,
        sigma2 = 0.05
      )
    ),
    
    run_loo_variety(
      variety = "V2",
      noise_kind = "gaussian",
      noise_args = list(
        sigma1 = 0.10,
        sigma2 = 0.05
      )
    ),
    
    run_loo_variety(
      variety = "V1",
      noise_kind = "levy",
      noise_args = list(
        alpha = 1.5,
        beta = 0,
        scale = 0.20
      )
    ),
    
    run_loo_variety(
      variety = "V2",
      noise_kind = "levy",
      noise_args = list(
        alpha = 1.5,
        beta = 0,
        scale = 0.20
      )
    )
  )
}


# ============================================================
# Load or compute LOO results
# ============================================================

if (
  file.exists(LOO_CACHE_FILE) &&
  !FORCE_LOO_RECALCULATION
) {
  
  message(
    "Loading cached LOO results: ",
    LOO_CACHE_FILE
  )
  
  loo_results <- readRDS(
    LOO_CACHE_FILE
  )
  
} else {
  
  message(
    "Running complete Leave-One-Out validation."
  )
  
  loo_results <- run_full_loo()
  
  saveRDS(
    loo_results,
    LOO_CACHE_FILE
  )
}


# ============================================================
# LOO performance summary
# ============================================================

loo_summary <- loo_results |>
  dplyr::group_by(
    VAR,
    Noise
  ) |>
  dplyr::summarise(
    
    RMSE_LOO_B =
      sqrt(
        mean(
          (obs_B - pred_B)^2,
          na.rm = TRUE
        )
      ),
    
    MAE_LOO_B =
      mean(
        abs(
          obs_B - pred_B
        ),
        na.rm = TRUE
      ),
    
    RMSE_LOO_LAI =
      sqrt(
        mean(
          (obs_LAI - pred_LAI)^2,
          na.rm = TRUE
        )
      ),
    
    MAE_LOO_LAI =
      mean(
        abs(
          obs_LAI - pred_LAI
        ),
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  dplyr::arrange(
    VAR,
    Noise
  )


# ============================================================
# Save summary
# ============================================================

saveRDS(
  loo_summary,
  file.path(
    CACHE_DIR,
    "loo_summary.rds"
  )
)


# ============================================================
# Console output
# ============================================================

print(
  loo_summary
)