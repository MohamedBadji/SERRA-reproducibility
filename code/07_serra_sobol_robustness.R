# ============================================================
# SERRA Sobol robustness analysis
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#   03_serra_run.R
#   06_serra_sobol.R
#
# This file performs:
# - Sobol robustness analysis using ±25% parameter domains;
# - reuse of existing robustness-analysis caches;
# - optional recalculation when explicitly requested;
# - assembly of the final robustness results table.
# ============================================================


# ============================================================
# Source main Sobol workflow
# ============================================================

source(file.path("code", "06_serra_sobol.R"))

# ============================================================
# Robustness-analysis settings
# ============================================================

SOBOL_ROBUSTNESS_MARGIN <- 0.25

FORCE_SOBOL_ROBUSTNESS_RECALCULATION <- FALSE


# ============================================================
# Existing robustness cache files
# ============================================================

SOBOL_V1_G_25_FILE <- file.path(
  CACHE_DIR,
  "sob_V1_G_25.rds"
)

SOBOL_V1_L_25_FILE <- file.path(
  CACHE_DIR,
  "sob_V1_L_25.rds"
)

SOBOL_V2_G_25_FILE <- file.path(
  CACHE_DIR,
  "sob_V2_G_25.rds"
)

SOBOL_V2_L_25_FILE <- file.path(
  CACHE_DIR,
  "sob_V2_L_25.rds"
)


# ============================================================
# Robustness-analysis cache helper
# ============================================================

load_or_run_sobol_robustness <- function(
    cache_file,
    variety,
    noise_kind,
    force_recalculation =
      FORCE_SOBOL_ROBUSTNESS_RECALCULATION
) {
  
  if (
    file.exists(cache_file) &&
    !force_recalculation
  ) {
    
    message(
      "Loading cached Sobol robustness analysis: ",
      cache_file
    )
    
    return(
      readRDS(
        cache_file
      )
    )
  }
  
  
  message(
    "Running Sobol robustness analysis for ",
    variety,
    " / ",
    noise_kind,
    " with ±25% parameter domains."
  )
  
  
  result <- run_sobol_analysis(
    variety = variety,
    noise_kind = noise_kind,
    n = SOBOL_N,
    margin = SOBOL_ROBUSTNESS_MARGIN,
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
# Load or run the four robustness analyses
# ============================================================

sobol_V1_gaussian_25 <- load_or_run_sobol_robustness(
  cache_file =
    SOBOL_V1_G_25_FILE,
  variety =
    "V1",
  noise_kind =
    "gaussian"
)


sobol_V1_levy_25 <- load_or_run_sobol_robustness(
  cache_file =
    SOBOL_V1_L_25_FILE,
  variety =
    "V1",
  noise_kind =
    "levy"
)


sobol_V2_gaussian_25 <- load_or_run_sobol_robustness(
  cache_file =
    SOBOL_V2_G_25_FILE,
  variety =
    "V2",
  noise_kind =
    "gaussian"
)


sobol_V2_levy_25 <- load_or_run_sobol_robustness(
  cache_file =
    SOBOL_V2_L_25_FILE,
  variety =
    "V2",
  noise_kind =
    "levy"
)


# ============================================================
# Combined robustness results
# ============================================================

sobol_robustness_results <- dplyr::bind_rows(
  sobol_V1_gaussian_25$sobol,
  sobol_V1_levy_25$sobol,
  sobol_V2_gaussian_25$sobol,
  sobol_V2_levy_25$sobol
) |>
  dplyr::arrange(
    VAR,
    Noise,
    param
  )


# ============================================================
# Combined robustness parameter domains
# ============================================================

sobol_robustness_bounds <- dplyr::bind_rows(
  
  sobol_V1_gaussian_25$bounds |>
    dplyr::mutate(
      VAR = "V1",
      Noise = "gaussian"
    ),
  
  sobol_V1_levy_25$bounds |>
    dplyr::mutate(
      VAR = "V1",
      Noise = "levy"
    ),
  
  sobol_V2_gaussian_25$bounds |>
    dplyr::mutate(
      VAR = "V2",
      Noise = "gaussian"
    ),
  
  sobol_V2_levy_25$bounds |>
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


# ============================================================
# Save reusable outputs
# ============================================================

saveRDS(
  sobol_robustness_results,
  file.path(
    CACHE_DIR,
    "sobol_robustness_results.rds"
  )
)


saveRDS(
  sobol_robustness_bounds,
  file.path(
    CACHE_DIR,
    "sobol_robustness_bounds.rds"
  )
)


# ============================================================
# Console output
# ============================================================

print(
  sobol_robustness_results
)

print(
  sobol_robustness_bounds
)