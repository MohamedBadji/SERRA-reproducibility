# ============================================================
# SERRA simulator
# Gaussian and alpha-stable Levy-type stochastic formulations
#
# This file defines:
# - stochastic-noise utilities;
# - positive multiplicative perturbations;
# - single-trajectory simulation;
# - Monte Carlo ensemble simulation.
# ============================================================


# ============================================================
# Stochastic-noise utilities
# ============================================================

winsorize <- function(x, q = 0.995) {
  
  qs <- stats::quantile(
    x,
    probs = c(1 - q, q),
    na.rm = TRUE
  )
  
  x <- pmax(x, qs[1])
  x <- pmin(x, qs[2])
  
  x
}


cms_alpha_stable <- function(
    alpha = 1.5,
    beta = 0,
    n = 1L
) {
  
  U <- stats::runif(
    n,
    min = -pi / 2,
    max = pi / 2
  )
  
  W <- stats::rexp(
    n,
    rate = 1
  )
  
  a <- alpha
  b <- beta
  
  if (abs(a - 1) < 1e-8) {
    a <- 1 + 1e-6
  }
  
  phi <- (1 / a) *
    atan(
      b * tan(pi * a / 2)
    )
  
  S <- (
    1 +
      b^2 *
      tan(pi * a / 2)^2
  )^(1 / (2 * a))
  
  S *
    sin(a * (U + phi)) /
    cos(U)^(1 / a) *
    (
      cos(U - a * (U + phi)) / W
    )^((1 - a) / a)
}


make_noise <- function(
    kind = c("gaussian", "levy"),
    sigma1 = 0.10,
    sigma2 = 0.05,
    alpha = 1.5,
    beta = 0,
    scale = 0.20
) {
  
  kind <- match.arg(kind)
  
  list(
    kind = kind,
    sigma1 = sigma1,
    sigma2 = sigma2,
    alpha = alpha,
    beta = beta,
    scale = scale
  )
}


# ============================================================
# Positive multiplicative perturbations
# ============================================================

pos_mult_gaussian <- function(
    sigma,
    n = 1L
) {
  
  exp(
    stats::rnorm(
      n,
      mean = -0.5 * sigma^2,
      sd = sigma
    )
  )
}


pos_mult_levy <- function(
    alpha = 1.5,
    beta = 0,
    scale = 0.20,
    n = 1L,
    q_winsor = 0.995
) {
  
  z <- cms_alpha_stable(
    alpha = alpha,
    beta = beta,
    n = n
  )
  
  z <- winsorize(
    z,
    q = q_winsor
  )
  
  exp(
    scale * z -
      0.5 * scale^2
  )
}


cap_mult <- function(
    g,
    g_min = 0.05,
    g_max = 20
) {
  
  pmin(
    pmax(g, g_min),
    g_max
  )
}


# ============================================================
# Deterministic helper functions
# ============================================================

smooth_step <- function(
    x,
    k = 0.08
) {
  
  1 / (
    1 +
      exp(-k * x)
  )
}


clip_projection <- function(
    TT,
    B,
    LAI,
    P5
) {
  
  TT <- max(TT, 0)
  B <- max(B, 0)
  
  LAI <- min(
    max(LAI, 0),
    P5
  )
  
  c(
    TT = TT,
    B = B,
    LAI = LAI
  )
}


# ============================================================
# Monte Carlo summaries
# ============================================================

ic95_from_mat <- function(
    mat_rows_are_runs
) {
  
  list(
    
    mean = apply(
      mat_rows_are_runs,
      2,
      mean,
      na.rm = TRUE
    ),
    
    lo = apply(
      mat_rows_are_runs,
      2,
      stats::quantile,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    hi = apply(
      mat_rows_are_runs,
      2,
      stats::quantile,
      probs = 0.975,
      na.rm = TRUE
    )
    
  )
}


# ============================================================
# Single-trajectory simulation
# ============================================================

simulate_once_v3_senB <- function(
    TMIN,
    TMAX,
    PAR,
    params,
    noise,
    X0 = c(
      TT = 0,
      B = 0,
      LAI = 0.05
    ),
    project = TRUE,
    seed = NULL,
    g_min = 0.05,
    g_max = 20,
    g_max_LAI = 2.0,
    k_step = 0.005,
    k_senB = 0.01
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- length(PAR)
  
  stopifnot(
    length(TMIN) == n,
    length(TMAX) == n
  )
  
  TT <- numeric(n)
  B <- numeric(n)
  LAI <- numeric(n)
  
  TT[1] <- X0["TT"]
  B[1] <- X0["B"]
  LAI[1] <- X0["LAI"]
  
  TT0 <- if (!is.null(params$TT0_LAI)) {
    params$TT0_LAI
  } else {
    0
  }
  
  P9 <- if (!is.null(params$P9)) {
    params$P9
  } else {
    0.001
  }
  
  noise_kind <- if (!is.null(noise$kind)) {
    noise$kind
  } else {
    noise$type
  }
  
  if (is.null(noise_kind)) {
    stop(
      "The noise specification must contain either 'kind' or 'type'."
    )
  }
  
  
  for (j in seq_len(n - 1)) {
    
    U1 <- TMIN[j]
    U2 <- TMAX[j]
    U3 <- PAR[j]
    
    
    # Thermal-time accumulation
    
    dTT <- max(
      (U1 + U2) / 2 -
        params$P4,
      0
    )
    
    TT_next <- TT[j] + dTT
    
    
    # Stochastic multiplicative factors
    
    if (noise_kind == "gaussian") {
      
      g1 <- pos_mult_gaussian(
        sigma = noise$sigma1,
        n = 1
      )
      
      g2 <- pos_mult_gaussian(
        sigma = noise$sigma2,
        n = 1
      )
      
    } else if (noise_kind == "levy") {
      
      g1 <- pos_mult_levy(
        alpha = noise$alpha,
        beta = noise$beta,
        scale = noise$scale,
        n = 1,
        q_winsor = 0.995
      )
      
      g2 <- pos_mult_levy(
        alpha = noise$alpha,
        beta = noise$beta,
        scale = noise$scale,
        n = 1,
        q_winsor = 0.995
      )
      
    } else {
      
      stop(
        "Unknown stochastic formulation."
      )
    }
    
    
    g1 <- cap_mult(
      as.numeric(g1),
      g_min = g_min,
      g_max = g_max
    )
    
    g2 <- cap_mult(
      as.numeric(g2),
      g_min = g_min,
      g_max = g_max_LAI
    )
    
    
    # Biomass production
    
    growthB <- if (
      TT[j] <= params$P6
    ) {
      
      params$P2 *
        U3 *
        (
          1 -
            exp(
              -params$P1 *
                LAI[j]
            )
        ) *
        g1
      
    } else {
      
      0
    }
    
    
    # Biomass senescence
    
    sB <- smooth_step(
      TT[j] - params$P6,
      k = k_senB
    )
    
    decayB <- sB *
      P9 *
      B[j] *
      dTT
    
    
    B_next <- max(
      B[j] +
        growthB -
        decayB,
      0
    )
    
    
    # Leaf-area development
    
    sTT <- smooth_step(
      TT[j] - TT0,
      k = k_step
    )
    
    
    dLAI <- if (
      TT[j] <= params$P7
    ) {
      
      sTT *
        params$P3 *
        dTT *
        LAI[j] *
        max(
          params$P5 -
            LAI[j],
          0
        ) *
        g2
      
    } else {
      
      0
    }
    
    
    LAI_next <- min(
      max(
        LAI[j] +
          dLAI,
        0
      ),
      params$P5
    )
    
    
    # State-space projection
    
    if (project) {
      
      projected_state <- clip_projection(
        TT = TT_next,
        B = B_next,
        LAI = LAI_next,
        P5 = params$P5
      )
      
      TT_next <- projected_state["TT"]
      B_next <- projected_state["B"]
      LAI_next <- projected_state["LAI"]
    }
    
    
    TT[j + 1] <- TT_next
    B[j + 1] <- B_next
    LAI[j + 1] <- LAI_next
  }
  
  
  list(
    TT = TT,
    B = B,
    LAI = LAI
  )
}


# ============================================================
# Monte Carlo ensemble simulation
# ============================================================

simulate_ensemble_v3_senB <- function(
    TMIN,
    TMAX,
    PAR,
    params,
    noise,
    m = 300,
    seed = 123,
    g_min = 0.05,
    g_max = 20,
    g_max_LAI = 2.0,
    k_step = 0.005,
    k_senB = 0.01
) {
  
  set.seed(seed)
  
  n <- length(PAR)
  
  Bmat <- matrix(
    NA_real_,
    nrow = m,
    ncol = n
  )
  
  LAImat <- matrix(
    NA_real_,
    nrow = m,
    ncol = n
  )
  
  
  for (r in seq_len(m)) {
    
    sim <- simulate_once_v3_senB(
      TMIN = TMIN,
      TMAX = TMAX,
      PAR = PAR,
      params = params,
      noise = noise,
      project = TRUE,
      seed = NULL,
      g_min = g_min,
      g_max = g_max,
      g_max_LAI = g_max_LAI,
      k_step = k_step,
      k_senB = k_senB
    )
    
    Bmat[r, ] <- sim$B
    LAImat[r, ] <- sim$LAI
  }
  
  
  list(
    
    B = ic95_from_mat(
      Bmat
    ),
    
    LAI = ic95_from_mat(
      LAImat
    )
    
  )
}