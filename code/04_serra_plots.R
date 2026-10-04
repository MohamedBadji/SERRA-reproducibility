# ============================================================
# SERRA publication figures
#
# Requires:
#   01_serra_simulator.R
#   02_serra_core.R
#   03_serra_run.R
#
# This file generates:
# - final model-fit figure;
# - residual diagnostics;
# - predictive-uncertainty summaries;
# - uncertainty-width figure;
# - combined residual figure;
# - optional sample-path illustrations.
# ============================================================


# ============================================================
# Packages
# ============================================================

required_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
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

FIGURE_DIR <- "figures"

if (!dir.exists(FIGURE_DIR)) {
  dir.create(
    FIGURE_DIR,
    recursive = TRUE
  )
}


# ============================================================
# Publication theme
# ============================================================

theme_serra <- function(
    base_size = 13
) {
  
  ggplot2::theme_classic(
    base_size = base_size
  ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
}


# ============================================================
# Publication labels
# ============================================================

prediction_plot_data <- prediction_results |>
  dplyr::mutate(
    Model = dplyr::case_when(
      Noise == "gaussian" ~ "Gaussian",
      Noise == "levy" ~ "Levy",
      TRUE ~ Noise
    )
  ) |>
  dplyr::left_join(
    observations |>
      dplyr::select(
        VAR,
        DAS,
        LAI_obs = LAI,
        B_obs = BiomassDW_g
      ),
    by = c(
      "VAR",
      "DAS"
    )
  )


# ============================================================
# Model colors
# ============================================================

model_colors <- c(
  "Gaussian" = "#0072B2",
  "Levy" = "#D55E00"
)


# ============================================================
# Biomass fit panel
# ============================================================

plot_biomass_fit <- function(
    data,
    variety,
    panel_label = NULL
) {
  
  plot_data <- data |>
    dplyr::filter(
      VAR == variety
    )
  
  panel_title <- if (is.null(panel_label)) {
    paste(
      "Biomass -",
      variety
    )
  } else {
    paste(
      panel_label,
      "Biomass -",
      variety
    )
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = DAS,
      y = B_mean,
      color = Model,
      fill = Model
    )
  ) +
    
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = B_lo,
        ymax = B_hi
      ),
      alpha = 0.20,
      colour = NA
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.1
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        y = B_obs
      ),
      shape = 21,
      size = 2.6,
      stroke = 0.8,
      fill = "white",
      colour = "black",
      na.rm = TRUE
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors
    ) +
    
    ggplot2::scale_fill_manual(
      values = model_colors
    ) +
    
    ggplot2::labs(
      title = panel_title,
      x = "Days after sowing (DAS)",
      y = expression(
        Biomass ~ (g ~ m^{-2})
      )
    ) +
    
    theme_serra()
}


# ============================================================
# LAI fit panel
# ============================================================

plot_lai_fit <- function(
    data,
    variety,
    panel_label = NULL
) {
  
  plot_data <- data |>
    dplyr::filter(
      VAR == variety
    )
  
  panel_title <- if (is.null(panel_label)) {
    paste(
      "LAI -",
      variety
    )
  } else {
    paste(
      panel_label,
      "LAI -",
      variety
    )
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = DAS,
      y = LAI_mean,
      color = Model,
      fill = Model
    )
  ) +
    
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = LAI_lo,
        ymax = LAI_hi
      ),
      alpha = 0.20,
      colour = NA
    ) +
    
    ggplot2::geom_line(
      linewidth = 1.1
    ) +
    
    ggplot2::geom_point(
      ggplot2::aes(
        y = LAI_obs
      ),
      shape = 21,
      size = 2.6,
      stroke = 0.8,
      fill = "white",
      colour = "black",
      na.rm = TRUE
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors
    ) +
    
    ggplot2::scale_fill_manual(
      values = model_colors
    ) +
    
    ggplot2::labs(
      title = panel_title,
      x = "Days after sowing (DAS)",
      y = "Leaf Area Index"
    ) +
    
    theme_serra()
}


# ============================================================
# Main 2 x 2 model-fit figure
# ============================================================

plot_model_fit <- function(
    data,
    varieties = c(
      "V1",
      "V2"
    )
) {
  
  if (length(varieties) != 2) {
    stop(
      "Exactly two varieties are required for the 2 x 2 publication figure.",
      call. = FALSE
    )
  }
  
  p1 <- plot_biomass_fit(
    data = data,
    variety = varieties[1],
    panel_label = "(a)"
  )
  
  p2 <- plot_lai_fit(
    data = data,
    variety = varieties[1],
    panel_label = "(b)"
  )
  
  p3 <- plot_biomass_fit(
    data = data,
    variety = varieties[2],
    panel_label = "(c)"
  )
  
  p4 <- plot_lai_fit(
    data = data,
    variety = varieties[2],
    panel_label = "(d)"
  )
  
  (
    p1 | p2
  ) /
    (
      p3 | p4
    )
}


# ============================================================
# Residual table
# ============================================================

build_residual_table <- function(
    data
) {
  
  data |>
    dplyr::mutate(
      residual_B =
        B_obs -
        B_mean,
      
      residual_LAI =
        LAI_obs -
        LAI_mean
    )
}


# ============================================================
# Residual plot
# ============================================================

plot_residuals <- function(
    residual_data,
    variety,
    variable = c(
      "biomass",
      "lai"
    ),
    show_title = TRUE
) {
  
  variable <- match.arg(
    variable
  )
  
  if (variable == "biomass") {
    
    residual_column <- "residual_B"
    y_label <- "Biomass residual"
    variable_label <- "Biomass"
    
  } else {
    
    residual_column <- "residual_LAI"
    y_label <- "LAI residual"
    variable_label <- "LAI"
  }
  
  
  plot_data <- residual_data |>
    dplyr::filter(
      VAR == variety,
      !is.na(
        .data[[residual_column]]
      )
    )
  
  
  plot_title <- if (show_title) {
    
    paste(
      variety,
      "-",
      variable_label
    )
    
  } else {
    
    NULL
  }
  
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = DAS,
      y = .data[[residual_column]],
      color = Model
    )
  ) +
    
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    
    ggplot2::geom_point(
      size = 2
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors
    ) +
    
    ggplot2::labs(
      title = plot_title,
      x = "Days after sowing (DAS)",
      y = y_label
    ) +
    
    theme_serra()
}


# ============================================================
# Predictive-band widths
# ============================================================

build_uncertainty_width_table <- function(
    data
) {
  
  data |>
    dplyr::mutate(
      width_B =
        B_hi -
        B_lo,
      
      width_LAI =
        LAI_hi -
        LAI_lo
    )
}


# ============================================================
# Uncertainty-width summary
# ============================================================

summarize_uncertainty_width <- function(
    width_data
) {
  
  width_data |>
    dplyr::group_by(
      VAR,
      Model
    ) |>
    dplyr::summarise(
      median_width_B =
        stats::median(
          width_B,
          na.rm = TRUE
        ),
      
      median_width_LAI =
        stats::median(
          width_LAI,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    )
}


# ============================================================
# Biomass uncertainty-width figure
# ============================================================

plot_uncertainty_width_biomass <- function(
    width_data
) {
  
  ggplot2::ggplot(
    width_data,
    ggplot2::aes(
      x = DAS,
      y = width_B,
      color = Model
    )
  ) +
    
    ggplot2::geom_line(
      linewidth = 1
    ) +
    
    ggplot2::facet_wrap(
      ~VAR
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors
    ) +
    
    ggplot2::labs(
      x = "Days after sowing (DAS)",
      y = expression(
        "Width of 95% biomass predictive band" ~
          (g ~ m^{-2})
      )
    ) +
    
    theme_serra()
}


# ============================================================
# Sample trajectories
# ============================================================

simulate_sample_paths <- function(
    variety,
    calibration_results,
    base_params,
    TMIN,
    TMAX,
    PAR,
    days,
    noise_kind,
    noise_args,
    n_paths = 20,
    seed = 123,
    k_step = 0.005,
    k_senB = 0.01
) {
  
  result_row <- calibration_results[
    calibration_results$VAR ==
      variety,
    ,
    drop = FALSE
  ]
  
  if (nrow(result_row) == 0) {
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
  
  
  noise <- build_noise_spec(
    noise_kind = noise_kind,
    noise_args = noise_args
  )
  
  
  set.seed(
    seed
  )
  
  
  purrr::map_dfr(
    seq_len(
      n_paths
    ),
    
    function(path_id) {
      
      simulation <- simulate_once_v3_senB(
        TMIN = TMIN,
        TMAX = TMAX,
        PAR = PAR,
        params = params,
        noise = noise,
        project = TRUE,
        seed = NULL,
        k_step = k_step,
        k_senB = k_senB
      )
      
      tibble::tibble(
        VAR = variety,
        Model = ifelse(
          noise_kind == "gaussian",
          "Gaussian",
          "Levy"
        ),
        path = path_id,
        DAS = days,
        B = simulation$B,
        LAI = simulation$LAI
      )
    }
  )
}


# ============================================================
# Sample-path plot
# ============================================================

plot_sample_paths <- function(
    path_data,
    variable = c(
      "biomass",
      "lai"
    )
) {
  
  variable <- match.arg(
    variable
  )
  
  if (variable == "biomass") {
    
    y_column <- "B"
    
    y_label <- expression(
      Biomass ~ (g ~ m^{-2})
    )
    
  } else {
    
    y_column <- "LAI"
    y_label <- "Leaf Area Index"
  }
  
  
  ggplot2::ggplot(
    path_data,
    ggplot2::aes(
      x = DAS,
      y = .data[[y_column]],
      group = interaction(
        Model,
        path
      ),
      color = Model
    )
  ) +
    
    ggplot2::geom_line(
      alpha = 0.5
    ) +
    
    ggplot2::scale_color_manual(
      values = model_colors
    ) +
    
    ggplot2::labs(
      x = "Days after sowing (DAS)",
      y = y_label
    ) +
    
    theme_serra()
}


# ============================================================
# Build publication objects
# ============================================================

main_fit_figure <- plot_model_fit(
  data = prediction_plot_data,
  varieties = c(
    "V1",
    "V2"
  )
)


residual_table <- build_residual_table(
  prediction_plot_data
)


uncertainty_width_table <-
  build_uncertainty_width_table(
    prediction_plot_data
  )


uncertainty_width_summary <-
  summarize_uncertainty_width(
    uncertainty_width_table
  )


uncertainty_width_figure <-
  plot_uncertainty_width_biomass(
    uncertainty_width_table
  )


# ============================================================
# Combined residual figure
# ============================================================

residual_v1_biomass <- plot_residuals(
  residual_data = residual_table,
  variety = "V1",
  variable = "biomass"
) +
  ggplot2::theme(
    legend.position = "none"
  )


residual_v2_biomass <- plot_residuals(
  residual_data = residual_table,
  variety = "V2",
  variable = "biomass"
) +
  ggplot2::theme(
    legend.position = "none"
  )


residual_v1_lai <- plot_residuals(
  residual_data = residual_table,
  variety = "V1",
  variable = "lai"
) +
  ggplot2::theme(
    legend.position = "none"
  )


residual_v2_lai <- plot_residuals(
  residual_data = residual_table,
  variety = "V2",
  variable = "lai"
) +
  ggplot2::theme(
    legend.position = "none"
  )


residual_figure <-
  (
    residual_v1_biomass |
      residual_v2_biomass
  ) /
  (
    residual_v1_lai |
      residual_v2_lai
  )


# ============================================================
# Save publication figures
# ============================================================

ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_serra_fit.png"
  ),
  plot = main_fit_figure,
  width = 11,
  height = 7,
  dpi = 600
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_serra_fit.pdf"
  ),
  plot = main_fit_figure,
  width = 11,
  height = 7
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_uncertainty_width.png"
  ),
  plot = uncertainty_width_figure,
  width = 8,
  height = 5,
  dpi = 300
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_uncertainty_width.pdf"
  ),
  plot = uncertainty_width_figure,
  width = 8,
  height = 5
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_residuals.png"
  ),
  plot = residual_figure,
  width = 10,
  height = 8,
  dpi = 600
)


ggplot2::ggsave(
  filename = file.path(
    FIGURE_DIR,
    "fig_residuals.pdf"
  ),
  plot = residual_figure,
  width = 10,
  height = 8
)


# ============================================================
# Console summaries
# ============================================================

print(
  uncertainty_width_summary
)

print(
  main_fit_figure
)

print(
  uncertainty_width_figure
)

print(
  residual_figure
)