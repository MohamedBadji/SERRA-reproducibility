# SERRA reproducibility repository

This repository contains the data, R code, cached numerical results, and supporting material used for the numerical analyses presented in:

**SERRA: A Stochastic Framework for Plant Biomass Growth Combining Itô--Lévy Analysis and Numerical Evaluation on Irrigated Rice**

Authors: Mohamed BADJI, Omar Ndaw FAYE, Mamadou CISS, El Hadji DEME, and Mamadou NDIAYE.

The study develops SERRA (*Stochastic Environmental Resource Response for Agriculture*), a stochastic modelling framework for plant biomass and leaf area dynamics under environmental variability.

The numerical implementation compares Gaussian and alpha-stable Lévy-type stochastic formulations and includes calibration, Monte Carlo prediction, Leave-One-Out cross-validation, global Sobol sensitivity analysis, robustness analysis, and complementary numerical diagnostics.

## Repository structure

SERRA-reproducibility/
├── code/
├── data/
├── cache/
├── results/
├── figures/
├── README.md
├── CITATION.cff
├── LICENSE
└── .gitignore

## Code

The `code/` directory contains the final R scripts used in the reproducibility workflow:

- `01_serra_simulator.R`
- `02_serra_core.R`
- `03_serra_run.R`
- `04_serra_plots.R`
- `05_serra_loo_validation.R`
- `06_serra_sobol.R`
- `07_serra_sobol_robustness.R`
- `08_serra_diagnostics.R`

## Data

The `data/` directory contains the model inputs used in the numerical study:

- `by_variety_stage_means_2023.rds`: aggregated field observations by variety and growth stage;
- `base_params_used.rds`: baseline parameter values;
- `climat_fanaye1.xlsx`: meteorological forcing used for the Fanaye simulations.

## Cache

The `cache/` directory contains computationally expensive intermediate results retained to ensure reproducibility without requiring recalculation of all calibration, cross-validation, and Sobol analyses.

## Results

The `results/` directory contains processed numerical results used for manuscript tables and diagnostics.

## Figures

The `figures/` directory contains figures generated from the final numerical workflow.

## Main numerical settings

The final numerical experiments use the following settings:

- calibration: 50 Monte Carlo trajectories for LAI and 50 for biomass;
- final goodness-of-fit evaluation: 500 Monte Carlo trajectories;
- Leave-One-Out recalibration: 200 Monte Carlo trajectories;
- Leave-One-Out held-out prediction: 300 Monte Carlo trajectories;
- predictive uncertainty: 500 Monte Carlo trajectories;
- Sobol analysis: N = 1000 with eight parameters;
- Monte Carlo sample size within each Sobol model evaluation: 150;
- Gaussian formulation: sigma_B = 0.10 and sigma_LAI = 0.05;
- Lévy-type formulation: alpha = 1.5, beta = 0, and scale = 0.20.

The main Sobol analysis uses ±20% parameter domains around the calibrated reference values. A ±25% analysis is additionally provided as a robustness assessment.

## Running the workflow

The working directory must be the repository root.

The workflow can then be executed sequentially using:

source(file.path("code", "03_serra_run.R"))
source(file.path("code", "04_serra_plots.R"))
source(file.path("code", "05_serra_loo_validation.R"))
source(file.path("code", "06_serra_sobol.R"))
source(file.path("code", "07_serra_sobol_robustness.R"))
source(file.path("code", "08_serra_diagnostics.R"))

Scripts `01_serra_simulator.R` and `02_serra_core.R` are sourced automatically by the main workflow.

Existing cached results are loaded by default when available. Expensive analyses are recalculated only when the corresponding recalculation flag is explicitly activated.

## Software

The analyses were implemented in R.

Principal R packages used by the workflow include:

- `dplyr`
- `purrr`
- `tibble`
- `readxl`
- `sensitivity`

## Reproducibility note

The continuous Itô--Lévy formulation developed in the mathematical analysis and the discrete stochastic numerical implementation serve distinct purposes.

The numerical implementation uses positive multiplicative stochastic factors generated from either Gaussian or alpha-stable distributions. The repository therefore reproduces the numerical experiments reported in the final manuscript and should not be interpreted as a direct Euler--Maruyama discretization of the analytical stochastic differential equation.

## Data availability

The repository includes the aggregated field observations and model inputs required to reproduce the reported numerical analyses.

Meteorological forcing used in the study is derived from publicly available climate information.

## Citation

Citation information is provided in `CITATION.cff`.

## License

See the `LICENSE` file for reuse conditions.

## Contact

Mohamed BADJI  
LERSTAD, Université Gaston Berger  
Saint-Louis, Senegal  
Email: badji.mohamed@ugb.edu.sn