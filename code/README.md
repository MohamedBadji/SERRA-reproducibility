# SERRA R code

This directory contains the final R scripts used to reproduce the numerical analyses associated with the SERRA framework.

## Execution order

The scripts are organized as follows:

1. `01_serra_simulator.R`  
   Core stochastic simulation functions.

2. `02_serra_core.R`  
   Thermal-time construction, calibration objectives, goodness-of-fit metrics, sequential calibration, prediction, and comparison utilities.

3. `03_serra_run.R`  
   Main analysis workflow. Loads model inputs, reuses cached calibrations when available, generates final Monte Carlo predictions, and assembles calibration and comparison results.

4. `04_serra_plots.R`  
   Generates final graphical outputs and predictive-uncertainty summaries.

5. `05_serra_loo_validation.R`  
   Leave-One-Out cross-validation with recalibration and held-out prediction.

6. `06_serra_sobol.R`  
   Main Sobol global sensitivity analysis using ±20% parameter domains.

7. `07_serra_sobol_robustness.R`  
   Robustness assessment of the Sobol analysis using ±25% parameter domains.

8. `08_serra_diagnostics.R`  
   Residual-based AIC/BIC comparison and numerical terminal-stability diagnostic.

## Running the workflow

Run the scripts from the repository root, not from inside the `code/` directory.

Execute:

source(file.path("code", "03_serra_run.R"))
source(file.path("code", "04_serra_plots.R"))
source(file.path("code", "05_serra_loo_validation.R"))
source(file.path("code", "06_serra_sobol.R"))
source(file.path("code", "07_serra_sobol_robustness.R"))
source(file.path("code", "08_serra_diagnostics.R"))

Scripts `01_serra_simulator.R` and `02_serra_core.R` are sourced automatically by the main workflow.

Cached results are loaded by default when available. Expensive recalculations require explicit activation of the corresponding recalculation flags.