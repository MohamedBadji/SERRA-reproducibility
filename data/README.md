# SERRA input data

This directory contains the input data used by the final numerical SERRA workflow.

## Files

### `by_variety_stage_means_2023.rds`

Aggregated field observations used for model calibration, validation, and comparison.

The file contains observations summarized by rice variety and growth stage, including:

- variety identifier;
- growth stage;
- days after sowing;
- leaf area index;
- above-ground dry biomass.

These are aggregated analytical data used in the manuscript and do not represent the complete raw field dataset.

### `base_params_used.rds`

Baseline model parameters used to initialize and support the numerical simulations.

Some parameters are subsequently calibrated, fixed from observations, or derived from thermal-time information according to the workflow implemented in `code/02_serra_core.R`.

### `climat_fanaye1.xlsx`

Meteorological forcing used for the Fanaye simulations.

The workflow uses:

- days after sowing;
- minimum temperature;
- maximum temperature;
- photosynthetically active radiation.

These variables are used to calculate thermal time and drive biomass and LAI simulations.

## Data use

The files in this directory are provided to support reproducibility of the numerical analyses reported in the associated SERRA study.

Users should cite the associated article and repository when reusing these materials.