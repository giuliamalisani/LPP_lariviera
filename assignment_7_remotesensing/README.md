# Assignment 7 — Remote sensing in Bregaglia

This folder contains the complete data, analysis code, model comparisons, and
prediction outputs for the Bregaglia land-cover Random Forest assignment.

## Main analysis

- `practical_7.Rmd` — primary workflow using Sentinel-2 spectral predictors,
  DHM25 terrain variables, spatial block validation, class weighting, and
  Random Forest parameter tuning.
- `bregaglia_rf_model_comparison.R` — compares local, Switzerland-wide, and
  hybrid national/local models on the same held-out Bregaglia spatial fold.
- `bregaglia_rf_hierarchy_analysis.R` — compares the detailed nine-class model,
  post-hoc aggregation to six classes, and a model retrained directly on the
  six coarse classes.
- `bregaglia_rf_experiments.R` — temporal and neighbourhood-context feature
  experiments.
- `bregaglia_rf_dem_experiment.R` — isolated DEM feature experiment.

## Principal results

All comparable model scores use the same 367-point spatial Bregaglia holdout.

| Model | Overall accuracy | Balanced accuracy | Macro-F1 |
|---|---:|---:|---:|
| Basic Bregaglia, 10 bands | 56.7% | 47.1% | 45.7% |
| Tuned Bregaglia terrain | 67.8% | 65.9% | 65.2% |
| Hybrid Switzerland + Bregaglia terrain | 68.9% | 66.2% | 65.3% |
| Fine model evaluated after six-class aggregation | 78.2% | 76.3% | 75.7% |
| Model retrained directly on six classes | 76.6% | 72.8% | 73.1% |

Artificial/amenity surfaces remain separate from the existing Static/bare
class in the six-class legend.

## Prediction outputs

- `landcover_RF_SSPCH_Bregaglia_2025.tif` — detailed nine-class prediction.
- `landcover_RF_Bregaglia_2025_posthoc_6class.tif` — detailed model recomposed
  into six higher-level classes.
- `landcover_RF_Bregaglia_2025_coarse_6class.tif` — prediction from a model
  trained directly on six higher-level classes.

## Result tables

- `bregaglia_rf_model_comparison.csv`
- `bregaglia_rf_model_comparison_by_class.csv`
- `bregaglia_rf_hierarchy_comparison.csv`
- `bregaglia_rf_hierarchy_by_class.csv`
- `confusion_fine_9class.csv`
- `confusion_posthoc_6class.csv`
- `confusion_retrained_6class.csv`

## Software

The R workflows require `terra`, `ranger`, `jsonlite`, and `knitr`. Large
rasters and GeoPackages are stored with Git LFS. Run `git lfs pull` after
cloning to download the full input data.

