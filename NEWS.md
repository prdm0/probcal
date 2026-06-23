# probcal 0.1.1

- Renamed package to probcal.

# probcal 0.1.0

- CRAN submission version, when the package was still named calibratr.

# probcal 0.0.0.9000

- Added multiclass calibration. `cal_temperature()` and `cal_cv()` accept a logit or probability matrix, and new constructors `cal_vector_scaling()`, `cal_dirichlet()`, and `cal_ovr()` cover vector scaling, Dirichlet calibration, and one-vs-rest calibration.
- `ece()`, `mce()`, `ace()`, and `reliability_diagram()` accept a probability matrix with a `type` argument for classwise or top-label confidence evaluation.
- Added `mmce()`, a binning-free Maximum Mean Calibration Error metric for binary and multiclass predictions.
- Added `inst/CITATION` so users can cite the package with `citation("probcal")`.
- Added applied workflow, calibrator selection, and numerical validation vignettes.
- `print()` and `summary()` respect `options(probcal.emoji = FALSE)` to suppress the decorative glyph in console output.
- `reliability_diagram()` now reports ECE in the subtitle by default and can use either count-scaled or fixed-size points.
- Added optional development validation tests against Python `netcal` and R `betacal`.
- Initial development version with binary calibration methods, calibration metrics, reliability diagrams, and out-of-fold calibration support.
