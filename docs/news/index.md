# Changelog

## probcal 0.1.1

- Renamed package to probcal.

## probcal 0.1.0

- CRAN submission version, when the package was still named calibratr.

## probcal 0.0.0.9000

- Added multiclass calibration.
  [`cal_temperature()`](https://prdm0.github.io/probcal/reference/cal_temperature.md)
  and [`cal_cv()`](https://prdm0.github.io/probcal/reference/cal_cv.md)
  accept a logit or probability matrix, and new constructors
  [`cal_vector_scaling()`](https://prdm0.github.io/probcal/reference/cal_vector_scaling.md),
  [`cal_dirichlet()`](https://prdm0.github.io/probcal/reference/cal_dirichlet.md),
  and
  [`cal_ovr()`](https://prdm0.github.io/probcal/reference/cal_ovr.md)
  cover vector scaling, Dirichlet calibration, and one-vs-rest
  calibration.
- [`ece()`](https://prdm0.github.io/probcal/reference/ece.md),
  [`mce()`](https://prdm0.github.io/probcal/reference/mce.md),
  [`ace()`](https://prdm0.github.io/probcal/reference/ace.md), and
  [`reliability_diagram()`](https://prdm0.github.io/probcal/reference/reliability_diagram.md)
  accept a probability matrix with a `type` argument for classwise or
  top-label confidence evaluation.
- Added [`mmce()`](https://prdm0.github.io/probcal/reference/mmce.md), a
  binning-free Maximum Mean Calibration Error metric for binary and
  multiclass predictions.
- Added `inst/CITATION` so users can cite the package with
  `citation("probcal")`.
- Added applied workflow, calibrator selection, and numerical validation
  vignettes.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) respect
  `options(probcal.emoji = FALSE)` to suppress the decorative glyph in
  console output.
- [`reliability_diagram()`](https://prdm0.github.io/probcal/reference/reliability_diagram.md)
  now reports ECE in the subtitle by default and can use either
  count-scaled or fixed-size points.
- Added optional development validation tests against Python `netcal`
  and R `betacal`.
- Initial development version with binary calibration methods,
  calibration metrics, reliability diagrams, and out-of-fold calibration
  support.
