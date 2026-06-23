# Isotonic calibration

`cal_isotonic()` fits a monotone calibration curve with
[`stats::isoreg()`](https://rdrr.io/r/stats/isoreg.html). New
probabilities are calibrated by linear interpolation. Predictions below
the training range use the leftmost fitted value; predictions above the
range use the rightmost fitted value.

## Usage

``` r
cal_isotonic(p, y)
```

## Arguments

- p:

  Numeric vector of uncalibrated probabilities in `[0, 1]`.

- y:

  Binary outcome vector coded as `0` and `1`.

## Value

A `cal_isotonic` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new
probabilities to obtain calibrated probabilities.

## Details

Ties in the training probabilities are ordered with positive labels
first before isotonic regression and then collapsed to a single fitted
value per unique probability.

Isotonic calibration estimates a nondecreasing function \\g\\ that maps
raw probabilities to calibrated event probabilities. Let \\\pi\\ be the
ordering that sorts observations by increasing \\p_i\\ and, for equal
\\p_i\\, decreasing \\y_i\\. Thus positive labels precede negative
labels within a tied probability value. The fitted values solve the
projection problem

\$\$\min\_{m_1 \le \cdots \le m_n} \sum\_{i = 1}^n (y\_{\pi(i)} -
m_i)^2.\$\$

The implementation uses
[`stats::isoreg()`](https://rdrr.io/r/stats/isoreg.html) for the
constrained least-squares problem and clips the fitted values to
`[0, 1]`. The label vector must contain at least one `0` and one `1`.

Prediction uses linear interpolation between the unique training
probabilities and their fitted values. If a new probability is below the
smallest training value, prediction returns the leftmost fitted value.
If it is above the largest training value, prediction returns the
rightmost fitted value. Training ties are collapsed to one fitted value
per unique probability after the isotonic fit by averaging the fitted
values within each tied group. If the training data contain a single
unique probability, prediction is the resulting constant fitted value.
The fitted object stores the unique probabilities in `x_thresholds`, the
collapsed fitted values in `y_calibrated`, the
[`stats::isoreg()`](https://rdrr.io/r/stats/isoreg.html) object in
`fit`, and the original call. Prediction uses
`stats::approx(method = "linear")` with constant extrapolation at the
two endpoints, so the package prediction rule is the interpolated
monotone curve rather than the unmodified PAVA step function.

## References

Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
accurate multiclass probability estimates. Proceedings of the Eighth ACM
SIGKDD International Conference on Knowledge Discovery and Data Mining.
<doi:10.1145/775047.775151>.

## Examples

``` r
set.seed(4)
calibration <- data.frame(raw_p = sort(stats::runif(120))) |>
  dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))

fit <- cal_isotonic(calibration$raw_p, calibration$y)

calibration |>
  dplyr::mutate(calibrated = predict(fit, raw_p)) |>
  dplyr::summarise(
    raw_ece = ece(raw_p, y, bins = 10),
    calibrated_ece = ece(calibrated, y, bins = 10)
  )
#>      raw_ece calibrated_ece
#> 1 0.06602408              0
```
