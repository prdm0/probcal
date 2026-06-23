# Histogram binning calibration

`cal_histogram()` partitions `[0, 1]` into bins and replaces each
probability with the empirical event frequency in its bin. Equal-width
bins use fixed intervals. Equal-frequency bins use sample quantiles as
break points.

## Usage

``` r
cal_histogram(p, y, bins = 10, strategy = c("equal_width", "equal_freq"))
```

## Arguments

- p:

  Numeric vector of uncalibrated probabilities in `[0, 1]`.

- y:

  Binary outcome vector coded as `0` and `1`.

- bins:

  Number of bins. Must be a single positive integer.

- strategy:

  Binning strategy. Use `"equal_width"` for fixed-width bins or
  `"equal_freq"` for quantile bins.

## Value

A `cal_histogram` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new
probabilities to obtain calibrated probabilities.

## Details

Empty training bins inherit the empirical rate from the nearest
non-empty bin. This makes prediction defined over the whole interval
`[0, 1]`.

Histogram binning estimates a piecewise constant calibration map. Given
distinct break points \\0 = b_0 \< b_1 \< \cdots \< b_J = 1\\, the
implementation uses left-closed bins. For \\j \< J\\,

\$\$I_j = \\i: b\_{j-1} \le p_i \< b_j\\,\$\$

and the last bin is

\$\$I_J = \\i: b\_{J-1} \le p_i \le b_J\\.\$\$

The fitted value for a non-empty bin is the empirical event frequency,

\$\$\hat q_j = \frac{1}{n_j}\sum\_{i \in I_j} y_i, \quad n_j =
\|I_j\|.\$\$

A new probability receives the fitted value of the bin into which it
falls. Values exactly on an internal break point are assigned to the bin
that starts at that break point; the value `1` is assigned to the last
bin.

With `strategy = "equal_width"`, the break points are equally spaced on
`[0, 1]`, so \\J = B\\ when `bins = B`. With `strategy = "equal_freq"`,
provisional break points are

\$\$b_j = Q_8(j / B), \quad j = 0, \ldots, B,\$\$

where \\Q_8\\ is the sample quantile computed by
`stats::quantile(type = 8)`. The first and last break points are then
forced to `0` and `1`. Duplicated break points are removed, so the
actual number of bins \\J\\ can be smaller than `bins`. Empty bins are
assigned the value of the nearest non-empty bin by bin index; if an
empty bin is equally close to two non-empty bins, the lower-index
non-empty bin is used. If no non-empty bin is available, the global
event rate is used as a fallback.

The returned object stores the requested `bins`, the realized
`actual_bins`, `strategy`, `breaks`, per-bin fitted values in
`bin_values`, training `counts`, `global_rate`, and the original call.

## References

Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
accurate multiclass probability estimates. Proceedings of the Eighth ACM
SIGKDD International Conference on Knowledge Discovery and Data Mining.
<doi:10.1145/775047.775151>.

## Examples

``` r
set.seed(5)
calibration <- data.frame(raw_p = stats::runif(120)) |>
  dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))

fit <- cal_histogram(calibration$raw_p, calibration$y, bins = 5)

calibration |>
  dplyr::mutate(calibrated = predict(fit, raw_p)) |>
  dplyr::summarise(
    raw_ece = ece(raw_p, y, bins = 5),
    calibrated_ece = ece(calibrated, y, bins = 5)
  )
#>      raw_ece calibrated_ece
#> 1 0.04634435              0
```
