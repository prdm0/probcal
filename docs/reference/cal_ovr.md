# One-vs-rest multiclass calibration

`cal_ovr()` extends any binary calibrator to a multiclass problem with
the one-vs-rest reduction. For each class it fits a binary calibrator
that separates that class from the others, applies the calibrators
column by column, and renormalizes each row to sum to one. This is the
default strategy that binning methods use for multiclass calibration.

## Usage

``` r
cal_ovr(
  x,
  y,
  method = c("platt", "beta", "isotonic", "histogram", "temperature"),
  ...
)
```

## Arguments

- x:

  Numeric matrix of uncalibrated values with one row per observation and
  one column per class. For `method = "platt"`, entries may be arbitrary
  finite scores. For `"beta"`, `"isotonic"`, and `"histogram"`, entries
  must be probabilities in `[0, 1]`. For `"temperature"`, entries are
  logits.

- y:

  A factor or a vector of integer class codes in `1:K`, where `K` is the
  number of columns of `x`.

- method:

  Binary calibrator applied to each one-vs-rest problem.

- ...:

  Additional arguments passed to the binary calibrator, such as `bins`
  for `method = "histogram"`.

## Value

A `cal_ovr` object that also inherits from `cal_multiclass`. The object
stores `calibrators`, `base_method`, `k`, `levels`, `input`, and the
original call. Use [`predict()`](https://rdrr.io/r/stats/predict.html)
with a new score matrix to obtain a numeric matrix of calibrated
probabilities whose rows sum to one.

## Details

The columns of `x` are the per-class uncalibrated values. Use scores or
probabilities for `method = "platt"`, probabilities for `"beta"`,
`"isotonic"`, and `"histogram"`, and binary one-vs-rest logits for
`"temperature"`. Rows of `x` are not required to sum to one. Every class
must appear at least once in `y`, because each one-vs-rest problem needs
both labels.

For \\K\\ classes, column \\k\\ of `x` corresponds to integer class code
\\k\\; if `y` is a factor, column \\k\\ corresponds to `levels(y)[k]`.
One-vs-rest calibration creates \\K\\ binary labels,

\$\$y_i^{(k)} = \mathbf{1}\\y_i = k\\, \quad k = 1, \ldots, K.\$\$

A separate binary calibrator \\f_k\\ is fitted to column \\k\\ of `x`
and the binary labels \\y_i^{(k)}\\. On new data, the classwise
calibrated scores are

\$\$r\_{ik} = f_k(x\_{ik}).\$\$

Because the \\K\\ binary calibrators are fitted independently, the row
sums of \\r\_{ik}\\ need not equal one. Let \\S_i = \sum\_{\ell = 1}^K
r\_{i\ell}\\. If \\S_i\\ is finite and positive, the final multiclass
probabilities are renormalized by row,

\$\$q\_{ik} = \frac{r\_{ik}}{\sum\_{\ell = 1}^K r\_{i\ell}}.\$\$

If \\S_i\\ is zero or non-finite, the prediction for that row is
replaced by the uniform distribution \\q\_{ik} = 1 / K\\. This fallback
keeps the output on the probability simplex. The renormalization changes
the individual \\r\_{ik}\\ values unless \\S_i = 1\\, so final columns
should not be interpreted as the raw outputs of the independently
calibrated binary problems. The renormalized probabilities are
simplex-valued, but the one-vs-rest reduction does not by itself
guarantee joint multiclass calibration.

## References

Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
accurate multiclass probability estimates. Proceedings of the Eighth ACM
SIGKDD International Conference on Knowledge Discovery and Data Mining.
<doi:10.1145/775047.775151>.

## Examples

``` r
set.seed(21)
raw <- matrix(stats::runif(150 * 3), ncol = 3)
raw <- raw / rowSums(raw)
labels <- max.col(raw)

fit <- cal_ovr(raw, labels, method = "isotonic")
calibrated <- predict(fit, raw)
head(calibrated)
#>              1         2         3
#> [1,] 1.0000000 0.0000000 0.0000000
#> [2,] 0.7368421 0.2631579 0.0000000
#> [3,] 0.8166667 0.1833333 0.0000000
#> [4,] 0.0000000 1.0000000 0.0000000
#> [5,] 0.7368421 0.2631579 0.0000000
#> [6,] 0.8888889 0.0000000 0.1111111
```
