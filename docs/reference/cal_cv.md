# Cross-validated calibration

`cal_cv()` fits a calibrator with out-of-fold predictions. The function
expects scores, probabilities, or logits that were already produced by a
model. It does not train the underlying classifier.

## Usage

``` r
cal_cv(
  x,
  y,
  method = c("platt", "temperature", "beta", "isotonic", "histogram", "vector",
    "dirichlet", "ovr"),
  folds = 5,
  seed = NULL,
  ...
)
```

## Arguments

- x:

  Numeric vector of uncalibrated values for binary calibration, or a
  numeric matrix with one column per class for multiclass calibration.
  Use logits for `method = "temperature"` and `"vector"`, probabilities
  for `"beta"`, `"isotonic"`, `"histogram"`, and `"dirichlet"`, and
  scores or probabilities for `"platt"`.

- y:

  Binary outcome vector coded as `0` and `1`, or a factor or vector of
  integer class codes in `1:K` for multiclass calibration.

- method:

  Calibration method.

- folds:

  Number of stratified folds. Must be a single integer at least `2` and
  no larger than the smallest class count.

- seed:

  Optional integer seed used only for fold assignment.

- ...:

  Additional arguments passed to the selected calibrator, such as `bins`
  for histogram binning or `base_method` for one-vs-rest calibration.

## Value

A `cal_cv` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) to apply the final
calibrator to new values. The object stores `fold_id`,
`oof_predictions`, `fold_calibrators`, and `final_calibrator`. For
binary calibration, `oof_predictions` is a numeric vector. For
multiclass calibration, it is a numeric matrix with one row per
observation and one column per class, with column names given by the
class levels.

## Details

Folds are stratified by the outcome. The returned object stores the
out-of-fold calibrated probabilities and a final calibrator fitted on
all observations for future prediction. Binary and multiclass problems
are handled through the type of `x`. A numeric vector triggers binary
calibration. A numeric matrix with one column per class triggers
multiclass calibration, the out-of-fold predictions become a matrix, and
the available methods are `"temperature"`, `"vector"`, `"dirichlet"`,
and `"ovr"`. For `method = "ovr"`, pass the binary method through
`base_method`.

Cross-validated calibration estimates how the calibration map behaves on
observations not used to fit that map. Let \\F_i \in \\1, \ldots, V\\\\
denote the fold assigned to observation \\i\\. For each fold \\v\\, a
calibrator \\\hat f^{(-v)}\\ is fitted using observations with \\F_i \ne
v\\. The out-of-fold calibrated prediction for an observation in fold
\\v\\ is then

\$\$\hat q_i^{\mathrm{oof}} = \hat f^{(-v)}(x_i), \quad F_i = v.\$\$

These out-of-fold predictions are stored in `oof_predictions` and are
useful for estimating calibration metrics without evaluating a
calibrator on the same observations used to fit it. In binary
calibration, \\\hat q_i^{\mathrm{oof}}\\ is a scalar event probability.
In multiclass calibration, it is the row vector \\(\hat
q\_{i1}^{\mathrm{oof}}, \ldots, \hat q\_{iK}^{\mathrm{oof}})\\ on the
probability simplex. After the out-of-fold predictions are computed, a
final calibrator \\\hat f\\ is fitted on all observations. The S3
[`predict()`](https://rdrr.io/r/stats/predict.html) method for a
`cal_cv` object uses this final calibrator for future data.

The folds are stratified by the observed labels. Setting `seed` affects
only the fold assignment and restores the previous random-number state
after the assignment is made. The function assumes that `x` already
contains model outputs from another classifier; it does not refit that
classifier inside each fold. Thus the predictions are out of fold for
the calibration map only, unless `x` itself was produced out of fold by
the underlying classifier.

`folds` must be at least `2` and no larger than the smallest class
count. Within each class, observations are randomly permuted and
assigned fold labels \\1, \ldots, V, 1, \ldots\\ in sequence. For
multiclass inputs, column \\k\\ corresponds to integer class code \\k\\;
if `y` is a factor, column \\k\\ corresponds to `levels(y)[k]`. For
`method = "ovr"`, `base_method` is read from `...`; if it is not
supplied, the default base method is `"platt"`.

## Examples

``` r
set.seed(7)
predictions <- data.frame(raw_p = stats::runif(120)) |>
  dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))

fit <- cal_cv(
  predictions$raw_p,
  predictions$y,
  method = "histogram",
  folds = 3,
  bins = 5,
  seed = 1
)

predictions |>
  dplyr::mutate(calibrated = fit$oof_predictions) |>
  dplyr::summarise(ece = ece(calibrated, y, bins = 5))
#>          ece
#> 1 0.06772856
```
