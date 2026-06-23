# Platt scaling

`cal_platt()` fits a logistic regression that maps an uncalibrated score
to a calibrated probability. The binary targets are adjusted with
Platt's target correction before fitting, which shrinks labels away from
exact `0` and `1`.

## Usage

``` r
cal_platt(x, y)
```

## Arguments

- x:

  Numeric vector of uncalibrated scores or raw probabilities.

- y:

  Binary outcome vector coded as `0` and `1`.

## Value

A `cal_platt` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new scores to
obtain calibrated probabilities.

## Details

Let \\(x_i, y_i), i = 1, \ldots, n\\ be the calibration sample, where
\\x_i\\ is the supplied score and \\y_i \in \\0, 1\\\\ is the observed
label. Write \\n\_+ = \sum_i y_i\\ and \\n\_- = n - n\_+\\. Platt's
correction replaces the binary labels by fractional targets. Positive
labels use

\$\$t\_+ = \frac{n\_+ + 1}{n\_+ + 2},\$\$

and negative labels use

\$\$t\_- = \frac{1}{n\_- + 2}.\$\$

Thus \\t_i = t\_+\\ when \\y_i = 1\\ and \\t_i = t\_-\\ when \\y_i =
0\\. The fitted logistic map is

\$\$q_i(\alpha, \beta) = \operatorname{logit}^{-1}(\alpha + \beta
x_i),\$\$

and \\(\alpha, \beta)\\ are estimated by minimizing the binomial
cross-entropy with the corrected fractional targets,

\$\$\ell(\alpha, \beta) = -\sum\_{i = 1}^n \\t_i \log q_i(\alpha,
\beta) + (1 - t_i) \log\[1 - q_i(\alpha, \beta)\]\\.\$\$

The implementation fits this model with
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) using the formula
`y_adj ~ x`. The label vector must contain at least one `0` and one `1`.
The returned object stores `coefficients`, where `(Intercept)` is
\\\hat\alpha\\ and `x` is \\\hat\beta\\, as well as the full `glm`
object in `fit` and the corrected targets `target_pos` and `target_neg`.
Prediction applies \\\operatorname{logit}^{-1}(\hat\alpha + \hat\beta
x\_{new})\\ to new scores. The argument `x` may be a score on any
real-valued scale or a raw probability, but the fitted map is always a
logistic function of the supplied values. The slope is unconstrained;
the fitted map is increasing in `x` only when \\\hat\beta \ge 0\\.

## References

Platt, J. (1999). Probabilistic outputs for support vector machines and
comparisons to regularized likelihood methods. In Advances in Large
Margin Classifiers.

## Examples

``` r
set.seed(1)
calibration <- data.frame(score = rnorm(120)) |>
  dplyr::mutate(
    truth = inv_logit(score),
    y = rbinom(dplyr::n(), 1, truth)
  )

fit <- cal_platt(calibration$score, calibration$y)

calibration |>
  dplyr::mutate(calibrated = predict(fit, score)) |>
  dplyr::summarise(ece = ece(calibrated, y, bins = 10))
#>          ece
#> 1 0.07781091
```
