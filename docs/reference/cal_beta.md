# Beta calibration

`cal_beta()` fits the beta calibration model
`inv_logit(a * log(p) - b * log(1 - p) + c)`. Probabilities are clipped
to to have lower bound `eps` and upper bound `1 - eps` before taking
logarithms.

## Usage

``` r
cal_beta(p, y, eps = 1e-15)
```

## Arguments

- p:

  Numeric vector of uncalibrated probabilities in `[0, 1]`.

- y:

  Binary outcome vector coded as `0` and `1`.

- eps:

  Clipping constant satisfying `0 < eps < 0.5`. Probabilities must first
  be valid values in `[0, 1]`; values below `eps` and above `1 - eps`
  are clipped before taking logarithms.

## Value

A `cal_beta` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new
probabilities to obtain calibrated probabilities.

## Details

Beta calibration treats the uncalibrated event probability \\p_i\\
through two log-transformed features. Before the transformation,
probabilities are clipped by

\$\$p_i^\* = C\_\epsilon(p_i) = \min\\\max(p_i, \epsilon), 1 -
\epsilon\\.\$\$

The calibrated probability is

\$\$q_i = \operatorname{logit}^{-1} \\a \log(p_i^\*) - b \log(1 -
p_i^\*) + c\\.\$\$

The implementation fits an ordinary unpenalized binomial
[`glm()`](https://rdrr.io/r/stats/glm.html) with the original binary
labels, without Platt target correction. Its linear predictor is

\$\$\eta_i = \gamma_0 + \gamma_1 \log(p_i^\*) + \gamma_2 \log(1 -
p_i^\*).\$\$

Equivalently, the fitted coefficients minimize the binomial
cross-entropy

\$\$-\sum\_{i = 1}^n \\y_i \log q_i + (1 - y_i) \log(1 - q_i)\\.\$\$

The beta-calibration parameters are the following reparameterization of
the fitted [`glm()`](https://rdrr.io/r/stats/glm.html) coefficients:

\$\$\hat a = \hat\gamma_1, \quad \hat b = -\hat\gamma_2, \quad \hat c =
\hat\gamma_0.\$\$

Thus prediction first computes \\p\_{new}^\* = C\_\epsilon(p\_{new})\\
and then evaluates

\$\$\hat q(p\_{new}) = \operatorname{logit}^{-1}\\ \hat a
\log(p\_{new}^\*) - \hat b \log(1 - p\_{new}^\*) + \hat c\\.\$\$

The object element `coefficients` contains \\(\hat\gamma_0,
\hat\gamma_1, \hat\gamma_2)\\ from
[`glm()`](https://rdrr.io/r/stats/glm.html), while `a`, `b`, and `c`
contain the reparameterized beta-calibration coefficients. Since
\\d\eta_i / dp_i = a / p_i + b / (1 - p_i)\\, monotone increase on
`(0, 1)` is guaranteed when \\a \ge 0\\ and \\b \ge 0\\. The
implementation does not impose these constraints.

## References

Kull, M., Silva Filho, T. M., & Flach, P. (2017). Beta calibration: A
well-founded and easily implemented improvement on logistic calibration
for binary classifiers. Electronic Journal of Statistics, 11(2),
5052-5080. <doi:10.1214/17-EJS1338SI>.

## Examples

``` r
set.seed(3)
calibration <- data.frame(raw_p = stats::rbeta(120, 2, 2)) |>
  dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))

fit <- cal_beta(calibration$raw_p, calibration$y)

calibration |>
  dplyr::mutate(calibrated = predict(fit, raw_p)) |>
  dplyr::summarise(
    raw_ece = ece(raw_p, y, bins = 10),
    calibrated_ece = ece(calibrated, y, bins = 10)
  )
#>     raw_ece calibrated_ece
#> 1 0.1203746     0.08519107
```
