# Temperature scaling

`cal_temperature()` estimates a single positive temperature parameter by
minimizing the negative log-likelihood. Inputs must be logits, not
probabilities. For binary probabilities,
[`logit()`](https://prdm0.github.io/probcal/reference/logit.md) gives
the corresponding logit. For strictly positive multiclass probability
rows, \\z\_{ik} = \log p\_{ik}\\ is a valid softmax logit
representation, up to row-wise additive constants. If probabilities have
zero entries, the user must choose and supply a transformed logit
matrix, such as clipped log-probabilities. `cal_temperature()` does not
accept or clip probability matrices.

## Usage

``` r
cal_temperature(logits, y)
```

## Arguments

- logits:

  For binary calibration, a numeric vector of uncalibrated logits. For
  multiclass calibration, a numeric matrix of logits with one row per
  observation and one column per class.

- y:

  Outcome labels. For binary calibration, a vector coded as `0` and `1`.
  For multiclass calibration, a factor or a vector of integer class
  codes in `1:K`, where `K` is the number of columns of `logits`.

## Value

A `cal_temperature` object. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new logits to
obtain calibrated probabilities. Multiclass objects also inherit from
`cal_multiclass`.

## Details

The function handles both binary and multiclass problems through the
type of `logits`. A numeric vector triggers binary temperature scaling
and the calibrated probability is `inv_logit(logits / T)`. A numeric
matrix with one column per class triggers multiclass temperature scaling
and the calibrated probabilities are `softmax(logits / T)`. Because
dividing every logit by the same positive scalar preserves the row
ordering and argmax, temperature scaling leaves the predicted class
unchanged apart from existing ties and only sharpens or softens the
probabilities.

In the binary case, let \\z_i\\ be an uncalibrated logit. For a positive
temperature \\T\\, the calibrated event probability is

\$\$q_i(T) = \operatorname{logit}^{-1}(z_i / T).\$\$

The fitted temperature is found by a bounded one-dimensional
optimization on \\\[10^{-3}, 10^3\]\\:

\$\$\hat T \in \arg\min\_{10^{-3} \le T \le 10^3} -\sum\_{i = 1}^n \\y_i
\log q_i(T) + (1 - y_i) \log\[1 - q_i(T)\]\\.\$\$

In the multiclass case, let \\z\_{ik}\\ be the logit for class \\k\\ and
observation \\i\\. The calibrated probabilities are

\$\$q\_{ik}(T) = \frac{\exp(z\_{ik} / T)} {\sum\_{\ell = 1}^K
\exp(z\_{i\ell} / T)},\$\$

and \\T\\ is chosen by minimizing the average multiclass negative
log-likelihood over the same interval,

\$\$L(T) = -\frac{1}{n}\sum\_{i = 1}^n \log q\_{i y_i}(T).\$\$

For multiclass labels, column \\k\\ of the logit matrix corresponds to
class code \\k\\. If `y` is a factor, the stored order of `levels(y)`
defines the column order. The numerical objective clips probabilities
that enter logarithms to `[1e-15, 1 - 1e-15]`. The optimization uses
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) with method
`"Brent"` and initial value `1` on the bounded interval above. The
returned object stores `temperature`, the optimizer `value`, and the
optimizer `convergence` code; multiclass fits also store `k` and
`levels`.

Values \\T \> 1\\ soften the probability vector, while values \\0 \< T
\< 1\\ make it more concentrated. Dividing all class logits by the same
positive constant preserves their order, so the predicted class is
unchanged apart from ties already present in the logits.

## References

Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration
of modern neural networks. Proceedings of the 34th International
Conference on Machine Learning.

## Examples

``` r
set.seed(2)
calibration <- data.frame(logits = rnorm(120)) |>
  dplyr::mutate(
    raw_p = inv_logit(logits),
    y = rbinom(dplyr::n(), 1, raw_p)
  )

fit <- cal_temperature(calibration$logits, calibration$y)

calibration |>
  dplyr::mutate(calibrated = predict(fit, logits)) |>
  dplyr::summarise(
    raw_ece = ece(raw_p, y, bins = 10),
    calibrated_ece = ece(calibrated, y, bins = 10)
  )
#>      raw_ece calibrated_ece
#> 1 0.07986217     0.08507891

# Multiclass temperature scaling with a logit matrix and integer labels.
set.seed(20)
logits <- matrix(rnorm(150 * 3), ncol = 3)
labels <- max.col(logits) # integer codes in 1:3
mc_fit <- cal_temperature(logits, labels)
head(predict(mc_fit, logits))
#>      1            2             3
#> [1,] 1 6.206673e-39  0.000000e+00
#> [2,] 0 0.000000e+00  1.000000e+00
#> [3,] 1 0.000000e+00  0.000000e+00
#> [4,] 0 0.000000e+00  1.000000e+00
#> [5,] 0 1.000000e+00 4.446591e-323
#> [6,] 1 3.281336e-45 3.287823e-303
```
