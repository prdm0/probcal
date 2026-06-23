# Vector scaling

`cal_vector_scaling()` is the multiclass generalization of temperature
scaling that gives each class its own scale and bias. It rescales a
logit matrix column by column and applies the softmax. With a single
shared scale and no bias it reduces to temperature scaling, so it is
more flexible while remaining cheap to fit.

## Usage

``` r
cal_vector_scaling(logits, y)
```

## Arguments

- logits:

  Numeric matrix of uncalibrated logits with one row per observation and
  one column per class.

- y:

  A factor or a vector of integer class codes in `1:K`, where `K` is the
  number of columns of `logits`.

## Value

A `cal_vector_scaling` object that also inherits from `cal_multiclass`.
Use [`predict()`](https://rdrr.io/r/stats/predict.html) with new logits
to obtain calibrated probabilities.

## Details

The calibrated probabilities are `softmax(s * logits + b)`, where `s` is
a length `K` vector of per-class scales applied column by column and `b`
is a length `K` vector of per-class biases. Parameters are estimated by
minimizing the average multiclass negative log-likelihood.

Let \\z\_{ik}\\ be the uncalibrated logit for observation \\i\\ and
class \\k\\. Vector scaling estimates class-specific scales \\s_k\\ and
intercepts \\b_k\\, then forms calibrated logits

\$\$\eta\_{ik} = s_k z\_{ik} + b_k.\$\$

The predicted probabilities are obtained with the softmax,

\$\$q\_{ik} = \frac{\exp(\eta\_{ik})} {\sum\_{\ell = 1}^K
\exp(\eta\_{i\ell})}.\$\$

Parameters are estimated by minimizing

\$\$L(s, b) = -\frac{1}{n}\sum\_{i = 1}^n \log q\_{i y_i}.\$\$

For multiclass labels, column \\k\\ of `logits` corresponds to class
code \\k\\; if `y` is a factor, column \\k\\ corresponds to
`levels(y)[k]`. The implementation uses
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) with method
`"BFGS"`, analytic gradients, initial scales \\s_k = 1\\, initial biases
\\b_k = 0\\, and `maxit = 500`. True-class probabilities entering
logarithms are clipped to `[1e-15, 1 - 1e-15]`. The returned object
stores `scale`, `bias`, the optimized average negative log-likelihood
`value`, and the optimizer `convergence` code.

The scales are unconstrained in the fitted optimization, so a negative
scale is possible when it improves the likelihood on the calibration
data. Unlike temperature scaling, vector scaling can change the
predicted class because scales and biases vary by class. As with any
softmax model, adding the same constant to every class bias does not
change the resulting probability vector, so the fitted bias vector is
identifiable only up to a common additive constant.

## References

Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration
of modern neural networks. Proceedings of the 34th International
Conference on Machine Learning.

## Examples

``` r
set.seed(22)
logits <- matrix(rnorm(200 * 3), ncol = 3)
labels <- max.col(logits)
fit <- cal_vector_scaling(logits, labels)
head(predict(fit, logits))
#>                 1             2            3
#> [1,] 9.344745e-29  1.000000e+00 3.648798e-26
#> [2,] 1.000000e+00 2.108993e-107 1.224984e-16
#> [3,] 1.000000e+00  6.005877e-82 2.511543e-08
#> [4,] 9.081856e-01  3.736703e-25 9.181436e-02
#> [5,] 5.769597e-14  3.027013e-48 1.000000e+00
#> [6,] 1.000000e+00  2.170890e-77 2.723746e-34
```
