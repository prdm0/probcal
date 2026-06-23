# Dirichlet calibration

`cal_dirichlet()` is the multiclass generalization of beta calibration.
It fits a linear map on the log of the predicted probabilities followed
by a softmax, which is equivalent to a multinomial logistic regression
with the log-probabilities as features. An off-diagonal and intercept
regularization (ODIR) penalty shrinks the off-diagonal weights and the
intercepts toward zero, which reduces overfitting risk when the number
of classes is large.

## Usage

``` r
cal_dirichlet(p, y, lambda = NULL, eps = 1e-12)
```

## Arguments

- p:

  Numeric matrix of uncalibrated probabilities with one row per
  observation and one column per class. Rows must sum to one within
  absolute tolerance `1e-6`.

- y:

  A factor or a vector of integer class codes in `1:K`, where `K` is the
  number of columns of `p`.

- lambda:

  Non-negative ODIR regularization strength. When `NULL` it is chosen by
  cross-validation.

- eps:

  Clipping constant satisfying `0 < eps < 0.5`. Probabilities must first
  be valid values in `[0, 1]`; values below `eps` and above `1 - eps`
  are clipped before taking logarithms.

## Value

A `cal_dirichlet` object that also inherits from `cal_multiclass`. Use
[`predict()`](https://rdrr.io/r/stats/predict.html) with new
probabilities to obtain calibrated probabilities.

## Details

The calibrated probabilities are computed row-wise as
`softmax(log(p) %*% t(W) + b)`, where `W` is a `K` by `K` weight matrix
and `b` is a length `K` intercept vector. Probabilities are clipped to
to have lower bound `eps` and upper bound `1 - eps` before taking
logarithms. When `lambda` is `NULL`, it is selected from a small
deterministic grid by cross-validated log-likelihood.

Let \\p\_{ik}\\ be the uncalibrated probability assigned to class \\k\\
for observation \\i\\. Each row of `p` must sum to one within absolute
tolerance `1e-6`. Column \\k\\ corresponds to integer class code \\k\\;
if `y` is a factor, column \\k\\ corresponds to `levels(y)[k]`. The
entries are clipped elementwise by

\$\$p\_{ik}^\* = \min\\\max(p\_{ik}, \epsilon), 1 - \epsilon\\,\$\$

and transformed to \\u\_{ik} = \log(p\_{ik}^\*)\\. The clipped feature
matrix is not renormalized; normalization occurs only after the linear
map, through the final softmax. Dirichlet calibration fits a multinomial
logistic regression on these log-probability features,

\$\$\eta\_{ik} = b_k + \sum\_{\ell = 1}^K W\_{k\ell} u\_{i\ell},\$\$

followed by

\$\$q\_{ik} = \frac{\exp(\eta\_{ik})}{\sum\_{m = 1}^K
\exp(\eta\_{im})}.\$\$

With fixed \\\lambda\\, the fitted parameters minimize

\$\$-\frac{1}{n}\sum_i \log q\_{i y_i} + \lambda\left(\sum\_{k \ne \ell}
W\_{k\ell}^2 + \sum_k b_k^2\right).\$\$

This is the off-diagonal and intercept regularization penalty. Diagonal
weights are not penalized. For fixed `lambda`, optimization uses BFGS
with analytic gradients, initial weight matrix \\W = I_K\\, initial bias
\\b = 0\\, and `maxit = 500`. True-class probabilities entering
logarithms are clipped to `[1e-15, 1 - 1e-15]`. The returned `weight` is
a \\K \times K\\ matrix whose row \\k\\ produces the logit for class
\\k\\; `bias` is a length-\\K\\ vector of intercepts. The object also
stores `lambda`, `value`, and the optimizer `convergence` code.

If `lambda = NULL`, the implementation evaluates the grid
`c(0, 1e-4, 1e-3, 1e-2, 1e-1)` with at most three deterministic
stratified folds. Class indices are assigned to folds in their existing
order. The selected value minimizes the unweighted average of the fold
mean held-out negative log-likelihoods; ties choose the first grid
value. If fewer than two observations are available in the smallest
class during selection, the fallback value is `1e-3`. With `lambda = 0`,
the multinomial softmax parameterization is not unique: adding the same
linear function of the features to every class logit leaves all
probabilities unchanged. The calibrated probabilities are the identified
output.

## References

Kull, M., Perello-Nieto, M., Kängsepp, M., Silva Filho, T., Song, H., &
Flach, P. (2019). Beyond temperature scaling: Obtaining well-calibrated
multi-class probabilities with Dirichlet calibration. Advances in Neural
Information Processing Systems 32.

## Examples

``` r
set.seed(23)
prob <- matrix(stats::runif(200 * 3), ncol = 3)
prob <- prob / rowSums(prob)
labels <- max.col(prob)
fit <- cal_dirichlet(prob, labels)
head(predict(fit, prob))
#>                  1            2             3
#> [1,]  2.749610e-39 1.678165e-17  1.000000e+00
#> [2,] 2.481962e-177 1.000000e+00  0.000000e+00
#> [3,] 7.895177e-175 1.000000e+00  0.000000e+00
#> [4,]  4.089151e-61 1.000000e+00  1.825486e-49
#> [5,]  1.000000e+00 0.000000e+00  1.662300e-56
#> [6,] 2.490598e-103 1.000000e+00 2.702539e-321
```
