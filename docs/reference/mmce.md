# Maximum Mean Calibration Error

`mmce()` is a binning-free empirical calibration statistic built from a
kernel mean embedding of the calibration error. Unlike
[`ece()`](https://prdm0.github.io/probcal/reference/ece.md), it does not
partition the probability space into bins, so it avoids sensitivity to
the number and placement of bins. It still depends on the kernel and
bandwidth. The returned value is an empirical kernel statistic, not a
population calibration parameter by itself.

## Usage

``` r
mmce(p, y, bandwidth = 0.2)
```

## Arguments

- p:

  Predicted probabilities. A numeric vector in `[0, 1]` for binary
  problems, or a numeric matrix with one column per class for multiclass
  problems. Matrix inputs must have finite entries in `[0, 1]`, at least
  two columns, and rows summing to one within absolute tolerance `1e-6`.

- y:

  Outcome labels. A vector coded as `0` and `1` for binary problems, or
  a factor or vector of integer class codes in `1:K` for multiclass
  problems.

- bandwidth:

  Positive finite scalar bandwidth of the Laplacian kernel.

## Value

A single numeric value.

## Details

For a binary input the residual compares the event indicator `y` with
the predicted event probability `p`. For a multiclass probability matrix
the confidence is the top-label probability and correctness indicates
whether the predicted class is right. For multiclass inputs, `mmce()`
implements only this top-label confidence form; there is no classwise
`type` argument. The statistic uses a Laplacian kernel \\k(a, b) =
\exp(-\|a - b\| / \text{bandwidth})\\. The computation builds an
observation by observation kernel matrix, so both time and memory scale
as \\O(n^2)\\.

Let \\r_i\\ be the scalar probability assigned to observation \\i\\ and
\\c_i\\ the corresponding binary target. In the binary case, \\r_i =
p_i\\ and \\c_i = y_i\\. In the multiclass case, ties are broken by the
first class, \\\hat y_i = \min\\k: p\_{ik} = \max\_\ell p\_{i\ell}\\\\,
\\r_i = p\_{i\hat y_i}\\, and \\c_i = \mathbf{1}\\\hat y_i = y_i\\\\.
The residual used by the statistic is

\$\$e_i = c_i - r_i.\$\$

With the Laplacian kernel

\$\$k(r_i, r_j) = \exp\left(-\frac{\|r_i - r_j\|}{h}\right),\$\$

where \\h\\ is `bandwidth`, the returned value is the V-statistic
plug-in estimate with diagonal terms,

\$\$\operatorname{MMCE} = \left\\\frac{1}{n^2}\sum\_{i = 1}^n\sum\_{j =
1}^n e_i e_j k(r_i, r_j)\right\\^{1/2}.\$\$

The square-root argument is truncated at zero after numerical
computation to avoid negative values caused only by floating-point
error, so the returned value is nonnegative.

## References

Kumar, A., Sarawagi, S., & Jain, U. (2018). Trainable calibration
measures for neural networks from kernel mean embeddings. Proceedings of
the 35th International Conference on Machine Learning.

## Examples

``` r
set.seed(31)
p <- stats::runif(200)
y <- rbinom(200, 1, p)
mmce(p, y)
#> [1] 0.03170179
```
