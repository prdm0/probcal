# Expected Calibration Error

`ece()` returns the empirical weighted average gap between mean
confidence and empirical event frequency across equal-width probability
bins. It is zero when confidence and accuracy match in every non-empty
bin of the chosen partition.

## Usage

``` r
ece(p, y, bins = 10, type = c("classwise", "confidence"))
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

- bins:

  Number of equal-width bins on `[0, 1]`. Must be a single positive
  integer.

- type:

  Multiclass aggregation, either `"classwise"` or `"confidence"`.
  Ignored for binary inputs.

## Value

A single numeric value.

## Details

For binary problems `p` is a probability vector. For multiclass problems
`p` is a probability matrix with one column per class and `type` selects
the multiclass definition. The `"classwise"` form averages the binary
ECE over the one-vs-rest columns, also known as the static calibration
error. The `"confidence"` form applies the binary ECE to the top-label
confidence and whether the predicted class is correct, which is the
definition used by Guo et al. (2017).

For binary calibration, the interval `[0, 1]` is split into \\B\\
equal-width bins. The package uses left-closed bins, \\I_b = \\i: (b -
1)/B \le p_i \< b/B\\\\ for \\b \< B\\, and \\I_B = \\i: (B - 1)/B \le
p_i \le 1\\\\ for the last bin. Let \\n_b = \|I_b\|\\ and \\n = \sum_b
n_b\\. For each non-empty bin,

\$\$\operatorname{conf}(b) = \frac{1}{n_b}\sum\_{i \in I_b} p_i,\$\$

and

\$\$\operatorname{acc}(b) = \frac{1}{n_b}\sum\_{i \in I_b} y_i.\$\$

The returned empirical ECE is

\$\$\operatorname{ECE} = \sum\_{b: n_b \> 0} \frac{n_b}{n}
\|\operatorname{acc}(b) - \operatorname{conf}(b)\|.\$\$

Empty bins have zero weight. The estimate depends on `bins`; changing
the number of bins changes the empirical partition and can change the
value. A value of zero means equality of sample bin means for this
partition, not full population calibration.

For a probability matrix, `type = "classwise"` computes the binary ECE
for each one-vs-rest column \\p\_{\cdot k}\\ against \\\mathbf{1}\\y_i =
k\\\\ and returns their arithmetic mean,

\$\$\operatorname{ECE}\_{\mathrm{cw}} = \frac{1}{K}\sum\_{k = 1}^K
\operatorname{ECE}(p\_{\cdot k}, \mathbf{1}\\y_i = k\\).\$\$

`type = "confidence"` uses the top-label rule \\\hat y_i = \min\\k:
p\_{ik} = \max\_\ell p\_{i\ell}\\\\, the confidence \\r_i = p\_{i\hat
y_i}\\, and the correctness indicator \\c_i = \mathbf{1}\\\hat y_i =
y_i\\\\, then applies the binary definition to \\(r_i, c_i)\\:
\\\operatorname{ECE}\_{\mathrm{conf}} = \operatorname{ECE}(r, c)\\. For
matrix inputs, column \\k\\ corresponds to integer class code \\k\\; if
`y` is a factor, column \\k\\ corresponds to `levels(y)[k]`.

Here "calibrated" refers to the output of a fitted calibration map. It
does not imply population calibration. Binary population calibration can
be stated as \\E(Y \mid Q) = Q\\ for the predicted probability random
variable \\Q\\. For top-label confidence \\R\\, the analogous condition
is \\E\[\mathbf{1}\\\hat Y = Y\\ \mid R\] = R\\.

## References

Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration
of modern neural networks. Proceedings of the 34th International
Conference on Machine Learning.

## Examples

``` r
predictions <- data.frame(
  p = c(0.10, 0.20, 0.80, 0.90),
  y = c(0, 0, 1, 1)
)

predictions |>
  dplyr::summarise(ece = ece(p, y, bins = 2))
#>    ece
#> 1 0.15

# Multiclass classwise ECE from a probability matrix.
set.seed(30)
prob <- matrix(stats::runif(150 * 3), ncol = 3)
prob <- prob / rowSums(prob)
labels <- max.col(prob)
ece(prob, labels, bins = 10, type = "classwise")
#> [1] 0.2264214
```
