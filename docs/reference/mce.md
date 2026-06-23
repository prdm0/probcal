# Maximum Calibration Error

`mce()` returns the largest empirical absolute gap between mean
confidence and empirical event frequency among non-empty equal-width
bins. For multiclass inputs the `"classwise"` form returns the largest
binary MCE across the one-vs-rest columns and the `"confidence"` form
uses the top-label confidence.

## Usage

``` r
mce(p, y, bins = 10, type = c("classwise", "confidence"))
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

Using the same bin notation and endpoint convention as
[`ece()`](https://prdm0.github.io/probcal/reference/ece.md), the binary
empirical maximum calibration error is

\$\$\operatorname{MCE} = \max\_{b: n_b \> 0} \|\operatorname{acc}(b) -
\operatorname{conf}(b)\|.\$\$

Empty bins are ignored. For a multiclass probability matrix,
`type = "classwise"` returns the maximum of the one-vs-rest binary MCE
values across classes,

\$\$\operatorname{MCE}\_{\mathrm{cw}} = \max\_{1 \le k \le K}
\operatorname{MCE}(p\_{\cdot k}, \mathbf{1}\\y_i = k\\).\$\$

`type = "confidence"` returns \\\operatorname{MCE}(r, c)\\ using the
top-label confidence and correctness variables defined in
[`ece()`](https://prdm0.github.io/probcal/reference/ece.md).

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
  dplyr::summarise(mce = mce(p, y, bins = 2))
#>    mce
#> 1 0.15
```
