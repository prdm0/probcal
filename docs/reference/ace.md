# Average Calibration Error

`ace()` returns the empirical unweighted mean absolute calibration gap
over non-empty equal-width bins. Unlike
[`ece()`](https://prdm0.github.io/probcal/reference/ece.md), each
non-empty bin contributes equally. For multiclass inputs the
`"classwise"` form averages the binary ACE over the one-vs-rest columns
and the `"confidence"` form uses the top-label confidence.

## Usage

``` r
ace(p, y, bins = 10, type = c("classwise", "confidence"))
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
[`ece()`](https://prdm0.github.io/probcal/reference/ece.md), let \\M\\
be the number of non-empty bins. The binary empirical average
calibration error is

\$\$\operatorname{ACE} = \frac{1}{M}\sum\_{b: n_b \> 0}
\|\operatorname{acc}(b) - \operatorname{conf}(b)\|.\$\$

Unlike ECE, ACE does not weight bins by their sample sizes. Sparse bins
and dense bins therefore contribute equally once they are non-empty.
This implementation uses equal-width bins on `[0, 1]`; it does not
construct adaptive or equal-frequency bins. For a multiclass probability
matrix, `type = "classwise"` returns the arithmetic mean of the
one-vs-rest binary ACE values,

\$\$\operatorname{ACE}\_{\mathrm{cw}} = \frac{1}{K}\sum\_{k = 1}^K
\operatorname{ACE}(p\_{\cdot k}, \mathbf{1}\\y_i = k\\).\$\$

`type = "confidence"` returns \\\operatorname{ACE}(r, c)\\ using
top-label confidence and correctness.

## References

Niculescu-Mizil, A., & Caruana, R. (2005). Predicting good probabilities
with supervised learning. Proceedings of the 22nd International
Conference on Machine Learning.

## Examples

``` r
predictions <- data.frame(
  p = c(0.10, 0.20, 0.80, 0.90),
  y = c(0, 0, 1, 1)
)

predictions |>
  dplyr::summarise(ace = ace(p, y, bins = 2))
#>    ace
#> 1 0.15
```
