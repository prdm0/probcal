# Logit transformation

`logit()` maps probabilities from `(0, 1)` to the real line. Inputs must
lie in `[0, 1]`; values outside this probability interval are rejected.
Valid probabilities below `eps` and above `1 - eps` are clipped before
the transformation, because the mathematical logit is infinite at the
boundary.

## Usage

``` r
logit(p, eps = .Machine$double.eps)
```

## Arguments

- p:

  Numeric vector of probabilities in `[0, 1]`.

- eps:

  Positive clipping constant in `(0, 0.5)` used before applying the
  logit.

## Value

A numeric vector on the logit scale with the same length as `p`.

## Details

For a probability \\p \in (0, 1)\\, the logit is

\$\$\operatorname{logit}(p) = \log\left(\frac{p}{1 - p}\right).\$\$

The transformation is monotone increasing and maps probabilities below
\\0.5\\ to negative values, \\0.5\\ to zero, and probabilities above
\\0.5\\ to positive values. Because the expression is not finite at \\p
= 0\\ or \\p = 1\\, the implementation first computes

\$\$p^\* = \min\\\max(p, \epsilon), 1 - \epsilon\\,\$\$

where \\\epsilon\\ is `eps`, and then returns
\\\operatorname{logit}(p^\*)\\. The returned vector has the same length
as `p`.

## Examples

``` r
probabilities <- data.frame(p = c(0.05, 0.25, 0.5, 0.75, 0.95)) |>
  dplyr::mutate(
    logit_p = logit(p),
    recovered = inv_logit(logit_p)
  )

probabilities
#>      p   logit_p recovered
#> 1 0.05 -2.944439      0.05
#> 2 0.25 -1.098612      0.25
#> 3 0.50  0.000000      0.50
#> 4 0.75  1.098612      0.75
#> 5 0.95  2.944439      0.95
```
