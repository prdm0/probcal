# Inverse logit transformation

`inv_logit()` maps finite real values to probabilities. Mathematically
the range is `(0, 1)`, although floating-point results can round to `0`
or `1` for extreme finite inputs. It is used by temperature scaling and
by the parametric calibrators fitted with logistic regression.

## Usage

``` r
inv_logit(x)
```

## Arguments

- x:

  Numeric vector on the logit scale.

## Value

A numeric vector of probabilities with the same length as `x`.

## Details

The inverse logit, also called the logistic function, is

\$\$\operatorname{logit}^{-1}(x) = \frac{1}{1 + \exp(-x)}.\$\$

It maps real-valued scores to probabilities, is monotone increasing, and
satisfies \\\operatorname{logit}^{-1}(0) = 0.5\\. The implementation
uses [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html), which
evaluates the same transformation with stable numerical handling for
large positive or negative inputs. The implementation accepts finite
numeric inputs only; infinite values are rejected even though the
mathematical limits of the logistic function are defined. The returned
vector has the same length as `x`.

## Examples

``` r
scores <- data.frame(logit_score = c(-2, -1, 0, 1, 2)) |>
  dplyr::mutate(probability = inv_logit(logit_score))

scores
#>   logit_score probability
#> 1          -2   0.1192029
#> 2          -1   0.2689414
#> 3           0   0.5000000
#> 4           1   0.7310586
#> 5           2   0.8807971
```
