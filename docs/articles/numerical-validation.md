# Numerical Validation

``` r

library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
```

## Purpose

`probcal` is implemented in R. It has no Python dependency at runtime.
During development, optional tests compare selected outputs against
external reference implementations. These tests provide evidence that
shared methods follow the same numerical definitions where the APIs
overlap.

## Optional checks

``` r

validation_targets <- data.frame(
  reference = c("Python netcal", "Python netcal", "Python netcal", "R betacal"),
  compared = c(
    "ece(), mce(), ace()",
    "multiclass confidence ECE and temperature scaling",
    "cal_histogram() with equal-width bins",
    "cal_beta() predictions"
  ),
  test_file = c(
    "test-netcal.R",
    "test-netcal-multiclass.R",
    "test-netcal.R",
    "test-betacal.R"
  )
) |>
  mutate(runtime_dependency = "no")

validation_targets
#>       reference                                          compared
#> 1 Python netcal                               ece(), mce(), ace()
#> 2 Python netcal multiclass confidence ECE and temperature scaling
#> 3 Python netcal             cal_histogram() with equal-width bins
#> 4     R betacal                            cal_beta() predictions
#>                  test_file runtime_dependency
#> 1            test-netcal.R                 no
#> 2 test-netcal-multiclass.R                 no
#> 3            test-netcal.R                 no
#> 4           test-betacal.R                 no
```

The tests skip when the optional dependency is unavailable. This is
intentional: users should be able to install and use the package without
Python.

## Why not compare every method

Some `netcal` methods expose broader behavior than the current scope of
`probcal`. For example, `netcal` includes detection calibration,
Bayesian fitting, and optimizer-specific constraints. Those features are
outside the current package scope.

The initial validation therefore focuses on functions where the
numerical contract is directly comparable: confidence calibration
metrics and equal-width histogram binning. Additional comparisons can be
added once each convention is matched explicitly.

## Running the optional tests

The ordinary test suite runs without Python. To run the Python-backed
checks, install `reticulate`, configure Python for `reticulate`, and
install `netcal` in that Python environment. Then run the package tests
in the usual way.

``` r

devtools::test()
```

The `betacal` check runs when the R package `betacal` is installed.
