test_that("beta calibration can be cross-checked with betacal", {
  testthat::skip_if_not_installed("betacal")

  ns <- asNamespace("betacal")
  if (!exists("beta_calibration", envir = ns, inherits = FALSE)) {
    testthat::skip("Installed betacal does not expose beta_calibration().")
  }

  set.seed(11)
  p <- stats::rbeta(180, 2.5, 2)
  y <- rbinom(180, 1, p)
  p_new <- seq(0.05, 0.95, length.out = 25)

  r_fit <- cal_beta(p, y)
  r_pred <- predict(r_fit, p_new)

  beta_calibration <- get("beta_calibration", envir = ns)
  b_fit <- beta_calibration(p, y, parameters = "abm")

  b_pred <- tryCatch(
    {
      if (is.function(b_fit)) {
        b_fit(p_new)
      } else {
        stats::predict(b_fit, p_new)
      }
    },
    error = function(e) {
      testthat::skip("Installed betacal prediction API was not recognized.")
    }
  )

  expect_equal(r_pred, as.numeric(b_pred), tolerance = 1e-6)
})
