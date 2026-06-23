test_that("cross-validated calibration returns out-of-fold predictions", {
  set.seed(7)
  p <- stats::runif(120)
  y <- rbinom(120, 1, p)

  fit <- cal_cv(p, y, method = "histogram", folds = 3, bins = 5, seed = 10)

  expect_s3_class(fit, "cal_cv")
  expect_equal(length(fit$oof_predictions), length(p))
  expect_equal(anyNA(fit$oof_predictions), FALSE)
  expect_equal(length(predict(fit, p)), length(p))
})
