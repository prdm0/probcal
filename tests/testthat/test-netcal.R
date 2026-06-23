test_that("confidence metrics match netcal on binary probabilities", {
  metrics <- import_netcal("netcal.metrics.confidence")

  bins <- 10L
  p <- (seq_len(200) - 0.37) / 201
  y <- as.integer(sin(seq_along(p) * 1.7) > 0)

  py_ece <- as.numeric(metrics$ECE(as.integer(bins))$measure(p, y))
  py_mce <- as.numeric(metrics$MCE(as.integer(bins))$measure(p, y))
  py_ace <- as.numeric(metrics$ACE(as.integer(bins))$measure(p, y))

  expect_equal(ece(p, y, bins = bins), py_ece, tolerance = 1e-10)
  expect_equal(mce(p, y, bins = bins), py_mce, tolerance = 1e-10)
  expect_equal(ace(p, y, bins = bins), py_ace, tolerance = 1e-10)
})

test_that("equal-width histogram binning matches netcal", {
  binning <- import_netcal("netcal.binning")

  bins <- 5L
  p_train <- (seq_len(200) - 0.5) / 200
  y_train <- as.integer((seq_len(200) %% 7) %in% c(0, 1, 2))
  p_new <- c(0.03, 0.12, 0.28, 0.41, 0.56, 0.73, 0.91)

  r_fit <- cal_histogram(
    p_train,
    y_train,
    bins = bins,
    strategy = "equal_width"
  )
  r_pred <- predict(r_fit, p_new)

  py_fit <- binning$HistogramBinning(
    as.integer(bins),
    equal_intervals = TRUE
  )
  py_fit$fit(p_train, y_train)
  py_pred <- as.numeric(py_fit$transform(p_new))

  expect_equal(r_pred, py_pred, tolerance = 1e-12)
})
