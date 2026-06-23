test_that("multiclass temperature scaling matches netcal", {
  scaling <- import_netcal("netcal.scaling")

  set.seed(123)
  k <- 3L
  n <- 300L
  logits <- matrix(stats::rnorm(n * k), ncol = k)
  prob <- t(apply(logits, 1L, function(row) {
    z <- exp(row - max(row))
    z / sum(z)
  }))
  labels <- apply(prob, 1L, function(row) sample.int(k, 1L, prob = row))

  py_fit <- scaling$TemperatureScaling()
  py_fit$fit(prob, as.integer(labels - 1L))
  py_pred <- py_fit$transform(prob)

  r_fit <- cal_temperature(logits, labels)
  r_pred <- predict(r_fit, logits)

  # Softmax is invariant to a constant shift, so scaling the log-probabilities
  # used by netcal and scaling the raw logits used here give the same optimum.
  expect_equal(unname(r_pred), unname(py_pred), tolerance = 1e-3)
})

test_that("multiclass confidence ECE matches netcal", {
  metrics <- import_netcal("netcal.metrics.confidence")

  set.seed(321)
  k <- 4L
  n <- 400L
  logits <- matrix(stats::rnorm(n * k), ncol = k)
  prob <- t(apply(logits, 1L, function(row) {
    z <- exp(row - max(row))
    z / sum(z)
  }))
  labels <- apply(prob, 1L, function(row) sample.int(k, 1L, prob = row))

  py_ece <- as.numeric(metrics$ECE(10L)$measure(prob, as.integer(labels - 1L)))
  r_ece <- ece(prob, labels, bins = 10, type = "confidence")

  expect_equal(r_ece, py_ece, tolerance = 1e-8)
})
