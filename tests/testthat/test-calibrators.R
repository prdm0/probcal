test_that("temperature scaling returns probabilities", {
  set.seed(1)
  logits <- rnorm(120)
  y <- rbinom(120, 1, inv_logit(logits))

  fit <- cal_temperature(logits, y)
  pred <- predict(fit, logits)

  expect_s3_class(fit, "cal_temperature")
  expect_gt(fit$temperature, 0)
  expect_equal(pmin(pmax(pred, 0), 1), pred)
})

test_that("Platt scaling returns probabilities", {
  set.seed(2)
  score <- rnorm(120)
  y <- rbinom(120, 1, inv_logit(score))

  fit <- cal_platt(score, y)
  pred <- predict(fit, score)

  expect_s3_class(fit, "cal_platt")
  expect_equal(pmin(pmax(pred, 0), 1), pred)
})

test_that("beta calibration handles boundary probabilities", {
  set.seed(3)
  p <- c(0, stats::rbeta(118, 2, 2), 1)
  y <- sample(rep(c(0, 1), 60))

  fit <- cal_beta(p, y)
  pred <- predict(fit, p)

  expect_s3_class(fit, "cal_beta")
  expect_equal(pmin(pmax(pred, 0), 1), pred)
})

test_that("isotonic calibration is monotone on a grid", {
  set.seed(4)
  p <- stats::runif(120)
  y <- rbinom(120, 1, p)

  fit <- cal_isotonic(p, y)
  grid <- seq(0, 1, length.out = 50)
  pred <- predict(fit, grid)

  expect_s3_class(fit, "cal_isotonic")
  expect_equal(pred, sort(pred), tolerance = 1e-12)
  expect_equal(pmin(pmax(pred, 0), 1), pred)
})

test_that("histogram calibration supports both binning strategies", {
  set.seed(5)
  p <- stats::runif(120)
  y <- rbinom(120, 1, p)

  equal_width <- cal_histogram(p, y, bins = 6, strategy = "equal_width")
  equal_freq <- cal_histogram(p, y, bins = 6, strategy = "equal_freq")

  expect_s3_class(equal_width, "cal_histogram")
  expect_s3_class(equal_freq, "cal_histogram")
  expect_equal(length(predict(equal_width, p)), length(p))
  expect_equal(length(predict(equal_freq, p)), length(p))
})
