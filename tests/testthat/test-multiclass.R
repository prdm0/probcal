test_that("multiclass temperature scaling returns a probability matrix", {
  data <- make_multiclass(n = 200, k = 3, seed = 1)
  fit <- cal_temperature(data$logits, data$labels)
  pred <- predict(fit, data$logits)

  expect_s3_class(fit, "cal_temperature")
  expect_s3_class(fit, "cal_multiclass")
  expect_gt(fit$temperature, 0)
  expect_equal(dim(pred), dim(data$logits))
  expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-10)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("temperature equal to one recovers the softmax of the logits", {
  data <- make_multiclass(n = 150, k = 4, seed = 2)
  fit <- cal_temperature(data$logits, data$labels)
  fit$temperature <- 1
  pred <- predict(fit, data$logits)

  expect_equal(unname(pred), unname(data$prob), tolerance = 1e-10)
})

test_that("temperature scaling preserves the predicted class", {
  data <- make_multiclass(n = 200, k = 3, seed = 3)
  fit <- cal_temperature(data$logits, data$labels)
  pred <- predict(fit, data$logits)

  expect_equal(max.col(pred), max.col(data$logits))
})

test_that("multiclass temperature scaling accepts factor labels and labels columns", {
  data <- make_multiclass(n = 120, k = 3, seed = 4)
  labels <- factor(data$labels, labels = c("a", "b", "c"))
  fit <- cal_temperature(data$logits, labels)
  pred <- predict(fit, data$logits)

  expect_equal(fit$levels, c("a", "b", "c"))
  expect_equal(colnames(pred), c("a", "b", "c"))
})

test_that("multiclass temperature scaling validates label and column agreement", {
  data <- make_multiclass(n = 80, k = 3, seed = 5)
  expect_error(cal_temperature(data$logits, data$labels[1:10]))
  expect_error(predict(cal_temperature(data$logits, data$labels), data$logits[, 1:2]))
})

test_that("one-vs-rest calibration returns rows that sum to one", {
  data <- make_multiclass(n = 200, k = 3, seed = 6)
  for (method in c("platt", "beta", "isotonic", "histogram")) {
    fit <- cal_ovr(data$prob, data$labels, method = method)
    pred <- predict(fit, data$prob)

    expect_s3_class(fit, "cal_ovr")
    expect_s3_class(fit, "cal_multiclass")
    expect_equal(dim(pred), dim(data$prob))
    expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-10)
    expect_true(all(pred >= 0 & pred <= 1))
  }
})

test_that("one-vs-rest Platt calibration records score input scale", {
  data <- make_multiclass(n = 200, k = 3, seed = 60)
  fit <- cal_ovr(data$logits, data$labels, method = "platt")

  expect_equal(fit$input, "scores or probabilities (matrix)")
  expect_equal(dim(predict(fit, data$logits)), dim(data$logits))
})

test_that("one-vs-rest delegation matches a manual per-class binary fit", {
  data <- make_multiclass(n = 200, k = 3, seed = 7)
  fit <- cal_ovr(data$prob, data$labels, method = "platt")

  manual <- vapply(seq_len(3), function(j) {
    binary <- cal_platt(data$prob[, j], as.integer(data$labels == j))
    predict(binary, data$prob[, j])
  }, numeric(nrow(data$prob)))
  expected <- manual / rowSums(manual)

  expect_equal(unname(predict(fit, data$prob)), unname(expected), tolerance = 1e-10)
})

test_that("one-vs-rest calibration passes extra arguments to the binary method", {
  data <- make_multiclass(n = 200, k = 3, seed = 8)
  fit <- cal_ovr(data$prob, data$labels, method = "histogram", bins = 5)
  expect_equal(fit$calibrators[[1L]]$bins, 5L)
})

test_that("one-vs-rest temperature scaling accepts a logit matrix", {
  data <- make_multiclass(n = 150, k = 4, seed = 9)
  fit <- cal_ovr(data$logits, data$labels, method = "temperature")
  pred <- predict(fit, data$logits)

  expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-10)
})

test_that("classwise ECE is zero for a perfectly calibrated one-hot matrix", {
  codes <- rep(1:3, each = 40)
  p <- one_hot(codes, 3)

  expect_equal(ece(p, codes, bins = 10, type = "classwise"), 0)
  expect_equal(mce(p, codes, bins = 10, type = "classwise"), 0)
  expect_equal(ace(p, codes, bins = 10, type = "classwise"), 0)
})

test_that("confidence ECE matches the binary ECE when there are two classes", {
  set.seed(40)
  p1 <- stats::runif(300)
  y <- rbinom(300, 1, p1)
  p <- cbind(1 - p1, p1)
  codes <- y + 1L

  confidence <- apply(p, 1, max)
  predicted <- max.col(p, ties.method = "first")
  correct <- as.integer(predicted == codes)

  expect_equal(
    ece(p, codes, bins = 10, type = "confidence"),
    ece(confidence, correct, bins = 10)
  )
})

test_that("multiclass confidence metrics validate label length", {
  data <- make_multiclass(n = 12, k = 3, seed = 400)

  expect_error(ece(data$prob, data$labels[1:3], type = "confidence"))
  expect_error(mmce(data$prob, data$labels[1:3]))
  expect_error(reliability_diagram(data$prob, data$labels[1:3], type = "confidence"))
})

test_that("mmce is non-negative and zero for calibrated predictions", {
  set.seed(41)
  p <- stats::runif(300)
  y <- rbinom(300, 1, p)

  expect_gte(mmce(p, y), 0)
  expect_lt(mmce(p, y), 0.1)

  data <- make_multiclass(n = 200, k = 3, seed = 42)
  expect_gte(mmce(data$prob, data$labels), 0)
})

test_that("mmce uses confidence distances in the kernel", {
  p <- c(0.2, 0.7)
  y <- c(0, 1)
  bandwidth <- 0.5
  residual <- y - p
  kernel <- exp(-abs(outer(p, p, "-")) / bandwidth)
  expected <- sqrt(
    as.numeric(crossprod(residual, kernel %*% residual)) / length(p)^2
  )

  expect_equal(mmce(p, y, bandwidth = bandwidth), expected)
})

test_that("mmce validates the bandwidth", {
  expect_error(mmce(stats::runif(10), rbinom(10, 1, 0.5), bandwidth = 0))
  expect_error(mmce(stats::runif(10), rbinom(10, 1, 0.5), bandwidth = -1))
})

test_that("vector scaling returns a probability matrix", {
  data <- make_multiclass(n = 200, k = 3, seed = 10)
  fit <- cal_vector_scaling(data$logits, data$labels)
  pred <- predict(fit, data$logits)

  expect_s3_class(fit, "cal_vector_scaling")
  expect_s3_class(fit, "cal_multiclass")
  expect_equal(length(fit$scale), 3L)
  expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-10)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("Dirichlet calibration returns a probability matrix", {
  data <- make_multiclass(n = 300, k = 3, seed = 11)
  fit <- cal_dirichlet(data$prob, data$labels)
  pred <- predict(fit, data$prob)

  expect_s3_class(fit, "cal_dirichlet")
  expect_s3_class(fit, "cal_multiclass")
  expect_equal(dim(fit$weight), c(3L, 3L))
  expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-10)
  expect_true(all(pred >= 0 & pred <= 1))
})

test_that("Dirichlet calibration reduces miscalibration on skewed data", {
  set.seed(12)
  n <- 600
  k <- 3
  truth <- matrix(stats::runif(n * k), ncol = k)
  truth <- truth / rowSums(truth)
  labels <- apply(truth, 1, function(row) sample.int(k, 1, prob = row))
  raw <- truth^1.5
  raw <- raw / rowSums(raw)

  fit <- cal_dirichlet(raw, labels)
  calibrated <- predict(fit, raw)

  expect_lt(
    ece(calibrated, labels, type = "classwise"),
    ece(raw, labels, type = "classwise")
  )
})

test_that("Dirichlet lambda selection is deterministic", {
  data <- make_multiclass(n = 300, k = 3, seed = 13)
  fit1 <- cal_dirichlet(data$prob, data$labels)
  fit2 <- cal_dirichlet(data$prob, data$labels)
  expect_equal(fit1$lambda, fit2$lambda)
  expect_gte(fit1$lambda, 0)
})

test_that("Dirichlet calibration validates the clipping constant", {
  data <- make_multiclass(n = 120, k = 3, seed = 130)

  expect_error(cal_dirichlet(data$prob, data$labels, eps = 0))
  expect_error(cal_dirichlet(data$prob, data$labels, eps = 0.5))
})

test_that("high lambda shrinks the off-diagonal Dirichlet weights", {
  data <- make_multiclass(n = 300, k = 3, seed = 14)
  low <- cal_dirichlet(data$prob, data$labels, lambda = 0)
  high <- cal_dirichlet(data$prob, data$labels, lambda = 5)

  off_low <- sum(low$weight[upper.tri(low$weight) | lower.tri(low$weight)]^2)
  off_high <- sum(high$weight[upper.tri(high$weight) | lower.tri(high$weight)]^2)
  expect_lt(off_high, off_low)
})

test_that("multiclass reliability diagram returns a ggplot for both layouts", {
  data <- make_multiclass(n = 200, k = 3, seed = 15)
  classwise <- reliability_diagram(data$prob, data$labels, bins = 8, type = "classwise")
  confidence <- reliability_diagram(data$prob, data$labels, bins = 8, type = "confidence")

  expect_s3_class(classwise, "ggplot")
  expect_s3_class(confidence, "ggplot")
})

test_that("multiclass cross-validation returns an out-of-fold matrix", {
  data <- make_multiclass(n = 240, k = 3, seed = 16)
  for (method in c("temperature", "vector", "dirichlet")) {
    input <- if (method == "dirichlet") data$prob else data$logits
    fit <- cal_cv(input, data$labels, method = method, folds = 3, seed = 1)
    expect_s3_class(fit, "cal_cv")
    expect_s3_class(fit, "cal_multiclass")
    expect_equal(dim(fit$oof_predictions), dim(data$prob))
    expect_false(anyNA(fit$oof_predictions))
    expect_equal(unname(rowSums(fit$oof_predictions)), rep(1, nrow(data$prob)), tolerance = 1e-8)
  }
})

test_that("multiclass cross-validation preserves factor levels for prediction", {
  data <- make_multiclass(n = 240, k = 3, seed = 160)
  labels <- factor(data$labels, labels = c("setosa", "versicolor", "virginica"))
  fit <- cal_cv(data$logits, labels, method = "temperature", folds = 3, seed = 1)
  pred <- predict(fit, data$logits)

  expect_equal(colnames(fit$oof_predictions), levels(labels))
  expect_equal(colnames(pred), levels(labels))
})

test_that("multiclass cross-validation supports one-vs-rest with a base method", {
  data <- make_multiclass(n = 240, k = 3, seed = 17)
  fit <- cal_cv(data$prob, data$labels, method = "ovr", base_method = "isotonic", folds = 3, seed = 1)
  pred <- predict(fit, data$prob)

  expect_s3_class(fit, "cal_cv")
  expect_equal(unname(rowSums(pred)), rep(1, nrow(pred)), tolerance = 1e-8)
})

test_that("multiclass cross-validation rejects vector-only methods", {
  data <- make_multiclass(n = 120, k = 3, seed = 18)
  expect_error(cal_cv(data$prob, data$labels, method = "beta", folds = 3))
})
