test_that("logit and inverse logit are inverse transformations", {
  p <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  expect_equal(inv_logit(logit(p)), p, tolerance = 1e-12)
})

test_that("logit clips boundary probabilities", {
  z <- logit(c(0, 1), eps = 1e-6)
  expect_equal(length(z), 2)
  expect_equal(is.finite(z), c(TRUE, TRUE))
})
