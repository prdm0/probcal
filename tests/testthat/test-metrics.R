test_that("calibration metrics use shared bin statistics", {
  p <- c(0.2, 0.4, 0.8, 0.9)
  y <- c(0, 1, 1, 1)

  expect_equal(ece(p, y, bins = 2), 0.175, tolerance = 1e-12)
  expect_equal(mce(p, y, bins = 2), 0.2, tolerance = 1e-12)
  expect_equal(ace(p, y, bins = 2), 0.175, tolerance = 1e-12)
})

test_that("perfect bin-level predictions have zero ECE", {
  p <- c(0, 0, 1, 1)
  y <- c(0, 0, 1, 1)

  expect_equal(ece(p, y, bins = 2), 0)
})
