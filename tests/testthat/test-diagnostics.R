test_that("reliability diagram returns a ggplot object", {
  set.seed(6)
  p <- stats::runif(80)
  y <- rbinom(80, 1, p)

  plot <- reliability_diagram(p, y, bins = 8)

  expect_s3_class(plot, "ggplot")
  expect_match(plot$labels$subtitle, "ECE =")
})

test_that("reliability diagram can use fixed point size", {
  set.seed(6)
  p <- stats::runif(80)
  y <- rbinom(80, 1, p)

  plot <- reliability_diagram(p, y, bins = 8, show_counts = FALSE, show_ece = FALSE)

  expect_s3_class(plot, "ggplot")
  expect_null(plot$labels$subtitle)
})
