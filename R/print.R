# The header glyph is decorative. Set `options(probcal.emoji = FALSE)` to drop
# it, which is used when capturing console output for plain-text contexts such as
# the manuscript.
#' @export
print.calibrator <- function(x, ...) {
  prefix <- if (isTRUE(getOption("probcal.emoji", TRUE))) {
    paste0(emoji_balance, " ")
  } else {
    ""
  }
  cli::cli_h1("{prefix}probcal calibrator")
  items <- list(
    Method = x$method,
    Observations = x$n,
    Input = x$input
  )
  if (!is.null(x$k)) {
    items$Classes <- x$k
  }
  cli::cli_dl(items)
  invisible(x)
}

#' @export
summary.calibrator <- function(object, ...) {
  print(object)
  cli::cli_h2("Parameters")
  params <- cal_summary_params(object)
  if (length(params) > 0L) {
    cli::cli_dl(params)
  }
  invisible(object)
}

# Returns the method-specific parameters shown by `summary()`. Dispatch is a
# single switch on the calibrator subclass, which keeps the per-method details
# in one place and lets new calibrators add a branch without touching
# `summary.calibrator()`.
cal_summary_params <- function(object) {
  switch(
    class(object)[1L],
    cal_temperature = list(Temperature = round(object$temperature, 4L)),
    cal_platt = lapply(as.list(object$coefficients), round, digits = 4L),
    cal_beta = list(
      a = round(object$a, 4L),
      b = round(object$b, 4L),
      c = round(object$c, 4L),
      eps = object$eps
    ),
    cal_isotonic = list(Thresholds = length(object$x_thresholds)),
    cal_histogram = list(
      RequestedBins = object$bins,
      ActualBins = object$actual_bins,
      Strategy = object$strategy
    ),
    cal_vector_scaling = list(
      Scales = paste(round(object$scale, 3L), collapse = ", "),
      Biases = paste(round(object$bias, 3L), collapse = ", ")
    ),
    cal_dirichlet = list(Lambda = object$lambda),
    cal_ovr = list(BaseMethod = object$base_method),
    cal_cv = list(
      CalibrationMethod = object$calibration_method,
      Folds = object$folds
    ),
    list()
  )
}
