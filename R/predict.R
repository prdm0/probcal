#' @export
predict.cal_platt <- function(object, newdata, ...) {
  check_numeric_vector(newdata, arg = "newdata")
  coefficients <- object$coefficients
  eta <- unname(coefficients[["(Intercept)"]] + coefficients[["x"]] * newdata)
  inv_logit(eta)
}

#' @export
predict.cal_temperature <- function(object, newdata, ...) {
  if (inherits(object, "cal_multiclass") || is.matrix(newdata)) {
    check_logit_matrix(newdata, arg = "newdata")
    if (ncol(newdata) != object$k) {
      cli::cli_abort("Argument {.arg newdata} must have {object$k} column{?s}, one per class.")
    }
    prob <- softmax(newdata / object$temperature)
    colnames(prob) <- object$levels
    return(prob)
  }

  check_numeric_vector(newdata, arg = "newdata")
  inv_logit(newdata / object$temperature)
}

#' @export
predict.cal_beta <- function(object, newdata, ...) {
  check_probability(newdata, arg = "newdata")
  p <- clip_prob(newdata, eps = object$eps, arg = "newdata")
  eta <- object$a * log(p) - object$b * log1p(-p) + object$c
  inv_logit(eta)
}

#' @export
predict.cal_isotonic <- function(object, newdata, ...) {
  check_probability(newdata, arg = "newdata")

  if (length(object$x_thresholds) == 1L) {
    return(rep(object$y_calibrated, length(newdata)))
  }

  out <- stats::approx(
    x = object$x_thresholds,
    y = object$y_calibrated,
    xout = newdata,
    method = "linear",
    yleft = object$y_calibrated[1L],
    yright = object$y_calibrated[length(object$y_calibrated)],
    ties = "ordered"
  )$y

  pmin(pmax(out, 0), 1)
}

#' @export
predict.cal_histogram <- function(object, newdata, ...) {
  check_probability(newdata, arg = "newdata")
  bin <- find_bins(newdata, object$breaks)
  object$bin_values[bin]
}

#' @export
predict.cal_vector_scaling <- function(object, newdata, ...) {
  check_logit_matrix(newdata, arg = "newdata")
  if (ncol(newdata) != object$k) {
    cli::cli_abort("Argument {.arg newdata} must have {object$k} column{?s}, one per class.")
  }
  logits <- sweep(newdata, 2L, object$scale, "*") +
    matrix(object$bias, nrow = nrow(newdata), ncol = object$k, byrow = TRUE)
  prob <- softmax(logits)
  colnames(prob) <- object$levels
  prob
}

#' @export
predict.cal_dirichlet <- function(object, newdata, ...) {
  check_prob_matrix(newdata, arg = "newdata")
  if (ncol(newdata) != object$k) {
    cli::cli_abort("Argument {.arg newdata} must have {object$k} column{?s}, one per class.")
  }
  features <- log(pmin(pmax(newdata, object$eps), 1 - object$eps))
  logits <- features %*% t(object$weight) +
    matrix(object$bias, nrow = nrow(newdata), ncol = object$k, byrow = TRUE)
  prob <- softmax(logits)
  colnames(prob) <- object$levels
  prob
}

#' @export
predict.cal_ovr <- function(object, newdata, ...) {
  if (!is.matrix(newdata) || !is.numeric(newdata)) {
    cli::cli_abort("Argument {.arg newdata} must be a numeric matrix with one column per class.")
  }
  if (ncol(newdata) != object$k) {
    cli::cli_abort("Argument {.arg newdata} must have {object$k} column{?s}, one per class.")
  }

  columns <- lapply(seq_len(object$k), function(j) {
    stats::predict(object$calibrators[[j]], newdata = newdata[, j])
  })
  prob <- renormalize_rows(do.call(cbind, columns))
  colnames(prob) <- object$levels
  prob
}

#' @export
predict.cal_cv <- function(object, newdata, ...) {
  stats::predict(object$final_calibrator, newdata = newdata, ...)
}
