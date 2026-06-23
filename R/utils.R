#' Logit transformation
#'
#' `logit()` maps probabilities from `(0, 1)` to the real line. Inputs must lie
#' in `[0, 1]`; values outside this probability interval are rejected. Valid
#' probabilities below `eps` and above `1 - eps` are clipped before the
#' transformation, because the mathematical logit is infinite at the boundary.
#'
#' @details
#' For a probability \eqn{p \in (0, 1)}, the logit is
#'
#' \deqn{\operatorname{logit}(p) = \log\left(\frac{p}{1 - p}\right).}{
#' logit(p) = log(p / (1 - p)).}
#'
#' The transformation is monotone increasing and maps probabilities below
#' \eqn{0.5} to negative values, \eqn{0.5} to zero, and probabilities above
#' \eqn{0.5} to positive values. Because the expression is not finite at
#' \eqn{p = 0} or \eqn{p = 1}, the implementation first computes
#'
#' \deqn{p^* = \min\{\max(p, \epsilon), 1 - \epsilon\},}{
#' p^* = min(max(p, eps), 1 - eps),}
#'
#' where \eqn{\epsilon} is `eps`, and then returns
#' \eqn{\operatorname{logit}(p^*)}. The returned vector has the same length as
#' `p`.
#'
#' @param p Numeric vector of probabilities in `[0, 1]`.
#' @param eps Positive clipping constant in `(0, 0.5)` used before applying the
#' logit.
#'
#' @return A numeric vector on the logit scale with the same length as `p`.
#' @export
#'
#' @examples
#' probabilities <- data.frame(p = c(0.05, 0.25, 0.5, 0.75, 0.95)) |>
#'   dplyr::mutate(
#'     logit_p = logit(p),
#'     recovered = inv_logit(logit_p)
#'   )
#'
#' probabilities
logit <- function(p, eps = .Machine$double.eps) {
  check_probability(p, arg = "p")
  p <- clip_prob(p, eps = eps, arg = "p")
  log(p / (1 - p))
}

#' Inverse logit transformation
#'
#' `inv_logit()` maps finite real values to probabilities. Mathematically the
#' range is `(0, 1)`, although floating-point results can round to `0` or `1`
#' for extreme finite inputs. It is used by temperature scaling and by the
#' parametric calibrators fitted with logistic regression.
#'
#' @details
#' The inverse logit, also called the logistic function, is
#'
#' \deqn{\operatorname{logit}^{-1}(x) = \frac{1}{1 + \exp(-x)}.}{
#' logit^{-1}(x) = 1 / (1 + exp(-x)).}
#'
#' It maps real-valued scores to probabilities, is monotone increasing, and
#' satisfies \eqn{\operatorname{logit}^{-1}(0) = 0.5}. The implementation uses
#' `stats::plogis()`, which evaluates the same transformation with stable
#' numerical handling for large positive or negative inputs. The implementation
#' accepts finite numeric inputs only; infinite values are rejected even though
#' the mathematical limits of the logistic function are defined. The returned
#' vector has the same length as `x`.
#'
#' @param x Numeric vector on the logit scale.
#'
#' @return A numeric vector of probabilities with the same length as `x`.
#' @export
#'
#' @examples
#' scores <- data.frame(logit_score = c(-2, -1, 0, 1, 2)) |>
#'   dplyr::mutate(probability = inv_logit(logit_score))
#'
#' scores
inv_logit <- function(x) {
  check_numeric_vector(x, arg = "x")
  stats::plogis(x)
}

clip_prob <- function(p, eps = .Machine$double.eps, arg = "p") {
  check_numeric_vector(p, arg = arg)

  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0 || eps >= 0.5) {
    cli::cli_abort("Argument {.arg eps} must be a single number in `(0, 0.5)`.")
  }

  pmin(pmax(p, eps), 1 - eps)
}

check_numeric_vector <- function(x, arg = "x", finite = TRUE) {
  if (!is.numeric(x) || is.null(x) || !is.atomic(x)) {
    cli::cli_abort("Argument {.arg {arg}} must be a numeric vector.")
  }

  if (length(x) == 0L) {
    cli::cli_abort("Argument {.arg {arg}} must not be empty.")
  }

  if (anyNA(x)) {
    cli::cli_abort("Argument {.arg {arg}} must not contain missing values.")
  }

  if (isTRUE(finite) && any(!is.finite(x))) {
    cli::cli_abort("Argument {.arg {arg}} must contain only finite values.")
  }

  invisible(x)
}

check_probability <- function(p, arg = "p") {
  check_numeric_vector(p, arg = arg)

  if (any(p < 0 | p > 1)) {
    cli::cli_abort("Argument {.arg {arg}} must contain probabilities in `[0, 1]`.")
  }

  invisible(p)
}

check_binary_y <- function(y, arg = "y", require_both = FALSE) {
  if (is.logical(y)) {
    y <- as.integer(y)
  } else if (is.numeric(y) && is.atomic(y)) {
    if (anyNA(y) || any(!is.finite(y))) {
      cli::cli_abort("Argument {.arg {arg}} must not contain missing or infinite values.")
    }
    if (any(!(y %in% c(0, 1)))) {
      cli::cli_abort("Argument {.arg {arg}} must contain only binary labels `0` and `1`.")
    }
    y <- as.integer(y)
  } else {
    cli::cli_abort("Argument {.arg {arg}} must contain binary labels `0` and `1`.")
  }

  if (length(y) == 0L) {
    cli::cli_abort("Argument {.arg {arg}} must not be empty.")
  }

  if (isTRUE(require_both) && length(unique(y)) < 2L) {
    cli::cli_abort("Argument {.arg {arg}} must contain at least one `0` and one `1`.")
  }

  y
}

check_same_length <- function(x, y, x_arg = "x", y_arg = "y") {
  if (length(x) != length(y)) {
    cli::cli_abort(
      "Arguments {.arg {x_arg}} and {.arg {y_arg}} must have the same length."
    )
  }

  invisible(NULL)
}

check_bins <- function(bins, arg = "bins") {
  if (!is.numeric(bins) || length(bins) != 1L || is.na(bins) || !is.finite(bins)) {
    cli::cli_abort("Argument {.arg {arg}} must be a single positive integer.")
  }

  if (bins < 1 || bins != as.integer(bins)) {
    cli::cli_abort("Argument {.arg {arg}} must be a single positive integer.")
  }

  as.integer(bins)
}

check_folds <- function(folds, y) {
  folds <- check_bins(folds, arg = "folds")

  if (folds < 2L) {
    cli::cli_abort("Argument {.arg folds} must be at least `2`.")
  }

  class_counts <- table(y)
  if (folds > min(class_counts)) {
    cli::cli_abort(
      "Argument {.arg folds} must not exceed the smallest class count."
    )
  }

  folds
}

check_prob_matrix <- function(p, arg = "p", tolerance = 1e-6) {
  if (!is.matrix(p) || !is.numeric(p)) {
    cli::cli_abort("Argument {.arg {arg}} must be a numeric matrix with one column per class.")
  }

  if (ncol(p) < 2L) {
    cli::cli_abort("Argument {.arg {arg}} must have at least two columns, one per class.")
  }

  if (nrow(p) == 0L) {
    cli::cli_abort("Argument {.arg {arg}} must not be empty.")
  }

  if (anyNA(p) || any(!is.finite(p))) {
    cli::cli_abort("Argument {.arg {arg}} must not contain missing or infinite values.")
  }

  if (any(p < 0 | p > 1)) {
    cli::cli_abort("Argument {.arg {arg}} must contain probabilities in `[0, 1]`.")
  }

  row_sums <- rowSums(p)
  if (any(abs(row_sums - 1) > tolerance)) {
    cli::cli_abort("Each row of {.arg {arg}} must sum to `1`.")
  }

  invisible(p)
}

check_logit_matrix <- function(x, arg = "logits") {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("Argument {.arg {arg}} must be a numeric matrix with one column per class.")
  }

  if (ncol(x) < 2L) {
    cli::cli_abort("Argument {.arg {arg}} must have at least two columns, one per class.")
  }

  if (nrow(x) == 0L) {
    cli::cli_abort("Argument {.arg {arg}} must not be empty.")
  }

  if (anyNA(x) || any(!is.finite(x))) {
    cli::cli_abort("Argument {.arg {arg}} must not contain missing or infinite values.")
  }

  invisible(x)
}

check_multiclass_y <- function(y, n_classes, arg = "y") {
  if (anyNA(y)) {
    cli::cli_abort("Argument {.arg {arg}} must not contain missing values.")
  }

  if (is.factor(y)) {
    levels_y <- levels(y)
    codes <- as.integer(y)
  } else if (is.numeric(y) && is.atomic(y)) {
    if (any(!is.finite(y)) || any(y != as.integer(y))) {
      cli::cli_abort("Argument {.arg {arg}} must contain integer class codes or a factor.")
    }
    codes <- as.integer(y)
    if (min(codes) < 1L) {
      cli::cli_abort("Integer class codes in {.arg {arg}} must start at `1`.")
    }
    levels_y <- as.character(seq_len(max(codes)))
  } else {
    cli::cli_abort("Argument {.arg {arg}} must be a factor or a vector of integer class codes.")
  }

  k <- length(levels_y)
  if (k < 2L) {
    cli::cli_abort("Argument {.arg {arg}} must have at least two classes.")
  }

  if (!missing(n_classes) && !is.null(n_classes) && k != n_classes) {
    cli::cli_abort(
      "Argument {.arg {arg}} has {k} class{?es} but the input has {n_classes} column{?s}."
    )
  }

  if (max(codes) > k) {
    cli::cli_abort("Class codes in {.arg {arg}} must not exceed the number of classes.")
  }

  list(codes = codes, levels = levels_y, k = k)
}

new_calibrator <- function(subclass, method, n, input, ..., call = NULL) {
  object <- list(method = method, n = n, input = input, ..., call = call)
  class(object) <- c(subclass, "calibrator")
  object
}

binary_log_loss <- function(p, y, eps = 1e-15) {
  p <- clip_prob(p, eps = eps, arg = "p")
  -sum(y * log(p) + (1 - y) * log1p(-p))
}

find_bins <- function(p, breaks) {
  n_bins <- length(breaks) - 1L
  bin <- findInterval(p, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  bin[bin < 1L] <- 1L
  bin[bin > n_bins] <- n_bins
  as.integer(bin)
}

fill_empty_bins <- function(values, counts, fallback) {
  empty <- which(counts == 0L)
  if (length(empty) == 0L) {
    return(values)
  }

  non_empty <- which(counts > 0L)
  if (length(non_empty) == 0L) {
    values[] <- fallback
    return(values)
  }

  for (i in empty) {
    nearest <- non_empty[which.min(abs(non_empty - i))]
    values[i] <- values[nearest]
  }

  values
}
