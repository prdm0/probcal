# Internal infrastructure shared by the multiclass calibrators. These helpers
# are not exported. They keep the softmax cross-entropy machinery in one place
# so that temperature scaling, vector scaling, and Dirichlet calibration reuse
# the same numerically stable building blocks.

# Row-wise softmax. The maximum of each row is subtracted before exponentiating
# to avoid overflow. Recycling works column-major: because the subtracted vector
# has length equal to the number of rows, entry [i, j] is offset by row i.
softmax <- function(x) {
  m <- apply(x, 1L, max)
  z <- exp(x - m)
  z / rowSums(z)
}

# One-hot encoding of integer class codes in `1:k`.
one_hot <- function(codes, k) {
  out <- matrix(0, nrow = length(codes), ncol = k)
  out[cbind(seq_along(codes), codes)] <- 1
  out
}

# Mean negative log-likelihood of the true class under a probability matrix.
multiclass_log_loss <- function(prob, codes, eps = 1e-15) {
  picked <- prob[cbind(seq_along(codes), codes)]
  picked <- pmin(pmax(picked, eps), 1 - eps)
  -mean(log(picked))
}

# Rescale each row to sum to one. Rows that sum to zero or are not finite become
# uniform. Used by the one-vs-rest meta-calibrator after stacking per-class
# binary predictions.
renormalize_rows <- function(m) {
  row_sums <- rowSums(m)
  bad <- !is.finite(row_sums) | row_sums == 0
  if (any(bad)) {
    m[bad, ] <- 1 / ncol(m)
    row_sums[bad] <- 1
  }
  m / row_sums
}

# Dense multinomial logistic regression fitted by minimizing the softmax
# cross-entropy. The model maps features to class logits through a `k x p`
# weight matrix `W` and a length-`k` bias `b`, then applies the softmax.
#
# With `lambda > 0` and square features (`p == k`), an off-diagonal and
# intercept regularization (ODIR) penalty is added: the squared off-diagonal
# entries of `W` and the squared entries of `b`. This is the regularization used
# by Dirichlet calibration (Kull et al. 2019). The diagonal of `W` is never
# penalized.
#
# An analytic gradient is supplied so the optimization is well behaved even with
# many classes. The objective is convex when `lambda >= 0`.
fit_multinomial <- function(features, codes, k, lambda = 0, init = NULL) {
  features <- as.matrix(features)
  n <- nrow(features)
  p <- ncol(features)
  targets <- one_hot(codes, k)

  offdiag <- NULL
  if (p == k) {
    offdiag <- matrix(1, nrow = k, ncol = k)
    diag(offdiag) <- 0
  }

  unpack <- function(par) {
    weight <- matrix(par[seq_len(k * p)], nrow = k, ncol = p)
    bias <- par[k * p + seq_len(k)]
    list(weight = weight, bias = bias)
  }

  logits_of <- function(par) {
    par <- unpack(par)
    features %*% t(par$weight) + matrix(par$bias, nrow = n, ncol = k, byrow = TRUE)
  }

  objective <- function(par) {
    prob <- softmax(logits_of(par))
    value <- multiclass_log_loss(prob, codes)
    if (lambda > 0) {
      fit <- unpack(par)
      penalized <- if (is.null(offdiag)) fit$weight else fit$weight * offdiag
      value <- value + lambda * (sum(penalized^2) + sum(fit$bias^2))
    }
    value
  }

  gradient <- function(par) {
    fit <- unpack(par)
    prob <- softmax(logits_of(par))
    resid <- (prob - targets) / n
    grad_weight <- t(resid) %*% features
    grad_bias <- colSums(resid)
    if (lambda > 0) {
      penalized <- if (is.null(offdiag)) fit$weight else fit$weight * offdiag
      grad_weight <- grad_weight + 2 * lambda * penalized
      grad_bias <- grad_bias + 2 * lambda * fit$bias
    }
    c(as.vector(grad_weight), grad_bias)
  }

  if (is.null(init)) {
    start_weight <- if (p == k) diag(k) else matrix(0, nrow = k, ncol = p)
    init <- c(as.vector(start_weight), rep(0, k))
  }

  opt <- stats::optim(
    init,
    objective,
    gradient,
    method = "BFGS",
    control = list(maxit = 500)
  )

  fit <- unpack(opt$par)
  list(
    weight = fit$weight,
    bias = fit$bias,
    value = opt$value,
    convergence = opt$convergence
  )
}

#' One-vs-rest multiclass calibration
#'
#' `cal_ovr()` extends any binary calibrator to a multiclass problem with the
#' one-vs-rest reduction. For each class it fits a binary calibrator that
#' separates that class from the others, applies the calibrators column by
#' column, and renormalizes each row to sum to one. This is the default strategy
#' that binning methods use for multiclass calibration.
#'
#' The columns of `x` are the per-class uncalibrated values. Use scores or
#' probabilities for `method = "platt"`, probabilities for `"beta"`,
#' `"isotonic"`, and `"histogram"`, and binary one-vs-rest logits for
#' `"temperature"`. Rows of `x` are not required to sum to one. Every class
#' must appear at least once in `y`, because each one-vs-rest problem needs both
#' labels.
#'
#' @details
#' For \eqn{K} classes, column \eqn{k} of `x` corresponds to integer class code
#' \eqn{k}; if `y` is a factor, column \eqn{k} corresponds to `levels(y)[k]`.
#' One-vs-rest calibration creates \eqn{K} binary labels,
#'
#' \deqn{y_i^{(k)} = \mathbf{1}\{y_i = k\},
#'   \quad k = 1, \ldots, K.}{
#' y_i^(k) = 1 if y_i = k and 0 otherwise, for k = 1, ..., K.}
#'
#' A separate binary calibrator \eqn{f_k} is fitted to column \eqn{k} of `x` and
#' the binary labels \eqn{y_i^{(k)}}. On new data, the classwise calibrated
#' scores are
#'
#' \deqn{r_{ik} = f_k(x_{ik}).}{r_ik = f_k(x_ik).}
#'
#' Because the \eqn{K} binary calibrators are fitted independently, the row sums
#' of \eqn{r_{ik}} need not equal one. Let
#' \eqn{S_i = \sum_{\ell = 1}^K r_{i\ell}}{S_i = sum_l r_il}. If
#' \eqn{S_i} is finite and positive, the final multiclass probabilities are
#' renormalized by row,
#'
#' \deqn{q_{ik} = \frac{r_{ik}}{\sum_{\ell = 1}^K r_{i\ell}}.}{
#' q_ik = r_ik / sum_l r_il.}
#'
#' If \eqn{S_i} is zero or non-finite, the prediction for that row is replaced
#' by the uniform distribution \eqn{q_{ik} = 1 / K}{q_ik = 1 / K}. This fallback
#' keeps the output on the probability simplex. The renormalization changes the
#' individual \eqn{r_{ik}} values unless \eqn{S_i = 1}{S_i = 1}, so final
#' columns should not be interpreted as the raw outputs of the independently
#' calibrated binary problems. The renormalized probabilities are
#' simplex-valued, but the one-vs-rest reduction does not by itself guarantee
#' joint multiclass calibration.
#'
#' @param x Numeric matrix of uncalibrated values with one row per observation
#' and one column per class. For `method = "platt"`, entries may be arbitrary
#' finite scores. For `"beta"`, `"isotonic"`, and `"histogram"`, entries must
#' be probabilities in `[0, 1]`. For `"temperature"`, entries are logits.
#' @param y A factor or a vector of integer class codes in `1:K`, where `K` is
#' the number of columns of `x`.
#' @param method Binary calibrator applied to each one-vs-rest problem.
#' @param ... Additional arguments passed to the binary calibrator, such as
#' `bins` for `method = "histogram"`.
#'
#' @return A `cal_ovr` object that also inherits from `cal_multiclass`. The
#' object stores `calibrators`, `base_method`, `k`, `levels`, `input`, and the
#' original call. Use `predict()` with a new score matrix to obtain a numeric
#' matrix of calibrated probabilities whose rows sum to one.
#' @references
#' Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
#' accurate multiclass probability estimates. Proceedings of the Eighth ACM
#' SIGKDD International Conference on Knowledge Discovery and Data Mining.
#' <doi:10.1145/775047.775151>.
#' @export
#'
#' @examples
#' set.seed(21)
#' raw <- matrix(stats::runif(150 * 3), ncol = 3)
#' raw <- raw / rowSums(raw)
#' labels <- max.col(raw)
#'
#' fit <- cal_ovr(raw, labels, method = "isotonic")
#' calibrated <- predict(fit, raw)
#' head(calibrated)
cal_ovr <- function(
  x,
  y,
  method = c("platt", "beta", "isotonic", "histogram", "temperature"),
  ...
) {
  method <- match.arg(method)

  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("Argument {.arg x} must be a numeric matrix with one column per class.")
  }
  if (ncol(x) < 2L) {
    cli::cli_abort("Argument {.arg x} must have at least two columns, one per class.")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    cli::cli_abort("Argument {.arg x} must not contain missing or infinite values.")
  }

  label <- check_multiclass_y(y, n_classes = ncol(x), arg = "y")
  if (nrow(x) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg x} and {.arg y} must have the same number of observations.")
  }

  probability_based <- method %in% c("platt", "beta", "isotonic", "histogram")
  if (probability_based && method != "platt" && any(x < 0 | x > 1)) {
    cli::cli_abort("Argument {.arg x} must contain probabilities in `[0, 1]` for {.val {method}}.")
  }

  present <- sort(unique(label$codes))
  if (length(present) != label$k) {
    cli::cli_abort("Every class must appear in {.arg y} for one-vs-rest calibration.")
  }

  calibrators <- vector("list", label$k)
  for (j in seq_len(label$k)) {
    target <- as.integer(label$codes == j)
    calibrators[[j]] <- fit_calibrator(method, x[, j], target, ...)
  }

  input <- switch(
    method,
    platt = "scores or probabilities (matrix)",
    temperature = "logits (matrix)",
    "probabilities (matrix)"
  )

  object <- new_calibrator(
    "cal_ovr",
    method = paste0("one-vs-rest ", calibrators[[1L]]$method),
    n = nrow(x),
    input = input,
    base_method = method,
    calibrators = calibrators,
    k = label$k,
    levels = label$levels,
    call = match.call()
  )
  class(object) <- c("cal_ovr", "cal_multiclass", "calibrator")
  object
}
