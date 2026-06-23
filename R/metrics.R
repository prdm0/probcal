# Build the (confidence, correctness) pair used by the confidence form of the
# multiclass metrics. The predicted class is the column with the highest
# probability, the confidence is that probability, and correctness indicates
# whether the predicted class matches the label.
multiclass_confidence <- function(p, y) {
  label <- check_multiclass_y(y, n_classes = ncol(p), arg = "y")
  if (nrow(p) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg p} and {.arg y} must have the same number of observations.")
  }
  predicted <- max.col(p, ties.method = "first")
  confidence <- p[cbind(seq_len(nrow(p)), predicted)]
  list(confidence = confidence, correct = as.integer(predicted == label$codes))
}

# Shared dispatcher for the matrix form of ece(), mce(), and ace(). The
# classwise form averages (or maximizes) the binary metric computed on each
# one-vs-rest column. The confidence form applies the binary metric to the
# top-label confidence and correctness.
multiclass_calibration_error <- function(metric, p, y, bins, type, aggregate) {
  check_prob_matrix(p, arg = "p")

  if (identical(type, "confidence")) {
    conf <- multiclass_confidence(p, y)
    return(metric(conf$confidence, conf$correct, bins = bins))
  }

  label <- check_multiclass_y(y, n_classes = ncol(p), arg = "y")
  per_class <- vapply(
    seq_len(ncol(p)),
    function(k) metric(p[, k], as.integer(label$codes == k), bins = bins),
    numeric(1)
  )
  aggregate(per_class)
}

bin_stats <- function(p, y, bins = 10) {
  check_probability(p, arg = "p")
  y <- check_binary_y(y, arg = "y")
  check_same_length(p, y, x_arg = "p", y_arg = "y")
  bins <- check_bins(bins)

  breaks <- seq(0, 1, length.out = bins + 1L)
  bin <- find_bins(p, breaks)

  out <- data.frame(
    bin = seq_len(bins),
    lower = breaks[-length(breaks)],
    upper = breaks[-1L],
    n = integer(bins),
    confidence = rep(NA_real_, bins),
    accuracy = rep(NA_real_, bins),
    gap = rep(NA_real_, bins)
  )

  for (i in seq_len(bins)) {
    idx <- bin == i
    out$n[i] <- sum(idx)
    if (out$n[i] > 0L) {
      out$confidence[i] <- mean(p[idx])
      out$accuracy[i] <- mean(y[idx])
      out$gap[i] <- abs(out$accuracy[i] - out$confidence[i])
    }
  }

  out
}

#' Expected Calibration Error
#'
#' `ece()` returns the empirical weighted average gap between mean confidence
#' and empirical event frequency across equal-width probability bins. It is zero
#' when confidence and accuracy match in every non-empty bin of the chosen
#' partition.
#'
#' For binary problems `p` is a probability vector. For multiclass problems `p`
#' is a probability matrix with one column per class and `type` selects the
#' multiclass definition. The `"classwise"` form averages the binary ECE over
#' the one-vs-rest columns, also known as the static calibration error. The
#' `"confidence"` form applies the binary ECE to the top-label confidence and
#' whether the predicted class is correct, which is the definition used by Guo
#' et al. (2017).
#'
#' @details
#' For binary calibration, the interval `[0, 1]` is split into \eqn{B}
#' equal-width bins. The package uses left-closed bins,
#' \eqn{I_b = \{i: (b - 1)/B \le p_i < b/B\}}{I_b is the set with (b - 1)/B <= p_i < b/B}
#' for \eqn{b < B}{b < B}, and
#' \eqn{I_B = \{i: (B - 1)/B \le p_i \le 1\}}{I_B is the set with (B - 1)/B <= p_i <= 1}
#' for the last bin. Let \eqn{n_b = |I_b|}{n_b = |I_b|} and
#' \eqn{n = \sum_b n_b}{n = sum_b n_b}. For each non-empty bin,
#'
#' \deqn{\operatorname{conf}(b) = \frac{1}{n_b}\sum_{i \in I_b} p_i,}{
#' conf(b) = (1 / n_b) sum_{i in I_b} p_i,}
#'
#' and
#'
#' \deqn{\operatorname{acc}(b) = \frac{1}{n_b}\sum_{i \in I_b} y_i.}{
#' acc(b) = (1 / n_b) sum_{i in I_b} y_i.}
#'
#' The returned empirical ECE is
#'
#' \deqn{\operatorname{ECE} =
#'   \sum_{b: n_b > 0} \frac{n_b}{n}
#'   |\operatorname{acc}(b) - \operatorname{conf}(b)|.}{
#' ECE = sum_{b: n_b > 0} (n_b / n) |acc(b) - conf(b)|.}
#'
#' Empty bins have zero weight. The estimate depends on `bins`; changing the
#' number of bins changes the empirical partition and can change the value. A
#' value of zero means equality of sample bin means for this partition, not full
#' population calibration.
#'
#' For a probability matrix, `type = "classwise"` computes the binary ECE for
#' each one-vs-rest column \eqn{p_{\cdot k}} against
#' \eqn{\mathbf{1}\{y_i = k\}}{indicator(y_i = k)} and returns their
#' arithmetic mean,
#'
#' \deqn{\operatorname{ECE}_{\mathrm{cw}} =
#'   \frac{1}{K}\sum_{k = 1}^K
#'   \operatorname{ECE}(p_{\cdot k}, \mathbf{1}\{y_i = k\}).}{
#' ECE_cw = K^{-1} sum_k ECE(p_.k, indicator(y_i = k)).}
#'
#' `type = "confidence"` uses the top-label rule
#' \eqn{\hat y_i = \min\{k: p_{ik} = \max_\ell p_{i\ell}\}}{hat y_i is the first maximal class},
#' the confidence \eqn{r_i = p_{i\hat y_i}}{r_i = p_i,hat_y_i}, and the
#' correctness indicator
#' \eqn{c_i = \mathbf{1}\{\hat y_i = y_i\}}{c_i = indicator(hat_y_i = y_i)}, then
#' applies the binary definition to \eqn{(r_i, c_i)}:
#' \eqn{\operatorname{ECE}_{\mathrm{conf}} = \operatorname{ECE}(r, c)}{ECE_conf = ECE(r, c)}.
#' For matrix inputs, column \eqn{k} corresponds to integer class code \eqn{k};
#' if `y` is a factor, column \eqn{k} corresponds to `levels(y)[k]`.
#'
#' Here "calibrated" refers to the output of a fitted calibration map. It does
#' not imply population calibration. Binary population calibration can be stated
#' as \eqn{E(Y \mid Q) = Q}{E(Y | Q) = Q} for the predicted probability random
#' variable \eqn{Q}. For top-label confidence \eqn{R}, the analogous condition
#' is \eqn{E[\mathbf{1}\{\hat Y = Y\} \mid R] = R}{E[indicator(Yhat = Y) | R] = R}.
#'
#' @param p Predicted probabilities. A numeric vector in `[0, 1]` for binary
#' problems, or a numeric matrix with one column per class for multiclass
#' problems. Matrix inputs must have finite entries in `[0, 1]`, at least two
#' columns, and rows summing to one within absolute tolerance `1e-6`.
#' @param y Outcome labels. A vector coded as `0` and `1` for binary problems,
#' or a factor or vector of integer class codes in `1:K` for multiclass
#' problems.
#' @param bins Number of equal-width bins on `[0, 1]`. Must be a single
#' positive integer.
#' @param type Multiclass aggregation, either `"classwise"` or `"confidence"`.
#' Ignored for binary inputs.
#'
#' @return A single numeric value.
#' @references
#' Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration of
#' modern neural networks. Proceedings of the 34th International Conference on
#' Machine Learning.
#' @export
#'
#' @examples
#' predictions <- data.frame(
#'   p = c(0.10, 0.20, 0.80, 0.90),
#'   y = c(0, 0, 1, 1)
#' )
#'
#' predictions |>
#'   dplyr::summarise(ece = ece(p, y, bins = 2))
#'
#' # Multiclass classwise ECE from a probability matrix.
#' set.seed(30)
#' prob <- matrix(stats::runif(150 * 3), ncol = 3)
#' prob <- prob / rowSums(prob)
#' labels <- max.col(prob)
#' ece(prob, labels, bins = 10, type = "classwise")
ece <- function(p, y, bins = 10, type = c("classwise", "confidence")) {
  if (is.matrix(p)) {
    return(multiclass_calibration_error(ece, p, y, bins, match.arg(type), mean))
  }

  stats <- bin_stats(p, y, bins)
  n_total <- sum(stats$n)
  sum((stats$n / n_total) * stats$gap, na.rm = TRUE)
}

#' Maximum Calibration Error
#'
#' `mce()` returns the largest empirical absolute gap between mean confidence
#' and empirical event frequency among non-empty equal-width bins. For
#' multiclass inputs the `"classwise"` form returns the largest binary MCE
#' across the one-vs-rest columns and the `"confidence"` form uses the
#' top-label confidence.
#'
#' @details
#' Using the same bin notation and endpoint convention as [ece()], the binary
#' empirical maximum calibration error is
#'
#' \deqn{\operatorname{MCE} =
#'   \max_{b: n_b > 0}
#'   |\operatorname{acc}(b) - \operatorname{conf}(b)|.}{
#' MCE = max_{b: n_b > 0} |acc(b) - conf(b)|.}
#'
#' Empty bins are ignored. For a multiclass probability matrix,
#' `type = "classwise"` returns the maximum of the one-vs-rest binary MCE values
#' across classes,
#'
#' \deqn{\operatorname{MCE}_{\mathrm{cw}} =
#'   \max_{1 \le k \le K}
#'   \operatorname{MCE}(p_{\cdot k}, \mathbf{1}\{y_i = k\}).}{
#' MCE_cw = max_k MCE(p_.k, indicator(y_i = k)).}
#'
#' `type = "confidence"` returns \eqn{\operatorname{MCE}(r, c)}{MCE(r, c)}
#' using the top-label confidence and correctness variables defined in [ece()].
#'
#' @inheritParams ece
#'
#' @return A single numeric value.
#' @references
#' Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration of
#' modern neural networks. Proceedings of the 34th International Conference on
#' Machine Learning.
#' @export
#'
#' @examples
#' predictions <- data.frame(
#'   p = c(0.10, 0.20, 0.80, 0.90),
#'   y = c(0, 0, 1, 1)
#' )
#'
#' predictions |>
#'   dplyr::summarise(mce = mce(p, y, bins = 2))
mce <- function(p, y, bins = 10, type = c("classwise", "confidence")) {
  if (is.matrix(p)) {
    return(multiclass_calibration_error(mce, p, y, bins, match.arg(type), max))
  }

  stats <- bin_stats(p, y, bins)
  max(stats$gap[stats$n > 0L], na.rm = TRUE)
}

#' Average Calibration Error
#'
#' `ace()` returns the empirical unweighted mean absolute calibration gap over
#' non-empty equal-width bins. Unlike `ece()`, each non-empty bin contributes
#' equally. For multiclass inputs the `"classwise"` form averages the binary ACE
#' over the one-vs-rest columns and the `"confidence"` form uses the top-label
#' confidence.
#'
#' @details
#' Using the same bin notation and endpoint convention as [ece()], let \eqn{M}
#' be the number of non-empty bins. The binary empirical average calibration
#' error is
#'
#' \deqn{\operatorname{ACE} =
#'   \frac{1}{M}\sum_{b: n_b > 0}
#'   |\operatorname{acc}(b) - \operatorname{conf}(b)|.}{
#' ACE = (1 / M) sum_{b: n_b > 0} |acc(b) - conf(b)|.}
#'
#' Unlike ECE, ACE does not weight bins by their sample sizes. Sparse bins and
#' dense bins therefore contribute equally once they are non-empty. This
#' implementation uses equal-width bins on `[0, 1]`; it does not construct
#' adaptive or equal-frequency bins. For a multiclass probability matrix,
#' `type = "classwise"` returns the arithmetic mean of the one-vs-rest binary
#' ACE values,
#'
#' \deqn{\operatorname{ACE}_{\mathrm{cw}} =
#'   \frac{1}{K}\sum_{k = 1}^K
#'   \operatorname{ACE}(p_{\cdot k}, \mathbf{1}\{y_i = k\}).}{
#' ACE_cw = K^{-1} sum_k ACE(p_.k, indicator(y_i = k)).}
#'
#' `type = "confidence"` returns \eqn{\operatorname{ACE}(r, c)}{ACE(r, c)}
#' using top-label confidence and correctness.
#'
#' @inheritParams ece
#'
#' @return A single numeric value.
#' @references
#' Niculescu-Mizil, A., & Caruana, R. (2005). Predicting good probabilities
#' with supervised learning. Proceedings of the 22nd International Conference
#' on Machine Learning.
#' @export
#'
#' @examples
#' predictions <- data.frame(
#'   p = c(0.10, 0.20, 0.80, 0.90),
#'   y = c(0, 0, 1, 1)
#' )
#'
#' predictions |>
#'   dplyr::summarise(ace = ace(p, y, bins = 2))
ace <- function(p, y, bins = 10, type = c("classwise", "confidence")) {
  if (is.matrix(p)) {
    return(multiclass_calibration_error(ace, p, y, bins, match.arg(type), mean))
  }

  stats <- bin_stats(p, y, bins)
  mean(stats$gap[stats$n > 0L], na.rm = TRUE)
}

#' Maximum Mean Calibration Error
#'
#' `mmce()` is a binning-free empirical calibration statistic built from a
#' kernel mean embedding of the calibration error. Unlike `ece()`, it does not
#' partition the probability space into bins, so it avoids sensitivity to the
#' number and placement of bins. It still depends on the kernel and bandwidth.
#' The returned value is an empirical kernel statistic, not a population
#' calibration parameter by itself.
#'
#' For a binary input the residual compares the event indicator `y` with the
#' predicted event probability `p`. For a multiclass probability matrix the
#' confidence is the top-label probability and correctness indicates whether the
#' predicted class is right. For multiclass inputs, `mmce()` implements only
#' this top-label confidence form; there is no classwise `type` argument. The
#' statistic uses a Laplacian kernel
#' \eqn{k(a, b) = \exp(-|a - b| / \text{bandwidth})}. The computation builds an
#' observation by observation kernel matrix, so both time and memory scale as
#' \eqn{O(n^2)}{O(n^2)}.
#'
#' @details
#' Let \eqn{r_i} be the scalar probability assigned to observation \eqn{i} and
#' \eqn{c_i} the corresponding binary target. In the binary case,
#' \eqn{r_i = p_i} and \eqn{c_i = y_i}. In the multiclass case, ties are broken
#' by the first class,
#' \eqn{\hat y_i = \min\{k: p_{ik} = \max_\ell p_{i\ell}\}}{hat y_i is the first maximal class},
#' \eqn{r_i = p_{i\hat y_i}}{r_i = p_i,hat_y_i}, and
#' \eqn{c_i = \mathbf{1}\{\hat y_i = y_i\}}{c_i = indicator(hat_y_i = y_i)}.
#' The residual used by the statistic is
#'
#' \deqn{e_i = c_i - r_i.}{e_i = c_i - r_i.}
#'
#' With the Laplacian kernel
#'
#' \deqn{k(r_i, r_j) = \exp\left(-\frac{|r_i - r_j|}{h}\right),}{
#' k(r_i, r_j) = exp(-|r_i - r_j| / h),}
#'
#' where \eqn{h} is `bandwidth`, the returned value is the V-statistic plug-in
#' estimate with diagonal terms,
#'
#' \deqn{\operatorname{MMCE} =
#'   \left\{\frac{1}{n^2}\sum_{i = 1}^n\sum_{j = 1}^n
#'   e_i e_j k(r_i, r_j)\right\}^{1/2}.}{
#' MMCE = sqrt(n^{-2} sum_i sum_j e_i e_j k(r_i, r_j)).}
#'
#' The square-root argument is truncated at zero after numerical computation to
#' avoid negative values caused only by floating-point error, so the returned
#' value is nonnegative.
#'
#' @param p Predicted probabilities. A numeric vector in `[0, 1]` for binary
#' problems, or a numeric matrix with one column per class for multiclass
#' problems. Matrix inputs must have finite entries in `[0, 1]`, at least two
#' columns, and rows summing to one within absolute tolerance `1e-6`.
#' @param y Outcome labels. A vector coded as `0` and `1` for binary problems,
#' or a factor or vector of integer class codes in `1:K` for multiclass
#' problems.
#' @param bandwidth Positive finite scalar bandwidth of the Laplacian kernel.
#'
#' @return A single numeric value.
#' @references
#' Kumar, A., Sarawagi, S., & Jain, U. (2018). Trainable calibration measures
#' for neural networks from kernel mean embeddings. Proceedings of the 35th
#' International Conference on Machine Learning.
#' @export
#'
#' @examples
#' set.seed(31)
#' p <- stats::runif(200)
#' y <- rbinom(200, 1, p)
#' mmce(p, y)
mmce <- function(p, y, bandwidth = 0.2) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    cli::cli_abort("Argument {.arg bandwidth} must be a single positive number.")
  }

  if (is.matrix(p)) {
    check_prob_matrix(p, arg = "p")
    confidence_pair <- multiclass_confidence(p, y)
    confidence <- confidence_pair$confidence
    correct <- confidence_pair$correct
  } else {
    check_probability(p, arg = "p")
    y <- check_binary_y(y, arg = "y")
    check_same_length(p, y, x_arg = "p", y_arg = "y")
    confidence <- p
    correct <- y
  }

  residual <- correct - confidence
  kernel <- exp(-abs(outer(confidence, confidence, "-")) / bandwidth)
  value <- as.numeric(crossprod(residual, kernel %*% residual)) / length(residual)^2
  sqrt(max(value, 0))
}
