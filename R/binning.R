#' Isotonic calibration
#'
#' `cal_isotonic()` fits a monotone calibration curve with `stats::isoreg()`.
#' New probabilities are calibrated by linear interpolation. Predictions below
#' the training range use the leftmost fitted value; predictions above the range
#' use the rightmost fitted value.
#'
#' Ties in the training probabilities are ordered with positive labels first
#' before isotonic regression and then collapsed to a single fitted value per
#' unique probability.
#'
#' @details
#' Isotonic calibration estimates a nondecreasing function \eqn{g} that maps raw
#' probabilities to calibrated event probabilities. Let \eqn{\pi}{pi} be the
#' ordering that sorts observations by increasing \eqn{p_i} and, for equal
#' \eqn{p_i}, decreasing \eqn{y_i}. Thus positive labels precede negative labels
#' within a tied probability value. The fitted values solve the projection
#' problem
#'
#' \deqn{\min_{m_1 \le \cdots \le m_n}
#'   \sum_{i = 1}^n (y_{\pi(i)} - m_i)^2.}{
#' minimize sum_i (y_pi(i) - m_i)^2 subject to m_1 <= ... <= m_n.}
#'
#' The implementation uses `stats::isoreg()` for the constrained least-squares
#' problem and clips the fitted values to `[0, 1]`. The label vector must
#' contain at least one `0` and one `1`.
#'
#' Prediction uses linear interpolation between the unique training
#' probabilities and their fitted values. If a new probability is below the
#' smallest training value, prediction returns the leftmost fitted value. If it
#' is above the largest training value, prediction returns the rightmost fitted
#' value. Training ties are collapsed to one fitted value per unique probability
#' after the isotonic fit by averaging the fitted values within each tied group.
#' If the training data contain a single unique probability, prediction is the
#' resulting constant fitted value. The fitted object stores the unique
#' probabilities in `x_thresholds`, the collapsed fitted values in
#' `y_calibrated`, the `stats::isoreg()` object in `fit`, and the original call.
#' Prediction uses `stats::approx(method = "linear")` with constant
#' extrapolation at the two endpoints, so the package prediction rule is the
#' interpolated monotone curve rather than the unmodified PAVA step function.
#'
#' @param p Numeric vector of uncalibrated probabilities in `[0, 1]`.
#' @param y Binary outcome vector coded as `0` and `1`.
#'
#' @return A `cal_isotonic` object. Use `predict()` with new probabilities to
#' obtain calibrated probabilities.
#' @references
#' Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
#' accurate multiclass probability estimates. Proceedings of the Eighth ACM
#' SIGKDD International Conference on Knowledge Discovery and Data Mining.
#' <doi:10.1145/775047.775151>.
#' @export
#'
#' @examples
#' set.seed(4)
#' calibration <- data.frame(raw_p = sort(stats::runif(120))) |>
#'   dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))
#'
#' fit <- cal_isotonic(calibration$raw_p, calibration$y)
#'
#' calibration |>
#'   dplyr::mutate(calibrated = predict(fit, raw_p)) |>
#'   dplyr::summarise(
#'     raw_ece = ece(raw_p, y, bins = 10),
#'     calibrated_ece = ece(calibrated, y, bins = 10)
#'   )
cal_isotonic <- function(p, y) {
  check_probability(p, arg = "p")
  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(p, y, x_arg = "p", y_arg = "y")

  order_id <- order(p, -y)
  p_ordered <- p[order_id]
  y_ordered <- y[order_id]
  fit <- stats::isoreg(p_ordered, y_ordered)
  fitted <- pmin(pmax(fit$yf, 0), 1)

  runs <- rle(p_ordered)
  ends <- cumsum(runs$lengths)
  starts <- c(1L, ends[-length(ends)] + 1L)
  y_fitted <- numeric(length(runs$values))
  for (i in seq_along(runs$values)) {
    y_fitted[i] <- mean(fitted[starts[i]:ends[i]])
  }

  new_calibrator(
    "cal_isotonic",
    method = "isotonic regression",
    n = length(y),
    input = "probabilities",
    x_thresholds = runs$values,
    y_calibrated = y_fitted,
    fit = fit,
    call = match.call()
  )
}

#' Histogram binning calibration
#'
#' `cal_histogram()` partitions `[0, 1]` into bins and replaces each probability
#' with the empirical event frequency in its bin. Equal-width bins use fixed
#' intervals. Equal-frequency bins use sample quantiles as break points.
#'
#' Empty training bins inherit the empirical rate from the nearest non-empty
#' bin. This makes prediction defined over the whole interval `[0, 1]`.
#'
#' @details
#' Histogram binning estimates a piecewise constant calibration map. Given
#' distinct break points
#' \eqn{0 = b_0 < b_1 < \cdots < b_J = 1}{0 = b_0 < ... < b_J = 1},
#' the implementation uses left-closed bins. For \eqn{j < J},
#'
#' \deqn{I_j = \{i: b_{j-1} \le p_i < b_j\},}{
#' I_j is the set with b_{j-1} <= p_i < b_j,}
#'
#' and the last bin is
#'
#' \deqn{I_J = \{i: b_{J-1} \le p_i \le b_J\}.}{
#' I_J is the set with b_{J-1} <= p_i <= b_J.}
#'
#' The fitted value for a non-empty bin is the empirical event frequency,
#'
#' \deqn{\hat q_j = \frac{1}{n_j}\sum_{i \in I_j} y_i,
#'   \quad n_j = |I_j|.}{q_hat_j = (1 / n_j) sum_{i in I_j} y_i.}
#'
#' A new probability receives the fitted value of the bin into which it falls.
#' Values exactly on an internal break point are assigned to the bin that starts
#' at that break point; the value `1` is assigned to the last bin.
#'
#' With `strategy = "equal_width"`, the break points are equally spaced on
#' `[0, 1]`, so \eqn{J = B}{J = B} when `bins = B`. With
#' `strategy = "equal_freq"`, provisional break points are
#'
#' \deqn{b_j = Q_8(j / B), \quad j = 0, \ldots, B,}{
#' b_j = Q_8(j / B), for j = 0, ..., B,}
#'
#' where \eqn{Q_8}{Q_8} is the sample quantile computed by
#' `stats::quantile(type = 8)`. The first and last break points are then forced
#' to `0` and `1`. Duplicated break points are removed, so the actual number of
#' bins \eqn{J} can be smaller than `bins`. Empty bins are assigned the value of
#' the nearest non-empty bin by bin index; if an empty bin is equally close to
#' two non-empty bins, the lower-index non-empty bin is used. If no non-empty
#' bin is available, the global event rate is used as a fallback.
#'
#' The returned object stores the requested `bins`, the realized `actual_bins`,
#' `strategy`, `breaks`, per-bin fitted values in `bin_values`, training
#' `counts`, `global_rate`, and the original call.
#'
#' @param p Numeric vector of uncalibrated probabilities in `[0, 1]`.
#' @param y Binary outcome vector coded as `0` and `1`.
#' @param bins Number of bins. Must be a single positive integer.
#' @param strategy Binning strategy. Use `"equal_width"` for fixed-width bins or
#' `"equal_freq"` for quantile bins.
#'
#' @return A `cal_histogram` object. Use `predict()` with new probabilities to
#' obtain calibrated probabilities.
#' @references
#' Zadrozny, B., & Elkan, C. (2002). Transforming classifier scores into
#' accurate multiclass probability estimates. Proceedings of the Eighth ACM
#' SIGKDD International Conference on Knowledge Discovery and Data Mining.
#' <doi:10.1145/775047.775151>.
#' @export
#'
#' @examples
#' set.seed(5)
#' calibration <- data.frame(raw_p = stats::runif(120)) |>
#'   dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))
#'
#' fit <- cal_histogram(calibration$raw_p, calibration$y, bins = 5)
#'
#' calibration |>
#'   dplyr::mutate(calibrated = predict(fit, raw_p)) |>
#'   dplyr::summarise(
#'     raw_ece = ece(raw_p, y, bins = 5),
#'     calibrated_ece = ece(calibrated, y, bins = 5)
#'   )
cal_histogram <- function(p, y, bins = 10, strategy = c("equal_width", "equal_freq")) {
  check_probability(p, arg = "p")
  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(p, y, x_arg = "p", y_arg = "y")
  bins <- check_bins(bins)
  strategy <- match.arg(strategy)

  breaks <- histogram_breaks(p, bins = bins, strategy = strategy)
  bin <- find_bins(p, breaks)
  n_bins <- length(breaks) - 1L
  counts <- tabulate(bin, nbins = n_bins)
  values <- rep(NA_real_, n_bins)

  for (i in seq_len(n_bins)) {
    idx <- bin == i
    if (any(idx)) {
      values[i] <- mean(y[idx])
    }
  }

  global_rate <- mean(y)
  values <- fill_empty_bins(values, counts, fallback = global_rate)

  new_calibrator(
    "cal_histogram",
    method = "histogram binning",
    n = length(y),
    input = "probabilities",
    bins = bins,
    actual_bins = n_bins,
    strategy = strategy,
    breaks = breaks,
    bin_values = values,
    counts = counts,
    global_rate = global_rate,
    call = match.call()
  )
}

histogram_breaks <- function(p, bins, strategy) {
  if (identical(strategy, "equal_width")) {
    return(seq(0, 1, length.out = bins + 1L))
  }

  probs <- seq(0, 1, length.out = bins + 1L)
  breaks <- stats::quantile(p, probs = probs, type = 8, names = FALSE)
  breaks[1L] <- 0
  breaks[length(breaks)] <- 1
  sort(unique(breaks))
}
