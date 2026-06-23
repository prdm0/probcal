#' Reliability diagram
#'
#' `reliability_diagram()` returns a `ggplot2` object comparing mean predicted
#' confidence with the observed event frequency in equal-width probability bins.
#' By default, points are sized by the number of observations in each non-empty
#' bin and the subtitle reports the ECE computed with the same bins.
#'
#' For a probability matrix the function builds a multiclass diagram. The
#' `"classwise"` form draws one panel per class from the one-vs-rest view. The
#' `"confidence"` form draws a single panel from the top-label confidence and
#' whether the predicted class is correct.
#'
#' @details
#' The diagram is a visual version of the binned summaries used by [ece()]. For
#' binary inputs, the package uses the same left-closed equal-width bins as
#' [ece()], with the last bin closed on the right. For each non-empty bin
#' \eqn{b}, the x-coordinate is the mean predicted probability,
#'
#' \deqn{\operatorname{conf}(b) = \frac{1}{n_b}\sum_{i \in I_b} p_i,}{
#' conf(b) = (1 / n_b) sum_{i in I_b} p_i,}
#'
#' and the y-coordinate is the observed event frequency,
#'
#' \deqn{\operatorname{acc}(b) = \frac{1}{n_b}\sum_{i \in I_b} y_i.}{
#' acc(b) = (1 / n_b) sum_{i in I_b} y_i.}
#'
#' Points near the diagonal line have similar average confidence and empirical
#' frequency within the bin. Points below the diagonal indicate over-confident
#' predictions in that bin, and points above the diagonal indicate
#' under-confident predictions. Empty bins are omitted from the plotted data.
#' The diagonal reference line is the set where the bin mean predicted
#' probability equals the empirical event frequency.
#'
#' For multiclass inputs, `type = "classwise"` builds these summaries separately
#' for each one-vs-rest class and displays them in facets. `type = "confidence"`
#' replaces \eqn{p_i} by the top-label probability and \eqn{y_i} by the
#' indicator that the top-label prediction is correct. Ties in the top-label
#' rule are broken by the first column, matching `max.col(..., ties.method =
#' "first")`. When `show_ece = TRUE`, the subtitle reports
#' `ece(p, y, bins = bins)` for binary inputs and
#' `ece(p, y, bins = bins, type = type)` for multiclass inputs.
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
#' @param show_ece Logical. If `TRUE`, include the ECE in the plot subtitle.
#' @param show_counts Logical. If `TRUE`, map point size to the number of
#' observations in each bin.
#' @param type Multiclass layout, either `"classwise"` or `"confidence"`.
#' Ignored for binary inputs.
#'
#' @return A `ggplot` object.
#' @references
#' Niculescu-Mizil, A., & Caruana, R. (2005). Predicting good probabilities
#' with supervised learning. Proceedings of the 22nd International Conference
#' on Machine Learning.
#' @export
#'
#' @examples
#' set.seed(6)
#' predictions <- data.frame(raw_p = stats::runif(120)) |>
#'   dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))
#'
#' reliability_diagram(predictions$raw_p, predictions$y, bins = 8)
#'
#' # Multiclass reliability diagram with one panel per class.
#' set.seed(60)
#' prob <- matrix(stats::runif(150 * 3), ncol = 3)
#' prob <- prob / rowSums(prob)
#' labels <- max.col(prob)
#' reliability_diagram(prob, labels, bins = 8, type = "classwise")
reliability_diagram <- function(
  p,
  y,
  bins = 10,
  show_ece = TRUE,
  show_counts = TRUE,
  type = c("classwise", "confidence")
) {
  if (is.matrix(p)) {
    return(reliability_diagram_multiclass(
      p, y, bins, show_ece, show_counts, match.arg(type)
    ))
  }

  stats <- bin_stats(p, y, bins)
  stats <- stats[stats$n > 0L, , drop = FALSE]
  subtitle <- NULL
  if (isTRUE(show_ece)) {
    subtitle <- sprintf("ECE = %.4f", ece(p, y, bins = bins))
  }

  plot <- ggplot2::ggplot(stats, ggplot2::aes(x = confidence, y = accuracy)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Mean predicted probability",
      y = "Observed event frequency",
      title = "Reliability diagram",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(show_counts)) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.8) +
      ggplot2::scale_size_continuous(name = "Count")
  } else {
    plot <- plot + ggplot2::geom_point(size = 2.5, alpha = 0.8)
  }

  plot
}

reliability_diagram_multiclass <- function(p, y, bins, show_ece, show_counts, type) {
  check_prob_matrix(p, arg = "p")

  if (identical(type, "confidence")) {
    conf <- multiclass_confidence(p, y)
    plot <- reliability_diagram(
      conf$confidence, conf$correct, bins = bins,
      show_ece = show_ece, show_counts = show_counts
    )
    return(plot + ggplot2::labs(
      title = "Reliability diagram (confidence)",
      x = "Top-label confidence",
      y = "Observed accuracy"
    ))
  }

  label <- check_multiclass_y(y, n_classes = ncol(p), arg = "y")
  panels <- lapply(seq_len(ncol(p)), function(k) {
    stats <- bin_stats(p[, k], as.integer(label$codes == k), bins)
    stats <- stats[stats$n > 0L, , drop = FALSE]
    stats$class_label <- factor(label$levels[k], levels = label$levels)
    stats
  })
  stats <- do.call(rbind, panels)

  subtitle <- NULL
  if (isTRUE(show_ece)) {
    subtitle <- sprintf("Classwise ECE = %.4f", ece(p, y, bins = bins, type = "classwise"))
  }

  plot <- ggplot2::ggplot(stats, ggplot2::aes(x = confidence, y = accuracy)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    ggplot2::facet_wrap(ggplot2::vars(class_label)) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Mean predicted probability",
      y = "Observed event frequency",
      title = "Reliability diagram (classwise)",
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(show_counts)) {
    plot <- plot +
      ggplot2::geom_point(ggplot2::aes(size = n), alpha = 0.8) +
      ggplot2::scale_size_continuous(name = "Count")
  } else {
    plot <- plot + ggplot2::geom_point(size = 2.5, alpha = 0.8)
  }

  plot
}
