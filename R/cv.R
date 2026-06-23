#' Cross-validated calibration
#'
#' `cal_cv()` fits a calibrator with out-of-fold predictions. The function
#' expects scores, probabilities, or logits that were already produced by a
#' model. It does not train the underlying classifier.
#'
#' Folds are stratified by the outcome. The returned object stores the
#' out-of-fold calibrated probabilities and a final calibrator fitted on all
#' observations for future prediction. Binary and multiclass problems are
#' handled through the type of `x`. A numeric vector triggers binary
#' calibration. A numeric matrix with one column per class triggers multiclass
#' calibration, the out-of-fold predictions become a matrix, and the available
#' methods are `"temperature"`, `"vector"`, `"dirichlet"`, and `"ovr"`. For
#' `method = "ovr"`, pass the binary method through `base_method`.
#'
#' @details
#' Cross-validated calibration estimates how the calibration map behaves on
#' observations not used to fit that map. Let \eqn{F_i \in \{1, \ldots, V\}}
#' denote the fold assigned to observation \eqn{i}. For each fold \eqn{v}, a
#' calibrator \eqn{\hat f^{(-v)}} is fitted using observations with
#' \eqn{F_i \ne v}. The out-of-fold calibrated prediction for an observation in
#' fold \eqn{v} is then
#'
#' \deqn{\hat q_i^{\mathrm{oof}} = \hat f^{(-v)}(x_i),
#'   \quad F_i = v.}{q_hat_i^oof = f_hat^(-v)(x_i), for F_i = v.}
#'
#' These out-of-fold predictions are stored in `oof_predictions` and are useful
#' for estimating calibration metrics without evaluating a calibrator on the
#' same observations used to fit it. In binary calibration,
#' \eqn{\hat q_i^{\mathrm{oof}}}{q_hat_i^oof} is a scalar event probability.
#' In multiclass calibration, it is the row vector
#' \eqn{(\hat q_{i1}^{\mathrm{oof}}, \ldots,
#' \hat q_{iK}^{\mathrm{oof}})}{(q_hat_i1^oof, ..., q_hat_iK^oof)} on the
#' probability simplex. After the out-of-fold predictions are computed, a final
#' calibrator \eqn{\hat f} is fitted on all observations. The S3 `predict()`
#' method for a `cal_cv` object uses this final calibrator for future data.
#'
#' The folds are stratified by the observed labels. Setting `seed` affects only
#' the fold assignment and restores the previous random-number state after the
#' assignment is made. The function assumes that `x` already contains model
#' outputs from another classifier; it does not refit that classifier inside
#' each fold. Thus the predictions are out of fold for the calibration map only,
#' unless `x` itself was produced out of fold by the underlying classifier.
#'
#' `folds` must be at least `2` and no larger than the smallest class count.
#' Within each class, observations are randomly permuted and assigned fold
#' labels \eqn{1, \ldots, V, 1, \ldots}{1, ..., V, 1, ...} in sequence. For
#' multiclass inputs, column \eqn{k} corresponds to integer class code \eqn{k};
#' if `y` is a factor, column \eqn{k} corresponds to `levels(y)[k]`. For
#' `method = "ovr"`, `base_method` is read from `...`; if it is not supplied,
#' the default base method is `"platt"`.
#'
#' @param x Numeric vector of uncalibrated values for binary calibration, or a
#' numeric matrix with one column per class for multiclass calibration. Use
#' logits for `method = "temperature"` and `"vector"`, probabilities for
#' `"beta"`, `"isotonic"`, `"histogram"`, and `"dirichlet"`, and scores or
#' probabilities for `"platt"`.
#' @param y Binary outcome vector coded as `0` and `1`, or a factor or vector of
#' integer class codes in `1:K` for multiclass calibration.
#' @param method Calibration method.
#' @param folds Number of stratified folds. Must be a single integer at least
#' `2` and no larger than the smallest class count.
#' @param seed Optional integer seed used only for fold assignment.
#' @param ... Additional arguments passed to the selected calibrator, such as
#' `bins` for histogram binning or `base_method` for one-vs-rest calibration.
#'
#' @return A `cal_cv` object. Use `predict()` to apply the final calibrator to
#' new values. The object stores `fold_id`, `oof_predictions`,
#' `fold_calibrators`, and `final_calibrator`. For binary calibration,
#' `oof_predictions` is a numeric vector. For multiclass calibration, it is a
#' numeric matrix with one row per observation and one column per class, with
#' column names given by the class levels.
#' @export
#'
#' @examples
#' set.seed(7)
#' predictions <- data.frame(raw_p = stats::runif(120)) |>
#'   dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))
#'
#' fit <- cal_cv(
#'   predictions$raw_p,
#'   predictions$y,
#'   method = "histogram",
#'   folds = 3,
#'   bins = 5,
#'   seed = 1
#' )
#'
#' predictions |>
#'   dplyr::mutate(calibrated = fit$oof_predictions) |>
#'   dplyr::summarise(ece = ece(calibrated, y, bins = 5))
cal_cv <- function(
  x,
  y,
  method = c("platt", "temperature", "beta", "isotonic", "histogram",
             "vector", "dirichlet", "ovr"),
  folds = 5,
  seed = NULL,
  ...
) {
  method <- match.arg(method)

  if (is.matrix(x)) {
    return(cal_cv_multiclass(x, y, method, folds, seed, ..., call = match.call()))
  }

  if (method %in% c("vector", "dirichlet", "ovr")) {
    cli::cli_abort("Method {.val {method}} requires a matrix {.arg x} with one column per class.")
  }

  check_numeric_vector(x, arg = "x")
  if (method %in% c("beta", "isotonic", "histogram")) {
    check_probability(x, arg = "x")
  }

  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(x, y, x_arg = "x", y_arg = "y")
  folds <- check_folds(folds, y)

  fold_id <- stratified_folds(y, folds = folds, seed = seed)
  oof_predictions <- numeric(length(y))
  fold_calibrators <- vector("list", folds)

  for (fold in seq_len(folds)) {
    assessment <- fold_id == fold
    analysis <- !assessment
    fold_calibrators[[fold]] <- fit_calibrator(method, x[analysis], y[analysis], ...)
    oof_predictions[assessment] <- stats::predict(
      fold_calibrators[[fold]],
      newdata = x[assessment]
    )
  }

  final_calibrator <- fit_calibrator(method, x, y, ...)

  new_calibrator(
    "cal_cv",
    method = paste("cross-validated", final_calibrator$method),
    n = length(y),
    input = final_calibrator$input,
    calibration_method = method,
    folds = folds,
    fold_id = fold_id,
    oof_predictions = oof_predictions,
    fold_calibrators = fold_calibrators,
    final_calibrator = final_calibrator,
    call = match.call()
  )
}

fit_calibrator <- function(method, x, y, ...) {
  switch(
    method,
    platt = cal_platt(x = x, y = y),
    temperature = cal_temperature(logits = x, y = y),
    beta = cal_beta(p = x, y = y, ...),
    isotonic = cal_isotonic(p = x, y = y),
    histogram = cal_histogram(p = x, y = y, ...)
  )
}

cal_cv_multiclass <- function(x, y, method, folds, seed, ..., call = NULL) {
  if (method %in% c("platt", "beta", "isotonic", "histogram")) {
    cli::cli_abort("Method {.val {method}} expects a vector {.arg x}. Use {.val ovr} for multiclass binning.")
  }

  if (method %in% c("temperature", "vector")) {
    check_logit_matrix(x, arg = "x")
  } else if (method == "dirichlet") {
    check_prob_matrix(x, arg = "x")
  } else if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    cli::cli_abort("Argument {.arg x} must be a finite numeric matrix.")
  }

  label <- check_multiclass_y(y, n_classes = ncol(x), arg = "y")
  if (nrow(x) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg x} and {.arg y} must have the same number of observations.")
  }
  folds <- check_folds(folds, label$codes)

  fold_id <- stratified_folds(label$codes, folds = folds, seed = seed)
  oof_predictions <- matrix(NA_real_, nrow = nrow(x), ncol = ncol(x))
  fold_calibrators <- vector("list", folds)

  for (fold in seq_len(folds)) {
    assessment <- fold_id == fold
    analysis <- !assessment
    fold_calibrators[[fold]] <- fit_calibrator_mc(
      method, x[analysis, , drop = FALSE], label$codes[analysis], ...
    )
    oof_predictions[assessment, ] <- stats::predict(
      fold_calibrators[[fold]],
      newdata = x[assessment, , drop = FALSE]
    )
  }

  final_calibrator <- fit_calibrator_mc(method, x, y, ...)
  colnames(oof_predictions) <- label$levels

  object <- new_calibrator(
    "cal_cv",
    method = paste("cross-validated", final_calibrator$method),
    n = nrow(x),
    input = final_calibrator$input,
    calibration_method = method,
    folds = folds,
    fold_id = fold_id,
    oof_predictions = oof_predictions,
    fold_calibrators = fold_calibrators,
    final_calibrator = final_calibrator,
    k = label$k,
    levels = label$levels,
    call = call
  )
  class(object) <- c("cal_cv", "cal_multiclass", "calibrator")
  object
}

fit_calibrator_mc <- function(method, x, y, ...) {
  dots <- list(...)
  switch(
    method,
    temperature = cal_temperature(x, y),
    vector = cal_vector_scaling(x, y),
    dirichlet = do.call(cal_dirichlet, c(list(x, y), dots)),
    ovr = {
      base_method <- if (is.null(dots$base_method)) "platt" else dots$base_method
      dots$base_method <- NULL
      do.call(cal_ovr, c(list(x, y, method = base_method), dots))
    }
  )
}

stratified_folds <- function(y, folds, seed = NULL) {
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  fold_id <- integer(length(y))
  for (class in sort(unique(y))) {
    idx <- sample(which(y == class))
    fold_id[idx] <- rep(seq_len(folds), length.out = length(idx))
  }

  fold_id
}
