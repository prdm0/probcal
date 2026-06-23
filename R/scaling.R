#' Platt scaling
#'
#' `cal_platt()` fits a logistic regression that maps an uncalibrated score to
#' a calibrated probability. The binary targets are adjusted with Platt's target
#' correction before fitting, which shrinks labels away from exact `0` and `1`.
#'
#' @details
#' Let \eqn{(x_i, y_i), i = 1, \ldots, n}{(x_i, y_i), i = 1, ..., n}
#' be the calibration sample, where \eqn{x_i} is the supplied score and
#' \eqn{y_i \in \{0, 1\}}{y_i is 0 or 1} is the observed label. Write
#' \eqn{n_+ = \sum_i y_i}{n_+ = sum_i y_i} and
#' \eqn{n_- = n - n_+}{n_- = n - n_+}. Platt's correction replaces the
#' binary labels by fractional targets. Positive labels use
#'
#' \deqn{t_+ = \frac{n_+ + 1}{n_+ + 2},}{t_+ = (n_+ + 1) / (n_+ + 2),}
#'
#' and negative labels use
#'
#' \deqn{t_- = \frac{1}{n_- + 2}.}{t_- = 1 / (n_- + 2).}
#'
#' Thus \eqn{t_i = t_+}{t_i = t_+} when \eqn{y_i = 1}{y_i = 1}
#' and \eqn{t_i = t_-}{t_i = t_-} when \eqn{y_i = 0}{y_i = 0}. The
#' fitted logistic map is
#'
#' \deqn{q_i(\alpha, \beta) =
#'   \operatorname{logit}^{-1}(\alpha + \beta x_i),}{
#' q_i(alpha, beta) = logit^{-1}(alpha + beta x_i),}
#'
#' and \eqn{(\alpha, \beta)}{(alpha, beta)} are estimated by minimizing the
#' binomial cross-entropy with the corrected fractional targets,
#'
#' \deqn{\ell(\alpha, \beta) =
#'   -\sum_{i = 1}^n \{t_i \log q_i(\alpha, \beta) +
#'   (1 - t_i) \log[1 - q_i(\alpha, \beta)]\}.}{
#' ell(alpha, beta) = -sum_i {t_i log q_i + (1 - t_i) log(1 - q_i)}.}
#'
#' The implementation fits this model with `stats::glm()` using the formula
#' `y_adj ~ x`. The label vector must contain at least one `0` and one `1`.
#' The returned object stores `coefficients`, where `(Intercept)` is
#' \eqn{\hat\alpha}{alpha_hat} and `x` is \eqn{\hat\beta}{beta_hat}, as
#' well as the full `glm` object in `fit` and the corrected targets
#' `target_pos` and `target_neg`. Prediction applies
#' \eqn{\operatorname{logit}^{-1}(\hat\alpha + \hat\beta x_{new})}{
#' logit^{-1}(alpha_hat + beta_hat x_new)} to new scores. The argument `x` may
#' be a score on any real-valued scale or a raw probability, but the fitted map
#' is always a logistic function of the supplied values. The slope is
#' unconstrained; the fitted map is increasing in `x` only when
#' \eqn{\hat\beta \ge 0}{beta_hat >= 0}.
#'
#' @param x Numeric vector of uncalibrated scores or raw probabilities.
#' @param y Binary outcome vector coded as `0` and `1`.
#'
#' @return A `cal_platt` object. Use `predict()` with new scores to obtain
#' calibrated probabilities.
#' @references
#' Platt, J. (1999). Probabilistic outputs for support vector machines and
#' comparisons to regularized likelihood methods. In Advances in Large Margin
#' Classifiers.
#' @export
#'
#' @examples
#' set.seed(1)
#' calibration <- data.frame(score = rnorm(120)) |>
#'   dplyr::mutate(
#'     truth = inv_logit(score),
#'     y = rbinom(dplyr::n(), 1, truth)
#'   )
#'
#' fit <- cal_platt(calibration$score, calibration$y)
#'
#' calibration |>
#'   dplyr::mutate(calibrated = predict(fit, score)) |>
#'   dplyr::summarise(ece = ece(calibrated, y, bins = 10))
cal_platt <- function(x, y) {
  check_numeric_vector(x, arg = "x")
  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(x, y, x_arg = "x", y_arg = "y")

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)
  target_pos <- (n_pos + 1) / (n_pos + 2)
  target_neg <- 1 / (n_neg + 2)
  y_adj <- ifelse(y == 1L, target_pos, target_neg)

  data <- data.frame(y_adj = y_adj, x = x)
  fit <- suppressWarnings(stats::glm(
    y_adj ~ x,
    data = data,
    family = stats::binomial()
  ))

  new_calibrator(
    "cal_platt",
    method = "Platt scaling",
    n = length(y),
    input = "scores or probabilities",
    coefficients = stats::coef(fit),
    fit = fit,
    target_pos = target_pos,
    target_neg = target_neg,
    call = match.call()
  )
}

#' Temperature scaling
#'
#' `cal_temperature()` estimates a single positive temperature parameter by
#' minimizing the negative log-likelihood. Inputs must be logits, not
#' probabilities. For binary probabilities, `logit()` gives the corresponding
#' logit. For strictly positive multiclass probability rows,
#' \eqn{z_{ik} = \log p_{ik}}{z_ik = log p_ik} is a valid softmax logit
#' representation, up to row-wise additive constants. If probabilities have
#' zero entries, the user must choose and supply a transformed logit matrix,
#' such as clipped log-probabilities. `cal_temperature()` does not accept or
#' clip probability matrices.
#'
#' The function handles both binary and multiclass problems through the type of
#' `logits`. A numeric vector triggers binary temperature scaling and the
#' calibrated probability is `inv_logit(logits / T)`. A numeric matrix with one
#' column per class triggers multiclass temperature scaling and the calibrated
#' probabilities are `softmax(logits / T)`. Because dividing every logit by the
#' same positive scalar preserves the row ordering and argmax, temperature
#' scaling leaves the predicted class unchanged apart from existing ties and
#' only sharpens or softens the probabilities.
#'
#' @details
#' In the binary case, let \eqn{z_i} be an uncalibrated logit. For a positive
#' temperature \eqn{T}, the calibrated event probability is
#'
#' \deqn{q_i(T) = \operatorname{logit}^{-1}(z_i / T).}{
#' q_i(T) = logit^{-1}(z_i / T).}
#'
#' The fitted temperature is found by a bounded one-dimensional optimization on
#' \eqn{[10^{-3}, 10^3]}{[1e-3, 1e3]}:
#'
#' \deqn{\hat T \in \arg\min_{10^{-3} \le T \le 10^3}
#'   -\sum_{i = 1}^n \{y_i \log q_i(T) +
#'   (1 - y_i) \log[1 - q_i(T)]\}.}{
#' T_hat minimizes the binary negative log-likelihood over [1e-3, 1e3].}
#'
#' In the multiclass case, let \eqn{z_{ik}} be the logit for class \eqn{k} and
#' observation \eqn{i}. The calibrated probabilities are
#'
#' \deqn{q_{ik}(T) =
#'   \frac{\exp(z_{ik} / T)}
#'        {\sum_{\ell = 1}^K \exp(z_{i\ell} / T)},}{
#' q_ik(T) = exp(z_ik / T) / sum_l exp(z_il / T),}
#'
#' and \eqn{T} is chosen by minimizing the average multiclass negative
#' log-likelihood over the same interval,
#'
#' \deqn{L(T) = -\frac{1}{n}\sum_{i = 1}^n \log q_{i y_i}(T).}{
#' L(T) = -(1 / n) sum_i log q_i,y_i(T).}
#'
#' For multiclass labels, column \eqn{k} of the logit matrix corresponds to
#' class code \eqn{k}. If `y` is a factor, the stored order of `levels(y)`
#' defines the column order. The numerical objective clips probabilities that
#' enter logarithms to `[1e-15, 1 - 1e-15]`. The optimization uses
#' `stats::optim()` with method `"Brent"` and initial value `1` on the bounded
#' interval above. The returned object stores `temperature`, the optimizer
#' `value`, and the optimizer `convergence` code; multiclass fits also store
#' `k` and `levels`.
#'
#' Values \eqn{T > 1} soften the probability vector, while values
#' \eqn{0 < T < 1} make it more concentrated. Dividing all class logits by the
#' same positive constant preserves their order, so the predicted class is
#' unchanged apart from ties already present in the logits.
#'
#' @param logits For binary calibration, a numeric vector of uncalibrated
#' logits. For multiclass calibration, a numeric matrix of logits with one row
#' per observation and one column per class.
#' @param y Outcome labels. For binary calibration, a vector coded as `0` and
#' `1`. For multiclass calibration, a factor or a vector of integer class codes
#' in `1:K`, where `K` is the number of columns of `logits`.
#'
#' @return A `cal_temperature` object. Use `predict()` with new logits to obtain
#' calibrated probabilities. Multiclass objects also inherit from
#' `cal_multiclass`.
#' @references
#' Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration of
#' modern neural networks. Proceedings of the 34th International Conference on
#' Machine Learning.
#' @export
#'
#' @examples
#' set.seed(2)
#' calibration <- data.frame(logits = rnorm(120)) |>
#'   dplyr::mutate(
#'     raw_p = inv_logit(logits),
#'     y = rbinom(dplyr::n(), 1, raw_p)
#'   )
#'
#' fit <- cal_temperature(calibration$logits, calibration$y)
#'
#' calibration |>
#'   dplyr::mutate(calibrated = predict(fit, logits)) |>
#'   dplyr::summarise(
#'     raw_ece = ece(raw_p, y, bins = 10),
#'     calibrated_ece = ece(calibrated, y, bins = 10)
#'   )
#'
#' # Multiclass temperature scaling with a logit matrix and integer labels.
#' set.seed(20)
#' logits <- matrix(rnorm(150 * 3), ncol = 3)
#' labels <- max.col(logits) # integer codes in 1:3
#' mc_fit <- cal_temperature(logits, labels)
#' head(predict(mc_fit, logits))
cal_temperature <- function(logits, y) {
  if (is.matrix(logits)) {
    return(cal_temperature_multiclass(logits, y, call = match.call()))
  }

  check_numeric_vector(logits, arg = "logits")
  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(logits, y, x_arg = "logits", y_arg = "y")

  objective <- function(temperature) {
    binary_log_loss(inv_logit(logits / temperature), y)
  }

  opt <- stats::optim(
    par = 1,
    fn = objective,
    method = "Brent",
    lower = 1e-3,
    upper = 1e3
  )

  temperature <- unname(opt$par)
  if (!is.finite(temperature) || temperature <= 0) {
    cli::cli_abort("The fitted temperature must be positive and finite.")
  }

  new_calibrator(
    "cal_temperature",
    method = "temperature scaling",
    n = length(y),
    input = "logits",
    temperature = temperature,
    value = opt$value,
    convergence = opt$convergence,
    call = match.call()
  )
}

cal_temperature_multiclass <- function(logits, y, call = NULL) {
  check_logit_matrix(logits, arg = "logits")
  label <- check_multiclass_y(y, n_classes = ncol(logits), arg = "y")
  if (nrow(logits) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg logits} and {.arg y} must have the same number of observations.")
  }

  objective <- function(temperature) {
    multiclass_log_loss(softmax(logits / temperature), label$codes)
  }

  opt <- stats::optim(
    par = 1,
    fn = objective,
    method = "Brent",
    lower = 1e-3,
    upper = 1e3
  )

  temperature <- unname(opt$par)
  if (!is.finite(temperature) || temperature <= 0) {
    cli::cli_abort("The fitted temperature must be positive and finite.")
  }

  object <- new_calibrator(
    "cal_temperature",
    method = "multiclass temperature scaling",
    n = nrow(logits),
    input = "logits (matrix)",
    temperature = temperature,
    value = opt$value,
    convergence = opt$convergence,
    k = label$k,
    levels = label$levels,
    call = call
  )
  class(object) <- c("cal_temperature", "cal_multiclass", "calibrator")
  object
}

#' Beta calibration
#'
#' `cal_beta()` fits the beta calibration model
#' `inv_logit(a * log(p) - b * log(1 - p) + c)`. Probabilities are clipped to
#' to have lower bound `eps` and upper bound `1 - eps` before taking logarithms.
#'
#' @details
#' Beta calibration treats the uncalibrated event probability \eqn{p_i} through
#' two log-transformed features. Before the transformation, probabilities are
#' clipped by
#'
#' \deqn{p_i^* = C_\epsilon(p_i) =
#'   \min\{\max(p_i, \epsilon), 1 - \epsilon\}.}{
#' p_i^* = C_eps(p_i) = min(max(p_i, eps), 1 - eps).}
#'
#' The calibrated probability is
#'
#' \deqn{q_i = \operatorname{logit}^{-1}
#'   \{a \log(p_i^*) - b \log(1 - p_i^*) + c\}.}{
#' q_i = logit^{-1}(a log(p_i^*) - b log(1 - p_i^*) + c).}
#'
#' The implementation fits an ordinary unpenalized binomial `glm()` with the
#' original binary labels, without Platt target correction. Its linear
#' predictor is
#'
#' \deqn{\eta_i = \gamma_0 + \gamma_1 \log(p_i^*) +
#'   \gamma_2 \log(1 - p_i^*).}{
#' eta_i = gamma_0 + gamma_1 log(p_i^*) + gamma_2 log(1 - p_i^*).}
#'
#' Equivalently, the fitted coefficients minimize the binomial cross-entropy
#'
#' \deqn{-\sum_{i = 1}^n \{y_i \log q_i +
#'   (1 - y_i) \log(1 - q_i)\}.}{
#' -sum_i {y_i log q_i + (1 - y_i) log(1 - q_i)}.}
#'
#' The beta-calibration parameters are the following reparameterization of the
#' fitted `glm()` coefficients:
#'
#' \deqn{\hat a = \hat\gamma_1, \quad
#'   \hat b = -\hat\gamma_2, \quad
#'   \hat c = \hat\gamma_0.}{
#' a_hat = gamma_1_hat, b_hat = -gamma_2_hat, c_hat = gamma_0_hat.}
#'
#' Thus prediction first computes
#' \eqn{p_{new}^* = C_\epsilon(p_{new})}{p_new^* = C_eps(p_new)} and then
#' evaluates
#'
#' \deqn{\hat q(p_{new}) = \operatorname{logit}^{-1}\{
#'   \hat a \log(p_{new}^*) - \hat b \log(1 - p_{new}^*) + \hat c\}.}{
#' q_hat(p_new) = logit^{-1}(a_hat log(p_new^*) -
#' b_hat log(1 - p_new^*) + c_hat).}
#'
#' The object element `coefficients` contains
#' \eqn{(\hat\gamma_0, \hat\gamma_1, \hat\gamma_2)}{gamma_0_hat,
#' gamma_1_hat, gamma_2_hat} from `glm()`, while `a`, `b`, and `c` contain the
#' reparameterized beta-calibration coefficients. Since
#' \eqn{d\eta_i / dp_i = a / p_i + b / (1 - p_i)}{d eta / dp = a / p + b / (1 - p)},
#' monotone increase on `(0, 1)` is guaranteed when
#' \eqn{a \ge 0}{a >= 0} and \eqn{b \ge 0}{b >= 0}. The implementation does
#' not impose these constraints.
#'
#' @param p Numeric vector of uncalibrated probabilities in `[0, 1]`.
#' @param y Binary outcome vector coded as `0` and `1`.
#' @param eps Clipping constant satisfying `0 < eps < 0.5`. Probabilities must
#' first be valid values in `[0, 1]`; values below `eps` and above `1 - eps`
#' are clipped before taking logarithms.
#'
#' @return A `cal_beta` object. Use `predict()` with new probabilities to obtain
#' calibrated probabilities.
#' @references
#' Kull, M., Silva Filho, T. M., & Flach, P. (2017). Beta calibration: A
#' well-founded and easily implemented improvement on logistic calibration for
#' binary classifiers. Electronic Journal of Statistics, 11(2), 5052-5080.
#' <doi:10.1214/17-EJS1338SI>.
#' @export
#'
#' @examples
#' set.seed(3)
#' calibration <- data.frame(raw_p = stats::rbeta(120, 2, 2)) |>
#'   dplyr::mutate(y = rbinom(dplyr::n(), 1, raw_p))
#'
#' fit <- cal_beta(calibration$raw_p, calibration$y)
#'
#' calibration |>
#'   dplyr::mutate(calibrated = predict(fit, raw_p)) |>
#'   dplyr::summarise(
#'     raw_ece = ece(raw_p, y, bins = 10),
#'     calibrated_ece = ece(calibrated, y, bins = 10)
#'   )
cal_beta <- function(p, y, eps = 1e-15) {
  check_probability(p, arg = "p")
  y <- check_binary_y(y, arg = "y", require_both = TRUE)
  check_same_length(p, y, x_arg = "p", y_arg = "y")

  p_clip <- clip_prob(p, eps = eps, arg = "p")
  data <- data.frame(
    y = y,
    log_p = log(p_clip),
    log_1mp = log1p(-p_clip)
  )

  fit <- stats::glm(
    y ~ log_p + log_1mp,
    data = data,
    family = stats::binomial()
  )

  coefficients <- stats::coef(fit)

  new_calibrator(
    "cal_beta",
    method = "beta calibration",
    n = length(y),
    input = "probabilities",
    coefficients = coefficients,
    a = unname(coefficients[["log_p"]]),
    b = -unname(coefficients[["log_1mp"]]),
    c = unname(coefficients[["(Intercept)"]]),
    eps = eps,
    fit = fit,
    call = match.call()
  )
}

#' Vector scaling
#'
#' `cal_vector_scaling()` is the multiclass generalization of temperature
#' scaling that gives each class its own scale and bias. It rescales a logit
#' matrix column by column and applies the softmax. With a single shared scale
#' and no bias it reduces to temperature scaling, so it is more flexible while
#' remaining cheap to fit.
#'
#' The calibrated probabilities are `softmax(s * logits + b)`, where `s` is a
#' length `K` vector of per-class scales applied column by column and `b` is a
#' length `K` vector of per-class biases. Parameters are estimated by minimizing
#' the average multiclass negative log-likelihood.
#'
#' @details
#' Let \eqn{z_{ik}} be the uncalibrated logit for observation \eqn{i} and class
#' \eqn{k}. Vector scaling estimates class-specific scales \eqn{s_k} and
#' intercepts \eqn{b_k}, then forms calibrated logits
#'
#' \deqn{\eta_{ik} = s_k z_{ik} + b_k.}{eta_ik = s_k z_ik + b_k.}
#'
#' The predicted probabilities are obtained with the softmax,
#'
#' \deqn{q_{ik} =
#'   \frac{\exp(\eta_{ik})}
#'        {\sum_{\ell = 1}^K \exp(\eta_{i\ell})}.}{
#' q_ik = exp(eta_ik) / sum_l exp(eta_il).}
#'
#' Parameters are estimated by minimizing
#'
#' \deqn{L(s, b) = -\frac{1}{n}\sum_{i = 1}^n \log q_{i y_i}.}{
#' L(s, b) = -(1 / n) sum_i log q_i,y_i.}
#'
#' For multiclass labels, column \eqn{k} of `logits` corresponds to class code
#' \eqn{k}; if `y` is a factor, column \eqn{k} corresponds to `levels(y)[k]`.
#' The implementation uses `stats::optim()` with method `"BFGS"`, analytic
#' gradients, initial scales \eqn{s_k = 1}{s_k = 1}, initial biases
#' \eqn{b_k = 0}{b_k = 0}, and `maxit = 500`. True-class probabilities entering
#' logarithms are clipped to `[1e-15, 1 - 1e-15]`. The returned object stores
#' `scale`, `bias`, the optimized average negative log-likelihood `value`, and
#' the optimizer `convergence` code.
#'
#' The scales are unconstrained in the fitted optimization, so a negative scale
#' is possible when it improves the likelihood on the calibration data. Unlike
#' temperature scaling, vector scaling can change the predicted class because
#' scales and biases vary by class. As with any softmax model, adding the same
#' constant to every class bias does not change the resulting probability
#' vector, so the fitted bias vector is identifiable only up to a common
#' additive constant.
#'
#' @param logits Numeric matrix of uncalibrated logits with one row per
#' observation and one column per class.
#' @param y A factor or a vector of integer class codes in `1:K`, where `K` is
#' the number of columns of `logits`.
#'
#' @return A `cal_vector_scaling` object that also inherits from
#' `cal_multiclass`. Use `predict()` with new logits to obtain calibrated
#' probabilities.
#' @references
#' Guo, C., Pleiss, G., Sun, Y., & Weinberger, K. Q. (2017). On calibration of
#' modern neural networks. Proceedings of the 34th International Conference on
#' Machine Learning.
#' @export
#'
#' @examples
#' set.seed(22)
#' logits <- matrix(rnorm(200 * 3), ncol = 3)
#' labels <- max.col(logits)
#' fit <- cal_vector_scaling(logits, labels)
#' head(predict(fit, logits))
cal_vector_scaling <- function(logits, y) {
  check_logit_matrix(logits, arg = "logits")
  label <- check_multiclass_y(y, n_classes = ncol(logits), arg = "y")
  if (nrow(logits) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg logits} and {.arg y} must have the same number of observations.")
  }

  n <- nrow(logits)
  k <- ncol(logits)
  targets <- one_hot(label$codes, k)

  calibrated_logits <- function(par) {
    scale <- par[seq_len(k)]
    bias <- par[k + seq_len(k)]
    sweep(logits, 2L, scale, "*") + matrix(bias, nrow = n, ncol = k, byrow = TRUE)
  }

  objective <- function(par) {
    multiclass_log_loss(softmax(calibrated_logits(par)), label$codes)
  }

  gradient <- function(par) {
    prob <- softmax(calibrated_logits(par))
    resid <- (prob - targets) / n
    grad_scale <- colSums(resid * logits)
    grad_bias <- colSums(resid)
    c(grad_scale, grad_bias)
  }

  opt <- stats::optim(
    c(rep(1, k), rep(0, k)),
    objective,
    gradient,
    method = "BFGS",
    control = list(maxit = 500)
  )

  object <- new_calibrator(
    "cal_vector_scaling",
    method = "vector scaling",
    n = n,
    input = "logits (matrix)",
    scale = opt$par[seq_len(k)],
    bias = opt$par[k + seq_len(k)],
    value = opt$value,
    convergence = opt$convergence,
    k = k,
    levels = label$levels,
    call = match.call()
  )
  class(object) <- c("cal_vector_scaling", "cal_multiclass", "calibrator")
  object
}

#' Dirichlet calibration
#'
#' `cal_dirichlet()` is the multiclass generalization of beta calibration. It
#' fits a linear map on the log of the predicted probabilities followed by a
#' softmax, which is equivalent to a multinomial logistic regression with the
#' log-probabilities as features. An off-diagonal and intercept regularization
#' (ODIR) penalty shrinks the off-diagonal weights and the intercepts toward
#' zero, which reduces overfitting risk when the number of classes is large.
#'
#' The calibrated probabilities are computed row-wise as
#' `softmax(log(p) %*% t(W) + b)`, where `W` is a `K` by `K` weight matrix and
#' `b` is a length `K` intercept vector. Probabilities are clipped to
#' to have lower bound `eps` and upper bound `1 - eps` before taking logarithms.
#' When `lambda` is `NULL`, it is selected from a small deterministic grid by
#' cross-validated log-likelihood.
#'
#' @details
#' Let \eqn{p_{ik}} be the uncalibrated probability assigned to class \eqn{k}
#' for observation \eqn{i}. Each row of `p` must sum to one within absolute
#' tolerance `1e-6`. Column \eqn{k} corresponds to integer class code \eqn{k};
#' if `y` is a factor, column \eqn{k} corresponds to `levels(y)[k]`. The
#' entries are clipped elementwise by
#'
#' \deqn{p_{ik}^* = \min\{\max(p_{ik}, \epsilon), 1 - \epsilon\},}{
#' p_ik^* = min(max(p_ik, eps), 1 - eps),}
#'
#' and transformed to \eqn{u_{ik} = \log(p_{ik}^*)}. The clipped feature matrix
#' is not renormalized; normalization occurs only after the linear map, through
#' the final softmax. Dirichlet calibration fits a multinomial logistic
#' regression on these log-probability features,
#'
#' \deqn{\eta_{ik} = b_k + \sum_{\ell = 1}^K W_{k\ell} u_{i\ell},}{
#' eta_ik = b_k + sum_l W_kl u_il,}
#'
#' followed by
#'
#' \deqn{q_{ik} =
#'   \frac{\exp(\eta_{ik})}{\sum_{m = 1}^K \exp(\eta_{im})}.}{
#' q_ik = exp(eta_ik) / sum_m exp(eta_im).}
#'
#' With fixed \eqn{\lambda}, the fitted parameters minimize
#'
#' \deqn{-\frac{1}{n}\sum_i \log q_{i y_i}
#'   + \lambda\left(\sum_{k \ne \ell} W_{k\ell}^2
#'   + \sum_k b_k^2\right).}{
#' -(1 / n) sum_i log q_i,y_i + lambda(sum_{k != l} W_kl^2 + sum_k b_k^2).}
#'
#' This is the off-diagonal and intercept regularization penalty. Diagonal
#' weights are not penalized. For fixed `lambda`, optimization uses BFGS with
#' analytic gradients, initial weight matrix \eqn{W = I_K}{W = I_K}, initial
#' bias \eqn{b = 0}{b = 0}, and `maxit = 500`. True-class probabilities
#' entering logarithms are clipped to `[1e-15, 1 - 1e-15]`. The returned
#' `weight` is a \eqn{K \times K}{K by K} matrix whose row \eqn{k} produces
#' the logit for class \eqn{k}; `bias` is a length-\eqn{K} vector of
#' intercepts. The object also stores `lambda`, `value`, and the optimizer
#' `convergence` code.
#'
#' If `lambda = NULL`, the implementation evaluates the grid
#' `c(0, 1e-4, 1e-3, 1e-2, 1e-1)` with at most three deterministic stratified
#' folds. Class indices are assigned to folds in their existing order. The
#' selected value minimizes the unweighted average of the fold mean held-out
#' negative log-likelihoods; ties choose the first grid value. If fewer than two
#' observations are available in the smallest class during selection, the
#' fallback value is `1e-3`. With `lambda = 0`, the multinomial softmax
#' parameterization is not unique: adding the same linear function of the
#' features to every class logit leaves all probabilities unchanged. The
#' calibrated probabilities are the identified output.
#'
#' @param p Numeric matrix of uncalibrated probabilities with one row per
#' observation and one column per class. Rows must sum to one within absolute
#' tolerance `1e-6`.
#' @param y A factor or a vector of integer class codes in `1:K`, where `K` is
#' the number of columns of `p`.
#' @param lambda Non-negative ODIR regularization strength. When `NULL` it is
#' chosen by cross-validation.
#' @param eps Clipping constant satisfying `0 < eps < 0.5`. Probabilities must
#' first be valid values in `[0, 1]`; values below `eps` and above `1 - eps`
#' are clipped before taking logarithms.
#'
#' @return A `cal_dirichlet` object that also inherits from `cal_multiclass`.
#' Use `predict()` with new probabilities to obtain calibrated probabilities.
#' @references
#' Kull, M., Perello-Nieto, M., Kängsepp, M., Silva Filho, T., Song, H., &
#' Flach, P. (2019). Beyond temperature scaling: Obtaining well-calibrated
#' multi-class probabilities with Dirichlet calibration. Advances in Neural
#' Information Processing Systems 32.
#' @export
#'
#' @examples
#' set.seed(23)
#' prob <- matrix(stats::runif(200 * 3), ncol = 3)
#' prob <- prob / rowSums(prob)
#' labels <- max.col(prob)
#' fit <- cal_dirichlet(prob, labels)
#' head(predict(fit, prob))
cal_dirichlet <- function(p, y, lambda = NULL, eps = 1e-12) {
  check_prob_matrix(p, arg = "p")
  label <- check_multiclass_y(y, n_classes = ncol(p), arg = "y")
  if (nrow(p) != length(label$codes)) {
    cli::cli_abort("Arguments {.arg p} and {.arg y} must have the same number of observations.")
  }
  if (!is.null(lambda) && (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda < 0)) {
    cli::cli_abort("Argument {.arg lambda} must be `NULL` or a single non-negative number.")
  }

  k <- ncol(p)
  features <- log(clip_prob(p, eps = eps, arg = "p"))

  if (is.null(lambda)) {
    lambda <- select_dirichlet_lambda(features, label$codes, k)
  }

  fit <- fit_multinomial(features, label$codes, k, lambda = lambda)

  object <- new_calibrator(
    "cal_dirichlet",
    method = "Dirichlet calibration",
    n = nrow(p),
    input = "probabilities (matrix)",
    weight = fit$weight,
    bias = fit$bias,
    lambda = lambda,
    eps = eps,
    value = fit$value,
    convergence = fit$convergence,
    k = k,
    levels = label$levels,
    call = match.call()
  )
  class(object) <- c("cal_dirichlet", "cal_multiclass", "calibrator")
  object
}

# Choose the ODIR strength by minimizing the average held-out log-likelihood
# over a small grid. Folds are assigned deterministically within each class so
# the selection does not depend on the random number generator.
select_dirichlet_lambda <- function(features, codes, k,
                                     grid = c(0, 1e-4, 1e-3, 1e-2, 1e-1)) {
  counts <- table(codes)
  folds <- min(3L, min(counts))
  if (folds < 2L) {
    return(1e-3)
  }

  fold_id <- integer(length(codes))
  for (cls in sort(unique(codes))) {
    idx <- which(codes == cls)
    fold_id[idx] <- rep_len(seq_len(folds), length(idx))
  }

  scores <- vapply(grid, function(lambda) {
    total <- 0
    for (f in seq_len(folds)) {
      train <- fold_id != f
      test <- fold_id == f
      fit <- fit_multinomial(features[train, , drop = FALSE], codes[train], k, lambda = lambda)
      logits <- features[test, , drop = FALSE] %*% t(fit$weight) +
        matrix(fit$bias, nrow = sum(test), ncol = k, byrow = TRUE)
      total <- total + multiclass_log_loss(softmax(logits), codes[test])
    }
    total / folds
  }, numeric(1))

  grid[which.min(scores)]
}
