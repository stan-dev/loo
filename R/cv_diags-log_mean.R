#' Truncated-GPD MCSE and Jensen bias of a log mean
#'
#' The estimate is \eqn{\log \bar z} with \eqn{\bar z = S^{-1} \sum_s z^s}.
#' The upper tail of \eqn{z^s} gets a generalized Pareto fit. The variance
#' of \eqn{z} combines the empirical bulk with the moments of the fitted
#' tail, truncated at the smaller of the bound `upper` and the
#' finite-sample bound \eqn{S^{\hat k}}.
#'
#' @noRd
#' @param x Draws of \eqn{z^s}: a vector, or an iterations-by-chains matrix.
#' @param are_log `TRUE` if `x` holds \eqn{\log z^s}.
#' @param upper Upper bound of \eqn{z^s}, on the scale of \eqn{z}. `Inf`
#'   leaves only the finite-sample bound.
#' @param r_eff Relative efficiency. It sets only the tail length. `NULL`
#'   computes it from the chains of `x`; 1 for a vector.
#' @return A list with `mcse`, `bias` (\eqn{-\mathrm{mcse}^2/2}), `pareto_k`,
#'   `sigma`, `threshold`, `tau_eff` and `truncation` (`"none"`,
#'   `"sample_size"` or `"physical"`). `sigma` and `threshold` are on the
#'   scale of \eqn{z / \max_s z^s}. All values are `NA` when \eqn{S < 50},
#'   when `x` is not finite, or when the tail cannot be fitted.
#'
.log_mean_diags <- function(x, are_log = FALSE, upper = Inf, r_eff = NULL) {
  log_z <- if (are_log) x else log(x)
  S <- length(log_z)
  if (S < 50 || !all(is.finite(log_z))) {
    return(.log_mean_diags_na())
  }
  if (is.null(r_eff)) {
    r_eff <- if (is.matrix(log_z)) {
      relative_eff(array(exp(log_z - max(log_z)), c(dim(log_z), 1)))
    } else {
      1
    }
    if (!is.finite(r_eff)) r_eff <- 1
  }

  log_max <- max(log_z)
  z <- sort(exp(as.vector(log_z) - log_max))
  fit <- .fit_right_tail(z, n_pareto(r_eff, S))
  if (is.null(fit)) {
    return(.log_mean_diags_na())
  }
  upper <- exp(log(upper) - log_max)
  if (upper <= fit$threshold) {
    upper <- Inf
  }
  tail_moments <- .truncated_gpd_moments(fit$k, fit$sigma, fit$threshold, upper, S)

  share <- fit$tail_len / S
  bulk <- z[seq_len(S - fit$tail_len)]
  mean_bulk <- mean(bulk)
  var_bulk <- mean((bulk - mean_bulk)^2)
  mean_z <- (1 - share) * mean_bulk + share * tail_moments$mean
  var_z <- (1 - share) * var_bulk + share * tail_moments$var +
    share * (1 - share) * (mean_bulk - tail_moments$mean)^2

  mcse <- sqrt(var_z / S) / mean_z
  if (!is.finite(mcse) || mean_z <= 0) {
    return(.log_mean_diags_na())
  }
  list(
    mcse = mcse,
    bias = -mcse^2 / 2,
    pareto_k = fit$k,
    sigma = fit$sigma,
    threshold = fit$threshold,
    tau_eff = tail_moments$tau_eff,
    truncation = tail_moments$truncation
  )
}

#' The `NA` result of `.log_mean_diags()`
#' @noRd
.log_mean_diags_na <- function() {
  list(mcse = NA_real_, bias = NA_real_, pareto_k = NA_real_, sigma = NA_real_,
       threshold = NA_real_, tau_eff = NA_real_, truncation = NA_character_)
}

#' Pareto \eqn{\hat k} of the right tail of `exp(lw)`
#'
#' The PSIS \eqn{\hat k} of `do_psis_i()`; `NA` when it cannot be estimated.
#'
#' @noRd
#' @param lw Log of the draws.
#' @param r_eff Relative efficiency; sets the tail length.
#'
.right_tail_khat <- function(lw, r_eff = 1) {
  if (!all(is.finite(lw))) {
    return(NA_real_)
  }
  khat <- do_psis_i(lw, n_pareto(r_eff, length(lw)))$pareto_k
  if (is.finite(khat)) khat else NA_real_
}

#' GPD fit to the largest draws
#'
#' @noRd
#' @param z Draws, sorted in increasing order.
#' @param tail_len Number of draws in the tail.
#' @return A list with `k`, `sigma`, `threshold` and `tail_len`; `NULL` when
#'   the tail is too short or the fit fails.
#'
.fit_right_tail <- function(z, tail_len) {
  S <- length(z)
  if (!enough_tail_samples(tail_len) || tail_len >= S) {
    return(NULL)
  }
  threshold <- z[S - tail_len]
  exceedances <- z[seq(S - tail_len + 1, S)] - threshold
  exceedances <- exceedances[exceedances > 0]
  if (!enough_tail_samples(length(exceedances))) {
    return(NULL)
  }
  fit <- tryCatch(posterior::gpdfit(exceedances, sort_x = FALSE),
                  error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$k) || !is.finite(fit$sigma) || fit$sigma <= 0) {
    return(NULL)
  }
  list(k = fit$k, sigma = fit$sigma, threshold = threshold, tail_len = tail_len)
}

#' Mean and variance of a GPD tail, truncated at an upper bound
#'
#' On the standardized scale \eqn{X = 1 + k (z - u) / \sigma}, the tail is
#' \eqn{\mathrm{Pareto}(1, 1/k)}. The truncation point is
#' \eqn{\tau_{\mathrm{eff}} = \min(1 + k (\mathrm{upper} - u) / \sigma,\ S^k)}.
#' Near \eqn{k = 0} the tail is exponential, truncated at `upper`. For a
#' clearly negative \eqn{k} the tail is bounded; the untruncated moments
#' apply.
#'
#' @noRd
#' @param k,sigma GPD shape and scale.
#' @param threshold The threshold \eqn{u}.
#' @param upper Upper bound of \eqn{z}, above the threshold, or `Inf`.
#' @param S Number of draws.
#' @return A list with `mean`, `var`, `tau_eff` and `truncation`.
#'
.truncated_gpd_moments <- function(k, sigma, threshold, upper, S) {
  if (abs(k) < 1e-3) {
    return(.truncated_exponential_moments(sigma, threshold, upper))
  }
  if (k < 0) {
    return(list(
      mean = threshold + sigma / (1 - k),
      var = sigma^2 / ((1 - k)^2 * (1 - 2 * k)),
      tau_eff = Inf,
      truncation = "none"
    ))
  }
  tau_physical <- 1 + k * (upper - threshold) / sigma
  tau_sample <- S^k
  tau_eff <- min(tau_physical, tau_sample)
  moments <- .truncated_pareto_moments(alpha = 1 / k, log_tau = log(tau_eff))
  list(
    mean = threshold + sigma / k * (moments$mean - 1),
    var = (sigma / k)^2 * moments$var,
    tau_eff = tau_eff,
    truncation = if (tau_physical < tau_sample) "physical" else "sample_size"
  )
}

#' Mean and variance of an exponential tail, truncated at `upper`
#'
#' The limit \eqn{k \to 0} of the GPD. With \eqn{a = (\mathrm{upper} - u) /
#' \sigma}, the excess has mean \eqn{\sigma (1 - (1 + a) e^{-a}) / (1 - e^{-a})}.
#' For \eqn{a > 700} the truncation has no effect.
#'
#' @noRd
#'
.truncated_exponential_moments <- function(sigma, threshold, upper) {
  a <- (upper - threshold) / sigma
  if (!is.finite(a) || a > 700) {
    return(list(mean = threshold + sigma, var = sigma^2, tau_eff = Inf,
                truncation = "none"))
  }
  mass <- -expm1(-a)
  m1 <- sigma * (1 - (1 + a) * exp(-a)) / mass
  m2 <- sigma^2 * (2 - (2 + 2 * a + a^2) * exp(-a)) / mass
  list(mean = threshold + m1, var = max(m2 - m1^2, 0), tau_eff = Inf,
       truncation = "physical")
}

#' Mean and variance of Pareto(1, alpha), truncated at tau
#'
#' \eqn{E[X^m \mid X \le \tau] = \alpha\, g(\alpha - m) / (1 - \tau^{-\alpha})}
#' with \eqn{g(a) = (1 - \tau^{-a}) / a}, and \eqn{g(0) = \log \tau}.
#'
#' @noRd
#' @param alpha Pareto shape, \eqn{1 / k}.
#' @param log_tau \eqn{\log \tau}, positive and finite.
#' @return A list with `mean` and `var` of \eqn{X}.
#'
.truncated_pareto_moments <- function(alpha, log_tau) {
  g <- function(a) {
    if (abs(a * log_tau) < 1e-8) log_tau else -expm1(-a * log_tau) / a
  }
  mass <- -expm1(-alpha * log_tau)
  m1 <- alpha * g(alpha - 1) / mass
  m2 <- alpha * g(alpha - 2) / mass
  list(mean = m1, var = max(m2 - m1^2, 0))
}
