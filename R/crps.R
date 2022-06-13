#' Continuously ranked probability score
#'
#' The `crps()` and `scrps()` and their `loo` counterparts can be used to
#' compute the continuously ranked probability score (CRPS) and scaled CRPS
#' (see Bolin and Wallin, 2022).
#'
#' @export
#' @param x1 A.default or an array of posterior predictive draws.
#' @param x2 Draws from the same distribution as draws in `x1`.
#'     Sholud be of identical dimension.
#' @param y A vector of observations
#' @param ... Passed on to `E_loo()` in `loo_*`-functions.
#'
#' @references
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules,
#' Prediction, and Estimation. Journal of the American Statistical Association,
#' 102(477), 359–378.
#'
#' Bolin, D., & Wallin, J. (2022). Local scale invariance and robustness of
#' proper scoring rules. arXiv.
#' <https://doi.org/10.48550/arXiv.1912.05642>
crps <- function(x1, ...) {
  UseMethod("crps")
}


#' @rdname crps
#' @export
scrps <- function(x1, ...) {
  UseMethod("scrps")
}

#' @rdname crps
#' @export
#' @param ll A log-likelihood.default or array (for `loo_*`-functions)
#' @param r_eff A vector of relative effective sample size estimates containing
#'     one element per observation (for `loo_*`-functions).
loo_crps <- function(x1, ...) {
  UseMethod("loo_crps")
}

#' @rdname crps
#' @export
#' @param ll A log-likelihood.default or array (for `loo_`-functions)
#' @param r_eff A vector of relative effective sample size estimates containing
#'     one element per observation (for `loo_`-functions).
loo_scrps <- function(x1, ...) {
  UseMethod("loo_scrps")
}


#' @rdname crps
#' @export
crps.default <- function(x1, x2, y) {
  validate_crps_input(x1, x2, y)
  EXX <- colMeans(abs(x1 - x2))
  EXy <- colMeans(abs(sweep(x1, 2, y)))
  return(.crps_fun(EXX, EXy))
}


#' @rdname crps
#' @export
loo_crps.default <- function(x1, x2, y, ll, r_eff = NULL, ...) {
  validate_crps_input(x1, x2, y, ll)
  psis_obj <- psis(-ll, r_eff = r_eff)
  EXX <- E_loo(abs(x1 - x2), psis_obj, log_ratios = -ll, ...)$value
  EXy <- E_loo(abs(sweep(x1, 2, y)), psis_obj, log_ratios = -ll, ...)$value
  return(.crps_fun(EXX, EXy))
}


#' @rdname crps
#' @export
scrps.default <- function(x1, x2, y) {
  validate_crps_input(x1, x2, y)
  EXy <- colMeans(abs(sweep(x1, 2, y)))
  EXX <- colMeans(abs(x1 - x2))
  return(.crps_fun(EXX, EXy, scale = TRUE))
}

#' @rdname crps
#' @export
loo_scrps.default <- function(x1, x2, y, ll, r_eff = NULL, ...) {
  validate_crps_input(x1, x2, y, ll)
  psis_obj <- psis(-ll, r_eff = r_eff)

  EXy <- E_loo(abs(sweep(x1, 2, y)), psis_obj, log_ratios = -ll, ...)$value
  EXX <- E_loo(abs(x1 - x2), psis_obj, log_ratios = -ll, ...)$value
  return(.crps_fun(EXX, EXy, scale = TRUE))
}


# ------------ Internals ----------------

#' Function to compute crps and scrps
#' @noRd
.crps_fun <- function(EXX, EXy, scale = FALSE) {
  if(scale) return(-EXy/EXX - 0.5 * log(EXX))
  return(0.5 * EXX - EXy)
}


#' Validate input of CRPS functions
#'
#' Check that predictive draws and observed data are of compatible shape
#' @noRd
validate_crps_input <- function(x1, x2, y, ll = NULL) {
  stopifnot(is.numeric(x1),
            is.numeric(x2),
            is.numeric(y),
            identical(dim(x1), dim(x2)),
            ncol(x1) == length(y),
            ifelse(is.null(ll), TRUE, identical(dim(ll), dim(x1)))
            )
}
