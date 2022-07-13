#' Continuously ranked probability score
#'
#' The `crps()` and `scrps()` and their `loo` counterparts can be used to
#' compute the continuously ranked probability score (CRPS) and scaled CRPS
#' (SCRPS) (see Bolin and Wallin, 2022). CRPS is a proper scoring rule, and
#' strictly proper when the first moment of the predictive distribution is
#' finite. Both can be expressed in terms of samples form the predictive
#' distribution. See e.g. Gneiting and Raftery (2007) for a comprehensive
#' discussion on CRPS.
#'
#' To compute (S)CRPS, user needs to provide two sets of draws from the
#' predictive distribution. I.e. one needs to call, for instance,
#' `posterior_predict()` two times and pass the returned values on to the CRPS
#' functions. This is due to the fact that formulas used to compute CRPS involve
#' an expectation of the absolute difference of `X1` and `X2`, both having the
#' same distribution. See Gneiting and Raftery (2007) for details.
#'
#' @export
#' @param x1 A S by N matrix, or a vector of length S when only single
#'     observation is provided.
#' @param x2 Independent draws from the same distribution as draws in `x1`.
#'     Should be of the identical dimension.
#' @param y A vector of observations or a single value.
#' @param permutations An integer, with default value of 1,  specifying how many
#'     times the expected value of  |X - X'| is computed. The row order of `x2`
#'     is shuffled to as elements `x1` and `x2` are typically drawn given the
#'     same value of parameters. This happens e.g. when one calls
#'     `posterior_predict()` twice for a fitted `rstanarm` model. Generating
#'     more permutations is expected to decrease the variance of the computed
#'     expected value.
#' @param ... Passed on to `E_loo()` in `loo_*`-functions.
#'
#' @return A list containing two elements: `estimates` and `pointwise`.
#'   The former reports estimator and standard error and latter the pointwise
#'   values.
#'
#' @examples
#' \dontrun{
#' # An example using rstanarm
#' library(rstanarm)
#' data("kidiq")
#' fit <- stan_glm(kid_score ~ mom_hs + mom_iq, data = kidiq)
#' ypred1 <- posterior_predict(fit)
#' ypred2 <- posterior_predict(fit)
#' y <- get_y(fit)
#' crps(ypred1, ypred2, y)
#' loo_crps(ypred1, ypred2, y, ll = log_lik(fit))
#' }
#'
#' @references
#' #' Bolin, D., & Wallin, J. (2022). Local scale invariance and robustness of
#' proper scoring rules. arXiv.
#' <https://doi.org/10.48550/arXiv.1912.05642>
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules,
#' Prediction, and Estimation. Journal of the American Statistical Association,
#' 102(477), 359–378.
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
crps.matrix <- function(x1, x2, y, permutations = 1) {
  validate_crps_input(x1, x2, y)
  repeats <- replicate(permutations, EXX_compute(x1, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x1, 2, y)))
  return(crps_output(.crps_fun(EXX, EXy)))
}


#' Method for a single data point
#' @rdname crps
#' @export
crps.numeric <- function(x1, x2, y, permutations = 1) {
  stopifnot(length(x1) == length(x2),
            length(y) == 1)
  crps.matrix(as.matrix(x1), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
loo_crps.matrix <-
  function(x1,
           x2,
           y,
           ll,
           permutations = 1,
           r_eff = NULL,
           ...) {
  validate_crps_input(x1, x2, y, ll)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x1, x2, ll, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-ll, r_eff = r_eff)
  EXy <- E_loo(abs(sweep(x1, 2, y)), psis_obj, log_ratios = -ll, ...)$value
  return(crps_output(.crps_fun(EXX, EXy)))
}


#' @rdname crps
#' @export
scrps.matrix <- function(x1, x2, y, permutations = 1) {
  validate_crps_input(x1, x2, y)
  repeats <- replicate(1, EXX_compute(x1, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x1, 2, y)))
  return(crps_output(.crps_fun(EXX, EXy, scale = TRUE)))
}

#' @rdname crps
#' @export
scrps.numeric <- function(x1, x2, y, permutations = 1) {
  stopifnot(length(x1) == length(x2),
            length(y) == 1)
  scrps.matrix(as.matrix(x1), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
loo_scrps.matrix <-
  function(
    x1,
    x2,
    y,
    ll,
    permutations = 1,
    r_eff = NULL,
    ...) {
  validate_crps_input(x1, x2, y, ll)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x1, x2, ll, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-ll, r_eff = r_eff)
  EXy <- E_loo(abs(sweep(x1, 2, y)), psis_obj, log_ratios = -ll, ...)$value
  return(crps_output(.crps_fun(EXX, EXy, scale = TRUE)))
}

# ------------ Internals ----------------


EXX_compute <- function(x1, x2) {
  S <- nrow(x1)
  EXX <- colMeans(abs(x1 - x2[sample(1:S),]))
  return(EXX)
}


EXX_loo_compute <- function(x1, x2, ll, r_eff = NULL, ...) {
  S <- nrow(x1)
  shuffle <- sample (1:S)
  x2 <- x2[shuffle,]
  ll2 <- ll[shuffle,]
  psis_obj_joint <- psis(-ll - ll2 , r_eff = r_eff)
  EXX <- E_loo(abs(x1 - x2), psis_obj_joint, log_ratios = -ll - ll2, ...)$value
  return(EXX)
}


#' Function to compute crps and scrps
#' @noRd
.crps_fun <- function(EXX, EXy, scale = FALSE) {
  if(scale) return(-EXy/EXX - 0.5 * log(EXX))
  return( 0.5 * EXX - EXy)
}

#' Compute output data for crps functions
#' @noRd
crps_output <- function(crps_pw) {
  n <- length(crps_pw)
  out <- list()
  out$estimates <- c(mean(crps_pw), sd(crps_pw) / sqrt(n))
  names(out$estimates) <- c('Estimate', 'SE')
  out$pointwise <- crps_pw
  return(out)
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
