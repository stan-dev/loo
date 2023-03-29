#' Continuously ranked probability score
#'
#' The `crps()` and `scrps()` functions and their `loo_*()` counterparts can be
#' used to compute the continuously ranked probability score (CRPS) and scaled
#' CRPS (SCRPS) (see Bolin and Wallin, 2022). CRPS is a proper scoring rule, and
#' strictly proper when the first moment of the predictive distribution is
#' finite. Both can be expressed in terms of samples form the predictive
#' distribution. See e.g. Gneiting and Raftery (2007) for a comprehensive
#' discussion on CRPS.
#'
#' To compute (S)CRPS, the user needs to provide two sets of draws, `x` and
#' `x2`, from the predictive distribution. This is due to the fact that formulas
#' used to compute CRPS involve an expectation of the absolute difference of `x`
#' and `x2`, both having the same distribution. See the `permutations` argument,
#' as well as Gneiting and Raftery (2007) for details.
#'
#' @export
#' @param x A `S` by `N` matrix (draws by observations), or a vector of length
#'   `S` when only single observation is provided in `y`.
#' @param x2 Independent draws from the same distribution as draws in `x`.
#'   Should be of the identical dimension.
#' @param y A vector of observations or a single value.
#' @param permutations An integer, with default value of 1,  specifying how many
#'   times the expected value of  |X - X'| (`|x - x2|`) is computed. The row
#'   order of `x2` is shuffled as elements `x` and `x2` are typically drawn
#'   given the same values of parameters. This happens, e.g., when one calls
#'   `posterior_predict()` twice for a fitted \pkg{rstanarm} or \pkg{brms}
#'   model. Generating more permutations is expected to decrease the variance of
#'   the computed expected value.
#' @param ... Passed on to [E_loo()] in the `loo_*()` version of these
#'   functions.
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
#' crps(ypred1, ypred2, y = fit$y)
#' loo_crps(ypred1, ypred2, y = fit$y, log_lik = log_lik(fit))
#' }
#'
#' @references
#' Bolin, D., & Wallin, J. (2022). Local scale invariance and robustness of
#' proper scoring rules. arXiv. \doi{10.48550/arXiv.1912.05642}
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules,
#' Prediction, and Estimation. Journal of the American Statistical Association,
#' 102(477), 359–378.
crps <- function(x, ...) {
  UseMethod("crps")
}


#' @rdname crps
#' @export
scrps <- function(x, ...) {
  UseMethod("scrps")
}

#' @rdname crps
#' @export
loo_crps <- function(x, ...) {
  UseMethod("loo_crps")
}

#' @rdname crps
#' @export
loo_scrps <- function(x, ...) {
  UseMethod("loo_scrps")
}


#' @rdname crps
#' @export
crps.matrix <- function(x, x2, y, ..., permutations = 1) {
  validate_crps_input(x, x2, y)
  repeats <- replicate(permutations, EXX_compute(x, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x, 2, y)))
  crps_output(.crps_fun(EXX, EXy))
}


#' Method for a single data point
#' @rdname crps
#' @export
crps.numeric <- function(x, x2, y, ..., permutations = 1) {
  stopifnot(length(x) == length(x2),
            length(y) == 1)
  crps.matrix(as.matrix(x), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
#' @param log_lik A log-likelihood matrix the same size as `x`.
#' @param r_eff An optional vector of relative effective sample size estimates
#'   containing one element per observation. See [psis()] for details.
#' @param cores The number of cores to use for parallelization of `[psis()]`.
#'   See [psis()] for details.
loo_crps.matrix <-
  function(x,
           x2,
           y,
           log_lik,
           ...,
           permutations = 1,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
  validate_crps_input(x, x2, y, log_lik)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x, x2, log_lik, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-log_lik, r_eff = r_eff, cores = cores)
  EXy <- E_loo(abs(sweep(x, 2, y)), psis_obj, log_ratios = -log_lik, ...)$value
  crps_output(.crps_fun(EXX, EXy))
}


#' @rdname crps
#' @export
scrps.matrix <- function(x, x2, y, ..., permutations = 1) {
  validate_crps_input(x, x2, y)
  repeats <- replicate(permutations, EXX_compute(x, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x, 2, y)))
  crps_output(.crps_fun(EXX, EXy, scale = TRUE))
}

#' @rdname crps
#' @export
scrps.numeric <- function(x, x2, y, ..., permutations = 1) {
  stopifnot(length(x) == length(x2),
            length(y) == 1)
  scrps.matrix(as.matrix(x), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
loo_scrps.matrix <-
  function(
    x,
    x2,
    y,
    log_lik,
    ...,
    permutations = 1,
    r_eff = NULL,
    cores = getOption("mc.cores", 1)) {
  validate_crps_input(x, x2, y, log_lik)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x, x2, log_lik, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-log_lik, r_eff = r_eff, cores = cores)
  EXy <- E_loo(abs(sweep(x, 2, y)), psis_obj, log_ratios = -log_lik, ...)$value
  crps_output(.crps_fun(EXX, EXy, scale = TRUE))
}

# ------------ Internals ----------------


EXX_compute <- function(x, x2) {
  S <- nrow(x)
  colMeans(abs(x - x2[sample(1:S),]))
}


EXX_loo_compute <- function(x, x2, log_lik, r_eff = NULL, ...) {
  S <- nrow(x)
  shuffle <- sample (1:S)
  x2 <- x2[shuffle,]
  log_lik2 <- log_lik[shuffle,]
  psis_obj_joint <- psis(-log_lik - log_lik2 , r_eff = r_eff)
  E_loo(abs(x - x2), psis_obj_joint, log_ratios = -log_lik - log_lik2, ...)$value
}


#' Function to compute crps and scrps
#' @noRd
.crps_fun <- function(EXX, EXy, scale = FALSE) {
  if (scale) return(-EXy/EXX - 0.5 * log(EXX))
  0.5 * EXX - EXy
}

#' Compute output data for crps functions
#' @noRd
crps_output <- function(crps_pw) {
  n <- length(crps_pw)
  out <- list()
  out$estimates <- c(mean(crps_pw), sd(crps_pw) / sqrt(n))
  names(out$estimates) <- c('Estimate', 'SE')
  out$pointwise <- crps_pw
  out
}

#' Validate input of CRPS functions
#'
#' Check that predictive draws and observed data are of compatible shape
#' @noRd
validate_crps_input <- function(x, x2, y, log_lik = NULL) {
  stopifnot(is.numeric(x),
            is.numeric(x2),
            is.numeric(y),
            identical(dim(x), dim(x2)),
            ncol(x) == length(y),
            ifelse(is.null(log_lik), TRUE, identical(dim(log_lik), dim(x)))
            )
}
