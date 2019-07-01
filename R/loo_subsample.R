#' @title
#' Efficient approximate leave-one-out cross-validation (LOO) for posterior approximations
#'
#'
#' @param x A log-likelihood array, matrix, or function. The **Methods (by class)**
#'   section, below, has detailed descriptions of how to specify the inputs for
#'   each method.
#' @export loo_subsample loo_subsample.function
#'
#' @param save_psis Should the `"psis"` object created internally by `loo_approximate_posterior()` be
#'   saved in the returned object? See \link{loo} for details.
#' @template cores
#'
#' @details The `loo_approximate_posterior()` function is an S3 generic and methods are provided for
#'   3-D pointwise log-likelihood arrays, pointwise log-likelihood matrices, and
#'   log-likelihood functions. The implementation work for posterior approximations
#'   where it is possible to compute the log density for the posterior approximation.
#'
#' @return The `loo_approximate_posterior()` methods return a named list with class
#'   `c("psis_loo_ap", "psis_loo", "loo")` with the additional slot:
#' \describe{
#'  \item{`posterior_approximation`}{
#'   A list with two vectors, `log_p` and `log_g` of the same length
#'   containing the posterior density and the approximation density
#'   for the individual dras.
#'  }
#' }
#'
#' @seealso loo, psis, loo_compare
#'
#' @template loo-large-data-references
#'
#' @export
loo_subsample <- function(x, ...) {
  if (!requireNamespace("checkmate", quietly=TRUE)) {
    stop("Please install the 'checkmate' package to use this function.", call. = FALSE)
  }
  UseMethod("loo_subsample")
}

#' @export
#' @templateVar fn loo_subsample
#' @template function
#' @param data,draws,... For the `loo_subsample.function()` method and the `loo_i()`
#'   function, these are the data, posterior draws, and other arguments to pass
#'   to the log-likelihood function. See the **Methods (by class)** section
#'   below for details on how to specify these arguments.
#' @param observations The subsample observations to use. The argument can take four (4) types of arguments.
#'   If \code{NULL} is supplied, all observations are used and the algorithm just use standard psis loo.
#'   If a single integer is supplied this is the number of integers used.
#'   If a vector of integers, this will be the indecies used to subset the data. Note, this should only be used
#'   when subsampling have been done previously.
#'   If a \code{psis_loo_ss} object is supplied, the observations used for that is used.
#'
#' @param log_p Should be supplied if approximate posterior draws are used. Default (NULL) posterior draws from true posterior (i.e. using MCMC). The log-posterior (target) evaluated at S samples from the proposal distribution (g). A vector of length S.
#' @param log_q Should be supplied if approximate posterior draws are used. Default (NULL) posterior draws from true posterior (i.e. using MCMC). The log-density (proposal) evaluated at S samples from the proposal distribution (g). A vector of length S.
#'
#' @param loo_approximation What type of approximation the PSIS-LOO:s should be used.
#'   Default is \code{plpd} or log predictive density using the posterior expectation.
#'   If \code{NULL}, no approximation is done.
#'   There are six different methods implemented to approximate loo_i:s.
#'   See references for more details.
#'   \describe{
#'     \item{\code{plpd}}{
#'     uses the lpd based on point estimates (ie. \eqn{p(y_i|\hat{\theta})})}
#'     \item{\code{lpd}}{
#'     uses the lpds (ie. \eqn{p(y_i|y)})}
#'     \item{\code{tis}}{
#'     uses truncated importance sampling to approximate PSIS-LOO}
#'     \item{\code{waic}}{
#'     uses waic (ie. \eqn{p(y_i|y) - p_{waic}})}
#'     \item{\code{waic_grad_marginal}}{
#'     uses waic approximation using first order delta method and
#'     posterior marginal variances to approximate \eqn{p_{waic}}
#'     (ie. \eqn{p(y_i|\hat{\theta})}-p_waic_grad_marginal).
#'     Require gradient of likelihood function.}
#'     \item{\code{waic_grad}}{
#'     uses waic approximation using first order delta method and
#'     posterior covariance to approximate \eqn{p_{waic}}
#'     (ie. \eqn{p(y_i|\hat{\theta})}-p_waic_grad).
#'     Require gradient of likelihood function.}
#'     \item{\code{waic_grad}}{
#'     uses waic approximation using second order delta method and
#'     posterior covariance to approximate \eqn{p_{waic}}
#'     (ie. \eqn{p(y_i|\hat{\theta})}-p_waic_grad).
#'     Require gradient and Hessian of likelihood function.}
#'  }
#'  As points estimates of \eqn{\hat{\theta}}, the expectation
#'  of the posterior parameters  ares used.
#'
#' @param loo_approximation_draws The number of posterior draws
#' used when integrating over the posterior.
#' Used by \code{lpd} and \code{waic}.
#'
#' @param estimator How should elpd_loo be estimated. Default is \code{diff_srs}.
#'   \describe{
#'     \item{\code{diff_srs}}{
#'     uses the difference estimator with simple random sampling (srs)}
#'     \item{\code{hh}}{
#'     uses the Hansen-Hurwitz estimator with sampling proportional
#'     to size, where abs of loo_approximation is used as size.}
#'     \item{\code{srs}}{
#'     uses simple random sampling and ordinary estimation.}
#'  }
#'
#' @param llgrad the gradient of the log-likelihood. This is used
#'        with Default is \code{NULL}.
#' @param llhess the hessian of the log-likelihood. Default is \code{NULL}.
#' @template cores
#'
loo_subsample.function <-
  function(x,
           ...,
           data = NULL,
           draws = NULL,
           observations = 400,
           log_p = NULL,
           log_g = NULL,
           r_eff = NULL,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           loo_approximation = "plpd",
           loo_approximation_draws = NULL,
           estimator = "diff_srs",
           llgrad = NULL,
           llhess = NULL) {

    # Asserting inputs
    .llfun <- validate_llfun(x)
    stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
    if (checkmate::test_class(observations, "psis_loo_ss")) {
      # TODO: Fix this
      stop("Check also that the estimator of the object is the same")
      observations <- obs_idx(observations)
    }
    observations <- checkmate::assert_integerish(observations, null.ok = TRUE,
                                                 lower = 1,
                                                 min.len = 1,
                                                 any.missing = FALSE,
                                                 coerce = TRUE)
    if(length(observations) > 1) checkmate::assert_integerish(observations, upper = dim(data)[1])

    checkmate::assert_numeric(log_p, len = length(log_g), null.ok = TRUE)
    checkmate::assert_null(dim(log_p))
    checkmate::assert_numeric(log_g, len = length(log_p), null.ok = TRUE)
    checkmate::assert_null(dim(log_g))

    if(is.null(log_p) & is.null(log_g)){
      if (is.null(r_eff)) {
        throw_loo_r_eff_warning()
      } else {
        r_eff <- prepare_psis_r_eff(r_eff, len = N)
      }
    }
    checkmate::assert_flag(save_psis)
    cores <- loo_cores(cores)

    checkmate::assert_choice(loo_approximation, choices = loo_approximation_choices(), null.ok = FALSE)
    checkmate::assert_int(loo_approximation_draws, lower = 2, null.ok = TRUE)
    checkmate::assert_choice(estimator, choices = estimator_choices())

    .llgrad <- .llhess <- NULL
    if(!is.null(llgrad)) .llgrad <- validate_llfun(llgrad)
    if(!is.null(llhess)) .llhess <- validate_llfun(llhess)

    # Fallbacks
    if(is.null(observations)){
      if(is.null(log_p) & is.null(log_g)){
        lobj <- loo.function(x, ...,
                             data = data,
                             draws = draws,
                             r_eff = r_eff,
                             save_psis = save_psis,
                             cores = cores)
      } else {
        lobj <- loo_approximate_posterior.function(
                             x, ...,
                             log_p = log_p,
                             log_g = log_g,
                             data = data,
                             draws = draws,
                             r_eff = r_eff,
                             save_psis = save_psis,
                             cores = cores)
      }
      return(lobj)
    }

    # Compute loo approximation
    elpd_loo_approx <-
      elpd_loo_approximation(.llfun = .llfun,
                             data = data, draws = draws,
                             cores = cores,
                             loo_approximation = loo_approximation,
                             loo_approximation_draws = loo_approximation_draws,
                             .llgrad = .llgrad, .llhess = .llhess)

    # Draw subsample of observations
    if(length(observations) == 1){
      # Compute idxs
      idxs <- subsample_idxs(estimator = estimator,
                             elpd_loo_approximation = elpd_loo_approx,
                             m = observations)
    } else {
      # Compute idxs
      idxs <- compute_idxs(observations)
    }
    data_subsample <- data[idxs$idx,, drop = FALSE]
    if(length(r_eff) > 1) r_eff <- r_eff[idxs$idx]

    # Compute elpd_loo
    # TODO: Add test with long r_eff, i.e. handling r_eff
    # TODO: Fix so r_eff is only computed for subsample
    if(!is.null(log_p) & !is.null(log_g)){
      plo <- loo_approximate_posterior.function(x = .llfun,
                                                data = data_subsample,
                                                draws = draws,
                                                log_p = log_p,
                                                log_q = log_q,
                                                r_eff = r_eff,
                                                save_psis = save_psis,
                                                cores = cores)
    } else {
      plo <- loo.function(x = .llfun,
                          data = data_subsample,
                          draws = draws,
                          r_eff = r_eff,
                          save_psis = save_psis,
                          cores = cores)
    }

    # Construct ss object and estimate
    loo_ss <- psis_loo_ss_object(x = plo,
                                 idxs = idxs,
                                 elpd_loo_approx = elpd_loo_approx,
                                 loo_approximation = loo_approximation,
                                 loo_approximation_draws = loo_approximation_draws,
                                 estimator = estimator,
                                 .llfun = .llfun,
                                 .llgrad = .llgrad,
                                 .llhess = .llhess)
    loo_ss
  }


#' Update \code{psis_loo_ss} objects
#'
#' @details
#' If \code{observation} is updated, two approaches can be used.
#' Either a vector of indecies or a \code{psis_loo_ss} object,
#' the the updated object will have exactly the same observations
#' as supplied.
#' If an integer is supplied, new observations will be sampled to
#' reach the supplied sample size.
#'
#' @export
update.psis_loo_ss <- function(object,
                               data,
                               draws,
                               observations = NULL,
                               r_eff = NULL,
                               cores = getOption("mc.cores", 1),
                               loo_approximation = NULL,
                               loo_approximation_draws = NULL,
                               llgrad = NULL,
                               llhess = NULL){
  # Fallback
  # TODO: if nothing is updated, just return the object

  stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
  cores <- loo_cores(cores)

  # TODO: Add r_eff in object, ie handle it. See relative_eff() function.


  # Update elpd approximations
  if(!is.null(loo_approximation) | !is.null(loo_approximation_draws)){
    if(is.null(loo_approximation)) loo_approximation <- object$subsamling_loo$loo_approximation
    if(is.null(loo_approximation_draws)) loo_approximation_draws <- object$subsamling_loo$loo_approximation_draws
    if(is.null(llgrad)) .llgrad <- object$subsamling_loo$.llgrad else .llgrad <- validate_llfun(llgrad)
    if(is.null(llhess)) .llhess <- object$subsamling_loo$.llhess else .llhess <- validate_llfun(llhess)

    # TODO: Add tests for dimensions of data and params
    # Compute loo approximation
    elpd_loo_approx <-
      elpd_loo_approximation(.llfun = object$subsamling_loo$.llfun,
                             data = data, draws = draws,
                             cores = cores,
                             loo_approximation = loo_approximation,
                             loo_approximation_draws = loo_approximation_draws,
                             .llgrad = .llgrad, .llhess = .llhess)
    # Update object
    object$subsamling_loo$elpd_loo_approx <- elpd_loo_approx
    object$subsamling_loo$loo_approximation <- loo_approximation
    object$subsamling_loo$loo_approximation_draws <- loo_approximation_draws
    object$subsamling_loo$.llgrad <- .llgrad
    object$subsamling_loo$.llhess <- .llhess
    object$pointwise[, "elpd_loo_approx"] <- object$subsamling_loo$elpd_loo_approx[object$pointwise[, "idx"]]
  }

  # TODO: test update loo observations
  # Update observations
  if(!is.null(observations)){
    if (checkmate::test_class(observations, "psis_loo_ss")) {
      # TODO: Fix this
      stop("Check also that the estimator of the object is the same")
      observations <- obs_idx(observations)
    }
    observations <- checkmate::assert_integerish(observations, null.ok = TRUE,
                                                 lower = 1,
                                                 min.len = 1,
                                                 any.missing = FALSE,
                                                 coerce = TRUE)
    if(length(observations) > 1){
      checkmate::assert_integerish(observations, upper = dim(data)[1])
    } else {
      # Assert we add new observations
      checkmate::assert_int(observations, lower = nobs(object) + 1)
    }

    # TODO: Only update observations if a new vector is supplied or new obs are more than previous (i.e. adding)
    # TODO: Add in estimators that remove points that are m_i == 0

    # Compute subsample indecies
    if(length(observations) > 1){
      idxs <- compute_idxs(observations)
    } else {
      current_obs <- nobs(object)

      # If sampling with replacement
      if(object$subsamling_loo$estimator %in% c("hh")){
        idxs <- subsample_idxs(estimator = object$subsamling_loo$estimator,
                               elpd_loo_approximation = object$subsamling_loo$elpd_loo_approx,
                               m = observations - current_obs)
      }
      # If sampling without replacement
      # TODO: test this part
      if(object$subsamling_loo$estimator %in% c("diff_srs", "srs")){
        current_idxs <- obs_idx(object, rep = FALSE)
        new_idx <- (1:length(object$subsamling_loo$elpd_loo_approx))[-current_idxs]
        idxs <- subsample_idxs(estimator = object$subsamling_loo$estimator,
                               elpd_loo_approximation = object$subsamling_loo$elpd_loo_approx[-current_idxs],
                               m = observations - current_obs)
        idxs$idx <- new_idx[idxs$idx]
      }
    }

    # Identify how to update object
    cidxs <- compare_idxs(idxs, object)
    data_new_subsample <- data[cidxs$new$idx,, drop = FALSE]
    if(length(r_eff) > 1) r_eff <- r_eff[cidxs$new$idx]

    # TODO: Assert that we cant update loo_approx with HH (since need sampling)
    # TODO: Create function for set of estimators with and withour replacement

    # Compute new observations
    # TODO: Add test with long r_eff, i.e. handling r_eff
    # TODO: Fix so r_eff is only computed for subsample

    if(!is.null(cidxs$new)){
      if(!is.null(object$approximate_posterior$log_p) & !is.null(object$approximate_posterior$log_g)){
        plo <- loo_approximate_posterior.function(x = .llfun,
                                                  data = data_new_subsample,
                                                  draws = draws,
                                                  log_p = object$approximate_posterior$log_p,
                                                  log_q = object$approximate_posterior$log_g,
                                                  r_eff = r_eff,
                                                  save_psis = !is.null(object$psis_object),
                                                  cores = cores)
      } else {
        plo <- loo.function(x = .llfun,
                            data = data_new_subsample,
                            draws = draws,
                            r_eff = r_eff,
                            save_psis = !is.null(object$psis_object),
                            cores = cores)
      }
      # Add stuff to pointwise
      plo$pointwise <-
        add_subsampling_vars_to_pointwise(plo$pointwise,
                                          cidxs$new,
                                          object$subsamling_loo$elpd_loo_approx)
    }


    # TODO: Combine psis objects from psis_save in rbind.psis_loo_ss
    # TODO: Test so all are null below and it works as expected

    # Update object (diagnostics and pointwise)
    if(is.null(cidxs$new)) plo <- NULL
    if(length(observations) == 1){
      # Add new samples pointwise and diagnostic
      object<- rbind.psis_loo_ss(object, plo)

      # Update m_i for current pointwise (diagnostic stay the same)
      object$pointwise <- update_m_i_in_pointwise(pointwise, cidxs$add, type = "add")
    } else {
      # Add new samples pointwise and diagnostic
      object<- rbind.psis_loo_ss(object, plo)

      # Replace m_i current pointwise and diagnostics
      object$pointwise <- update_m_i_in_pointwise(pointwise, cidxs$add, type = "replace")

      # Remove samples
      object <- remove_idx.psis_loo_ss(object, cidxs$remove)
    }
  }

  # Compute estimates
  if(object$subsamling_loo$estimator == "hh"){
    object <- loo_subsample_estimation_hh(object)
  } else if(object$subsamling_loo$estimator == "diff_srs"){
    object <- loo_subsample_estimation_diff_srs(object)
  } else if(object$subsamling_loo$estimator == "srs"){
    object <- loo_subsample_estimation_srs(object)
  } else {
    stop("No correct estimator used.")
  }
  # TODO: Fix assert_psis_loo_ss
  assert_psis_loo_ss(object)
  object
}

# TODO: Test nobs() and obs_idx()

#' The observations indecies in the original data
#'
#' @details
#' Returns a vector of the same length as the number of observations
#' (if \code{rep = TRUE})
#' If sampling with replacement is used, an observation can have multiples
#' samples, these are then repeated (ie a vector c(1,1,2) indicate that
#' observation 1 has been ssubampled two times.
#' If \code{rep = FALSE} only the unique indecies are returned.
#'
#' @param x a \code{psis_loo_ss} object.
#' @param rep should multiple samples of same observation be returned?
#'
#' @return an integer lector of length \code{nobs}.
#'
#' @export
obs_idx <- function(x, rep = TRUE){
  checkmate::assert_class(x, "psis_loo_ss")
  if(rep){
    idxs <- as.integer(rep(x$pointwise[,"idx"], x$pointwise[,"m_i"]))
  } else {
    idxs <- as.integer(x$pointwise[,"idx"])
  }
  idxs
}

# TODO: Test that nobs and obs_idx give the same result for HH

#' @export
nobs.psis_loo_ss <- function(x){
  as.integer(sum(x$pointwise[,"m_i"]))
}


# internal ----------------------------------------------------------------

loo_approximation_choices <- function() {
  c("plpd", "lpd", "waic", "waic_grad_marginal", "waic_grad", "waic_hess", "none")
}

estimator_choices <- function() {
  c("hh", "diff_srs", "srs")
}


## Approximate elpd -----

#' Compute approximation to loo_i:s
#'
#' @details
#' See \link{\code{loo_subsample.function()}} \code{loo_approximation} argument.
#'
#' @inheritParams loo_subsample.function
#'
#' @return a vector with approximations of elpd_{loo,i}s
#'
elpd_loo_approximation <- function(.llfun, data, draws, cores, loo_approximation, loo_approximation_draws = NULL, .llgrad = NULL, .llhess = NULL){
  checkmate::assert_function(.llfun, args = c("data_i", "draws"), ordered = TRUE)
  stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
  checkmate::assert_choice(loo_approximation, choices = loo_approximation_choices(), null.ok = TRUE)
  checkmate::assert_int(loo_approximation_draws, lower = 2, null.ok = TRUE)
  cores <- loo_cores(cores)
  if(!is.null(.llgrad)){
    checkmate::assert_function(.llgrad, args = c("data_i", "draws"), ordered = TRUE)
  }
  if(!is.null(.llhess)){
    checkmate::assert_function(.llhess, args = c("data_i", "draws"), ordered = TRUE)
  }

  N <- dim(data)[1]
  # TODO: Test this in testsuite
  if(is.null(loo_approximation)) return(NULL)

  if(loo_approximation == "waic"){
    # TODO: stop("make this more efficient, by extracting")
    draws <- thin_draws(draws, loo_approximation_draws)
    waic_full_obj <- waic.function(.llfun, data = data, draws = draws)
    return(waic_full_obj$pointwise[,"elpd_waic"])
  }

  # Compute the lpd or log p(y_i|y_{-i})
  if (loo_approximation == "lpd"){
    draws <- thin_draws(draws, loo_approximation_draws)
    lpds <- unlist(lapply(seq_len(N), FUN = function(i) {
      ll_i <- .llfun(data_i = data[i,, drop=FALSE], draws = draws)
      ll_i <- as.vector(ll_i)
      lpd_i <- logMeanExp(ll_i)
      lpd_i
    }))
    return(lpds) # Use only the lpd
  }

  # Compute the point lpd or log p(y_i|\hat{\theta}) - also used in waic_delta approaches
  if(loo_approximation == "plpd" |
     loo_approximation == "waic_grad" |
     loo_approximation == "waic_grad_marginal" |
     loo_approximation == "waic_hess"){

    point_est <- compute_point_estimate(draws)
    lpds <- unlist(lapply(seq_len(N), FUN = function(i) {
      ll_i <- .llfun(data_i = data[i,, drop=FALSE], draws = point_est)
      ll_i <- as.vector(ll_i)
      lpd_i <- logMeanExp(ll_i)
      lpd_i
    }))
    if(loo_approximation == "plpd") return(lpds) # Use only the lpd
  }

  if(loo_approximation == "waic_grad" |
     loo_approximation == "waic_grad_marginal" |
     loo_approximation == "waic_hess") {
    checkmate::assert_true(!is.null(.llgrad))

    point_est <- compute_point_estimate(draws)
    # Compute the lpds
    lpds <- unlist(lapply(seq_len(N), FUN = function(i) {
      ll_i <- .llfun(data_i = data[i,, drop=FALSE], draws = point_est)
      ll_i <- as.vector(ll_i)
      lpd_i <- logMeanExp(ll_i)
      lpd_i
    }))
    # TODO: Cleanup this function (when test suits pass)
    # Refactor each method to its own function instead of the if-else

    if(loo_approximation == "waic_grad" |
       loo_approximation == "waic_hess") {
      cov_est <- cov(draws)
    }

    if(loo_approximation == "waic_grad_marginal") {
      marg_vars <- apply(draws, MARGIN = 2, var)
    }

    p_eff_approx <- numeric(N)
    if(cores>1) warning("Multicore is not implemented for waic_delta") # TODO: Look at this

    if(loo_approximation == "waic_grad"){
      for(i in 1:nrow(data)){
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        local_cov <- cov_est[rownames(grad_i), rownames(grad_i)]
        p_eff_approx[i] <-  t(grad_i) %*% local_cov %*% grad_i
      }
    } else if(loo_approximation == "waic_grad_marginal") {
      for(i in 1:nrow(data)){
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        p_eff_approx[i] <- sum(grad_i * marg_vars[rownames(grad_i)] * grad_i)
      }
    } else if(loo_approximation == "waic_hess") {
      checkmate::assert_true(!is.null(.llhess))
      for(i in 1:nrow(data)){
        # TODO: Check with Michael on efficient implementation
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        hess_i <- .llhess(data_i = data[i,,drop = FALSE], draws = point_est[,rownames(grad_i), drop = FALSE])[,,1]
        local_cov <- cov_est[rownames(grad_i), rownames(grad_i)]
        p_eff_approx[i] <- t(grad_i) %*% local_cov %*% grad_i +
          0.5 * sum(diag(local_cov %*% hess_i %*% local_cov %*% hess_i))
      }
    } else {stop(loo_approximation, " is not implemented!")}

    return(lpds - p_eff_approx)
  }

}

# TODO: Implement TIS

### TODO: MAKE THIS GENERIC AND REMOVE STANREG
#' Compute \eqn{E(\theta)} as point estimate.
#' @rdname thin_draws
compute_point_estimate <- function(draws){
  if(is.matrix(draws)){
    draws <- t(as.matrix(colMeans(draws)))
  } else if(is.stanreg.draws(draws)){
    for(i in seq_along(draws)){
      if(is.matrix(draws[[i]])){
        draws[[i]] <- t(as.matrix(colMeans(draws[[i]])))
      }
    }
  }
  draws
}

### TODO: MAKE THIS GENERIC AND REMOVE STANREG
#' Thin draws for a matrix with draws and for rstanarm stan_reg objects.
#' @param draws a draws object with draws from the posterior.
#' @param loo_approximation_draws the number of posterior draws to return (ie after thinning)
thin_draws <- function(draws, loo_approximation_draws){
  if(!is.null(loo_approximation_draws)) {
    if(is.matrix(draws)){
      draws <- thin_draws_matrix(draws, loo_approximation_draws)
    } else if(is.stanreg.draws(draws)){
      for(i in seq_along(draws)){
        if(is.matrix(draws[[i]])){
          draws[[i]] <- thin_draws_matrix(draws[[i]], loo_approximation_draws)
        }
      }
    }
  }
  return(draws)
}
#' @rdname thin_draws
thin_draws_matrix <- function(draws, loo_approximation_draws){
  S <- nrow(draws)
  idx <- 1:loo_approximation_draws * S %/% loo_approximation_draws
  draws <- draws[idx, , drop = FALSE]
  draws
}



## Subsampling -----

#' Subsampling strategy used
#'
#' @param estimator The estimator to use (\code{hh},\code{srs_diff})
#'
subsample_idxs <- function(estimator, elpd_loo_approximation, m){
  checkmate::assert_choice(estimator, choices = estimator_choices())
  checkmate::assert_numeric(elpd_loo_approximation)
  checkmate::assert_int(m, lower = 1, upper = length(elpd_loo_approximation))

  if(estimator == "hh"){
    pi_values <- hh_elpd_approx_to_pis(elpd_loo_approximation)
    idxs_df <- pps_sample(m, pis = pi_values)
  }

  if(estimator == "diff_srs" | estimator == "srs"){
    idx <- 1:length(elpd_loo_approximation)
    idx_m <- idx[order(runif(length(elpd_loo_approximation)))][1:m]
    idx_m <- idx_m[order(idx_m)]
    idxs_df <- data.frame(idx=as.integer(idx_m), m_i=1L)
  }
  assert_subsample_idxs(x = idxs_df)
  idxs_df
}

#' Compute pis from approximation for use in HH estimation
hh_elpd_loo_approximation_to_pis <- function(elpd_loo_approximation){
  checkmate::assert_numeric(elpd_loo_approximation)
  pi_values <- abs(elpd_loo_approximation)
  pi_values <- pi_values/sum(pi_values) # \tilde{\pi}
  pi_values
}


#' Compute subsampling indecies from an observation vector
compute_idxs <- function(observations){
  checkmate::assert_integer(observations, lower = 1, min.len = 2, any.missing = FALSE)
  tab <- table(observations)
  idxs_df <- data.frame(idx = as.integer(names(tab)), m_i = as.integer(unname(tab)))
  assert_subsample_idxs(idxs_df)
  idxs_df
}


compare_idxs  <- function(idxs, object){
  assert_subsample_idxs(idxs)
  current_idx <- compute_idxs(obs_idx(object))
  result <- list()
  new_idx <- !(idxs$idx %in% current_idx$idx)
  remove_idx <- !(current_idx$idx %in% idxs$idx)

  result$new <- idxs[new_idx, ]
  if(nrow(result$new) == 0) {
    result["new"] <- NULL
  } else {
    assert_subsample_idxs(result$new)
  }

  result$add <- idxs[!new_idx, ]
  if(nrow(result$add) == 0) {
    result["add"] <- NULL
  } else {
    assert_subsample_idxs(result$add)
  }

  result$remove <- current_idx[remove_idx, ]
  if(nrow(result$remove) == 0) {
    result["remove"] <- NULL
  } else {
    assert_subsample_idxs(result$remove)
  }

  result
}

assert_subsample_idxs <- function(x){
  checkmate::assert_data_frame(x,
                               types = c("integer", "integer"),
                               any.missing = FALSE,
                               min.rows = 1,
                               col.names = "named")
  checkmate::assert_names(names(x), identical.to = c("idx", "m_i"))
  checkmate::assert_integer(x$idx, lower = 1, any.missing = FALSE, unique = TRUE)
  checkmate::assert_integer(x$m_i, lower = 1, any.missing = FALSE)
}

assert_psis_loo_ss <- function(x){
  # Test
  warning("implement these assertions")
  # Add test that the elpd_loo_approx in pointwise is correct
  x
}



#' Draw a sample and return a idx_df
#' @details
#' We are sampling with replacement, hence we only want to compute elpd
#' for each observation once.
pps_sample <- function(m, pis){
  checkmate::assert_int(m)
  checkmate::assert_numeric(pis, min.len = 2, lower = 0, upper = 1)
  idx <- sample(1:length(pis), size = m, replace = TRUE, prob = pis)
  idxs_df <- as.data.frame(table(idx), stringsAsFactors = FALSE)
  colnames(idxs_df) <- c("idx", "m_i")
  idxs_df$idx <- as.integer(idxs_df$idx)
  idxs_df$m_i <- as.integer(idxs_df$m_i)
  idxs_df
}

## Constructor ---

#' Add components needed for subsampling loo
#'
#' @param x either a psis_loo object or a quick_psis_loo object.
#' @param approx_corr has approximation corrections been done
psis_loo_ss_object <- function(x,
                               idxs,
                               elpd_loo_approx,
                               loo_approximation, loo_approximation_draws,
                               estimator,
                               .llfun, .llgrad, .llhess){
  # Assertions
  checkmate::assert_class(x, "psis_loo")
  assert_subsample_idxs(idxs)
  checkmate::assert_numeric(elpd_loo_approx, any.missing = FALSE)
  checkmate::assert_choice(loo_approximation, loo_approximation_choices())
  checkmate::assert_int(loo_approximation_draws, null.ok = TRUE)
  checkmate::assert_choice(estimator, estimator_choices())
  checkmate::assert_function(.llfun, args = c("data_i", "draws"), ordered = TRUE)
  checkmate::assert_function(.llgrad, args = c("data_i", "draws"), ordered = TRUE, null.ok = TRUE)
  checkmate::assert_function(.llhess, args = c("data_i", "draws"), ordered = TRUE, null.ok = TRUE)
  # TODO: Assert that the approx loo is the same in pointwise as the idx in the approx_loo vector

  # Construct object
  class(x) <- c("psis_loo_ss", class(x))
  x$pointwise <- add_subsampling_vars_to_pointwise(pointwise = x$pointwise, idxs, elpd_loo_approx)
  x$estimates <- cbind(x$estimates, matrix(0, nrow = nrow(x$estimates)))
  colnames(x$estimates)[ncol(x$estimates)] <- "subsampling SE"

  x$subsamling_loo <- list()
  x$subsamling_loo$elpd_loo_approx <- elpd_loo_approx
  x$subsamling_loo$loo_approximation <- loo_approximation
  x$subsamling_loo["loo_approximation_draws"] <- list(loo_approximation_draws)
  x$subsamling_loo$estimator <- estimator
  x$subsamling_loo$.llfun <- .llfun
  x$subsamling_loo[".llgrad"] <- list(.llgrad)
  x$subsamling_loo[".llhess"] <- list(.llhess)

  # Compute estimates
  if(estimator == "hh"){
    x <- loo_subsample_estimation_hh(x)
  } else if(estimator == "diff_srs"){
    x <- loo_subsample_estimation_diff_srs(x)
  } else if(estimator == "srs"){
    x <- loo_subsample_estimation_srs(x)
  } else {
    stop("No correct estimator used.")
  }
  assert_psis_loo_ss(x)
  x
}

add_subsampling_vars_to_pointwise <- function(pointwise, idxs, elpd_loo_approx){
  checkmate::assert_matrix(pointwise,
                           any.missing = FALSE,
                           min.cols = 4)
  checkmate::assert_names(colnames(pointwise), identical.to = c("elpd_loo","mcse_elpd_loo","p_loo","looic"))
  assert_subsample_idxs(idxs)
  checkmate::assert_numeric(elpd_loo_approx)

  pw <- cbind(as.data.frame(pointwise), idxs)
  pw$elpd_loo_approx <- elpd_loo_approx[idxs$idx]
  pw <- as.matrix(pw)
  rownames(pw) <- NULL
  assert_subsampling_pointwise(pw)
  pw
}

rbind.psis_loo_ss <- function(object, x){
  checkmate::assert_class(object, "psis_loo_ss")
  if(is.null(x)) return(object) # Fallback
  checkmate::assert_class(x, "loo_ss")
  assert_subsampling_pointwise(object$pointwise)
  assert_subsampling_pointwise(x$pointwise)
  checkmate::assert_disjunct(object$pointwise[, "idx"], x$pointwise[, "idx"])

  object$pointwise <- rbind(object$pointwise,
                            x$pointwise)
  object$diagnostics$pareto_k <- c(object$diagnostics$pareto_k,
                                   plo$diagnostics$pareto_k)
  object$diagnostics$n_eff <- c(object$diagnostics$n_eff,
                                plo$diagnostics$n_eff)
  object
}

remove_idx.psis_loo_ss <- function(object, idxs = cidxs$remove){
  checkmate::assert_class(object, "psis_loo_ss")
  if(is.null(idxs)) return(object) # Fallback
  assert_subsample_idxs(idxs)
  checkmate::assert_choice(type, choices = c("replace", "add"))

  row_map <- data.frame(row_no = 1:nrow(object$pointwise), idx = object$pointwise[, "idx"])
  row_map <- merge(row_map, idxs, by = "idx", all.y = TRUE)

  object$pointwise <- object$pointwise[-row_map$row_no,,drop = FALSE]
  object$diagnostics$pareto_k <- object$diagnostics$pareto_k[-row_map$row_no]
  object$diagnostics$n_eff <- object$diagnostics$n_eff[-row_map$row_no]
  object
}

update_m_i_in_pointwise <- function(pointwise, idxs, type = "replace"){
  assert_subsampling_pointwise(pointwise)
  if(is.null(idxs)) return(pointwise) # Fallback
  assert_subsample_idxs(idxs)
  checkmate::assert_choice(type, choices = c("replace", "add"))

  row_map <- data.frame(row_no = 1:nrow(pointwise), idx = pointwise[, "idx"])
  row_map <- merge(row_map, idxs, by = "idx", all.y = TRUE)

  if(type == "replace"){
    pointwise[row_map$row_no, "m_i"] <- row_map$m_i
  }
  if(type == "add"){
    pointwise[row_map$row_no, "m_i"] <- pointwise[row_map$row_no, "m_i"] + row_map$m_i
  }
  pointwise
}

remove_idx.psis_loo_ss <- function(pointwise, idxs = cidxs$add, type = "replace"){
  assert_subsampling_pointwise(pointwise)
  assert_subsample_idxs(idxs)
  checkmate::assert_choice(type, choices = c("replace", "add"))

  row_map <- data.frame(row_no = 1:nrow(pointwise), idx = pointwise[, "idx"])
  row_map <- merge(row_map, cidxs$add, by = "idx", all.y = TRUE)

  if(type == "replace"){
    pointwise[row_map$row_no, "m_i"] <- row_map$m_i
  }
  if(type == "add"){
    pointwise[row_map$row_no, "m_i"] <- pointwise[row_map$row_no, "m_i"] + row_map$m_i
  }
  pointwise
}



assert_subsampling_pointwise <- function(pointwise){
  checkmate::assert_matrix(pointwise,
                           any.missing = FALSE,
                           ncols = 7)
  checkmate::assert_names(colnames(pointwise), identical.to = c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic", "idx", "m_i", "elpd_loo_approx"))
}


## Estimation ---

#' Estimate elpd using the Hansen-Hurwitz estimator
#' @param x a quick_psis_loo object
loo_subsample_estimation_hh <- function(x){
  checkmate::assert_class(x, "psis_loo_ss")
  N <- length(x$subsamling_loo$elpd_loo_approx)
  pis <- hh_elpd_loo_approximation_to_pis(x$subsamling_loo$elpd_loo_approx)
  pis_sample <- pis[x$pointwise$idx]

  hh_elpd_loo <- whhest(z = pis_sample, m_i = x$pointwise$m_i, y = x$pointwise$elpd_loo, N)
  x$estimates["elpd_loo", "Estimate"]  <- hh_elpd_loo$y_hat_ppz
  x$estimates["elpd_loo", "SE"]  <- sqrt(hh_elpd_loo$hat_v_y_ppz)
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(hh_elpd_loo$v_hat_y_ppz)

  hh_p_loo <- whhest(z = pis_sample, m_i = x$pointwise$m_i, y = x$pointwise$p_loo, N)
  x$estimates["p_loo", "Estimate"] <- hh_p_loo$y_hat_ppz
  x$estimates["p_loo", "SE"] <- sqrt(hh_p_loo$hat_v_y_ppz)
  x$estimates["p_loo", "subsampling SE"] <- sqrt(hh_p_loo$v_hat_y_ppz)

  x <- update_psis_loo_ss_object(x)
  x
}

#' Update object with estimates
update_psis_loo_ss_object <- function(x){
  checkmate::assert_class(x, "psis_loo_ss")

  warning("update_psis_loo_ss_object not implemented")
  # 1. Compute looic
  # 2. Update x$p_loo,  x$se_p_loo, x$looic,  x$se_looic, x$elpd_loo, x$se_elpd_loo
#  x$estimates["looic", "Estimate"] <- x$looic <- looic_loo_est$y_hat
#  x$estimates["looic", "SE"] <- x$se_looic <- sqrt(looic_loo_est$hat_v_y)
#  x$quick_psis_loo$se_estimates["looic"] <- sqrt(looic_loo_est$v_y_hat)
  x
}

#' Weighted Hansen-Hurwitz estimator
#' @param z Normalized probabilities for the observation
#' @param m_i The number of times obs i was selected
#' @param y the values observed
#' @param N total number of observations in finite population
whhest <- function(z, m_i, y, N){
  checkmate::assert_numeric(z, lower = 0, upper = 1)
  checkmate::assert_numeric(y, len = length(z))
  checkmate::assert_integerish(m_i, len = length(z))
  est_list <- list(m = sum(m_i))
  est_list$y_hat_ppz <- sum(m_i*(y/z))/est_list$m
  est_list$v_hat_y_ppz <- (sum(m_i*((y/z - est_list$y_hat_ppz)^2))/est_list$m)/(est_list$m-1)

  # See unbiadness proof in supplementary material to the article
  est_list$hat_v_y_ppz <-
    (sum(m_i*(y^2/z)) / est_list$m) +
    est_list$v_hat_y_ppz / N - est_list$y_hat_ppz^2 / N
  est_list
}


#' Estimate elpd using the difference estimator and srs wor
#' @param x a quick_psis_loo object
loo_subsample_estimation_diff_srs <- function(x){
  checkmate::assert_class(x, "psis_loo_ss")

  elpd_loo_est <- srs_diff_est(x$subsamling_loo$elpd_loo_approx, y = x$pointwise[, "elpd_loo"], y_idx = x$pointwise[, "idx"])
  x$estimates["elpd_loo", "Estimate"] <- elpd_loo_est$y_hat
  x$estimates["elpd_loo", "SE"] <- sqrt(elpd_loo_est$hat_v_y)
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(elpd_loo_est$v_y_hat)

  #  elpd_p_est <- srs_diff_est(x$quick_psis_loo$approx_loo, y = x$pointwise$elpd_loo, y_idx = x$pointwise$idx)
  x$estimates["p_loo", "Estimate"] <- NA
  x$estimates["p_loo", "SE"] <-  NA
  x$estimates["p_loo", "subsampling SE"] <- NA
  warning("not implemented p_loo computation")

  x <- update_psis_loo_ss_object(x)

  x
}

#' Difference estimation using SRS-WOR sampling
#' @param y_approx Approximated values of all observations
#' @param y the values observed
#' @param y_idx the index of y in y_approx
srs_diff_est <- function(y_approx, y, y_idx){
  checkmate::assert_numeric(y_approx)
  checkmate::assert_numeric(y, max.len = length(y_approx))
  checkmate::assert_integerish(y_idx, len = length(y))

  N <- length(y_approx)
  m <- length(y)
  y_approx_m <- y_approx[y_idx]

  e_i <- y - y_approx_m
  t_pi_tilde <- sum(y_approx)
  t_pi2_tilde <- sum(y_approx^2)
  t_e <- N * mean(e_i)
  t_hat_epsilon <- N * mean(y^2 - y_approx_m^2)

  est_list <- list(m = length(y), N = N)
  est_list$y_hat <- t_pi_tilde + t_e
  est_list$v_y_hat <- N^2 * (1 - m / N) * var(e_i) / m
  est_list$hat_v_y <- (t_pi2_tilde + t_hat_epsilon) - # a (has been checked)
    (1/N) * (t_e^2 - est_list$v_y_hat + 2 * t_pi_tilde * est_list$y_hat - t_pi_tilde^2) # b
  est_list
}


#' Estimate elpd using the standard SRS estimator and SRS WOR
#' @param x a quick_psis_loo object
loo_estimation_srs <- function(x){
  checkmate::assert_class(x, "psis_loo_ss")

  elpd_loo_est <- srs_est(y = x$pointwise[, "elpd_loo"], y_approx = x$subsamling_loo$elpd_loo_approx)
  x$estimates["elpd_loo", "Estimate"] <- elpd_loo_est$y_hat
  x$estimates["elpd_loo", "SE"] <- sqrt(elpd_loo_est$hat_v_y)
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(elpd_loo_est$v_y_hat)

  p_loo_est <- srs_est(y = x$pointwise[, "p_loo"], y_approx = x$subsamling_loo$elpd_loo_approx)
  x$estimates["p_loo", "Estimate"] <- p_loo_est$y_hat
  x$estimates["p_loo", "SE"] <- sqrt(p_loo_est$hat_v_y)
  x$estimates["p_loo", "subsampling SE"] <- sqrt(p_loo_est$v_y_hat)


  x <- update_psis_loo_ss_object(x)

  x
}

#' Simple SRS-WOR estimation
#' @param y the values observed
#' @param y_approx a vector of length N
srs_est <- function(y, y_approx){
  checkmate::assert_numeric(y)
  checkmate::assert_numeric(y_approx, min.len = length(y))
  N <- length(y_approx)
  m <- length(y)
  est_list <- list(m = m)
  est_list$y_hat <- N * mean(y)
  est_list$v_y_hat <- N^2 * (1-m/N) * var(y)/m
  est_list$hat_v_y <- N * var(y)

  est_list
}


## Print methods ---

# TODO: Fix  this


