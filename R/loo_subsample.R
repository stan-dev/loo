#' Efficient approximate leave-one-out cross-validation (LOO) using subsampling,
#' so that less costly and more approximate computation is made for all LOO-fold,
#' and more costly and accurate computations are made only for m<N LOO-folds.
#'
#' @param x A function. The **Methods (by class)** section, below, has detailed
#'   descriptions of how to specify the inputs.
#'
#' @inheritParams loo
#' @param save_psis Should the `"psis"` object created internally by
#'   `loo_subsample()` be saved in the returned object? See [loo()] for details.
#' @template cores
#'
#' @details The `loo_subsample()` function is an S3 generic and a methods is
#'   currently provided for log-likelihood functions. The implementation works
#'   for both MCMC and for posterior approximations where it is possible to
#'   compute the log density for the approximation.
#'
#' @return `loo_subsample()` returns a named list with class `c("psis_loo_ss",
#'   "psis_loo", "loo")`. This has the same structure as objects returned by
#'   [loo()] but with the additional slot:
#'   * `loo_subsampling`: A list with two vectors, `log_p` and `log_g`, of the
#'   same length containing the posterior density and the approximation density
#'   for the individual draws.
#'
#' @seealso [loo()], [psis()], [loo_compare()]
#' @template loo-large-data-references
#'
#' @export loo_subsample loo_subsample.function
#'
loo_subsample <- function(x, ...) {
  UseMethod("loo_subsample")
}

#' @export
#' @templateVar fn loo_subsample
#' @template function
#' @param data,draws,... For `loo_subsample.function()`, these are the data,
#'   posterior draws, and other arguments to pass to the log-likelihood
#'   function. Note that for some `loo_approximation`s, the draws will be replaced
#'   by the posteriors summary statistics to compute loo approximations. See
#'   argument `loo_approximation` for details.
#' @param observations The subsample observations to use. The argument can take
#'   four (4) types of arguments:
#'   * `NULL` to use all observations. The algorithm then just uses
#'     standard `loo()` or `loo_approximate_posterior()`.
#'   * A single integer to specify the number of observations to be subsampled.
#'   * A vector of integers to provide the indices used to subset the data.
#'     _These observations need to be subsampled with the same scheme as given by
#'     the `estimator` argument_.
#'   * A `psis_loo_ss` object to use the same observations that were used in a
#'     previous call to `loo_subsample()`.
#'
#' @param log_p,log_g Should be supplied only if approximate posterior draws are
#'   used. The default (`NULL`) indicates draws are from "true" posterior (i.e.
#'   using MCMC). If not `NULL` then they should be specified as described in
#'   [loo_approximate_posterior()].
#'
#' @param loo_approximation What type of approximation of the loo_i's should be used?
#'   The default is `"plpd"` (the log predictive density using the posterior expectation).
#'   There are six different methods implemented to approximate loo_i's
#'   (see the references for more details):
#'   * `"plpd"`: uses the lpd based on point estimates (i.e., \eqn{p(y_i|\hat{\theta})}).
#'   * `"lpd"`: uses the lpds (i,e., \eqn{p(y_i|y)}).
#'   * `"tis"`: uses truncated importance sampling to approximate PSIS-LOO.
#'   * `"waic"`: uses waic (i.e., \eqn{p(y_i|y) - p_{waic}}).
#'   * `"waic_grad_marginal"`: uses waic approximation using first order delta
#'     method and posterior marginal variances to approximate \eqn{p_{waic}} (ie.
#'     \eqn{p(y_i|\hat{\theta})}-p_waic_grad_marginal). Requires gradient of
#'     likelihood function.
#'   * `"waic_grad"`: uses waic approximation using first order delta method and
#'     posterior covariance to approximate \eqn{p_{waic}} (ie.
#'     \eqn{p(y_i|\hat{\theta})}-p_waic_grad). Requires gradient of likelihood
#'     function.
#'   * `"waic_hess"`: uses waic approximation using second order delta method and
#'     posterior covariance to approximate \eqn{p_{waic}} (ie.
#'     \eqn{p(y_i|\hat{\theta})}-p_waic_grad). Requires gradient and Hessian of
#'     likelihood function.
#'
#'  As point estimates of \eqn{\hat{\theta}}, the posterior expectations
#'  of the parameters are used.
#'
#' @param loo_approximation_draws The number of posterior draws used when
#'   integrating over the posterior. This is used if `loo_approximation` is set
#'   to `"lpd"`, `"waic"`, or `"tis"`.
#'
#' @param estimator How should `elpd_loo`, `p_loo` and `looic` be estimated?
#'  The default is `"diff_srs"`.
#'  * `"diff_srs"`: uses the difference estimator with simple random sampling
#'    without replacement (srs). `p_loo` is estimated using standard srs.
#'    (Magnusson et al., 2020)
#'  * `"hh"`: uses the Hansen-Hurwitz estimator with sampling with replacement
#'    proportional to size, where `abs` of loo_approximation is used as size.
#'    (Magnusson et al., 2019)
#'  * `"srs"`: uses simple random sampling and ordinary estimation.
#'
#' @param llgrad The gradient of the log-likelihood. This
#'   is only used when `loo_approximation` is `"waic_grad"`,
#'   `"waic_grad_marginal"`, or `"waic_hess"`. The default is `NULL`.
#' @param llhess The Hessian of the log-likelihood. This is only used
#'        with `loo_approximation = "waic_hess"`. The default is `NULL`.
#'
loo_subsample.function <-
  function(x,
           ...,
           data = NULL,
           draws = NULL,
           observations = 400,
           log_p = NULL,
           log_g = NULL,
           r_eff = 1,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           loo_approximation = "plpd",
           loo_approximation_draws = NULL,
           estimator = "diff_srs",
           llgrad = NULL,
           llhess = NULL) {
    cores <- loo_cores(cores)
    # Asserting inputs
    .llfun <- validate_llfun(x)
    stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
    observations <- assert_observations(observations,
                                        N = dim(data)[1],
                                        estimator)
    checkmate::assert_numeric(log_p, len = length(log_g), null.ok = TRUE)
    checkmate::assert_null(dim(log_p))
    checkmate::assert_numeric(log_g, len = length(log_p), null.ok = TRUE)
    checkmate::assert_null(dim(log_g))

    if (is.null(log_p) && is.null(log_g)) {
        r_eff <- prepare_psis_r_eff(r_eff, len = dim(data)[1])
    }
    checkmate::assert_flag(save_psis)
    cores <- loo_cores(cores)

    checkmate::assert_choice(loo_approximation, choices = loo_approximation_choices(), null.ok = FALSE)
    checkmate::assert_int(loo_approximation_draws, lower = 1, upper = .ndraws(draws), null.ok = TRUE)
    checkmate::assert_choice(estimator, choices = estimator_choices())

    .llgrad <- .llhess <- NULL
    if (!is.null(llgrad)) .llgrad <- validate_llfun(llgrad)
    if (!is.null(llhess)) .llhess <- validate_llfun(llhess)

    # Fallbacks
    if (is.null(observations)) {
      if (is.null(log_p) && is.null(log_g)) {
        loo_obj <- loo.function(
          .llfun,
          ...,
          data = data,
          draws = draws,
          r_eff = r_eff,
          save_psis = save_psis,
          cores = cores
        )
      } else {
        loo_obj <- loo_approximate_posterior.function(
          .llfun,
          ...,
          log_p = log_p,
          log_g = log_g,
          data = data,
          draws = draws,
          save_psis = save_psis,
          cores = cores
        )
      }
      return(loo_obj)
    }

    # Compute loo approximation
    elpd_loo_approx <- elpd_loo_approximation(
      .llfun = .llfun,
      data = data,
      draws = draws,
      cores = cores,
      loo_approximation = loo_approximation,
      loo_approximation_draws = loo_approximation_draws,
      .llgrad = .llgrad,
      .llhess = .llhess
    )

    # Draw subsample of observations
    if (length(observations) == 1) {
      # Compute idxs
      idxs <- subsample_idxs(
        estimator = estimator,
        elpd_loo_approximation = elpd_loo_approx,
        observations = observations
      )
    } else {
      # Compute idxs
      idxs <- compute_idxs(observations)
    }
    data_subsample <- data[idxs$idx,, drop = FALSE]
    if (length(r_eff) > 1) {
      r_eff <- r_eff[idxs$idx]
    }

    # Compute elpd_loo
    if (!is.null(log_p) && !is.null(log_g)) {
      loo_obj <- loo_approximate_posterior.function(
        x = .llfun,
        data = data_subsample,
        draws = draws,
        log_p = log_p,
        log_g = log_g,
        save_psis = save_psis,
        cores = cores
      )
    } else {
      loo_obj <- loo.function(
        x = .llfun,
        data = data_subsample,
        draws = draws,
        r_eff = r_eff,
        save_psis = save_psis,
        cores = cores
      )
    }

    # Construct ss object and estimate
    loo_ss <- psis_loo_ss_object(x = loo_obj,
                                 idxs = idxs,
                                 elpd_loo_approx = elpd_loo_approx,
                                 loo_approximation = loo_approximation,
                                 loo_approximation_draws = loo_approximation_draws,
                                 estimator = estimator,
                                 .llfun = .llfun,
                                 .llgrad = .llgrad,
                                 .llhess = .llhess,
                                 data_dim = dim(data),
                                 ndraws = .ndraws(draws))
    loo_ss
  }


#' Update `psis_loo_ss` objects
#'
#' @details
#' If `observations` is updated then if a vector of indices or a `psis_loo_ss`
#' object is supplied the updated object will have exactly the observations
#' indicated by the vector or `psis_loo_ss` object. If a single integer is
#' supplied, new observations will be sampled to reach the supplied sample size.
#'
#' @export
#' @inheritParams loo_subsample.function
#' @param data,draws See [loo_subsample.function()].
#' @param object A `psis_loo_ss` object to update.
#' @param ... Currently not used.
#' @return A `psis_loo_ss` object.
#' @importFrom stats update
update.psis_loo_ss <- function(object, ...,
                               data = NULL,
                               draws = NULL,
                               observations = NULL,
                               r_eff = 1,
                               cores = getOption("mc.cores", 1),
                               loo_approximation = NULL,
                               loo_approximation_draws = NULL,
                               llgrad = NULL,
                               llhess = NULL) {
  # Fallback
  if (is.null(observations) &
     is.null(loo_approximation) &
     is.null(loo_approximation_draws) &
     is.null(llgrad) &
     is.null(llhess)) return(object)

  if (!is.null(data)) {
    stopifnot(is.data.frame(data) || is.matrix(data))
    checkmate::assert_true(all(dim(data) == object$loo_subsampling$data_dim))
  }
  if (!is.null(draws)) {
    # No current checks
  }
  cores <- loo_cores(cores)

  # Update elpd approximations
  if (!is.null(loo_approximation) | !is.null(loo_approximation_draws)) {
    stopifnot(is.data.frame(data) || is.matrix(data) & !is.null(draws))
    if (object$loo_subsampling$estimator %in% "hh_pps") {
      # HH estimation uses elpd_loo approx to sample,
      # so updating it will lead to incorrect results
      stop("Can not update loo_approximation when using PPS sampling.", call. = FALSE)
    }
    if (is.null(loo_approximation)) loo_approximation <- object$loo_subsampling$loo_approximation
    if (is.null(loo_approximation_draws)) loo_approximation_draws <- object$loo_subsampling$loo_approximation_draws
    if (is.null(llgrad)) .llgrad <- object$loo_subsampling$.llgrad else .llgrad <- validate_llfun(llgrad)
    if (is.null(llhess)) .llhess <- object$loo_subsampling$.llhess else .llhess <- validate_llfun(llhess)

    # Compute loo approximation
    elpd_loo_approx <-
      elpd_loo_approximation(.llfun = object$loo_subsampling$.llfun,
                             data = data, draws = draws,
                             cores = cores,
                             loo_approximation = loo_approximation,
                             loo_approximation_draws = loo_approximation_draws,
                             .llgrad = .llgrad, .llhess = .llhess)
    # Update object
    object$loo_subsampling$elpd_loo_approx <- elpd_loo_approx
    object$loo_subsampling$loo_approximation <- loo_approximation
    object$loo_subsampling["loo_approximation_draws"] <- list(loo_approximation_draws)
    object$loo_subsampling$.llgrad <- .llgrad
    object$loo_subsampling$.llhess <- .llhess
    object$pointwise[, "elpd_loo_approx"] <- object$loo_subsampling$elpd_loo_approx[object$pointwise[, "idx"]]
  }

  # Update observations
  if (!is.null(observations)) {
    observations <- assert_observations(observations,
                                        N = object$loo_subsampling$data_dim[1],
                                        object$loo_subsampling$estimator)
    if (length(observations) == 1) {
      checkmate::assert_int(observations, lower = nobs(object) + 1)
      stopifnot(is.data.frame(data) || is.matrix(data) & !is.null(draws))
    }

    # Compute subsample indices
    if (length(observations) > 1) {
      idxs <- compute_idxs(observations)
    } else {
      current_obs <- nobs(object)

      # If sampling with replacement
      if (object$loo_subsampling$estimator %in% c("hh_pps")) {
        idxs <- subsample_idxs(estimator = object$loo_subsampling$estimator,
                               elpd_loo_approximation = object$loo_subsampling$elpd_loo_approx,
                               observations = observations - current_obs)
      }
      # If sampling without replacement
      if (object$loo_subsampling$estimator %in% c("diff_srs", "srs")) {
        current_idxs <- obs_idx(object, rep = FALSE)
        new_idx <- (1:length(object$loo_subsampling$elpd_loo_approx))[-current_idxs]
        idxs <- subsample_idxs(estimator = object$loo_subsampling$estimator,
                               elpd_loo_approximation = object$loo_subsampling$elpd_loo_approx[-current_idxs],
                               observations = observations - current_obs)
        idxs$idx <- new_idx[idxs$idx]
      }
    }

    # Identify how to update object
    cidxs <- compare_idxs(idxs, object)

    # Compute new observations
    if (!is.null(cidxs$new)) {
      stopifnot(is.data.frame(data) || is.matrix(data) & !is.null(draws))
      data_new_subsample <- data[cidxs$new$idx,, drop = FALSE]
      if (length(r_eff) > 1) r_eff <- r_eff[cidxs$new$idx]

      if (!is.null(object$approximate_posterior$log_p) & !is.null(object$approximate_posterior$log_g)) {
        loo_obj <- loo_approximate_posterior.function(x = object$loo_subsampling$.llfun,
                                                  data = data_new_subsample,
                                                  draws = draws,
                                                  log_p = object$approximate_posterior$log_p,
                                                  log_g = object$approximate_posterior$log_g,
                                                  save_psis = !is.null(object$psis_object),
                                                  cores = cores)
      } else {
        loo_obj <- loo.function(x = object$loo_subsampling$.llfun,
                            data = data_new_subsample,
                            draws = draws,
                            r_eff = r_eff,
                            save_psis = !is.null(object$psis_object),
                            cores = cores)
      }
      # Add stuff to pointwise
      loo_obj$pointwise <-
        add_subsampling_vars_to_pointwise(loo_obj$pointwise,
                                          cidxs$new,
                                          object$loo_subsampling$elpd_loo_approx)
    } else {
      loo_obj <- NULL
    }

    if (length(observations) == 1) {
      # Add new samples pointwise and diagnostic
      object <- rbind_psis_loo_ss(object, x = loo_obj)

      # Update m_i for current pointwise (diagnostic stay the same)
      object$pointwise <- update_m_i_in_pointwise(object$pointwise, cidxs$add, type = "add")
    } else {
      # Add new samples pointwise and diagnostic
      object <- rbind_psis_loo_ss(object, loo_obj)

      # Replace m_i current pointwise and diagnostics
      object$pointwise <- update_m_i_in_pointwise(object$pointwise, cidxs$add, type = "replace")

      # Remove samples
      object <- remove_idx.psis_loo_ss(object, idxs = cidxs$remove)
      stopifnot(setequal(obs_idx(object), observations))

      # Order object as in observations
      object <- order.psis_loo_ss(object, observations)
    }
  }


  # Compute estimates
  if (object$loo_subsampling$estimator == "hh_pps") {
    object <- loo_subsample_estimation_hh(object)
  } else if (object$loo_subsampling$estimator == "diff_srs") {
    object <- loo_subsample_estimation_diff_srs(object)
  } else if (object$loo_subsampling$estimator == "srs") {
    object <- loo_subsample_estimation_srs(object)
  } else {
    stop("No correct estimator used.")
  }

  assert_psis_loo_ss(object)
  object
}

#' Get observation indices used in subsampling
#'
#' @param x A `psis_loo_ss` object.
#' @param rep If sampling with replacement is used, an observation can have
#'   multiple samples and these are then repeated in the returned object if
#'   `rep=TRUE` (e.g., a vector `c(1,1,2)` indicates that observation 1 has been
#'   subampled two times). If `rep=FALSE` only the unique indices are returned.
#'
#' @return An integer vector.
#'
#' @export
obs_idx <- function(x, rep = TRUE) {
  checkmate::assert_class(x, "psis_loo_ss")
  if (rep) {
    idxs <- as.integer(rep(x$pointwise[,"idx"], x$pointwise[,"m_i"]))
  } else {
    idxs <- as.integer(x$pointwise[,"idx"])
  }
  idxs
}

#' The number of observations in a `psis_loo_ss` object.
#' @importFrom stats nobs
#' @param object a `psis_loo_ss` object.
#' @param ... Currently unused.
#' @export
nobs.psis_loo_ss <- function(object, ...) {
  as.integer(sum(object$pointwise[,"m_i"]))
}

# internal ----------------------------------------------------------------

#' The possible choices of loo_approximations implemented
#'
#' @details
#' The choice `psis` is returned if a `psis_loo` object
#' is converted to a `psis_loo_ss` object with `as.psis_loo_ss()`.
#' But `psis` cannot be chosen in the API of `loo_subsample()`.
#'
#' @noRd
#' @param api The choices available in the loo API or all possible choices.
#' @return A character vector of allowed choices.
loo_approximation_choices <- function(api = TRUE) {
  lac <- c("plpd", "lpd", "waic", "waic_grad_marginal", "waic_grad", "waic_hess", "tis", "sis", "none")
  if (!api) lac <- c(lac, "psis")
  lac
}

#' The estimators implemented
#'
#' @noRd
#' @return A character vector of allowed choices.
estimator_choices <- function() {
  c("hh_pps", "diff_srs", "srs")
}

## Approximate elpd -----

#' Utility function to apply user-specified log-likelihood to a single data point
#' @details
#' See [elpd_loo_approximation] and [compute_lpds] for usage examples
#' @noRd
#'
#' @return lpd value for a single data point i
lpd_i <- function(i, llfun, data, draws) {
  ll_i <- llfun(data_i = data[i,, drop=FALSE], draws = draws)
  ll_i <- as.vector(ll_i)
  lpd_i <- logMeanExp(ll_i)
  lpd_i
}


#' Utility function to compute lpd using user-defined likelihood function
#' using platform-dependent parallel backends when cores > 1
#'
#' @details
#' See [elpd_loo_approximation] for usage examples
#'
#' @noRd
#' @return a vector of computed log probability densities
compute_lpds <- function(N, data, draws, llfun, cores) {
  if (cores == 1) {
    lpds <- lapply(X = seq_len(N), FUN = lpd_i, llfun, data, draws)
  } else {
    if (.Platform$OS.type != "windows") {
      lpds <- mclapply(X = seq_len(N), mc.cores = cores, FUN = lpd_i, llfun, data, draws)
    } else {
      cl <- makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      lpds <- parLapply(cl, X = seq_len(N), fun = lpd_i, llfun, data, draws)
    }
  }

  unlist(lpds)
}

#' Compute approximation to loo_i:s
#'
#' @details
#' See [loo_subsample.function()] and the `loo_approximation` argument.
#' @noRd
#' @inheritParams loo_subsample.function
#'
#' @return a vector with approximations of elpd_{loo,i}s
elpd_loo_approximation <- function(.llfun, data, draws, cores, loo_approximation, loo_approximation_draws = NULL, .llgrad = NULL, .llhess = NULL) {
  checkmate::assert_function(.llfun, args = c("data_i", "draws"), ordered = TRUE)
  stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
  checkmate::assert_choice(loo_approximation, choices = loo_approximation_choices(), null.ok = FALSE)
  checkmate::assert_int(loo_approximation_draws, lower = 2, null.ok = TRUE)
  if (!is.null(.llgrad)) {
    checkmate::assert_function(.llgrad, args = c("data_i", "draws"), ordered = TRUE)
  }
  if (!is.null(.llhess)) {
    checkmate::assert_function(.llhess, args = c("data_i", "draws"), ordered = TRUE)
  }

  cores <- loo_cores(cores)
  N <- dim(data)[1]

  if (loo_approximation == "none") return(rep(1L,N))

  if (loo_approximation %in% c("tis", "sis")) {
    draws <- .thin_draws(draws, loo_approximation_draws)
    is_values <- suppressWarnings(loo.function(.llfun, data = data, draws = draws, is_method = loo_approximation))
    return(is_values$pointwise[, "elpd_loo"])
  }

  if (loo_approximation == "waic") {
    draws <- .thin_draws(draws, loo_approximation_draws)
    waic_full_obj <- waic.function(.llfun, data = data, draws = draws)
    return(waic_full_obj$pointwise[,"elpd_waic"])
  }

  # Compute the lpd or log p(y_i|y_{-i})
  if (loo_approximation == "lpd") {
    draws <- .thin_draws(draws, loo_approximation_draws)
    lpds <- compute_lpds(N, data, draws, .llfun, cores)
    return(lpds) # Use only the lpd
  }

  # Compute the point lpd or log p(y_i|\hat{\theta}) - also used in waic_delta approaches
  if (loo_approximation == "plpd" |
      loo_approximation == "waic_grad" |
      loo_approximation == "waic_grad_marginal" |
      loo_approximation == "waic_hess") {

    draws <- .thin_draws(draws, loo_approximation_draws)
    point_est <- .compute_point_estimate(draws)
    lpds <- compute_lpds(N, data, point_est, .llfun, cores)
    if (loo_approximation == "plpd") return(lpds) # Use only the lpd
  }

  if (loo_approximation == "waic_grad" |
      loo_approximation == "waic_grad_marginal" |
      loo_approximation == "waic_hess") {
    checkmate::assert_true(!is.null(.llgrad))

    point_est <- .compute_point_estimate(draws)
    # Compute the lpds
    lpds <- compute_lpds(N, data, point_est, .llfun, cores)

    if (loo_approximation == "waic_grad" |
        loo_approximation == "waic_hess") {
      cov_est <- stats::cov(draws)
    }

    if (loo_approximation == "waic_grad_marginal") {
      marg_vars <- apply(draws, MARGIN = 2, var)
    }

    p_eff_approx <- numeric(N)
    if (cores>1) warning("Multicore is not implemented for waic_delta",
                         call. = FALSE)

    if (loo_approximation == "waic_grad") {
      for(i in 1:nrow(data)) {
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        local_cov <- cov_est[rownames(grad_i), rownames(grad_i)]
        p_eff_approx[i] <-  t(grad_i) %*% local_cov %*% grad_i
      }
    } else if (loo_approximation == "waic_grad_marginal") {
      for(i in 1:nrow(data)) {
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        p_eff_approx[i] <- sum(grad_i * marg_vars[rownames(grad_i)] * grad_i)
      }
    } else if (loo_approximation == "waic_hess") {
      checkmate::assert_true(!is.null(.llhess))
      for(i in 1:nrow(data)) {
        grad_i <- t(.llgrad(data[i,,drop = FALSE], point_est))
        hess_i <- .llhess(data_i = data[i,,drop = FALSE], draws = point_est[,rownames(grad_i), drop = FALSE])[,,1]
        local_cov <- cov_est[rownames(grad_i), rownames(grad_i)]
        p_eff_approx[i] <- t(grad_i) %*% local_cov %*% grad_i +
          0.5 * sum(diag(local_cov %*% hess_i %*% local_cov %*% hess_i))
      }
    } else {
      stop(loo_approximation, " is not implemented!", call. = FALSE)
    }
    return(lpds - p_eff_approx)
  }

}


#' Compute a point estimate from a draws object
#'
#' @keywords internal
#' @export
#' @details This is a generic function to compute point estimates from draws
#'   objects. The function is internal and should only be used by developers to
#'   enable [loo_subsample()] for arbitrary draws objects.
#'
#' @param draws A draws object with draws from the posterior.
#' @return A 1 by P matrix with point estimates from a draws object.
.compute_point_estimate <- function(draws) {
  UseMethod(".compute_point_estimate")
}
#' @rdname dot-compute_point_estimate
#' @export
.compute_point_estimate.matrix <- function(draws) {
  t(as.matrix(colMeans(draws)))
}
#' @rdname dot-compute_point_estimate
#' @export
.compute_point_estimate.default <- function(draws) {
  stop(".compute_point_estimate() has not been implemented for objects of class '", class(draws), "'")
}

#' Thin a draws object
#'
#' @keywords internal
#' @export
#' @details This is a generic function to thin draws from arbitrary draws
#'   objects. The function is internal and should only be used by developers to
#'   enable [loo_subsample()] for arbitrary draws objects.
#'
#' @param draws A draws object with posterior draws.
#' @param loo_approximation_draws The number of posterior draws to return (ie after thinning).
#' @return A thinned draws object.
.thin_draws <- function(draws, loo_approximation_draws) {
  UseMethod(".thin_draws")
}
#' @rdname dot-thin_draws
#' @export
.thin_draws.matrix <- function(draws, loo_approximation_draws) {
  if (is.null(loo_approximation_draws)) return(draws)
  checkmate::assert_int(loo_approximation_draws, lower = 1, upper = .ndraws(draws), null.ok = TRUE)
  S <- .ndraws(draws)
  idx <- 1:loo_approximation_draws * S %/% loo_approximation_draws
  draws <- draws[idx, , drop = FALSE]
  draws
}
#' @rdname dot-thin_draws
#' @export
.thin_draws.numeric <- function(draws, loo_approximation_draws) {
  .thin_draws.matrix(as.matrix(draws), loo_approximation_draws)
}
#' @rdname dot-thin_draws
#' @export
.thin_draws.default <- function(draws, loo_approximation_draws) {
  stop(".thin_draws() has not been implemented for objects of class '", class(draws), "'")
}


#' The number of posterior draws in a draws object.
#'
#' @keywords internal
#' @export
#' @details This is a generic function to return the total number of draws from
#'   an arbitrary draws objects. The function is internal and should only be
#'   used by developers to enable [loo_subsample()] for arbitrary draws objects.
#'
#' @param x A draws object with posterior draws.
#' @return An integer with the number of draws.
.ndraws <- function(x) {
  UseMethod(".ndraws")
}
#' @rdname dot-ndraws
#' @export
.ndraws.matrix <- function(x) {
  nrow(x)
}
#' @rdname dot-ndraws
#' @export
.ndraws.default <- function(x) {
  stop(".ndraws() has not been implemented for objects of class '", class(x), "'")
}

## Subsampling -----

#' Subsampling strategy
#'
#' @noRd
#' @param estimator The estimator to use, see `estimator_choices()`.
#' @param elpd_loo_approximation A vector of loo approximations, see `elpd_loo_approximation()`.
#' @param observations The total number of subsample observations to sample.
#' @return A `subsample_idxs` data frame.
subsample_idxs <- function(estimator, elpd_loo_approximation, observations) {
  checkmate::assert_choice(estimator, choices = estimator_choices())
  checkmate::assert_numeric(elpd_loo_approximation)
  checkmate::assert_int(observations)

  if (estimator == "hh_pps") {
    pi_values <- pps_elpd_loo_approximation_to_pis(elpd_loo_approximation)
    idxs_df <- pps_sample(observations, pis = pi_values)
  }

  if (estimator == "diff_srs" | estimator == "srs") {
    if (observations > length(elpd_loo_approximation)) {
      stop("'observations' is larger than the total sample size in 'data'.", call. = FALSE)
    }
    idx <- 1:length(elpd_loo_approximation)
    idx_m <- idx[order(stats::runif(length(elpd_loo_approximation)))][1:observations]
    idx_m <- idx_m[order(idx_m)]
    idxs_df <- data.frame(idx=as.integer(idx_m), m_i=1L)
  }
  assert_subsample_idxs(x = idxs_df)
  idxs_df
}

#' Compute pis from approximation for use in pps sampling.
#' @noRd
#' @details pis are the sampling probabilities and sum to 1.
#' @inheritParams subsample_idxs
#' @return A vector of pis.
pps_elpd_loo_approximation_to_pis <- function(elpd_loo_approximation) {
  checkmate::assert_numeric(elpd_loo_approximation)
  pi_values <- abs(elpd_loo_approximation)
  pi_values <- pi_values/sum(pi_values) # \tilde{\pi}
  pi_values
}


#' Compute subsampling indices from an observation vector
#' @noRd
#' @param observation A vector of indices.
#' @return A `subsample_idxs` data frame.
compute_idxs <- function(observations) {
  checkmate::assert_integer(observations, lower = 1, min.len = 2, any.missing = FALSE)
  tab <- table(observations)
  idxs_df <- data.frame(idx = as.integer(names(tab)), m_i = as.integer(unname(tab)))
  assert_subsample_idxs(idxs_df)
  idxs_df
}


#' Compare the indices to prepare handling
#'
#' @details
#' The function compares the object and sampled indices into `new`
#' (observations not in `object`), `add` (observations in `object`), and
#' `remove` (observations in `object` but not in idxs).
#' @noRd
#' @param idxs A `subsample_idxs` data frame.
#' @param object A `psis_loo_ss` object.
#' @return A list of three `subsample_idxs` data frames. Elements without any
#'   observations return `NULL`.
compare_idxs  <- function(idxs, object) {
  assert_subsample_idxs(idxs)
  current_idx <- compute_idxs(obs_idx(object))
  result <- list()
  new_idx <- !(idxs$idx %in% current_idx$idx)
  remove_idx <- !(current_idx$idx %in% idxs$idx)

  result$new <- idxs[new_idx, ]
  if (nrow(result$new) == 0) {
    result["new"] <- NULL
  } else {
    assert_subsample_idxs(result$new)
  }

  result$add <- idxs[!new_idx, ]
  if (nrow(result$add) == 0) {
    result["add"] <- NULL
  } else {
    assert_subsample_idxs(result$add)
  }

  result$remove <- current_idx[remove_idx, ]
  if (nrow(result$remove) == 0) {
    result["remove"] <- NULL
  } else {
    assert_subsample_idxs(result$remove)
  }

  result
}


#' Draw a PPS sample with replacement and return a idx_df
#' @noRd
#' @details
#' We are sampling with replacement, hence we only want to compute elpd
#' for each observation once.
#' @param m The total sampling size.
#' @param pis The probability of selecting each observation.
#' @return a `subsample_idxs` data frame.
pps_sample <- function(m, pis) {
  checkmate::assert_int(m)
  checkmate::assert_numeric(pis, min.len = 2, lower = 0, upper = 1)
  idx <- sample(1:length(pis), size = m, replace = TRUE, prob = pis)
  idxs_df <- as.data.frame(table(idx), stringsAsFactors = FALSE)
  colnames(idxs_df) <- c("idx", "m_i")
  idxs_df$idx <- as.integer(idxs_df$idx)
  idxs_df$m_i <- as.integer(idxs_df$m_i)
  assert_subsample_idxs(idxs_df)
  idxs_df
}

## Constructor ---

#' Construct a `psis_loo_ss` object
#'
#' @noRd
#' @param x A `psis_loo` object.
#' @param idxs a `subsample_idxs` data frame.
#' @param elpd_loo_approximation A vector of loo approximations, see
#'   `elpd_loo_approximation()`.
#' @inheritParams loo_subsample
#' @param .llfun,.llgrad,.llhess  See llfun, llgrad and llhess in `loo_subsample()`.
#' @param data_dim Dimension of the data object.
#' @param ndraws Dimension of the draws object.
#' @return A `psis_loo_ss` object.
psis_loo_ss_object <- function(x,
                               idxs,
                               elpd_loo_approx,
                               loo_approximation, loo_approximation_draws,
                               estimator,
                               .llfun, .llgrad, .llhess,
                               data_dim, ndraws) {
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
  checkmate::assert_integer(data_dim, len = 2, lower = 1, any.missing = FALSE)
  checkmate::assert_int(ndraws, lower = 1)

  # Construct object
  class(x) <- c("psis_loo_ss", class(x))
  x$pointwise <- add_subsampling_vars_to_pointwise(pointwise = x$pointwise, idxs, elpd_loo_approx)
  x$estimates <- cbind(x$estimates, matrix(0, nrow = nrow(x$estimates)))
  colnames(x$estimates)[ncol(x$estimates)] <- "subsampling SE"

  x$loo_subsampling <- list()
  x$loo_subsampling$elpd_loo_approx <- elpd_loo_approx
  x$loo_subsampling$loo_approximation <- loo_approximation
  x$loo_subsampling["loo_approximation_draws"] <- list(loo_approximation_draws)
  x$loo_subsampling$estimator <- estimator
  x$loo_subsampling$.llfun <- .llfun
  x$loo_subsampling[".llgrad"] <- list(.llgrad)
  x$loo_subsampling[".llhess"] <- list(.llhess)
  x$loo_subsampling$data_dim <- data_dim
  x$loo_subsampling$ndraws <- ndraws

  # Compute estimates
  if (estimator == "hh_pps") {
    x <- loo_subsample_estimation_hh(x)
  } else if (estimator == "diff_srs") {
    x <- loo_subsample_estimation_diff_srs(x)
  } else if (estimator == "srs") {
    x <- loo_subsample_estimation_srs(x)
  } else {
    stop("No correct estimator used.")
  }
  assert_psis_loo_ss(x)
  x
}

as.psis_loo_ss <- function(x) {
  UseMethod("as.psis_loo_ss")
}
#' @export
as.psis_loo_ss.psis_loo_ss <- function(x) {
  x
}
#' @export
as.psis_loo_ss.psis_loo <- function(x) {
  class(x) <- c("psis_loo_ss", class(x))
  x$estimates <- cbind(x$estimates, matrix(0, nrow = nrow(x$estimates)))
  colnames(x$estimates)[ncol(x$estimates)] <- "subsampling SE"
  x$pointwise <- cbind(x$pointwise,
                       matrix(1:nrow(x$pointwise), byrow = FALSE, ncol = 1),
                       matrix(rep(1,nrow(x$pointwise)), byrow = FALSE, ncol = 1),
                       x$pointwise[, "elpd_loo"])
  ncp <- ncol(x$pointwise)
  colnames(x$pointwise)[(ncp-2):ncp] <- c("idx", "m_i", "elpd_loo_approx")
  x$loo_subsampling <- list(elpd_loo_approx=x$pointwise[, "elpd_loo"],
                           loo_approximation = "psis",
                           loo_approximation_draws = NULL,
                           estimator = "diff_srs",
                           data_dim = c(nrow(x$pointwise), NA),
                           ndraws = NA)
  assert_psis_loo_ss(x)
  x
}

as.psis_loo <- function(x) {
  UseMethod("as.psis_loo")
}

#' @export
as.psis_loo.psis_loo <- function(x) {
  x
}
#' @export
as.psis_loo.psis_loo_ss <- function(x) {
  if (x$loo_subsampling$data_dim[1] == nrow(x$pointwise)) {
    x$estimates <- x$estimates[, 1:2]
    x$pointwise <- x$pointwise[, 1:5]
    x$loo_subsampling <- NULL
    loo_obj <- importance_sampling_loo_object(pointwise = x$pointwise[, 1:5],
                           diagnostics = x$diagnostics,
                           dims = attr(x, "dims"),
                           is_method = "psis",
                           is_object = x$psis_object)
    if (inherits(x, "psis_loo_ap")) {
      loo_obj$approximate_posterior <- list(log_p = x$approximate_posterior$log_p,
                                        log_g = x$approximate_posterior$log_g)
      class(loo_obj) <- c("psis_loo_ap", class(loo_obj))
      assert_psis_loo_ap(loo_obj)
    }
  } else {
    stop("A subsampling loo object can only be coerced to a loo object ",
         "if all observations in data have been subsampled.", call. = FALSE)
  }

  loo_obj
}

#' Add subsampling information to the pointwise element in a `psis_loo` object.
#' @noRd
#' @param pointwise The `pointwise` element in a `psis_loo` object.
#' @param idxs A `subsample_idxs` data frame.
#' @param elpd_loo_approximation A vector of loo approximations, see `elpd_loo_approximation()`.
#' @return A `pointwise` matrix with subsampling information.
add_subsampling_vars_to_pointwise <- function(pointwise, idxs, elpd_loo_approx) {
  checkmate::assert_matrix(pointwise,
                           any.missing = FALSE,
                           min.cols = 5)
  checkmate::assert_names(colnames(pointwise), identical.to = c("elpd_loo","mcse_elpd_loo","p_loo","looic", "influence_pareto_k"))
  assert_subsample_idxs(idxs)
  checkmate::assert_numeric(elpd_loo_approx)

  pw <- cbind(as.data.frame(pointwise), idxs)
  pw$elpd_loo_approx <- elpd_loo_approx[idxs$idx]
  pw <- as.matrix(pw)
  rownames(pw) <- NULL
  assert_subsampling_pointwise(pw)
  pw
}

#' Add `psis_loo` object to a `psis_loo_ss` object
#' @noRd
#' @param object A `psis_loo_ss` object.
#' @param x A `psis_loo` object.
#' @return An updated `psis_loo_ss` object.
rbind_psis_loo_ss <- function(object, x) {
  checkmate::assert_class(object, "psis_loo_ss")
  if (is.null(x)) return(object) # Fallback
  checkmate::assert_class(x, "psis_loo")
  assert_subsampling_pointwise(object$pointwise)
  assert_subsampling_pointwise(x$pointwise)
  checkmate::assert_disjunct(object$pointwise[, "idx"], x$pointwise[, "idx"])

  object$pointwise <- rbind(object$pointwise, x$pointwise)
  object$diagnostics$pareto_k <-
    c(object$diagnostics$pareto_k, x$diagnostics$pareto_k)
  object$diagnostics$n_eff <- c(object$diagnostics$n_eff, x$diagnostics$n_eff)
  object$diagnostics$r_eff <- c(object$diagnostics$r_eff, x$diagnostics$r_eff)
  attr(object, "dims")[2] <- nrow(object$pointwise)
  object
}

#' Remove observations in `idxs` from object
#' @noRd
#' @param object A `psis_loo_ss` object.
#' @param idxs A `subsample_idxs` data frame.
#' @return A `psis_loo_ss` object.
remove_idx.psis_loo_ss <- function(object, idxs) {
  checkmate::assert_class(object, "psis_loo_ss")
  if (is.null(idxs)) return(object) # Fallback
  assert_subsample_idxs(idxs)

  row_map <- data.frame(
    row_no = 1:nrow(object$pointwise),
    idx = object$pointwise[, "idx"]
  )
  row_map <- merge(row_map, idxs, by = "idx", all.y = TRUE)

  object$pointwise <- object$pointwise[-row_map$row_no,,drop = FALSE]
  object$diagnostics$pareto_k <- object$diagnostics$pareto_k[-row_map$row_no]
  object$diagnostics$n_eff <- object$diagnostics$n_eff[-row_map$row_no]
  object$diagnostics$r_eff <- object$diagnostics$r_eff[-row_map$row_no]
  attr(object, "dims")[2] <- nrow(object$pointwise)
  object
}

#' Order object by `observations`.
#' @noRd
#' @param x A `psis_loo_ss` object.
#' @param observations A vector with indices.
#' @return An ordered `psis_loo_ss` object.
order.psis_loo_ss <- function(x, observations) {
  checkmate::assert_class(x, "psis_loo_ss")
  checkmate::assert_integer(observations, len = nobs(x))
  if (identical(obs_idx(x), observations)) return(x) # Fallback
  checkmate::assert_set_equal(obs_idx(x), observations)

  row_map_x <- data.frame(row_no_x = 1:nrow(x$pointwise), idx = x$pointwise[, "idx"])
  row_map_obs <- data.frame(row_no_obs = 1:length(observations), idx = observations)
  row_map <- merge(row_map_obs, row_map_x, by = "idx", sort = FALSE)
  x$pointwise <- x$pointwise[row_map$row_no_x,,drop = FALSE]
  x$diagnostics$pareto_k <- x$diagnostics$pareto_k[row_map$row_no_x]
  x$diagnostics$n_eff <- x$diagnostics$n_eff[row_map$row_no_x]
  x$diagnostics$r_eff <- x$diagnostics$r_eff[row_map$row_no_x]
  x
}

#' Update m_i in a `pointwise` element.
#' @noRd
#' @param x A `psis_loo_ss` `pointwise` data frame.
#' @param idxs A `subsample_idxs` data frame.
#' @param type should the m_i:s in `idxs` `"replace"` the current m_i:s or
#'   `"add"` to them.
#' @return An ordered `psis_loo_ss` object.
update_m_i_in_pointwise <- function(pointwise, idxs, type = "replace") {
  assert_subsampling_pointwise(pointwise)
  if (is.null(idxs)) return(pointwise) # Fallback
  assert_subsample_idxs(idxs)
  checkmate::assert_choice(type, choices = c("replace", "add"))

  row_map <- data.frame(row_no = 1:nrow(pointwise), idx = pointwise[, "idx"])
  row_map <- merge(row_map, idxs, by = "idx", all.y = TRUE)

  if (type == "replace") {
    pointwise[row_map$row_no, "m_i"] <- row_map$m_i
  }
  if (type == "add") {
    pointwise[row_map$row_no, "m_i"] <- pointwise[row_map$row_no, "m_i"] + row_map$m_i
  }
  pointwise
}



## Estimation ---

#' Estimate the elpd using the Hansen-Hurwitz estimator (Magnusson et al., 2019)
#' @noRd
#' @param x A `psis_loo_ss` object.
#' @return A `psis_loo_ss` object.
loo_subsample_estimation_hh <- function(x) {
  checkmate::assert_class(x, "psis_loo_ss")
  N <- length(x$loo_subsampling$elpd_loo_approx)
  pis <- pps_elpd_loo_approximation_to_pis(x$loo_subsampling$elpd_loo_approx)
  pis_sample <- pis[x$pointwise[,"idx"]]

  hh_elpd_loo <- whhest(z = pis_sample, m_i = x$pointwise[, "m_i"], y = x$pointwise[, "elpd_loo"], N)
  srs_elpd_loo <- srs_est(y = x$pointwise[, "elpd_loo"], y_approx = pis_sample)
  x$estimates["elpd_loo", "Estimate"]  <- hh_elpd_loo$y_hat_ppz
  if (hh_elpd_loo$hat_v_y_ppz > 0) {
    x$estimates["elpd_loo", "SE"]  <- sqrt(hh_elpd_loo$hat_v_y_ppz)
  } else {
    warning("Negative estimate of SE, more subsampling obs. needed.", call. = FALSE)
    x$estimates["elpd_loo", "SE"]  <- NaN
  }
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(hh_elpd_loo$v_hat_y_ppz)

  hh_p_loo <- whhest(z = pis_sample, m_i = x$pointwise[,"m_i"], y = x$pointwise[,"p_loo"], N)
  x$estimates["p_loo", "Estimate"] <- hh_p_loo$y_hat_ppz
  if (hh_p_loo$hat_v_y_ppz > 0) {
    x$estimates["p_loo", "SE"]  <- sqrt(hh_p_loo$hat_v_y_ppz)
  } else {
    warning("Negative estimate of SE, more subsampling obs. needed.", call. = FALSE)
    x$estimates["elpd_loo", "SE"]  <- NaN
  }
  x$estimates["p_loo", "subsampling SE"] <- sqrt(hh_p_loo$v_hat_y_ppz)
  update_psis_loo_ss_estimates(x)
}

#' Update a `psis_loo_ss` object with generic estimates
#'
#' @noRd
#' @details
#' Updates a `psis_loo_ss` with generic estimates (looic)
#' and updates components in the object based on x$estimate.
#' @param x A `psis_loo_ss` object.
#' @return x A `psis_loo_ss` object.
update_psis_loo_ss_estimates <- function(x) {
  checkmate::assert_class(x, "psis_loo_ss")

  x$estimates["looic", "Estimate"] <- (-2) * x$estimates["elpd_loo", "Estimate"]
  x$estimates["looic", "SE"] <- 2 * x$estimates["elpd_loo", "SE"]
  x$estimates["looic", "subsampling SE"] <- 2 * x$estimates["elpd_loo", "subsampling SE"]

  x$elpd_loo <- x$estimates["elpd_loo", "Estimate"]
  x$p_loo <- x$estimates["p_loo", "Estimate"]
  x$looic <- x$estimates["looic", "Estimate"]
  x$se_elpd_loo <- x$estimates["elpd_loo", "SE"]
  x$se_p_loo <- x$estimates["p_loo", "SE"]
  x$se_looic <- x$estimates["looic", "SE"]

  x
}

#' Weighted Hansen-Hurwitz estimator (Magnusson et al., 2019)
#' @noRd
#' @param z Normalized probabilities for the observation.
#' @param m_i The number of times obs i was selected.
#' @param y The values observed.
#' @param N The total number of observations in finite population.
#' @return A list with estimates.
whhest <- function(z, m_i, y, N) {
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


#' Estimate elpd using the difference estimator and SRS-WOR (Magnusson et al., 2020)
#' @noRd
#' @param x A `psis_loo_ss` object.
#' @return A `psis_loo_ss` object.
loo_subsample_estimation_diff_srs <- function(x) {
  checkmate::assert_class(x, "psis_loo_ss")

  elpd_loo_est <- srs_diff_est(y_approx = x$loo_subsampling$elpd_loo_approx, y = x$pointwise[, "elpd_loo"], y_idx = x$pointwise[, "idx"])
  x$estimates["elpd_loo", "Estimate"] <- elpd_loo_est$y_hat
  x$estimates["elpd_loo", "SE"] <- sqrt(elpd_loo_est$hat_v_y)
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(elpd_loo_est$v_y_hat)

  p_loo_est <- srs_est(y = x$pointwise[, "p_loo"], y_approx = x$loo_subsampling$elpd_loo_approx)
  x$estimates["p_loo", "Estimate"] <- p_loo_est$y_hat
  x$estimates["p_loo", "SE"] <- sqrt(p_loo_est$hat_v_y)
  x$estimates["p_loo", "subsampling SE"] <- sqrt(p_loo_est$v_y_hat)

  update_psis_loo_ss_estimates(x)
}

#' Difference estimation using SRS-WOR sampling (Magnusson et al., 2020)
#' @noRd
#' @param y_approx Approximated values of all observations.
#' @param y The values observed.
#' @param y_idx The index of `y` in `y_approx`.
#' @return A list with estimates.
srs_diff_est <- function(y_approx, y, y_idx) {
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
  # eq (7)
  est_list$y_hat <- t_pi_tilde + t_e
  # eq (8)
  est_list$v_y_hat <- N^2 * (1 - m / N) * var(e_i) / m
  # eq (9) first row second `+` should be `-`
  # Supplementary material eq (6) has this correct
  # Here the variance is for sum, while in the paper the variance is for mean
  # which explains the proportional difference of 1/N
  est_list$hat_v_y <- (t_pi2_tilde + t_hat_epsilon) - # a (has been checked)
    (1/N) * (t_e^2 - est_list$v_y_hat + 2 * t_pi_tilde * est_list$y_hat - t_pi_tilde^2) # b
  est_list
}


#' Estimate elpd using the standard simple-re-sample without
#' resampling (SRS-WOR) estimator
#' @noRd
#' @param x A `psis_loo_ss` object.
#' @return A `psis_loo_ss` object.
loo_subsample_estimation_srs <- function(x) {
  checkmate::assert_class(x, "psis_loo_ss")

  elpd_loo_est <- srs_est(y = x$pointwise[, "elpd_loo"], y_approx = x$loo_subsampling$elpd_loo_approx)
  x$estimates["elpd_loo", "Estimate"] <- elpd_loo_est$y_hat
  x$estimates["elpd_loo", "SE"] <- sqrt(elpd_loo_est$hat_v_y)
  x$estimates["elpd_loo", "subsampling SE"] <- sqrt(elpd_loo_est$v_y_hat)

  p_loo_est <- srs_est(y = x$pointwise[, "p_loo"], y_approx = x$loo_subsampling$elpd_loo_approx)
  x$estimates["p_loo", "Estimate"] <- p_loo_est$y_hat
  x$estimates["p_loo", "SE"] <- sqrt(p_loo_est$hat_v_y)
  x$estimates["p_loo", "subsampling SE"] <- sqrt(p_loo_est$v_y_hat)

  update_psis_loo_ss_estimates(x)
}

#' Simple-re-sample without resampling (SRS-WOR) estimation
#' @noRd
#' @param y The values observed.
#' @param y_approx A vector of length N.
#' @return A list of estimates.
srs_est <- function(y, y_approx) {
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



## Specialized assertions of objects ---

#' Assert that the object has the expected properties
#' @noRd
#' @param x An object to assert.
#' @param N The total number of data points in data.
#' @param estimator The estimator used.
#' @return An asserted object of `x`.
assert_observations <- function(x, N, estimator) {
  checkmate::assert_int(N)
  checkmate::assert_choice(estimator, choices = estimator_choices())
  if (is.null(x)) return(x)
  if (checkmate::test_class(x, "psis_loo_ss")) {
    x <- obs_idx(x)
    checkmate::assert_integer(x, lower = 1, upper = N, any.missing = FALSE)
    return(x)
  }
  x <- as.integer(x)
  if (length(x) > 1) {
    checkmate::assert_integer(x, lower = 1, upper = N, any.missing = FALSE)
    if (estimator %in% "hh_pps") {
      message("Sampling proportional to elpd approximation and with replacement assumed.")
    }
    if (estimator %in% c("diff_srs", "srs")) {
      message("Simple random sampling with replacement assumed.")
    }
  } else {
    checkmate::assert_integer(x, lower = 1, any.missing = FALSE)
  }
  x
}

#' Assert that the object has the expected properties
#' @noRd
#' @inheritParams assert_observations
#' @return An asserted object of `x`.
assert_subsample_idxs <- function(x) {
  checkmate::assert_data_frame(x,
                               types = c("integer", "integer"),
                               any.missing = FALSE,
                               min.rows = 1,
                               col.names = "named")
  checkmate::assert_names(names(x), identical.to = c("idx", "m_i"))
  checkmate::assert_integer(x$idx, lower = 1, any.missing = FALSE, unique = TRUE)
  checkmate::assert_integer(x$m_i, lower = 1, any.missing = FALSE)
  x
}

#' Assert that the object has the expected properties
#' @noRd
#' @inheritParams assert_observations
#' @return An asserted object of `x`.
assert_psis_loo_ss <- function(x) {
  checkmate::assert_class(x, "psis_loo_ss")
  checkmate::assert_names(names(x), must.include = c("estimates", "pointwise", "diagnostics", "psis_object", "loo_subsampling"))
  checkmate::assert_names(rownames(x$estimates), must.include = c("elpd_loo", "p_loo", "looic"))
  checkmate::assert_names(colnames(x$estimates), must.include = c("Estimate", "SE", "subsampling SE"))
  assert_subsampling_pointwise(x$pointwise)
  checkmate::assert_names(names(x$loo_subsampling),
                          must.include = c("elpd_loo_approx",
                                           "loo_approximation", "loo_approximation_draws",
                                           "estimator",
                                           "data_dim", "ndraws"))
  checkmate::assert_numeric(x$loo_subsampling$elpd_loo_approx, any.missing = FALSE, len = x$loo_subsampling$data_dim[1])
  checkmate::assert_choice(x$loo_subsampling$loo_approximation, choices = loo_approximation_choices(api = FALSE))
  checkmate::assert_int(x$loo_subsampling$loo_approximation_draws, null.ok = TRUE)
  checkmate::assert_choice(x$loo_subsampling$estimator, choices = estimator_choices())
  checkmate::assert_integer(x$loo_subsampling$data_dim, any.missing = TRUE, len = 2)
  checkmate::assert_int(x$loo_subsampling$data_dim[1], na.ok = FALSE)
  checkmate::assert_integer(x$loo_subsampling$ndraws, len = 1, any.missing = TRUE)
  x
}

#' Assert that the object has the expected properties
#' @noRd
#' @inheritParams assert_observations
#' @return An asserted object of `x`.
assert_subsampling_pointwise <- function(x) {
  checkmate::assert_matrix(x,
                           any.missing = FALSE,
                           ncols = 8)
  checkmate::assert_names(colnames(x), identical.to = c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic", "influence_pareto_k", "idx", "m_i", "elpd_loo_approx"))
  x
}
