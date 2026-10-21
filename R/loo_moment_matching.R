#' Moment matching for efficient approximate leave-one-out cross-validation (LOO)
#'
#' Moment matching algorithm for updating a loo object when Pareto k estimates
#' are large.
#'
#' @export loo_moment_match loo_moment_match.default
#' @param x A fitted model object.
#' @param loo A loo object to be modified.
#' @param post_draws A function the takes `x` as the first argument and returns
#'   a matrix of posterior draws of the model parameters.
#' @param log_lik_i A function that takes `x` and `i` and returns a matrix (one
#'   column per chain) or a vector (all chains stacked) of log-likelihood draws
#'   of the `i`th observation based on the model `x`. If the draws are obtained
#'   using MCMC, the matrix with MCMC chains separated is preferred.
#' @param unconstrain_pars A function that takes arguments `x`, and `pars` and
#'   returns posterior draws on the unconstrained space based on the posterior
#'   draws on the constrained space passed via `pars`.
#' @param log_prob_upars A function that takes arguments `x` and `upars` and
#'   returns a matrix of log-posterior density values of the unconstrained
#'   posterior draws passed via `upars`.
#' @param log_lik_i_upars A function that takes arguments `x`, `upars`, and `i`
#'   and returns a vector of log-likelihood draws of the `i`th observation based
#'   on the unconstrained posterior draws passed via `upars`.
#' @param max_iters Maximum number of moment matching iterations. Usually this
#'   does not need to be modified. If the maximum number of iterations is
#'   reached, there will be a warning, and increasing `max_iters` may improve
#'   accuracy.
#' @param k_threshold Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is `min(1 - 1/log10(S), 0.7)`,
#'   where `S` is the sample size.
#' @param split Logical; Indicate whether to do the split transformation or not
#'   at the end of moment matching for each LOO fold.
#' @param cov Logical; Indicate whether to match the covariance matrix of the
#'   samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched.
#' @template cores
#' @param ... Further arguments passed to the custom functions documented above.
#'
#' @return The `loo_moment_match()` methods return an updated `loo` object. The
#'   structure of the updated `loo` object is similar, but the method also
#'   stores the original Pareto k diagnostic values in the diagnostics field.
#'
#' @details The `loo_moment_match()` function is an S3 generic and we provide a
#'   default method that takes as arguments user-specified functions
#'   `post_draws`, `log_lik_i`, `unconstrain_pars`, `log_prob_upars`, and
#'   `log_lik_i_upars`. All of these functions should take `...`. as an argument
#'   in addition to those specified for each function.
#'
#' @seealso [loo()], [loo_moment_match_split()]
#' @template moment-matching-references
#'
#' @examples
#' # See the vignette for loo_moment_match()

#' @export
loo_moment_match <- function(x, ...) {
  UseMethod("loo_moment_match")
}


#' @describeIn loo_moment_match A default method that takes as arguments a
#'   user-specified model object `x`, a `loo` object and user-specified
#'   functions `post_draws`, `log_lik_i`, `unconstrain_pars`, `log_prob_upars`,
#'   and `log_lik_i_upars`.
#' @export
loo_moment_match.default <- function(x, loo, post_draws, log_lik_i,
                          unconstrain_pars, log_prob_upars,
                          log_lik_i_upars, max_iters = 30L,
                          k_threshold = NULL, split = TRUE,
                          cov = TRUE, cores = getOption("mc.cores", 1),
                          ...) {

  # input checks
  checkmate::assertClass(loo,classes = "loo")
  checkmate::assertFunction(post_draws)
  checkmate::assertFunction(log_lik_i)
  checkmate::assertFunction(unconstrain_pars)
  checkmate::assertFunction(log_prob_upars)
  checkmate::assertFunction(log_lik_i_upars)
  checkmate::assertNumber(max_iters)
  checkmate::assertNumber(k_threshold, null.ok=TRUE)
  checkmate::assertLogical(split)
  checkmate::assertLogical(cov)
  checkmate::assertNumber(cores)


  if ("psis_loo" %in% class(loo)) {
    is_method <- "psis"
  } else {
    stop("loo_moment_match currently supports only the \"psis\" importance sampling class.")
  }


  S <- dim(loo)[1]
  N <- dim(loo)[2]
  if (is.null(k_threshold)) {
    k_threshold <- ps_khat_threshold(S)
  }
  pars <- post_draws(x, ...)
  # transform the model parameters to unconstrained space
  upars <- unconstrain_pars(x, pars = pars, ...)
  # number of parameters in the **parameters** block only
  npars <- dim(upars)[2]
  # if more parameters than samples, do not do Cholesky transformation
  cov <- cov && S >= 10 * npars
  # compute log-probabilities of the original parameter values
  orig_log_prob <- log_prob_upars(x, upars = upars, ...)

  # loop over all observations whose Pareto k is high
  ks <- loo$diagnostics$pareto_k
  kfs <- rep(0,N)
  I <- which(ks > k_threshold)

  loo_moment_match_i_fun <- function(i) {
    loo_moment_match_i(i = i, x = x, log_lik_i = log_lik_i,
                       unconstrain_pars = unconstrain_pars,
                       log_prob_upars = log_prob_upars,
                       log_lik_i_upars = log_lik_i_upars,
                       max_iters = max_iters, k_threshold = k_threshold,
                       split = split, cov = cov, N = N, S = S, upars = upars,
                       orig_log_prob = orig_log_prob, k = ks[i],
                       is_method = is_method, npars = npars, ...)
  }

  if (cores == 1) {
    mm_list <- lapply(X = I, FUN = function(i) loo_moment_match_i_fun(i))
  }
  else {
    if (!os_is_windows()) {
      mm_list <- parallel::mclapply(X = I, mc.cores = cores,
                                    FUN = function(i) loo_moment_match_i_fun(i))
    }
    else {
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      mm_list <- parallel::parLapply(cl = cl, X = I,
                                    fun = function(i) loo_moment_match_i_fun(i))
    }
  }

  # update results
  failed <- integer(0)
  stopped <- integer(0)
  split_failed <- integer(0)
  error_msg <- NULL
  for (ii in seq_along(I)) {
    i <- I[ii]
    # parallel::mclapply() returns NULL if a worker dies, for example when it
    # runs out of memory, and a "try-error" string if a worker stops with an
    # error. Errors of the user functions are caught in loo_moment_match_i(),
    # so a worker error is an error in loo itself, as it is with cores = 1.
    if (is.null(mm_list[[ii]])) {
      mm_list[[ii]] <- list(i = i, failed = TRUE,
                            error = "a parallel worker did not return a result")
    } else if (inherits(mm_list[[ii]], "try-error")) {
      stop(try_error_message(mm_list[[ii]]), call. = FALSE)
    }
    if (!is.null(mm_list[[ii]]$error)) {
      error_msg <- c(error_msg, mm_list[[ii]]$error)
    }
    if (isTRUE(mm_list[[ii]]$failed)) {
      failed <- c(failed, i)
      next
    }
    if (!is.null(mm_list[[ii]]$error)) {
      stopped <- c(stopped, i)
    }
    if (isTRUE(mm_list[[ii]]$split_failed)) {
      split_failed <- c(split_failed, i)
    }
    loo$pointwise[i, "elpd_loo"] <- mm_list[[ii]]$elpd_loo_i
    loo$pointwise[i, "p_loo"] <- mm_list[[ii]]$p_loo
    loo$pointwise[i, "mcse_elpd_loo"] <- mm_list[[ii]]$mcse_elpd_loo
    loo$pointwise[i, "looic"] <- mm_list[[ii]]$looic

    loo$diagnostics$pareto_k[i] <- mm_list[[ii]]$k
    loo$diagnostics$n_eff[i] <- mm_list[[ii]]$n_eff
    kfs[i] <- mm_list[[ii]]$kf

    if (!is.null(loo$psis_object)) {
      loo$psis_object$log_weights[, i] <- mm_list[[ii]]$lwi
    }
  }
  if (length(I) > 0 && length(failed) == length(I)) {
    stop("Moment matching failed for all observations. The first error ",
         "message:\n", error_msg[1], call. = FALSE)
  }
  if (!is.null(loo$psis_object)) {
    attr(loo$psis_object, "norm_const_log") <- matrixStats::colLogSumExps(loo$psis_object$log_weights)
    loo$psis_object$diagnostics <- loo$diagnostics
  }

  # combined estimates
  cols_to_summarize <- !(colnames(loo$pointwise) %in% c("mcse_elpd_loo",
                                                        "influence_pareto_k"))
  loo$estimates <- table_of_estimates(loo$pointwise[, cols_to_summarize,
                                                    drop = FALSE])

  # these will be deprecated at some point
  loo$elpd_loo <- loo$estimates["elpd_loo","Estimate"]
  loo$p_loo <- loo$estimates["p_loo","Estimate"]
  loo$looic <- loo$estimates["looic","Estimate"]
  loo$se_elpd_loo <- loo$estimates["elpd_loo","SE"]
  loo$se_p_loo <- loo$estimates["p_loo","SE"]
  loo$se_looic <- loo$estimates["looic","SE"]

  throw_moment_match_error_warning(failed, stopped, error_msg)
  # Warn if some Pareto ks are still high
  throw_pareto_warnings(loo$diagnostics$pareto_k, k_threshold)
  # if we don't split, accuracy may be compromised
  if (!split) {
    throw_large_kf_warning(kfs, k_threshold)
  } else if (length(split_failed) > 0) {
    throw_large_kf_warning(kfs[split_failed], k_threshold, split_failed = TRUE)
  }

  loo
}




# Internal functions ---------------


#' Do moment matching for a single observation.
#'
#' @noRd
#' @param i observation number.
#' @param x A fitted model object.
#' @param log_lik_i A function that takes `x` and `i` and returns a matrix (one
#'   column per chain) or a vector (all chains stacked) of log-likelihood draws
#'   of the `i`th observation based on the model `x`. If the draws are obtained
#'   using MCMC, the matrix with MCMC chains separated is preferred.
#' @param unconstrain_pars A function that takes arguments `x`, and `pars` and
#'   returns posterior draws on the unconstrained space based on the posterior
#'   draws on the constrained space passed via `pars`.
#' @param log_prob_upars A function that takes arguments `x` and `upars` and
#'   returns a matrix of log-posterior density values of the unconstrained
#'   posterior draws passed via `upars`.
#' @param log_lik_i_upars A function that takes arguments `x`, `upars`, and `i`
#'   and returns a vector of log-likelihood draws of the `i`th observation based
#'   on the unconstrained posterior draws passed via `upars`.
#' @param max_iters Maximum number of moment matching iterations. Usually this
#'   does not need to be modified. If the maximum number of iterations is
#'   reached, there will be a warning, and increasing `max_iters` may improve
#'   accuracy.
#' @param k_threshold Threshold value for Pareto k values above which the moment
#'   matching algorithm is used. The default value is 0.5.
#' @param split Logical; Indicate whether to do the split transformation or not
#'   at the end of moment matching for each LOO fold.
#' @param cov Logical; Indicate whether to match the covariance matrix of the
#'   samples or not. If `FALSE`, only the mean and marginal variances are
#'   matched.
#' @param N Number of observations.
#' @param S number of MCMC draws.
#' @param upars A matrix representing a sample of vector-valued parameters in
#'   the unconstrained space.
#' @param orig_log_prob log probability densities of the original draws from the
#'   model `x`.
#' @param k Pareto k value before moment matching
#' @template is_method
#' @param npars Number of parameters in the model
#' @param ... Further arguments passed to the custom functions documented above.
#' @return List with the updated elpd values and diagnostics
#'
loo_moment_match_i <- function(i,
                               x,
                               log_lik_i,
                               unconstrain_pars,
                               log_prob_upars,
                               log_lik_i_upars,
                               max_iters,
                               k_threshold,
                               split,
                               cov,
                               N,
                               S,
                               upars,
                               orig_log_prob,
                               k,
                               is_method,
                               npars,
                               ...) {
  # initialize values for this LOO-fold
  uparsi <- upars
  ki <- k
  kfi <- 0
  log_liki <- tryCatch(log_lik_i(x, i, ...), error = function(e) e)
  if (inherits(log_liki, "error")) {
    # the user function failed for this observation, which keeps its PSIS
    # results
    return(list(i = i, failed = TRUE, error = conditionMessage(log_liki)))
  }
  S_per_chain <- NROW(log_liki)
  N_chains <- NCOL(log_liki)
  dim(log_liki) <- c(S_per_chain, N_chains, 1)
  r_eff_i <- loo::relative_eff(exp(log_liki), cores = 1)
  dim(log_liki) <- NULL
  lpd <- matrixStats::logSumExp(log_liki) - log(length(log_liki))

  is_obj <- suppressWarnings(importance_sampling.default(-log_liki,
                                                         method = is_method,
                                                         r_eff = r_eff_i,
                                                         cores = 1))
  lwi <- as.vector(weights(is_obj))
  lwfi <- rep(-matrixStats::logSumExp(rep(0, S)),S)

  # initialize objects that keep track of the total transformation
  total <- list(
    shift = rep(0, npars),
    scaling = rep(1, npars),
    mapping = diag(npars)
  )

  # message of an error that stopped moment matching for this observation
  error_msg <- NULL

  # try several transformations one by one
  # if one does not work, do not apply it and try another one
  # to accept the transformation, Pareto k needs to improve
  # when transformation succeeds, start again from the first one
  iterind <- 1
  while (iterind <= max_iters && ki > k_threshold) {

    if (iterind == max_iters) {
      throw_moment_match_max_iters_warning()
    }

    # 1. match means
    trans <- shift(x, uparsi, lwi)
    # gather updated quantities
    quantities_i <- try(update_quantities_i(x, trans$upars,  i = i,
                                        orig_log_prob = orig_log_prob,
                                        log_prob_upars = log_prob_upars,
                                        log_lik_i_upars = log_lik_i_upars,
                                        r_eff_i = r_eff_i,
                                        cores = 1,
                                        is_method = is_method,
                                        ...),
                        silent = TRUE
                        )
    if (inherits(quantities_i, "try-error")) {
      # log_prob_upars() or log_lik_i_upars() failed, for example because the
      # transformed draws are extreme
      error_msg <- try_error_message(quantities_i)
      break
    }
    if (quantities_i$ki < ki) {
      uparsi <- trans$upars
      total <- add_step_to_total(total, trans, cov)

      lwi <- quantities_i$lwi
      lwfi <- quantities_i$lwfi
      ki <- quantities_i$ki
      kfi <- quantities_i$kfi
      log_liki <- quantities_i$log_liki
      iterind <- iterind + 1
      next
    }

    # 2. match means and marginal variances
    trans <- shift_and_scale(x, uparsi, lwi)
    # gather updated quantities
    quantities_i <- try(update_quantities_i(x, trans$upars,  i = i,
                                        orig_log_prob = orig_log_prob,
                                        log_prob_upars = log_prob_upars,
                                        log_lik_i_upars = log_lik_i_upars,
                                        r_eff_i = r_eff_i,
                                        cores = 1,
                                        is_method = is_method,
                                        ...),
                        silent = TRUE
                        )
    if (inherits(quantities_i, "try-error")) {
      # log_prob_upars() or log_lik_i_upars() failed, for example because the
      # transformed draws are extreme
      error_msg <- try_error_message(quantities_i)
      break
    }
    if (quantities_i$ki < ki) {
      uparsi <- trans$upars
      total <- add_step_to_total(total, trans, cov)

      lwi <- quantities_i$lwi
      lwfi <- quantities_i$lwfi
      ki <- quantities_i$ki
      kfi <- quantities_i$kfi
      log_liki <- quantities_i$log_liki
      iterind <- iterind + 1
      next
    }

    # 3. match means and covariances
    if (cov) {
      trans <- shift_and_cov(x, uparsi, lwi)
      # gather updated quantities
      quantities_i <- try(update_quantities_i(x, trans$upars,  i = i,
                                          orig_log_prob = orig_log_prob,
                                          log_prob_upars = log_prob_upars,
                                          log_lik_i_upars = log_lik_i_upars,
                                          r_eff_i = r_eff_i,
                                          cores = 1,
                                          is_method = is_method,
                                          ...),
                          silent = TRUE
                          )

      if (inherits(quantities_i, "try-error")) {
        # log_prob_upars() or log_lik_i_upars() failed, for example because the
        # transformed draws are extreme
        error_msg <- try_error_message(quantities_i)
        break
      }
      if (quantities_i$ki < ki) {
        uparsi <- trans$upars
        total <- add_step_to_total(total, trans, cov)

        lwi <- quantities_i$lwi
        lwfi <- quantities_i$lwfi
        ki <- quantities_i$ki
        kfi <- quantities_i$kfi
        log_liki <- quantities_i$log_liki
        iterind <- iterind + 1
        next
      }
    }
    # none of the transformations improved khat
    # so there is no need to try further
    break
  }

  # transformations are now done
  # if we don't do split transform, or
  # if no transformations were successful
  # stop and collect values
  split_done <- FALSE
  if (split && (iterind > 1)) {
    # compute split transformation
    split_obj <- try(loo_moment_match_split(
      x, upars, cov, total$shift, total$scaling, total$mapping, i,
      log_prob_upars = log_prob_upars, log_lik_i_upars = log_lik_i_upars,
      cores = 1, r_eff_i = r_eff_i, is_method = is_method, ...
    ), silent = TRUE)
    if (inherits(split_obj, "try-error")) {
      # use the result without the split transformation
      error_msg <- try_error_message(split_obj)
    } else {
      log_liki <- split_obj$log_liki
      lwi <- split_obj$lwi
      lwfi <- split_obj$lwfi
      r_eff_i <- split_obj$r_eff_i
      split_done <- TRUE
    }
  }
  if (!split_done) {
    dim(log_liki) <- c(S_per_chain, N_chains, 1)
    r_eff_i <- loo::relative_eff(exp(log_liki), cores = 1)
    dim(log_liki) <- NULL
  }

  # pointwise estimates
  elpd_loo_i <- matrixStats::logSumExp(log_liki + lwi)
  mcse_elpd_loo <- mcse_elpd(
    ll = as.matrix(log_liki), lw = as.matrix(lwi),
    E_elpd = elpd_loo_i, r_eff = r_eff_i
  )

  list(elpd_loo_i = elpd_loo_i,
       p_loo = lpd - elpd_loo_i,
       mcse_elpd_loo = mcse_elpd_loo,
       looic = -2 * elpd_loo_i,
       k = ki,
       kf = kfi,
       n_eff = min(1.0 / sum(exp(2 * lwi)),
                   1.0 / sum(exp(2 * lwfi))) * r_eff_i,
       lwi = lwi,
       i = i,
       error = error_msg,
       # an error before any transformation was accepted leaves the PSIS
       # results
       failed = !is.null(error_msg) && iterind == 1,
       split_failed = split && (iterind > 1) && !split_done)
}




#' Update the importance weights, Pareto diagnostic and log-likelihood
#' for observation `i` based on model `x`.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in
#' the unconstrained space.
#' @param i observation number.
#' @param orig_log_prob log probability densities of the original draws from
#' the model `x`.
#' @param log_prob_upars A function that takes arguments `x` and
#'   `upars` and returns a matrix of log-posterior density values of the
#'   unconstrained posterior draws passed via `upars`.
#' @param log_lik_i_upars A function that takes arguments `x`, `upars`,
#'   and `i` and returns a vector of log-likelihood draws of the `i`th
#'   observation based on the unconstrained posterior draws passed via
#'   `upars`.
#' @param r_eff_i MCMC effective sample size divided by the total sample size
#' for 1/exp(log_ratios) for observation i.
#' @template is_method
#' @return List with the updated importance weights, Pareto diagnostics and
#' log-likelihood values.
#'
update_quantities_i <- function(x, upars, i, orig_log_prob,
                                log_prob_upars, log_lik_i_upars,
                                r_eff_i, is_method, ...) {
  log_prob_new <- log_prob_upars(x, upars = upars, ...)
  log_liki_new <- log_lik_i_upars(x, upars = upars, i = i, ...)
  # compute new log importance weights

  # If log_liki_new and log_prob_new both have same element as Inf,
  # replace the log ratio with -Inf
  lr <- -log_liki_new + log_prob_new - orig_log_prob
  lr[is.na(lr)] <- -Inf
  is_obj_new <- suppressWarnings(importance_sampling.default(lr,
                                                             method = is_method,
                                                             r_eff = r_eff_i,
                                                             cores = 1))
  lwi_new <- as.vector(weights(is_obj_new))
  ki_new <- is_obj_new$diagnostics$pareto_k

  is_obj_f_new <- suppressWarnings(importance_sampling.default(log_prob_new - orig_log_prob,
                                                               method = is_method,
                                                               r_eff = r_eff_i,
                                                               cores = 1))
  lwfi_new <- as.vector(weights(is_obj_f_new))
  kfi_new <- is_obj_f_new$diagnostics$pareto_k

  # gather results
  list(
    lwi = lwi_new,
    lwfi = lwfi_new,
    ki = ki_new,
    kfi = kfi_new,
    log_liki = log_liki_new
  )
}



#' Add one accepted moment matching step to the total transformation
#'
#' Each step maps the current draws `u` to `L (u - mean(u)) + mean(u) + shift`,
#' where `L` is the identity, `diag(scaling)` or `mapping`. The shifts add up,
#' and the linear parts compose in the order in which the steps were accepted.
#' After any accepted step the loop in `loo_moment_match_i()` starts again with
#' the shift step, so a scaling step can follow a covariance step. With
#' `cov = TRUE` the scaling is therefore applied to the rows of `mapping`
#' (`diag(scaling) %*% mapping`), and `scaling` stays one. With `cov = FALSE`
#' there is no covariance step, and the scalings are kept in `scaling`.
#'
#' @noRd
#' @param total A list with `shift`, `scaling` and `mapping`.
#' @param trans The output of `shift()`, `shift_and_scale()` or
#'   `shift_and_cov()`.
#' @param cov Logical; the `cov` argument of the moment matching.
#' @return The updated `total`.
#'
add_step_to_total <- function(total, trans, cov) {
  total$shift <- total$shift + trans$shift
  if (!is.null(trans$scaling)) {
    if (cov) {
      total$mapping <- trans$scaling * total$mapping
    } else {
      total$scaling <- total$scaling * trans$scaling
    }
  }
  if (!is.null(trans$mapping)) {
    total$mapping <- trans$mapping %*% total$mapping
  }
  total
}

#' Shift a matrix of parameters to their weighted mean.
#' Also calls update_quantities_i which updates the importance weights based on
#' the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in
#' the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift that was performed, and the new parameter matrix.
#'
shift <- function(x, upars, lwi) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  # transform posterior draws
  upars_new <- sweep(upars, 2, shift, "+")
  list(
    upars = upars_new,
    shift = shift
  )
}




#' Shift a matrix of parameters to their weighted mean and scale the marginal
#' variances to match the weighted marginal variances. Also calls
#' update_quantities_i which updates the importance weights based on
#' the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in
#' the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift and scaling that were performed, and the new
#' parameter matrix.
#'
#'
shift_and_scale <- function(x, upars, lwi) {
  # compute moments using log weights
  S <- dim(upars)[1]
  mean_original <- colMeans(upars)
  weights <- exp(lwi)
  mean_weighted <- colSums(weights * upars)
  shift <- mean_weighted - mean_original
  # The two-pass form avoids the cancellation in E[x^2] - E[x]^2 and is
  # equivalent to it only because `weights` sums to one.
  centered <- sweep(upars, 2, mean_weighted)
  mii <- colSums(weights * centered^2)
  mii <- mii * S / (S - 1)
  scaling <- sqrt(mii / matrixStats::colVars(upars))
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- sweep(upars_new, 2, scaling, "*")
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")

  list(
    upars = upars_new,
    shift = shift,
    scaling = scaling
  )
}

#' Shift a matrix of parameters to their weighted mean and scale the covariance
#' to match the weighted covariance.
#' Also calls update_quantities_i which updates the importance weights based on
#' the supplied model object.
#'
#' @noRd
#' @param x A fitted model object.
#' @param upars A matrix representing a sample of vector-valued parameters in
#' the unconstrained space
#' @param lwi A vector representing the log-weight of each parameter
#' @return List with the shift and mapping that were performed, and the new
#' parameter matrix.
#'
shift_and_cov <- function(x, upars, lwi, ...) {
  # compute moments using log weights
  mean_original <- colMeans(upars)
  mean_weighted <- colSums(exp(lwi) * upars)
  shift <- mean_weighted - mean_original
  covv <- stats::cov(upars)
  wcovv <- stats::cov.wt(upars, wt = exp(lwi))$cov
  chol1 <- tryCatch(
    {
      chol(wcovv)
    },
    error = function(cond)
    {
      return(NULL)
    }
  )
  if (is.null(chol1)) {
    mapping <- diag(length(mean_original))
  }
  else {
    chol2 <- chol(covv)
    mapping <- t(chol1) %*% solve(t(chol2))
  }
  # transform posterior draws
  upars_new <- sweep(upars, 2, mean_original, "-")
  upars_new <- tcrossprod(upars_new, mapping)
  upars_new <- sweep(upars_new, 2, mean_weighted, "+")
  colnames(upars_new) <- colnames(upars)

  list(
    upars = upars_new,
    shift = shift,
    mapping = mapping
  )
}


#' Message of a "try-error" object without the call and the newline
#' @noRd
try_error_message <- function(x) {
  cond <- attr(x, "condition")
  if (!is.null(cond)) conditionMessage(cond) else trimws(as.character(x))
}

#' Warning message if a function call failed during moment matching
#' @noRd
#' @param failed Observations for which moment matching failed completely.
#' @param stopped Observations for which moment matching stopped early.
#' @param error_msg Error messages, the first one is shown.
throw_moment_match_error_warning <- function(failed, stopped, error_msg) {
  if (length(failed) == 0 && length(stopped) == 0) {
    return(invisible(NULL))
  }
  msg <- NULL
  if (length(stopped) > 0) {
    msg <- c(msg, paste0(
      "Moment matching stopped early for observation(s) ",
      paste(stopped, collapse = ", "), ". The results for these ",
      "observations use the last transformation that succeeded.\n"
    ))
  }
  if (length(failed) > 0) {
    msg <- c(msg, paste0(
      "Moment matching failed for observation(s) ",
      paste(failed, collapse = ", "), ". The results for these ",
      "observations are the PSIS results without moment matching.\n"
    ))
  }
  warning(
    msg,
    "The cause was an error in 'log_prob_upars', 'log_lik_i_upars' or ",
    "'log_lik_i', or a parallel worker that did not return a result. The ",
    "first error message:\n", error_msg[1],
    call. = FALSE
  )
}

#' Warning message if max_iters is reached
#' @noRd
throw_moment_match_max_iters_warning <- function() {
  warning(
    "The maximum number of moment matching iterations ('max_iters' argument)
    was reached.\n",
    "Increasing the value may improve accuracy.",
    call. = FALSE
  )
}

#' Warning message if not using split transformation and accuracy is
#' compromised
#' @noRd
throw_large_kf_warning <- function(kf, k_threshold, split_failed = FALSE) {
  if (any(kf > k_threshold)) {
    if (split_failed) {
      advice <- paste0(
        "The split transformation failed for these observations, ",
        "see the warning about the error."
      )
    } else {
      advice <- "Setting the argument 'split' to 'TRUE' will likely improve accuracy."
    }
    warning(
      "The accuracy of self-normalized importance sampling may be bad.\n",
      advice,
      call. = FALSE
    )
  }
}

