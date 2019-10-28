
#' @rdname psis
importance_sampling <- function(log_ratios, ...) UseMethod("importance_sampling")

#' @noRd
#' @description
#' Currently implemented importance sampling methods
implemented_is_methods <- function() c("psis", "tis", "sis")

#' @rdname psis
importance_sampling.array <-
  function(log_ratios, ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1),
           is_method) {
    cores <- loo_cores(cores)
    stopifnot(length(dim(log_ratios)) == 3)
    stopifnot(is_method %in% implemented_is_methods())
    log_ratios <- validate_ll(log_ratios)
    log_ratios <- llarray_to_matrix(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, cores = cores, is_method = is_method)
  }

#' @rdname psis
importance_sampling.matrix <-
  function(log_ratios,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1),
           is_method) {
    cores <- loo_cores(cores)
    stopifnot(is_method %in% implemented_is_methods())
    log_ratios <- validate_ll(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, cores = cores, is_method = is_method)
  }

#' @rdname psis
importance_sampling.default <-
  function(log_ratios, ..., r_eff = NULL,
           is_method) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1)
    stopifnot(is_method %in% implemented_is_methods())
    dim(log_ratios) <- c(length(log_ratios), 1)
    r_eff <- prepare_psis_r_eff(r_eff, len = 1)
    importance_sampling.matrix(log_ratios, r_eff = r_eff, cores = 1, is_method = is_method)
  }


#' @export
dim.importance_sampling <- function(x) {
  attr(x, "dims")
}

#' @rdname psis
#' @export
#' @export weights.importance_sampling
#' @method weights importance_sampling
#' @param object For the `weights()` method, an object returned by `psis()` (a
#'   list with class `"psis"`).
#' @param log For the `weights()` method, should the weights be returned on
#'   the log scale? Defaults to `TRUE`.
#' @param normalize For the `weights()` method, should the weights be
#'   normalized? Defaults to `TRUE`.
#'
#' @return The `weights()` method returns an object with the same dimensions as
#'   the `log_weights` component of the `"psis"` object. The `normalize` and
#'   `log` arguments control whether the returned weights are normalized and
#'   whether or not to return them on the log scale.
#'
weights.importance_sampling <-
  function(object,
           ...,
           log = TRUE,
           normalize = TRUE) {
    out <- object[["log_weights"]] # smoothed but unnormalized log weights
    if (normalize) {
      out <- sweep(out,
                   MARGIN = 2,
                   STATS = attr(object, "norm_const_log"), # colLogSumExp(log_weights)
                   check.margin = FALSE)
    }
    if (!log) {
      out <- exp(out)
    }

    return(out)
  }

# internal ----------------------------------------------------------------

#' Structure the object returned by the importance_sampling methods
#'
#' @noRd
#' @param unnormalized_log_weights Smoothed and possibly truncated log weights,
#'   but unnormalized.
#' @param pareto_k Vector of GPD k estimates.
#' @param tail_len Vector of tail lengths used to fit GPD.
#' @param r_eff Vector of relative MCMC n_eff for `exp(log lik)`
#' @param is_method
#'   Importance sampling method to use. The following approaches are implemented:
#'   * Pareto-Smoothed Importance Sampling (\code{"psis"}).
#'   * Truncated Importance Sampling with truncation at \code{sqrt(S)} (\code{"tis"}).
#'   * Standard Importance Sampling  (\code{"sis"}).
#'   Defaults to \code{"psis"}.
#' @return A list of class `"psis"` with structure described in the main doc at
#'   the top of this file.
#'
importance_sampling_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff,
           is_method) {
    stopifnot(is.matrix(unnormalized_log_weights))
    is_methods <- unique(is_method)
    stopifnot(all(is_methods %in% implemented_is_methods()))
    if(length(is_methods) == 1) {
      is_method <- is_methods
      classes <- c(tolower(is_methods), "importance_sampling", "list")
    } else {
      classes <- c("importance_sampling", "list")
    }

    norm_const_log <- matrixStats::colLogSumExps(unnormalized_log_weights)
    out <- structure(
      list(
        log_weights = unnormalized_log_weights,
        diagnostics = list(pareto_k = pareto_k, n_eff = NULL)
      ),
      # attributes
      norm_const_log = norm_const_log,
      tail_len = tail_len,
      r_eff = r_eff,
      dims = dim(unnormalized_log_weights),
      is_method = is_method,
      class = classes
    )

    # need normalized weights (not on log scale) for psis_n_eff
    w <- weights(out, normalize = TRUE, log = FALSE)
    out$diagnostics[["n_eff"]] <- psis_n_eff(w, r_eff)
    return(out)
  }

#' Do importance sampling given matrix of log weights
#'
#' @noRd
#' @param lr Matrix of log ratios (`-loglik`)
#' @param r_eff Vector of relative effective sample sizes
#' @param cores User's integer `cores` argument
#' @return A list with class `"psis"` and structure described in the main doc at
#'   the top of this file.
#'
do_importance_sampling <- function(log_ratios, r_eff, cores, is_method) {
  stopifnot(cores == as.integer(cores))
  stopifnot(is_method %in% implemented_is_methods())
  N <- ncol(log_ratios)
  S <- nrow(log_ratios)
  tail_len <- n_pareto(r_eff, S)

  if(is_method == "psis") {
    is_fun <- do_psis_i
    throw_tail_length_warnings(tail_len)
  } else if(is_method == "tis") {
    is_fun <- do_tis_i
  } else if(is_method == "sis") {
    is_fun <- do_sis_i
  } else {
    stop("Incorrect IS method.")
  }

  if (cores == 1) {
    lw_list <- lapply(seq_len(N), function(i)
      is_fun(log_ratios[, i], tail_len[i]))
  } else {
    if (.Platform$OS.type != "windows") {
      lw_list <- parallel::mclapply(
        X = seq_len(N),
        mc.cores = cores,
        FUN = function(i)
          is_fun(log_ratios[, i], tail_len[i])
      )
    } else {
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      lw_list <-
        parallel::parLapply(
          cl = cl,
          X = seq_len(N),
          fun = function(i)
            is_fun(log_ratios[, i], tail_len[i])
        )
    }
  }

  log_weights <- psis_apply(lw_list, "log_weights", fun_val = numeric(S))
  pareto_k <- psis_apply(lw_list, "pareto_k")
  throw_pareto_warnings(pareto_k)

  importance_sampling_object(
    unnormalized_log_weights = log_weights,
    pareto_k = pareto_k,
    tail_len = tail_len,
    r_eff = r_eff,
    is_method = rep(is_method, length(pareto_k)) # Conform to other attr that exist per obs.
  )
}
