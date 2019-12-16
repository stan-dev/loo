#' A parent class for different importance sampling methods.
#' @keywords internal
#' @inheritParams psis
#' @param method The importance sampling method to use. The following methods
#'   are implemented:
#' * [`"psis"`][psis]: Pareto-Smoothed Importance Sampling (PSIS). Default method.
#' * [`"tis"`][tis]: Truncated Importance Sampling (TIS) with truncation at
#'   `sqrt(S)`, where `S` is the number of posterior draws.
#' * [`"sis"`][sis]: Standard Importance Sampling (SIS).
#'
importance_sampling <- function(log_ratios, method, ...) UseMethod("importance_sampling")

#' @noRd
#' @keywords internal
#' @description
#' Currently implemented importance sampling methods
assert_importance_sampling_method_is_implemented <- function(x){
  if (!x %in% implemented_is_methods()) {
    stop("Importance sampling method '",
         x,
         "' is not implemented. Implemented methods: '",
         paste0(implemented_is_methods, collapse = "', '"),
         "'")
  }
}
implemented_is_methods <- function() c("psis", "tis", "sis")

#' Importance sampling of array
#' @inheritParams psis
#' @template is_method
#' @keywords internal
#' @export
importance_sampling.array <-
  function(log_ratios, method,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    stopifnot(length(dim(log_ratios)) == 3)
    assert_importance_sampling_method_is_implemented(method)
    log_ratios <- validate_ll(log_ratios)
    log_ratios <- llarray_to_matrix(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, cores = cores, method = method)
  }

#' Importance sampling of matrices
#' @inheritParams psis
#' @template is_method
#' @keywords internal
#' @export
importance_sampling.matrix <-
  function(log_ratios, method,
           ...,
           r_eff = NULL,
           cores = getOption("mc.cores", 1)) {
    cores <- loo_cores(cores)
    assert_importance_sampling_method_is_implemented(method)
    log_ratios <- validate_ll(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, cores = cores, method = method)
  }

#' Importance sampling (default)
#' @inheritParams psis
#' @template is_method
#' @keywords internal
#' @export
importance_sampling.default <-
  function(log_ratios, method, ..., r_eff = NULL) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1)
    assert_importance_sampling_method_is_implemented(method)
    dim(log_ratios) <- c(length(log_ratios), 1)
    r_eff <- prepare_psis_r_eff(r_eff, len = 1)
    importance_sampling.matrix(log_ratios, r_eff = r_eff, cores = 1, method = method)
  }


#' @export
dim.importance_sampling <- function(x) {
  attr(x, "dims")
}

#' @rdname psis
#' @export
#' @method weights importance_sampling
#' @param object For the `weights()` method, an object returned by `psis()/tis()/sis()`
#' (a list with same class name `"psis"/"tis"/"sis"`.).
#' @param log For the `weights()` method, should the weights be returned on
#'   the log scale? Defaults to `TRUE`.
#' @param normalize For the `weights()` method, should the weights be
#'   normalized? Defaults to `TRUE`.
#'
#' @return The `weights()` method returns an object with the same dimensions as
#'   the `log_weights` component of the `"psis"/"tis"/"sis"` object. The `normalize` and
#'   `log` arguments control whether the returned weights are normalized and
#'   whether or not to return them on the log scale.
#' @export weights.importance_sampling
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
#' @template is_method
#' @return A list of class `"psis"` with structure described in the main doc at
#'   the top of this file.
#'
importance_sampling_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff,
           method) {
    stopifnot(is.matrix(unnormalized_log_weights))
    methods <- unique(method)
    stopifnot(all(methods %in% implemented_is_methods()))
    if (length(methods) == 1) {
      method <- methods
      classes <- c(tolower(method), "importance_sampling", "list")
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
      method = method,
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
do_importance_sampling <- function(log_ratios, r_eff, cores, method) {
  stopifnot(cores == as.integer(cores))
  assert_importance_sampling_method_is_implemented(method)
  N <- ncol(log_ratios)
  S <- nrow(log_ratios)
  tail_len <- n_pareto(r_eff, S)

  if (method == "psis") {
    is_fun <- do_psis_i
    throw_tail_length_warnings(tail_len)
  } else if (method == "tis") {
    is_fun <- do_tis_i
  } else if (method == "sis") {
    is_fun <- do_sis_i
  } else {
    stop("Incorrect IS method.")
  }

  if (cores == 1) {
    lw_list <- lapply(seq_len(N), function(i)
      is_fun(log_ratios_i = log_ratios[, i], tail_len_i = tail_len[i]))
  } else {
    if (!os_is_windows()) {
      lw_list <- parallel::mclapply(
        X = seq_len(N),
        mc.cores = cores,
        FUN = function(i)
          is_fun(log_ratios_i = log_ratios[, i], tail_len_i = tail_len[i])
      )
    } else {
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      lw_list <-
        parallel::parLapply(
          cl = cl,
          X = seq_len(N),
          fun = function(i)
            is_fun(log_ratios_i = log_ratios[, i], tail_len_i = tail_len[i])
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
    method = rep(method, length(pareto_k)) # Conform to other attr that exist per obs.
  )
}
