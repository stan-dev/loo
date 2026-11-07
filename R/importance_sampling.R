#' A parent class for different importance sampling methods.
#'
#' @inheritParams psis
#' @param method The importance sampling method to use. The following methods
#'   are implemented:
#' * [`"psis"`][psis]: Pareto-Smoothed Importance Sampling (PSIS). Default method.
#' * [`"tis"`][tis]: Truncated Importance Sampling (TIS) with truncation at
#'   `sqrt(S)`, where `S` is the number of posterior draws.
#' * [`"sis"`][sis]: Standard Importance Sampling (SIS).
#'
importance_sampling <- function(log_ratios, method, ...) {
  UseMethod("importance_sampling")
}

validate_log_ratios <- function(x) {
  validate_ll(x)
  # validate_ll() has already ruled out NA and +Inf, so a column without a
  # finite value is exactly a column whose maximum is -Inf
  if (any(matrixStats::colMaxs(x) == -Inf)) {
    stop("Each column of log ratios must contain at least one finite value.")
  }
  invisible(x)
}


#' @rdname importance_sampling
#' @inheritParams psis
#' @export
importance_sampling.array <-
  function(log_ratios, method,
           ...,
           r_eff = 1,
           parallel = getOption("loo.parallel", TRUE)) {
    stopifnot(length(dim(log_ratios)) == 3)
    assert_importance_sampling_method_is_implemented(method)
    log_ratios <- llarray_to_matrix(log_ratios)
    log_ratios <- validate_log_ratios(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, parallel = parallel, method = method)
  }

#' @rdname importance_sampling
#' @inheritParams psis
#' @export
importance_sampling.matrix <-
  function(log_ratios, method,
           ...,
           r_eff = 1,
           parallel = getOption("loo.parallel", TRUE)) {
    assert_importance_sampling_method_is_implemented(method)
    log_ratios <- validate_log_ratios(log_ratios)
    r_eff <- prepare_psis_r_eff(r_eff, len = ncol(log_ratios))
    do_importance_sampling(log_ratios, r_eff = r_eff, parallel = parallel, method = method)
  }

#' @rdname importance_sampling
#' @inheritParams psis
#' @export
importance_sampling.default <-
  function(log_ratios, method, ..., r_eff = 1) {
    stopifnot(is.null(dim(log_ratios)) || length(dim(log_ratios)) == 1)
    assert_importance_sampling_method_is_implemented(method)
    dim(log_ratios) <- c(length(log_ratios), 1)
    r_eff <- prepare_psis_r_eff(r_eff, len = 1)
    importance_sampling.matrix(log_ratios, r_eff = r_eff, parallel = FALSE, method = method)
  }


#' @export
dim.importance_sampling <- function(x) {
  attr(x, "dims")
}


#' Extract importance sampling weights
#'
#' @export
#' @export weights.importance_sampling
#' @method weights importance_sampling
#' @param object An object returned by [psis()], [tis()], or [sis()].
#' @param log Should the weights be returned on the log scale? Defaults to
#'   `TRUE`.
#' @param normalize Should the weights be normalized? Defaults to `TRUE`.
#' @param ... Ignored.
#'
#' @return The `weights()` method returns an object with the same dimensions as
#'   the `log_weights` component of `object`. The `normalize` and `log`
#'   arguments control whether the returned weights are normalized and whether
#'   or not to return them on the log scale.
#'
#' @examples
#' # See the examples at help("psis")
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

#' Validate selected importance sampling method
#' @noRd
#' @keywords internal
#' @description
#' Currently implemented importance sampling methods
assert_importance_sampling_method_is_implemented <- function(x){
  if (!x %in% implemented_is_methods()) {
    stop("Importance sampling method '",
         x,
         "' is not implemented. Implemented methods: '",
         paste0(implemented_is_methods(), collapse = "', '"),
         "'")
  }
}
implemented_is_methods <- function() c("psis", "tis", "sis")


#' Structure the object returned by the importance_sampling methods
#'
#' @noRd
#' @param unnormalized_log_weights Smoothed and possibly truncated log weights,
#'   but unnormalized.
#' @param pareto_k Vector of GPD k estimates.
#' @param tail_len Vector of tail lengths used to fit GPD.
#' @param r_eff Vector of relative MCMC ESS (n_eff) for `exp(log lik)`
#' @param norm_const_log,n_eff Optional per-observation normalizing constants
#'   and effective sample sizes that were already computed, one per column of
#'   `unnormalized_log_weights`. Computed here if `NULL`.
#' @template is_method
#' @return A list of class `"psis"` with structure described in the main doc at
#'   the top of this file.
#'
importance_sampling_object <-
  function(unnormalized_log_weights,
           pareto_k,
           tail_len,
           r_eff,
           method,
           norm_const_log = NULL,
           n_eff = NULL) {
    stopifnot(is.matrix(unnormalized_log_weights))
    methods <- unique(method)
    stopifnot(all(methods %in% implemented_is_methods()))
    if (length(methods) == 1) {
      method <- methods
      classes <- c(tolower(method), "importance_sampling", "list")
    } else {
      classes <- c("importance_sampling", "list")
    }

    if (is.null(norm_const_log)) {
      norm_const_log <- matrixStats::colLogSumExps(unnormalized_log_weights)
    }
    out <- structure(
      list(
        log_weights = unnormalized_log_weights,
        diagnostics = list(pareto_k = pareto_k, n_eff = NULL, r_eff = r_eff)
      ),
      # attributes
      norm_const_log = norm_const_log,
      tail_len = tail_len,
      r_eff = r_eff,
      dims = dim(unnormalized_log_weights),
      method = method,
      class = classes
    )

    if (is.null(n_eff)) {
      w <- weights(out, normalize = TRUE, log = FALSE)
      n_eff <- psis_n_eff(w, r_eff)
    }
    out$diagnostics[["n_eff"]] <- n_eff
    return(out)
  }

#' Do importance sampling given matrix of log weights
#'
#' @noRd
#' @param lr Matrix of log ratios (`-loglik`)
#' @param r_eff Vector of relative effective sample sizes
#' @param parallel Logical
#' @return A list with class `"psis"` and structure described in the main doc at
#'   the top of this file.
#'
do_importance_sampling <- function(log_ratios, r_eff, parallel, method) {
  assert_importance_sampling_method_is_implemented(method)
  N <- ncol(log_ratios)
  S <- nrow(log_ratios)
  k_threshold <- ps_khat_threshold(S)
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

  lw_list <- loo_map(
    seq_len(N),
    importance_sampling_i,
    is_fun = is_fun,
    log_ratios = slices(log_ratios, parallel),
    tail_len = tail_len,
    r_eff = r_eff,
    parallel = parallel
  )

  log_weights <- psis_apply(lw_list, "log_weights", fun_val = numeric(S))
  pareto_k <- psis_apply(lw_list, "pareto_k")
  throw_pareto_warnings(pareto_k, k_threshold)

  importance_sampling_object(
    unnormalized_log_weights = log_weights,
    pareto_k = pareto_k,
    tail_len = tail_len,
    r_eff = r_eff,
    method = rep(method, length(pareto_k)), # Conform to other attr that exist per obs.
    norm_const_log = psis_apply(lw_list, "norm_const_log"),
    n_eff = psis_apply(lw_list, "n_eff")
  )
}

importance_sampling_i <- function(i, is_fun, log_ratios, tail_len, r_eff) {
  out <- is_fun(log_ratios_i = slice(log_ratios, i), tail_len_i = tail_len[i])
  out$norm_const_log <- matrixStats::logSumExp(out$log_weights)
  out$n_eff <- psis_n_eff(exp(out$log_weights - out$norm_const_log), r_eff[i])
  out
}
