#' Efficient approximate leave-one-out cross-validation (LOO) for posterior
#' approximations
#'
#' @param x A log-likelihood array, matrix, or function.
#'   The **Methods (by class)** section, below, has detailed descriptions of how
#'   to specify the inputs for each method.
#' @param save_psis Should the `"psis"` object created internally by
#'   `loo_approximate_posterior()` be saved in the returned object? See
#'   [loo()] for details.
#' @template cores
#' @inheritParams psis_approximate_posterior
#'
#' @details The `loo_approximate_posterior()` function is an S3 generic and
#'   methods are provided for 3-D pointwise log-likelihood arrays, pointwise
#'   log-likelihood matrices, and log-likelihood functions. The implementation
#'   works for posterior approximations where it is possible to compute the log
#'   density for the posterior approximation.
#'
#' @return The `loo_approximate_posterior()` methods return a named list with
#'   class `c("psis_loo_ap", "psis_loo", "loo")`. It has the same structure
#'   as the objects returned by [loo()] but with the additional slot:
#' \describe{
#'  \item{`posterior_approximation`}{
#'   A list with two vectors, `log_p` and `log_g` of the same length
#'   containing the posterior density and the approximation density
#'   for the individual draws.
#'  }
#' }
#'
#' @seealso [loo()], [psis()], [loo_compare()]
#' @template loo-large-data-references
#'
#' @export loo_approximate_posterior
#' @export loo_approximate_posterior.array loo_approximate_posterior.matrix loo_approximate_posterior.function
#'
loo_approximate_posterior <- function(x, log_p, log_g, ...) {
  UseMethod("loo_approximate_posterior")
}

#' @export
#' @templateVar fn loo_approximate_posterior
#' @template array
loo_approximate_posterior.array <-
  function(x,
           log_p,
           log_g,
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    checkmate::assert_flag(save_psis)
    checkmate::assert_int(cores)
    checkmate::assert_matrix(log_p, mode = "numeric", nrows = dim(x)[1], ncols = dim(x)[2])
    checkmate::assert_matrix(log_g, mode = "numeric", nrows = nrow(log_p), ncols = ncol(log_p))

    ll <- llarray_to_matrix(x)
    log_p <- as.vector(log_p)
    log_g <- as.vector(log_g)
    loo_approximate_posterior.matrix(
      ll,
      log_p = log_p,
      log_g = log_g,
      ...,
      save_psis = save_psis,
      cores = cores
    )
  }

#' @export
#' @templateVar fn loo_approximate_posterior
#' @template matrix
loo_approximate_posterior.matrix <-
  function(x,
           log_p,
           log_g,
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    checkmate::assert_flag(save_psis)
    checkmate::assert_int(cores)
    checkmate::assert_numeric(log_p, len = nrow(x))
    checkmate::assert_null(dim(log_p))
    checkmate::assert_numeric(log_g, len = length(log_p))
    checkmate::assert_null(dim(log_g))

    ap_psis <-
      psis_approximate_posterior(
        log_p = log_p,
        log_g = log_g,
        log_liks = x,
        ...,
        cores = cores,
        save_psis = save_psis
      )
    ap_psis$approximate_posterior <- list(log_p = log_p, log_g = log_g)
    class(ap_psis) <- c("psis_loo_ap", class(ap_psis))
    assert_psis_loo_ap(ap_psis)
    ap_psis
  }

#' @export
#' @templateVar fn loo_approximate_posterior
#' @template function
#' @param data,draws,... For the `loo_approximate_posterior.function()` method,
#'   these are the data, posterior draws, and other arguments to pass to the
#'   log-likelihood function. See the **Methods (by class)** section below for
#'   details on how to specify these arguments.
#'
loo_approximate_posterior.function <-
  function(x,
           ...,
           data = NULL,
           draws = NULL,
           log_p = NULL,
           log_g = NULL,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {

    checkmate::assert_numeric(log_p, len = length(log_g))
    checkmate::assert_numeric(log_g, len = length(log_p))
    cores <- loo_cores(cores)
    stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
    .llfun <- validate_llfun(x)
    N <- dim(data)[1]

    psis_list <- parallel_psis_list(N = N,
                                    .loo_i = .loo_ap_i,
                                    .llfun = .llfun,
                                    data = data,
                                    draws = draws,
                                    r_eff = 1, # r_eff is ignored
                                    save_psis = save_psis,
                                    log_p = log_p,
                                    log_g = log_g,
                                    cores = cores,
                                    ...)

    pointwise <- lapply(psis_list, "[[", "pointwise")
    if (save_psis) {
      psis_object_list <- lapply(psis_list, "[[", "psis_object")
      psis_out <- list2importance_sampling(psis_object_list)
      diagnostics <- psis_out$diagnostics
    } else {
      diagnostics_list <- lapply(psis_list, "[[", "diagnostics")
      diagnostics <- list(
        pareto_k = psis_apply(diagnostics_list, "pareto_k"),
        n_eff = psis_apply(diagnostics_list, "n_eff")
      )
    }

    ap_psis <- importance_sampling_loo_object(
      pointwise = do.call(rbind, pointwise),
      diagnostics = diagnostics,
      dims = c(attr(psis_list[[1]], "S"), N),
      is_method = "psis",
      is_object = if (save_psis) psis_out else NULL
    )

    ap_psis$approximate_posterior <- list(log_p = log_p, log_g = log_g)
    class(ap_psis) <- c("psis_loo_ap", class(ap_psis))
    assert_psis_loo_ap(ap_psis)
    ap_psis
  }

# Function that is passed to the FUN argument of lapply, mclapply, or parLapply
# for the loo_approximate_posterior.function method. The arguments and return
# value are the same as the ones documented for the user-facing loo_i function.
.loo_ap_i <-
  function(i,
           llfun,
           ...,
           data,
           draws,
           log_p,
           log_g,
           r_eff = 1,
           save_psis = FALSE,
           is_method) {

    if (is_method != "psis") stop(is_method, " not implemented for aploo.")
    d_i <- data[i, , drop = FALSE]
    ll_i <- llfun(data_i = d_i, draws = draws, ...)
    if (!is.matrix(ll_i)) {
      ll_i <- as.matrix(ll_i)
    }
    psis_out <- ap_psis(log_ratios = -ll_i, log_p = log_p, log_g = log_g, cores = 1)

    structure(
      list(
        pointwise = pointwise_loo_calcs(ll_i, psis_out),
        diagnostics = psis_out$diagnostics,
        psis_object = if (save_psis) psis_out else NULL
      ),
      S = dim(psis_out)[1],
      N = 1
    )
  }


assert_psis_loo_ap <- function(x) {
  checkmate::assert_class(x, "psis_loo_ap")
  checkmate::assert_names(names(x), must.include = c("estimates", "pointwise", "diagnostics", "psis_object", "approximate_posterior"))
  checkmate::assert_names(names(x$approximate_posterior), must.include = c("log_p", "log_g"))
  checkmate::assert_numeric(x$approximate_posterior$log_p, len = length(x$approximate_posterior$log_g), any.missing = FALSE)
  checkmate::assert_numeric(x$approximate_posterior$log_g, len = length(x$approximate_posterior$log_p), any.missing = FALSE)
}

