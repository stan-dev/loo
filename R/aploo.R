#' Efficient approximate leave-one-out cross-validation (LOO) for posterior approximations
#'
#'
#' @param x A log-likelihood array, matrix, or function. The **Methods (by class)**
#'   section, below, has detailed descriptions of how to specify the inputs for
#'   each method.
#' @export aploo aploo.array aploo.matrix aploo.function
#'
#' @param save_psis Should the `"psis"` object created internally by `aploo()` be
#'   saved in the returned object? See \link{loo} for details.
#' @template cores
#'
#' @details The `aploo()` function is an S3 generic and methods are provided for
#'   3-D pointwise log-likelihood arrays, pointwise log-likelihood matrices, and
#'   log-likelihood functions. The implementation work for posterior approximations
#'   where it is possible to compute the log density for the posterior approximation.
#'
#' @return The `aploo()` methods return a named list with class
#'   `c("psis_aploo", "psis_loo", "loo")` with the additional slot:
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
#' @template large-data-references
#'
#' @examples
#' ### Array and matrix methods (using example objects included with loo package)
#' # Array method
#' LLarr <- example_loglik_array()
#' rel_n_eff <- relative_eff(exp(LLarr))
#' loo(LLarr, r_eff = rel_n_eff, cores = 2)
#'
#' # Matrix method
#' LLmat <- example_loglik_matrix()
#' rel_n_eff <- relative_eff(exp(LLmat), chain_id = rep(1:2, each = 500))
#' loo(LLmat, r_eff = rel_n_eff, cores = 2)
#'
#'
#'
aploo <- function(x, log_p, log_g, ...) {
  UseMethod("aploo")
}

#' @rdname aploo
#' @export
approximate_posterior_loo <- function(x, log_p, log_g, ...) {
  UseMethod("aploo")
}

#' @export
#' @templateVar fn aploo
#' @template array
#'
aploo.array <-
  function(x,
           log_p,
           log_g,
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    if (!requireNamespace("checkmate", quietly=TRUE)) {
      stop("Please install the 'checkmate' package to use this function.", call. = FALSE)
    }
    checkmate::assert_flag(save_psis)
    checkmate::assert_int(cores)
    checkmate::assert_matrix(log_p, mode = "numeric", nrows = dim(x)[1], ncols = dim(x)[2])
    checkmate::assert_matrix(log_g, mode = "numeric", nrows = nrow(log_p), ncols = ncol(log_p))

    ll <- llarray_to_matrix(x)
    log_p <- as.vector(log_p)
    log_g <- as.vector(log_g)
    aploo.matrix(ll, log_p = log_p, log_g = log_g, ...,  save_psis = save_psis, cores = cores)
  }

#' @export
#' @templateVar fn aploo
#' @template matrix
aploo.matrix <-
  function(x,
           log_p,
           log_g,
           ...,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1)) {
    if (!requireNamespace("checkmate", quietly=TRUE)) {
      stop("Please install the 'checkmate' package to use this function.", call. = FALSE)
    }
    checkmate::assert_flag(save_psis)
    checkmate::assert_int(cores)
    checkmate::assert_numeric(log_p, len = nrow(x))
    checkmate::assert_null(dim(log_p))
    checkmate::assert_numeric(log_g, len = length(log_p))
    checkmate::assert_null(dim(log_g))

    ap_psis <- psis_approximate_posterior(log_p = log_p, log_q = log_g, log_liks = x, ..., cores = cores, save_psis = save_psis)
    class(ap_psis) <- c("psis_aploo", class(ap_psis))
    ap_psis
  }

#' @export
#' @templateVar fn aploo
#' @template function
#' @param data,draws,... For the `aploo.function()` method and the `loo_i()`
#'   function, these are the data, posterior draws, and other arguments to pass
#'   to the log-likelihood function. See the **Methods (by class)** section
#'   below for details on how to specify these arguments.
#'
aploo.function <-
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
                                    .loo_i = .aploo_i,
                                    .llfun = .llfun,
                                    data = data,
                                    draws = draws,
                                    r_eff = NULL, # r_eff is ignored
                                    save_psis = save_psis,
                                    log_p = log_p,
                                    log_g = log_g,
                                    cores = cores,
                                    ...)

    pointwise <- lapply(psis_list, "[[", "pointwise")
    if (save_psis) {
      psis_object_list <- lapply(psis_list, "[[", "psis_object")
      psis_out <- list2psis(psis_object_list)
      diagnostics <- psis_out$diagnostics
    } else {
      diagnostics_list <- lapply(psis_list, "[[", "diagnostics")
      diagnostics <- list(
        pareto_k = psis_apply(diagnostics_list, "pareto_k"),
        n_eff = psis_apply(diagnostics_list, "n_eff")
      )
    }

    psis_loo_object(
      pointwise = do.call(rbind, pointwise),
      diagnostics = diagnostics,
      dims = c(attr(psis_list[[1]], "S"), N),
      psis_object = if (save_psis) psis_out else NULL
    )
  }


.aploo_i <-
  function(i,
           llfun,
           ...,
           data,
           draws,
           log_p,
           log_g,
           r_eff = NULL,
           save_psis = FALSE) {

    if(!is.null(r_eff)) warning("r_eff not implemented for aploo")
    d_i <- data[i, , drop = FALSE]
    ll_i <- llfun(data_i = d_i, draws = draws, ...)
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

