#' Efficient approximate leave-one-out cross-validation (LOO)
#'
#' The `loo()` methods for arrays, matrices, and functions compute PSIS-LOO
#' CV, efficient approximate leave-one-out (LOO) cross-validation for Bayesian
#' models using Pareto smoothed importance sampling ([PSIS][psis()]). This is an
#' implementation of the methods described in Vehtari, Gelman, and Gabry (2017a,
#' 2017b).
#'
#' @export loo loo.array loo.matrix loo.function
#' @param x A log-likelihood array, matrix, or function. The **Methods (by class)**
#'   section, below, has detailed descriptions of how to specify the inputs for
#'   each method.
#' @param r_eff Vector of relative effective sample size estimates for the
#'   likelihood (`exp(log_lik)`) of each observation. This is related to
#'   the relative efficiency of estimating the normalizing term in
#'   self-normalizing importance sampling when using posterior draws obtained
#'   with MCMC. If MCMC draws are used and `r_eff` is not provided then
#'   the reported PSIS effective sample sizes and Monte Carlo error estimates
#'   will be over-optimistic. If the posterior draws are independent then
#'   `r_eff=1` and can be omitted. See the [relative_eff()]
#'   helper functions for computing `r_eff`.
#'
#' @param save_psis Should the `"psis"` object created internally by `loo()` be
#'   saved in the returned object? The `loo()` function calls [psis()]
#'   internally but by default discards the (potentially large) `"psis"` object
#'   after using it to compute the LOO-CV summaries. Setting `save_psis=TRUE`
#'   will add a `psis_object` component to the list returned by `loo`. Currently
#'   this is only needed if you plan to use the [E_loo()] function to compute
#'   weighted expectations after running `loo`.
#' @template cores
#' @template is_method
#'
#' @details The `loo()` function is an S3 generic and methods are provided for
#'   3-D pointwise log-likelihood arrays, pointwise log-likelihood matrices, and
#'   log-likelihood functions. The array and matrix methods are the most
#'   convenient, but for models fit to very large datasets the `loo.function()`
#'   method is more memory efficient and may be preferable.
#'
#' @section Defining `loo()` methods in a package: Package developers can define
#'   `loo()` methods for fitted models objects. See the example `loo.stanfit()`
#'   method in the **Examples** section below for an example of defining a
#'   method that calls `loo.array()`. The `loo.stanreg()` method in the
#'   **rstanarm** package is an example of defining a method that calls
#'   `loo.function()`.
#'
#' @return The `loo()` methods return a named list with class
#'   `c("psis_loo", "loo")` and components:
#' \describe{
#'  \item{`estimates`}{
#'   A matrix with two columns (`Estimate`, `SE`) and four rows
#'   (`elpd_loo`, `mcse_elpd_loo`, `p_loo`, `looic`). This
#'   contains point estimates and standard errors of the expected log pointwise
#'   predictive density (`elpd_loo`), the Monte Carlo standard error of
#'   `elpd_loo` (`mcse_elpd_loo`), the effective number of parameters
#'   (`p_loo`) and the LOO information criterion `looic` (which is
#'   just `-2 * elpd_loo`, i.e., converted to deviance scale).
#'  }
#'
#'  \item{`pointwise`}{
#'   A matrix with four columns (and number of rows equal to the number of
#'   observations) containing the pointwise contributions of each of the above
#'   measures (`elpd_loo`, `mcse_elpd_loo`, `p_loo`, `looic`).
#'  }
#'
#'  \item{`diagnostics`}{
#'  A named list containing two vectors:
#'    * `pareto_k`: Estimates of the shape parameter \eqn{k} of the
#'      generalized Pareto fit to the importance ratios for each leave-one-out
#'      distribution. See the [pareto-k-diagnostic] page for details.
#'    * `n_eff`: PSIS effective sample size estimates.
#'  }
#'
#'  \item{`psis_object`}{
#'  This component will be `NULL` unless the `save_psis` argument is set to
#'  `TRUE` when calling `loo()`. In that case `psis_object` will be the object
#'  of class `"psis"` that is created when the `loo()` function calls [psis()]
#'  internally to do the PSIS procedure.
#'  }
#' }
#'
#' @seealso
#'  * The **loo** package [vignettes](https://mc-stan.org/loo/articles/index.html)
#'    for demonstrations.
#'  * [psis()] for the underlying Pareto Smoothed Importance Sampling (PSIS)
#'    procedure used in the LOO-CV approximation.
#'  * [pareto-k-diagnostic] for convenience functions for looking at diagnostics.
#'  * [loo_compare()] for model comparison.
#'
#' @template loo-and-psis-references
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
#' \dontrun{
#' ### Usage with stanfit objects
#' # see ?extract_log_lik
#' log_lik1 <- extract_log_lik(stanfit1, merge_chains = FALSE)
#' rel_n_eff <- relative_eff(exp(log_lik1))
#' loo(log_lik1, r_eff = rel_n_eff, cores = 2)
#' }
#'
#' ### Using log-likelihood function instead of array or matrix
#' set.seed(124)
#'
#' # Simulate data and draw from posterior
#' N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
#' p <- rbeta(1, a0, b0)
#' y <- rbinom(N, size = K, prob = p)
#' a <- a0 + sum(y); b <- b0 + N * K - sum(y)
#' fake_posterior <- as.matrix(rbeta(S, a, b))
#' dim(fake_posterior) # S x 1
#' fake_data <- data.frame(y,K)
#' dim(fake_data) # N x 2
#'
#' llfun <- function(data_i, draws) {
#'   # each time called internally within loo the arguments will be equal to:
#'   # data_i: ith row of fake_data (fake_data[i,, drop=FALSE])
#'   # draws: entire fake_posterior matrix
#'   dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
#' }
#'
#' # Use the loo_i function to check that llfun works on a single observation
#' # before running on all obs. For example, using the 3rd obs in the data:
#' loo_3 <- loo_i(i = 3, llfun = llfun, data = fake_data, draws = fake_posterior, r_eff = NA)
#' print(loo_3$pointwise[, "elpd_loo"])
#'
#' # Use loo.function method (setting r_eff=NA since this posterior not obtained via MCMC)
#' loo_with_fn <- loo(llfun, draws = fake_posterior, data = fake_data, r_eff = NA)
#'
#' # If we look at the elpd_loo contribution from the 3rd obs it should be the
#' # same as what we got above with the loo_i function and i=3:
#' print(loo_with_fn$pointwise[3, "elpd_loo"])
#' print(loo_3$pointwise[, "elpd_loo"])
#'
#' # Check that the loo.matrix method gives same answer as loo.function method
#' log_lik_matrix <- sapply(1:N, function(i) {
#'   llfun(data_i = fake_data[i,, drop=FALSE], draws = fake_posterior)
#' })
#' loo_with_mat <- loo(log_lik_matrix, r_eff = NA)
#' all.equal(loo_with_mat$estimates, loo_with_fn$estimates) # should be TRUE!
#'
#'
#' \dontrun{
#' ### For package developers: defining loo methods
#'
#' # An example of a possible loo method for 'stanfit' objects (rstan package).
#' # A similar method is planned for a future release of rstan (or is already
#' # released, depending on when you are reading this). In order for users
#' # to be able to call loo(stanfit) instead of loo.stanfit(stanfit) the
#' # NAMESPACE needs to be handled appropriately (roxygen2 and devtools packages
#' # are good for that).
#' #
#' loo.stanfit <-
#'  function(x,
#'          pars = "log_lik",
#'          ...,
#'          save_psis = FALSE,
#'          cores = getOption("mc.cores", 1)) {
#'   stopifnot(length(pars) == 1L)
#'   LLarray <- loo::extract_log_lik(stanfit = x,
#'                                   parameter_name = pars,
#'                                   merge_chains = FALSE)
#'   r_eff <- loo::relative_eff(x = exp(LLarray), cores = cores)
#'   loo::loo.array(LLarray,
#'                  r_eff = r_eff,
#'                  cores = cores,
#'                  save_psis = save_psis)
#' }
#' }
#'
#'

loo <- function(x, ...) {
  UseMethod("loo")
}

#' @export
#' @templateVar fn loo
#' @template array
#'
loo.array <-
  function(x,
           ...,
           r_eff = NULL,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           is_method = c("psis", "tis", "sis")) {
    if (is.null(r_eff)) throw_loo_r_eff_warning()
    is_method <- match.arg(is_method)
    psis_out <- importance_sampling.array(log_ratios = -x, r_eff = r_eff, cores = cores, method = is_method)
    ll <- llarray_to_matrix(x)
    pointwise <- pointwise_loo_calcs(ll, psis_out)
    importance_sampling_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out$diagnostics,
      dims = dim(psis_out),
      is_method = is_method,
      is_object = if (save_psis) psis_out else NULL
    )
  }

#' @export
#' @templateVar fn loo
#' @template matrix
#'
loo.matrix <-
  function(x,
           ...,
           r_eff = NULL,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           is_method = c("psis", "tis", "sis")) {
    is_method <- match.arg(is_method)
    if (is.null(r_eff)) {
      throw_loo_r_eff_warning()
    }
    psis_out <-
      importance_sampling.matrix(
        log_ratios = -x,
        r_eff = r_eff,
        cores = cores,
        method = is_method
      )
    pointwise <- pointwise_loo_calcs(x, psis_out)
    importance_sampling_loo_object(
      pointwise = pointwise,
      diagnostics = psis_out$diagnostics,
      dims = dim(psis_out),
      is_method = is_method,
      is_object = if (save_psis) psis_out else NULL
    )
  }

#' @export
#' @templateVar fn loo
#' @template function
#' @param data,draws,... For the `loo.function()` method and the `loo_i()`
#'   function, these are the data, posterior draws, and other arguments to pass
#'   to the log-likelihood function. See the **Methods (by class)** section
#'   below for details on how to specify these arguments.
#'
loo.function <-
  function(x,
           ...,
           data = NULL,
           draws = NULL,
           r_eff = NULL,
           save_psis = FALSE,
           cores = getOption("mc.cores", 1),
           is_method = c("psis", "tis", "sis")) {
    is_method <- match.arg(is_method)
    cores <- loo_cores(cores)
    stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))
    assert_importance_sampling_method_is_implemented(is_method)
    .llfun <- validate_llfun(x)
    N <- dim(data)[1]

    if (is.null(r_eff)) {
      throw_loo_r_eff_warning()
    } else {
      r_eff <- prepare_psis_r_eff(r_eff, len = N)
    }

    psis_list <-
      parallel_importance_sampling_list(
        N = N,
        .loo_i = .loo_i,
        .llfun = .llfun,
        data = data,
        draws = draws,
        r_eff = r_eff,
        save_psis = save_psis,
        cores = cores,
        method = is_method,
        ...
      )

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

    importance_sampling_loo_object(
      pointwise = do.call(rbind, pointwise),
      diagnostics = diagnostics,
      dims = c(attr(psis_list[[1]], "S"), N),
      is_method = is_method,
      is_object = if (save_psis) psis_out else NULL
    )
  }


#' @description The `loo_i()` function enables testing log-likelihood
#'   functions for use with the `loo.function()` method.
#'
#' @rdname loo
#' @export
#'
#' @param i For `loo_i()`, an integer in `1:N`.
#' @param llfun For `loo_i()`, the same as `x` for the
#'   `loo.function()` method. A log-likelihood function as described in the
#'   **Methods (by class)** section.
#'
#' @return The `loo_i()` function returns a named list with components
#'   `pointwise` and `diagnostics`. These components have the same
#'   structure as the `pointwise` and `diagnostics` components of the
#'   object returned by `loo()` except they contain results for only a single
#'   observation.
#'
loo_i <-
  function(i,
           llfun,
           ...,
           data = NULL,
           draws = NULL,
           r_eff = NULL,
           is_method = "psis"
           ) {
    stopifnot(
      i == as.integer(i),
      is.function(llfun) || is.character(llfun),
      is.data.frame(data) || is.matrix(data),
      i <= dim(data)[1],
      !is.null(draws),
      is_method %in% implemented_is_methods()
    )
    .loo_i(
      i = as.integer(i),
      llfun = match.fun(llfun),
      data = data,
      draws = draws,
      r_eff = r_eff[i],
      save_psis = FALSE,
      is_method = is_method,
      ...
    )
  }


# Function that is passed to the FUN argument of lapply, mclapply, or parLapply
# for the loo.function method. The arguments and return value are the same as
# the ones documented above for the user-facing loo_i function.
.loo_i <-
  function(i,
           llfun,
           ...,
           data,
           draws,
           r_eff = NULL,
           save_psis = FALSE,
           is_method) {

    if (!is.null(r_eff)) {
      r_eff <- r_eff[i]
    }
    d_i <- data[i, , drop = FALSE]
    ll_i <- llfun(data_i = d_i, draws = draws, ...)
    if (!is.matrix(ll_i)) {
      ll_i <- as.matrix(ll_i)
    }
    psis_out <-
      importance_sampling.matrix(
        log_ratios = -ll_i,
        r_eff = r_eff,
        cores = 1,
        method = is_method
      )
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


#' @export
dim.loo <- function(x) {
  attr(x, "dims")
}

#' @rdname loo
#' @export
is.loo <- function(x) {
  inherits(x, "loo")
}

#' @export
dim.psis_loo <- function(x) {
  attr(x, "dims")
}

#' @rdname loo
#' @export
is.psis_loo <- function(x) {
  inherits(x, "psis_loo") && is.loo(x)
}


# internal ----------------------------------------------------------------

#' Compute pointwise elpd_loo, p_loo, looic from log lik matrix and
#' psis log weights
#'
#' @noRd
#' @param ll Log-likelihood matrix.
#' @param psis_object The object returned by `psis()`.
#' @return Named list with pointwise elpd_loo, p_loo, and looic.
#'
pointwise_loo_calcs <- function(ll, psis_object) {
  if (!is.matrix(ll)) {
    ll <- as.matrix(ll)
  }
  lw <- weights(psis_object, normalize = TRUE, log = TRUE)
  elpd_loo <- matrixStats::colLogSumExps(ll + lw)
  lpd <- matrixStats::colLogSumExps(ll) - log(nrow(ll)) # colLogMeanExps
  p_loo <- lpd - elpd_loo
  mcse_elpd_loo <- mcse_elpd(ll, lw, E_elpd = elpd_loo, r_eff = relative_eff(psis_object))
  looic <- -2 * elpd_loo
  cbind(elpd_loo, mcse_elpd_loo, p_loo, looic)
}

#' Structure the object returned by the loo methods
#'
#' @noRd
#' @param pointwise Matrix containing columns elpd_loo, mcse_elpd_loo, p_loo,
#'   looic.
#' @param diagnostics Named list containing vector `pareto_k` and vector `n_eff`.
#' @param dims Log likelihood matrix dimensions (attribute of `"psis"` object).
#' @template is_method
#' @param is_object An object of class `"psis"/"tis"/"sis"`, as returned by the [psis()/tis()/sis()] function.
#' @return A `'importance_sampling_loo'` object as described in the Value section of the [loo()]
#'   function documentation.
#'
importance_sampling_loo_object <- function(pointwise, diagnostics, dims,
                                           is_method, is_object = NULL) {
  if (!is.matrix(pointwise)) stop("Internal error ('pointwise' must be a matrix)")
  if (!is.list(diagnostics)) stop("Internal error ('diagnositcs' must be a list)")
  assert_importance_sampling_method_is_implemented(is_method)

  cols_to_summarize <- !(colnames(pointwise) %in% "mcse_elpd_loo")
  estimates <- table_of_estimates(pointwise[, cols_to_summarize, drop=FALSE])

  out <- nlist(estimates, pointwise, diagnostics)
  if (is.null(is_object)) {
    out[paste0(is_method, "_object")] <- list(NULL)
  } else {
    out[[paste0(is_method, "_object")]] <- is_object
  }

  # maintain backwards compatibility
  old_nms <- c("elpd_loo", "p_loo", "looic", "se_elpd_loo", "se_p_loo", "se_looic")
  out <- c(out, setNames(as.list(estimates), old_nms))

  structure(
    out,
    dims = dims,
    class = c(paste0(is_method, "_loo"), "importance_sampling_loo", "loo")
  )
}


#' Compute Monte Carlo standard error for ELPD
#'
#' @noRd
#' @param ll Log-likelihood matrix.
#' @param E_elpd elpd_loo column of pointwise matrix.
#' @param psis_object Object returned by [psis()].
#' @param n_samples Number of draws to take from `Normal(E[epd_i], SD[epd_i])`.
#' @return Vector of standard error estimates.
#'
mcse_elpd <- function(ll, lw, E_elpd, r_eff, n_samples = 1000) {
  lik <- exp(ll)
  w2 <- exp(lw)^2
  E_epd <- exp(E_elpd)
  var_elpd <-
    vapply(
      seq_len(ncol(w2)),
      FUN.VALUE = numeric(1),
      FUN = function(i) {
        var_epd_i <- sum(w2[, i] * (lik[, i] - E_epd[i]) ^ 2)
        sd_epd_i <- sqrt(var_epd_i)
        z <- rnorm(n_samples, mean = E_epd[i], sd = sd_epd_i)
        var(log(z[z > 0]))
      }
    )
  sqrt(var_elpd / r_eff)
}


#' Warning message if r_eff not specified
#' @noRd
throw_loo_r_eff_warning <- function() {
  warning(
    "Relative effective sample sizes ('r_eff' argument) not specified.\n",
    "For models fit with MCMC, the reported PSIS effective sample sizes and \n",
    "MCSE estimates will be over-optimistic.",
    call. = FALSE
  )
}

#' Combine many psis objects into a single psis object
#'
#' @noRd
#' @param objects List of `"psis"` objects, each for a single observation.
#' @return A single `"psis"` object.
#'
list2importance_sampling <- function(objects) {
  log_weights <- sapply(objects, "[[", "log_weights")
  diagnostics <- lapply(objects, "[[", "diagnostics")

  method <- psis_apply(objects, "method", fun = "attr", fun_val = character(1))
  methods <- unique(method)
  if (length(methods) == 1) {
    method <- methods
    classes <- c(methods, "importance_sampling", "list")
  } else {
    classes <- c("importance_sampling", "list")
  }

  structure(
    list(
      log_weights = log_weights,
      diagnostics = list(
        pareto_k = psis_apply(diagnostics, item = "pareto_k"),
        n_eff = psis_apply(diagnostics, item = "n_eff")
      )
    ),
    norm_const_log = psis_apply(objects, "norm_const_log", fun = "attr"),
    tail_len = psis_apply(objects, "tail_len", fun = "attr"),
    r_eff = psis_apply(objects, "r_eff", fun = "attr"),
    dims = dim(log_weights),
    method = method,
    class = classes
  )
}

#' Extractor methods
#'
#' These are only defined in order to deprecate with a warning (rather than
#' remove and break backwards compatibility) the old way of accessing the point
#' estimates in a `"psis_loo"` or `"psis"` object. The new way as of
#' v2.0.0 is to get them from the `"estimates"` component of the object.
#'
#' @name old-extractors
#' @keywords internal
#' @param x,i,exact,name See \link{Extract}.
#'
NULL

#' @rdname old-extractors
#' @keywords internal
#' @export
`[.loo` <- function(x, i) {
  flags <- c("elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic",
            "elpd_waic", "se_elpd_waic", "p_waic", "se_p_waic", "waic", "se_waic")

  if (is.character(i)) {
    needs_warning <- which(flags == i)
    if (length(needs_warning)) {
      warning(
        "Accessing ", flags[needs_warning], " using '[' is deprecated ",
        "and will be removed in a future release. ",
        "Please extract the ", flags[needs_warning],
        " estimate from the 'estimates' component instead.",
        call. = FALSE
      )
    }
  }
  NextMethod()
}

#' @rdname old-extractors
#' @keywords internal
#' @export
`[[.loo` <- function(x, i, exact=TRUE) {
  flags <- c("elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic",
             "elpd_waic", "se_elpd_waic", "p_waic", "se_p_waic", "waic", "se_waic")

  if (is.character(i)) {
    needs_warning <- which(flags == i)
    if (length(needs_warning)) {
      warning(
        "Accessing ", flags[needs_warning], " using '[[' is deprecated ",
        "and will be removed in a future release. ",
        "Please extract the ", flags[needs_warning],
        " estimate from the 'estimates' component instead.",
        call. = FALSE
      )
    }
  }
  NextMethod()
}

#' @rdname old-extractors
#' @keywords internal
#' @export
#'
`$.loo` <- function(x, name) {
  flags <- c("elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo", "looic", "se_looic",
             "elpd_waic", "se_elpd_waic", "p_waic", "se_p_waic", "waic", "se_waic")
  needs_warning <- which(flags == name)
  if (length(needs_warning)) {
    warning(
      "Accessing ", flags[needs_warning], " using '$' is deprecated ",
      "and will be removed in a future release. ",
      "Please extract the ", flags[needs_warning],
      " estimate from the 'estimates' component instead.",
      call. = FALSE
    )
  }
  NextMethod()
}


#' Parallel psis list computations
#'
#' @details Refactored function to handle parallel computations
#' for psis_list
#'
#' @keywords internal
#' @inheritParams loo.function
#' @param .loo_i The function used to compute individual loo contributions.
#' @param .llfun See `llfun` in [loo.function()].
#' @param N The total number of observations (i.e. `nrow(data)`).
#' @param method See `is_method` for [loo()]
#'
parallel_psis_list <- function(N, .loo_i, .llfun,
                               data, draws, r_eff,
                               save_psis, cores,
                               ...){
  parallel_importance_sampling_list(N, .loo_i, .llfun,
                                    data, draws, r_eff,
                                    save_psis, cores,
                                    method = "psis", ...)
}

#' @rdname parallel_psis_list
parallel_importance_sampling_list <- function(N, .loo_i, .llfun,
                                              data, draws, r_eff,
                                              save_psis, cores,
                                              method, ...){
  if (cores == 1) {
    psis_list <-
      lapply(
        X = seq_len(N),
        FUN = .loo_i,
        llfun = .llfun,
        data = data,
        draws = draws,
        r_eff = r_eff,
        save_psis = save_psis,
        is_method = method,
        ...
      )
  } else {
    if (!os_is_windows()) {
      # On Mac or Linux use mclapply() for multiple cores
      psis_list <-
        parallel::mclapply(
          mc.cores = cores,
          X = seq_len(N),
          FUN = .loo_i,
          llfun = .llfun,
          data = data,
          draws = draws,
          r_eff = r_eff,
          save_psis = save_psis,
          is_method = method,
          ...
        )
    } else {
      # On Windows use makePSOCKcluster() and parLapply() for multiple cores
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(parallel::stopCluster(cl))
      psis_list <-
        parallel::parLapply(
          cl = cl,
          X = seq_len(N),
          fun = .loo_i,
          llfun = .llfun,
          data = data,
          draws = draws,
          r_eff = r_eff,
          save_psis = save_psis,
          is_method = method,
          ...
        )
    }
  }
}
