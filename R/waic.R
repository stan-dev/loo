#' Widely applicable information criterion (WAIC)
#'
#' The `waic()` methods can be used to compute WAIC from the pointwise
#' log-likelihood. However, we recommend LOO-CV using PSIS (as implemented by
#' the [loo()] function) because PSIS provides useful diagnostics as well as
#' effective sample size and Monte Carlo estimates.
#'
#' @export waic waic.array waic.matrix waic.function
#' @inheritParams loo
#'
#' @return A named list (of class `c("waic", "loo")`) with components:
#'
#' \describe{
#'  \item{`estimates`}{
#'  A matrix with two columns (`"Estimate"`, `"SE"`) and three
#'  rows (`"elpd_waic"`, `"p_waic"`, `"waic"`). This contains
#'  point estimates and standard errors of the expected log pointwise predictive
#'  density (`elpd_waic`), the effective number of parameters
#'  (`p_waic`) and the information criterion `waic` (which is just
#'  `-2 * elpd_waic`, i.e., converted to deviance scale).
#'  }
#'  \item{`pointwise`}{
#'  A matrix with three columns (and number of rows equal to the number of
#'  observations) containing the pointwise contributions of each of the above
#'  measures (`elpd_waic`, `p_waic`, `waic`).
#'  }
#' }
#'
#' @seealso
#' * The __loo__ package [vignettes](https://mc-stan.org/loo/articles/) and
#'   Vehtari, Gelman, and Gabry (2017) and Vehtari, Simpson, Gelman, Yao,
#'   and Gabry (2019) for more details on why we prefer `loo()` to `waic()`.
#' * [loo_compare()] for comparing models on approximate LOO-CV or WAIC.
#'
#' @references
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
#' widely application information criterion in singular learning theory.
#' *Journal of Machine Learning Research* **11**, 3571-3594.
#'
#' @template loo-and-psis-references
#'
#' @examples
#' ### Array and matrix methods
#' LLarr <- example_loglik_array()
#' dim(LLarr)
#'
#' LLmat <- example_loglik_matrix()
#' dim(LLmat)
#'
#' waic_arr <- waic(LLarr)
#' waic_mat <- waic(LLmat)
#' identical(waic_arr, waic_mat)
#'
#'
#' \dontrun{
#' log_lik1 <- extract_log_lik(stanfit1)
#' log_lik2 <- extract_log_lik(stanfit2)
#' (waic1 <- waic(log_lik1))
#' (waic2 <- waic(log_lik2))
#' print(compare(waic1, waic2), digits = 2)
#' }
#'
waic <- function(x, ...) {
  UseMethod("waic")
}

#' @export
#' @templateVar fn waic
#' @template array
#'
waic.array <- function(x, ...) {
  waic.matrix(llarray_to_matrix(x), ...)
}

#' @export
#' @templateVar fn waic
#' @template matrix
#'
waic.matrix <- function(x, ...) {
  ll <- validate_ll(x)
  lldim <- dim(ll)
  lpd <- matrixStats::colLogSumExps(ll) - log(nrow(ll)) # colLogMeanExps
  p_waic <- matrixStats::colVars(ll)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  pointwise <- cbind(elpd_waic, p_waic, waic)

  throw_pwaic_warnings(pointwise[, "p_waic"], digits = 1)
  return(waic_object(pointwise, dims = lldim))
}


#' @export
#' @templateVar fn waic
#' @template function
#' @param draws,data,... For the function method only. See the
#' **Methods (by class)** section below for details on these arguments.
#'
waic.function <-
  function(x,
           ...,
           data = NULL,
           draws = NULL) {
    stopifnot(is.data.frame(data) || is.matrix(data), !is.null(draws))

    .llfun <- validate_llfun(x)
    N <- dim(data)[1]
    S <- length(as.vector(.llfun(data_i = data[1,, drop=FALSE], draws = draws, ...)))
    waic_list <- lapply(seq_len(N), FUN = function(i) {
      ll_i <- .llfun(data_i = data[i,, drop=FALSE], draws = draws, ...)
      ll_i <- as.vector(ll_i)
      lpd_i <- logMeanExp(ll_i)
      p_waic_i <- var(ll_i)
      elpd_waic_i <- lpd_i - p_waic_i
      c(elpd_waic = elpd_waic_i, p_waic = p_waic_i)
    })
    pointwise <- do.call(rbind, waic_list)
    pointwise <- cbind(pointwise, waic = -2 * pointwise[, "elpd_waic"])

    throw_pwaic_warnings(pointwise[, "p_waic"], digits = 1)
    waic_object(pointwise, dims = c(S, N))
  }


#' @export
dim.waic <- function(x) {
  attr(x, "dims")
}

#' @rdname waic
#' @export
is.waic <- function(x) {
  inherits(x, "waic") && is.loo(x)
}


# internal ----------------------------------------------------------------

# structure the object returned by the waic methods
waic_object <- function(pointwise, dims) {
  estimates <- table_of_estimates(pointwise)
  out <- nlist(estimates, pointwise)
  # maintain backwards compatibility
  old_nms <- c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")
  out <- c(out, setNames(as.list(estimates), old_nms))
  structure(
    out,
    dims = dims,
    class = c("waic", "loo")
  )
}

# waic warnings
# @param p 'p_waic' estimates
throw_pwaic_warnings <- function(p, digits = 1, warn = TRUE) {
  badp <- p > 0.4
  if (any(badp)) {
    count <- sum(badp)
    prop <- count / length(badp)
    msg <- paste0("\n", count, " (", .fr(100 * prop, digits),
                  "%) p_waic estimates greater than 0.4. ",
                  "We recommend trying loo instead.")
    if (warn) .warn(msg) else cat(msg, "\n")
  }
  invisible(NULL)
}

