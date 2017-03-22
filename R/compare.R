#' Model comparison
#'
#' Compare fitted models on LOO or WAIC
#'
#' @export
#' @param ... At least two objects returned by \code{\link{loo}} (or
#'   \code{\link{waic}}).
#' @param x A list of at least two objects returned by \code{\link{loo}} (or
#'   \code{\link{waic}}). This argument can be used as an alternative to
#'   specifying the models in \code{...}.
#'
#' @return A vector or matrix with class \code{'compare.loo'} that has its own
#'   print method. If exactly two objects are provided in \code{...} or
#'   \code{x}, then the difference in expected predictive accuracy and the
#'   standard error of the difference are returned (see Details). \emph{The
#'   difference will be positive if the expected predictive accuracy for the
#'   second model is higher.} If more than two objects are provided then a
#'   matrix of summary information is returned.
#'
#' @details When comparing two fitted models, we can estimate the difference in
#'   their expected predictive accuracy by the difference in \code{elpd_loo} or
#'   \code{elpd_waic} (multiplied by \eqn{-2}, if desired, to be on the deviance
#'   scale). To compute the standard error of this difference we can use a
#'   paired estimate to take advantage of the fact that the same set of \eqn{N}
#'   data points was used to fit both models. These calculations should be most
#'   useful when \eqn{N} is large, because then non-normality of the
#'   distribution is not such an issue when estimating the uncertainty in these
#'   sums. These standard errors, for all their flaws, should give a better
#'   sense of uncertainty than what is obtained using the current standard
#'   approach of comparing differences of deviances to a Chi-squared
#'   distribution, a practice derived for Gaussian linear models or
#'   asymptotically, and which only applies to nested models in any case.
#'
#' @template loo-and-psis-references
#'
#' @examples
#' \dontrun{
#' loo1 <- loo(log_lik1)
#' loo2 <- loo(log_lik2)
#' print(compare(loo1, loo2), digits = 3)
#'
#' waic1 <- waic(log_lik1)
#' waic2 <- waic(log_lik2)
#' compare(waic1, waic2)
#' }
#'
compare <- function(..., x) {
  dots <- list(...)
  if (length(dots)) {
    if (!missing(x))
      stop("If 'x' is specified then '...' should not be specified.")
    nms <- as.character(match.call(expand.dots = TRUE))[-1L]
  } else {
    if (!is.list(x) || !length(x))
      stop("'x' must be a list.")
    dots <- x
    nms <- names(dots)
    if (!length(nms))
      nms <- paste0("model", seq_along(dots))
  }
  if (!all(sapply(dots, is.loo)))
    stop("All inputs should have class 'loo'.")
  if (length(dots) <= 1L) {
    stop("'compare' requires at least two models.")
  } else if (length(dots) == 2L) {
    a <- dots[[1L]]
    b <- dots[[2L]]
    pa <- a$pointwise
    pb <- b$pointwise
    Na <- nrow(pa)
    Nb <- nrow(pb)
    if (Na != Nb)
      stop(paste("Models don't have the same number of data points.",
                 "\nFound N_1 =", Na, "and N_2 =", Nb), call. = FALSE)
    sqrtN <- sqrt(Na)
    elpd <- grep("^elpd", colnames(pa))
    diff <- pb[, elpd] - pa[, elpd]
    comp <- c(elpd_diff = sum(diff), se = sqrtN * sd(diff))
    structure(comp, class = "compare.loo")
  } else {
    Ns <- sapply(dots, function(x) nrow(x$pointwise))
    if (!all(Ns == Ns[1L]))
      stop("Not all models have the same number of data points.", call. = FALSE)
    sel <- grep("pointwise|pareto_k", names(dots[[1L]]), invert = TRUE)
    x <- sapply(dots, function(x) unlist(x[sel]))
    colnames(x) <- nms
    rnms <- rownames(x)
    comp <- x
    patts <- c("^waic$|^looic$", "^se_waic$|^se_looic$", "elpd", "p_")
    row_ord <- unlist(sapply(patts, function(p) grep(p, rownames(comp))),
                      use.names = FALSE)
    col_ord <- order(x[grep("^elpd", rnms), ], decreasing = TRUE)
    comp <- t(comp[row_ord, col_ord])
    class(comp) <- c("compare.loo", class(comp))
    comp
  }
}

#' @export
print.compare.loo <- function(x, ..., digits = 1) {
  print(.fr(x, digits), quote = FALSE)
  invisible(x)
}
