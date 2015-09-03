#' Model comparison
#'
#' Compare fitted models on LOO or WAIC
#'
#' @export
#' @param ... At least two objects returned by \code{\link{loo}} or
#'   \code{\link{waic}}.
#' @return A vector or matrix with class \code{'compare.loo'}. If \code{...} has
#'   more than two objects then a matrix is returned. This matrix summarizes the
#'   objects and also reports model weights (the posterior probability that each
#'   model has the best expected out-of-sample predictive accuracy). If
#'   \code{...} contains exactly two objects then the difference in expected
#'   predictive accuracy and the standard error of the difference are returned
#'   (see Details) in addition to model weights.
#'
#' @details When comparing two fitted models, we can estimate the difference in
#'   their expected predictive accuracy by the difference in \code{elpd_waic} or
#'   \code{elpd_loo} (multiplied by -2, if desired, to be on the deviance
#'   scale). To compute the standard error of this difference we can use a
#'   paired estimate to take advantage of the fact that the same set of \eqn{N}
#'   data points is being used to fit both models. We would think that these
#'   calculations would be most useful when \eqn{N} is large, because then
#'   non-normality of the distribution is not such an issue when estimating the
#'   uncertainty of these sums. In any case, we suspect that these standard
#'   errors, for all their flaws, should give a better sense of uncertainty than
#'   what is obtained using the current standard approach of comparing
#'   differences of deviances to a Chi-squared distribution, a practice derived
#'   for Gaussian linear models or asymptotically and which only applies to
#'   nested models in any case.
#'
#' @seealso \code{\link{loo}}, \code{\link{waic}},
#'   \code{\link{print.compare.loo}}
#' @examples
#' \dontrun{
#' loo1 <- loo(log_lik1)
#' loo2 <- loo(log_lik2)
#' diff <- compare(loo1, loo2)
#' print(diff, digits = 1)
#'
#' waic1 <- waic(log_lik1)
#' waic2 <- waic(log_lik2)
#' compare(waic1, waic2)
#' }
#'

compare <- function(...) {
  dots <- list(...)
  nms <- as.character(match.call())[-1L]
  if (!all(sapply(dots, is.loo)))
    stop("All inputs should have class 'loo'", call. = FALSE)

  if (length(dots) <= 1L)
    stop("'compare' requires at least two models.", call. = FALSE)
  else if (length(dots) == 2L) {
    a <- dots[[1L]]
    b <- dots[[2L]]
    pa <- a$pointwise
    pb <- b$pointwise
    Na <- nrow(pa)
    Nb <- nrow(pb)
    if (Na != Nb)
      stop(paste("Models a and b should have the same number of data points.",
                 "\nFound N_a =", Na, "and N_b =", Nb), call. = FALSE)
    sqrtN <- sqrt(Na)
    elpd <- grep("^elpd", colnames(pa))
    diff <- pb[, elpd] - pa[, elpd]
    uwts <- c(sum(pa[, elpd]), sum(pb[, elpd]))
    uwts <- exp(uwts - max(uwts))
    wts <- uwts / sum(uwts)
    comp <- c(elpd_diff = sum(diff), se = sqrtN * sd(diff),
              weight1 = wts[1L], weight2 = wts[2L])
    class(comp) <- "compare.loo"
    comp
  }
  else {
    Ns <- sapply(dots, function(x) nrow(x$pointwise))
    if (!all(Ns == Ns[1L]))
      stop("Not all models have the same number of data points.", call. = FALSE)
    sel <- grep("pointwise|pareto_k", names(dots[[1L]]), invert = TRUE)
    x <- sapply(dots, function(x) unlist(x[sel]))
    colnames(x) <- nms
    rnms <- rownames(x)
    uwts <- x[grep("^elpd", rnms), ]
    uwts <- exp(uwts - max(uwts))
    comp <- rbind(x, weights = uwts / sum(uwts))
    col_ord <- order(uwts, decreasing = TRUE)

    patts <- c("^waic$|^looic$", "^se_waic$|^se_looic$", "elpd", "p_", "weights")
    row_ord <- unlist(sapply(patts, function(p) grep(p, rownames(comp))),
                      use.names = FALSE)
    comp <- t(comp[row_ord, col_ord])
    class(comp) <- c("compare.loo", class(comp))
    comp
  }
}
