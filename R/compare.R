#' Model comparison
#'
#' Compare two fitted models on LOO or WAIC
#'
#' @export
#' @param a,b objects returned by \code{\link{loo}} or \code{\link{waic}}.
#' @return A named list with class \code{'compare.loo'}.
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

compare <- function(a, b) {
  namea <- deparse(substitute(a))
  nameb <- deparse(substitute(b))
  if (!is.loo(a))
    stop(paste(namea, "does not have class 'loo'"), call. = FALSE)
  if (!is.loo(b))
    stop(paste(nameb, "does not have class 'loo'"), call. = FALSE)

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
  comp <- list(elpd_diff = sum(diff), se = sqrtN * sd(diff))
  class(comp) <- "compare.loo"
  comp
}
