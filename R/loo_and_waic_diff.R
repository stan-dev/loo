#' Model comparison
#'
#' Compare two fitted models
#'
#' @export
#' @param loo1,loo2 lists (of class \code{'loo'}) returned by
#' \code{\link{loo_and_waic}}.
#' @return a named list with class \code{'compare.loo'}.
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
#' @seealso \code{\link{loo_and_waic}}, \code{\link{loo-package}}
#' @examples
#' \dontrun{
#' loo1 <- loo_and_waic(log_lik1)
#' loo2 <- loo_and_waic(log_lik2)
#' diff <- loo_and_waic_diff(loo1, loo2)
#' print(diff, digits = 1)
#' }
#'
loo_and_waic_diff <- function(loo1, loo2) {
  N1 <- nrow(loo1$pointwise)
  N2 <- nrow(loo2$pointwise)
  if (N1 != N2) {
    stop(paste("Models being compared should have the same number of data points.",
               "Found N1 =", N1, "and N2 =", N2))
  }
  sqrtN <- sqrt(N1)
  loo_diff <- loo2$pointwise[, "elpd_loo"] - loo1$pointwise[, "elpd_loo"]
  waic_diff <- loo2$pointwise[, "elpd_waic"] - loo1$pointwise[, "elpd_waic"]

  out <- list(elpd_loo_diff = sum(loo_diff),
              lpd_loo_diff = sqrtN * sd(loo_diff),
              elpd_waic_diff = sum(waic_diff),
              lpd_waic_diff = sqrtN * sd(waic_diff))
  class(out) <- "compare.loo"
  out
}
