#' Model comparison (deprecated, old version)
#'
#' **This function will be deprecated in a future release**. Please
#' use the new [loo_compare()] function instead.
#'
#' @export
#' @param ... At least two objects returned by [loo()] (or [waic()]).
#' @param x A list of at least two objects returned by [loo()] (or
#'   [waic()]). This argument can be used as an alternative to
#'   specifying the objects in `...`.
#'
#' @return A vector or matrix with class `'compare.loo'` that has its own
#'   print method. If exactly two objects are provided in `...` or
#'   `x`, then the difference in expected predictive accuracy and the
#'   standard error of the difference are returned. If more than two objects are
#'   provided then a matrix of summary information is returned (see **Details**).
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in `elpd_loo` or
#'   `elpd_waic` (or multiplied by -2, if desired, to be on the
#'   deviance scale).
#'
#'   *When that difference, `elpd_diff`, is positive then the expected
#'   predictive accuracy for the second model is higher. A negative
#'   `elpd_diff` favors the first model.*
#'
#'   When using `compare()` with more than two models, the values in the
#'   `elpd_diff` and `se_diff` columns of the returned matrix are
#'   computed by making pairwise comparisons between each model and the model
#'   with the best ELPD (i.e., the model in the first row).
#'   Although the `elpd_diff` column is equal to the difference in
#'   `elpd_loo`, do not expect the `se_diff` column to be equal to the
#'   the difference in `se_elpd_loo`.
#'
#'   To compute the standard error of the difference in ELPD we use a
#'   paired estimate to take advantage of the fact that the same set of _N_
#'   data points was used to fit both models. These calculations should be most
#'   useful when _N_ is large, because then non-normality of the
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
#' print(compare(x = list(loo1, loo2)))
#'
#' waic1 <- waic(log_lik1)
#' waic2 <- waic(log_lik2)
#' compare(waic1, waic2)
#' }
#'
compare <- function(..., x = list()) {
  .Deprecated("loo_compare")
  dots <- list(...)
  if (length(dots)) {
    if (length(x)) {
      stop("If 'x' is specified then '...' should not be specified.",
           call. = FALSE)
    }
    nms <- as.character(match.call(expand.dots = TRUE))[-1L]
  } else {
    if (!is.list(x) || !length(x)) {
      stop("'x' must be a list.", call. = FALSE)
    }
    dots <- x
    nms <- names(dots)
    if (!length(nms)) {
      nms <- paste0("model", seq_along(dots))
    }
  }

  if (!all(sapply(dots, is.loo))) {
    stop("All inputs should have class 'loo'.")
  }
  if (length(dots) <= 1L) {
    stop("'compare' requires at least two models.")
  } else if (length(dots) == 2L) {
    loo1 <- dots[[1]]
    loo2 <- dots[[2]]
    comp <- compare_two_models(loo1, loo2)
    class(comp) <- c(class(comp), "old_compare.loo")
    return(comp)
  } else {
    Ns <- sapply(dots, function(x) nrow(x$pointwise))
    if (!all(Ns == Ns[1L])) {
      stop("Not all models have the same number of data points.", call. = FALSE)
    }

    x <- sapply(dots, function(x) {
      est <- x$estimates
      setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))) )
    })
    colnames(x) <- nms
    rnms <- rownames(x)
    comp <- x
    ord <- order(x[grep("^elpd", rnms), ], decreasing = TRUE)
    comp <- t(comp)[ord, ]
    patts <- c("elpd", "p_", "^waic$|^looic$", "^se_waic$|^se_looic$")
    col_ord <- unlist(sapply(patts, function(p) grep(p, colnames(comp))),
                      use.names = FALSE)
    comp <- comp[, col_ord]

    # compute elpd_diff and se_elpd_diff relative to best model
    rnms <- rownames(comp)
    diffs <- mapply(elpd_diffs, dots[ord[1]], dots[ord])
    elpd_diff <- apply(diffs, 2, sum)
    se_diff <- apply(diffs, 2, se_elpd_diff)
    comp <- cbind(elpd_diff = elpd_diff, se_diff = se_diff, comp)
    rownames(comp) <- rnms
    class(comp) <- c("compare.loo", class(comp), "old_compare.loo")
    comp
  }
}



# internal ----------------------------------------------------------------
compare_two_models <- function(loo_a, loo_b, return = c("elpd_diff", "se"), check_dims = TRUE) {
  if (check_dims) {
    if (dim(loo_a$pointwise)[1] != dim(loo_b$pointwise)[1]) {
      stop(paste("Models don't have the same number of data points.",
                 "\nFound N_1 =", dim(loo_a$pointwise)[1], "and N_2 =", dim(loo_b$pointwise)[1]), call. = FALSE)
    }
  }

  diffs <- elpd_diffs(loo_a, loo_b)
  comp <- c(elpd_diff = sum(diffs), se = se_elpd_diff(diffs))
  structure(comp, class = "compare.loo")
}
