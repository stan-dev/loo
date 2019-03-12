#' Model comparison
#'
#' Compare fitted models on LOO or WAIC.
#'
#' @export
#' @param x An object of class \code{"loo"} or a list of such objects.
#' @param ... Additional objects of class \code{"loo"}.
#'
#' @return A matrix with class \code{"compare.loo"} that has its own
#'   print method. See the \strong{Details} section for more .
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in \code{elpd_loo} or
#'   \code{elpd_waic} (or multiplied by \eqn{-2}, if desired, to be on the
#'   deviance scale).
#'
#'   When using \code{loo_compare()}, the returned matrix will have one row per
#'   model and several columns of estimates. The values in the \code{elpd_diff}
#'   and \code{se_diff} columns of the returned matrix are computed by making
#'   pairwise comparisons between each model and the model with the largest ELPD
#'   (the model in the first row). For this reason the \code{elpd_diff} column
#'   will always have the value \code{0} in the first row (i.e., the difference
#'   between the preferred model and itself) and negative values in subsequent
#'   rows for the remaining models.
#'
#'   To compute the standard error of the difference in ELPD --- which should
#'   not be expected to equal the difference of the standard errors --- we use a
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
#' print(loo_compare(loo1, loo2), digits = 3)
#' print(loo_compare(x = list(loo1, loo2)))
#'
#' waic1 <- waic(log_lik1)
#' waic2 <- waic(log_lik2)
#' loo_compare(waic1, waic2)
#' }
#'
loo_compare <- function(x, ...) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ...) {
  if (is.loo(x)) {
    dots <- list(...)
    loos <- c(list(x), dots)
  } else {
    if (!is.list(x) || !length(x)) {
      stop("'x' must be a list if not a 'loo' object.")
    }
    if (length(list(...))) {
      stop("If 'x' is a list then '...' should not be specified.")
    }
    loos <- x
  }

  if (!all(sapply(loos, is.loo))) {
    stop("All inputs should have class 'loo'.")
  }
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.")
  }

  Ns <- sapply(loos, function(x) nrow(x$pointwise))
  if (!all(Ns == Ns[1L])) {
    stop("Not all models have the same number of data points.")
  }

  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))))
  })

  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  comp <- tmp
  ord <- order(tmp[grep("^elpd", rnms), ], decreasing = TRUE)
  comp <- t(comp)[ord, ]
  patts <- c("elpd", "p_", "^waic$|^looic$", "^se_waic$|^se_looic$")
  col_ord <- unlist(sapply(patts, function(p) grep(p, colnames(comp))),
                    use.names = FALSE)
  comp <- comp[, col_ord]

  # compute elpd_diff and se_elpd_diff relative to best model
  rnms <- rownames(comp)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[1]], loos[ord])
  elpd_diff <- apply(diffs, 2, sum)
  se_diff <- apply(diffs, 2, se_elpd_diff)
  comp <- cbind(elpd_diff = elpd_diff, se_diff = se_diff, comp)
  rownames(comp) <- rnms

  class(comp) <- c("compare.loo", class(comp))
  return(comp)
}

#' @rdname loo_compare
#' @export
#' @param digits For the print method only, the number of digits to use when
#'   printing.
#' @param simplify For the print method only, should only the essential columns
#'   of the summary matrix be printed? The entire matrix is always returned, but
#'   by default only the most important columns are printed.
print.compare.loo <- function(x, ..., digits = 1, simplify = TRUE) {
  xcopy <- x
  if (inherits(xcopy, "old_compare.loo")) {
    if (NCOL(xcopy) >= 2 && simplify) {
      patts <- "^elpd_|^se_diff|^p_|^waic$|^looic$"
      xcopy <- xcopy[, grepl(patts, colnames(xcopy))]
    }
  } else if (NCOL(xcopy) >= 2 && simplify) {
     xcopy <- xcopy[, c("elpd_diff", "se_diff")]
  }
  print(.fr(xcopy, digits), quote = FALSE)
  invisible(x)
}



# internal ----------------------------------------------------------------

#' Compute pointwise elpd differences
#' @noRd
#' @param loo_a,loo_b Two loo objects.
elpd_diffs <- function(loo_a, loo_b) {
  pt_a <- loo_a$pointwise
  pt_b <- loo_b$pointwise
  elpd <- grep("^elpd", colnames(pt_a))
  pt_b[, elpd] - pt_a[, elpd]
}

#' Compute standard error of the elpd difference
#' @noRd
#' @param diffs Vector of pointwise elpd differences
se_elpd_diff <- function(diffs) {
  N <- length(diffs)
  sqrt(N) * sd(diffs)
}



#' Find the model names associated with loo objects
#'
#' @export
#' @keywords internal
#' @param x List of loo objects.
#' @return Character vector of model names the same length as x.
#'
find_model_names <- function(x) {
  stopifnot(is.list(x))
  out_names <- character(length(x))

  names1 <- names(x)
  names2 <- lapply(x, "attr", "model_name", exact = TRUE)
  names3 <- lapply(x, "[[", "model_name")
  names4 <- paste0("model", seq_along(x))

  for (j in seq_along(x)) {
    if (isTRUE(nzchar(names1[j]))) {
      out_names[j] <- names1[j]
    } else if (length(names2[[j]])) {
      out_names[j] <- names2[[j]]
    } else if (length(names3[[j]])) {
      out_names[j] <- names3[[j]]
    } else {
      out_names[j] <- names4[j]
    }
  }

  return(out_names)
}
