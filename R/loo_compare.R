#' Model comparison
#'
#' @description Compare fitted models based on [ELPD][loo-glossary].
#'
#'   By default the print method shows only the most important information. Use
#'   `print(..., simplify=FALSE)` to print a more detailed summary.
#'
#' @export
#' @param x An object of class `"loo"` or a list of such objects. If a list is
#'   used then the list names will be used as the model names in the output. See
#'   **Examples**.
#' @param ... Additional objects of class `"loo"`, if not passed in as a single
#'   list.
#'
#' @return A matrix with class `"compare.loo"` that has its own
#'   print method. See the **Details** section.
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in
#'   [`elpd_loo`][loo-glossary] or `elpd_waic` (or multiplied by \eqn{-2}, if
#'   desired, to be on the deviance scale).
#'
#'   When using `loo_compare()`, the returned matrix will have one row per model
#'   and several columns of estimates. The values in the
#'   [`elpd_diff`][loo-glossary] and [`se_diff`][loo-glossary] columns of the
#'   returned matrix are computed by making pairwise comparisons between each
#'   model and the model with the largest ELPD (the model in the first row). For
#'   this reason the `elpd_diff` column will always have the value `0` in the
#'   first row (i.e., the difference between the preferred model and itself) and
#'   negative values in subsequent rows for the remaining models.
#'
#'   To compute the standard error of the difference in [ELPD][loo-glossary] ---
#'   which should not be expected to equal the difference of the standard errors
#'   --- we use a paired estimate to take advantage of the fact that the same
#'   set of \eqn{N} data points was used to fit both models. These calculations
#'   should be most useful when \eqn{N} is large, because then non-normality of
#'   the distribution is not such an issue when estimating the uncertainty in
#'   these sums. These standard errors, for all their flaws, should give a
#'   better sense of uncertainty than what is obtained using the current
#'   standard approach of comparing differences of deviances to a Chi-squared
#'   distribution, a practice derived for Gaussian linear models or
#'   asymptotically, and which only applies to nested models in any case.
#'
#' @seealso
#' * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
#' @template loo-and-psis-references
#'
#' @examples
#' # very artificial example, just for demonstration!
#' LL <- example_loglik_array()
#' loo1 <- loo(LL, r_eff = NA)     # should be worst model when compared
#' loo2 <- loo(LL + 1, r_eff = NA) # should be second best model when compared
#' loo3 <- loo(LL + 2, r_eff = NA) # should be best model when compared
#'
#' comp <- loo_compare(loo1, loo2, loo3)
#' print(comp, digits = 2)
#'
#' # show more details with simplify=FALSE
#' # (will be the same for all models in this artificial example)
#' print(comp, simplify = FALSE, digits = 3)
#'
#' # can use a list of objects with custom names
#' # will use apple, banana, and cherry, as the names in the output
#' loo_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#'
#' \dontrun{
#' # works for waic (and kfold) too
#' loo_compare(waic(LL), waic(LL - 10))
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

  # If subsampling is used
  if (any(sapply(loos, inherits, "psis_loo_ss"))) {
    return(loo_compare.psis_loo_ss_list(loos))
  }

  loo_compare_checks(loos)

  comp <- loo_compare_matrix(loos)
  ord <- loo_compare_order(loos)

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
#' @param loo_a,loo_b Two `"loo"` objects.
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
  # As `elpd_diff` is defined as the sum of N independent components,
  # we can compute the standard error by using the standard deviation
  # of the N components and multiplying by `sqrt(N)`.
  sqrt(N) * sd(diffs)
}


#' Perform checks on `"loo"` objects before comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @return Nothing, just possibly throws errors/warnings.
loo_compare_checks <- function(loos) {
  ## errors
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.", call.=FALSE)
  }
  if (!all(sapply(loos, is.loo))) {
    stop("All inputs should have class 'loo'.", call.=FALSE)
  }

  Ns <- sapply(loos, function(x) nrow(x$pointwise))
  if (!all(Ns == Ns[1L])) {
    stop("Not all models have the same number of data points.", call.=FALSE)
  }

  ## warnings

  yhash <- lapply(loos, attr, which = "yhash")
  yhash_ok <- sapply(yhash, function(x) { # ok only if all yhash are same (all NULL is ok)
    isTRUE(all.equal(x, yhash[[1]]))
  })
  if (!all(yhash_ok)) {
    warning("Not all models have the same y variable. ('yhash' attributes do not match)",
            call. = FALSE)
  }

  if (all(sapply(loos, is.kfold))) {
    Ks <- unlist(lapply(loos, attr, which = "K"))
    if (!all(Ks == Ks[1])) {
      warning("Not all kfold objects have the same K value. ",
              "For a more accurate comparison use the same number of folds. ",
              call. = FALSE)
    }
  } else if (any(sapply(loos, is.kfold)) && any(sapply(loos, is.psis_loo))) {
    warning("Comparing LOO-CV to K-fold-CV. ",
            "For a more accurate comparison use the same number of folds ",
            "or loo for all models compared.",
            call. = FALSE)
  }
}


#' Find the model names associated with `"loo"` objects
#'
#' @export
#' @keywords internal
#' @param x List of `"loo"` objects.
#' @return Character vector of model names the same length as `x.`
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
  out_names
}


#' Compute the loo_compare matrix
#' @keywords internal
#' @noRd
#' @param loos List of `"loo"` objects.
loo_compare_matrix <- function(loos){
  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))))
  })
  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  comp <- tmp
  ord <- loo_compare_order(loos)
  comp <- t(comp)[ord, ]
  patts <- c("elpd", "p_", "^waic$|^looic$", "^se_waic$|^se_looic$")
  col_ord <- unlist(sapply(patts, function(p) grep(p, colnames(comp))),
                    use.names = FALSE)
  comp <- comp[, col_ord]
  comp
}

#' Computes the order of loos for comparison
#' @noRd
#' @keywords internal
#' @param loos List of `"loo"` objects.
loo_compare_order <- function(loos){
  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est), paste0("se_", rownames(est))))
  })
  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  ord <- order(tmp[grep("^elpd", rnms), ], decreasing = TRUE)
  ord
}
