#' Model comparison (previous name)
#'
#' @description
#' `loo_compare()` is the previous name of [model_compare()] and is kept as a
#' working alias. It is a generic so that packages registering their own
#' `loo_compare` methods keep dispatching as before; the default method simply
#' forwards to [model_compare()].
#'
#' New code should call [model_compare()], which additionally compares
#' [`kfold_pred_measure()`][kfold_pred_measure],
#' [`test_pred_measure()`][test_pred_measure], and
#' [`insample_pred_measure()`][insample_pred_measure] results.
#'
#' @export
#' @inheritParams model_compare
#' @return See [model_compare()].
#'
#' @seealso [model_compare()]
#'
#' @examples
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)
#' loo2 <- loo(LL + 1)
#'
#' comp <- loo_compare(loo1, loo2, loo3)
#' print(comp, digits = 2)
#' print(comp, simplify = FALSE) # full table
#'
loo_compare <- function(x, ..., rank_by = NULL, custom_se_fn) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ..., rank_by = NULL, custom_se_fn) {
  # Forward `custom_se_fn` only when it was actually supplied, so that
  # `model_compare()` can still tell an omitted argument from an explicit NULL.
  if (missing(custom_se_fn)) {
    model_compare(x, ..., rank_by = rank_by)
  } else {
    model_compare(x, ..., rank_by = rank_by, custom_se_fn = custom_se_fn)
  }
}

#' @rdname loo_compare
#' @export
#' @param digits For the print method only, the number of digits to use when
#'   printing.
#' @param p_worse For the print method only, should we include the normal
#'   approximation based probability of each model having worse performance than
#'   the best model? The default is `TRUE`.
#' @param simplify For the print method only, should the output be simplified to
#'   only include the model names, ELPD differences, and (when `p_worse = TRUE`)
#'   diagnostic columns? The default is `TRUE`. Set to `FALSE` to also print the
#'   available estimate columns (pointwise ELPD, LOOIC/WAIC, and their standard
#'   errors).
print.compare.loo <- function(x, ..., digits = 1, p_worse = TRUE,
                              simplify = TRUE) {
  if (inherits(x, "old_compare.loo")) {
    return(unclass(x))
  }
  if (!inherits(x, "data.frame")) {
    class(x) <- c(class(x), "data.frame")
  }
  if (!all(c("model", "elpd_diff", "se_diff") %in% colnames(x))) {
    print(as.data.frame(x))
    return(x)
  }
  base_cols <- c("model", "elpd_diff", "se_diff")
  diag_cols <- c("p_worse", "diag_diff", "diag_elpd")
  show_diag <- p_worse && "p_worse" %in% colnames(x)

  estimate_cols <- setdiff(colnames(x), c(base_cols, diag_cols))
  estimate_cols <- estimate_cols[vapply(x[estimate_cols], is.numeric, logical(1))]

  cols <- c(
    base_cols,
    if (show_diag) diag_cols,
    if (!simplify) estimate_cols
  )
  cols <- intersect(cols, colnames(x))

  x2 <- x[, cols, drop = FALSE]

  fmt_cols <- setdiff(cols, c("model", "diag_diff", "diag_elpd"))
  if (length(fmt_cols)) {
    if ("p_worse" %in% fmt_cols) {
      x2$p_worse <- .fr(x2$p_worse, digits = 2)
      fmt_cols <- setdiff(fmt_cols, "p_worse")
    }
    if (length(fmt_cols)) {
      x2[fmt_cols] <- .fr(x2[fmt_cols], digits)
    }
  }
  # Use `as.data.frame(x2)` here to drop "compare.loo"
  # so print() uses print.data.frame.
  print(as.data.frame(x2), quote = FALSE, row.names = FALSE)

  # show glossary for diagnostic flags
  has_diag <- any(nzchar(x[["diag_diff"]], keepNA = FALSE), na.rm = TRUE) ||
              any(nzchar(x[["diag_elpd"]], keepNA = FALSE), na.rm = TRUE)
  if (has_diag && p_worse) {
    message(
      "\nDiagnostic flags present.\n",
      "See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)\n",
      "or https://mc-stan.org/loo/reference/loo-glossary.html."
    )
  }
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

  Ns <- vapply(loos, function(x) nrow(x$pointwise), integer(1))
  if (any(Ns != Ns[1L])) {
    stop(
      paste0(
        "All models must have the same number of observations, but models have inconsistent observation counts: ",
        paste(paste0("'", find_model_names(loos), "' (", Ns, ")"), collapse = ", ")
      ),
      call. = FALSE
    )
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

#' Perform checks on `"loo"` objects __after__ comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param ord List of `"loo"` object orderings.
#' @return Nothing, just possibly throws errors/warnings.
loo_order_stat_check <- function(loos, ord) {

  ## breaks

  if (length(loos) <= 11L) {
    # procedure cannot be diagnosed for fewer than ten candidate models
    # (total models = worst model + ten candidates)
    # break from function
    return(NULL)
  }

  ## warnings

  # compute the elpd differences from the median model
  baseline_idx <- middle_idx(ord)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[baseline_idx]], loos[ord])
  elpd_diff <- apply(diffs, 2, sum)

  # estimate the standard deviation of the upper-half-normal
  diff_median <- stats::median(elpd_diff)
  elpd_diff_trunc <- elpd_diff[elpd_diff >= diff_median]
  n_models <- sum(!is.na(elpd_diff_trunc))
  candidate_sd <- sqrt(1 / n_models * sum(elpd_diff_trunc^2, na.rm = TRUE))

  # estimate expected best diff under null hypothesis
  K <- length(loos) - 1
  order_stat <- order_stat_heuristic(K, candidate_sd)

  if (max(elpd_diff) <= order_stat) {
    # flag warning if we suspect no model is theoretically better than the baseline
    warning("Difference in performance potentially due to chance. ",
            "See McLatchie and Vehtari (2023) for details.",
            call. = FALSE)
  }
}

#' Returns the middle index of a vector
#' @noRd
#' @param vec A vector.
#' @return Integer index value.
middle_idx <- function(vec) floor(length(vec) / 2)

#' Computes maximum order statistic from K Gaussians
#' @noRd
#' @param K Number of Gaussians.
#' @param c Scaling of the order statistic.
#' @return Numeric expected maximum from K samples from a Gaussian with mean
#' zero and scale `"c"`
order_stat_heuristic <- function(K, c) {
  qnorm(p = 1 - 1 / (K * 2), mean = 0, sd = c)
}

#' Count number of high Pareto k values in PSIS-LOO and create diagnostic message
#' @noRd
#' @param loos Ordered list of loo objects.
#' @return Character vector of diagnostic messages.
diag_elpd <- function(loos) {
  sapply(loos, function(loo) {
    k <- loo$diagnostics[["pareto_k"]]
    if (is.null(k)) {
      out <- ""
    } else {
      S <- dim(loo)[1]
      khat_threshold <- ps_khat_threshold(S)
      K <- sum(k > khat_threshold)
      out <- ifelse(K == 0, "", paste0(K, " k_psis > ", round(khat_threshold, 2)))
    }
    out
  })
}

#' Create diagnostic for elpd differences
#' @noRd
#' @param N Number of data points.
#' @param elpd_diff Vector of elpd differences.
#' @return Character vector of diagnostic messages.
diag_diff <- function(N, elpd_diff) {
  if (N < 100) {
    diag_diff <- rep("N < 100", length(elpd_diff))
    diag_diff[elpd_diff == 0] <- ""
  } else {
    diag_diff <- rep("", length(elpd_diff))
    diag_diff[elpd_diff > -4 & elpd_diff != 0] <- "|elpd_diff| < 4"
  }
  diag_diff
}
