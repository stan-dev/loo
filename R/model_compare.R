#' Model comparison
#'
#' @description Compare fitted models on [ELPD][loo-glossary] or, for
#'   [`pred_measure`][pred_measure] results, on several predictive performance
#'   measures at once.
#'
#'   `model_compare()` accepts two families of input:
#'
#'   * **Classic results** --- `"loo"`, `"waic"`, and `"kfold"` objects, compared
#'     on ELPD alone.
#'   * **Predictive measure results** --- objects from
#'     [`loo_pred_measure()`][loo_pred_measure],
#'     [`kfold_pred_measure()`][kfold_pred_measure],
#'     [`test_pred_measure()`][test_pred_measure], or
#'     [`insample_pred_measure()`][insample_pred_measure], compared on every
#'     measure the models share.
#'
#'   All models in one call must be evaluated the same way. Differences between,
#'   say, a LOO and a k-fold result would contrast held-out schemes rather than
#'   models, so mixed inputs are an error.
#'
#' @export
#' @param x An object of class `"loo"` or `"pred_measure"`, or a list of such
#'   objects. List names are used as the model names in the output. See
#'   **Examples**.
#' @param ... Additional objects of class `"loo"` or `"pred_measure"`, if not
#'   passed in as a single list. Naming every model here, as in
#'   `model_compare(A = m1, B = m2)`, names the models in the output, exactly as
#'   the list form does.
#' @return A data frame of class `"compare.loo"` with one row per model and its
#'   own print method.
#'
#'   For classic `"loo"` / `"waic"` / `"kfold"` comparisons the columns are
#'   unchanged from previous versions: `model`, `elpd_diff`, `se_diff`,
#'   `p_worse`, `diag_diff`, `diag_elpd`, and the estimate columns of the input
#'   objects.
#'
#'   For [`pred_measure`][pred_measure] comparisons there is a `{measure}_diff`
#'   and a `{measure}_se_diff` column for every measure shared by all models
#'   (e.g. `rmse_diff`, `rmse_se_diff`). ELPD-family measures use `elpd_diff`
#'   and `se_diff` instead. `p_worse` and `diag_diff` are computed for ELPD
#'   only. `diag_elpd` holds per-model Pareto \eqn{\hat{k}} diagnostics and is
#'   present only for [`loo_pred_measure()`][loo_pred_measure] comparisons, the
#'   only source with Pareto \eqn{\hat{k}} values.
#'
#'   The object also carries the following attributes:
#'   \describe{
#'     \item{`compare_reference`}{
#'       A named character vector giving, for each measure, the model its
#'       differences were computed against, which is that measure's own best
#'       model.
#'     }
#'     \item{`compare_measures`}{
#'       Bare names of all measures that were compared.
#'     }
#'     \item{`sign_converted_measures`}{
#'       Bare names of the loss measures whose sign was flipped onto the utility
#'       scale.
#'     }
#'     \item{`compare_source`}{
#'       The shared evaluation source: `"loo"`, `"kfold"`, `"test"`, or
#'       `"insample"`.
#'     }
#'   }
#'   `compare_reference` is set for every comparison; the last three are set
#'   for [`pred_measure`][pred_measure] comparisons only.
#'
#' @details
#' ## Differences and their standard errors
#'   Differences are pairwise: every model is compared with one reference model,
#'   whose own `{measure}_diff` is therefore `0`. The reference is the best
#'   model on that measure, so `mse_diff` may use a different reference than
#'   `elpd_diff`, and the remaining differences for a measure are all negative.
#'   Rows are ordered by `"elpd"` when all models share it. Otherwise, rows are
#'   ordered by the first shared measure in alphabetical order.
#'
#'   The standard error of a difference is a paired estimate, which uses the
#'   fact that the same \eqn{N} data points were used for both models. It should
#'   not be expected to equal the difference of the two models' standard errors.
#'
#' ## `p_worse`, `diag_diff`, and `diag_elpd`
#'   `p_worse` is the probability that a model has worse ELPD than the reference
#'   model, computed with a normal approximation from `elpd_diff` and `se_diff`.
#'   Sivula et al. (2025) give the conditions under which that approximation is
#'   good; `diag_diff` reports the two that fail most often:
#'
#'   * `N < 100` (small data)
#'   * `|elpd_diff| < 4` (models make similar predictions)
#'
#'   Either message means the error distribution is skewed or thick tailed, the
#'   normal approximation is not well calibrated, and `p_worse` is likely too
#'   large. If `|elpd_diff|` is many times `se_diff` the difference is
#'   quite certain. Model misspecification and outliers also skew the error
#'   distribution, and can be diagnosed with the usual predictive checks.
#'
#'   `diag_elpd` reports the PSIS-LOO Pareto \eqn{\hat{k}} diagnostic for each
#'   model's pointwise ELPD. An entry `K k_psis > 0.7`, where `K` counts the
#'   high Pareto \eqn{\hat{k}} values, warns of possible bias in `elpd_diff`
#'   favoring models with many such values. Pareto \eqn{\hat{k}} describes a
#'   model's PSIS-LOO approximation rather than any one measure or pair of
#'   models, and every LOO measure uses the same importance weights, so for
#'   `pred_measure` comparisons `print()` reports it once per model in a block
#'   above the difference tables instead of as a column inside one of them. The
#'   `diag_elpd` column is still returned on the object.
#'
#' ## Comparing `pred_measure` objects
#'   When all inputs are predictive measure results sharing one evaluation
#'   source, paired differences are computed for every measure present in all
#'   models. Measures are matched on their bare names, so the source suffix
#'   (`_loo`, `_kfold`, `_test`, or none for in-sample) is handled
#'   transparently. When the models were evaluated on different `measure` sets,
#'   only the shared measures are compared and a warning lists the omitted ones.
#'
#'   The data frame carries one row order for all measures, but each *printed*
#'   measure table is sorted by its own difference, so the best model on that
#'   measure always leads its table and the differences run in decreasing order.
#'   Use `print(x, measures = "all")` to display a table for every compared
#'   measure; see [loo-glossary] for column definitions.
#'
#' ## Utility scale and sign conversion
#'   Measures differ in orientation in their raw form: ELPD and SRPS/SCRPS are
#'   utilities (higher is better), while MSE, RPS/CRPS and the Brier score are
#'   losses (lower is better). All `{measure}_diff` values are reported on a
#'   common utility scale, so loss measures have their sign flipped and a
#'   negative `{measure}_diff` always means worse performance than the
#'   reference. Which measures are losses is recorded in the `loss` element of
#'   each measure's entry in the `measure_info` attribute of an
#'   `*_pred_measure()` result. The flipped measures are named in the
#'   `sign_converted_measures` attribute. `print()` marks them with
#'   "sign flipped" in the table header and names them below the tables.
#'
#'   A custom measure is treated as a utility unless it declares otherwise with
#'   `attr(my_fun, "measure_loss") <- TRUE`. The declaration also determines the
#'   ranking direction, so an undeclared loss is both flipped and ranked in the
#'   wrong direction; see [insample_pred_measure()].
#'
#' ## Standard error of a measure difference
#'   How `{measure}_se_diff` is obtained is recorded in the `diff_method`
#'   element of the measure's entry in `measure_info`:
#'
#'   * `"sum"` or `"mean"`: the overall estimate is the sum (`elpd`, `ic`) or the
#'     mean (`mlpd`, `mae`, `mse`, `acc`, `rps`, `srps`, `brier`) of its
#'     pointwise contributions, so the standard error is computed from paired
#'     pointwise differences (the same formula as `se_diff`).
#'   * `"measure_specific"`: the overall estimate is not a sum or mean of
#'     pointwise contributions (`r2`, `rmse`, `bacc`), so the measure supplies
#'     its own standard error of the difference.
#'   * `"custom"`: the standard error comes from the measure's own
#'     `attr(my_fun, "measure_se_diff")` declaration, set with
#'     [custom_measure()]. `{measure}_se_diff` is `NA` when the measure
#'     declares nothing.
#'
#' ## Source-specific behavior
#'   Comparisons behave the same way across sources, with three exceptions:
#'
#'   * **`diag_elpd`** is produced only for
#'     [`loo_pred_measure()`][loo_pred_measure] comparisons, since Pareto
#'     \eqn{\hat{k}} diagnostics exist only for PSIS-LOO.
#'   * **K-fold** comparisons warn when the models do not share the same number
#'     of folds, matching the behavior for plain `"kfold"` objects.
#'   * **In-sample** comparisons warn that in-sample scores are optimistically
#'     biased and favor more complex models. They are supported for
#'     completeness, but out-of-sample sources should be preferred for model
#'     selection.
#'
#' ## Warnings for many model comparisons
#'   If more than \eqn{11} models are compared, we internally recompute the model
#'   differences using the median model (by ELPD, or by the ranking measure
#'   for `pred_measure` comparisons) as the baseline, and estimate whether the
#'   differences in predictive performance are potentially due to chance as
#'   described by McLatchie and Vehtari (2023). This flags a warning if there is
#'   a risk of over-fitting due to the selection process. In that case users are
#'   recommended to avoid model selection based on LOO-CV, and instead to favor
#'   model averaging/stacking or projection predictive inference.
#'
#' @seealso
#' * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
#' * The article
#'   [Model comparison: Explanation of `model_compare()`](https://mc-stan.org/loo/articles/articles-online-only/model-comparison.html)
#'   on the __loo__ website, for how the differences and their standard errors
#'   are computed for each measure and when the normal approximation behind
#'   `p_worse` can be trusted.
#' @template loo-and-compare-references
#'
#' @examples
#' # very artificial example, just for demonstration!
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)     # should be worst model when compared
#' loo2 <- loo(LL + 1) # should be second best model when compared
#' loo3 <- loo(LL + 2) # should be best model when compared
#'
#' comp <- model_compare(loo1, loo2, loo3)
#' print(comp, digits = 2)
#'
#' # can use a list of objects with custom names
#' # the names will be used in the output
#' model_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#'
#' \dontrun{
#' # works for waic (and kfold) too
#' model_compare(waic(LL), waic(LL - 10))
#'
#' # compare multiple predictive measures from loo_pred_measure()
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   fit1 <- brms::brm(
#'     Reaction ~ Days, data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   fit2 <- brms::brm(
#'     Reaction ~ poly(Days, 2), data = lme4::sleepstudy,
#'     refresh = 0, chains = 2, iter = 1000
#'   )
#'   pm1 <- loo_pred_measure(
#'     loo = loo(fit1, save_psis = TRUE),
#'     y = fit1$data$Reaction,
#'     mupred = brms::posterior_epred(fit1),
#'     measures = c("rmse", "r2")
#'   )
#'   pm2 <- loo_pred_measure(
#'     loo = loo(fit2, save_psis = TRUE),
#'     y = fit2$data$Reaction,
#'     mupred = brms::posterior_epred(fit2),
#'     measures = c("rmse", "r2")
#'   )
#'   comp <- model_compare(pm1, pm2)
#'   print(comp)
#'   print(comp, measures = "all")
#'
#'   # the same works for k-fold CV
#'   folds <- kfold_split_random(K = 5, N = nrow(lme4::sleepstudy))
#'   kf1 <- brms::kfold(fit1, folds = folds, save_fits = TRUE)
#'   kf2 <- brms::kfold(fit2, folds = folds, save_fits = TRUE)
#'   kpm1 <- kfold_pred_measure(
#'     y = fit1$data$Reaction,
#'     mupred = brms::kfold_predict(kf1, method = "fitted")$yrep,
#'     kfold = kf1,
#'     measures = "rmse"
#'   )
#'   kpm2 <- kfold_pred_measure(
#'     y = fit2$data$Reaction,
#'     mupred = brms::kfold_predict(kf2, method = "fitted")$yrep,
#'     kfold = kf2,
#'     measures = "rmse"
#'   )
#'   model_compare(kpm1, kpm2)
#'
#'   # mixing evaluation sources is an error
#'   try(model_compare(pm1, kpm2))
#' }
#' }
#'
model_compare <- function(x, ...) {
  if (missing(x)) {
    dots <- list(...)
    if (!length(dots)) {
      stop("No models supplied.", call. = FALSE)
    }
    return(model_compare(dots))
  }
  UseMethod("model_compare")
}

#' @rdname model_compare
#' @export
model_compare.default <- function(x, ...) {
  loos <- .model_compare_inputs(x, ...)

  # if subsampling is used
  if (any(sapply(loos, inherits, "psis_loo_ss"))) {
    return(model_compare.psis_loo_ss_list(loos))
  }

  # `pred_measure` objects must be tested before any `is.loo()` check: results
  # from `loo_pred_measure()` and `kfold_pred_measure()` inherit the classes of
  # the `loo`/`kfold` object they were built from.
  is_pm <- vapply(loos, is.pred_measure, logical(1))

  if (all(is_pm)) {
    return(compare_pred_measure(loos))
  }

  if (any(is_pm)) {
    stop(
      "Cannot mix 'pred_measure' objects with plain 'loo' objects. ",
      "Compare models using the same *_pred_measure() function for each model.",
      call. = FALSE
    )
  }

  # run pre-comparison checks
  model_compare_checks(loos)

  # compute elpd_diff and se_elpd_diff relative to best model
  ord <- model_compare_order(loos)
  comp <- model_compare_matrix(loos, ord = ord)
  rnms <- rownames(comp)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[1L]], loos[ord])
  colnames(diffs) <- rnms
  elpd_diff <- apply(diffs, 2, sum)
  se_diff <- apply(diffs, 2, se_elpd_diff)

  # compute probabilities that a model has worse elpd than the best model
  # using a normal approximation
  # (Sivula et al., 2025)
  p_worse <- stats::pnorm(0, elpd_diff, se_diff)
  p_worse[elpd_diff == 0] <- NA

  comp <- cbind(
    data.frame(
      model = rnms,
      elpd_diff = elpd_diff,
      se_diff = se_diff,
      p_worse = p_worse,
      diag_diff = diag_diff(nrow(diffs), elpd_diff),
      diag_elpd = diag_elpd(loos[ord])
    ),
    as.data.frame(comp)
  )
  rownames(comp) <- NULL

  # run order statistics-based checks for many model comparisons
  model_order_stat_check(loos, ord)

  # Same attribute contract as the `pred_measure` path, with the single
  # measure `"elpd"`.
  attr(comp, "compare_reference") <- c(elpd = rnms[[1L]])
  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' Reference model a measure's differences were computed against
#'
#' Each measure has its own best model as reference, recorded in attribute `compare_reference`. Falls back to the first row for objects
#' created before that attribute existed.
#' @noRd
.measure_ref_model <- function(x, measure) {
  refs <- attr(x, "compare_reference")
  if (!is.null(refs) && measure %in% names(refs)) {
    return(refs[[measure]])
  }
  x$model[[1L]]
}

#' Normalize `model_compare()` inputs to a list of model results
#' @noRd
.model_compare_inputs <- function(x, ...) {
  if (is.loo(x) || inherits(x, "pred_measure")) {
    dots <- list(...)
    return(c(list(x), dots))
  }
  if (!is.list(x) || !length(x)) {
    stop(
      "'x' must be a list if not a 'loo' or 'pred_measure' object.",
      call. = FALSE
    )
  }
  if (length(list(...))) {
    stop("If 'x' is a list then '...' should not be specified.", call. = FALSE)
  }
  x
}

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

#' Warn when k-fold results do not share the same number of folds
#' @noRd
#' @param loos List of `"kfold"` or `"kfold_pred_measure"` objects.
throw_kfold_K_mismatch_warning <- function(loos) {
  Ks <- unlist(lapply(loos, attr, which = "K"))
  if (length(Ks) == length(loos) && !all(Ks == Ks[1])) {
    warning(
      "Not all kfold objects have the same K value. ",
      "For a more accurate comparison use the same number of folds. ",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Warn when k-fold results do not share the same fold assignment
#' @noRd
#' @param loos List of `"kfold"` or `"kfold_pred_measure"` objects.
#' @details The fold labels are arbitrary. The check therefore relabels each
#'   vector by first appearance. Two runs that split the data in the same way
#'   then agree, whatever the labels are. A `NULL` `folds` attribute means the
#'   object does not record the split. The check is then not possible.
throw_kfold_folds_mismatch_warning <- function(loos) {
  folds <- lapply(loos, attr, which = "folds")
  if (any(vapply(folds, is.null, logical(1)))) {
    return(invisible(NULL))
  }
  canonical <- lapply(folds, function(f) {
    as.integer(factor(f, levels = unique(f)))
  })
  same <- vapply(
    canonical,
    function(f) identical(f, canonical[[1L]]),
    logical(1)
  )
  if (!all(same)) {
    warning(
      "Not all kfold objects use the same fold assignment.", call. = FALSE
    )
  }
  invisible(NULL)
}

#' Perform checks on `"loo"` objects before comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param class_check Function returning `TRUE` for valid input objects.
#' @param class_msg Error message when `class_check` fails.
#' @param kfold_checks If `TRUE`, run k-fold comparison warnings.
#' @param n_fun Function returning one model's number of observations. A
#'   `"psis_loo_ss"` object subsamples its `pointwise` matrix, so it reports the
#'   size of the full data instead.
#' @return Nothing, just possibly throws errors/warnings.
model_compare_checks <- function(
  loos,
  class_check = is.loo,
  class_msg = "All inputs should have class 'loo'.",
  kfold_checks = TRUE,
  n_fun = function(x) nrow(x$pointwise)
) {
  ## errors
  if (length(loos) <= 1L) {
    stop("At least two models are required for comparison.", call. = FALSE)
  }
  if (!all(vapply(loos, class_check, logical(1)))) {
    stop(class_msg, call. = FALSE)
  }

  Ns <- vapply(loos, function(x) as.integer(n_fun(x)), integer(1))
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
  yhash_ok <- vapply(yhash, function(x) {
    isTRUE(all.equal(x, yhash[[1]]))
  }, logical(1))
  if (!all(yhash_ok)) {
    warning(
      "Not all models have the same y variable. ('yhash' attributes do not match)",
      call. = FALSE
    )
  }

  if (!kfold_checks) {
    return(invisible(NULL))
  }

  if (all(vapply(loos, is.kfold, logical(1)))) {
    throw_kfold_K_mismatch_warning(loos)
    throw_kfold_folds_mismatch_warning(loos)
  } else if (any(vapply(loos, is.kfold, logical(1))) &&
      any(vapply(loos, is.psis_loo, logical(1)))) {
    warning(
      "Comparing LOO-CV to K-fold-CV. ",
      "For a more accurate comparison use the same number of folds ",
      "or loo for all models compared.",
      call. = FALSE
    )
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

#' Build estimates table for `model_compare()` ordering and matrix output
#' @noRd
.model_compare_estimates_table <- function(loos, bare_names = FALSE,
                                          subsampling = FALSE) {
  sapply(loos, function(x) {
    est <- x$estimates
    rows <- if (bare_names) .display_name(rownames(est), loos) else rownames(est)
    nms <- c(rows, paste0("se_", rows))
    # A `psis_loo_ss` object carries a third estimate column, the subsampling
    # standard error, so its table needs a third name set.
    if (subsampling) {
      nms <- c(nms, paste0("subsampling_se_", rows))
    }
    setNames(c(est), nm = nms)
  })
}

#' Compute the model_compare matrix
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param bare_names If `TRUE`, strip `_loo` suffixes from estimate row names.
#' @param ord Optional model ordering indices; computed from ELPD when `NULL`.
model_compare_matrix <- function(loos, bare_names = FALSE, ord = NULL,
                                 subsampling = FALSE) {
  tmp <- .model_compare_estimates_table(
    loos,
    bare_names = bare_names,
    subsampling = subsampling
  )
  colnames(tmp) <- find_model_names(loos)
  comp <- t(tmp)

  if (is.null(ord)) {
    ord <- model_compare_order(loos)
  }
  comp <- comp[ord, , drop = FALSE]

  patts <- if (bare_names) {
    c("^elpd$", "^p$", "^se_elpd$", "^se_p$")
  } else if (subsampling) {
    # Left unanchored, so each `subsampling_se_*` column is picked up beside
    # its `se_*` counterpart.
    c("elpd", "p_", "^waic$|^looic$", "se_waic$|se_looic$")
  } else {
    c("elpd", "p_", "^waic$|^looic$", "^se_waic$|^se_looic$")
  }
  col_ord <- unique(unlist(
    lapply(patts, function(p) grep(p, colnames(comp))),
    use.names = FALSE
  ))
  if (bare_names) {
    other <- setdiff(seq_len(ncol(comp)), col_ord)
    comp <- comp[, c(col_ord, other), drop = FALSE]
  } else {
    comp <- comp[, col_ord, drop = FALSE]
  }
  comp
}

#' Computes the order of loos for comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param rank_col Optional internal `pointwise` column name used for ranking.
model_compare_order <- function(loos, rank_col = NULL) {
  if (is.null(rank_col)) {
    tmp <- .model_compare_estimates_table(loos, bare_names = FALSE)
    colnames(tmp) <- find_model_names(loos)
    rnms <- rownames(tmp)
    return(order(tmp[grep("^elpd", rnms), ], decreasing = TRUE))
  }

  est_row <- vapply(loos, function(x) {
    val <- x$estimates[rank_col, "Estimate"]
    if (.measure_is_loss(rank_col, loos)) -val else val
  }, numeric(1))
  order(est_row, decreasing = TRUE)
}

#' Perform checks on `"loo"` objects __after__ comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param ord List of `"loo"` object orderings.
#' @param measure_diff Optional precomputed model differences for the rank
#'   measure; computed from the median model when `NULL`.
#' @param rank_col Optional internal `pointwise` column name used for the
#'   median-baseline differences when `measure_diff` is `NULL` and inputs are not
#'   classic `"loo"` objects.
#' @return Nothing, just possibly throws errors/warnings.
model_order_stat_check <- function(loos, ord, measure_diff = NULL, rank_col = NULL) {

  ## breaks

  if (length(loos) <= 11L) {
    # procedure cannot be diagnosed for fewer than ten candidate models
    # (total models = worst model + ten candidates)
    # break from function
    return(invisible(NULL))
  }

  ## warnings

  if (is.null(measure_diff)) {
    # compute differences from the median model
    baseline_idx <- middle_idx(ord)
    ref_loo <- loos[[ord[baseline_idx]]]
    if (is.null(rank_col)) {
      diffs <- mapply(FUN = elpd_diffs, loos[ord[baseline_idx]], loos[ord])
      measure_diff <- apply(diffs, 2, sum)
    } else {
      method <- .measure_pointwise_diff_method(loos, rank_col)
      measure_diff <- vapply(
        loos[ord],
        .pair_measure_stats,
        FUN.VALUE = c(diff = 0, se = 0),
        ref = ref_loo,
        col = rank_col,
        method = method,
        loos = loos
      )["diff", ]
    }
  }

  # estimate the standard deviation of the upper-half-normal
  diff_median <- stats::median(measure_diff)
  measure_diff_trunc <- measure_diff[measure_diff >= diff_median]
  n_models <- sum(!is.na(measure_diff_trunc))
  candidate_sd <- sqrt(1 / n_models * sum(measure_diff_trunc^2, na.rm = TRUE))

  # estimate expected best diff under null hypothesis
  K <- length(loos) - 1
  order_stat <- order_stat_heuristic(K, candidate_sd)

  if (max(measure_diff) <= order_stat) {
    # flag warning if we suspect no model is theoretically better than the baseline
    warning("Difference in performance potentially due to chance. ",
            "See McLatchie and Vehtari (2023) for details.",
            call. = FALSE)
  }
  invisible(NULL)
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
    # The reference model need not be the best one, so a difference can be
    # positive: the flag is about the magnitude, not the sign.
    diag_diff[abs(elpd_diff) < 4 & elpd_diff != 0] <- "|elpd_diff| < 4"
  }
  diag_diff
}
