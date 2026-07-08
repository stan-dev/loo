#' Model comparison
#'
#' @description Compare fitted models based on [ELPD][loo-glossary] or, for
#'   [`loo_pred_measure()`][loo_pred_measure] results, multiple predictive
#'   performance measures at once.
#'
#' @export
#' @param x An object of class `"loo"` or `"loo_pred_measure"`, or a list of
#'   such objects. If a list is used then the list names will be used as the
#'   model names in the output. See **Examples**.
#' @param ... Additional objects of class `"loo"` or `"loo_pred_measure"`, if not
#'   passed in as a single list.
#' @param rank_by For [`loo_pred_measure()`][loo_pred_measure] comparisons only,
#'   the bare measure name used to rank models and define the reference model
#'   for all pairwise differences (default `"elpd"`). For example,
#'   `rank_by = "mse"` ranks models by predictive MSE (best/lowest MSE first) and
#'   computes all measure differences relative to that model on a utility scale
#'   (higher is better; loss measures such as MSE have their sign flipped).
#'
#' @return A data frame with class `"compare.loo"` that has its own
#'   print method. See the **Details** and **Examples** sections.
#'
#'   For classic `"loo"` / `"waic"` / `"kfold"` comparisons, the returned
#'   columns are unchanged from previous versions.
#'
#'   For [`loo_pred_measure()`][loo_pred_measure] comparisons, the data frame
#'   additionally contains `{measure}_diff` and `{measure}_se_diff` columns for
#'   every predictive measure common to all models (e.g. `rmse_diff`,
#'   `rmse_se_diff`). ELPD-family measures use `elpd_diff` and `se_diff`.
#'   `p_worse` and `diag_diff` are computed for ELPD only. Per-model PSIS
#'   diagnostics appear in `diag_elpd`. Attributes `compare_measures` and
#'   `sign_converted_measures` record which measures were compared and which
#'   loss measures had their sign flipped for comparison. Attribute `rank_by` is
#'   set when `rank_by` was passed explicitly (default ranking is by `"elpd"`).
#'   Attribute `measures_no_pointwise_se` lists measures without pointwise-based
#'   `{measure}_se_diff` values.
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in
#'   [`elpd_loo`][loo-glossary] or `elpd_waic` (or multiplied by \eqn{-2}, if
#'   desired, to be on the deviance scale).
#'
#' ## `elpd_diff` and `se_diff`
#'   When using `loo_compare()`, the returned data frame will have one row per
#'   model and several columns of estimates. The values of
#'   [`elpd_diff`][loo-glossary] and [`se_diff`][loo-glossary] are computed by
#'   making pairwise comparisons between each model and the model with the
#'   largest ELPD (the model listed first). Therefore, the first `elpd_diff`
#'   value will always be `0` (i.e., the difference between the preferred model
#'   and itself) and the rest of the values will be negative.
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
#' ## `p_worse`, `diag_diff`, and `diag_elpd`
#'   The values in the `p_worse` column show the probability of each model
#'   having worse ELPD than the best model. These probabilities are computed
#'   with a normal approximation using the values from `elpd_diff` and
#'   `se_diff`. Sivula et al. (2025) present the conditions when the normal
#'   approximation used for SE and `se_diff` is good, and the column
#'   `diag_diff` contains possible diagnostic messages:
#'
#'   * `N < 100` (small data)
#'   * `|elpd_diff| < 4` (models make similar predictions)
#'
#'   If either of these diagnostic messages is shown, the error distribution is
#'   skewed or thick tailed and the normal approximation based on `elpd_diff`
#'   and `se_diff` is not well calibrated. In that case, the probabilities
#'   `p_worse` are likely to be too large. However, `elpd_diff` and `se_diff`
#'   will still be indicative of the differences and uncertainties (for example,
#'   if `|elpd_diff|` is many times larger than `se_diff` the difference is quite
#'   certain). In addition, if the model is not well specificed and there are
#'   outliers, the error distribution can also be skewed or thick tailed and the
#'   normal approximation is not well calibrated. Possible model misspecification
#'   and outliers can be diagnosed with usual predictive checking methods.
#'
#'   The column `diag_elpd` shows the PSIS-LOO Pareto k diagnostic for the
#'   pointwise ELPD computations for each model. If `K k_psis > 0.7` is shown,
#'   where `K` is the number of high Pareto k values in the PSIS
#'   computation, then there may be significant bias in `elpd_diff` favoring
#'   models with a large number of high Pareto k values.
#'
#' ## Comparing `loo_pred_measure` objects
#'   When all inputs are [`loo_pred_measure()`][loo_pred_measure] objects,
#'   `loo_compare()` computes paired differences for every predictive measure
#'   present in all models. Models are ranked by `rank_by` (default `"elpd"`);
#'   the top-ranked model is the reference for all `{measure}_diff` columns.
#'   Measures may use different orientations in their raw form (e.g. ELPD and
#'   CRPS/RPS are returned on a utility scale where higher is better, while MSE
#'   and Brier score are loss measures where lower is better). For comparison,
#'   all `{measure}_diff` values are reported on a common utility scale (higher
#'   is better). Loss measures have their sign flipped from the raw loss
#'   orientation so that negative `{measure}_diff` values indicate worse
#'   performance than the reference model. Each `*_pred_measure()` result stores
#'   attribute `measure_higher_is_better`, a named list recording the
#'   `higher_is_better` setting used when each measure was computed. When loss measures are compared
#'   on a utility scale, `loo_compare()` emits a short message naming the affected
#'   measures, for example:
#'   "For model comparison, differences for mse are reported on a utility scale
#'   (higher is better)."
#'
#'   `p_worse` and `diag_diff` are computed for ELPD-family measures only. Other
#'   measures receive `{measure}_diff` and `{measure}_se_diff` from paired
#'   pointwise contributions when the overall estimate is a sum or mean of those
#'   contributions (using the same standard error formula as `se_diff`). For
#'   measures where pointwise values do not define the overall estimate (e.g.
#'   `r2`, `mse`, `rmse`), `{measure}_diff` is the difference between overall
#'   estimates and `{measure}_se_diff` is `NA`. When models were fit with
#'   different `measure` sets, only measures common to all models are compared; a
#'   warning lists omitted measures. Use `print(x, measures = "all")` to display
#'   diff tables for every compared measure; see [loo-glossary] for column
#'   definitions.
#'
#' ## Warnings for many model comparisons
#'   If more than \eqn{11} models are compared, we internally recompute the model
#'   differences using the median model (by ELPD, or by `rank_by` for
#'   `loo_pred_measure` comparisons) as the baseline model. We then
#'   estimate whether the differences in predictive performance are potentially
#'   due to chance as described by McLatchie and Vehtari (2023). This will flag
#'   a warning if it is deemed that there is a risk of over-fitting due to the
#'   selection process. In that case users are recommended to avoid model
#'   selection based on LOO-CV, and instead to favor model averaging/stacking or
#'   projection predictive inference.
#'
#' @seealso
#' * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'   the __loo__ website for answers to frequently asked questions.
#' @template loo-and-compare-references
#'
#' @examples
#' # very artificial example, just for demonstration!
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)     # should be worst model when compared
#' loo2 <- loo(LL + 1) # should be second best model when compared
#' loo3 <- loo(LL + 2) # should be best model when compared
#'
#' comp <- loo_compare(loo1, loo2, loo3)
#' print(comp, digits = 2)
#'
#' # can use a list of objects with custom names
#' # the names will be used in the output
#' loo_compare(list("apple" = loo1, "banana" = loo2, "cherry" = loo3))
#'
#' \dontrun{
#' # works for waic (and kfold) too
#' loo_compare(waic(LL), waic(LL - 10))
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
#'     measure = c("rmse", "r2")
#'   )
#'   pm2 <- loo_pred_measure(
#'     loo = loo(fit2, save_psis = TRUE),
#'     y = fit2$data$Reaction,
#'     mupred = brms::posterior_epred(fit2),
#'     measure = c("rmse", "r2")
#'   )
#'   comp <- loo_compare(pm1, pm2)
#'   print(comp)                      # ranked by elpd (default)
#'   print(comp, measures = "all")    # all measure diff tables
#'   loo_compare(pm1, pm2, rank_by = "rmse")
#' }
#' }
#'
loo_compare <- function(x, ..., rank_by = NULL) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ..., rank_by = NULL) {
  loos <- .loo_compare_inputs(x, ...)

  # if subsampling is used
  if (any(sapply(loos, inherits, "psis_loo_ss"))) {
    return(loo_compare.psis_loo_ss_list(loos))
  }

  if (all(vapply(loos, is.loo_pred_measure, logical(1)))) {
    return(compare_loo_pred_measure(loos, rank_by = rank_by))
  }

  if (any(vapply(loos, is.loo_pred_measure, logical(1)))) {
    stop(
      "Cannot mix 'loo_pred_measure' objects with other 'loo' objects. ",
      "Compare models using 'loo_pred_measure()' for each model.",
      call. = FALSE
    )
  }

  if (any(vapply(loos, inherits, what = "pred_measure", logical(1)))) {
    stop(
      "'loo_compare' for predictive measures requires 'loo_pred_measure' objects. ",
      "Use loo_pred_measure() instead of insample_pred_measure(), ",
      "kfold_pred_measure(), or test_pred_measure().",
      call. = FALSE
    )
  }

  if (!is.null(rank_by)) {
    warning(
      "`rank_by` is only used for `loo_pred_measure` comparisons and will be ignored.",
      call. = FALSE
    )
  }

  # run pre-comparison checks
  loo_compare_checks(loos)

  # compute elpd_diff and se_elpd_diff relative to best model
  ord <- loo_compare_order(loos)
  comp <- loo_compare_matrix(loos, ord = ord)
  rnms <- rownames(comp)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[1]], loos[ord])
  colnames(diffs) <- rnms
  elpd_diff <- apply(diffs, 2, sum)
  se_diff <- apply(diffs, 2, se_elpd_diff)

  # compute probabilities that a model has worse elpd than the best model
  # using a normal approximation (Sivula et al., 2025)
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
  loo_order_stat_check(loos, ord)

  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' @rdname loo_compare
#' @export
#' @param digits For the print method only, the number of digits to use when
#'   printing.
#' @param p_worse For the print method only, should we include the normal
#'   approximation based probability of each model having worse performance than
#'   the best model? The default is `TRUE`.
#' @param measures For `loo_pred_measure` comparisons only, which measures to
#'   print diff tables for. `NULL` (default) prints only the ranking measure
#'   (`"elpd"` when `rank_by` was not set, otherwise `rank_by`);
#'   `"all"` prints all compared measures; or a character vector of measure
#'   names (e.g. `c("elpd", "mse")`). Printed tables use the column name
#'   `se_diff` for the standard error of the difference even when the data
#'   frame column is `{measure}_se_diff`.
print.compare.loo <- function(x, ..., digits = 1, p_worse = TRUE, measures = NULL) {
  if (inherits(x, "old_compare.loo")) {
    return(unclass(x))
  }
  if (!inherits(x, "data.frame")) {
    class(x) <- c(class(x), "data.frame")
  }

  compare_measures <- attr(x, "compare_measures")
  if (!is.null(compare_measures)) {
    return(.print_compare_loo_pred_measure(
      x,
      digits = digits,
      p_worse = p_worse,
      measures = measures
    ))
  }

  if (!all(c("model", "elpd_diff", "se_diff") %in% colnames(x))) {
    print(as.data.frame(x))
    return(x)
  }
  x2 <- cbind(
    model = x$model,
    .fr(x[, c("elpd_diff", "se_diff")], digits)
  )
  if (p_worse && "p_worse" %in% colnames(x)) {
    x2 <- cbind(
      x2,
      p_worse = .fr(x[, "p_worse"], digits = 2),
      diag_diff = x[, "diag_diff"],
      diag_elpd = x[, "diag_elpd"]
    )
  }
  print(x2, quote = FALSE, row.names = FALSE)

  .print_compare_diag_message(x, p_worse = p_worse)
  invisible(x)
}


# internal ----------------------------------------------------------------

#' Print `compare.loo` results from `loo_pred_measure` comparisons
#' @noRd
.print_compare_loo_pred_measure <- function(x, digits, p_worse, measures) {
  rank_by <- attr(x, "rank_by")
  compare_measures <- attr(x, "compare_measures")
  ref_model <- x$model[[1L]]
  primary_measure <- if (is.null(rank_by)) "elpd" else rank_by

  measures_to_print <- if (is.null(measures)) {
    primary_measure
  } else if (identical(measures, "all")) {
    compare_measures
  } else {
    measures
  }

  unknown <- setdiff(measures_to_print, compare_measures)
  if (length(unknown)) {
    stop(
      paste0(
        "Unknown measure(s) in `measures`: ",
        paste(unknown, collapse = ", "),
        ". Available measures: ",
        paste(compare_measures, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (identical(measures, "all") && length(compare_measures) > 4L) {
    message(
      "Printing ", length(compare_measures), " measure comparisons; ",
      "consider `measures = c(...)`."
    )
  }

  if (is.null(measures) && !is.null(rank_by)) {
    message(
      "Models ranked by ", rank_by, " (reference: ", ref_model, ")."
    )
  }

  show_diag_elpd_primary <- is.null(measures) || identical(measures, "all")
  for (i in seq_along(measures_to_print)) {
    measure <- measures_to_print[[i]]
    if (!is.null(measures)) {
      cat("\n-- ", measure, " (vs ", ref_model, ") --\n", sep = "")
    }
    .print_compare_measure_table(
      x,
      measure = measure,
      digits = digits,
      p_worse = p_worse,
      show_diag_elpd = show_diag_elpd_primary &&
        identical(measure, primary_measure) &&
        i == match(primary_measure, measures_to_print)
    )
  }

  has_diag_msg <- .print_compare_diag_message(
    x,
    p_worse = p_worse,
    measures = measures_to_print
  )

  if (is.null(measures)) {
    other <- setdiff(compare_measures, primary_measure)
    if (length(other)) {
      message(
        if (has_diag_msg) "\n",
        "Other measures compared: ",
        paste(other, collapse = ", "),
        ". Use print(x, measures = \"all\")."
      )
    }
  }

  .warn_measures_no_pointwise_se(attr(x, "measures_no_pointwise_se"))

  invisible(x)
}

#' Print one measure's comparison table
#' @noRd
.print_compare_measure_table <- function(x, measure, digits, p_worse, show_diag_elpd) {
  if (.is_elpd_measure(measure)) {
    diff_col <- "elpd_diff"
    se_col <- "se_diff"
    diff_name <- "elpd_diff"
    se_name <- "se_diff"
  } else {
    diff_col <- paste0(measure, "_diff")
    se_col <- paste0(measure, "_se_diff")
    diff_name <- diff_col
    se_name <- "se_diff"
  }

  if (!all(c(diff_col, se_col) %in% colnames(x))) {
    stop(
      "Comparison columns for measure '", measure, "' are missing.",
      call. = FALSE
    )
  }

  x2 <- data.frame(
    model = x$model,
    diff = unname(.fr(x[[diff_col]], digits)),
    se_diff = unname(.fr(x[[se_col]], digits)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(x2)[2:3] <- c(diff_name, se_name)

  if (.is_elpd_measure(measure) && p_worse && "p_worse" %in% colnames(x)) {
    x2$p_worse <- unname(.fr(x[["p_worse"]], digits = 2))
    x2$diag_diff <- x[["diag_diff"]]
  }
  if (show_diag_elpd && "diag_elpd" %in% colnames(x)) {
    x2$diag_elpd <- x[["diag_elpd"]]
  }

  print(x2, quote = FALSE, row.names = FALSE)
}

#' Print diagnostic glossary message for compare output
#' @noRd
.print_compare_diag_message <- function(x, p_worse, measures = NULL) {
  diag_cols <- c("diag_elpd")
  if (is.null(measures) || "elpd" %in% measures) {
    diag_cols <- c("diag_diff", diag_cols)
  } else if (!is.null(measures)) {
    elpd_measures <- measures[vapply(measures, .is_elpd_measure, logical(1))]
    if (length(elpd_measures)) {
      diag_cols <- c("diag_diff", diag_cols)
    }
  }

  has_diag <- any(
    vapply(
      intersect(diag_cols, colnames(x)),
      function(col) any(nzchar(x[[col]], keepNA = FALSE), na.rm = TRUE),
      logical(1)
    )
  )
  if (has_diag && p_worse) {
    message(
      "\nDiagnostic flags present.\n",
      "See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)\n",
      "or https://mc-stan.org/loo/reference/loo-glossary.html."
    )
  }
  invisible(has_diag && p_worse)
}

#' Is an object a PSIS-LOO predictive measure result?
#' @noRd
is.loo_pred_measure <- function(x) {
  inherits(x, "loo_pred_measure")
}

#' Normalize `loo_compare()` inputs to a list of model results
#' @noRd
.loo_compare_inputs <- function(x, ...) {
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

#' Compare `loo_pred_measure` objects (multi-measure path)
#' @noRd
#' @param loos List of `loo_pred_measure` objects.
#' @param rank_by Bare measure name used to order models.
compare_loo_pred_measure <- function(loos, rank_by = NULL) {
  loo_compare_checks(
    loos,
    class_check = is.loo_pred_measure,
    class_msg = "All inputs must have class 'loo_pred_measure'.",
    kfold_checks = FALSE
  )
  .compare_metadata_check(loos)
  .warn_omitted_compare_measures(loos)

  rank_measure <- .resolve_rank_measure(loos, rank_by)
  compare_cols <- .compare_pointwise_cols(loos)
  .inform_compare_sign_conversion(compare_cols, loos)
  ord <- loo_compare_order(loos, rank_measure$internal)
  loos_ord <- loos[ord]
  ref_loo <- loos_ord[[1L]]

  comp <- loo_compare_matrix(
    loos_ord,
    bare_names = TRUE,
    ord = seq_along(loos_ord)
  )
  rnms <- rownames(comp)
  n_obs <- nrow(loos_ord[[1L]]$pointwise)

  diff_cols <- list()
  measures_no_pointwise_se <- list()
  for (col in compare_cols) {
    bare <- .display_name(col)
    method <- .measure_pointwise_diff_method(loos_ord, col)
    if (method == "estimates_only") {
      measures_no_pointwise_se[[bare]] <- bare
    }
    pair_stats <- vapply(
      loos_ord,
      .pair_measure_stats,
      FUN.VALUE = c(diff = 0, se = 0),
      ref = ref_loo,
      col = col,
      method = method,
      loos = loos_ord
    )
    measure_diff <- pair_stats["diff", ]
    measure_se <- pair_stats["se", ]

    if (.is_elpd_measure(col)) {
      diff_cols$elpd_diff <- measure_diff
      diff_cols$se_diff <- measure_se
      p_worse <- stats::pnorm(0, measure_diff, measure_se)
      p_worse[measure_diff == 0] <- NA_real_
      diff_cols$p_worse <- p_worse
      diff_cols$diag_diff <- diag_diff(n_obs, measure_diff)
    } else {
      diff_cols[[paste0(bare, "_diff")]] <- measure_diff
      diff_cols[[paste0(bare, "_se_diff")]] <- measure_se
    }
  }

  comp <- cbind(
    data.frame(
      model = rnms,
      diff_cols,
      diag_elpd = diag_elpd(loos_ord),
      stringsAsFactors = FALSE
    ),
    as.data.frame(comp)
  )
  rownames(comp) <- NULL

  loo_order_stat_check(
    loos_ord,
    seq_along(loos_ord),
    rank_col = rank_measure$internal
  )

  attr(comp, "measures_no_pointwise_se") <- unique(unlist(measures_no_pointwise_se))
  if (!is.null(rank_by)) {
    attr(comp, "rank_by") <- rank_measure$bare
  }
  attr(comp, "compare_measures") <- .compare_measures(loos)
  attr(comp, "sign_converted_measures") <- .compare_sign_converted_measures(
    compare_cols,
    loos
  )
  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' Strip `_loo` suffix for `loo_compare` display names
#' @noRd
.display_name <- function(col) {
  sub("_loo$", "", col)
}

#' Map bare measure name to `pointwise` column name
#' @noRd
.pointwise_col <- function(name, cols) {
  if (name %in% cols) {
    return(name)
  }
  internal <- paste0(name, "_loo")
  if (internal %in% cols) {
    return(internal)
  }
  stop(
    paste0(
      "Measure '", name, "' not found in all models. ",
      "Available measures: ",
      paste(vapply(cols, .display_name, character(1)), collapse = ", ")
    ),
    call. = FALSE
  )
}

#' Common `pointwise` columns across models, excluding complexity terms
#' @noRd
.compare_pointwise_cols <- function(loos) {
  cols <- Reduce(
    intersect,
    lapply(loos, function(x) colnames(x$pointwise))
  )
  cols[!grepl("^p_", cols)]
}

#' Check that comparison metadata is consistent across models
#' @noRd
.compare_metadata_check <- function(loos) {
  bare_measures <- .compare_measures(loos)
  if (!length(bare_measures)) {
    return(invisible(NULL))
  }

  for (bare in bare_measures) {
    metas <- lapply(loos, function(x) {
      compare_meta <- attr(x, "measure_compare_meta")
      if (is.null(compare_meta)) {
        return(NULL)
      }
      compare_meta[[bare]]
    })
    has_meta <- !vapply(metas, is.null, logical(1))
    if (any(has_meta) && !all(has_meta)) {
      stop(
        "Not all models provide comparison metadata for measure '",
        bare,
        "'. Recompute all inputs with the current version of `loo_pred_measure()`.",
        call. = FALSE
      )
    }
    non_null <- metas[has_meta]
    if (length(non_null) > 1L) {
      ref <- non_null[[1L]]
      inconsistent <- vapply(
        non_null[-1L],
        function(meta) !identical(meta, ref),
        logical(1)
      )
      if (any(inconsistent)) {
        stop(
          "Models disagree on comparison metadata for measure '",
          bare,
          "'. Ensure all models use the same `higher_is_better` settings for each measure.",
          call. = FALSE
        )
      }
    }
  }

  invisible(NULL)
}

#' Warn when models do not share the same predictive measures
#' @noRd
.warn_omitted_compare_measures <- function(loos) {
  model_names <- find_model_names(loos)
  if (anyDuplicated(model_names)) {
    model_names <- make.unique(model_names, sep = "_")
  }
  by_model <- stats::setNames(
    lapply(loos, function(x) {
      cols <- colnames(x$pointwise)
      cols <- cols[!grepl("^p_", cols)]
      unname(vapply(cols, .display_name, character(1)))
    }),
    model_names
  )
  common <- Reduce(intersect, by_model)
  omitted <- setdiff(unique(unlist(by_model)), common)
  if (!length(omitted)) {
    return(invisible(NULL))
  }
  omitted <- sort(omitted)

  omitted_detail <- vapply(
    omitted,
    function(measure) {
      present <- names(by_model)[vapply(
        by_model,
        function(measures) measure %in% measures,
        logical(1)
      )]
      paste0(measure, " (", paste(present, collapse = ", "), ")")
    },
    character(1)
  )

  warning(
    paste0(
      "Omitted measures: ",
      paste(omitted_detail, collapse = ", "),
      ". Compared: ",
      paste(common, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

#' Warn when `se_diff` is unavailable for compared measures
#' @noRd
.warn_measures_no_pointwise_se <- function(measures) {
  if (!length(measures)) {
    return(invisible(NULL))
  }
  warning(
    paste0(
      "se_diff unavailable for: ",
      paste(measures, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
  invisible(NULL)
}

#' Bare measure names available for comparison across models
#' @noRd
.compare_measures <- function(loos) {
  cols <- .compare_pointwise_cols(loos)
  unname(vapply(cols, .display_name, character(1)))
}

#' Resolve `rank_by` to bare and internal `pointwise` column names
#' @noRd
.resolve_rank_measure <- function(loos, rank_by = NULL) {
  cols <- .compare_pointwise_cols(loos)
  bare <- if (is.null(rank_by)) "elpd" else rank_by
  internal <- .pointwise_col(bare, cols)
  list(
    bare = .display_name(internal),
    internal = internal
  )
}

#' Is a measure an ELPD-family measure (for `p_worse` / `diag_diff`)?
#' @noRd
.is_elpd_measure <- function(name) {
  grepl("^elpd", .display_name(name))
}

#' Look up per-measure comparison metadata on a result object
#' @noRd
.get_measure_compare_meta <- function(loos, bare) {
  compare_meta <- attr(loos[[1L]], "measure_compare_meta")
  if (is.null(compare_meta)) {
    return(NULL)
  }
  compare_meta[[bare]]
}

#' Whether stored values are on a loss scale (lower is better)
#' @noRd
.measure_lower_is_better <- function(name, loos = NULL) {
  bare <- .display_name(name)
  higher_is_better <- NULL
  loss <- NULL

  if (!is.null(loos)) {
    meta <- .get_measure_compare_meta(loos, bare)
    if (!is.null(meta)) {
      higher_is_better <- meta$higher_is_better
      loss <- meta$loss
    } else {
      hib_attr <- attr(loos[[1L]], "measure_higher_is_better")
      if (!is.null(hib_attr) && bare %in% names(hib_attr)) {
        higher_is_better <- hib_attr[[bare]]
      }
    }
  }

  if (!is.null(higher_is_better)) {
    return(!isTRUE(higher_is_better))
  }

  if (is.null(loss)) {
    spec <- .measure_spec[[bare]]
    loss <- if (!is.null(spec)) {
      isTRUE(spec$loss)
    } else {
      bare %in% c("ic", "mae", "mse", "rmse", "brier", "srps")
    }
  }

  isTRUE(loss)
}

#' Bare names of measures whose sign is flipped for `loo_compare()`
#' @noRd
.compare_sign_converted_measures <- function(cols, loos) {
  bare <- vapply(cols, .display_name, character(1))
  unique(bare[vapply(
    cols,
    function(col) .measure_lower_is_better(col, loos),
    logical(1)
  )])
}

#' Inform when measure signs are flipped for comparison
#' @noRd
.inform_compare_sign_conversion <- function(cols, loos) {
  converted <- .compare_sign_converted_measures(cols, loos)
  if (!length(converted)) {
    return(invisible(NULL))
  }
  message(
    "For model comparison, differences for ",
    paste(converted, collapse = ", "),
    " ",
    if (length(converted) == 1L) "is" else "are",
    " reported on a utility scale (higher is better)."
  )
  invisible(NULL)
}

#' How to aggregate paired pointwise differences for a measure
#'
#' Returns `"sum"` when the overall estimate equals the sum of pointwise
#' contributions, `"mean"` when it equals the mean, and `"estimates_only"`
#' when pointwise values do not define the overall estimate.
#' @noRd
.measure_pointwise_diff_method <- function(loos, col) {
  bare <- .display_name(col)
  meta <- .get_measure_compare_meta(loos, bare)
  if (!is.null(meta) && !identical(meta$diff_method, "auto")) {
    return(meta$diff_method)
  }

  spec <- .measure_spec[[bare]]
  if (!is.null(spec) && identical(spec$diff_method, "estimates_only")) {
    return("estimates_only")
  }
  if (.is_elpd_measure(col) || bare == "ic") {
    return("sum")
  }

  ref <- loos[[1L]]
  est <- ref$estimates[col, "Estimate"]
  pw <- ref$pointwise[, col, drop = TRUE]
  if (!length(pw) || !is.finite(est)) {
    return("estimates_only")
  }

  tol <- sqrt(.Machine$double.eps) * max(abs(c(est, pw)), na.rm = TRUE)
  if (isTRUE(all.equal(est, sum(pw), tolerance = tol, check.attributes = FALSE))) {
    return("sum")
  }
  if (isTRUE(all.equal(est, mean(pw), tolerance = tol, check.attributes = FALSE))) {
    return("mean")
  }
  "estimates_only"
}

#' Paired measure difference and SE for one model vs a reference
#' @noRd
.pair_measure_stats <- function(cmp, ref, col, method = NULL, loos = list(ref)) {
  if (is.null(method)) {
    method <- .measure_pointwise_diff_method(c(list(ref, cmp)), col)
  }

  flip <- .measure_lower_is_better(col, loos)

  if (method == "estimates_only") {
    est_utility <- function(estimates) {
      val <- estimates[col, "Estimate"]
      if (flip) -val else val
    }
    return(c(
      diff = est_utility(cmp$estimates) - est_utility(ref$estimates),
      se = NA_real_
    ))
  }

  to_utility <- function(pointwise) {
    x <- pointwise[, col, drop = TRUE]
    if (flip) -x else x
  }
  diffs <- to_utility(cmp$pointwise) - to_utility(ref$pointwise)

  diff <- if (method == "sum") sum(diffs) else mean(diffs)
  se <- if (method == "sum") {
    se_elpd_diff(diffs)
  } else {
    N <- length(diffs)
    if (N <= 1L) 0 else stats::sd(diffs) / sqrt(N)
  }
  c(diff = diff, se = se)
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

#' Perform checks on `"loo"` objects before comparison
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param class_check Function returning `TRUE` for valid input objects.
#' @param class_msg Error message when `class_check` fails.
#' @param kfold_checks If `TRUE`, run k-fold comparison warnings.
#' @return Nothing, just possibly throws errors/warnings.
loo_compare_checks <- function(
  loos,
  class_check = is.loo,
  class_msg = "All inputs should have class 'loo'.",
  kfold_checks = TRUE
) {
  ## errors
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.", call. = FALSE)
  }
  if (!all(vapply(loos, class_check, logical(1)))) {
    stop(class_msg, call. = FALSE)
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
    Ks <- unlist(lapply(loos, attr, which = "K"))
    if (!all(Ks == Ks[1])) {
      warning(
        "Not all kfold objects have the same K value. ",
        "For a more accurate comparison use the same number of folds. ",
        call. = FALSE
      )
    }
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


#' Build estimates table for `loo_compare()` ordering and matrix output
#' @noRd
.loo_compare_estimates_table <- function(loos, bare_names = FALSE) {
  sapply(loos, function(x) {
    est <- x$estimates
    rows <- if (bare_names) .display_name(rownames(est)) else rownames(est)
    setNames(c(est), nm = c(rows, paste0("se_", rows)))
  })
}

#' Compute the loo_compare matrix
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param bare_names If `TRUE`, strip `_loo` suffixes from estimate row names.
#' @param ord Optional model ordering indices; computed from ELPD when `NULL`.
loo_compare_matrix <- function(loos, bare_names = FALSE, ord = NULL) {
  tmp <- .loo_compare_estimates_table(loos, bare_names = bare_names)
  colnames(tmp) <- find_model_names(loos)
  comp <- t(tmp)

  if (is.null(ord)) {
    ord <- loo_compare_order(loos)
  }
  comp <- comp[ord, , drop = FALSE]

  patts <- if (bare_names) {
    c("^elpd$", "^p$", "^se_elpd$", "^se_p$")
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
loo_compare_order <- function(loos, rank_col = NULL) {
  if (is.null(rank_col)) {
    tmp <- .loo_compare_estimates_table(loos, bare_names = FALSE)
    colnames(tmp) <- find_model_names(loos)
    rnms <- rownames(tmp)
    return(order(tmp[grep("^elpd", rnms), ], decreasing = TRUE))
  }

  est_row <- vapply(loos, function(x) {
    val <- x$estimates[rank_col, "Estimate"]
    if (.measure_lower_is_better(rank_col, loos)) -val else val
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
loo_order_stat_check <- function(loos, ord, measure_diff = NULL, rank_col = NULL) {

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
    diag_diff[elpd_diff > -4 & elpd_diff != 0] <- "|elpd_diff| < 4"
  }
  diag_diff
}
