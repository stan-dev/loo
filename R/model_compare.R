#' Model comparison
#'
#' @description Compare fitted models based on [ELPD][loo-glossary] or, for
#'   [`pred_measure`][pred_measure] results, multiple predictive performance
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
#'   All models in one call must be evaluated the same way: every input must
#'   come from the same `*_pred_measure()` function, since paired differences
#'   between, say, a LOO and a k-fold result would contrast held-out schemes
#'   rather than models. Mixed inputs are an error.
#'
#' @export
#' @param x An object of class `"loo"` or `"pred_measure"`, or a list of
#'   such objects. If a list is used then the list names will be used as the
#'   model names in the output. See **Examples**.
#' @param ... Additional objects of class `"loo"` or `"pred_measure"`, if not
#'   passed in as a single list.
#' @param rank_by A single string naming either a **measure** or a **model**,
#'   used to define one reference model for all pairwise differences.
#'
#'   A **measure name** ([`pred_measure`][pred_measure] comparisons only) ranks
#'   models by that measure and makes the top-ranked model the reference. Bare
#'   names are used regardless of source, so `rank_by = "rmse"` selects
#'   `rmse_loo`, `rmse_kfold`, or `rmse_test` as appropriate. For example,
#'   `rank_by = "mse"` ranks models by predictive MSE (best/lowest MSE first)
#'   and computes *all* measure differences relative to that one model on a
#'   utility scale (higher is better; loss measures such as MSE have their sign
#'   flipped).
#'
#'   A **model name** (one of the names shown in the `model` column, i.e. the
#'   list names or `model1`, `model2`, ...) pins that model as the reference for
#'   all differences, whichever model performs best. Rows stay ordered by
#'   `"elpd"`. This form also works for plain `"loo"` comparisons, where
#'   `elpd_diff` is then relative to the named model rather than to the best
#'   one. If a name matches both a measure and a model, the measure wins and a
#'   warning is issued.
#'
#'   When `rank_by` is `NULL` (the default), rows are ordered by `"elpd"` but
#'   each measure is compared against *its own* best model, so `mse_diff` may be
#'   relative to a different model than `elpd_diff`. Each `{measure}_diff`
#'   column then has exactly one `0` entry, at that measure's best model.
#' @param custom_se_fn How to compute the standard error of the difference
#'   between two models for a **custom** measure. Required whenever a custom
#'   measure is compared; nothing is inferred from the measure's values. One of:
#'   \itemize{
#'     \item a **function** called as `custom_se_fn(ref, cmp)` (see
#'       **Custom measure standard errors** below);
#'     \item `"sum"`, for a measure whose estimate is the sum of its pointwise
#'       values, giving `sqrt(N) * sd(d_i)` as for `elpd`;
#'     \item `"mean"`, for a measure whose estimate is the mean of its pointwise
#'       values, giving `sd(d_i) / sqrt(N)` as for `mae`;
#'     \item `NULL`, to report the difference with an `NA` standard error.
#'   }
#'   When two or more custom measures are compared, pass a list named by bare
#'   measure name, e.g. `list(huber = "mean", nrmse = my_se_fn)`. Ignored, with
#'   a warning, when no custom measure is present.
#'
#' @section Custom measure standard errors:
#'   A function passed as `custom_se_fn` is called once per comparison as
#'   `custom_se_fn(ref = <list>, cmp = <list>)`, with **named** arguments. Each
#'   argument describes one model and has elements `estimate` (scalar), `se`
#'   (that model's own standard error), `pointwise` (a plain numeric vector, not
#'   a matrix), and `extra` (whatever the measure returned as `extra`, or
#'   `NULL`). All values are on the measure's natural scale, so the function
#'   does not need to account for the utility-scale conversion applied to the
#'   reported differences. It must return the standard error of the difference
#'   as a numeric scalar. For example:
#'
#'   ```
#'   my_se_fn <- function(ref, cmp) {
#'     d <- cmp$pointwise - ref$pointwise
#'     sd(d) / sqrt(length(d))
#'   }
#'   ```
#'
#' @return A data frame with class `"compare.loo"` that has its own
#'   print method. See the **Details** and **Examples** sections.
#'
#'   For classic `"loo"` / `"waic"` / `"kfold"` comparisons, the returned
#'   columns are unchanged from previous versions.
#'
#'   For [`pred_measure`][pred_measure] comparisons, the data frame
#'   additionally contains `{measure}_diff` and `{measure}_se_diff` columns for
#'   every predictive measure common to all models (e.g. `rmse_diff`,
#'   `rmse_se_diff`). ELPD-family measures use `elpd_diff` and `se_diff`.
#'   `p_worse` and `diag_diff` are computed for ELPD only. `diag_elpd` holds
#'   per-model PSIS diagnostics and is present only for
#'   [`loo_pred_measure()`][loo_pred_measure] comparisons, the only source with
#'   Pareto \eqn{\hat{k}} values. Attributes `compare_measures` and
#'   `sign_converted_measures` record which measures were compared and which
#'   loss measures had their sign flipped for comparison. Attribute
#'   `compare_source` records the shared evaluation source (`"loo"`,
#'   `"kfold"`, `"test"`, or `"insample"`). Attribute `rank_by` is
#'   set when `rank_by` named a measure (default ranking is by `"elpd"`), and
#'   attribute `compare_ref_model` is set when it named a model.
#'   Attribute `compare_reference` is a named character vector giving the
#'   reference model each measure's differences were computed against; all
#'   entries are that single reference model when `rank_by` was supplied.
#'
#' @details
#'   When comparing two fitted models, we can estimate the difference in their
#'   expected predictive accuracy by the difference in
#'   [`elpd_loo`][loo-glossary] or `elpd_waic` (or multiplied by \eqn{-2}, if
#'   desired, to be on the deviance scale).
#'
#' ## `elpd_diff` and `se_diff`
#'   When using `model_compare()`, the returned data frame will have one row per
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
#'   Pareto \eqn{\hat{k}} describes a model's PSIS-LOO approximation rather than
#'   any one measure or any one pair of models, and every LOO measure is computed
#'   from the same importance weights. For `pred_measure` comparisons `print()`
#'   therefore reports it once per model, in a block above the difference tables,
#'   instead of as a column inside one of them. The `diag_elpd` column is still
#'   returned on the object.
#'
#' ## Comparing `pred_measure` objects
#'   When all inputs are predictive measure results sharing one evaluation
#'   source,
#'   `model_compare()` computes paired differences for every predictive measure
#'   present in all models. Measures are matched on their bare names, so the
#'   source suffix (`_loo`, `_kfold`, `_test`, or none for in-sample) is
#'   handled transparently. Rows are ordered by `rank_by` (default `"elpd"`).
#'   By default each measure is compared against the model that is best on that
#'   measure, so `mse_diff` can use a different reference model than
#'   `elpd_diff`; the reference used for each measure is recorded in attribute
#'   `compare_reference` and shown by `print(x, measures = "all")`. Supplying
#'   `rank_by` instead pins a single reference --- the top-ranked model --- for
#'   every `{measure}_diff` column. The returned data frame carries one row
#'   order for all measures, but each *printed* measure table is sorted by its
#'   own difference, so the best model on that measure is always the first row
#'   and the differences run in decreasing order.
#'   Measures may use different orientations in their raw form (e.g. ELPD and
#'   SRPS/SCRPS are returned on a utility scale where higher is better, while
#'   MSE, RPS/CRPS and the Brier score are loss measures where lower is better).
#'   For comparison, all `{measure}_diff` values are reported on a common
#'   utility scale (higher is better). Loss measures have their sign flipped
#'   from the raw loss orientation so that negative `{measure}_diff` values
#'   indicate worse performance than the reference model. Which measures are
#'   losses is recorded in the `loss` element of each measure's entry in the
#'   `measure_info` attribute of an `*_pred_measure()` result. When loss
#'   measures are compared on a utility scale, `model_compare()` emits a short
#'   message naming the affected measures, for example:
#'   "For model comparison, differences for mse are reported on a utility scale
#'   (higher is better)."
#'
#'   A custom measure is treated as a utility unless it declares otherwise with
#'   `attr(my_fun, "measure_loss") <- TRUE`. The declaration also determines the
#'   direction of `rank_by`, so an undeclared loss is both flipped and ranked in
#'   the wrong direction; see [insample_pred_measure()].
#'
#'   `p_worse` and `diag_diff` are computed for ELPD-family measures only. Other
#'   measures receive `{measure}_diff` and `{measure}_se_diff`. How
#'   `{measure}_se_diff` is obtained is recorded in the `diff_method` element of
#'   the measure's entry in `measure_info`:
#'
#'   * `"sum"` or `"mean"`: the overall estimate is the sum (`elpd`, `ic`) or the
#'     mean (`mlpd`, `mae`, `mse`, `acc`, `rps`, `srps`, `brier`) of its
#'     pointwise contributions, so the standard error is computed from paired
#'     pointwise differences (the same formula as `se_diff`).
#'   * `"measure_specific"`: the overall estimate is not a sum or mean of
#'     pointwise contributions (`r2`, `rmse`, `bacc`), so the measure supplies
#'     its own standard error of the difference.
#'   * `"custom"`: a custom measure declares nothing, so `custom_se_fn` must be
#'     supplied. `{measure}_se_diff` is `NA` only when `custom_se_fn` is an
#'     explicit `NULL` for that measure.
#'
#'   Objects carrying no `measure_info` at all fall back to the difference
#'   between overall estimates with `{measure}_se_diff` set to `NA`.
#'   When models were fit with
#'   different `measure` sets, only measures common to all models are compared; a
#'   warning lists omitted measures. Use `print(x, measures = "all")` to display
#'   diff tables for every compared measure; see [loo-glossary] for column
#'   definitions.
#'
#' ## Source-specific behavior
#'   Comparisons behave the same way across sources, with three exceptions:
#'
#'   * **`diag_elpd`** is only produced for
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
#'     measure = c("rmse", "r2")
#'   )
#'   pm2 <- loo_pred_measure(
#'     loo = loo(fit2, save_psis = TRUE),
#'     y = fit2$data$Reaction,
#'     mupred = brms::posterior_epred(fit2),
#'     measure = c("rmse", "r2")
#'   )
#'   comp <- model_compare(pm1, pm2)
#'   print(comp)                      # ranked by elpd (default)
#'   print(comp, measures = "all")    # all measure diff tables
#'   model_compare(pm1, pm2, rank_by = "rmse")
#'
#'   # `rank_by` also takes a model name: every measure is then compared
#'   # against that model, whether or not it is the best one
#'   model_compare(list(m1 = pm1, m2 = pm2), rank_by = "m1")
#'
#'   # the same works for k-fold CV; `rank_by` still takes the bare name
#'   # even though the measures are stored as `elpd_kfold`, `rmse_kfold`, ...
#'   kf1 <- brms::kfold(fit1, K = 5, save_fits = TRUE)
#'   kf2 <- brms::kfold(fit2, K = 5, save_fits = TRUE)
#'   kpm1 <- kfold_pred_measure(
#'     y = fit1$data$Reaction,
#'     mupred = brms::kfold_predict(kf1, method = "fitted")$yrep,
#'     kfold = kf1,
#'     measure = "rmse"
#'   )
#'   kpm2 <- kfold_pred_measure(
#'     y = fit2$data$Reaction,
#'     mupred = brms::kfold_predict(kf2, method = "fitted")$yrep,
#'     kfold = kf2,
#'     measure = "rmse"
#'   )
#'   model_compare(kpm1, kpm2, rank_by = "rmse")
#'
#'   # mixing evaluation sources is an error
#'   try(model_compare(pm1, kpm2))
#' }
#' }
#'
model_compare <- function(x, ..., rank_by = NULL, custom_se_fn) {
  UseMethod("model_compare")
}

#' @rdname model_compare
#' @export
model_compare.default <- function(x, ..., rank_by = NULL, custom_se_fn) {
  # `custom_se_fn` is deliberately given no default: an omitted argument and an
  # explicit `NULL` mean different things (error vs. "report an NA se_diff").
  custom_se_fn_supplied <- !missing(custom_se_fn)
  if (!custom_se_fn_supplied) {
    custom_se_fn <- NULL
  }

  loos <- .model_compare_inputs(x, ...)

  # if subsampling is used
  if (any(sapply(loos, inherits, "psis_loo_ss"))) {
    if (custom_se_fn_supplied) {
      stop(
        "`custom_se_fn` is not supported for subsampled loo objects, which ",
        "are compared on elpd only.",
        call. = FALSE
      )
    }
    return(model_compare.psis_loo_ss_list(loos))
  }

  # `pred_measure` objects must be tested before any `is.loo()` check: results
  # from `loo_pred_measure()` and `kfold_pred_measure()` inherit the classes of
  # the `loo`/`kfold` object they were built from.
  is_pm <- vapply(loos, is.pred_measure, logical(1))

  if (all(is_pm)) {
    return(compare_pred_measure(
      loos,
      rank_by = rank_by,
      custom_se_fn = custom_se_fn,
      custom_se_fn_supplied = custom_se_fn_supplied
    ))
  }

  if (any(is_pm)) {
    stop(
      "Cannot mix 'pred_measure' objects with plain 'loo' objects. ",
      "Compare models using the same *_pred_measure() function for each model.",
      call. = FALSE
    )
  }

  # For plain `loo` objects only the model-name form of `rank_by` applies:
  # there is a single measure (elpd), so there is nothing to rank by.
  ref_model <- NULL
  if (!is.null(rank_by)) {
    if (is.character(rank_by) && length(rank_by) == 1L &&
        !is.na(rank_by) && rank_by %in% find_model_names(loos)) {
      ref_model <- rank_by
    } else {
      warning(
        "`rank_by` is only used for `pred_measure` comparisons, or to name the ",
        "reference model, and will be ignored.",
        call. = FALSE
      )
    }
  }
  if (custom_se_fn_supplied) {
    warning(
      "`custom_se_fn` is only used for `pred_measure` comparisons and will be ignored.",
      call. = FALSE
    )
  }

  # run pre-comparison checks
  model_compare_checks(loos)

  # compute elpd_diff and se_elpd_diff relative to best model
  ord <- model_compare_order(loos)
  comp <- model_compare_matrix(loos, ord = ord)
  rnms <- rownames(comp)
  ref_idx <- if (is.null(ref_model)) 1L else match(ref_model, rnms)
  diffs <- mapply(FUN = elpd_diffs, loos[ord[ref_idx]], loos[ord])
  colnames(diffs) <- rnms
  elpd_diff <- apply(diffs, 2, sum)
  se_diff <- apply(diffs, 2, se_elpd_diff)

  # compute probabilities that a model has worse elpd than the reference model
  # (the best model unless `rank_by` named one) using a normal approximation
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

  if (!is.null(ref_model)) {
    attr(comp, "compare_ref_model") <- ref_model
  }
  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' @rdname model_compare
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
#'   names (e.g. `c("elpd", "mse")`). Each table is sorted by its own measure,
#'   best model first, so the same model need not lead every table.
print.compare.loo <- function(x, ..., digits = 1, p_worse = TRUE, measures = NULL) {
  if (inherits(x, "old_compare.loo")) {
    return(unclass(x))
  }
  if (!inherits(x, "data.frame")) {
    class(x) <- c(class(x), "data.frame")
  }

  compare_measures <- attr(x, "compare_measures")
  if (!is.null(compare_measures)) {
    return(.print_compare_pred_measure(
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

  ref_model_attr <- attr(x, "compare_ref_model")
  if (!is.null(ref_model_attr)) {
    message("Differences computed against model ", ref_model_attr, ".")
  }
  .print_compare_diag_message(x, p_worse = p_worse)
  invisible(x)
}


# internal ----------------------------------------------------------------

#' Print `compare.loo` results from `pred_measure` comparisons
#' @noRd
.print_compare_pred_measure <- function(x, digits, p_worse, measures) {
  rank_by <- attr(x, "rank_by")
  ref_model_attr <- attr(x, "compare_ref_model")
  compare_measures <- attr(x, "compare_measures")
  compare_source <- attr(x, "compare_source")
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

  # The reference is always named, whichever way it was chosen: without
  # `rank_by` each measure keeps its own best model, which is the case most in
  # need of saying so and the only one that used to say nothing.
  # Printed rather than messaged: it labels the tables below, and `message()`
  # output is suppressed wholesale by knitr chunks and `suppressMessages()`.
  .cat_wrapped(
    .compare_reference_line(x, rank_by, ref_model_attr, compare_measures)
  )

  # LOO is the familiar default, so only name the source when it is not LOO.
  if (!is.null(compare_source) && !identical(compare_source, "loo")) {
    .cat_wrapped(
      "Predictive measures evaluated on ",
      .compare_source_label(compare_source),
      "."
    )
  }

  # Pareto k is a property of a model's PSIS-LOO approximation, not of any one
  # measure or of the comparison, so it is reported once for all models rather
  # than as a column inside a per-measure difference table.
  psis_shown <- .print_psis_diag_block(x)
  # A per-measure header already opens with a blank line; only the single-table
  # form needs one inserted here.
  if (psis_shown && is.null(measures)) {
    cat("\n")
  }

  for (measure in measures_to_print) {
    if (!is.null(measures)) {
      cat(
        "\n-- ", measure, " (vs ", .compare_ref_model(x, measure), ") --\n",
        sep = ""
      )
    }
    .print_compare_measure_table(
      x,
      measure = measure,
      digits = digits,
      p_worse = p_worse
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
      # The per-measure references are named in the header line, so this only
      # has to say which measures exist and how to see them.
      message(
        if (has_diag_msg) "\n",
        .wrap(
          "Other measures compared: ",
          paste(other, collapse = ", "),
          ". Use print(x, measures = \"all\")."
        )
      )
    }
  }

  invisible(x)
}

#' Reference model a measure's differences were computed against
#'
#' Without `rank_by` each measure has its own best model as reference, recorded
#' in attribute `compare_reference`. Falls back to the first row for objects
#' created before that attribute existed.
#' @noRd
.compare_ref_model <- function(x, measure) {
  refs <- attr(x, "compare_reference")
  if (!is.null(refs) && measure %in% names(refs)) {
    return(refs[[measure]])
  }
  x$model[[1L]]
}

#' Header line naming the reference each difference is computed against
#'
#' Always printed, so the reference is never left implicit. `rank_by` pins one
#' reference for every measure; without it each measure keeps its own best
#' model, which is the case most in need of being spelled out.
#' @noRd
#' @param x A `"compare.loo"` data frame.
#' @param rank_by Attribute `rank_by`, or `NULL`.
#' @param ref_model_attr Attribute `compare_ref_model`, or `NULL`.
#' @param compare_measures Bare names of all compared measures.
#' @return A single string.
.compare_reference_line <- function(x, rank_by, ref_model_attr,
                                    compare_measures) {
  if (!is.null(rank_by)) {
    return(paste0(
      "Models ranked by ", rank_by,
      " (reference: ", .compare_ref_model(x, rank_by), ")."
    ))
  }
  if (!is.null(ref_model_attr)) {
    return(paste0("All measures compared against model ", ref_model_attr, "."))
  }

  refs <- vapply(compare_measures, .compare_ref_model, character(1), x = x)
  if (length(compare_measures) == 1L) {
    return(paste0(
      "Models ranked by ", compare_measures, " (reference: ", refs, ")."
    ))
  }
  # Naming every reference stops being readable once there are many measures.
  if (length(compare_measures) > 4L) {
    return("Each measure compared against its own best model.")
  }
  paste0(
    "Each measure compared against its own best model (",
    paste0(compare_measures, ": ", refs, collapse = ", "),
    ")."
  )
}

#' Split `diag_elpd` entries into their count and threshold
#' @noRd
#' @param flags Character vector of `diag_elpd` values, such as
#'   `"25 k_psis > 0.62"`.
#' @return Data frame with numeric `bad_k` and `threshold`, both `NA` for an
#'   entry that does not match, so an unrecognised value is passed through
#'   rather than silently dropped.
.parse_diag_psis <- function(flags) {
  parts <- regmatches(flags, regexec("^([0-9]+) k_psis > ([0-9.]+)$", flags))
  field <- function(i) {
    vapply(
      parts,
      function(p) if (length(p) == 3L) as.numeric(p[[i]]) else NA_real_,
      numeric(1)
    )
  }
  data.frame(bad_k = field(2L), threshold = field(3L))
}

#' Wrap printed prose to the conventional 80-column terminal width
#'
#' Tables are wrapped by `print.data.frame()` at `getOption("width")`; this does
#' the same for the sentences around them, capped at 80 so the output stays
#' within a standard terminal however wide the option is set.
#' @noRd
#' @param ... Pieces of a single line, pasted together.
#' @return A single string with embedded newlines.
.wrap <- function(...) {
  width <- min(getOption("width", 80L), 80L)
  paste(strwrap(paste0(...), width = width), collapse = "\n")
}

#' Print prose wrapped to 80 columns
#' @noRd
#' @param ... Pieces of a single line, pasted together.
.cat_wrapped <- function(...) {
  cat(.wrap(...), "\n", sep = "")
}

#' Describe how many of the compared models a flag applies to
#' @noRd
#' @param n Number of flagged models.
#' @param total Number of compared models.
#' @return A string such as `"2 of 3 models"`, `"all 3 models"`, or
#'   `"both models"`.
.n_of_models <- function(n, total) {
  if (n < total) {
    return(paste0(n, " of ", total, " models"))
  }
  if (total == 2L) "both models" else paste0("all ", total, " models")
}

#' Print the PSIS-LOO diagnostics block
#'
#' Pareto \eqn{\hat{k}} describes a model's PSIS-LOO approximation, not any one
#' measure and not the comparison, so it is reported once per model above the
#' per-measure difference tables. Nothing is printed when no model is flagged,
#' or for sources other than LOO, which carry no `diag_elpd` column.
#' @noRd
#' @param x A `"compare.loo"` data frame.
#' @return `TRUE` invisibly when a block was printed.
.print_psis_diag_block <- function(x) {
  col <- x[["diag_elpd"]]
  if (is.null(col)) {
    return(invisible(FALSE))
  }
  flagged <- !is.na(col) & nzchar(col, keepNA = FALSE)
  if (!any(flagged)) {
    return(invisible(FALSE))
  }

  n_models <- length(col)
  models <- x$model[flagged]
  parsed <- .parse_diag_psis(col[flagged])

  # An unparsed entry has no count to sort or tabulate by, so fall back to the
  # stored strings rather than inventing numbers for them.
  if (anyNA(parsed$bad_k)) {
    .cat_wrapped(
      "PSIS-LOO unreliable for ", .n_of_models(length(models), n_models),
      "; measures may be biased."
    )
    print(
      data.frame(
        model = models,
        diag_elpd = col[flagged],
        stringsAsFactors = FALSE
      ),
      quote = FALSE,
      row.names = FALSE
    )
    return(invisible(TRUE))
  }

  ord <- order(parsed$bad_k, decreasing = TRUE)
  models <- models[ord]
  parsed <- parsed[ord, , drop = FALSE]
  # Thresholds depend on the number of draws, so models need not share one.
  common <- length(unique(parsed$threshold)) == 1L

  if (length(models) == 1L) {
    .cat_wrapped(
      "PSIS-LOO unreliable for ", models, " (", parsed$bad_k,
      " obs, k_psis > ", parsed$threshold, "); measures may be biased."
    )
    return(invisible(TRUE))
  }

  .cat_wrapped(
    "PSIS-LOO unreliable for ", .n_of_models(length(models), n_models),
    if (common) paste0(" (k_psis > ", parsed$threshold[[1L]], ")") else "",
    "; measures may be biased."
  )
  block <- data.frame(
    model = models,
    bad_k = parsed$bad_k,
    stringsAsFactors = FALSE
  )
  if (!common) {
    block$k_psis_threshold <- parsed$threshold
  }
  print(block, quote = FALSE, row.names = FALSE)
  invisible(TRUE)
}

#' Print one measure's comparison table
#' @noRd
.print_compare_measure_table <- function(x, measure, digits, p_worse) {
  if (.is_elpd_measure(measure)) {
    diff_col <- "elpd_diff"
    se_col <- "se_diff"
    diff_name <- "elpd_diff"
    se_name <- "se_diff"
  } else {
    diff_col <- paste0(measure, "_diff")
    se_col <- paste0(measure, "_se_diff")
    diff_name <- diff_col
    # Print the column name the object actually carries, so the header matches
    # `comp$mae_se_diff` and names the measure an NA belongs to.
    se_name <- se_col
  }

  if (!all(c(diff_col, se_col) %in% colnames(x))) {
    stop(
      "Comparison columns for measure '", measure, "' are missing.",
      call. = FALSE
    )
  }

  # The data frame carries one row order for all measures (by `rank_by`), but a
  # measure's own best model need not be first in it. Sort each printed table by
  # its own difference so the best model is always the first row and the
  # differences run in decreasing order.
  ord <- order(x[[diff_col]], decreasing = TRUE, na.last = TRUE)

  x2 <- data.frame(
    model = x$model[ord],
    diff = unname(.fr(x[[diff_col]][ord], digits)),
    se_diff = unname(.fr(x[[se_col]][ord], digits)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(x2)[2:3] <- c(diff_name, se_name)

  if (.is_elpd_measure(measure) && p_worse && "p_worse" %in% colnames(x)) {
    x2$p_worse <- unname(.fr(x[["p_worse"]][ord], digits = 2))
    x2$diag_diff <- x[["diag_diff"]][ord]
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

#' Is an object a predictive measure result?
#' @noRd
is.pred_measure <- function(x) {
  inherits(x, "pred_measure")
}

#' Is an object a PSIS-LOO predictive measure result?
#' @noRd
is.loo_pred_measure <- function(x) {
  inherits(x, "loo_pred_measure")
}

#' Resolve the shared evaluation source of `pred_measure` inputs
#'
#' All models in one comparison must be evaluated the same way: paired
#' differences between, say, a LOO and a k-fold result would contrast different
#' held-out schemes rather than the models themselves. Measure names are also
#' suffixed by source, so mixed inputs share no comparable `pointwise` columns.
#' @noRd
#' @param loos List of `"pred_measure"` objects.
#' @return The shared `source` string: `"loo"`, `"kfold"`, `"test"`, or
#'   `"insample"`.
.compare_source <- function(loos) {
  sources <- vapply(loos, function(x) {
    source <- attr(x, "source")
    if (is.null(source)) NA_character_ else source
  }, character(1))

  if (anyNA(sources)) {
    stop(
      "All inputs must be results of insample_pred_measure(), ",
      "loo_pred_measure(), kfold_pred_measure(), or test_pred_measure().",
      call. = FALSE
    )
  }
  if (any(sources != sources[1L])) {
    labels <- unique(vapply(loos, .pred_measure_source_label, character(1)))
    stop(
      paste0(
        "All models must be evaluated on the same source, but got: ",
        paste(labels, collapse = ", "),
        ". Recompute all models with the same *_pred_measure() function."
      ),
      call. = FALSE
    )
  }
  # `loos` may be a named list, which would make `vapply()` return a named
  # vector and break the `identical()` checks against a bare string.
  unname(sources[1L])
}

#' Warn that in-sample comparisons are optimistically biased
#' @noRd
.warn_insample_compare <- function(source) {
  if (!identical(source, "insample")) {
    return(invisible(NULL))
  }
  warning(
    "Comparing in-sample predictive measures. In-sample scores are ",
    "optimistically biased and favor more complex models. For out-of-sample ",
    "comparison use loo_pred_measure(), kfold_pred_measure(), or ",
    "test_pred_measure().",
    call. = FALSE
  )
  invisible(NULL)
}

#' Warn when k-fold results do not share the same number of folds
#' @noRd
#' @param loos List of `"kfold"` or `"kfold_pred_measure"` objects.
.warn_kfold_K_mismatch <- function(loos) {
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

#' Compare `loo_pred_measure` objects (multi-measure path)
#' @noRd
#' @param loos List of `pred_measure` objects, all sharing one evaluation
#'   source.
#' @param rank_by Bare measure name used to order models.
#' @param custom_se_fn How to compute the standard error of the difference for
#'   custom measures, see `.resolve_custom_se_fns()`.
#' @param custom_se_fn_supplied Whether `custom_se_fn` was given at all, as
#'   opposed to being an explicit `NULL`.
compare_pred_measure <- function(loos, rank_by = NULL, custom_se_fn = NULL,
                                 custom_se_fn_supplied = FALSE) {
  # Resolve the source before the generic checks: mixed sources usually also
  # differ in their number of observations, and "you mixed LOO with k-fold" is
  # far more actionable than "your models have different N".
  source <- .compare_source(loos)
  model_compare_checks(
    loos,
    class_check = is.pred_measure,
    class_msg = "All inputs must have class 'pred_measure'.",
    kfold_checks = FALSE
  )
  .warn_insample_compare(source)
  if (identical(source, "kfold")) {
    .warn_kfold_K_mismatch(loos)
  }
  .compare_metadata_check(loos)
  .warn_omitted_compare_measures(loos)

  rank_spec <- .resolve_rank_by(loos, rank_by)
  rank_measure <- rank_spec$measure
  compare_cols <- .compare_pointwise_cols(loos)
  custom_se_fns <- .resolve_custom_se_fns(
    loos,
    compare_cols,
    custom_se_fn,
    custom_se_fn_supplied
  )
  .inform_compare_sign_conversion(compare_cols, loos)
  ord <- model_compare_order(loos, rank_measure$internal)
  loos_ord <- loos[ord]
  # With an explicit `rank_by` a single model is the reference for every
  # measure: the top-ranked one when `rank_by` names a measure, the named one
  # when it names a model. Without `rank_by`, each measure gets its own best
  # model as reference, so e.g. `mse_diff` may be relative to a different model
  # than `elpd_diff`.
  per_measure_ref <- identical(rank_spec$kind, "default")

  comp <- model_compare_matrix(
    loos_ord,
    bare_names = TRUE,
    ord = seq_along(loos_ord)
  )
  rnms <- rownames(comp)
  n_obs <- nrow(loos_ord[[1L]]$pointwise)
  pinned_ref_idx <- if (identical(rank_spec$kind, "model")) {
    match(rank_spec$model, rnms)
  } else {
    1L
  }

  diff_cols <- list()
  ref_models <- character(0)
  for (col in compare_cols) {
    bare <- .display_name(col, loos_ord)
    ref_idx <- if (per_measure_ref) {
      model_compare_order(loos_ord, col)[[1L]]
    } else {
      pinned_ref_idx
    }
    ref_loo <- loos_ord[[ref_idx]]
    ref_models[[bare]] <- rnms[[ref_idx]]
    method <- .measure_pointwise_diff_method(loos_ord, col)
    se_fn <- if (identical(method, "custom")) custom_se_fns[[bare]] else NULL
    if (is.character(se_fn)) {
      .check_declared_aggregation(loos_ord, col, se_fn)
    }
    pair_stats <- vapply(
      loos_ord,
      .pair_measure_stats,
      FUN.VALUE = c(diff = 0, se = 0),
      ref = ref_loo,
      col = col,
      method = method,
      loos = loos_ord,
      se_fn = se_fn
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

  # `diag_elpd` reports PSIS Pareto k, which only exists for the LOO source;
  # for the others it would be an all-blank column.
  model_cols <- data.frame(
    model = rnms,
    diff_cols,
    stringsAsFactors = FALSE
  )
  if (identical(source, "loo")) {
    model_cols$diag_elpd <- diag_elpd(loos_ord)
  }

  comp <- cbind(model_cols, as.data.frame(comp))
  rownames(comp) <- NULL

  model_order_stat_check(
    loos_ord,
    seq_along(loos_ord),
    rank_col = rank_measure$internal
  )

  if (identical(rank_spec$kind, "measure")) {
    attr(comp, "rank_by") <- rank_measure$bare
  } else if (identical(rank_spec$kind, "model")) {
    attr(comp, "compare_ref_model") <- rank_spec$model
  }
  attr(comp, "compare_reference") <- ref_models
  attr(comp, "compare_source") <- source
  attr(comp, "compare_measures") <- .compare_measures(loos)
  attr(comp, "sign_converted_measures") <- .compare_sign_converted_measures(
    compare_cols,
    loos
  )
  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' Measure-name suffix used by a comparison's evaluation source
#'
#' `.measure_result_name()` suffixes measure names by source (`elpd_loo`,
#' `elpd_kfold`, `elpd_test`, and a bare `elpd` for in-sample). This is the
#' inverse, so display names can be recovered for any source.
#' @noRd
#' @param loos List of `"pred_measure"` objects, or `NULL` when the source is
#'   unknown (falls back to the LOO suffix).
.compare_suffix <- function(loos = NULL) {
  if (is.null(loos)) {
    return("_loo")
  }
  source <- attr(loos[[1L]], "source")
  if (is.null(source) || identical(source, "insample")) "" else paste0("_", source)
}

#' Strip the source suffix for `model_compare` display names
#' @noRd
#' @param col Measure column name.
#' @param loos List of model results the column came from; determines which
#'   suffix to strip. Deriving the suffix from the source (rather than matching
#'   any of `_loo|_kfold|_test`) keeps a custom measure named e.g. `my_test`
#'   intact outside a test-set comparison.
.display_name <- function(col, loos = NULL) {
  suffix <- .compare_suffix(loos)
  if (!nzchar(suffix)) {
    return(col)
  }
  sub(paste0(suffix, "$"), "", col)
}

#' Map bare measure name to `pointwise` column name
#' @noRd
.pointwise_col <- function(name, cols, loos = NULL) {
  if (name %in% cols) {
    return(name)
  }
  internal <- paste0(name, .compare_suffix(loos))
  if (internal %in% cols) {
    return(internal)
  }
  stop(
    paste0(
      "Measure '", name, "' not found in all models. ",
      "Available measures: ",
      paste(vapply(cols, .display_name, character(1), loos = loos), collapse = ", ")
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

#' Check that `measure_info` is consistent across models
#' @noRd
.compare_metadata_check <- function(loos) {
  bare_measures <- .compare_measures(loos)
  if (!length(bare_measures)) {
    return(invisible(NULL))
  }

  for (bare in bare_measures) {
    infos <- lapply(loos, function(x) {
      measure_info <- attr(x, "measure_info")
      if (is.null(measure_info)) {
        return(NULL)
      }
      measure_info[[bare]]
    })
    has_info <- !vapply(infos, is.null, logical(1))
    if (any(has_info) && !all(has_info)) {
      stop(
        "Not all models provide `measure_info` for measure '",
        bare,
        "'. Recompute all inputs with the current version of `loo_pred_measure()`.",
        call. = FALSE
      )
    }
    # `extra` holds per-measure auxiliary data (for `r2`, the pointwise
    # baseline derived from `y`), which legitimately differs when models are
    # fitted to different data. That case is already reported by the `yhash`
    # warning, so comparing `extra` here would only mislabel it as a
    # disagreement about the measure itself.
    non_null <- lapply(infos[has_info], function(info) {
      info$extra <- NULL
      info
    })
    if (length(non_null) > 1L) {
      ref <- non_null[[1L]]
      inconsistent <- vapply(
        non_null[-1L],
        function(info) !identical(info, ref),
        logical(1)
      )
      if (any(inconsistent)) {
        stop(
          "Models disagree on `measure_info` for measure '",
          bare,
          "'. For a custom measure, ensure all models use the same ",
          "`measure_loss` declaration.",
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
      unname(vapply(cols, .display_name, character(1), loos = loos))
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

#' Bare measure names available for comparison across models
#' @noRd
.compare_measures <- function(loos) {
  cols <- .compare_pointwise_cols(loos)
  unname(vapply(cols, .display_name, character(1), loos = loos))
}

#' Resolve `rank_by` to bare and internal `pointwise` column names
#' @noRd
.resolve_rank_measure <- function(loos, rank_by = NULL) {
  cols <- .compare_pointwise_cols(loos)
  bare <- if (is.null(rank_by)) "elpd" else rank_by
  internal <- .pointwise_col(bare, cols, loos)
  list(
    bare = .display_name(internal, loos),
    internal = internal
  )
}

#' Match `rank_by` against the measures shared by all models
#'
#' Like `.pointwise_col()` but returns `NULL` instead of erroring, so callers
#' can fall back to interpreting `rank_by` as a model name.
#' @noRd
.match_rank_measure <- function(loos, rank_by, cols) {
  if (rank_by %in% cols) {
    return(rank_by)
  }
  internal <- paste0(rank_by, .compare_suffix(loos))
  if (internal %in% cols) {
    return(internal)
  }
  NULL
}

#' Resolve `rank_by` to either a measure or a reference model
#'
#' `rank_by` accepts a bare measure name (rank models by that measure and use
#' the top-ranked model as reference) or a model name (keep the default `elpd`
#' ordering but pin that model as the reference for every measure).
#' @noRd
#' @return A list with `kind` (`"default"`, `"measure"`, or `"model"`),
#'   `measure` (the resolved ranking measure, as `.resolve_rank_measure()`
#'   returns it) and `model` (the pinned reference model name, or `NULL`).
.resolve_rank_by <- function(loos, rank_by = NULL) {
  if (is.null(rank_by)) {
    return(list(
      kind = "default",
      measure = .resolve_rank_measure(loos),
      model = NULL
    ))
  }

  if (!is.character(rank_by) || length(rank_by) != 1L || is.na(rank_by)) {
    stop(
      "`rank_by` must be a single measure name or model name.",
      call. = FALSE
    )
  }

  cols <- .compare_pointwise_cols(loos)
  internal <- .match_rank_measure(loos, rank_by, cols)
  model_names <- find_model_names(loos)

  if (!is.null(internal) && rank_by %in% model_names) {
    warning(
      "`rank_by = \"", rank_by, "\"` matches both a measure and a model name; ",
      "ranking by the measure. Rename the model to rank by the model instead.",
      call. = FALSE
    )
  }

  if (!is.null(internal)) {
    return(list(
      kind = "measure",
      measure = list(
        bare = .display_name(internal, loos),
        internal = internal
      ),
      model = NULL
    ))
  }

  if (rank_by %in% model_names) {
    return(list(
      kind = "model",
      measure = .resolve_rank_measure(loos),
      model = rank_by
    ))
  }

  stop(
    paste0(
      "`rank_by` value '", rank_by, "' is neither a measure nor a model name. ",
      "Available measures: ",
      paste(vapply(cols, .display_name, character(1), loos = loos), collapse = ", "),
      ". Available models: ",
      paste(model_names, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

#' Is a measure an ELPD-family measure (for `p_worse` / `diag_diff`)?
#'
#' Matches on the raw column name: every source suffix (`elpd_loo`,
#' `elpd_kfold`, `elpd_test`, bare `elpd`) shares the `elpd` prefix, so no
#' suffix stripping is needed here.
#' @noRd
.is_elpd_measure <- function(name) {
  grepl("^elpd", name)
}

#' Look up the per-measure information recorded on a result object
#' @noRd
.get_measure_info <- function(loos, bare) {
  measure_info <- attr(loos[[1L]], "measure_info")
  if (is.null(measure_info)) {
    return(NULL)
  }
  measure_info[[bare]]
}

#' Names of the built-in measures that are losses (lower is better)
#' @noRd
.builtin_loss_measures <- function() {
  names(Filter(function(spec) isTRUE(spec$loss), .measure_spec))
}

#' Whether a measure is a loss (lower is better)
#'
#' Measure values are always stored on the measure's own scale, so this equally
#' describes the measure and the values recorded for it.
#' @noRd
.measure_is_loss <- function(name, loos = NULL) {
  bare <- .display_name(name, loos)

  if (!is.null(loos)) {
    info <- .get_measure_info(loos, bare)
    if (!is.null(info) && !is.null(info$loss)) {
      return(isTRUE(info$loss))
    }
  }

  spec <- .measure_spec[[bare]]
  if (!is.null(spec)) {
    return(isTRUE(spec$loss))
  }
  bare %in% .builtin_loss_measures()
}

#' Bare names of measures whose sign is flipped for `model_compare()`
#' @noRd
.compare_sign_converted_measures <- function(cols, loos) {
  bare <- vapply(cols, .display_name, character(1), loos = loos)
  unique(bare[vapply(
    cols,
    function(col) .measure_is_loss(col, loos),
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
    "\nreported on a utility scale (higher is better)."
  )
  invisible(NULL)
}

#' How to aggregate paired pointwise differences for a measure
#'
#' Taken from the measure's stored `measure_info`: `"sum"` or `"mean"` when
#' the overall estimate is the sum or the mean of its pointwise contributions,
#' `"measure_specific"` when the built-in measure supplies its own
#' `se_diff_fun`, and `"custom"` for custom measures, whose standard error is
#' supplied at comparison time through `model_compare(custom_se_fn = )`.
#' Nothing is inferred. `"estimates_only"` is only reached by legacy objects
#' carrying no `measure_info` at all.
#' @noRd
.measure_pointwise_diff_method <- function(loos, col) {
  bare <- .display_name(col, loos)
  info <- .get_measure_info(loos, bare)
  if (!is.null(info) && !is.null(info$diff_method)) {
    return(info$diff_method)
  }

  if (.is_elpd_measure(col) || bare == "ic") {
    return("sum")
  }

  "estimates_only"
}

#' Check that a declared `"sum"`/`"mean"` aggregation matches the estimate
#'
#' Only called when the user declares `custom_se_fn = "sum"` or `"mean"` for a
#' custom measure. This is the computation `.measure_pointwise_diff_method()`
#' used to run as autodetection, inverted: rather than guessing the aggregation,
#' it verifies the one the user asserted.
#' @noRd
.check_declared_aggregation <- function(loos, col, method) {
  ref <- loos[[1L]]
  est <- ref$estimates[col, "Estimate"]
  pw <- ref$pointwise[, col, drop = TRUE]
  agg <- if (identical(method, "sum")) sum(pw) else mean(pw)

  ok <- length(pw) > 0L && is.finite(est) && is.finite(agg)
  if (ok) {
    tol <- sqrt(.Machine$double.eps) * max(abs(c(est, pw)), na.rm = TRUE)
    ok <- isTRUE(all.equal(est, agg, tolerance = tol, check.attributes = FALSE))
  }
  if (!ok) {
    warning(
      "`custom_se_fn = \"", method, "\"` was declared for measure '",
      .display_name(col, loos), "', but ", method,
      "(pointwise) does not reproduce its estimate.\n",
      "The reported standard error may be wrong.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Resolve a built-in measure's `se_diff_fun`
#'
#' Built-in measures with `diff_method = "measure_specific"` name an entry of
#' `.se_diff_funs`. Custom measures never reach this; their standard error comes
#' from `model_compare(custom_se_fn = )`, see `.resolve_custom_se_fns()`.
#' @noRd
.measure_se_diff_fun <- function(loos, col) {
  bare <- .display_name(col, loos)
  info <- .get_measure_info(loos, bare)

  fun <- info$se_diff_fun
  if (is.null(fun)) {
    fun <- .measure_spec[[bare]]$se_diff_fun
  }
  if (is.character(fun)) {
    fun <- .se_diff_funs[[fun]]
  }
  if (!is.function(fun)) {
    stop(
      paste0(
        "No 'se_diff_fun' available for measure '", bare, "'."
      ),
      call. = FALSE
    )
  }
  fun
}

#' Accepted string shorthands for `custom_se_fn`
#' @noRd
.custom_se_fn_keywords <- c("sum", "mean")

#' Validate one `custom_se_fn` value
#' @noRd
#' @return The value itself, or `NULL`.
.check_custom_se_fn_value <- function(value, bare) {
  if (is.null(value) || is.function(value)) {
    return(value)
  }
  if (is.character(value) && length(value) == 1L &&
      value %in% .custom_se_fn_keywords) {
    return(value)
  }
  stop(
    "Invalid `custom_se_fn` for measure '", bare,
    "'. It must be a function, ",
    paste0("\"", .custom_se_fn_keywords, "\"", collapse = " or "),
    ", or NULL.",
    call. = FALSE
  )
}

#' Message listing what `custom_se_fn` accepts
#' @noRd
.custom_se_fn_help <- function() {
  paste0(
    "Pass a function computing the standard error of the difference, ",
    paste0("\"", .custom_se_fn_keywords, "\"", collapse = " or "),
    "\nto use the paired pointwise formula, or NULL to report the difference ",
    "with an NA standard error."
  )
}

#' Resolve `custom_se_fn` to a per-measure lookup
#'
#' Custom measures carry `diff_method = "custom"` and declare nothing about
#' their standard error, so the person running the comparison supplies it. A
#' bare value is only unambiguous when exactly one custom measure is compared;
#' otherwise a list keyed by bare measure name is required.
#' @noRd
#' @param custom_se_fn The user's `custom_se_fn` argument, already normalised to
#'   `NULL` when it was not supplied.
#' @param supplied Whether the argument was given at all, as opposed to being
#'   an explicit `NULL`.
#' @return Named list keyed by bare measure name; each element is a function,
#'   `"sum"`, `"mean"`, or `NULL`.
.resolve_custom_se_fns <- function(loos, compare_cols, custom_se_fn, supplied) {
  is_custom <- vapply(
    compare_cols,
    function(col) {
      identical(.measure_pointwise_diff_method(loos, col), "custom")
    },
    logical(1)
  )
  custom_bare <- unname(vapply(
    compare_cols[is_custom],
    .display_name,
    character(1),
    loos = loos
  ))

  if (!length(custom_bare)) {
    if (supplied) {
      warning(
        "`custom_se_fn` is only used for custom measures and will be ignored.",
        call. = FALSE
      )
    }
    return(list())
  }

  if (!supplied) {
    stop(
      if (length(custom_bare) == 1L) "Measure '" else "Measures '",
      paste(custom_bare, collapse = "', '"),
      if (length(custom_bare) == 1L) {
        "' is a custom measure, so `custom_se_fn` must be supplied.\n"
      } else {
        "' are custom measures, so `custom_se_fn` must be supplied.\n"
      },
      .custom_se_fn_help(),
      call. = FALSE
    )
  }

  # A bare function or keyword: unambiguous only for a single custom measure.
  if (is.null(custom_se_fn) || is.function(custom_se_fn) ||
      is.character(custom_se_fn)) {
    if (!is.null(custom_se_fn) && length(custom_bare) > 1L) {
      stop(
        "`custom_se_fn` must be a named list when more than one custom measure ",
        "is compared.\nName an entry for each of: ",
        paste(custom_bare, collapse = ", "), ".",
        call. = FALSE
      )
    }
    value <- .check_custom_se_fn_value(custom_se_fn, custom_bare[[1L]])
    return(stats::setNames(rep(list(value), length(custom_bare)), custom_bare))
  }

  if (!is.list(custom_se_fn)) {
    stop(
      "`custom_se_fn` must be a function, ",
      paste0("\"", .custom_se_fn_keywords, "\"", collapse = " or "),
      ", NULL, or a named list of those.",
      call. = FALSE
    )
  }

  nms <- names(custom_se_fn)
  if (length(custom_se_fn) && (is.null(nms) || any(!nzchar(nms)))) {
    stop(
      "Every element of `custom_se_fn` must be named after a custom measure. ",
      "Expected name(s): ",
      paste(custom_bare, collapse = ", "), ".",
      call. = FALSE
    )
  }
  unknown <- setdiff(nms, custom_bare)
  if (length(unknown)) {
    stop(
      "Unknown measure(s) in `custom_se_fn`: ",
      paste(unknown, collapse = ", "),
      ". Custom measure(s) compared: ",
      paste(custom_bare, collapse = ", "), ".",
      call. = FALSE
    )
  }
  missing_measures <- setdiff(custom_bare, nms)
  if (length(missing_measures)) {
    stop(
      "`custom_se_fn` has no entry for custom measure(s): ",
      paste(missing_measures, collapse = ", "), ".\n",
      .custom_se_fn_help(),
      call. = FALSE
    )
  }

  stats::setNames(
    lapply(custom_bare, function(bare) {
      .check_custom_se_fn_value(custom_se_fn[[bare]], bare)
    }),
    custom_bare
  )
}

#' Assemble one model's inputs for an `se_diff_fun`
#'
#' Every element describes the single model `x`, on the measure's natural scale,
#' including `extra`, which is read from that model's own `measure_info` rather
#' than the reference model's.
#' @noRd
.se_diff_input <- function(x, col) {
  list(
    estimate = x$estimates[col, "Estimate"],
    se = x$estimates[col, "SE"],
    pointwise = x$pointwise[, col, drop = TRUE],
    extra = .get_measure_info(list(x), .display_name(col, list(x)))$extra
  )
}

#' Validate the value returned by an `se_diff_fun` or `custom_se_fn`
#' @noRd
#' @param what Name of the argument or attribute the function came from, used
#'   only to make the error point at what the user can change.
.validate_se_diff <- function(se, col, loos = NULL, what = "se_diff_fun") {
  if (!is.numeric(se) || length(se) != 1L) {
    stop(
      paste0(
        "The `", what, "` for measure '", .display_name(col, loos),
        "' must return a numeric scalar."
      ),
      call. = FALSE
    )
  }
  unname(se)
}

#' Paired measure difference and SE for one model vs a reference
#' @noRd
#' @param se_fn For `method = "custom"` only: the value resolved from
#'   `model_compare(custom_se_fn = )` for this measure. A function, the string
#'   `"sum"` or `"mean"`, or `NULL` for an `NA` standard error.
.pair_measure_stats <- function(cmp, ref, col, method = NULL, loos = list(ref),
                                se_fn = NULL) {
  if (is.null(method)) {
    method <- .measure_pointwise_diff_method(c(list(ref, cmp)), col)
  }

  flip <- .measure_is_loss(col, loos)
  est_utility <- function(estimates) {
    val <- estimates[col, "Estimate"]
    if (flip) -val else val
  }

  if (method == "custom") {
    if (is.character(se_fn)) {
      # "sum"/"mean" reuse the paired pointwise branch below, so a custom
      # measure declaring "mean" behaves exactly like the built-in `mae`.
      method <- se_fn
    } else {
      diff <- est_utility(cmp$estimates) - est_utility(ref$estimates)
      if (is.null(se_fn)) {
        return(c(diff = diff, se = NA_real_))
      }
      se <- se_fn(
        ref = .se_diff_input(ref, col),
        cmp = .se_diff_input(cmp, col)
      )
      return(c(
        diff = diff,
        se = .validate_se_diff(se, col, loos, what = "custom_se_fn")
      ))
    }
  }

  if (method == "estimates_only") {
    return(c(
      diff = est_utility(cmp$estimates) - est_utility(ref$estimates),
      se = NA_real_
    ))
  }

  if (method == "measure_specific") {
    se_diff_fun <- .measure_se_diff_fun(loos, col)
    se <- se_diff_fun(
      ref = .se_diff_input(ref, col),
      cmp = .se_diff_input(cmp, col)
    )
    return(c(
      diff = est_utility(cmp$estimates) - est_utility(ref$estimates),
      se = .validate_se_diff(se, col, loos)
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
model_compare_checks <- function(
  loos,
  class_check = is.loo,
  class_msg = "All inputs should have class 'loo'.",
  kfold_checks = TRUE
) {
  ## errors
  if (length(loos) <= 1L) {
    stop("At least two models are required for comparison.", call. = FALSE)
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


#' Build estimates table for `model_compare()` ordering and matrix output
#' @noRd
.model_compare_estimates_table <- function(loos, bare_names = FALSE) {
  sapply(loos, function(x) {
    est <- x$estimates
    rows <- if (bare_names) .display_name(rownames(est), loos) else rownames(est)
    setNames(c(est), nm = c(rows, paste0("se_", rows)))
  })
}

#' Compute the model_compare matrix
#' @noRd
#' @param loos List of `"loo"` objects.
#' @param bare_names If `TRUE`, strip `_loo` suffixes from estimate row names.
#' @param ord Optional model ordering indices; computed from ELPD when `NULL`.
model_compare_matrix <- function(loos, bare_names = FALSE, ord = NULL) {
  tmp <- .model_compare_estimates_table(loos, bare_names = bare_names)
  colnames(tmp) <- find_model_names(loos)
  comp <- t(tmp)

  if (is.null(ord)) {
    ord <- model_compare_order(loos)
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
