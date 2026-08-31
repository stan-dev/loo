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
#' @param rank_by A single string naming either a **measure** or a **model**,
#'   used to pin one reference model for all pairwise differences.
#'
#'   A **measure name** ([`pred_measure`][pred_measure] comparisons only) orders
#'   models by that measure and makes the top-ranked model the reference. Bare
#'   names are matched regardless of source, so `rank_by = "rmse"` selects
#'   `rmse_loo`, `rmse_kfold`, or `rmse_test` as appropriate.
#'
#'   A **model name** (one of the names in the `model` column, i.e. the list
#'   names or `model1`, `model2`, ...) pins that model as the reference,
#'   whichever model performs best, and leaves rows ordered by `"elpd"`. This
#'   form also works for classic comparisons, where `elpd_diff` is then relative
#'   to the named model rather than to the best one. A name matching both a
#'   measure and a model is treated as the measure, with a warning.
#'
#'   With `rank_by = NULL` (the default) rows are ordered by `"elpd"` and each
#'   measure is compared against *its own* best model, so `mse_diff` may use a
#'   different reference than `elpd_diff`. Each `{measure}_diff` column then has
#'   exactly one `0` entry, at that measure's best model.
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
#'   For two or more custom measures, pass a list named by bare measure name,
#'   e.g. `list(huber = "mean", nrmse = my_se_fn)`. Ignored, with a warning,
#'   when no custom measure is present.
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
#'     \item{`rank_by`}{
#'       How the reference model was chosen, as a list with elements `kind`
#'       (`"default"`, `"measure"`, or `"model"`, for the three cases described
#'       under `rank_by` above), `measure` (bare name of the measure the rows
#'       are ordered by, always set, `"elpd"` by default) and `model` (the
#'       pinned reference model, or `NULL` unless `kind` is `"model"`).
#'     }
#'     \item{`compare_reference`}{
#'       A named character vector giving, for each measure, the model its
#'       differences were computed against. All entries name the same model
#'       unless `kind` is `"default"`.
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
#'   `rank_by` and `compare_reference` are set for every comparison; the last
#'   three are set for [`pred_measure`][pred_measure] comparisons only.
#'
#' @details
#' ## Differences and their standard errors
#'   Differences are pairwise: every model is compared with one reference model,
#'   whose own `{measure}_diff` is therefore `0`. See `rank_by` for how that
#'   reference is chosen. When it is the best model on a measure, as in classic
#'   comparisons, the remaining differences for that measure are all negative.
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
#'   `sign_converted_measures` attribute and in a message, for example:
#'   "For model comparison, differences for mse are reported on a utility scale
#'   (higher is better)."
#'
#'   A custom measure is treated as a utility unless it declares otherwise with
#'   `attr(my_fun, "measure_loss") <- TRUE`. The declaration also determines the
#'   direction of `rank_by`, so an undeclared loss is both flipped and ranked in
#'   the wrong direction; see [insample_pred_measure()].
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
#'   * `"custom"`: a custom measure declares nothing, so `custom_se_fn` must be
#'     supplied. `{measure}_se_diff` is `NA` only when `custom_se_fn` is an
#'     explicit `NULL` for that measure.
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
#'   differences using the median model (by ELPD, or by `rank_by` for
#'   `pred_measure` comparisons) as the baseline, and estimate whether the
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
#'   [Differences and their standard errors in model comparison](https://mc-stan.org/loo/articles/articles-online-only/model-comparison.html)
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
  if (missing(x)) {
    dots <- list(...)
    if (!length(dots)) {
      stop("No models supplied.", call. = FALSE)
    }
    # `custom_se_fn` has no default: omitted and explicit `NULL` differ, so it
    # is forwarded only when the caller supplied it.
    args <- list(dots, rank_by = rank_by)
    if (!missing(custom_se_fn)) {
      args$custom_se_fn <- custom_se_fn
    }
    return(do.call(model_compare, args))
  }
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

  # Same attribute contract as the `pred_measure` path, with the single
  # measure `"elpd"`: `rank_by` records how the reference was chosen and
  # `compare_reference` names the model it resolved to.
  attr(comp, "rank_by") <- list(
    kind = if (is.null(ref_model)) "default" else "model",
    measure = "elpd",
    model = ref_model
  )
  attr(comp, "compare_reference") <- c(elpd = rnms[[ref_idx]])
  class(comp) <- c("compare.loo", class(comp))
  comp
}

#' Reference model a measure's differences were computed against
#'
#' Without `rank_by` each measure has its own best model as reference, recorded
#' in attribute `compare_reference`. Falls back to the first row for objects
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
