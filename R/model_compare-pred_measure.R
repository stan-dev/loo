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
throw_insample_compare_warning <- function(source) {
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
  throw_insample_compare_warning(source)
  if (identical(source, "kfold")) {
    throw_kfold_K_mismatch_warning(loos)
  }
  .compare_metadata_check(loos)
  throw_omitted_compare_measures_warning(loos)

  rank_spec <- .resolve_rank_by(loos, rank_by)
  rank_measure <- rank_spec$measure
  compare_cols <- .compare_pointwise_cols(loos)
  custom_se_fns <- .resolve_custom_se_fns(
    loos,
    compare_cols,
    custom_se_fn,
    custom_se_fn_supplied
  )
  inform_compare_sign_conversion(compare_cols, loos)
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

  # `rank_by` records how the reference was chosen, as a tagged record: `kind`
  # is the branch taken, `measure` is the measure rows are ordered by (always
  # set), and `model` is the pinned reference model, or `NULL`. `kind` is kept
  # because a name can match both a measure and a model.
  attr(comp, "rank_by") <- list(
    kind = rank_spec$kind,
    measure = rank_measure$bare,
    model = rank_spec$model
  )
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
    if (!all(has_info)) {
      stop(
        if (any(has_info)) "Not all models provide" else "No model provides",
        " `measure_info` for measure '",
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
throw_omitted_compare_measures_warning <- function(loos) {
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
inform_compare_sign_conversion <- function(cols, loos) {
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
#' Nothing is inferred. Every compared measure carries a `diff_method`;
#' `.compare_metadata_check()` has already rejected the inputs otherwise, so the
#' `elpd`/`ic` branch below only covers direct internal calls.
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

  stop(
    "No `measure_info` for measure '",
    bare,
    "'. Recompute all inputs with the current version of `loo_pred_measure()`.",
    call. = FALSE
  )
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
