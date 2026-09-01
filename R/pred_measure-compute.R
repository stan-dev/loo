#' Shared parameters for predictive measure wrappers
#'
#' @description
#' Parameter definitions shared by the user-facing entry points and the
#' internal engine [do_pred_measure()].
#'
#' @param y Vector of observed values (`n`). Required for distributional and
#'   point-prediction measures such as `crps`, `mae`, and `acc`.
#' @param ypred Matrix of posterior predictive draws (`S` draws × `n`
#'   observations), typically from [brms::posterior_predict()]. Required for
#'   distributional measures such as `crps`, `rps`, and `scrps`.
#' @param mupred Matrix of posterior expected values (`S` × `n`), typically from
#'   [brms::posterior_epred()]. Required for point-prediction measures such as
#'   `mae`, `rmse`, `r2`, and `acc`.
#' @param ylp Matrix of pointwise log predictive densities or probabilities
#'   (`S` × `n`), typically from [brms::log_lik()]. Required for density-based
#'   summaries (`elpd`, `mlpd`, `ic`).
#' @param ylp_test Matrix of pointwise log predictive densities for holdout
#'   observations (`S` × `n_test`), typically from
#'   `brms::log_lik(fit, newdata = test_data)`. Used with `ylp` (from the
#'   training fit) in [test_pred_measure()] to score genuinely new data.
#' @param predperf An existing predictive measure object (class
#'   `"pred_measure"`) to update. When supplied, base density summaries and
#'   (for LOO) PSIS weights are reused instead of recomputed.
#' @param measure Additional measures beyond the base summary `elpd` (always
#'   included when `ylp` is available). Can be:
#'   \itemize{
#'     \item A **character vector** of built-in names; see
#'       [supported_measures_list].
#'     \item A **function** with attribute `"measure_name"` for one custom measure.
#'     \item A **list** mixing character scalars (built-in names) and named
#'       functions (custom measures), e.g. `list("rps", my_metric = my_fun)`.
#'   }
#'   Custom functions are called with any of `y`, `ypred`, `mupred`, `ylp`, and
#'   `log_weights` that appear in their formals, plus arguments from `control`.
#'   They must return a list with  `estimates` and `pointwise`.
#'
#'   A custom measure declares whether it is a loss (lower is better) or a
#'   utility (higher is better) with attribute `"measure_loss"`:
#'   `attr(my_fun, "measure_loss") <- TRUE` for a loss. Without it a custom
#'   measure is taken to be a utility. [model_compare()] uses the declaration to
#'   put all measures on a common utility scale and to rank models, so an
#'   undeclared loss is compared and ranked in the wrong direction.
#'
#'   A custom measure declares nothing about the standard error of a difference
#'   between two models. That is supplied at comparison time through the
#'   `custom_se_fn` argument of [model_compare()], which accepts a function
#'   `function(ref, cmp) ...`, the shorthands `"sum"` and `"mean"` for the
#'   paired pointwise formulas, or `NULL` to report the difference with an `NA`
#'   standard error. A function receives one list per model with elements
#'   `estimate`, `se`, `pointwise`, and `extra`, always on the measure's natural
#'   scale, and must return the standard error of the difference as a numeric
#'   scalar.
#'
#'   `extra` is for anything the standard error needs that the pointwise values
#'   do not carry. Return it as an additional list element `extra` from the
#'   measure function and it is stored alongside the estimates and passed on to
#'   `custom_se_fn`; the built-in `r2` uses it for the baseline
#'   `(y_i - mean(y))^2`, which cannot be recovered once `y` is out of scope.
#' @param measure_name Only needed when `measure` is a single custom function
#'   passed directly (not inside a list) — it sets the name that function is
#'   reported under. Set it with `attr(my_fun, "measure_name") <- "my_metric"`
#'   before passing `my_fun` to `measure`. If you pass the same function inside
#'   a list instead (e.g. `list(my_metric = my_fun)`), it takes its name from
#'   the list element and `measure_name` isn't needed — this also lets the
#'   same function be reused under several names. If both are set and disagree,
#'   the list name wins and a warning is issued. Either way,
#'   `attr(my_fun, "measure_loss")` is still read from the function itself.
#' @param group_ids Optional vector of group identifiers for grouped summaries
#'   (reserved; not yet implemented).
#' @param loo A [loo::loo()] result, computed with
#'   `save_psis = TRUE` so that PSIS weights are available for additional
#'   measures. See [loo_pred_measure()].
#' @param kfold A `kfold` object from [brms::kfold()]. Supplies ELPD summaries
#'   and fold structure for [kfold_pred_measure()]; pass `y`, `ypred`, and/or
#'   `mupred` when requesting additional measures.
#' @param psis_object A `psis` object with LOO importance weights. An
#'   alternative to passing a full `loo` object; must be supplied together with
#'   `ylp` when computing `elpd`.
#' @param save_psis Logical. If `TRUE`, store the `psis` object in the result
#'   so that additional measures can be added later with [pred_measure()] without
#'   recomputing PSIS weights.
#' @param control Named list of per-measure settings. Each name must match an
#'   element of `measure`; the value is a list of arguments passed to that
#'   measure's summary function (e.g. `list(new_measure = list(add_arg = 10))`).
#' @param source Character string indicating the evaluation mode: `"insample"`,
#'   `"loo"`, `"kfold"`, or `"test"`. Set automatically by the wrapper
#'   functions; required when calling [do_pred_measure()] directly.
#'
#' @keywords internal
#' @name pred_measure_params
do_pred_measure <- function(
  y = NULL,
  ypred = NULL,
  mupred = NULL,
  ylp = NULL,
  ylp_test = NULL,
  measure = NULL,
  predperf = NULL,
  loo = NULL,
  kfold = NULL,
  group_ids = NULL,
  psis_object = NULL,
  save_psis = FALSE,
  source = NULL,
  control = list()
) {
  # input validation ---------------------------------------------------
  if (!is.null(group_ids)) {
    cli::cli_abort(
      "`group_ids` is reserved for future feature but is not yet implemented."
    )
  }

  measures <- .prepare_measures(
    measure, predperf, supported_measures_list, source
  )
  # validated against every requested measure, including the ones
  # `.prepare_measures()` dropped as already present: those are reported by their
  # own warning, and a control entry for them is not a mistake
  .validate_control(control, .normalize_measure(measure))

  if (source == "loo") {
    if (is.null(predperf)) {
      if (!is.null(loo) && is.null(loo$psis_object)) {
        cli::cli_abort(c(
          "No `psis_object` found in `loo` object. Did you run loo(..., save_psis = 'TRUE')."
        ))
      }
      if (is.null(ylp) && !is.null(psis_object)) {
        cli::cli_abort(c(
          "For computation of `elpd` it is required to pass `ylp` besides `psis_object`."
        ))
      }
    } else {
      if (is.null(psis_object) && !is.null(predperf$psis_object)) {
        cli::cli_inform("Using psis_object for LOO CV from `predperf`")
        psis_object <- predperf$psis_object
      }
    }
    checkmate::assert_null(ylp_test)
  } else if (source == "insample") {
    checkmate::assert_null(ylp_test)
    checkmate::assert_null(psis_object)
  } else if (source == "test") {
    checkmate::assert_null(psis_object)
  } else { # kfold
    checkmate::assert_null(psis_object)
    if (is.null(predperf) && is.null(kfold)) {
      cli::cli_abort(c(
        "{.arg kfold} is required for {.fn kfold_pred_measure}.",
        "i" = "Pass a {.cls kfold} object from {.fn brms::kfold}."
      ))
    }
  }

  # core logic ---------------------------------------------------
  if (source == "loo") {
    psis_object <- .get_psis_object(
      loo = loo,
      predperf = predperf,
      psis_object = psis_object,
      ylp = ylp,
      r_eff = 1
    )
  }
  
  base_measure <- .compute_base_measure(
    ylp = ylp,
    ylp_test = ylp_test,
    loo = loo,
    kfold = kfold,
    predperf = predperf,
    psis_object = psis_object,
    source = source
  )

  estimates <- base_measure$estimates
  pointwise <- base_measure$pointwise

  log_weights <- if (!is.null(psis_object)) psis_object$log_weights else NULL

  for (entry in measures) {
    sel_measure <- .compute_measure(
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      measure_entry = entry,
      log_weights = log_weights,
      control = control,
      base_measure = base_measure
    )
    result_name <- attr(sel_measure, "measure")
    if (is.null(result_name)) {
      result_name <- entry$name
    }
    # A measure may rename its own result: `rps` with `scaled = TRUE` returns
    # `srps`. Read the spec under that name, or the requested measure's
    # orientation leaks into the renamed row and inverts the ranking.
    info_entry <- entry
    if (entry$type == "builtin" && !is.null(.measure_spec[[result_name]])) {
      info_entry$key <- result_name
    }
    # add new measures to existing pred_measure results
    name_updated <- .measure_result_name(source, result_name)
    if (!is.null(estimates) && name_updated %in% rownames(estimates)) {
      cli::cli_warn(c(
        "{.field {name_updated}} already present in results. Skipping the update."
      ))
      next
    }
    estimates <- .merge_matrix(
      source = source,
      mat = estimates,
      name = result_name,
      values = .measure_estimate_se(sel_measure),
      margin = 1,
      measure_entry = info_entry,
      extra = sel_measure$extra
    )
    pointwise <- .merge_matrix(
      source = source,
      mat = pointwise,
      name = result_name,
      values = sel_measure$pointwise,
      margin = 2
    )
  }

  predperf_res <- .build_pred_measure(
    estimates = estimates,
    pointwise = pointwise,
    diagnostics = base_measure$diagnostics,
    psis_object = psis_object,
    save_psis = save_psis
  )

  predperf_res <- .add_attributes(
    save_psis = save_psis,
    predperf_res = predperf_res,
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    ylp_test = ylp_test,
    kfold = kfold,
    loo = loo,
    predperf = predperf,
    source = source
  )
  predperf_res
}

# internal helper functions ---------------------------------------------------

#' Resolve or compute the PSIS object for LOO scoring
#'
#' @description
#' Selects a PSIS object from the available inputs, or computes one from `ylp`
#' when no precomputed weights are supplied (assuming `r_eff = 1`).
#'
#' Resolution order:
#' \enumerate{
#'   \item If both `psis_object` and `loo` are provided, return `psis_object`
#'     after verifying it matches `loo$psis_object`.
#'   \item Extract from `loo$psis_object` when `loo` is provided.
#'   \item Use the supplied `psis_object`.
#'   \item Reuse `predperf$psis_object` when accumulating measures.
#'   \item Reuse `predperf$log_weights` when only the weights were stored
#'     (`save_psis = FALSE`); the result carries no `diagnostics`.
#'   \item Compute from `ylp` via [loo::psis()] on `-ylp` log ratios.
#' }
#'
#' @param ylp Matrix of pointwise log predictive densities (`S` × `n`).
#' @param loo Optional [loo::loo()] result containing a `psis_object`.
#' @param predperf Optional existing measure object with a stored `psis_object`.
#' @param psis_object Optional precomputed PSIS object.
#' @param r_eff Relative effective sample size passed to [loo::psis()];
#'   default `1`.
#'
#' @return A `psis` object with `log_weights` and `diagnostics`.
#'
#' @note See developer notes on computation of the `psis_object` for details.
#' @noRd
.get_psis_object <- function(
  ylp,
  loo,
  predperf,
  psis_object,
  r_eff = 1
) {
  # psis_object + loo are both provided -> return psis_object
  if (!is.null(psis_object) && !is.null(loo)) {
    psis_equal_loo <- isTRUE(all.equal(
      psis_object$log_weights,
      loo$psis_object$log_weights
    ))
    if (!psis_equal_loo) {
      cli::cli_abort(
        "Provided `psis_object` and `loo$psis_object` are not identical."
      )
    } 
    return(psis_object)
  # loo is provided
  } else if (!is.null(loo)) {
    return(loo$psis_object)
  # psis_object is provided
  } else if (!is.null(psis_object)) {
    return(psis_object)
  # predperf with psis_object is provided
  } else if (!is.null(predperf$psis_object)) {
    return(predperf$psis_object)
  # predperf carries log_weights only (save_psis = FALSE)
  } else if (!is.null(predperf$log_weights)) {
    return(list(log_weights = predperf$log_weights))
  # ylp is provided
  } else if (is.null(loo) && is.null(psis_object) && !is.null(ylp)) {
    cli::cli_inform(
      "Compute `psis_object` internally from `ylp` assuming `r_eff = 1`."
    )
    log_ratios <- -1 * ylp
    return(psis(log_ratios, r_eff = r_eff))
  # nothing is provided
  } else {
    cli::cli_abort(c(
      "psis_object can not be computed, either of `psis_object`, `loo`, or",
      "`ylp` needs to be provided."
    ))
  }}

#' Extract estimate and SE from a measure result
#'
#' `.create_measure_structure()` gives every builtin measure an `estimates`
#' matrix. A custom measure returns either `estimates` or `estimate` and `se`.
#'
#' @noRd
.measure_estimate_se <- function(res) {
  if (!is.null(res$estimates)) {
    return(res$estimates)
  }
  c(res$estimate, res$se)
}

#' Pointwise ELPD from the base measure block
#'
#' `.compute_base_measure()` appends the source suffix to the base column names,
#' so the column is one of `elpd`, `elpd_loo`, `elpd_kfold`, `elpd_test`.
#' Measures that declare `needs_elpd` in `.measure_spec` read the log predictive
#' density here instead of recomputing it from `ylp`. `ylp` is absent when only
#' a `loo` or `kfold` object is supplied, and on the `test` source it holds the
#' training data while the base block holds the holdout data.
#'
#' @noRd
.base_elpd_pointwise <- function(base_measure) {
  known <- c("elpd", "elpd_loo", "elpd_kfold", "elpd_test")
  hit <- intersect(known, colnames(base_measure$pointwise))
  if (length(hit) == 0L) {
    cli::cli_abort(c(
      "No {.field elpd} column found in the base measure.",
      "i" = "{.val mlpd} and {.val ic} are derived from {.val elpd}."
    ))
  }
  base_measure$pointwise[, hit[1L]]
}

#' Compute a single predictive measure
#'
#' @description
#' Dispatches one requested predictive measure to the appropriate summary
#' function. No measure has a fixed input list. The function builds a pool of
#' candidate arguments, then keeps only the ones the measure function declares:
#'
#'   `args <- pool[intersect(names(formals(measure_fun)), names(pool))]`
#'
#' The pool holds `y`, `ypred`, `mupred`, `ylp` and `log_weights`, plus the
#' measure's slice of `control`. A measure that sets `needs_elpd` in
#' `.measure_spec` gets a different pool: `pointwise` holds the ELPD column
#' taken from `base_measure`, and `ylp` and `log_weights` are `NULL`.
#'
#' @param y Vector of observed values (n).
#' @param ypred Matrix of posterior predictive draws (S × n).
#' @param mupred Matrix of posterior point predictions (S × n).
#' @param ylp Matrix of pointwise log predictive densities (S × n).
#' @param measure_entry A normalized measure entry with elements 
#'   `name`, `type` (`"builtin"` or `"custom"`), and `key`.
#' @param log_weights Matrix of log-weights (S × n), as returned by 
#'   `.compute_log_weights()`.
#' @param control Named list of per-measure settings passed from
#'   [pred_measure()]; the active slice is `control[[measure_entry$name]]`.
#' @param base_measure The base measure block from `.compute_base_measure()`.
#'   Read only when the measure sets `needs_elpd` in `.measure_spec`.
#'
#' @return The result of the measure function, in one of two shapes.
#'   \describe{
#'     \item{builtin}{A `"measure"` object from `.create_measure_structure()`:
#'       `estimates`, a 1 by 2 matrix with columns `Estimate` and `SE`, and
#'       `pointwise`, an n by 1 matrix.}
#'     \item{custom}{Either `estimates`, a length-2 numeric vector, or
#'       `estimate` and `se` as scalars. `pointwise` is a numeric vector.
#'       `.validate_measure_result()` accepts both.}
#'   }
#'
#' @noRd
.compute_measure <- function(
  y,
  ypred,
  mupred,
  ylp,
  measure_entry,
  log_weights,
  control = list(),
  base_measure
) {
  if (measure_entry$type == "builtin") {
    spec <- .measure_spec[[measure_entry$key]]
    if (is.null(spec)) {
      cli::cli_abort("Unknown built-in measure {.val {measure_entry$key}}.")
    }
    measure_fun <- spec$fun
  } else {
    spec <- NULL
    measure_fun <- measure_entry$key
  }

  measure_control <- control[[measure_entry$name]]
  if (is.null(measure_control)) {
    measure_control <- list()
  }

  pool <- if (isTRUE(spec$needs_elpd)) {
    lppd_i <- .base_elpd_pointwise(base_measure)
    if (is.function(spec$elpd_transform)) {
      lppd_i <- spec$elpd_transform(lppd_i)
    }
    # `ylp` has no default, and `.inform_ignored_inputs()` forces it. Pass it as
    # NULL rather than leaving it out.
    list(ylp = NULL, log_weights = NULL, pointwise = lppd_i)
  } else {
    list(
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      log_weights = log_weights
    )
  }

  pool <- c(pool, measure_control)
  args <- pool[intersect(names(formals(measure_fun)), names(pool))]
  res <- do.call(measure_fun, args)
  if (measure_entry$type == "custom") {
    n_obs <- .measure_n_obs(y, ypred, mupred, ylp)
    res <- .validate_measure_result(res, measure_entry$name, n_obs = n_obs)
  }
  res
}

#' Compute base density summaries for a predictive measure object
#'
#' @description
#' Forms the default block of log-density summaries (`elpd` and for loo and kfold
#' the effective number of parameters `p_loo`\`p_kfold`) that underlie every 
#' [pred_measure()] result. Additional measures
#' requested via `measure` are merged into the returned matrices later.
#'
#' When `predperf` is provided (incremental update), returns `predperf` unchanged.
#' For `source = "kfold"` or `"loo"` with a precomputed object, extracts existing
#' summaries from `kfold` or `loo`.
#'
#' Otherwise ELPD is computed from:
#' \itemize{
#'   \item `insample`: `ylp` directly (in-sample log predictive density).
#'   \item `loo`: `ylp` reweighted with PSIS log weights from `psis_object`.
#'   \item `test`: `ylp_test` on holdout observations.
#' }
#'
#' Effective number of parameters is computed by \code{.compute_effective_param()} 
#' as the difference between in-sample and LOO log predictive density; see 
#' `p_loo` in [loo::loo()] and the [CV-FAQ on p_loo](https://users.aalto.fi/~ave/CV-FAQ.html#p_loo).
#'
#' @param ylp Matrix of pointwise log predictive densities for training data
#'   (`S` × `n`).
#' @param ylp_test Matrix of pointwise log predictive densities for holdout
#'   data (`S` × `n_test`; `test` source only).
#' @param loo Optional [loo::loo()] result. When supplied for `source = "loo"`,
#'   base summaries are taken from the `loo` object.
#' @param kfold Optional `kfold` object. When supplied for `source = "kfold"`,
#'   base summaries are taken from the `kfold` object.
#' @param predperf Existing [pred_measure()] object. When not `NULL`, returned
#'   unchanged so base summaries are not recomputed.
#' @param psis_object PSIS object with LOO log weights and diagnostics.
#' @param source Character string; one of `"insample"`, `"loo"`, `"kfold"`, or
#'   `"test"`.
#'
#' @return A named list with:
#' \describe{
#'   \item{`estimates`}{Matrix with rows `elpd` and optionally `p` (number of 
#'   effective parameters, for LOO and kfold), plus a source suffix when 
#'   applicable (`_loo`, `_kfold`, `_test`).}
#'   \item{`pointwise`}{Matrix of observation-level contributions.}
#'   \item{`diagnostics`}{From `psis_object$diagnostics` when LOO weights are
#'     used; otherwise `NULL` or taken from the input `loo`/`kfold` object.}
#' }
#'
#' @noRd
.compute_base_measure <- function(
  ylp,
  ylp_test,
  loo,
  kfold,
  predperf,
  psis_object,
  source
) {
  if (!is.null(predperf)) return(predperf)
  
  if (source == "kfold") {
    components <- if ("diagnostics" %in% names(kfold)) {
      c("estimates", "pointwise", "diagnostics")
    } else {
      c("estimates", "pointwise")
    }
    return(subset_measures(
      kfold, 
      measures = c("elpd_kfold", "p_kfold"),
      components = components
    ))
  }
  
  if (source == "loo" && !is.null(loo)) {
    return(subset_measures(
      loo, 
      measures = c("elpd_loo", "p_loo"),
      components = c("estimates", "pointwise", "diagnostics")
    ))
  }
  
  elpd_res <- switch(source,
    insample = measure_elpd(ylp = ylp),
    loo = measure_elpd(ylp = ylp, log_weights = psis_object$log_weights),
    test = measure_elpd(ylp_test)
  )
  
  est_se <- .measure_estimate_se(elpd_res)
  elpd_res <- list(
    estimate = unname(est_se[1]),
    se = unname(est_se[2]),
    pointwise = elpd_res$pointwise
  )
  
  suffix <- if (source == "insample") "" else paste0("_", source)
  add_p_eff <- source == "loo"
  
  estimates <- rbind(elpd = c(elpd_res$estimate, elpd_res$se))
  pointwise <- cbind(elpd = elpd_res$pointwise)
  
  if (add_p_eff) {
    p_loo <- .compute_effective_param(ylp, elpd_res$pointwise)
    estimates <- rbind(estimates, p = c(p_loo$estimate, p_loo$se))
    pointwise <- cbind(pointwise, p = p_loo$pointwise)
  }
  
  rownames(estimates) <- paste0(rownames(estimates), suffix)
  colnames(estimates) <- c("Estimate", "SE")
  colnames(pointwise) <- paste0(colnames(pointwise), suffix)
  
  list(
    estimates = estimates, 
    pointwise = pointwise,
    diagnostics = psis_object$diagnostics
  )
}

# compute the effective number of parameters
#' Compute effective number of parameters (`p_loo`)
#'
#' @description
#' Per-observation effective number of parameters as the difference between
#' the log posterior predictive density (`lpd`) and the cross-validated log
#' predictive density (`elpd`). Summed across observations, this matches
#' `p_loo` from the **loo** package: it describes how much harder it is to
#' predict held-out data than the data used to fit the model.
#'
#' @param ylp Matrix of pointwise log predictive densities (`S` × `n`).
#' @param elpd_cv_i Numeric vector of length `n` with LOO pointwise ELPD
#'   contributions.
#'
#' @return A named list with `estimate` (total `p_loo`), `se`, and `pointwise`
#'   (`p_loo` per observation).
#'
#' @references See the **loo** package glossary (`vignette("loo2", package = "loo")`)
#'   and \url{https://users.aalto.fi/~ave/CV-FAQ.html#p_loo}.
#'
#' @noRd
.compute_effective_param <- function(ylp, elpd_cv_i) {
  lpd_i <- matrixStats::colLogSumExps(ylp) - log(nrow(ylp))
  p_eff_i <- lpd_i - elpd_cv_i[ ,"elpd"]
  
  list(
      estimate = sum(p_eff_i),
      se = sqrt(ncol(ylp) * var(p_eff_i)),
      pointwise = p_eff_i
    )
}

#' Add or update a row or column in a summary matrix
#'
#' @description
#' Builds the `estimates` and `pointwise` matrices used in
#' [pred_measure()] results. When `margin = 1`, appends or updates a **row**
#' (for estimates). When `margin = 2`, appends or updates a **column** 
#' (for pointwise).
#'
#' If `mat` is `NULL`, returns a new one-row or one-column matrix for `name`.
#' If `name` is already present along that margin, the update is skipped with a
#' warning. Otherwise the new slice is bound with [rbind()] or [cbind()].
#'
#' @param source Character string; one of `"loo"`, `"insample"`, `"kfold"`, or
#'   `"test"`.
#' @param mat Existing matrix, or `NULL` when adding the first measure.
#' @param name Character label for the row (`margin = 1`) or column
#'   (`margin = 2`).
#' @param values Numeric vector to store. For `margin = 1`, length-2 vector
#'   `(estimate, se)`; for `margin = 2`, length-`n` pointwise vector.
#' @param margin `1` to merge along rows (estimates table), `2` along columns
#'   (pointwise table).
#' @param measure_entry Optional normalized measure entry; when merging an
#'   estimates row (`margin = 1`), the `measure_info` used by [model_compare()]
#'   is recorded from this entry.
#' @param extra Optional list of auxiliary data the measure stores for its
#'   `se_diff_fun` (the measure result's `extra` element); recorded in `measure_info`
#'   when merging an estimates row (`margin = 1`).
#'
#' @return Updated matrix with `name` as a row or column name.
#'
#' @noRd
.measure_result_name <- function(source, name) {
  switch(
    source,
    kfold = paste0(name, "_kfold"),
    loo = paste0(name, "_loo"),
    test = paste0(name, "_test"),
    insample = name
  )
}

#' @noRd
.merge_matrix <- function(
  source,
  mat,
  name,
  values,
  margin,
  measure_entry = NULL,
  extra = NULL
) {
  is_row <- margin == 1
  bind_fn <- if (is_row) rbind else cbind
  name_updated <- .measure_result_name(source, name)

  new_slice <- if (is_row) {
    matrix(values, nrow = 1, dimnames = list(name_updated, c("Estimate", "SE")))
  } else {
    matrix(values, ncol = 1, dimnames = list(NULL, name_updated))
  }

  info <- if (is_row && !is.null(measure_entry)) {
    .measure_info(measure_entry)
  }
  if (!is.null(info) && !is.null(extra)) {
    info$extra <- extra
  }

  old_info <- if (is_row && !is.null(mat)) {
    attr(mat, "measure_info")
  }

  mat <- if (is.null(mat)) new_slice else bind_fn(mat, new_slice)

  if (is_row && (!is.null(info) || !is.null(old_info))) {
    measure_info <- old_info
    if (is.null(measure_info)) {
      measure_info <- list()
    }
    if (!is.null(info)) {
      measure_info[[name]] <- info
    }
    attr(mat, "measure_info") <- measure_info
  }

  mat
}

#' Construct the S3 predictive measure result object
#'
#' @description
#' Wraps computed summaries into the list structure returned by the
#' [pred_measure()] pipeline. S3 classes and attributes are attached later by
#' \code{.add_attributes()}.
#'
#' When `save_psis = TRUE`, the `psis_object` is stored in the result; otherwise
#' that slot is omitted. When a `psis_object` is available, its `log_weights` are
#' copied to the result.
#'
#' @param estimates Matrix of overall estimates and standard errors (rows =
#'   measures, columns = `Estimate` and `SE`).
#' @param pointwise Matrix of observation-level contributions (columns =
#'   measures).
#' @param diagnostics Optional PSIS or other diagnostic information, or `NULL`.
#' @param psis_object PSIS object with LOO weights and diagnostics, or `NULL`.
#' @param save_psis Logical; if `TRUE`, include `psis_object` in the result.
#'
#' @return A list with elements `estimates`, `pointwise`, and optionally
#'   `diagnostics`, `psis_object`, and `log_weights`. Attribute `measure_info`
#'   records per-measure metadata for measures added in the current call. Class
#'   attributes are added by \code{.add_attributes()}.
#'
#' @noRd
.build_pred_measure <- function(
  estimates,
  pointwise,
  diagnostics = NULL,
  psis_object,
  save_psis
) {
  measure_info <- attr(estimates, "measure_info")
  if (is.null(measure_info)) {
    measure_info <- list()
  }
  attr(estimates, "measure_info") <- NULL

  output_list <- list(
    estimates = estimates,
    pointwise = pointwise
  )
  if (!is.null(diagnostics)) {
    output_list$diagnostics <- diagnostics
  }
  if (isTRUE(save_psis)) {
    output_list$psis_object <- psis_object
  }
  if (!is.null(psis_object)) {
    output_list$log_weights <- psis_object$log_weights
  }

  structure(
    output_list,
    measure_info = measure_info
  )
}

#' Attach S3 classes and metadata attributes to a result
#'
#' @description
#' Sets `class`, `source`, `dims`, and `measure_info` attributes on a predictive
#' measure object.
#'
#' When updating an existing result (`predperf` is not `NULL`), copies attributes
#' from `predperf` and refreshes `dims` from newly supplied input matrices.
#' Merges `measure_info` from the prior result with any new entries supplied on
#' `predperf_res` (from \code{.build_pred_measure()}).
#' When `save_psis = FALSE`, clears any stored `psis_object` from the prior
#' result.
#'
#' For new objects, copies relevant attributes from `loo` or `kfold` inputs
#' (e.g. `yhash`, `model_name`, fold structure) and assigns a source-specific
#' subclass (`"insample_pred_measure"`, `"loo_pred_measure"`, etc.). Sets
#' `measure_info`, seeding the `elpd` entry.
#'
#' @param save_psis Logical; when `FALSE` and accumulating, clears stored
#'   `psis_object` from the prior result.
#' @param predperf_res List returned by \code{.build_pred_measure()}.
#' @param y Vector of observed values; used indirectly via matrix `dims`.
#' @param ypred Matrix of posterior predictive draws; used to set `dims`.
#' @param mupred Matrix of posterior expected values; used to set `dims`.
#' @param ylp Matrix of pointwise log predictive densities; used to set `dims`.
#' @param ylp_test Matrix of holdout log predictive densities; sets `dims` for
#'   `source = "test"`.
#' @param kfold Optional `kfold` object whose attributes are inherited.
#' @param loo Optional [loo::loo()] object whose attributes are inherited.
#' @param predperf Existing object when accumulating measures, or `NULL`.
#' @param source Character evaluation mode (`"insample"`, `"loo"`, `"kfold"`,
#'   or `"test"`).
#'
#' @return The updated `predperf_res` with class and attributes set.
#'
#' @noRd
.add_attributes <- function(
  save_psis,
  predperf_res,
  y,
  ypred,
  mupred,
  ylp,
  ylp_test,
  kfold,
  loo,
  predperf,
  source
) {
  new_info <- attr(predperf_res, "measure_info")
  if (is.null(new_info)) {
    new_info <- list()
  }

  if (!is.null(predperf)) {
    if (isFALSE(save_psis)) {
      predperf$psis_object <- NULL
    }
    attributes(predperf_res) <- attributes(predperf)

    dims <- if (!is.null(ypred)) {
      dim(ypred)
    } else if (!is.null(mupred)) {
      dim(mupred)
    } else if (!is.null(ylp)) {
      dim(ylp)
    } else {
      attr(predperf, "dims")
    }
    attr(predperf_res, "dims") <- dims
    measure_info <- attr(predperf, "measure_info")
    if (is.null(measure_info)) {
      measure_info <- list()
    }
    if (is.null(measure_info$elpd)) {
      measure_info$elpd <- .measure_info("elpd")
    }
    if (length(new_info)) {
      measure_info[names(new_info)] <- new_info
    }
    attr(predperf_res, "measure_info") <- measure_info

    return(predperf_res)
  }

  predperf_res <- switch(
    source,
    kfold = .copy_attrs(
      predperf_res,
      kfold,
      setdiff(names(attributes(kfold)), "names")
    ),
    loo = .copy_attrs(
      predperf_res,
      loo,
      setdiff(names(attributes(loo)), "names")
    ),
    test = , # fall through (same as insample)
    insample = predperf_res
  )
  
  if (source %in% c("insample", "test") || 
    (is.null(attr(predperf_res, "dims")) && !is.null(ylp))){
    # make attribute names consistent between pred_measure classes
    if (source == "test") {
      attr(predperf_res, "dims") <- attr(ylp_test, "dim")
    } else {
      attr(predperf_res, "dims") <- attr(ylp, "dim")  
    } 
  }
  
  classes <- c(
    switch(
      source,
      loo = "loo_pred_measure",
      insample = "insample_pred_measure",
      kfold = "kfold_pred_measure",
      test = "test_pred_measure",
      "pred_measure"
    ),
    "pred_measure",
    attr(predperf_res, "class")
  )
  if (source == "loo" && !"loo" %in% classes) {
    classes <- c(classes, "loo")
  }
  attr(predperf_res, "class") <- classes
  attr(predperf_res, "source") <- source
  measure_info <- list(
    elpd = .measure_info("elpd")
  )
  if (length(new_info)) {
    measure_info[names(new_info)] <- new_info
  }
  attr(predperf_res, "measure_info") <- measure_info

  predperf_res
}
