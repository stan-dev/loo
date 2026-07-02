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
#' @param measure Additional measures beyond the base summaries (`elpd` and
#'   `ic`, which are always included). Can be:
#'   \itemize{
#'     \item A **character vector** of built-in names; see
#'       [supported_measures_list()].
#'     \item A **function** with attribute `"measure_name"` for one custom measure.
#'     \item A **list** mixing character scalars (built-in names) and named
#'       functions (custom measures), e.g. `list("rps", my_metric = my_fun)`.
#'   }
#'   Custom functions are called with any of `y`, `ypred`, `mupred`, `ylp`, and
#'   `log_weights` that appear in their formals, plus arguments from `control`.
#'   They must return a list with  `estimates` and `pointwise`.
#' @param measure_name For a single custom function, set
#'   `attr(my_fun, "measure_name") <- "my_metric"` before passing `my_fun` to
#'   `measure`.
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
  .validate_control(control)
  
  measures <- .prepare_measures(measure, predperf, supported_measures_list)
  
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
    estimates <- .merge_matrix(
      source = source,
      mat = estimates,
      name = entry$name,
      values = .measure_estimate_se(sel_measure),
      margin = 1
    )
    pointwise <- .merge_matrix(
      source = source,
      mat = pointwise,
      name = entry$name,
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
  
  .add_attributes(
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
  )
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
    psis_equal_loo <- all(psis_object == loo$psis_object)
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

#' Compute a single predictive measure
#'
#' @description
#' Dispatches one requested predictive measure to the appropriate summary
#' function. The measure's `family` in `.measure_spec` determines which inputs
#' are passed through:
#'
#' \describe{
#'   \item{`metrics`}{`y` and `mupred` (e.g. MAE, RMSE, accuracy).}
#'   \item{`rank_scores`}{`y` and `ypred` (e.g. RPS, CRPS).}
#'   \item{`density_scores`}{`ylp` (e.g. ELPD, MLPD).}
#' }
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
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{`estimate`}{Scalar point estimate for the measure.}
#'     \item{`se`}{Scalar standard error for the estimate.}
#'     \item{`pointwise`}{Length-`n` vector of observation-level contributions.}
#'   }
#'
#' Extract estimate and SE from a measure result
#'
#' Supports `estimate`/`se` (pred_measure functions) and `estimates` (CRPS).
#'
#' @noRd
.measure_estimate_se <- function(res) {
  if (!is.null(res$estimates)) {
    return(res$estimates)
  }
  c(res$estimate, res$se)
}

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
    measure_fun <- .measure_spec[[measure_entry$key]]
    if (is.null(measure_fun)) {
      cli::cli_abort("Unknown built-in measure {.val {measure_entry$key}}.")
    }
  } else {
    measure_fun <- measure_entry$key
  }
  
  if (identical(measure_fun, measure_mlpd) || identical(measure_fun, measure_ic)) {
    elpd_i <- base_measure$pointwise[, 1]
    # ensure elpd_i is in correct orientation (negative log predictive density)
    if (any(elpd_i > 0)) elpd_i <- -elpd_i

    if (identical(measure_fun, measure_mlpd)) {
      return(measure_mlpd(ylp = NULL, pointwise = elpd_i))
    }
    return(measure_ic(ylp = NULL, pointwise = -2 * elpd_i))
  }

  measure_control <- control[[measure_entry$name]]
  if (is.null(measure_control)) {
    measure_control <- list()
  }

  pool <- c(
    list(
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      log_weights = log_weights
    ),
    measure_control
  )
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
    components <- if ("diagnostics" %in% kfold1) {
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
    insample = elpd(ylp = ylp),
    loo = elpd(ylp = ylp, log_weights = psis_object$log_weights),
    test = elpd(ylp_test)
  )
  
  if (!is.null(elpd_res$estimates)) {
    elpd_res <- list(
      estimate = unname(elpd_res$estimates[1]),
      se = unname(elpd_res$estimates[2]),
      pointwise = elpd_res$pointwise
    )
  }
  
  suffix <- if (source == "insample") "" else paste0("_", source)
  add_p_eff <- source == "loo"
  
  estimates <- rbind(elpd = c(elpd_res$estimate, elpd_res$se))
  pointwise <- cbind(elpd = elpd_res$pointwise)
  
  if (add_p_eff) {
    p_eff <- .compute_effective_param(ylp, elpd_res$pointwise)
    estimates <- rbind(estimates, p_eff = c(p_eff$estimate, p_eff$se))
    pointwise <- cbind(pointwise, p_eff = p_eff$pointwise)
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
#'
#' @return Updated matrix with `name` as a row or column name.
#'
#' @noRd
.merge_matrix <- function(source, mat, name, values, margin) {
  is_row <- margin == 1
  along <- if (is_row) rownames else colnames
  bind_fn <- if (is_row) rbind else cbind
  name_updated <- switch(
    source,
    kfold = paste0(name, "_kfold"),
    loo = paste0(name, "_loo"),
    test = paste0(name, "_test"),
    insample = name
  )
  
  new_slice <- if (is_row) {
    matrix(values, nrow = 1, dimnames = list(name_updated, c("Estimate", "SE")))
  } else {
    matrix(values, ncol = 1, dimnames = list(NULL, name_updated))
  }

  if (is.null(mat)) return(new_slice)
  # && margin == 1 condition ensures that warning is shown only once
  if (name_updated %in% along(mat) && margin == 1) {
    # TODO: Check whether this behavior is wanted
    cli::cli_warn(
      c(
        "{.field {name_updated}} already present in results. Skipping the update."
      )
    )
    return(mat)
  }
  bind_fn(mat, new_slice)
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
#'   `diagnostics`, `psis_object`, and `log_weights`. Class attributes are added
#'   by \code{.add_attributes()}.
#'
#' @noRd
.build_pred_measure <- function(
  estimates,
  pointwise,
  diagnostics = NULL,
  psis_object,
  save_psis
) {
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
  
  structure(output_list)
}

#' Attach S3 classes and metadata attributes to a result
#'
#' @description
#' Sets `class`, `source`, and `dims` attributes on a predictive measure object.
#'
#' When updating an existing result (`predperf` is not `NULL`), copies attributes
#' from `predperf` and refreshes `dims` from newly supplied input matrices.
#' When `save_psis = FALSE`, clears any stored `psis_object` from the prior
#' result.
#'
#' For new objects, copies relevant attributes from `loo` or `kfold` inputs
#' (e.g. `yhash`, `model_name`, fold structure) and assigns a source-specific
#' subclass (`"insample_pred_measure"`, `"loo_pred_measure"`, etc.).
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
.add_attributes <- function(save_psis, predperf_res, y, ypred, mupred, ylp, ylp_test, kfold, loo, predperf, source) {
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
    }
    attr(predperf_res, "dims") <- dims
    
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
  
  attr(predperf_res, "class") <- c(
    "pred_measure",
    switch(
      source,
      loo = "loo_pred_measure",
      insample = "insample_pred_measure",
      kfold = "kfold_pred_measure",
      test = "test_pred_measure",
      "pred_measure"
    ),
    attr(predperf_res, "class")
  )
  attr(predperf_res, "source") <- source
  
  return(predperf_res)
}