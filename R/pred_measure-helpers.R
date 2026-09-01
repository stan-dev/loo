#' Normalize the `measure` argument to an internal list
#'
#' @description
#' Converts `measure` (character, function, list, or `NULL`) into a list of
#' entries with elements `name`, `type` (`"builtin"` or `"custom"`), and `key`
#' (built-in name or function). Custom entries also carry `loss`, taken from
#' `attr(fun, "measure_loss")`.
#'
#' @param measure User-supplied `measure` argument.
#'
#' @return A list of normalized measure entries, or an empty list when
#'   `measure` is `NULL`.
#'
#' @noRd
.normalize_measure <- function(measure) {
  if (is.null(measure)) {
    return(list())
  }

  if (is.function(measure)) {
    entries <- list(.measure_entry_custom(measure))
    .check_measure_entry_names(entries)
    return(entries)
  }

  if (is.character(measure)) {
    entries <- lapply(measure, function(nm) {
      list(name = nm, type = "builtin", key = nm)
    })
    .check_measure_entry_names(entries)
    return(entries)
  }

  if (is.list(measure)) {
    if (length(measure) == 0L) {
      return(list())
    }
    entries <- lapply(seq_along(measure), function(i) {
      el <- measure[[i]]
      nm <- names(measure)[i]
      if (is.character(el) && length(el) == 1L) {
        list(name = el, type = "builtin", key = el)
      } else if (is.function(el)) {
        if (is.null(nm) || !nzchar(nm)) {
          cli::cli_abort(c(
            "Each custom function in {.arg measure} must be named.",
            "i" = "Use {.code measure = list(my_metric = my_fun)}."
          ))
        }
        .measure_entry_custom(el, name = nm)
      } else {
        cli::cli_abort(c(
          "Each element of {.arg measure} must be a character scalar (built-in",
          "name) or a function (custom measure).",
          "i" = "Element {i} has type {.cls {class(el)[1]}}."
        ))
      }
    })
    .check_measure_entry_names(entries)
    return(entries)
  }

  cli::cli_abort(c(
    "{.arg measure} must be a character vector, a function, a list, or",
    "{.code NULL}.",
    "i" = "Got an object of class {.cls {class(measure)[1]}}."
  ))
}

#' Build a custom measure entry from a function
#'
#' The name and the orientation of a measure are properties of its definition,
#' so both are declared once on the function via `attr(fun, "measure_name")`
#' and `attr(fun, "measure_loss")` rather than at every call. A custom measure
#' is a utility (higher is better) unless it declares itself a loss.
#'
#' @param fun Function implementing a custom measure.
#' @param name Measure name, when it comes from the name of a `measure` list
#'   element rather than from `attr(fun, "measure_name")`.
#' @noRd
.measure_entry_custom <- function(fun, name = NULL) {
  if (is.null(name)) {
    name <- attr(fun, "measure_name", exact = TRUE)
    if (is.null(name) || length(name) != 1L || !nzchar(name)) {
      stop(
        "A custom function passed to 'measure' must have attribute ",
        "'measure_name', e.g. attr(my_fun, \"measure_name\") <- \"my_metric\".",
        call. = FALSE
      )
    }
  } else {
    attr_name <- attr(fun, "measure_name", exact = TRUE)
    if (!is.null(attr_name) && length(attr_name) == 1L && nzchar(attr_name) &&
        !identical(attr_name, name)) {
      cli::cli_warn(c(
        "Custom measure named {.val {name}} in {.arg measure} also has",
        "{.code attr(fun, \"measure_name\") = {.val {attr_name}}}.",
        "i" = "Using the list name {.val {name}}; the attribute is ignored here."
      ))
    }
  }

  loss <- attr(fun, "measure_loss", exact = TRUE)
  if (is.null(loss)) {
    loss <- FALSE
  } else if (!is.logical(loss) || length(loss) != 1L || is.na(loss)) {
    stop(
      "Attribute 'measure_loss' of a custom measure must be TRUE or FALSE, ",
      "e.g. attr(my_fun, \"measure_loss\") <- TRUE for a measure where lower ",
      "values are better.",
      call. = FALSE
    )
  }

  list(name = name, type = "custom", key = fun, loss = loss)
}

#' Check duplicate measure names
#'
#' @param entries List of normalized measure entries from `.normalize_measure()`.
#' @noRd
.check_measure_entry_names <- function(entries) {
  names <- vapply(entries, `[[`, "", "name")
  dups <- names[duplicated(names)]
  if (length(dups) > 0L) {
    cli::cli_abort(c(
      "Duplicate measure names in {.arg measure}: {.val {unique(dups)}}",
      "i" = "Each measure may appear only once."
    ))
  }
  invisible(NULL)
}

#' Normalize, validate, and filter the `measure` argument
#'
#' @description
#' Converts `measure` via `.normalize_measure()`, validates built-in names and
#' custom functions, and drops entries already present in `predperf`.
#'
#' @param measure User-supplied `measure` argument (see `.normalize_measure()`).
#' @param predperf Existing pred_measure object used when accumulating measures.
#' @param supported_measures_list Character vector of allowed built-in names.
#' @param source Character evaluation mode (`"insample"`, `"loo"`, `"kfold"`,
#'   or `"test"`); selects the row-name suffix used to match `predperf`.
#'
#' @return A list of normalized measure entries ready for computation.
#'
#' @noRd
.prepare_measures <- function(
  measure, predperf, supported_measures_list, source
) {
  entries <- .normalize_measure(measure)
  if (length(entries) == 0L) {
    return(entries)
  }

  is_builtin <- vapply(entries, function(e) e$type == "builtin", logical(1L))
  builtin_keys <- vapply(entries[is_builtin], function(e) e$key, character(1L))
  invalid <- setdiff(builtin_keys, supported_measures_list)
  if (length(invalid) > 0L) {
    cli::cli_abort(c(
      "Invalid measure{?s}: {paste(shQuote(invalid), collapse = ', ')}",
      "i" = "Built-in measures must be one of:",
      " " = "{paste(shQuote(supported_measures_list), collapse = ', ')}"
    ))
  }

  is_custom <- !is_builtin
  if (any(is_custom)) {
    for (entry in entries[is_custom]) {
      checkmate::assert_function(entry$key, .var.name = "measure")
    }
  }

  if (!is.null(predperf)) {
    existing_measures <- rownames(predperf$estimates)
    entry_names <- vapply(
      entries,
      function(e) .measure_result_name(source, e$name),
      character(1L)
    )
    dups <- intersect(entry_names, existing_measures)
    if (length(dups) > 0L) {
      cli::cli_warn(c(
        "!" = "{cli::qty(length(dups))} Measure{?s} {.val {dups}} {?is/are}",
        "already present in {.arg predperf} and will be skipped."
      ))
    }
    keep <- !(entry_names %in% existing_measures)
    entries <- entries[keep]
  }

  entries
}

#' Infer number of observations from measure inputs
#'
#' @noRd
.measure_n_obs <- function(y, ypred, mupred, ylp) {
  if (!is.null(y)) {
    return(length(y))
  }
  if (!is.null(ypred)) {
    return(ncol(ypred))
  }
  if (!is.null(mupred)) {
    return(if (length(dim(mupred)) == 3L) dim(mupred)[2L] else ncol(mupred))
  }
  if (!is.null(ylp)) {
    return(ncol(ylp))
  }
  NULL
}

#' Validate the return value of a custom measure function
#'
#' @param res Object returned by a custom measure function.
#' @param measure_name Label used in error messages.
#' @param n_obs Expected length of `pointwise`, or `NULL` to skip.
#'
#' @return `res`, invisibly.
#'
#' @noRd
.validate_measure_result <- function(res, measure_name, n_obs = NULL) {
  if (!is.list(res)) {
    cli::cli_abort(c(
      "Custom measure {.val {measure_name}} must return a list.",
      "i" = "Got an object of class {.cls {class(res)[1]}}."
    ))
  }

  if (!is.null(res$estimates)) {
    if (!is.numeric(res$estimates) || length(res$estimates) != 2L) {
      cli::cli_abort(c(
        "{.field estimates} from custom measure {.val {measure_name}} must be",
        "a numeric vector of length 2 (estimate and SE)."
      ))
    }
  } else {
    missing <- setdiff(c("estimate", "se", "pointwise"), names(res))
    if (length(missing) > 0L) {
      cli::cli_abort(c(
        "Custom measure {.val {measure_name}} must return a list with",
        "{.field estimate}, {.field se}, and {.field pointwise}.",
        "x" = "Missing: {.field {missing}}"
      ))
    }
    if (!is.numeric(res$estimate) || length(res$estimate) != 1L) {
      cli::cli_abort(
        "{.field estimate} from custom measure {.val {measure_name}} must be a numeric scalar."
      )
    }
    if (!is.numeric(res$se) || length(res$se) != 1L) {
      cli::cli_abort(
        "{.field se} from custom measure {.val {measure_name}} must be a numeric scalar."
      )
    }
  }

  if (is.null(res$pointwise)) {
    cli::cli_abort(
      "Custom measure {.val {measure_name}} must return {.field pointwise}."
    )
  }
  if (!is.numeric(res$pointwise) || length(res$pointwise) < 1L) {
    cli::cli_abort(
      "{.field pointwise} from custom measure {.val {measure_name}} must be a numeric vector."
    )
  }
  if (!is.null(n_obs) && length(res$pointwise) != n_obs) {
    cli::cli_abort(c(
      "{.field pointwise} from custom measure {.val {measure_name}} must have",
      "length {.val {n_obs}}, not {.val {length(res$pointwise)}}."
    ))
  }
  # pass measure name if user set it as attribute
  attr(res, "measure") <- measure_name
  if (!is.null(res$extra) && !is.list(res$extra)) {
    cli::cli_abort(c(
      "{.field extra} from custom measure {.val {measure_name}} must be a list.",
      "i" = "It is handed to {.code custom_se_fn(ref, cmp)} as
             {.code ref$extra} and {.code cmp$extra}."
    ))
  }

  invisible(res)
}

#' Validate a numeric matrix argument
#'
#' @description
#' Checks that `x` is a numeric matrix with at least one row and column.
#' Optionally enforces expected `nrow` and/or `ncol`. Aborts via
#' [cli::cli_abort()] when validation fails.
#'
#' @param x Object to validate.
#' @param arg Name of the argument (used in error messages).
#' @param nrow Expected number of rows, or `NULL` to skip this check.
#' @param ncol Expected number of columns, or `NULL` to skip this check.
#'
#' @return `NULL`, invisibly, on success.
#'
#' @noRd
.validate_numeric_matrix <- function(x, arg, nrow = NULL, ncol = NULL) {
  if (!is.numeric(x) || (!is.matrix(x) && !is.array(x))) {
    cli::cli_abort(
      "{.arg {arg}} must be a numeric matrix or array, not {.obj_type_friendly {x}}."
    )
  }
  if (!is.null(nrow) && nrow(x) != nrow) {
    cli::cli_abort(
      "{.arg {arg}} must have {.val {nrow}} row{?s}, not {.val {nrow(x)}}."
    )
  }
  if (!is.null(ncol) && ncol(x) != ncol) {
    cli::cli_abort(
      "{.arg {arg}} must have {.val {ncol}} column{?s}, not {.val {ncol(x)}}."
    )
  }
  if (nrow(x) < 1 || ncol(x) < 1) {
    cli::cli_abort("{.arg {arg}} must have at least 1 row and 1 column.")
  }
}

#' Validate a numeric vector argument
#'
#' @description
#' Checks that `x` is a numeric atomic vector (not a matrix or array) with
#' length at least one. Optionally enforces an expected `len`. Aborts via
#' [cli::cli_abort()] when validation fails.
#'
#' @param x Object to validate.
#' @param arg Name of the argument (used in error messages).
#' @param len Expected length, or `NULL` to skip this check.
#'
#' @return `NULL`, invisibly, on success.
#'
#' @noRd
.validate_numeric_vector <- function(x, arg, len = NULL) {
  if (!is.atomic(x) || !is.numeric(x) || is.matrix(x) || is.array(x)) {
    cli::cli_abort("{.arg {arg}} must be a numeric vector.")
  }
  if (!is.null(len) && length(x) != len) {
    cli::cli_abort(
      "{.arg {arg}} must have length {.val {len}}, not {.val {length(x)}}."
    )
  }
  if (length(x) < 1) {
    cli::cli_abort("{.arg {arg}} must not be empty.")
  }
}

#' Validate and normalize log weights
#'
#' @description
#' Validates that `log_weights` is a numeric matrix of size `n_draws` by
#' `n_obs`, then column-normalizes it via `.normalize_log_weights()`.
#'
#' @param log_weights Numeric matrix of log weights (`n_draws` \eqn{\times}
#'   `n_obs`).
#' @param n_draws Expected number of rows (posterior draws).
#' @param n_obs Expected number of columns (observations).
#'
#' @return Numeric matrix of the same dimensions as `log_weights` with
#'   column-normalized log weights.
#'
#' @noRd
.normalize_and_validate_log_weights <- function(log_weights, n_draws, n_obs) {
  .validate_numeric_matrix(
    log_weights,
    arg = "log_weights",
    nrow = n_draws,
    ncol = n_obs
  )
  .normalize_log_weights(log_weights)
}

#' Inform about ignored inputs when pointwise is supplied
#'
#' @description
#' When `pointwise` is not `NULL`, emits an informative message via
#' [cli::cli_inform()] listing non-`NULL` entries in `ignored_args` that are
#' not used.
#'
#' @param pointwise Optional precomputed pointwise contributions. When
#'   `NULL`, no message is emitted.
#' @param ignored_args Named list of arguments that may be ignored (e.g.
#'   `ylp`, `log_weights`).
#' @param fun_name Name of the calling function (shown in the message).
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.inform_ignored_inputs <- function(pointwise, ignored_args, fun_name) {
  if (is.null(pointwise)) {
    return(invisible(NULL))
  }
  supplied <- names(ignored_args)[vapply(ignored_args, Negate(is.null), logical(1))]
  if (length(supplied) > 0L) {
    cli::cli_inform(
      "In {.fn {fun_name}}, {.arg pointwise} is provided; ignoring {.arg {supplied}}."
    )
  }
  invisible(NULL)
}

#' Validate probability values
#'
#' @description
#' Checks that all elements of `x` lie in the closed interval \eqn{[0, 1]}.
#' Aborts via [cli::cli_abort()] when any value is out of range.
#'
#' @param x Numeric vector or matrix of probabilities.
#' @param arg Name of the argument (used in error messages).
#'
#' @return `NULL`, invisibly, on success.
#'
#' @noRd
.validate_probs <- function(x, arg) {
  if (!all(x >= 0 & x <= 1)) {
    cli::cli_abort("{.arg {arg}} must contain values in [0, 1].")
  }
}

#' Pointwise log predictive density from measure inputs
#'
#' @description
#' `measure_elpd()`, `measure_mlpd()` and `measure_ic()` take the same inputs:
#' precomputed `pointwise` values, or an `ylp` matrix with optional
#' `log_weights`. This helper holds their shared validation.
#'
#' @param ylp A draws x observations matrix of log predictive densities, or a
#'   3-D array.
#' @param log_weights Optional log weights, normalized before use.
#' @param pointwise Optional numeric vector of precomputed pointwise values.
#' @param fun_name Name of the calling measure. Used in the message that
#'   reports ignored inputs.
#'
#' @return A list with `lppd_i`, `n_draws` (`NULL` when `pointwise` is
#'   supplied) and `n_obs`.
#'
#' @noRd
.lppd_from_inputs <- function(ylp, log_weights, pointwise, fun_name) {
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = fun_name
    )
    return(list(lppd_i = pointwise, n_draws = NULL, n_obs = length(pointwise)))
  }

  .validate_numeric_matrix(ylp, arg = "ylp")
  ylp <- if (is.array(ylp) && length(dim(ylp)) == 3) {
    llarray_to_matrix(ylp)
  } else {
    ylp
  }
  n_draws <- nrow(ylp)
  n_obs <- ncol(ylp)
  if (!is.null(log_weights)) {
    log_weights <- .normalize_and_validate_log_weights(
      log_weights = log_weights, n_draws = n_draws, n_obs = n_obs
    )
  }
  list(
    lppd_i = ptw_log_pred_density(ylp, log_weights),
    n_draws = n_draws,
    n_obs = n_obs
  )
}

#' Weighted pointwise classification accuracy
#'
#' @description
#' Shared by `measure_acc()` and `measure_bacc()`. Maps each observation to a
#' predicted class, then compares it with `y`. Binary `mupred` is thresholded
#' at 0.5. A 3-D `mupred` takes the argmax over categories.
#'
#' @param y An integer vector of observed class labels.
#' @param mupred A draws x observations matrix, or a draws x observations x
#'   categories array, of predicted probabilities.
#' @param log_weights Optional log weights. Draws are equally weighted when
#'   `NULL`.
#'
#' @return An integer vector of 0/1 accuracy contributions.
#'
#' @noRd
.acc_pointwise <- function(y, mupred, log_weights) {
  if (!is.numeric(mupred) || (length(dim(mupred)) != 2 && length(dim(mupred)) != 3)) {
    cli::cli_abort(
      "{.arg mupred} must be a numeric matrix or 3D numeric array."
    )
  }
  .validate_probs(mupred, arg = "mupred")

  if (!is.null(log_weights)) {
    weights <- exp(.normalize_and_validate_log_weights(
      log_weights = log_weights,
      n_draws = nrow(mupred),
      n_obs = dim(mupred)[2]
    ))
  } else {
    weights <- rep(1 / nrow(mupred), nrow(mupred))
  }

  if (length(dim(mupred)) == 3) {
    # Multiclass: (draws × obs × categories) > argmax over categories
    weighted_mupred <- apply(array(weights, dim(mupred)) * mupred, c(2, 3), sum)
    mupred_hat <- apply(weighted_mupred, 1, which.max)
  } else {
    .validate_numeric_matrix(mupred, arg = "mupred")
    weighted_mupred <- colSums(mupred * weights)
    mupred_hat <- (weighted_mupred > 0.5) * 1L
  }

  (mupred_hat == y) * 1L
}

#' Pointwise prediction error from measure inputs
#'
#' @description
#' Shared by `measure_mae()` and `measure_mse()`. Forms a point prediction for
#' each observation, then applies `transform` to the residual. The point
#' prediction is the mean of the `mupred` draws, or their weighted mean when
#' `log_weights` is supplied.
#'
#' @param y A numeric vector of observed outcomes.
#' @param mupred A draws x observations matrix of point predictions. A vector
#'   is coerced to a 1 x n matrix.
#' @param log_weights Optional log weights, normalized before use.
#' @param pointwise Optional numeric vector of precomputed pointwise errors.
#' @param fun_name Name of the calling measure. Used in the messages that
#'   report ignored inputs and the coercion of `mupred`.
#' @param transform Function applied to the residual `y - yhat`.
#'
#' @return A list with `err_i`, `n_draws` (`NULL` when `pointwise` is
#'   supplied) and `n_obs`.
#'
#' @noRd
.point_error_from_inputs <- function(y, mupred, log_weights, pointwise,
                                     fun_name, transform) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = fun_name
    )
    return(list(err_i = pointwise, n_draws = NULL, n_obs = length(pointwise)))
  }

  n_draws <- nrow(mupred)
  n_obs <- ncol(mupred)
  .validate_numeric_vector(y, arg = "y")
  if (!is.null(mupred) && !is.matrix(mupred)) {
    .validate_numeric_vector(mupred, arg = "mupred", len = length(y))
    cli::cli_inform(
      "Coercing {.arg mupred} from vector to 1 x n matrix for {.fn {fun_name}}."
    )
    mupred <- matrix(mupred, nrow = 1, ncol = length(mupred))
  }
  .validate_numeric_matrix(mupred, arg = "mupred", ncol = length(y))
  if (is.null(log_weights)) {
    yhat <- colMeans(mupred)
  } else {
    weights <- exp(.normalize_and_validate_log_weights(
      log_weights = log_weights, n_draws = n_draws, n_obs = n_obs
    ))
    yhat <- colSums(weights * mupred)
  }
  list(err_i = transform(y - yhat), n_draws = n_draws, n_obs = n_obs)
}

#' Copy selected attributes between objects
#'
#' @description
#' Copies attributes named in `which` from `from` onto `to`, overwriting any
#' existing attributes with the same names.
#'
#' @param to Object receiving attributes.
#' @param from Object supplying attributes.
#' @param which Character vector of attribute names to copy.
#'
#' @return `to`, with updated attributes.
#'
#' @noRd
.copy_attrs <- function(to, from, which) {
  for (nm in which) {
    attr(to, nm) <- attr(from, nm)
  }
  to
}

#' Normalize log weights
#'
#' @description
#' Normalizes a matrix of log weights column-wise so that the weights in each
#' column sum to one on the probability scale. Normalization is performed by
#' subtracting the log-sum-exp of each column from its elements, equivalent to
#' dividing each column's weights by their sum on the probability scale.
#'
#' @param log_weights Numeric matrix of log weights (`n_draws` \eqn{\times}
#'   `n_obs`), where rows are draws and columns are observations.
#'
#' @return Numeric matrix of the same dimensions as `log_weights` with
#'   column-normalized log weights.
#'
#' @noRd
.normalize_log_weights <- function(log_weights) {
  sweep(
    log_weights,
    2,
    matrixStats::colLogSumExps(log_weights),
    FUN = "-",
    check.margin = FALSE
  )
}


#' Probability-weighted moment estimator of E|X - X'|
#'
#' @description
#' Estimates \eqn{E[|X - X'|]}, the expected absolute difference between two
#' independent draws from the predictive distribution, which is the term the
#' (S)RPS and (S)CRPS scores in [measure_rps()] are built from. Written as a
#' weighted U-statistic over draw pairs, the estimator is
#' \deqn{E[|X - X'|] = \frac{\sum_i \sum_{j \neq i} w_i w_j |x_i - x_j|}{1 -
#'   \sum_i w_i^2} = \frac{2 \sum_s w_{(s)} x_{(s)} (C_s + C_{s-1} - 1)}{1 -
#'   \sum_s w_{(s)}^2},}
#' where \eqn{x_{(s)}} are the draws sorted in ascending order, \eqn{w_{(s)}}
#' their weights, and \eqn{C_s = \sum_{k \le s} w_{(k)}}. The second form is the
#' one computed here and needs only a single sort and cumulative sum per column.
#'
#' With equal weights \eqn{w_{(s)} = 1/S} this reduces exactly to the classic
#' unbiased pairwise estimator with the \eqn{1 / (S (S - 1))} normalization,
#' i.e. `colMeans(x_sorted * 2 * (2 * (1:S) - S - 1) / (S - 1))`, which is the
#' probability-weighted moment estimator of Taillardat et al. (2016). Because
#' the estimate is a convex combination of \eqn{|x_i - x_j|} it is always
#' non-negative, so `log()` of it in the scaled scores is always defined, and
#' the coefficients sum to zero, so it is invariant to shifts of `ypred`.
#'
#' This is the bias-corrected weighted Gini mean difference, not the estimator
#' derived in `notes/crps_pwm.pdf`; see decision D5 in `notes/developer-notes.md`
#' for why that derivation is not used here.
#'
#' @param ypred Numeric matrix of posterior predictive draws (`n_draws`
#'   \eqn{\times} `n_obs`), where rows are draws and columns are observations.
#' @param w Optional numeric matrix of column-normalized weights on the
#'   probability scale, of the same dimensions as `ypred`. `NULL` (the default)
#'   uses equal weights and takes the faster unweighted path.
#'
#' @return Numeric vector of length `ncol(ypred)` with one estimate of
#'   \eqn{E[|X - X'|]} per observation.
#'
#' @noRd
.exx_pwm <- function(ypred, w = NULL) {
  n_draws <- nrow(ypred)
  if (n_draws < 2) {
    stop(
      "`ypred` must have at least 2 draws (rows) to estimate E|X - X'|, ",
      "which the RPS/CRPS measures require.",
      call. = FALSE
    )
  }

  if (is.null(w)) {
    ypred_sorted <- apply(ypred, 2, sort)
    coefs <- 2 * (2 * seq_len(n_draws) - n_draws - 1) / (n_draws - 1)
    return(colMeans(ypred_sorted * coefs))
  }

  vapply(
    seq_len(ncol(ypred)),
    function(j) {
      ord <- order(ypred[, j])
      x_sorted <- ypred[ord, j]
      w_sorted <- w[ord, j]
      denominator <- 1 - sum(w_sorted^2)
      # All the weight sits on a single draw: the pair (X, X') is degenerate
      # and E|X - X'| is 0.
      if (denominator <= 0) {
        return(0)
      }
      # C_s + C_{s-1} = 2 * C_s - w_(s)
      coefs <- 2 * cumsum(w_sorted) - w_sorted - 1
      2 * sum(w_sorted * x_sorted * coefs) / denominator
    },
    numeric(1)
  )
}

#' Validate control argument
#'
#' @description
#' Validates that the arguments passed to the control list are valid
#' arguments for the corresponding function. If not, a warning is issued that
#' corresponding invalid argument is ignored.
#'
#' @param control Named list of per-measure settings.
#' @param measures Optional list of normalized measure entries from
#'   `.prepare_measures()`. When supplied, control names are resolved against
#'   the requested measures, so custom measures are validated against their own
#'   formals; without it only built-in names can be checked.
#'
#' @keywords internal
#' @noRd
.validate_control <- function(control, measures = NULL) {
  res <- checkmate::check_list(control, types = "list", names = "named")
  if (!isTRUE(res)) {
    cli::cli_abort(c(
      "{.arg control} must be a named list of named lists.",
      "i" = "Expected format: {.code list(fun_name = list(arg1 = val1, arg2 = val2))}"
    ))
  }

  # without `measures` the requested measures are unknown, so a control name is
  # only checked against the built-in registry
  known_measures <- !is.null(measures)
  if (is.null(measures)) {
    measures <- list()
  }
  entries <- stats::setNames(
    measures,
    vapply(measures, function(e) e$name, character(1L))
  )

  for (func_name in names(control)) {
    entry <- if (func_name %in% names(entries)) entries[[func_name]] else NULL
    # custom measures are validated against their own formals, built-ins
    # against the registry; a name matching neither accepts nothing
    valid_args <- if (!is.null(entry) && identical(entry$type, "custom")) {
      names(formals(entry$key))
    } else if (is.null(entry) && known_measures) {
      NULL
    } else {
      spec <- .measure_spec[[if (is.null(entry)) func_name else entry$key]]
      if (is.null(spec)) NULL else names(formals(spec$fun))
    }
    if (is.null(valid_args)) {
      cli::cli_warn(c(
        "Ignoring {.arg control} entry {.val {func_name}}, which matches no",
        "measure being computed."
      ))
      next
    }
    invalid_args <- setdiff(names(control[[func_name]]), valid_args)
    if (length(invalid_args) > 0) {
      cli::cli_warn(
        "Ignoring {.arg {invalid_args}} as it is not a valid argument of {.fn {func_name}}."
      )
    }
  }
  invisible(NULL)
}

#' Subset measure results
#' 
#' @description
#' Subsets the measure results to the specified measures and components.
#' 
#' @param x Measure results object.
#' @param measures Character vector of measures to subset.
#' @param components Character vector of components to subset.
#' 
#' @return Subsetted measure results object.
#' 
#' @noRd
subset_measures <- function(x, measures, components) {
  invalid_components <- setdiff(components, names(x))
  if (length(invalid_components) > 0) {
    cli::cli_abort(c(
      "{.arg components} contains invalid value{?s}: {.val {invalid_components}}.",
      "i" = "Valid components: {.val {names(x)}}."
    ))
  }
  components <- intersect(components, names(x))
  
  available_measures <- if ("estimates" %in% components) {
    rownames(x$estimates)
  } else if ("pointwise" %in% components) {
    colnames(x$pointwise)
  }

  if (!is.null(available_measures)) {
    invalid_measures <- setdiff(measures, available_measures)
    if (length(invalid_measures) > 0) {
      cli::cli_abort(c(
        "{.arg measures} contains invalid value{?s}: {.val {invalid_measures}}.",
        "i" = "Valid measures: {.val {available_measures}}."
      ))
    }
  }
  
  result <- x[components]

  if ("estimates" %in% components) {
    rows <- intersect(measures, rownames(result$estimates))
    result$estimates <- result$estimates[rows, , drop = FALSE]
  }

  if ("pointwise" %in% components) {
    cols <- intersect(measures, colnames(result$pointwise))
    result$pointwise <- result$pointwise[, cols, drop = FALSE]
  }

  result
}
