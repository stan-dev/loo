#' Normalize the `measure` argument to an internal list
#'
#' @description
#' Converts `measure` (character, function, list, or `NULL`) into a list of
#' entries with elements `name`, `type` (`"builtin"` or `"custom"`), and `key`
#' (built-in name or function).
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
        list(name = nm, type = "custom", key = el)
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
#' @param fun Function implementing a custom measure.
#' @noRd
.measure_entry_custom <- function(fun) {
  name <- attr(fun, "measure_name", exact = TRUE)
  if (is.null(name) || length(name) != 1L || !nzchar(name)) {
    cli::cli_abort(c(
      "A custom function passed to {.arg measure} must have attribute",
      "{.code measure_name}.",
      "i" = "Set {.code attr(my_fun, \"measure_name\") <- \"my_metric\"}."
    ))
  }
  list(name = name, type = "custom", key = fun)
}

#' Check duplicate and reserved measure names
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
#'
#' @return A list of normalized measure entries ready for computation.
#'
#' @noRd
.prepare_measures <- function(measure, predperf, supported_measures_list) {
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
    entry_names <- vapply(entries, function(e) e$name, character(1L))
    dups <- intersect(entry_names, existing_measures)
    if (length(dups) > 0L) {
      cli::cli_warn(c(
        "!" = "{cli::qty(length(dups))} Measure{?s} {.val {dups}} {?is/are}",
        "already present in {.arg predperf} and will be skipped."
      ))
    }
    keep <- !vapply(
      entries,
      function(e) e$name %in% existing_measures,
      logical(1L)
    )
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
  if (!is.matrix(x) && !is.array(x)) {
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
  if (!all(x >= 0 | x <= 1)) {
    cli::cli_abort("{.arg {arg}} must contain values in [0, 1].")
  }
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


#' Validate control argument
#' 
#' @description
#' Validates that the arguments passed to the control list are valid 
#' arguments for the corresponding function. If not, a warning is issued that 
#' corresponding invalid argument is ignored.
#' 
#' @param control Named list of per-measure settings.
#' 
.validate_control <- function(control) {
  res <- checkmate::check_list(control, types = "list", names = "named")
  if (!isTRUE(res)) {
    cli::cli_abort(c(
      "{.arg control} must be a named list of named lists.",
      "i" = "Expected format: {.code list(fun_name = list(arg1 = val1, arg2 = val2))}"
    ))
  }
  
  for (func_name in names(control)) {
    invalid_args <- names(control[[func_name]])[
      !names(control[[func_name]]) %in% names(formals(match.fun(func_name)))
    ]
    if (length(invalid_args) > 0) {
      cli::cli_warn(
        "Ignoring {.arg {invalid_args}} as it is not a valid argument of {.fn {func_name}}."
      )
    }
  }
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

  if ("diagnostics" %in% components) {
    result$diagnostics <- result$diagnostics
  }

  result
}