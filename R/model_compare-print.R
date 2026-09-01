#' @rdname model_compare
#' @export
#' @param digits For the print method only, the number of digits to use when
#'   printing.
#' @param p_worse For the print method only, should we include the normal
#'   approximation based probability of each model having worse performance than
#'   the reference model? The default is `TRUE`.
#' @param simplify For the print method only, should the output be simplified to
#'   only include the model names, differences, and (when `p_worse = TRUE`)
#'   diagnostic columns? The default is `TRUE`. Set to `FALSE` to also print the
#'   available estimate columns: pointwise ELPD, LOOIC/WAIC and their standard
#'   errors for classic comparisons. For [`pred_measure`][pred_measure]
#'   comparisons each printed table gains its own measure's estimate and
#'   standard error, and an ELPD table also gains `p` and `se_p`.
#' @param measures For `loo_pred_measure` comparisons only, which measures to
#'   print diff tables for. `NULL` (default) prints only the ranking measure
#'   (`"elpd"` when `rank_by` was not set, otherwise `rank_by`);
#'   `"all"` prints all compared measures; or a character vector of measure
#'   names (e.g. `c("elpd", "mse")`). Each table is sorted by its own measure,
#'   best model first, so the same model need not lead every table.
print.compare.loo <- function(x, ..., digits = 1, p_worse = TRUE,
                              simplify = TRUE, measures = NULL) {
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
      simplify = simplify,
      measures = measures
    ))
  }

  if (!all(c("model", "elpd_diff", "se_diff") %in% colnames(x))) {
    print(as.data.frame(x))
    return(x)
  }
  base_cols <- c("model", "elpd_diff", "se_diff", "subsampling_se_diff")
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

  rank_spec <- attr(x, "rank_by")
  if (identical(rank_spec$kind, "model")) {
    message("Differences computed against model ", rank_spec$model, ".")
  }
  .print_compare_diag_message(x, p_worse = p_worse)
  invisible(x)
}

#' Print `compare.loo` results from `pred_measure` comparisons
#' @noRd
.print_compare_pred_measure <- function(x, digits, p_worse, simplify,
                                       measures) {
  rank_spec <- attr(x, "rank_by")
  compare_measures <- attr(x, "compare_measures")
  compare_source <- attr(x, "compare_source")
  primary_measure <- if (is.null(rank_spec)) "elpd" else rank_spec$measure

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
    .compare_reference_line(x, rank_spec, compare_measures)
  )

  # LOO is the familiar default, so only name the source when it is not LOO.
  # Spelled out for a full sentence, unlike the short per-object tag from
  # `.pred_measure_source_label()`.
  if (!is.null(compare_source) && !identical(compare_source, "loo")) {
    .cat_wrapped(
      "Predictive measures evaluated on ",
      switch(
        compare_source,
        kfold = "K-fold cross-validation",
        test = "held-out test data",
        insample = "in-sample (training) data",
        compare_source
      ),
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
        "\n-- ", measure, " (vs ", .measure_ref_model(x, measure), ") --\n",
        sep = ""
      )
    }
    .print_compare_measure_table(
      x,
      measure = measure,
      digits = digits,
      p_worse = p_worse,
      simplify = simplify
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
        "Use print(x, measures = \"all\") to see all measures."
      )
    }
  }

  invisible(x)
}

#' Header line naming the reference each difference is computed against
#'
#' Always printed, so the reference is never left implicit. `rank_by` pins one
#' reference for every measure; without it each measure keeps its own best
#' model, which is the case most in need of being spelled out.
#' @noRd
#' @param x A `"compare.loo"` data frame.
#' @param rank_spec Attribute `rank_by`, or `NULL` for an object created before
#'   it was set.
#' @param compare_measures Bare names of all compared measures.
#' @return A single string.
.compare_reference_line <- function(x, rank_spec, compare_measures) {
  if (identical(rank_spec$kind, "measure")) {
    return(paste0(
      "Models ranked by ", rank_spec$measure,
      " (reference: ", .measure_ref_model(x, rank_spec$measure), ")."
    ))
  }
  if (identical(rank_spec$kind, "model")) {
    return(paste0("All measures compared against model ", rank_spec$model, "."))
  }

  refs <- vapply(compare_measures, .measure_ref_model, character(1), x = x)
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

#' Print prose wrapped to the conventional 80-column terminal width
#'
#' Tables are wrapped by `print.data.frame()` at `getOption("width")`; this does
#' the same for the sentences around them, capped at 80 so the output stays
#' within a standard terminal however wide the option is set.
#' @noRd
#' @param ... Pieces of a single line, pasted together.
.cat_wrapped <- function(...) {
  width <- min(getOption("width", 80L), 80L)
  cat(paste(strwrap(paste0(...), width = width), collapse = "\n"), "\n", sep = "")
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
.print_compare_measure_table <- function(x, measure, digits, p_worse,
                                        simplify = TRUE) {
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

  # The frame carries every measure's per-model estimate and SE under its bare
  # name (`model_compare_matrix(bare_names = TRUE)`). `simplify = FALSE` adds
  # the pair this table's own measure owns, after the diagnostic columns, as
  # the classic path does. `p` is an ELPD companion, so it rides with it.
  if (!simplify) {
    est_cols <- c(measure, paste0("se_", measure))
    if (.is_elpd_measure(measure)) {
      est_cols <- c(est_cols, "p", "se_p")
    }
    for (col in intersect(est_cols, colnames(x))) {
      x2[[col]] <- unname(.fr(x[[col]][ord], digits))
    }
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
