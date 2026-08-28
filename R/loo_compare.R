#' Model comparison (deprecated)
#'
#' @description
#' **This function is deprecated**. Please use the new [model_compare()] function
#' instead. See `vignette("migration-guide", package = "loo")` for details.
#'
#' `loo_compare()` remains an exported generic so that `loo_compare` methods
#' registered by other packages keep dispatching as before, but it is frozen at
#' its previous behavior: it compares `"loo"`, `"waic"`, and `"kfold"` objects on
#' [ELPD][loo-glossary] only. Comparing
#' [`pred_measure`][pred_measure] results, or using the `rank_by` and
#' `custom_se_fn` arguments, requires [model_compare()].
#'
#' @export
#' @param x An object of class `"loo"` or a list of such objects. If a list is
#'   used then the list names will be used as the model names in the output.
#' @param ... Additional objects of class `"loo"`, if not passed in as a single
#'   list.
#' @return See [model_compare()]. For the inputs `loo_compare()` still accepts,
#'   the result is identical to what [model_compare()] returns.
#'
#' @seealso [model_compare()]
#'
#' @examples
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)
#' loo2 <- loo(LL + 1)
#'
#' # deprecated; identical to model_compare(loo1, loo2)
#' suppressWarnings(loo_compare(loo1, loo2))
#'
#' # use this instead
#' model_compare(loo1, loo2)
#'
loo_compare <- function(x, ...) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ...) {
  # `old` is set explicitly: without it `.Deprecated()` reports the method name
  # (e.g. "loo_compare.default") rather than the function users called.
  .Deprecated("model_compare", old = "loo_compare")

  # `loo_compare()` keeps its old signature, so the arguments added to
  # `model_compare()` would arrive through `...` and be mistaken for models.
  dots <- list(...)
  new_args <- intersect(names(dots), c("rank_by", "custom_se_fn"))
  if (length(new_args)) {
    stop(
      "`", new_args[1L], "` is not supported by the deprecated `loo_compare()`. ",
      "Use `model_compare()` instead.",
      call. = FALSE
    )
  }

  loos <- .model_compare_inputs(x, ...)
  if (any(vapply(loos, is.pred_measure, logical(1)))) {
    stop(
      "`loo_compare()` compares only 'loo', 'waic', and 'kfold' objects. ",
      "Use `model_compare()` to compare 'pred_measure' results.",
      call. = FALSE
    )
  }

  model_compare(loos)
}

#' @rdname loo_compare
#' @export
loo_compare.psis_loo_ss_list <- function(x, ...) {
  # `old` is set explicitly: without it `.Deprecated()` reports the method name
  # (e.g. "loo_compare.default") rather than the function users called.
  .Deprecated("model_compare", old = "loo_compare")
  model_compare.psis_loo_ss_list(x, ...)
}
