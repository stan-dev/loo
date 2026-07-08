#' @section `measure_higher_is_better` attribute:
#' All `*_pred_measure()` and [pred_measure()] results include attribute
#' `measure_higher_is_better`: a named list recording the `higher_is_better`
#' setting used for each measure (`NULL`, `TRUE`, or `FALSE` per bare measure
#' name; `elpd` is always `NULL`). Built-in loss measures such as MSE are stored
#' on a loss scale by default (`NULL`); pass
#' `control = list(mse = list(higher_is_better = TRUE))` to store values on a
#' utility scale. When measures are added incrementally with [pred_measure()],
#' the attribute is updated for newly computed measures.
#'
#' Attribute `measure_compare_meta` records per-measure comparison metadata
#' (`higher_is_better`, `loss`, and `diff_method`) used by [loo_compare()].
#' Built-in measures take `loss` and `diff_method` from the package measure
#' registry; custom measures default to `diff_method = "auto"` and `loss = FALSE`
#' (utility scale). For custom loss measures, use `higher_is_better = TRUE` in
#' `control` so [loo_compare()] interprets the sign correctly. See [loo-glossary].
