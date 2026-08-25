#' @section `measure_higher_is_better` attribute:
#' All `*_pred_measure()` and [pred_measure()] results include attribute
#' `measure_higher_is_better`: a named list recording the `higher_is_better`
#' setting used for each measure (`TRUE` or `FALSE` per bare measure name).
#' Measures left at their natural orientation have no entry, which reads as
#' `NULL`; `elpd` is always `NULL`. Built-in loss measures such as MSE are
#' stored on a loss scale by default; pass
#' `control = list(mse = list(higher_is_better = TRUE))` to store values on a
#' utility scale. The same works for a custom measure, using the name it was
#' given in `measure`. When measures are added incrementally with
#' [pred_measure()], the attribute is updated for newly computed measures.
#'
#' Attribute `measure_compare_meta` records per-measure comparison metadata
#' (`higher_is_better`, `loss`, and `diff_method`) used by [model_compare()].
#' Built-in measures take `loss` and `diff_method` from the package measure
#' registry. Custom measures always get `diff_method = "custom"` and take the
#' standard error of their difference from the `custom_se_fn` argument of
#' [model_compare()]; their `loss` comes from
#' `attr(my_fun, "measure_loss") <- TRUE`, which declares that lower values are
#' better. Without that declaration a custom measure is treated as a utility, so
#' an undeclared loss is compared and ranked in the wrong direction. See
#' [loo-glossary].
