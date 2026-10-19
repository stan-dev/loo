#' @section `measure_info` attribute:
#' All `*_pred_measure()` and [pred_measure()] results include attribute
#' `measure_info`: a named list, keyed by bare measure name, recording what
#' [model_compare()] needs to know about each measure — `loss` (whether lower
#' values are better), `diff_method` (how the standard error of a difference is
#' obtained), and, where applicable, `se_diff_fun` and `extra`. When measures
#' are added incrementally with [pred_measure()], the attribute is extended for
#' the newly computed measures.
#'
#' Built-in measures take `loss`, `diff_method`, and `se_diff_fun` from the
#' package measure registry. Custom measures always get `diff_method = "custom"`
#' and take `loss` and `se_diff_fun` from [custom_measure()]. Without
#' `loss = TRUE` a custom measure is treated as a utility, so an undeclared
#' loss is compared and ranked in the wrong direction. See [loo-glossary].
