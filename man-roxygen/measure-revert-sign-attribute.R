#' @section `measure_revert_sign` attribute:
#' All `*_pred_measure()` and [pred_measure()] results include attribute
#' `measure_revert_sign`: a named list recording whether `revert_sign` was
#' applied to each measure when it was computed (logical per bare measure name;
#' `elpd` is always `FALSE`). Built-in loss measures such as MSE are stored on
#' a loss scale by default (`FALSE`); pass
#' `control = list(mse = list(revert_sign = TRUE))` to store negated values.
#' When measures are added incrementally with [pred_measure()], the attribute
#' is updated for newly computed measures. [loo_compare()] uses
#' `measure_revert_sign` when converting paired differences to a common utility
#' scale; see [loo-glossary].
