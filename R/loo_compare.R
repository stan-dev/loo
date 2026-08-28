#' Model comparison (previous name)
#'
#' @description
#' `loo_compare()` is the previous name of [model_compare()] and is kept as a
#' working alias. It is a generic so that packages registering their own
#' `loo_compare` methods keep dispatching as before; the default method simply
#' forwards to [model_compare()].
#'
#' New code should call [model_compare()], which additionally compares
#' [`kfold_pred_measure()`][kfold_pred_measure],
#' [`test_pred_measure()`][test_pred_measure], and
#' [`insample_pred_measure()`][insample_pred_measure] results.
#'
#' @export
#' @inheritParams model_compare
#' @return See [model_compare()].
#'
#' @seealso [model_compare()]
#'
#' @examples
#' LL <- example_loglik_array()
#' loo1 <- loo(LL)
#' loo2 <- loo(LL + 1)
#'
#' # identical to model_compare(loo1, loo2)
#' loo_compare(loo1, loo2)
#'
loo_compare <- function(x, ..., rank_by = NULL, custom_se_fn) {
  UseMethod("loo_compare")
}

#' @rdname loo_compare
#' @export
loo_compare.default <- function(x, ..., rank_by = NULL, custom_se_fn) {
  # Forward `custom_se_fn` only when it was actually supplied, so that
  # `model_compare()` can still tell an omitted argument from an explicit NULL.
  if (missing(custom_se_fn)) {
    model_compare(x, ..., rank_by = rank_by)
  } else {
    model_compare(x, ..., rank_by = rank_by, custom_se_fn = custom_se_fn)
  }
}

#' @rdname loo_compare
#' @export
loo_compare.psis_loo_ss_list <- function(x, ...) {
  model_compare.psis_loo_ss_list(x, ...)
}
