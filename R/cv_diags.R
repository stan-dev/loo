#' Build a validated `cv_diags` object from the per-fold results
#'
#' @noRd
#' @param pointwise Data frame with `elpd`, `mcse_elpd`, `bias_elpd`.
#' @param mcse_total Total MCSE of the method.
#' @param diagnostics Method-specific list.
#' @param method One of `"psis"`, `"mixis"`, `"brute_force"`, `"cv_bridge"`.
#' @param S Number of draws.
#'
.cv_diags <- function(pointwise, mcse_total, diagnostics, method, S) {
  rownames(pointwise) <- NULL
  mc_error <- .mc_error_summary(mcse_total, sum(pointwise$bias_elpd))
  .validate_cv_diags(.new_cv_diags(pointwise, mc_error, diagnostics, method,
                                   c(S, nrow(pointwise))))
}

#' Low-level constructor: no checks
#' @noRd
.new_cv_diags <- function(pointwise, mc_error, diagnostics, method, dims) {
  structure(
    list(pointwise = pointwise, mc_error = mc_error, diagnostics = diagnostics,
         method = method, dims = as.integer(dims)),
    class = "cv_diags"
  )
}

#' Check the invariants of a `cv_diags` object
#' @noRd
.validate_cv_diags <- function(x) {
  required <- c("elpd", "mcse_elpd", "bias_elpd")
  missing <- setdiff(required, names(x$pointwise))
  if (length(missing) > 0) {
    stop("Internal error: 'pointwise' lacks ", toString(missing), ".", call. = FALSE)
  }
  checkmate::assert_choice(x$method, c("psis", "mixis", "brute_force", "cv_bridge"))
  if (nrow(x$pointwise) != x$dims[2]) {
    stop("Internal error: 'pointwise' must have one row per fold.", call. = FALSE)
  }
  x
}

#' Totals of the Monte Carlo error
#'
#' The interval is for the exact elpd, as offsets from the estimate. A
#' positive bias widens it downward.
#'
#' @noRd
#'
.mc_error_summary <- function(mcse_total, bias_total) {
  list(
    mcse_total = mcse_total,
    bias_total = bias_total,
    rmse_pred = sqrt(bias_total^2 + mcse_total^2),
    interval = c(lower = -1.96 * mcse_total - max(0, bias_total),
                 upper = 1.96 * mcse_total - min(0, bias_total))
  )
}

#' @export
print.cv_diags <- function(x, digits = 3, n = 5, ...) {
  cat(sprintf("Monte Carlo diagnostics (%s): %d draws, %d folds\n",
              x$method, x$dims[1], x$dims[2]))
  e <- x$mc_error
  cat(sprintf("  elpd total  %10.*f\n", digits, sum(x$pointwise$elpd)))
  cat(sprintf("  mcse_total  %10.*f\n", digits, e$mcse_total))
  cat(sprintf("  bias_total  %10.*f\n", digits, e$bias_total))
  cat(sprintf("  rmse_pred   %10.*f\n", digits, e$rmse_pred))
  cat(sprintf("  interval    [%.*f, %.*f]\n", digits, e$interval[["lower"]],
              digits, e$interval[["upper"]]))
  cat(sprintf("First %d folds:\n", min(n, x$dims[2])))
  pointwise <- utils::head(x$pointwise, n)
  is_num <- vapply(pointwise, is.double, logical(1))
  pointwise[is_num] <- round(pointwise[is_num], digits)
  print(format(pointwise, nsmall = digits))
  invisible(x)
}

#' @export
as.data.frame.cv_diags <- function(x, ...) {
  x$pointwise
}

#' @noRd
.assert_log_lik <- function(x, name = deparse(substitute(x))) {
  checkmate::assert_matrix(x, mode = "numeric", min.rows = 1, min.cols = 1,
                           .var.name = name)
  validate_log_lik(x)
}

#' A valid `chain_id`: `NULL`, or chains 1 to C of equal length
#' @noRd
.assert_chain_id <- function(chain_id, S, name = deparse(substitute(chain_id))) {
  if (is.null(chain_id)) {
    return(invisible(NULL))
  }
  checkmate::assert_integerish(chain_id, len = S, any.missing = FALSE,
                               .var.name = name)
  counts <- table(chain_id)
  if (any(counts != counts[[1]])) {
    stop("All chains in '", name, "' must have the same number of draws.",
         call. = FALSE)
  }
  if (!setequal(chain_id, seq_along(counts))) {
    stop("'", name, "' must hold the values 1 to the number of chains.",
         call. = FALSE)
  }
  invisible(chain_id)
}

#' ESS of the mean of a draw sequence
#'
#' Without `chain_id`, the draws are independent (loo's `r_eff = 1`) and
#' get the number of draws. A constant finite sequence gets the number of
#' draws; its variance is 0 for any ESS. A failed sequence gets `NA`, so
#' that the MCSE is `NA` and not too small.
#'
#' @noRd
.n_eff_of_mean <- function(x, chain_id) {
  if (is.null(chain_id)) {
    return(length(x))
  }
  n_eff <- tryCatch(relative_eff(x, chain_id) * length(x),
                    error = function(e) NA_real_)
  if (is.finite(n_eff) && n_eff > 0) {
    n_eff
  } else if (all(is.finite(x)) && diff(range(x)) < .Machine$double.eps) {
    length(x)
  } else {
    NA_real_
  }
}

#' Variance of an estimate from the influence of each draw
#'
#' \eqn{\mathrm{Var} \approx (S / \mathrm{ESS}(\psi)) \sum_s \psi_s^2}.
#' The ESS accounts for the autocorrelation within the chains.
#'
#' @noRd
#' @param psi Influence of each draw on the estimate; sums to 0.
#' @param chain_id `NULL` or a valid `chain_id`.
#'
.ij_variance <- function(psi, chain_id) {
  sum(psi^2) * length(psi) / .n_eff_of_mean(psi, chain_id)
}
