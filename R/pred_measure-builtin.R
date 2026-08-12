#' Shared parameters for all measures
#'
#' @param log_weights Optional numeric matrix of unnormalized log-importance
#'   weights with dimensions \eqn{S \times n}. Weights are column-normalized
#'   before computing each per-observation contribution.
#' @param pointwise Optional numeric vector of precomputed per-observation
#'   contributions. When supplied, `ylp` and `log_weights` are ignored.
#' @param higher_is_better Logical or `NULL`; whether larger values indicate
#'   better predictive performance. `NULL` (default) keeps each measure's
#'   natural convention (e.g. MSE on a loss scale, ELPD on a utility scale).
#'   Set `TRUE` for a utility scale (higher is better) or `FALSE` for a loss
#'   scale (lower is better).
#'
#' @return An object of class `"measure"`: a list with:
#'   \describe{
#'     \item{`estimates`}{Named numeric vector with elements `Estimate`
#'       and `SE` (standard error).}
#'     \item{`pointwise`}{Numeric vector of length \eqn{n} with per-observation
#'       values.}
#'   }
#'   Attributes `measure` (i.e., measure name) and `dims` (draws \eqn{\times} 
#'   observations) are also set. Use [print()] for a readable summary.
#' 
#' @keywords internal
#' @name measure_params
NULL

#' Shared parameters for density scores
#'
#' @param ylp A numeric matrix or three-dimensional array of log predictive
#'   densities or probabilities:
#'   \itemize{
#'     \item **Matrix** (\eqn{S \times n}): \eqn{S} posterior draws (chains
#'       merged) by \eqn{n} observations.
#'     \item **Array** (\eqn{I \times C \times n}): \eqn{I} MCMC iterations per
#'       chain, \eqn{C} chains, and \eqn{n} observations. Converted to an
#'       \eqn{S \times n} matrix internally.
#'   }
#' 
#' @keywords internal
#' @name measure_density_params
NULL

#' Shared parameters for metrics
#'
#' @param y A vector of observed values.
#' @param mupred A numeric array of posterior predictive means. For binary
#'   outcomes use a draws x observations matrix. For multiclass outcomes use a
#'   draws x observations x categories array.
#' 
#' @keywords internal
#' @name measure_metric_params
NULL

#' Shared parameters for scores
#'
#' @param y A vector of observed values.
#' @param mupred A numeric array of posterior predictive means. For binary
#'   outcomes use a draws x observations matrix. For multiclass outcomes use a
#'   draws x observations x categories array.
#' 
#' @keywords internal
#' @name measure_score_params
NULL

#' Pointwise log predictive density (`lppd_i`)
#'
#' Computes pointwise log predictive density contributions from a matrix of
#' log predictive densities/probabilities for posterior draws. When PSIS
#' log-weights are supplied, they are used to form a weighted log-sum-exp per
#' observation.
#'
#' @param ylp A numeric matrix of log predictive densities/probabilities with
#'   dimensions draws x observations.
#' @param psis_log_weights Optional numeric matrix of normalized PSIS log
#'   weights with the same dimensions as `ylp`. Each column must sum to 1 on
#'   the probability scale.
#'
#' @returns A numeric vector of length `ncol(ylp)` with pointwise log
#'   predictive density values.
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' ptw_log_pred_density(ylp)
#'
#' lw <- matrix(log(c(0.7, 0.3, 0.6, 0.4)), nrow = 2)
#' ptw_log_pred_density(ylp, lw)
#' @export
ptw_log_pred_density <- function(ylp, psis_log_weights = NULL) {
  .validate_numeric_matrix(ylp, arg = "ylp")
  n_draws <- dim(ylp)[1]
  
  if (is.null(psis_log_weights)) {
    return(matrixStats::colLogSumExps(ylp) - log(n_draws))
  }
  .validate_numeric_matrix(
    psis_log_weights,
    arg = "psis_log_weights",
    nrow = n_draws,
    ncol = ncol(ylp)
  )
  
  col_sums <- colSums(exp(psis_log_weights))
  if (!isTRUE(all.equal(col_sums, rep(1, ncol(psis_log_weights))))) {
    cli::cli_abort(c(
      "'psis_log_weights' must be normalized (column sums equal to 1).",
      "i" = "Range of current column sums: [{min(col_sums)}, {max(col_sums)}]."
    ))
  }
  
  matrixStats::colLogSumExps(ylp + psis_log_weights)
}

#' Expected log pointwise predictive density (`elpd`)
#'
#' Computes the expected log pointwise predictive density (ELPD) as the sum of
#' pointwise log predictive density contributions (\eqn{\mathrm{lppd}_i}), using
#' [ptw_log_pred_density()]. ELPD is returned on the utility scale (higher is
#' better), consistent with the sign convention used throughout this package.
#' Manual change of sign convention is possible via `higher_is_better`.
#'
#' @inheritParams measure_density_params
#' @inheritParams measure_params
#' @param pointwise Optional numeric vector of precomputed \eqn{\mathrm{lppd}_i}
#'   values. When supplied, `ylp` and `log_weights` are ignored.
#'
#' @details
#' \deqn{\mathrm{elpd} = \sum_{i=1}^{n} \mathrm{lppd}_i,}
#' where each \eqn{\mathrm{lppd}_i} is computed by [ptw_log_pred_density()].
#' The standard error is \eqn{\sqrt{n}\,\mathrm{sd}(\mathrm{lppd}_i)}.
#'
#' @seealso [ptw_log_pred_density()], [measure_mlpd()], [measure_ic()]
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' measure_elpd(ylp)
#'
#' # With unnormalized importance weights (e.g., PSIS-LOO)
#' lw <- matrix(log(c(0.7, 0.3, 0.6, 0.4)), nrow = 2)
#' measure_elpd(ylp, log_weights = lw)
#'
#' # From a draws x chains x observations array
#' LLarr <- example_loglik_array()
#' measure_elpd(LLarr)
#' @export
measure_elpd <- function(
  ylp, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {  
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = "measure_elpd"
    )
    lppd_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    .validate_numeric_matrix(ylp, arg = "ylp")
    ylp <- if (is.array(ylp) && length(dim(ylp)) == 3) llarray_to_matrix(ylp) else ylp
    n_draws <- nrow(ylp)
    n_obs <- ncol(ylp)
    if (!is.null(log_weights)) {
      log_weights <- .normalize_and_validate_log_weights(
        log_weights = log_weights, n_draws = n_draws, n_obs = n_obs
      )
    }
    lppd_i <- ptw_log_pred_density(ylp, log_weights)
  }
  
  if (length(lppd_i) == 1L) {
    cli::cli_warn("Only one pointwise value supplied; standard error is set to 0.")
  }
  
  res <- list(
    estimate = sum(lppd_i),
    se = if (length(lppd_i) == 1L) 0 else sqrt(length(lppd_i) * var(lppd_i)),
    pointwise = lppd_i
  )
  
  .create_measure_structure(
    res, higher_is_better, "elpd", n_draws = n_draws, n_obs = n_obs
  )
}

#' Mean log pointwise predictive density (`mlpd`)
#'
#' Computes MLPD as the average of pointwise log predictive density (lppd_i) 
#' values. Inputs follow the same conventions as [measure_elpd()].
#'
#' @inheritParams measure_density_params
#' @inheritParams measure_params
#' @param pointwise Optional numeric vector of precomputed \eqn{\mathrm{lppd}_i}
#'   values. When supplied, `ylp` and `log_weights` are ignored.
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' measure_mlpd(ylp)
#' @export
measure_mlpd <- function(
  ylp, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = "measure_mlpd"
    )
    lppd_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(ylp)
    n_obs <- ncol(ylp)
    .validate_numeric_matrix(ylp, arg = "ylp")
    ylp <- if (is.array(ylp) && length(dim(ylp)) == 3) llarray_to_matrix(ylp) else ylp
    
    if (!is.null(log_weights)) {
      log_weights <- .normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = nrow(ylp),
        n_obs = ncol(ylp)
      )
    }
    lppd_i <- ptw_log_pred_density(ylp, log_weights)
  }
  
  n_obs <- length(lppd_i)
  if (n_obs == 1L) {
    cli::cli_warn("Only one pointwise value supplied; standard error is set to 0.")
  }
  
  res <- list(
    estimate = sum(lppd_i) / n_obs,
    se = if (n_obs == 1L) 0 else sqrt(n_obs * var(lppd_i)) / n_obs,
    pointwise = lppd_i
  )
  .create_measure_structure(
    res, higher_is_better, "mlpd", n_draws = n_draws, n_obs = n_obs
  )
}

#' Information Criteria (`ic`)
#'
#' Computes the information criteria as -2 x log predictive density (lppd_i) 
#' values. Inputs follow the same conventions as [measure_elpd()] and
#' [measure_mlpd()].
#'
#' @inheritParams measure_density_params
#' @inheritParams measure_params
#' @param pointwise Optional numeric vector of precomputed pointwise
#'   contributions \eqn{\mathrm{ic}_i = -1 \cdot \mathrm{lppd}_i}. If provided, 
#'   `ylp` and `log_weights` are ignored.
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' measure_ic(ylp)
#' @export
measure_ic <- function(
  ylp, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = "measure_ic"
    )
    ic_i <- pointwise
    n_draws = NULL
    n_obs = length(pointwise)
  } else {
    n_draws <- nrow(ylp)
    n_obs <- ncol(ylp)
    .validate_numeric_matrix(ylp, arg = "ylp")
    ylp <- if (is.array(ylp) && length(dim(ylp)) == 3) llarray_to_matrix(ylp) else ylp
    if (!is.null(log_weights)) {
      log_weights <- .normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      )
    }
    lppd_i <- ptw_log_pred_density(ylp, log_weights)
    ic_i <- -2 * lppd_i
  }
  
  n_obs <- length(ic_i)
  if (n_obs == 1L) {
    cli::cli_warn("Only one pointwise value supplied; standard error is set to 0.")
  }
  
  res <- list(
    estimate = sum(ic_i),
    se = if (n_obs == 1L) 0 else 2 * sqrt(n_obs * var(ic_i / (-2))),
    pointwise = ic_i
  )
  .create_measure_structure(
    res, higher_is_better, "ic", n_draws = n_draws, n_obs = n_obs
  )
}

#' Classification accuracy (`acc`)
#'
#' Computes pointwise and average classification accuracy for binary or
#' multiclass outcomes from posterior predictive class assignments. For binary
#' outcomes, each draw is thresholded at 0.5. For multiclass outcomes, each
#' draw is mapped to the most likely category via `which.max()`.
#' 
#' @inheritParams measure_score_params
#' @inheritParams measure_params
#' @param y An integer vector of observed class labels.
#' @param pointwise Optional numeric vector of precomputed pointwise accuracy
#'   contributions. If provided, `y`, `mupred`, and `log_weights` are ignored.
#'
#' @examples
#' y <- c(1L, 0L, 1L)
#' mupred <- matrix(c(0.8, 0.3, 0.7, 0.6, 0.4, 0.9), nrow = 2)
#' measure_acc(y, mupred)
#' @export
measure_acc <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, mupred = mupred, log_weights = log_weights),
      fun_name = "acc"
    )
    acc_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(mupred)
    n_obs <- dim(mupred)[2]
    .validate_numeric_vector(y, arg = "y")
    if (!is.numeric(mupred) || (length(dim(mupred)) != 2 && length(dim(mupred)) != 3)) {
      cli::cli_abort("{.arg mupred} must be a numeric matrix or 3D numeric array.")
    }
    .validate_probs(mupred, arg = "mupred")
    
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
    } else {
      weights <- rep(1 / nrow(mupred), nrow(mupred))
    }
    
    if (length(dim(mupred)) == 3) {
      # Multiclass: (draws × obs × categories) > argmax over categories
      weighted_mupred <- apply(
        sweep(mupred, 1, weights, `*`),
        c(2, 3),
        sum
      )
      mupred_hat <- apply(weighted_mupred, 1, which.max)
    } else {
      weighted_mupred <- colSums(mupred * weights)
      mupred_hat <- (weighted_mupred > 0.5) * 1L
    }
    
    acc_i <- (mupred_hat == y) * 1L
  }
  
  res <- list(
    estimate = mean(acc_i),
    se = sqrt(mean(acc_i) * (1 - mean(acc_i)) / n_obs),
    pointwise = acc_i
  )
  .create_measure_structure(
    res, higher_is_better, "acc", n_draws = n_draws, n_obs = n_obs
  )
}

#' Balanced classification accuracy (`bacc`)
#'
#' Computes balanced accuracy by averaging class-specific mean accuracy, giving
#' each observed class equal weight regardless of class frequency.
#'
#' @inheritParams measure_acc
#'
#' @examples
#' y <- c(1L, 1L, 2L, 2L)
#' mupred <- array(
#'   c(0.8, 0.2, 0.7, 0.3, 0.3, 0.7, 0.2, 0.8),
#'   dim = c(1, 4, 2)
#' )
#' measure_bacc(y, mupred)
#' @export
measure_bacc <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "bacc"
    )
    acc_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(mupred)
    n_obs <- ncol(mupred)
    .validate_numeric_vector(y, arg = "y")
    classes <- sort(unique(y))
    K <- length(classes)
    if (K < 2) {
      cli::cli_abort("{.fn bacc} requires at least two outcome classes.")
    }
    if (!is.numeric(mupred) || (length(dim(mupred)) != 2 && length(dim(mupred)) != 3)) {
      cli::cli_abort("{.arg mupred} must be a numeric matrix or 3D numeric array.")
    }
    .validate_probs(mupred, arg = "mupred")
    
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
    } else {
      weights <- rep(1 / nrow(mupred), nrow(mupred))
    }
    if (length(dim(mupred)) == 3) {
      # Multiclass: (draws × obs × categories) > argmax over categories
      weighted_mupred <- apply(
        sweep(mupred, 1, weights, `*`),
        c(2, 3),
        sum
      )
      mupred_hat <- apply(weighted_mupred, 1, which.max)
    } else {
      .validate_numeric_matrix(mupred, arg = "mupred")
      weighted_mupred <- colSums(mupred * weights)
      mupred_hat <- (weighted_mupred > 0.5) * 1L
    }
    acc_i <- (mupred_hat == y) * 1L
  }
  
  # recomputed rather than reused from the branch above, which the
  # precomputed-`pointwise` branch never enters
  classes <- sort(unique(y))
  K <- length(classes)
  class_id <- match(y, classes)
  n_c <- tabulate(class_id, nbins = K)

  acc_c <- vapply(classes, function(c) mean(acc_i[y == c]), numeric(1))
  bacc_i <- acc_i / (K * n_c[class_id])

  res <- list(
    estimate = mean(acc_c),
    se = sqrt(sum(acc_c * (1 - acc_c) / n_c)) / K,
    pointwise = bacc_i,
    # the class strata are what makes a difference of balanced accuracies more
    # than a mean of pointwise differences, and `.se_diff_bacc()` cannot
    # recover them from `bacc_i`: every misclassified observation stores a
    # zero whatever its class
    extra = list(class_id = class_id)
  )
  .create_measure_structure(
    res, higher_is_better, "bacc", n_draws = n_draws, n_obs = n_obs
  )
}

#' Brier score (`brier`)
#'
#' Computes the Brier score for binary outcomes as squared error between the
#' observed label and predicted event probability.
#'
#' @param y A numeric vector of binary outcomes coded as 0 or 1.
#' @param ypred A numeric matrix of posterior predictive probabilities with
#'   dimensions draws x observations.
#' @inheritParams measure_params
#' @param pointwise Optional numeric vector of precomputed pointwise Brier
#'   scores. If provided, `y`, `ypred`, and `log_weights` are ignored.
#'
#' @examples
#' y <- c(1, 0, 1)
#' ypred <- matrix(c(0.8, 0.2, 0.7, 0.9, 0.4, 0.6), nrow = 2)
#' measure_brier(y, ypred)
#' @export
measure_brier <- function(
  y, ypred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, ypred = ypred, log_weights = log_weights),
      fun_name = "brier"
    )
    bs_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(ypred)
    n_obs <- ncol(ypred)
    .validate_numeric_vector(y, arg = "y")
    if (any(y != 0 & y != 1)) {
      cli::cli_abort(c(
        "The brier score expects binary data 'y'.",
        "i" = "Observed range: [{min(y)}, {max(y)}]",
        "x" = "All elements of {.arg y} must be 0 or 1."
      ))
    }
    .validate_numeric_matrix(ypred, arg = "ypred", ncol = length(y))
    .validate_probs(ypred, arg = "ypred")
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
      prob_i <- colSums(ypred * weights)
    } else {
      prob_i <- colMeans(ypred)
    }
    bs_i <- (prob_i - y)^2  
  }
  
  res <- list(
    estimate = mean(bs_i),
    se = sqrt(var(bs_i) / length(bs_i)),
    pointwise = bs_i
  )
  .create_measure_structure(
    res, higher_is_better, "brier", n_draws = n_draws, n_obs = n_obs
  )
}

#' Mean absolute error (`mae`)
#'
#' Computes MAE between observed outcomes and posterior predictive point
#' predictions. Point predictions are obtained by averaging `mupred` draws, or
#' by PSIS-weighted averaging when `log_weights` is provided.
#'
#' @param y A numeric vector of observed outcomes.
#' @param mupred A numeric matrix of posterior expected predictions with
#'   dimensions draws x observations. A length-`n` vector is also accepted.
#' @param pointwise Optional numeric vector of precomputed pointwise absolute
#'   errors. If provided, `y`, `mupred`, and `log_weights` are ignored.
#' @inheritParams measure_params
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' measure_mae(y, mupred)
#' @export
measure_mae <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "mae"
    )
    mae_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(mupred)
    n_obs <- ncol(mupred)
    .validate_numeric_vector(y, arg = "y")
    if (!is.null(mupred) && !is.matrix(mupred)){
      .validate_numeric_vector(mupred, arg = "mupred", len = length(y))
      cli::cli_inform(
        "Coercing {.arg mupred} from vector to 1 x n matrix for {.fn mae}."
      )
      mupred <- matrix(mupred, nrow = 1, ncol = length(mupred))
    }
    .validate_numeric_matrix(mupred, arg = "mupred", ncol = length(y))
    if (is.null(log_weights)) {
      mae_i <- abs(y - colMeans(mupred))
    } else {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
      mae_i <- abs(y - colSums(weights * mupred))
    }
  }
  
  res <- list(
    estimate = mean(mae_i),
    se = sqrt(var(mae_i) / length(mae_i)),
    pointwise = mae_i
  )
  .create_measure_structure(
    res, higher_is_better, "mae", n_draws = n_draws, n_obs = n_obs
  )
}

#' Mean squared error (`mse`)
#'
#' Computes MSE between observed outcomes and posterior predictive point
#' predictions. Point predictions are obtained by averaging `mupred` draws, or
#' by PSIS-weighted averaging when `log_weights` is provided.
#'
#' @inheritParams measure_mae
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' measure_mse(y, mupred)
#' @export
measure_mse <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {  
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "mse"
    )
    sqe_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  } else {
    n_draws <- nrow(mupred)
    n_obs <- ncol(mupred)
    .validate_numeric_vector(y, arg = "y")
    if (!is.null(mupred) && !is.matrix(mupred)){
      .validate_numeric_vector(mupred, arg = "mupred", len = length(y))
      cli::cli_inform(
        "Coercing {.arg mupred} from vector to 1 x n matrix for {.fn mse}."
      )
      mupred <- matrix(mupred, nrow = 1, ncol = length(mupred))
    }
    .validate_numeric_matrix(mupred, arg = "mupred", ncol = length(y))
    if (is.null(log_weights)) {
      sqe_i <- (y - colMeans(mupred))^2
    } else {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
      sqe_i <- (y - colSums(weights * mupred))^2
    }
  }

  res <- list(
    estimate = mean(sqe_i),
    se = sqrt(var(sqe_i) / length(sqe_i)),
    pointwise = sqe_i
  )
  .create_measure_structure(
    res, higher_is_better, "mse", n_draws = n_draws, n_obs = n_obs
  )
}

#' Root mean squared error (`rmse`)
#'
#' Computes RMSE as the square root of MSE and propagates uncertainty via a
#' first-order delta-method approximation.
#'
#' @inheritParams measure_mae
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' measure_rmse(y, mupred)
#' @export
measure_rmse <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  mse_res <- measure_mse(
    y = y, mupred = mupred, log_weights = log_weights, 
    pointwise = pointwise
  )
  n_draws <- if (is.null(pointwise)) nrow(mupred) else NULL
  n_obs <- length(mse_res$pointwise)
  
  sqe_i <- mse_res$pointwise
  rmse_est <- sqrt(mse_res$estimates[1])
  if (rmse_est == 0) {
    rmse_se <- 0
  } else {
    rmse_se <- mse_res$estimates[2] / (2 * rmse_est)
  }

  res <- list(
    estimate = rmse_est,
    se = rmse_se,
    pointwise = sqe_i
  )
  .create_measure_structure(
    res, higher_is_better, "rmse", n_draws = n_draws, n_obs = n_obs
  )
}

#' Delta-method standard error of an R-squared quantity
#'
#' The R-squared of one model and the difference in R-squared between two
#' models have the same form: a mean of squared-error contributions divided by
#' the model-independent baseline `MSE(y)`. Their standard errors are therefore
#' the same first-order Taylor approximation, evaluated either at one model's
#' pointwise squared errors or at the pointwise differences between two models.
#'
#' Writing `c = mean(sqe) / MSE(y)`, the three-term expansion
#' `Var[MSE] - 2 * c * Cov[MSE, MSE(y)] + c^2 * Var[MSE(y)]`, scaled by
#' `MSE(y)^-1`, is exactly the standard error of the mean of
#' `sqe_i - c * mse_y_i`. That is the form used here: it needs one variance
#' rather than three moments, it cannot go negative under the square root, and
#' it is exactly `0` when `sqe` is identically zero, as it is when a model is
#' compared against itself.
#'
#' @noRd
#' @param sqe Pointwise squared errors of one model, or pointwise differences
#'   in squared error between two models.
#' @param mse_y_i Pointwise baseline `(y_i - mean(y))^2`.
#' @return Numeric scalar standard error.
.se_r2_delta <- function(sqe, mse_y_i) {
  mse_y_hat <- mean(mse_y_i)
  scaled <- sqe - (mean(sqe) / mse_y_hat) * mse_y_i
  sqrt(var(scaled) / length(sqe)) / mse_y_hat
}

#' Predictive R-squared (`r2`)
#'
#' Computes predictive R-squared as one minus the ratio of prediction MSE to
#' the empirical variance of `y`. The standard error is computed with a
#' first-order delta-method approximation.
#'
#' @inheritParams measure_mae
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' measure_r2(y, mupred)
#' @export
measure_r2 <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, higher_is_better = NULL
) {
  .validate_numeric_vector(y, arg = "y")
  if (var(y) == 0) {
    cli::cli_abort(
      "{.fn r2} is undefined when {.arg y} has zero variance."
    )
  }
  
  mse_res <- measure_mse(
    y = y, 
    mupred = mupred, 
    log_weights = log_weights,
    pointwise = pointwise
  )
  mse_hat <- mse_res$estimate[1]
  sqe_i <- mse_res$pointwise
  n_obs <- length(sqe_i)
  n_draws <- if (is.null(pointwise)) nrow(mupred) else NULL
  
  mse_y_i <- (y - mean(y))^2
  mse_y_hat <- mean(mse_y_i)

  est_r2 <- 1 - mse_hat / mse_y_hat
  se_r2 <- .se_r2_delta(sqe_i, mse_y_i)

  res <- list(
    estimate = est_r2,
    se = se_r2,
    pointwise = sqe_i,
    # `loo_compare()` needs the baseline to propagate uncertainty into the
    # standard error of an r2 difference; `y` is gone by then. See
    # `.se_diff_r2()`.
    extra = list(mse_y_i = mse_y_i)
  )
  .create_measure_structure(
    res, higher_is_better, "r2", n_draws = n_draws, n_obs = n_obs
  )
}

#' Ranked Probability Score (RPS, SRPS, CRPS, SCRPS)
#'
#' Computes proper scoring rules based on the ranked probability score family,
#' covering both discrete and continuous outcomes, with optional scaling.
#' The specific score computed depends on the type of `y` and `ypred` and the
#' value of `scaled`:
#'
#' | `y`/`ypred` type | `scaled = FALSE` | `scaled = TRUE` |
#' |------------------|-----------------|-----------------|
#' | Discrete         | RPS             | SRPS            |
#' | Continuous       | CRPS            | SCRPS           |
#'
#' **Scoring rules:**
#'
#' - **RPS** (Epstein, 1969): Compares predictive and observed cumulative
#'   distributions for ordered discrete outcomes.
#' - **CRPS** (Matheson & Winkler, 1976; Gneiting & Raftery, 2007): Generalizes
#'   RPS to continuous outcomes. Defined as
#'   \deqn{\mathrm{CRPS}(X; y) = \frac{1}{2} E[|X - X'|] - E[|X - y|],}
#'   where \eqn{X, X'} are independent draws from the predictive distribution.
#' - **SRPS/SCRPS** (Bolin & Wallin, 2023): Scaled variants that are invariant
#'   to the scale of the predictive distribution. Defined as
#'   \deqn{\mathrm{SCRPS}(X; y) = -\frac{E[|X - y|]}{E[|X - X'|]} -
#'     \frac{1}{2} \log E[|X - X'|].}
#'
#' **Estimation:**
#'
#' Scores are estimated using the probability-weighted moment (PWM) estimator
#' (Taillardat et al., 2016; Zamo & Naveau, 2018), which requires only a single
#' set of predictive draws unlike permutation-based estimators, which require
#' two independent draw sets. The PWM estimator is unbiased and generally more
#' accurate than single-permutation estimators. The same estimator is used for
#' both discrete and continuous outcomes; see Hosking (1990, 1996) for
#' theoretical justification in the discrete case.
#'
#' If log-weights (`log_weights`) are provided (e.g., PSIS weights
#' for LOO cross-validation), a weighted PWM estimator is used instead, which
#' accounts for the importance weights when estimating expectations.
#'
#' **Sign convention:**
#'
#' Unscaled scores are returned as utilities (higher is better), consistent 
#' with the sign convention of log score / ELPD. Scaled scores are negatively
#' oriented (lower is better). Use `higher_is_better = NULL` to obtain the
#' loss convention (lower is better) used in some references for unscaled scores,
#' or `higher_is_better = TRUE` for a utility scale when using scaled scores.
#'
#' @param y A numeric vector of \eqn{n} observed outcomes. May be integer-valued
#'   (for RPS/SRPS) or continuous (for CRPS/SCRPS).
#' @param ypred A numeric matrix of posterior predictive draws with dimensions
#'   \eqn{S \times n} (draws × observations).
#' @param pointwise Optional numeric vector of precomputed pointwise rps values.
#'   If provided, `y`, `ypred`, and `log_weights` are ignored.
#' @param scaled Logical; if `TRUE`, computes the scaled variant (SRPS for
#'   discrete outcomes, SCRPS for continuous outcomes). Default is `FALSE`.
#' @inheritParams measure_params
#'
#' @examples
#' # Discrete outcomes: RPS
#' y <- c(2L, 1L, 3L)
#' ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
#' measure_rps(y, ypred)
#'
#' # Discrete outcomes: SRPS (scaled)
#' measure_rps(y, ypred, scaled = TRUE)
#'
#' # Continuous outcomes: CRPS
#' y_cont <- c(0.5, -1.2, 2.3)
#' ypred_cont <- matrix(rnorm(200), nrow = 100, ncol = 3)
#' measure_rps(y_cont, ypred_cont)
#'
#' # With importance weights: LOO-CRPS
#' log_weights <- matrix(rnorm(200), nrow = 100, ncol = 3)
#' measure_rps(y_cont, ypred_cont, log_weights = log_weights)
#'
#' @references
#' Bolin, D. and Wallin, J. (2023). Local scale invariance and robustness of
#' proper scoring rules. *Statistical Science*, 38(1):140–159.
#'
#' Epstein, E. S. (1969). A scoring system for probability forecasts of ranked
#' categories. *Journal of Applied Meteorology*, 8(6):985–987.
#'
#' Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring rules,
#' prediction, and estimation. *Journal of the American Statistical
#' Association*, 102(477):359–378.
#'
#' Hosking, J. R. M. (1990). L-moments: Analysis and estimation of
#' distributions using linear combinations of order statistics. *Journal of
#' the Royal Statistical Society Series B*, 52(1):105–124.
#'
#' Hosking, J. R. M. (1996). Some theoretical results concerning L-moments.
#' Research report RC 14492. IBM Thomas J. Watson Research Division.
#'
#' Matheson, J. E. and Winkler, R. L. (1976). Scoring rules for continuous
#' probability distributions. *Management Science*, 22(10):1087–1096.
#'
#' Taillardat, M., Mestre, O., Zamo, M., and Naveau, P. (2016). Calibrated
#' ensemble forecasts using quantile regression forests and ensemble model
#' output statistics. *Monthly Weather Review*, 144(6):2375–2393.
#'
#' Zamo, M. and Naveau, P. (2018). Estimation of the continuous ranked
#' probability score with limited information and applications to ensemble
#' weather forecasts. *Mathematical Geosciences*, 50:209–234.
#'
#' @export
measure_rps <- function(y, ypred, log_weights = NULL, pointwise = NULL, scaled = FALSE, 
  higher_is_better = NULL) {
  if (is.null(pointwise)) {
    n_draws <- nrow(ypred)
    n_obs <- ncol(ypred)
    
    if (is.null(log_weights)) {
      EXy <- colMeans(abs(sweep(ypred, 2, y)))
      ypred_sorted <- apply(ypred, 2, sort)
      EXX  <- colMeans(ypred_sorted * ((1:n_draws) * (4 / (n_draws - 1)) - 2))
      if (scaled) {
        # Scaled version by Bolin & Wallin (2023)
        rps_i <- -EXy/EXX - 0.5 * log(EXX)
      } else {
        # Gneiting & Raftery (2007)
        rps_i <- - 0.5 * EXX + EXy
      }
    } else {
      w <- exp(.normalize_log_weights(log_weights))
      w_csum <- sapply(1:n_obs, function(j) {
        perm <- order(ypred[,j])
        result <- numeric(n_draws)
        result[perm] <- cumsum(w[perm, j])
        result
      })
      
      EXX <- 2 * colSums(ypred * w * (1 - (2*(1 - w_csum)) / (1 - w)))
      EXy <- colSums(w * abs(sweep(ypred, 2, y)))
      if (scaled) {
        # Scaled version by Bolin & Wallin (2023)
        rps_i <- -EXy/EXX - 0.5 * log(EXX)
      } else {
        # Gneiting & Raftery (2007)
        rps_i <- EXy - 0.5 * EXX
      }
    }
  } else {
    rps_i <- pointwise
    n_draws <- NULL
    n_obs <- length(pointwise)
  }
  
  res <- list(
    estimate = mean(rps_i),
    se = sqrt(var(rps_i) / n_obs),
    pointwise = rps_i
  )
  name <- if(isTRUE(scaled)) "srps" else "rps"
  .create_measure_structure(
    res, higher_is_better, name, n_draws = n_draws, n_obs = n_obs
  )
}

#' Scaled Ranked Probability Score (SRPS, SCRPS)
#'
#' A convenience wrapper around [measure_rps()] with `scaled = TRUE`. Computes the
#' scaled ranked probability score (SRPS) for discrete outcomes or the scaled
#' continuously ranked probability score (SCRPS) for continuous outcomes.
#' Scaling makes the score invariant to the spread of the predictive
#' distribution, which can be useful when comparing models across observations
#' with very different predictive uncertainties.
#'
#' See [measure_rps()] for full details on arguments, estimation, and references.
#'
#' @inheritParams measure_rps
#'
#' @examples
#' y <- c(2L, 1L, 3L)
#' ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
#' measure_srps(y, ypred)
#'
#' @export
measure_srps <- function(y, ypred, log_weights = NULL, pointwise = NULL,
  higher_is_better = NULL) {
  measure_rps(
    y = y, ypred = ypred, log_weights = log_weights, 
    scaled = TRUE, higher_is_better = higher_is_better
  )
}

# measure overview -----------------------------
#
# Built-in measures are registered in `.measure_spec`. Users can also pass
# custom functions to `measure` in the pred_measure family; those functions
# must set `attr(fun, "measure_name")` and return `estimate`, `se`, and
# `pointwise` (see `?insample_pred_measure`).

# internal function to get the measure specification
# @noRd
# @param measure The measure used.
# @return The measure specification.
.measure_spec <- list(
  elpd = list(fun = measure_elpd, loss = FALSE, diff_method = "sum"),
  ic = list(fun = measure_ic, loss = TRUE, diff_method = "sum"),
  mlpd = list(fun = measure_mlpd, loss = FALSE, diff_method = "sum"),
  mae = list(fun = measure_mae, loss = TRUE, diff_method = "mean"),
  r2 = list(
    fun = measure_r2,
    loss = FALSE,
    diff_method = "pairwise",
    se_diff_fun = "r2"
  ),
  rmse = list(
    fun = measure_rmse,
    loss = TRUE,
    diff_method = "pairwise",
    se_diff_fun = "rmse"
  ),
  mse = list(fun = measure_mse, loss = TRUE, diff_method = "mean"),
  acc = list(fun = measure_acc, loss = FALSE, diff_method = "mean"),
  bacc = list(
    fun = measure_bacc,
    loss = FALSE,
    diff_method = "pairwise",
    se_diff_fun = "bacc"
  ),
  rps = list(fun = measure_rps, loss = FALSE, diff_method = "mean"),
  srps = list(fun = measure_srps, loss = TRUE, diff_method = "mean"),
  brier = list(fun = measure_brier, loss = TRUE, diff_method = "mean")
)

# pairwise standard errors -----------------------------
#
# Measures whose overall estimate is not a sum or mean of pointwise
# contributions cannot use the paired pointwise standard error. They register a
# `se_diff_fun` in `.measure_spec`, naming an entry of `.se_diff_funs` below.
# Custom measures may instead attach a function directly via
# `attr(my_fun, "se_diff_fun")`.
#
# Such a function receives `ref` and `cmp`, each a list with the elements
# `estimate`, `se`, `pointwise`, and `extra` for one model, always on the
# measure's natural scale (`higher_is_better` sign flips are undone before the
# call and reapplied to the difference afterwards), and returns the standard
# error of the difference as a numeric scalar. The difference itself is always
# `estimate_cmp - estimate_ref` and is computed by `loo_compare()`.

#' Standard error of an RMSE difference
#'
#' First-order bivariate Taylor (delta method) approximation of the standard
#' error of \eqn{RMSE(M_cmp) - RMSE(M_ref)}, propagated from the MSE scale on
#' which the pointwise squared errors live.
#'
#' @noRd
#' @param ref,cmp Per-model inputs; `pointwise` holds squared errors.
#' @return Numeric scalar standard error.
.se_diff_rmse <- function(ref, cmp) {
  sqe_ref <- ref$pointwise
  sqe_cmp <- cmp$pointwise
  n <- length(sqe_ref)
  mse_ref <- mean(sqe_ref)
  mse_cmp <- mean(sqe_cmp)

  # a perfect predictor leaves the ratios below undefined; `measure_rmse()`
  # reports a zero standard error in that case, so do the same here
  if (n <= 1L || mse_ref <= 0 || mse_cmp <= 0) {
    return(0)
  }

  se_mse_ref <- sqrt(stats::var(sqe_ref) / n)
  se_mse_cmp <- sqrt(stats::var(sqe_cmp) / n)
  cov_mse <- sum((sqe_cmp - mse_cmp) * (sqe_ref - mse_ref)) / (n * (n - 1))

  # standard errors of the two MSEs relative to their own scale, and their
  # correlation
  rel_ref <- se_mse_ref / sqrt(mse_ref)
  rel_cmp <- se_mse_cmp / sqrt(mse_cmp)
  rho <- if (se_mse_ref <= 0 || se_mse_cmp <= 0) {
    0
  } else {
    cov_mse / (se_mse_cmp * se_mse_ref)
  }

  # algebraically `rel_cmp^2 + rel_ref^2 - 2 * rho * rel_cmp * rel_ref`, but
  # written so that two models with the same squared errors cancel exactly
  # rather than leaving rounding noise behind. Both terms are non-negative,
  # because the correlation cannot exceed one.
  v <- (rel_cmp - rel_ref)^2 + 2 * rel_cmp * rel_ref * max(1 - rho, 0)

  0.5 * sqrt(v)
}

#' Standard error of an R-squared difference
#'
#' First-order trivariate Taylor (delta method) approximation of the standard
#' error of \eqn{R^2(M_cmp) - R^2(M_ref)}. The difference equals
#' \eqn{-MSE(M_cmp, M_ref) / MSE(y)}, so it is the same expansion as the
#' single-model standard error in `measure_r2()` with the pointwise squared
#' errors replaced by their pointwise differences; both go through
#' `.se_r2_delta()`.
#'
#' @noRd
#' @param ref,cmp Per-model inputs; `pointwise` holds squared errors and
#'   `extra$mse_y_i` the baseline `(y_i - mean(y))^2` stored by `measure_r2()`.
#' @return Numeric scalar standard error, or `NA_real_` when the baseline is
#'   unavailable.
.se_diff_r2 <- function(ref, cmp) {
  # the baseline is a property of `y`, so either model's copy will do; models
  # fitted to different `y` are already reported by the `yhash` warning
  mse_y_i <- ref$extra$mse_y_i
  if (is.null(mse_y_i)) {
    mse_y_i <- cmp$extra$mse_y_i
  }

  # objects computed before the baseline was stored cannot support the
  # covariance terms; report the difference without a standard error rather
  # than refusing the whole comparison
  if (!is.numeric(mse_y_i) || length(mse_y_i) != length(ref$pointwise)) {
    return(NA_real_)
  }

  .se_r2_delta(cmp$pointwise - ref$pointwise, mse_y_i)
}

#' Standard error of a balanced-accuracy difference
#'
#' Balanced accuracy averages class-wise accuracies, so a difference of two
#' balanced accuracies is a difference of two stratified means, not a mean of
#' pointwise differences. The two sources of dependence separate: the class
#' strata are disjoint sets of observations and so contribute independent
#' variances, while within a stratum both models score the *same* `n_c`
#' observations and are therefore paired. Writing
#' \eqn{d_i = acc_i(M_cmp) - acc_i(M_ref)},
#'
#' \deqn{SE = \frac{1}{K} \sqrt{\sum_c Var(d_i : i \in c) / n_c}}
#'
#' which is the difference-analogue of the single-model
#' \eqn{\sqrt{\sum_c acc_c (1 - acc_c) / n_c} / K} in `measure_bacc()`: the
#' per-stratum binomial variance replaced by the paired-difference variance.
#' Within a stratum this is the McNemar variance of a paired difference of
#' proportions, \eqn{(b + c)/n_c^2 - (b - c)^2/n_c^3} in discordant-pair form,
#' up to the \eqn{n_c/(n_c - 1)} of the sample variance.
#'
#' For binary outcomes balanced accuracy is \eqn{(sens + spec)/2}, so this is
#' the estimand of Newcombe (2001) at a mixing parameter of one half, and half
#' the difference of two Youden indices in a paired design (Chen et al., 2015).
#'
#' @references
#' Newcombe, R. G. (2001). Simultaneous comparison of sensitivity and
#' specificity of two tests in the paired design: a straightforward graphical
#' approach. *Statistics in Medicine*, 20(6):907--915.
#'
#' Chen, F., Xue, Y., Tan, M. T., and Chen, P. (2015). Efficient statistical
#' tests to compare Youden index: accounting for contingency correlation.
#' *Statistics in Medicine*, 34(9):1560--1576.
#'
#' @noRd
#' @param ref,cmp Per-model inputs; `pointwise` holds `acc_i / (K * n_c)` and
#'   `extra$class_id` the class index stored by `measure_bacc()`.
#' @return Numeric scalar standard error, or `NA_real_` when the class strata
#'   are unavailable.
.se_diff_bacc <- function(ref, cmp) {
  # the strata are a property of `y`, so either model's copy will do; models
  # fitted to different `y` are already reported by the `yhash` warning
  class_id <- ref$extra$class_id
  if (is.null(class_id)) {
    class_id <- cmp$extra$class_id
  }

  n <- length(ref$pointwise)
  # objects computed before the strata were stored cannot be stratified;
  # report the difference without a standard error rather than refusing the
  # whole comparison
  if (!is.numeric(class_id) || length(class_id) != n) {
    return(NA_real_)
  }

  n_c <- tabulate(class_id)
  K <- length(n_c)
  # undo the `acc_i / (K * n_c)` scaling to recover the 0/1 accuracies, so the
  # variances below are on the natural per-observation scale
  d <- (cmp$pointwise - ref$pointwise) * (K * n_c[class_id])

  # a stratum holding a single observation supports no variance estimate and
  # contributes nothing, matching `measure_bacc()`, where `acc_c` is then 0 or
  # 1 and its binomial variance vanishes
  var_c <- vapply(
    seq_len(K),
    function(k) {
      if (n_c[k] < 2L) 0 else stats::var(d[class_id == k]) / n_c[k]
    },
    numeric(1)
  )

  sqrt(sum(var_c)) / K
}

# registry of built-in `se_diff_fun` implementations, referenced by name from
# `.measure_spec` so that stored objects carry a string rather than a closure
.se_diff_funs <- list(
  rmse = .se_diff_rmse,
  r2 = .se_diff_r2,
  bacc = .se_diff_bacc
)

#' Return comparison metadata for a measure
#' @noRd
#' @param measure_entry Normalized measure entry, or a built-in measure name.
#' @param higher_is_better Value of `higher_is_better` used for this measure.
.measure_compare_meta <- function(measure_entry, higher_is_better = NULL) {
  if (is.character(measure_entry)) {
    measure_entry <- list(
      name = measure_entry,
      type = "builtin",
      key = measure_entry
    )
  }

  if (measure_entry$type == "builtin") {
    entry <- .measure_spec[[measure_entry$key]]
    if (is.null(entry)) {
      return(NULL)
    }
    return(list(
      higher_is_better = higher_is_better,
      loss = isTRUE(entry$loss),
      diff_method = entry$diff_method,
      se_diff_fun = entry$se_diff_fun
    ))
  }

  se_diff_fun <- attr(measure_entry$key, "se_diff_fun", exact = TRUE)
  if (!is.null(se_diff_fun) && !is.function(se_diff_fun)) {
    cli::cli_abort(c(
      "Attribute {.code se_diff_fun} of custom measure",
      "{.val {measure_entry$name}} must be a function.",
      "i" = "It is called as {.code se_diff_fun(ref, cmp)} and must return a",
      " " = "numeric scalar."
    ))
  }

  list(
    higher_is_better = higher_is_better,
    loss = FALSE,
    diff_method = if (is.null(se_diff_fun)) "auto" else "pairwise",
    se_diff_fun = se_diff_fun
  )
}

#' Supported predictive measure names
#'
#' A character vector of measure names that can be passed to the `measure`
#' argument of [insample_pred_measure()], [loo_pred_measure()],
#' [kfold_pred_measure()], [test_pred_measure()], and [pred_measure()].
#'
#' @export
supported_measures_list <- names(.measure_spec)

# internal function that produces output format for measures
.create_measure_structure <- function(
  res, higher_is_better, measure_name, n_draws, n_obs
) {
  if (!is.null(higher_is_better)) {
    spec <- .measure_spec[[measure_name]]
    natural_higher <- is.null(spec) || !isTRUE(spec$loss)
    if (xor(natural_higher, isTRUE(higher_is_better))) {
      res$estimate <- -res$estimate
      res$pointwise <- -res$pointwise
    }
  }
  out <- list()
  out$estimates <- matrix(
    c(res$estimate, res$se),
    nrow = 1,
    dimnames = list(measure_name, c("Estimate", "SE"))
  )
  out$pointwise <- matrix(
    res$pointwise,
    ncol     = 1,
    dimnames = list(NULL, measure_name)
  )

  structure(
    out,
    class = c("measure", "loo"),
    measure = measure_name,
    dims = c(n_draws, n_obs),
    higher_is_better = higher_is_better,
    # not affected by the sign flip above: `extra` carries auxiliary data for
    # `se_diff_fun()`, which always works on the measure's natural scale
    compare_extra = res$extra
  )
}

#' Auxiliary data a measure stores for its `se_diff_fun`
#'
#' Built-in measures carry it as the `compare_extra` attribute added by
#' `.create_measure_structure()`; custom measures return it as an `extra`
#' element of their result list.
#' @noRd
#' @param res A measure result.
#' @return A list, or `NULL` when the measure stores nothing.
.measure_compare_extra <- function(res) {
  extra <- attr(res, "compare_extra", exact = TRUE)
  if (is.null(extra)) {
    extra <- res$extra
  }
  extra
}
