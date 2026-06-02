## Density scores -------------------------------------------------------------

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
#' Computes ELPD as the sum of pointwise log predictive density contributions
#' (lppd_i).
#'
#' @param ylp A numeric matrix of log predictive densities/probabilities with
#'   dimensions draws x observations.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise
#'   contributions. If provided, `ylp` and `log_weights` are ignored.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (sum over lppd_i), `se`
#'   (standard error across lppd_i), and `pointwise` (lppd_i).
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' elpd(ylp)
#' @export
elpd <- function(
  ylp, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {  
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = "elpd"
    )
    lppd_i <- pointwise
  } else {
    .validate_numeric_matrix(ylp, arg = "ylp")
    if (!is.null(log_weights)) {
      log_weights <- .normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = nrow(ylp),
        n_obs = ncol(ylp)
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
  
  .measure_output(res, revert_sign)
}

#' Mean log pointwise predictive density (`mlpd`)
#'
#' Computes MLPD as the average of pointwise log predictive density (lppd_i) 
#' values. Inputs follow the same conventions as [elpd()].
#'
#' @param ylp A numeric matrix of log predictive densities/probabilities with
#'   dimensions draws x observations.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise
#'   contributions. If provided, `ylp` and `log_weights` are ignored.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean over lppd_i), `se`
#'   (standard error), and `pointwise` (lppd_i).
#'
#' @examples
#' ylp <- matrix(log(c(0.2, 0.4, 0.3, 0.8)), nrow = 2)
#' mlpd(ylp)
#' @export
mlpd <- function(
  ylp, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .validate_numeric_vector(pointwise, arg = "pointwise")
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(ylp = ylp, log_weights = log_weights),
      fun_name = "mlpd"
    )
    lppd_i <- pointwise
  } else {
    .validate_numeric_matrix(ylp, arg = "ylp")
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
  .measure_output(res, revert_sign)
}

#' Classification accuracy (`acc`)
#'
#' Computes pointwise and average classification accuracy for binary or
#' multiclass outcomes from posterior predictive class assignments. For binary
#' outcomes, each draw is thresholded at 0.5. For multiclass outcomes, each
#' draw is mapped to the most likely category via `which.max()`.
#'
#' @param y An integer vector of observed class labels.
#' @param mupred A numeric array of posterior predictive means. For binary
#'   outcomes use a draws x observations matrix. For multiclass outcomes use a
#'   draws x observations x categories array.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise accuracy
#'   contributions. If provided, `y`, `mupred`, and `log_weights` are ignored.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean accuracy), `se`
#'   (standard error), and `pointwise` (pointwise accuracy).
#'
#' @examples
#' y <- c(1L, 0L, 1L)
#' mupred <- matrix(c(0.8, 0.3, 0.7, 0.6, 0.4, 0.9), nrow = 2)
#' loo::acc(y, mupred)
#' @export
acc <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, mupred = mupred, log_weights = log_weights),
      fun_name = "acc"
    )
    acc_i <- pointwise
  } else {
    .validate_numeric_vector(y, arg = "y")
    if (!is.numeric(mupred) || (length(dim(mupred)) != 2 && length(dim(mupred)) != 3)) {
      cli::cli_abort("{.arg mupred} must be a numeric matrix or 3D numeric array.")
    }
    .validate_probs(mupred, arg = "mupred")
    
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = nrow(mupred),
        n_obs = dim(mupred)[2]
      ))
    } else {
      weights <- rep(1 / nrow(mupred), nrow(mupred))
    }
    
    if (length(dim(mupred)) == 3) {
      # Multiclass: (draws × obs × categories) > argmax over categories
      ypred <- apply(mupred, c(1, 2), which.max)
    } else {
      .validate_numeric_matrix(mupred, arg = "mupred")
      # Binary: (draws × obs) > threshold at 0.5
      ypred <- (mupred > 0.5) * 1L
    }
    acc_i <- colSums(sweep(ypred, 2, y, `==`) * weights)
  }
  
  res <- list(
    estimate = mean(acc_i),
    se = sqrt(var(acc_i) / length(acc_i)),
    pointwise = acc_i
  )
  .measure_output(res, revert_sign)
}

#' Balanced classification accuracy (`bacc`)
#'
#' Computes balanced accuracy by averaging class-specific mean accuracy, giving
#' each observed class equal weight regardless of class frequency.
#'
#' @inheritParams acc
#'
#' @returns A named list with elements `estimate` (balanced accuracy), `se`
#'   (standard error based on balanced pointwise contributions), and
#'   `pointwise` (pointwise balanced accuracy terms).
#'
#' @examples
#' y <- c(1L, 1L, 2L, 2L)
#' mupred <- array(
#'   c(0.8, 0.2, 0.7, 0.3, 0.3, 0.7, 0.2, 0.8),
#'   dim = c(1, 4, 2)
#' )
#' bacc(y, mupred)
#' @export
bacc <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "bacc"
    )
    acc_i <- pointwise
  } else {
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
        n_draws = nrow(mupred),
        n_obs = dim(mupred)[2]
      ))
    } else {
      weights <- rep(1 / nrow(mupred), nrow(mupred))
    }
    if (length(dim(mupred)) == 3) {
      # Multiclass: (draws × obs × categories) > argmax over categories
      ypred <- apply(mupred, c(1, 2), which.max)
    } else {
      .validate_numeric_matrix(mupred, arg = "mupred")
      # Binary: (draws × obs) > threshold at 0.5
      ypred <- (mupred > 0.5) * 1L
    }
    acc_i <- colSums(sweep(ypred, 2, y, `==`) * weights)
  }
  
  acc_c <- vapply(classes, function(c) mean(acc_i[y == c]), numeric(1))
  n_c <- tabulate(match(y, classes))
  bacc_i <- acc_i / (K * n_c[match(y, classes)])
  
  res <- list(
    estimate = mean(acc_c),
    se = sqrt(var(bacc_i) / length(bacc_i)),
    pointwise = acc_i
  )
  .measure_output(res, revert_sign)
}

#' Brier score (`brier`)
#'
#' Computes the Brier score for binary outcomes as squared error between the
#' observed label and predicted event probability.
#'
#' @param y A numeric vector of binary outcomes coded as 0 or 1.
#' @param ypred A numeric matrix of posterior predictive probabilities with
#'   dimensions draws x observations.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise Brier
#'   scores. If provided, `y`, `ypred`, and `log_weights` are ignored.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean Brier score), `se`
#'   (standard error), and `pointwise` (pointwise brier score).
#'
#' @examples
#' y <- c(1, 0, 1)
#' ypred <- matrix(c(0.8, 0.2, 0.7, 0.9, 0.4, 0.6), nrow = 2)
#' brier(y, ypred)
#' @export
brier <- function(
  y, ypred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, ypred = ypred, log_weights = log_weights),
      fun_name = "brier"
    )
    bs_i <- pointwise
  } else {
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
        n_draws = nrow(ypred),
        n_obs = ncol(ypred)
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
  .measure_output(res, revert_sign)
}

#' Ranked probability score (`rps`)
#'
#' Computes the ranked probability score (RPS) for ordered discrete outcomes by
#' comparing predictive and observed cumulative distribution functions.
#'
#' @param y An integer vector of observed outcomes.
#' @param ypred A numeric matrix of posterior predictive draws with dimensions
#'   draws x observations.
#' @param size Maximum support value for `y`. The function assumes support
#'   values from 0 to `size`.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise scores.
#'   If provided, `y`, `ypred`, and `log_weights` are ignored.
#' @param scale Logical; if `TRUE`, computes the scaled variant used by
#'   [srps()].
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean RPS score), `se`
#'   (standard error), and `pointwise` (pointwise RPS scores).
#'
#' @examples
#' y <- c(2, 1, 3)
#' ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
#' rps(y, ypred, size = 3)
#' @export
rps <- function(
  y, ypred, size, log_weights = NULL, pointwise = NULL, scale = FALSE,
  revert_sign = FALSE
) {
  if (!is.numeric(size) || length(size) != 1 || !is.finite(size) ||
      size < 1 || size != as.integer(size)) {
    cli::cli_abort("{.arg size} must be a single positive integer.")
  }
  
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, ypred = ypred, log_weights = log_weights, 
        size = size),
      fun_name = "rps"
    )
    s_rps_i <- pointwise
  } else {
    .validate_numeric_vector(y, arg = "y")
    .validate_numeric_matrix(ypred, arg = "ypred", ncol = length(y))
    if (any(y < 0 | y > size) || any(y != as.integer(y))) {
      cli::cli_abort("{.arg y} must contain integers between 0 and {.val {size}}.")
    }
    if (any(ypred < 0 | ypred > size) || any(ypred != as.integer(ypred))) {
      cli::cli_abort(
        "{.arg ypred} must contain integer predictive draws between 0 and {.val {size}}."
      )
    }
    
    n_draws <- nrow(ypred)
    p_mat <- matrix(0, nrow = ncol(ypred), ncol = size + 1)
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = ncol(ypred)
      ))
      for (k in 1:size) {
        p_mat[, k] <- colSums((ypred == k) * weights)
      }
    } else {
      weights <- rep(1 / n_draws, n_draws)
      for (k in 1:size) {
        p_mat[, k] <- colSums(ypred == k) / n_draws
      }
    }

    predicted <- matrixStats::rowCumsums(p_mat)[, 1:size, drop = FALSE]
    observed  <- outer(y, 0:(size - 1), `<=`) * 1L

    s_rps_i <- rowSums((predicted - observed)^2)
    
    if (scale) {
      # E[|Y - y_i|]
      EYy <- colSums(abs(sweep(ypred, 2, y)) * weights)
      # E[|Y - Y'|] = 2 * sum_k F_k * (1 - F_k)
      EYY <- 2 * rowSums(predicted * (1 - predicted))
      if (any(EYY <= 0)) {
        cli::cli_abort(
          "Cannot compute scaled RPS because predictive spread is zero for at least one observation."
        )
      }
      s_rps_i <- -EYy / EYY - 0.5 * log(EYY)
    }
  }

  res <- list(
    estimate = mean(s_rps_i),
    se = sqrt(var(s_rps_i) / length(s_rps_i)),
    pointwise = s_rps_i
  )
  .measure_output(res, revert_sign)
}

#' Scaled ranked probability score (`srps`)
#'
#' Wrapper around [rps()] that computes the scaled ranked probability score
#' (SRPS), a positively oriented score based on normalized absolute error and
#' predictive spread.
#'
#' @inheritParams rps
#'
#' @returns A named list with elements `estimate`, `se`, and `pointwise`.
#'
#' @examples
#' y <- c(2, 1, 3)
#' ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
#' srps(y, ypred, size = 3)
#' @export
srps <- function(
  y, ypred, size, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  rps(
    y = y, 
    ypred = ypred, 
    size = size, 
    log_weights = log_weights, 
    scale = TRUE,
    pointwise = pointwise,
    revert_sign = revert_sign
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
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise absolute
#'   errors. If provided, `y`, `mupred`, and `log_weights` are ignored.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean absolute error), `se`
#'   (standard error), and `pointwise` (pointwise absolute errors).
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' mae(y, mupred)
#' @export
mae <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "mae"
    )
    mae_i <- pointwise
  } else {
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
        n_draws = nrow(mupred),
        n_obs = ncol(mupred)
      ))
      mae_i <- abs(y - colSums(weights * mupred))
    }
  }
  
  res <- list(
    estimate = mean(mae_i),
    se = sqrt(var(mae_i) / length(mae_i)),
    pointwise = mae_i
  )
  .measure_output(res, revert_sign)
}

#' Mean squared error (`mse`)
#'
#' Computes MSE between observed outcomes and posterior predictive point
#' predictions. Point predictions are obtained by averaging `mupred` draws, or
#' by PSIS-weighted averaging when `log_weights` is provided.
#'
#' @inheritParams mae
#'
#' @returns A named list with elements `estimate` (mean squared error), `se`
#'   (standard error), and `pointwise` (pointwise squared errors).
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' mse(y, mupred)
#' @export
mse <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {  
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(mupred = mupred, log_weights = log_weights),
      fun_name = "mse"
    )
    sqe_i <- pointwise
  } else {
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
        n_draws = nrow(mupred),
        n_obs = ncol(mupred)
      ))
      sqe_i <- (y - colSums(weights * mupred))^2
    }
  }

  res <- list(
    estimate = mean(sqe_i),
    se = sqrt(var(sqe_i) / length(sqe_i)),
    pointwise = sqe_i
  )
  .measure_output(res, revert_sign)
}

#' Root mean squared error (`rmse`)
#'
#' Computes RMSE as the square root of MSE and propagates uncertainty via a
#' first-order delta-method approximation.
#'
#' @inheritParams mae
#'
#' @returns A named list with elements `estimate` (root mean squared error),
#'   `se` (standard error), and `pointwise` (pointwise squared errors).
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' rmse(y, mupred)
#' @export
rmse <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  mse_res <- mse(
    y = y, mupred = mupred, log_weights = log_weights, 
    pointwise = pointwise
  )
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
  .measure_output(res, revert_sign)
}

#' Predictive R-squared (`r2`)
#'
#' Computes predictive R-squared as one minus the ratio of prediction MSE to
#' the empirical variance of `y`. The standard error is computed with a
#' first-order delta-method approximation.
#'
#' @inheritParams mae
#'
#' @returns A named list with elements `estimate` (predictive R2), `se`
#'   (standard error), and `pointwise` (pointwise squared errors).
#'
#' @examples
#' y <- c(1, 2, 3)
#' mupred <- matrix(c(0.9, 2.1, 2.8, 1.2, 1.9, 3.1), nrow = 2)
#' r2(y, mupred)
#' @export
r2 <- function(
  y, mupred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  .validate_numeric_vector(y, arg = "y")
  if (var(y) == 0) {
    cli::cli_abort(
      "{.fn r2} is undefined when {.arg y} has zero variance."
    )
  }
  
  mse_res <- mse(
    y = y, 
    mupred = mupred, 
    log_weights = log_weights,
    pointwise = pointwise
  )
  mse_hat <- mse_res$estimate
  sqe_i <- mse_res$pointwise
  n_obs <- length(sqe_i)
  
  mse_y_i <- (y - mean(y))^2
  mse_y_hat <- mean(mse_y_i)
   
  var_mse_hat <- mse_res$se^2     
  cov_mse_msey <- stats::cov(sqe_i, mse_y_i) / n_obs              
  var_mse_y_hat <- var(mse_y_i) / n_obs 
  
  t1 <- var_mse_hat
  t2 <- -2 * (mse_hat / mse_y_hat) * cov_mse_msey
  t3 <- (mse_hat^2 / mse_y_hat^2) * var_mse_y_hat
  se_r2 <- sqrt(t1 + t2 + t3) * (1 / mse_y_hat)
  
  est_r2 <- 1 - mse_hat / mse_y_hat
  
  res <- list(
    estimate = est_r2,
    se = se_r2,
    pointwise = sqe_i
  )
  .measure_output(res, revert_sign)
}

# crps ----------------------------------------

#' Continuously ranked probability score
#'
#' The `crps()` and `scrps()` functions and their `loo_*()` counterparts can be
#' used to compute the continuously ranked probability score (CRPS) and scaled
#' CRPS (SCRPS) (as defined by Bolin and Wallin, 2023). CRPS is a proper scoring rule, and
#' strictly proper when the first moment of the predictive distribution is
#' finite. Both can be expressed in terms of samples form the predictive
#' distribution. See, for example, a paper by Gneiting and Raftery (2007)
#' for a comprehensive discussion on CRPS.
#'
#' To compute (S)CRPS, the user needs to provide two sets of draws, `x` and
#' `x2`, from the predictive distribution. This is due to the fact that formulas
#' used to compute CRPS involve an expectation of the absolute difference of `x`
#' and `x2`, both having the same distribution. See the `permutations` argument,
#' as well as Gneiting and Raftery (2007) for details.
#'
#' @export
#' @param x A `S` by `N` matrix (draws by observations), or a vector of length
#'   `S` when only single observation is provided in `y`.
#' @param x2 Independent draws from the same distribution as draws in `x`.
#'   Should be of the identical dimension.
#' @param y A vector of observations or a single value.
#' @param permutations An integer, with default value of 1,  specifying how many
#'   times the expected value of  |X - X'| (`|x - x2|`) is computed. The row
#'   order of `x2` is shuffled as elements `x` and `x2` are typically drawn
#'   given the same values of parameters. This happens, e.g., when one calls
#'   `posterior_predict()` twice for a fitted \pkg{rstanarm} or \pkg{brms}
#'   model. Generating more permutations is expected to decrease the variance of
#'   the computed expected value.
#' @param ... Passed on to [E_loo()] in the `loo_*()` version of these
#'   functions.
#'
#' @return A list containing two elements: `estimates` and `pointwise`.
#'   The former reports estimator and standard error and latter the pointwise
#'   values. Following Bolin & Wallin (2023), a larger value is better.
#'
#' @examples
#' \dontrun{
#' # An example using rstanarm
#' library(rstanarm)
#' data("kidiq")
#' fit <- stan_glm(kid_score ~ mom_hs + mom_iq, data = kidiq)
#' ypred1 <- posterior_predict(fit)
#' ypred2 <- posterior_predict(fit)
#' crps(ypred1, ypred2, y = fit$y)
#' loo_crps(ypred1, ypred2, y = fit$y, log_lik = log_lik(fit))
#' }
#'
#' @references
#' Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
#' proper scoring rules. Statistical Science, 38(1):140-159.
#'
#' Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules,
#' Prediction, and Estimation. Journal of the American Statistical Association,
#' 102(477), 359–378.
crps <- function(x, ...) {
  UseMethod("crps")
}


#' @rdname crps
#' @export
scrps <- function(x, ...) {
  UseMethod("scrps")
}

#' @rdname crps
#' @export
loo_crps <- function(x, ...) {
  UseMethod("loo_crps")
}

#' @rdname crps
#' @export
loo_scrps <- function(x, ...) {
  UseMethod("loo_scrps")
}


#' @rdname crps
#' @export
crps.matrix <- function(x, x2, y, ..., permutations = 1) {
  validate_crps_input(x, x2, y)
  repeats <- replicate(permutations, EXX_compute(x, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x, 2, y)))
  crps_output(.crps_fun(EXX, EXy))
}


#' Method for a single data point
#' @rdname crps
#' @export
crps.numeric <- function(x, x2, y, ..., permutations = 1) {
  stopifnot(length(x) == length(x2),
            length(y) == 1)
  crps.matrix(as.matrix(x), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
#' @param log_lik A log-likelihood matrix the same size as `x`.
#' @param r_eff An optional vector of relative effective sample size estimates
#'   containing one element per observation. See [psis()] for details.
#' @param cores The number of cores to use for parallelization of `[psis()]`.
#'   See [psis()] for details.
loo_crps.matrix <-
  function(x,
           x2,
           y,
           log_lik,
           ...,
           permutations = 1,
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
  validate_crps_input(x, x2, y, log_lik)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x, x2, log_lik, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-log_lik, r_eff = r_eff, cores = cores)
  EXy <- E_loo(abs(sweep(x, 2, y)), psis_obj, log_ratios = -log_lik, ...)$value
  crps_output(.crps_fun(EXX, EXy))
}


#' @rdname crps
#' @export
scrps.matrix <- function(x, x2, y, ..., permutations = 1) {
  validate_crps_input(x, x2, y)
  repeats <- replicate(permutations, EXX_compute(x, x2), simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  EXy <- colMeans(abs(sweep(x, 2, y)))
  crps_output(.crps_fun(EXX, EXy, scale = TRUE))
}

#' @rdname crps
#' @export
scrps.numeric <- function(x, x2, y, ..., permutations = 1) {
  stopifnot(length(x) == length(x2),
            length(y) == 1)
  scrps.matrix(as.matrix(x), as.matrix(x2), y, permutations)
}


#' @rdname crps
#' @export
loo_scrps.matrix <-
  function(
    x,
    x2,
    y,
    log_lik,
    ...,
    permutations = 1,
    r_eff = 1,
    cores = getOption("mc.cores", 1)) {
  validate_crps_input(x, x2, y, log_lik)
  repeats <- replicate(permutations,
                       EXX_loo_compute(x, x2, log_lik, r_eff = r_eff, ...),
                       simplify = F)
  EXX <- Reduce(`+`, repeats) / permutations
  psis_obj <- psis(-log_lik, r_eff = r_eff, cores = cores)
  EXy <- E_loo(abs(sweep(x, 2, y)), psis_obj, log_ratios = -log_lik, ...)$value
  crps_output(.crps_fun(EXX, EXy, scale = TRUE))
}

# ------------ Internals ----------------


EXX_compute <- function(x, x2) {
  S <- nrow(x)
  colMeans(abs(x - x2[sample(1:S),]))
}


EXX_loo_compute <- function(x, x2, log_lik, r_eff = 1, ...) {
  S <- nrow(x)
  shuffle <- sample (1:S)
  x2 <- x2[shuffle,]
  log_lik2 <- log_lik[shuffle,]
  psis_obj_joint <- psis(-log_lik - log_lik2 , r_eff = r_eff)
  E_loo(abs(x - x2), psis_obj_joint, log_ratios = -log_lik - log_lik2, ...)$value
}


#' Function to compute crps and scrps
#' @noRd
.crps_fun <- function(EXX, EXy, scale = FALSE) {
  if (scale) return(-EXy/EXX - 0.5 * log(EXX))
  0.5 * EXX - EXy
}

#' Compute output data for crps functions
#' @noRd
crps_output <- function(crps_pw) {
  n <- length(crps_pw)
  out <- list()
  out$estimates <- c(mean(crps_pw), sd(crps_pw) / sqrt(n))
  names(out$estimates) <- c('Estimate', 'SE')
  out$pointwise <- crps_pw
  out
}

#' Validate input of CRPS functions
#'
#' Check that predictive draws and observed data are of compatible shape
#' @noRd
validate_crps_input <- function(x, x2, y, log_lik = NULL) {
  stopifnot(is.numeric(x),
            is.numeric(x2),
            is.numeric(y),
            identical(dim(x), dim(x2)),
            ncol(x) == length(y),
            ifelse(is.null(log_lik), TRUE, identical(dim(log_lik), dim(x)))
            )
}


#' Continuous ranked probability score (`crps`)
#'
#' Computes CRPS for continuous outcomes from posterior predictive draws using
#' weighted empirical distribution formulas. With `scale = TRUE`, returns the
#' scaled variant used by [scrps()].
#'
#' @param y A numeric vector of observed outcomes.
#' @param ypred A numeric matrix of posterior predictive draws with dimensions
#'   draws x observations.
#' @param log_weights Optional numeric matrix of unnormalized log-weights
#'   with dimensions draws x observations.
#' @param pointwise Optional numeric vector of precomputed pointwise scores.
#'   If provided, `y`, `ypred`, and `log_weights` are ignored.
#' @param scale Logical; if `TRUE`, computes the scaled CRPS.
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate` (mean CRPS score), `se`
#'   (standard error), and `pointwise` (pointwise CRPS scores).
#'
#' @examples
#' y <- c(0.2, -0.1, 1.4)
#' ypred <- matrix(c(0.1, -0.2, 1.3, 0.3, 0.0, 1.7), nrow = 2)
#' crps2(y, ypred)
#' @export
crps2 <- function(
  y, ypred, log_weights = NULL, pointwise = NULL, scale = FALSE, 
  revert_sign = FALSE
) {
  if (!is.null(pointwise)) {
    .inform_ignored_inputs(
      pointwise,
      ignored_args = list(y = y, ypred = ypred, log_weights = log_weights),
      fun_name = "crps"
    )
    crps_i <- pointwise
  } else {
    .validate_numeric_vector(y, arg = "y")
    .validate_numeric_matrix(ypred, arg = "ypred", ncol = length(y))
    n_draws <- nrow(ypred)
    n_obs <- ncol(ypred)
    if (!is.null(log_weights)) {
      weights <- exp(.normalize_and_validate_log_weights(
        log_weights = log_weights,
        n_draws = n_draws,
        n_obs = n_obs
      ))
    } else {
      weights <- rep(1 / n_draws, n_draws)
    }
    
    crps_i <- vapply(seq_len(n_obs), function(i) {
      x <- ypred[, i]
      ord <- order(x)
      x <- x[ord]
      w <- weights[ord]
      
      F_x <- cumsum(w)
      yi <- y[i]

      EXy <- sum(abs(x - yi) * w)  # E|X - y|
      EX <- sum(x * w)             # E[X]
      EXFx <- sum(x * F_x * w)     # E[X * F(X)]

      # E∣X-X'∣ = 4E[X⋅F(X)] - 2E[X]
      EXX <- 4 * EXFx - 2 * EX

      if (scale) {
        # SCRPS(F, y) = -E|X-y| / E|X-X'| - 0.5 * log(E|X-X'|)
        if (EXX <= 0) {
          cli::cli_abort(
            "Cannot compute scaled CRPS because predictive spread is zero for at least one observation."
          )
        }
        -EXy / EXX - 0.5 * log(EXX)
      } else {
        # CRPS = E|X-y| - 0.5*E|X-X'| = E|X-y| + E[X] - 2*E[X*F(X)]
        EXy + EX - 2 * EXFx
      }
    }, numeric(1))
  }
  
  res <- list(
    estimate = mean(crps_i),
    se = sqrt(var(crps_i) / length(crps_i)),
    pointwise = crps_i
  )
  .measure_output(res, revert_sign)
}

#' Scaled continuous ranked probability score (`scrps`)
#'
#' Wrapper around [crps2()] that computes SCRPS, a scale-normalized
#' transformation of CRPS with an additional predictive spread penalty.
#'
#' @inheritParams crps2
#' @param revert_sign Logical; if `TRUE`, multiply the estimate and pointwise
#'   values by -1 before returning.
#'
#' @returns A named list with elements `estimate`, `se`, and `pointwise`.
#'
#' @examples
#' y <- c(0.2, -0.1, 1.4)
#' ypred <- matrix(c(0.1, -0.2, 1.3, 0.3, 0.0, 1.7), nrow = 2)
#' scrps2(y, ypred)
#' @export
scrps2 <- function(
  y, ypred, log_weights = NULL, pointwise = NULL, revert_sign = FALSE
) {
  crps2(
    y = y, ypred = ypred, log_weights = log_weights,
    scale = TRUE, pointwise = pointwise, revert_sign = revert_sign
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
  elpd = elpd,
  mlpd = mlpd,
  mae = mae,
  r2 = r2,
  rmse = rmse,
  mse = mse,
  acc = acc,
  bacc = bacc,
  rps = rps,
  crps2 = crps2,
  srps = srps,
  scrps2 = scrps2
)

#' Supported predictive measure names
#'
#' A character vector of measure names that can be passed to the `measure`
#' argument of [insample_pred_measure()], [loo_pred_measure()],
#' [kfold_pred_measure()], [test_pred_measure()], and [pred_measure()].
#'
#' @export
supported_measures_list <- names(.measure_spec)

.measure_output <- function(res, revert_sign) {
  if (isTRUE(revert_sign)) {
    res$estimate <- -res$estimate
    res$pointwise <- -res$pointwise
  }
  out <- list()
  out$estimates <- c(res$estimate, res$se)
  names(out$estimates) <- c('Estimate', 'SE')
  out$pointwise <- res$pointwise
  out
}