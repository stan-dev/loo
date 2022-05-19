#' Estimate leave-one-out predictive errors.
#'
#' The `loo_predictive_error()` function computes estimates of leave-one-out
#' predictive errors given a set of predictions and observations. Curreantly
#' supported error metrics are mean absolute error, mean squared error and root
#' mean squared error for continuous predictions and accuracy and balanced
#' accuracy for binary classification. Predictions are passed on to the `E_loo()`
#' function, so this function assumes that the PSIS approximation is working well.
#'
#' @param x A numeric matrix of predictions.
#' @param y A numeric vector of observations. Length should be equal to the
#'     number of rows in `x`.
#' @param ll A matrix of pointwise log-likelihoods. Should be of same dimension
#'     as `x`.
#' @param pred_error The type of predictive error to be used. Currently
#'     supported options are `"mae"`, `"rmse"` and `"mse"` for regression and
#'     for binary classification `"acc"` and `"balanced_acc"`.
#'     \describe{
#'       \item{`"mae"`}{
#'          Mean absolute error.
#'       }
#'       \item{`"mse"`}{
#'          Mean squared error.
#'       }
#'       \item{`"rmse"`}{
#'          Root mean squared error, given by as the square root of `MSE`.
#'       }
#'       \item{`"acc"`}{
#'          The proportion of predictions indicating the correct outcome.
#'       }
#'       \item{`"balanced_acc"`}{
#'          Balanced accuracy is given by the average of true positive and true
#'          negative rates.
#'       }
#'     }
#' @param r_eff A Vector of relative effective sample size estimates containing
#'     one element per observation. See [psis()] for more details.
#' @param ... Additional arguments passed on to [E_loo()]
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`estimate`}{
#'   Estimate of the given predictive error measure.
#'   }
#'  \item{`se`}{
#'   Standard error of the estimate.
#'   }
#'  }
#' @export
#'
#' @examples
loo_predictive_error <- function(x, ...) {
  UseMethod("loo_predictive_error")
}

#' @rdname loo_predictive_error
#' @export
loo_predictive_error.default <-
  function(x,
           y,
           ll,
           pred_error = c("mae", "rmse", "mse", "acc", "balanced_acc"),
           r_eff = NULL,
           ...) {
    stopifnot(
      is.numeric(x),
      is.numeric(y),
      identical(ncol(x), length(y)),
      identical(dim(x), dim(ll))
    )
    pred_error <- match.arg(pred_error)
    psis_object <- psis(-ll, r_eff = r_eff)
    pred_loo <- E_loo(x, psis_object = psis_object, log_ratios = -ll, ...)$value

    predictive_error_fun <- .loo_predictive_error_fun(pred_error)

    return(predictive_error_fun(y, pred_loo))
  }


# ----------------------------- Internals -----------------------------

#' Select predictive error function based on user's `pred_error` argument
#'
#' @noRd
#' @param pred_error The predictive error used.
#' @return The function used to compute predictive error specified by `pred_error`.
.loo_predictive_error_fun <- function(pred_error) {
  switch(
    pred_error,
    'mae' = .mae,
    'rmse' = .rmse,
    'mse' = .mse,
    'acc' = .accuracy,
    'balanced_acc' = .balanced_accuracy
  )
}

#' Mean absolute error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.mae <-function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  n <- length(y)
  e <- abs(y - yhat)
  return(list(estimate = mean(e), se = sd(e) / sqrt(n)))
}

#' Mean squared error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.mse <-function(y, yhat) {
  stopifnot(length(y) == length(yhat))
  n <- length(y)
  e <- (y - yhat)^2
  return(list(estimate = mean(e), se = sd(e) / sqrt(n)))
}

#' Root mean squared error
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.rmse <-function(y, yhat) {
  est <- .mse(y, yhat)
  mean_mse <- est$estimate
  var_mse <- est$se^2
  var_rmse <- var_mse / mean_mse / 4 # Comes from the first order Taylor approx.
  return(list(estimate = sqrt(mean_mse), se = sqrt(var_rmse)))
}

#' Classification accuracy
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.accuracy <- function(y, yhat) {
  stopifnot(length(y) == length(yhat),
            all(y <= 1 & y >= 0),
            all(yhat <= 1 & yhat >= 0))
  n <- length(y)
  yhat <- as.integer(yhat > 0.5)
  acc <- as.integer(yhat == y)
  est <- mean(acc)
  return(list(estimate = est, se = sqrt(est * (1-est) / n) ))
}

#' Balanced classification accuracy
#'
#' @noRd
#' @param y A vector of observed values
#' @param yhat A vector of predictions
.balanced_accuracy <- function(y, yhat) {
  stopifnot(length(y) == length(yhat),
            all(y <= 1 & y >= 0),
            all(yhat <= 1 & yhat >= 0))
  n <- length(y)
  yhat <- as.integer(yhat > 0.5)
  mask <- y == 0

  tn <- mean(yhat[mask] == y[mask]) # True negatives
  tp <- mean(yhat[!mask] == y[!mask]) # True positives

  bls_acc <- (tp + tn) / 2
  # This approximation has quite large bias for small samples
  bls_acc_var <- (tp * (1 - tp) + tn * (1 - tn)) / 4
  return(list(estimate = bls_acc, se = sqrt(bls_acc_var / n)))
}

