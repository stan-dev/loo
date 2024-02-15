#' Estimate leave-one-out predictive performance..
#'
#' The `loo_predictive_metric()` function computes estimates of leave-one-out
#' predictive metrics given a set of predictions and observations. Currently
#' supported metrics are mean absolute error, mean squared error and root mean
#' squared error for continuous predictions and accuracy and balanced accuracy
#' for binary classification. Predictions are passed on to the [E_loo()]
#' function, so this function assumes that the PSIS approximation is working
#' well.
#'
#' @param x A numeric matrix of predictions.
#' @param y A numeric vector of observations. Length should be equal to the
#'     number of rows in `x`.
#' @param log_lik A matrix of pointwise log-likelihoods. Should be of same
#'     dimension as `x`.
#' @param metric The type of predictive metric to be used. Currently
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
#' @param cores The number of cores to use for parallelization of `[psis()]`.
#'   See [psis()] for details.
#' @param ... Additional arguments passed on to [E_loo()]
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`estimate`}{
#'   Estimate of the given metric.
#'   }
#'  \item{`se`}{
#'   Standard error of the estimate.
#'   }
#'  }
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#' # Use rstanarm package to quickly fit a model and get both a log-likelihood
#' # matrix and draws from the posterior predictive distribution
#' library("rstanarm")
#'
#' # data from help("lm")
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' d <- data.frame(
#'   weight = c(ctl, trt),
#'   group = gl(2, 10, 20, labels = c("Ctl","Trt"))
#' )
#' fit <- stan_glm(weight ~ group, data = d, refresh = 0)
#' ll <- log_lik(fit)
#' r_eff <- relative_eff(exp(-ll), chain_id = rep(1:4, each = 1000))
#'
#' mu_pred <- posterior_epred(fit)
#' # Leave-one-out mean absolute error of predictions
#' mae <- loo_predictive_metric(x = mu_pred, y = d$weight, log_lik = ll,
#'                             pred_error = 'mae', r_eff = r_eff)
#' # Leave-one-out 90%-quantile of mean absolute error
#' mae_90q <- loo_predictive_metric(x = mu_pred, y = d$weight, log_lik = ll,
#'                                 pred_error = 'mae', r_eff = r_eff,
#'                                 type = 'quantile', probs = 0.9)
#' }
#' }
loo_predictive_metric <- function(x, ...) {
  UseMethod("loo_predictive_metric")
}

#' @rdname loo_predictive_metric
#' @export
loo_predictive_metric.matrix <-
  function(x,
           y,
           log_lik,
           ...,
           metric = c("mae", "rmse", "mse", "acc", "balanced_acc"),
           r_eff = 1,
           cores = getOption("mc.cores", 1)) {
    stopifnot(
      is.numeric(x),
      is.numeric(y),
      identical(ncol(x), length(y)),
      identical(dim(x), dim(log_lik))
    )
    metric <- match.arg(metric)
    psis_object <- psis(-log_lik, r_eff = r_eff, cores = cores)
    pred_loo <- E_loo(x,
                      psis_object = psis_object,
                      log_ratios = -log_lik,
                      ...)$value

    predictive_metric_fun <- .loo_predictive_metric_fun(metric)

    predictive_metric_fun(y, pred_loo)
  }


# ----------------------------- Internals -----------------------------

#' Select predictive metric function based on user's `metric` argument
#'
#' @noRd
#' @param metric The metric used.
#' @return The function used to compute predictive error or accuracy specified
#' by the argument `metric`.
.loo_predictive_metric_fun <- function(metric) {
  switch(
    metric,
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
  list(estimate = mean(e), se = sd(e) / sqrt(n))
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
  list(estimate = mean(e), se = sd(e) / sqrt(n))
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
  list(estimate = est, se = sqrt(est * (1-est) / n) )
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
  list(estimate = bls_acc, se = sqrt(bls_acc_var / n))
}

