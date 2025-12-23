# Estimate leave-one-out predictive performance..

The `loo_predictive_metric()` function computes estimates of
leave-one-out predictive metrics given a set of predictions and
observations. Currently supported metrics are mean absolute error, mean
squared error and root mean squared error for continuous predictions and
accuracy and balanced accuracy for binary classification. Predictions
are passed on to the
[`E_loo()`](https://mc-stan.org/loo/reference/E_loo.md) function, so
this function assumes that the PSIS approximation is working well.

## Usage

``` r
loo_predictive_metric(x, ...)

# S3 method for class 'matrix'
loo_predictive_metric(
  x,
  y,
  log_lik,
  ...,
  metric = c("mae", "rmse", "mse", "acc", "balanced_acc"),
  r_eff = 1,
  cores = getOption("mc.cores", 1)
)
```

## Arguments

- x:

  A numeric matrix of predictions.

- ...:

  Additional arguments passed on to
  [`E_loo()`](https://mc-stan.org/loo/reference/E_loo.md)

- y:

  A numeric vector of observations. Length should be equal to the number
  of rows in `x`.

- log_lik:

  A matrix of pointwise log-likelihoods. Should be of same dimension as
  `x`.

- metric:

  The type of predictive metric to be used. Currently supported options
  are `"mae"`, `"rmse"` and `"mse"` for regression and for binary
  classification `"acc"` and `"balanced_acc"`.

  `"mae"`

  :   Mean absolute error.

  `"mse"`

  :   Mean squared error.

  `"rmse"`

  :   Root mean squared error, given by as the square root of `MSE`.

  `"acc"`

  :   The proportion of predictions indicating the correct outcome.

  `"balanced_acc"`

  :   Balanced accuracy is given by the average of true positive and
      true negative rates.

- r_eff:

  A Vector of relative effective sample size estimates containing one
  element per observation. See
  [`psis()`](https://mc-stan.org/loo/reference/psis.md) for more
  details.

- cores:

  The number of cores to use for parallelization of `[psis()]`. See
  [`psis()`](https://mc-stan.org/loo/reference/psis.md) for details.

## Value

A list with the following components:

- `estimate`:

  Estimate of the given metric.

- `se`:

  Standard error of the estimate.

## Examples

``` r
# \donttest{
if (requireNamespace("rstanarm", quietly = TRUE)) {
# Use rstanarm package to quickly fit a model and get both a log-likelihood
# matrix and draws from the posterior predictive distribution
library("rstanarm")

# data from help("lm")
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
d <- data.frame(
  weight = c(ctl, trt),
  group = gl(2, 10, 20, labels = c("Ctl","Trt"))
)
fit <- stan_glm(weight ~ group, data = d, refresh = 0)
ll <- log_lik(fit)
r_eff <- relative_eff(exp(-ll), chain_id = rep(1:4, each = 1000))

mu_pred <- posterior_epred(fit)
# Leave-one-out mean absolute error of predictions
mae <- loo_predictive_metric(x = mu_pred, y = d$weight, log_lik = ll,
                            pred_error = 'mae', r_eff = r_eff)
# Leave-one-out 90%-quantile of mean absolute error
mae_90q <- loo_predictive_metric(x = mu_pred, y = d$weight, log_lik = ll,
                                pred_error = 'mae', r_eff = r_eff,
                                type = 'quantile', probs = 0.9)
}
# }
```
