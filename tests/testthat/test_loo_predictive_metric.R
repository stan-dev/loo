LL <- example_loglik_matrix()
chain_id <- rep(1:2, each = dim(LL)[1] / 2)
r_eff <- relative_eff(exp(LL), chain_id)
psis_obj <- psis(-LL, r_eff = r_eff, cores = 2)

set.seed(123)
x <- matrix(rnorm(length(LL)), nrow = nrow(LL), ncol = ncol(LL))
x_prob <- 1 / (1 + exp(-x))
y <- rnorm(ncol(LL))
y_binary <- rbinom(ncol(LL), 1, 0.5)

mae_mean <- loo_predictive_metric(x, y, LL, metric = 'mae', r_eff = r_eff)
mae_quant <- loo_predictive_metric(x, y, LL, metric = 'mae', r_eff = r_eff,
                                  type = 'quantile', probs = 0.9)

rmse_mean <- loo_predictive_metric(x, y, LL, metric = 'rmse', r_eff = r_eff)
rmse_quant <- loo_predictive_metric(x, y, LL, metric = 'rmse', r_eff = r_eff,
                                   type = 'quantile', probs = 0.9)

mse_mean <- loo_predictive_metric(x, y, LL, metric = 'mse', r_eff = r_eff)
mse_quant <- loo_predictive_metric(x, y, LL, metric = 'mse', r_eff = r_eff,
                                  type = 'quantile', probs = 0.9)

acc_mean <- loo_predictive_metric(x_prob, y_binary, LL, metric = 'acc', r_eff = r_eff)
acc_quant <- loo_predictive_metric(x_prob, y_binary, LL, metric = 'acc', r_eff = r_eff,
                                  type = 'quantile', probs = 0.9)

bacc_mean <- loo_predictive_metric(x_prob, y_binary, LL, metric = 'balanced_acc', r_eff = r_eff)
bacc_quant <- loo_predictive_metric(x_prob, y_binary, LL, metric = 'balanced_acc', r_eff = r_eff,
                                  type = 'quantile', probs = 0.9)

test_that('loo_predictive_metric stops with incorrect inputs', {
  expect_error(loo_predictive_metric(as.character(x), y, LL, r_eff = r_eff),
               'no applicable method',
               fixed = TRUE)

  expect_error(loo_predictive_metric(x, as.character(y), LL, r_eff = r_eff),
               'is.numeric(y) is not TRUE',
               fixed = TRUE)

  x_invalid <- matrix(rnorm(9), nrow = 3)
  expect_error(loo_predictive_metric(x_invalid, y, LL, r_eff = r_eff),
               'identical(ncol(x), length(y)) is not TRUE',
               fixed = TRUE)

  x_invalid <- matrix(rnorm(64), nrow = 2)
  expect_error(loo_predictive_metric(x_invalid, y, LL, r_eff = r_eff),
               'identical(dim(x), dim(log_lik)) is not TRUE',
               fixed = TRUE)
})


test_that('loo_predictive_metric return types are correct', {
  # MAE
  expect_type(mae_mean, 'list')
  expect_type(mae_quant, 'list')
  expect_named(mae_mean, c('estimate', 'se'))
  expect_named(mae_quant, c('estimate', 'se'))
  # RMSE
  expect_type(rmse_mean, 'list')
  expect_type(rmse_quant, 'list')
  expect_named(rmse_mean, c('estimate', 'se'))
  expect_named(rmse_quant, c('estimate', 'se'))
  # MSE
  expect_type(mse_mean, 'list')
  expect_type(mse_quant, 'list')
  expect_named(mse_mean, c('estimate', 'se'))
  expect_named(mse_quant, c('estimate', 'se'))
  # Accuracy
  expect_type(acc_mean, 'list')
  expect_type(acc_quant, 'list')
  expect_named(acc_mean, c('estimate', 'se'))
  expect_named(acc_quant, c('estimate', 'se'))
  # Balanced accuracy
  expect_type(bacc_mean, 'list')
  expect_type(bacc_quant, 'list')
  expect_named(bacc_mean, c('estimate', 'se'))
  expect_named(bacc_quant, c('estimate', 'se'))
})

test_that('loo_predictive_metric is equal to reference', {
  expect_equal_to_reference(mae_mean, 'reference-results/loo_predictive_metric_mae_mean.rds')
  expect_equal_to_reference(mae_quant, 'reference-results/loo_predictive_metric_mae_quant.rds')
  expect_equal_to_reference(rmse_mean, 'reference-results/loo_predictive_metric_rmse_mean.rds')
  expect_equal_to_reference(rmse_quant, 'reference-results/loo_predictive_metric_rmse_quant.rds')
  expect_equal_to_reference(mse_mean, 'reference-results/loo_predictive_metric_mse_mean.rds')
  expect_equal_to_reference(mse_quant, 'reference-results/loo_predictive_metric_mse_quant.rds')
  expect_equal_to_reference(acc_mean, 'reference-results/loo_predictive_metric_acc_mean.rds')
  expect_equal_to_reference(acc_quant, 'reference-results/loo_predictive_metric_acc_quant.rds')
  expect_equal_to_reference(bacc_mean, 'reference-results/loo_predictive_metric_bacc_mean.rds')
  expect_equal_to_reference(bacc_quant, 'reference-results/loo_predictive_metric_bacc_quant.rds')
})

test_that('MAE computation is correct', {
  expect_equal(
    .mae(rep(0.5, 5), rep(1, 5))$estimate,
    0.5)
  expect_equal(
    .mae(rep(0.5, 5), rep(1, 5))$se,
    0.0)
  expect_error(
    .mae(rep(0.5, 5), rep(1, 3)),
    'length(y) == length(yhat) is not TRUE',
    fixed = TRUE)
})

test_that('MSE computation is correct', {
  expect_equal(
    .mse(rep(0.5, 5), rep(1, 5))$estimate,
    0.25)
  expect_equal(
    .mse(rep(0.5, 5), rep(1, 5))$se,
    0.0)
  expect_error(
    .mse(rep(0.5, 5), rep(1, 3)),
    'length(y) == length(yhat) is not TRUE',
    fixed = TRUE)
})

test_that('RMSE computation is correct', {
  expect_equal(
    .rmse(rep(0.5, 5), rep(1, 5))$estimate,
    sqrt(0.25))
  expect_equal(
    .mse(rep(0.5, 5), rep(1, 5))$se,
    0.0)
  expect_error(
    .mse(rep(0.5, 5), rep(1, 3)),
    'length(y) == length(yhat) is not TRUE',
    fixed = TRUE)
})

test_that('Accuracy computation is correct', {
  expect_equal(
    .accuracy(c(0, 0, 0, 1, 1, 1), c(0.2, 0.2, 0.2, 0.7, 0.7, 0.7))$estimate,
    1.0
  )
  expect_error(
    .accuracy(c(1, 0), c(0.5)),
    'length(y) == length(yhat) is not TRUE',
    fixed = TRUE)
  expect_error(
    .accuracy(c(2, 1), c(0.5, 0.5)),
    'all(y <= 1 & y >= 0) is not TRUE',
    fixed = TRUE
  )
  expect_error(
    .accuracy(c(1, 0), c(1.1, 0.5)),
    'all(yhat <= 1 & yhat >= 0) is not TRUE',
    fixed = TRUE
  )
})

test_that('Balanced accuracy computation is correct', {
  expect_equal(
    .balanced_accuracy(c(0, 0, 1, 1, 1, 1), c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9))$estimate,
    0.5
  )
  expect_error(
    .balanced_accuracy(c(1, 0), c(0.5)),
    'length(y) == length(yhat) is not TRUE',
    fixed = TRUE)
  expect_error(
    .balanced_accuracy(c(2, 1), c(0.5, 0.5)),
    'all(y <= 1 & y >= 0) is not TRUE',
    fixed = TRUE
  )
  expect_error(
    .balanced_accuracy(c(1, 0), c(1.1, 0.5)),
    'all(yhat <= 1 & yhat >= 0) is not TRUE',
    fixed = TRUE
  )
})
