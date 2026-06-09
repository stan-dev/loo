# load test data --------------------------------------
res_roaches <- readRDS("tests/testthat/data-for-tests/res_roaches.Rds")
res_sleep <- readRDS("data-for-tests/res_sleep.Rds")
res_binom <- readRDS("data-for-tests/res_binomial.Rds")
res_binary <- readRDS("data-for-tests/res_binary.Rds")
res_cat <- readRDS("data-for-tests/res_penguins.Rds")

# ptw_log_pred_density ------------------------
testthat::test_that("ptw_log_pred_density() works as expected", {
  res <- ptw_log_pred_density(ylp = res_roaches$ylp, psis_log_weights = NULL)

  expect_equal(length(res), dim(res_roaches$ylp)[2])
  expect_equal(res, matrixStats::colLogSumExps(res_roaches$ylp) - log(dim(res_roaches$ylp)[1]))
})

testthat::test_that("ptw_log_pred_density() with psis_log_weights works as expected", {
  norm_log_weights <- .normalize_log_weights(res_roaches$log_weights)
  res <- ptw_log_pred_density(
    ylp = res_roaches$ylp,
    psis_log_weights = norm_log_weights
  )

  expect_equal(length(res), dim(res_roaches$ylp)[2])
  expect_equal(res, matrixStats::colLogSumExps(res_roaches$ylp + norm_log_weights)
  )
})

testthat::test_that("ptw_log_pred_density() returns error when weights are not normalized", {
  expect_error(
    ptw_log_pred_density(
      ylp = res_roaches$ylp,
      psis_log_weights = res_roaches$log_weights
    ),
    regexp = "Range of current column sums"
  )
})

# elpd() -----------------------------------

testthat::test_that("elpd() works as expected", {
  res <- elpd(ylp = res_roaches$ylp, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
  
  expect_snapshot_output(
    elpd(ylp = res_roaches$ylp, log_weights = NULL)
  )
})

testthat::test_that("elpd() with unnormalized log-weights works as expected", {
  log_weights <- res_roaches$log_weights
  res <- elpd(ylp = res_roaches$ylp, log_weights = log_weights)

  expect_false(all(.normalize_log_weights(log_weights) == log_weights))
  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
})

testthat::test_that("elpd() with normalized log-weights works as expected", {
  res <- elpd(ylp = res_roaches$ylp, log_weights = res_roaches$predperf_loo$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
})

# ic() -----------------------------------

testthat::test_that("ic() works as expected", {
  res <- ic(ylp = res_roaches$ylp)
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_snapshot_output(ic(ylp = res_roaches$ylp))
})

# mlpd() -----------------------------------

testthat::test_that("mlpd() works as expected", {
  res <- mlpd(ylp = res_roaches$ylp, log_weights = NULL)
  res_elpd <- elpd(ylp = res_roaches$ylp, log_weights = NULL)
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), n_obs)
  expect_equal(res$estimates, res_elpd$estimates / n_obs)
  expect_equal(res$pointwise, res_elpd$pointwise)
  
  expect_snapshot_output(mlpd(ylp = res_roaches$ylp, log_weights = NULL))
})

testthat::test_that("mlpd() with pointwise works as expected", {
  res_elpd <- elpd(ylp = res_roaches$ylp, log_weights = NULL)
  res <- mlpd(ylp = NULL, pointwise = res_elpd$pointwise)
  
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), n_obs)
  expect_equal(res$estimates, res_elpd$estimates / n_obs)
  expect_equal(res$pointwise, res_elpd$pointwise)
  
  expect_snapshot_output(mlpd(ylp = NULL, pointwise = res_elpd$pointwise))
})

# rps() -------------------------------------

testthat::test_that("rps() with ordered categorial data works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binom$y))
  
  expect_snapshot_output(rps(y = res_binom$y, ypred = res_binom$ypred))
})

testthat::test_that("rps() scaled version with categorical data works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    scaled = TRUE
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binom$y))
  expect_true(all(res$pointwise < 0))
  
  expect_snapshot_output(srps(y = res_binom$y, ypred = res_binom$ypred))
})

testthat::test_that("rps() for categorical data with log-weights works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    log_weights = res_binom$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binary$y))
  expect_true(all(res$pointwise >= 0))
})

testthat::test_that("rps() with continuous data works as expected", {
  res <- rps(res_sleep$y, res_sleep$ypred)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
  
  expect_snapshot_output(rps(res_sleep$y, res_sleep$ypred))
})

testthat::test_that("rps() with continuous data and log-weights works as expected", {
  res <- rps(
    res_sleep$y, 
    res_sleep$ypred, 
    log_weights = res_sleep$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
})

testthat::test_that("rps() with continuous data and scaled version works as expected", {
  res <- rps(res_sleep$y, res_sleep$ypred, scaled = TRUE)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
  expect_true(all(res$pointwise < 0))
})


# brier() ---------------------------------------

testthat::test_that("brier() works as expected", {
  res_brier <- brier(y = res_binary$y, ypred = res_binary$ypred, log_weights = NULL)

  expect_equal(names(res_brier), c("estimates", "pointwise"))
  expect_equal(length(res_brier$estimates), 2)
  expect_equal(length(res_brier$pointwise), length(res_binary$y))
  expect_true(all(res_brier$pointwise >= 0 & res_brier$pointwise <= 1))
  
  expect_snapshot_output(brier(y = res_binary$y, ypred = res_binary$ypred))
})

testthat::test_that("brier() with log-weights works as expected", {
  res_brier <- brier(
    y = res_binary$y,
    ypred = res_binary$ypred,
    log_weights = res_binary$log_weights
  )

  expect_equal(names(res_brier), c("estimates", "pointwise"))
  expect_equal(length(res_brier$estimates), 2)
  expect_equal(length(res_brier$pointwise), length(res_binary$y))
  expect_true(all(res_brier$pointwise >= 0 & res_brier$pointwise <= 1))

  res_brier2 <- brier(
    y = res_binary$y,
    ypred = res_binary$ypred,
    log_weights = .normalize_log_weights(res_binary$log_weights)
  )
  expect_equal(res_brier2$pointwise, res_brier$pointwise)
})

testthat::test_that("brier() errors when y not binary", {
  expect_error(
    brier(
      y = res_binom$y,
      ypred = res_binary$ypred,
      log_weights = res_binary$log_weights
    ),
    regexp = "The brier score expects binary data 'y'."
  )
})

# mae() ------------------------------------------------
testthat::test_that("mae() works as expected", {
  res <- mae(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  
  expect_snapshot_output(mae(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("mae() with log_weights works as expected", {
  res <- mae(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = res_roaches$loo1$psis_object$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
})

# rmse() / mse() -----------------------------------------
testthat::test_that("mse() and rmse() work as expected", {
  res_mse <- mse(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  res_rmse <- rmse(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  expect_equal(names(res_mse), c("estimates", "pointwise"))
  expect_equal(length(res_mse$estimates), 2)
  expect_equal(length(res_mse$pointwise), length(res_roaches$y))

  expect_equal(sqrt(abs(res_mse$estimates[1])), res_rmse$estimates[1])
  expect_equal(res_mse$estimates[2]/(2*sqrt(abs(res_mse$estimates[1]))), 
  res_rmse$estimates[2])
  
  expect_snapshot_output(mse(y = res_roaches$y, mupred = res_roaches$mupred))
  expect_snapshot_output(rmse(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("rmse() works with se=0", {
  mupred0 <- t(replicate(4000, res_roaches$y))

  res_rmse <- rmse(y = res_roaches$y, mupred = mupred0,
    log_weights = NULL)

  expect_equal(names(res_rmse), c("estimates", "pointwise"))
  expect_equal(length(res_rmse$estimates), 2)
  expect_equal(length(res_rmse$pointwise), length(res_roaches$y))
  expect_equal(unname(res_rmse$estimates[2]), 0)
})

# r2() ---------------------------------------------------
testthat::test_that("r2() works as expected", {
  res <- r2(y = res_roaches$y, mupred = res_roaches$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  expect_true(all(res$estimates[1] >= 0 & res$estimates[1] <= 1))
  
  expect_snapshot_output(r2(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("r2() with log_weights works as expected", {
  res <- r2(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = res_roaches$loo1$psis_object$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  expect_true(all(res$estimates[1] >= 0 & res$estimates[1] <= 1))
})

# acc() / bacc() -----------------------------------------------------
testthat::skip_if_not_installed("brms")

testthat::test_that("acc() works as expected", {
  res <- acc(y = as.integer(res_cat$y), mupred = res_cat$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(all(res$pointwise >= 0 & res$pointwise <= 1))
  
  expect_snapshot_output(acc(y = as.integer(res_cat$y), mupred = res_cat$mupred))
})

testthat::test_that("acc() with log-weights works as expected", {
  res <- acc(
    y = as.integer(res_cat$y),
    mupred = res_cat$mupred,
    log_weights = res_cat$loo$psis_object$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})

testthat::test_that("bacc() works as expected", {
  res <- bacc(y = as.integer(res_cat$y), mupred = res_cat$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
  
  expect_snapshot_output(bacc(y = as.integer(res_cat$y), mupred = res_cat$mupred))
})

testthat::test_that("bacc() with log-weights works as expected", {
  res <- bacc(
    y = as.integer(res_cat$y),
    mupred = res_cat$mupred,
    log_weights = res_cat$loo$psis_object$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})
