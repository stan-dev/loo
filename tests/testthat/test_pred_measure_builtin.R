# load test data --------------------------------------
path <- ""#"tests/testthat/"
res_roaches <- readRDS(paste0(path, "data-for-tests/test_data_roaches.Rds"))
res_sleep <- readRDS(paste0(path, "data-for-tests/test_data_sleep.Rds"))
res_binom <- readRDS(paste0(path, "data-for-tests/test_data_binomial.Rds"))
res_binary <- readRDS(paste0(path, "data-for-tests/test_data_binary.Rds"))
res_cat <- readRDS(paste0(path, "data-for-tests/test_data_penguins.Rds"))

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

# measure_elpd() -----------------------------------

testthat::test_that("measure_elpd() works as expected", {
  res <- measure_elpd(ylp = res_roaches$ylp, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
  
  expect_snapshot_output(
    measure_elpd(ylp = res_roaches$ylp, log_weights = NULL)
  )
})

testthat::test_that("measure_elpd() with unnormalized log-weights works as expected", {
  log_weights <- res_roaches$log_weights
  res <- measure_elpd(ylp = res_roaches$ylp, log_weights = log_weights)

  expect_false(all(.normalize_log_weights(log_weights) == log_weights))
  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
})

testthat::test_that("measure_elpd() with normalized log-weights works as expected", {
  res <- measure_elpd(ylp = res_roaches$ylp, log_weights = res_roaches$predperf_loo$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), dim(res_roaches$ylp)[2])
})

# measure_ic() -----------------------------------

testthat::test_that("measure_ic() works as expected", {
  res <- measure_ic(ylp = res_roaches$ylp)
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_snapshot_output(measure_ic(ylp = res_roaches$ylp))
})

# measure_mlpd() -----------------------------------

testthat::test_that("measure_mlpd() works as expected", {
  res <- measure_mlpd(ylp = res_roaches$ylp, log_weights = NULL)
  res_elpd <- measure_elpd(ylp = res_roaches$ylp, log_weights = NULL)
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), n_obs)
  expect_equal(unname(res$estimates), unname(res_elpd$estimates) / n_obs)
  expect_equal(unname(res$pointwise), unname(res_elpd$pointwise))
  
  expect_snapshot_output(measure_mlpd(ylp = res_roaches$ylp, log_weights = NULL))
})

testthat::test_that("measure_mlpd() with pointwise works as expected", {
  res_elpd <- measure_elpd(ylp = res_roaches$ylp, log_weights = NULL)
  res <- measure_mlpd(ylp = NULL, pointwise = res_elpd$pointwise[ ,"elpd"])
  
  n_obs <- dim(res_roaches$ylp)[2]

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates[1]), 1)
  expect_equal(length(res$estimates[2]), 1)
  expect_equal(length(res$pointwise), n_obs)
  expect_equal(unname(res$estimates), unname(res_elpd$estimates) / n_obs)
  expect_equal(unname(res$pointwise), unname(res_elpd$pointwise))
  
  expect_snapshot_output(measure_mlpd(ylp = NULL, pointwise = res_elpd$pointwise[ ,"elpd"]))
})

# measure_rps() -------------------------------------

testthat::test_that("measure_rps() with ordered categorial data works as expected", {
  res <- measure_rps(
    y = res_binom$y,
    ypred = res_binom$ypred
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binom$y))
  
  expect_snapshot_output(measure_rps(y = res_binom$y, ypred = res_binom$ypred))
})

testthat::test_that("measure_rps() scaled version with categorical data works as expected", {
  res <- measure_rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    scaled = TRUE
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binom$y))
  expect_true(all(res$pointwise < 0))
  
  expect_snapshot_output(measure_srps(y = res_binom$y, ypred = res_binom$ypred))
})

testthat::test_that("measure_rps() for categorical data with log-weights works as expected", {
  res <- measure_rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    log_weights = res_binom$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binary$y))
  expect_true(all(res$pointwise >= 0))
})

testthat::test_that("measure_rps() with continuous data works as expected", {
  res <- measure_rps(res_sleep$y, res_sleep$ypred)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
  
  expect_snapshot_output(measure_rps(res_sleep$y, res_sleep$ypred))
})

testthat::test_that("measure_rps() with continuous data and log-weights works as expected", {
  res <- measure_rps(
    res_sleep$y, 
    res_sleep$ypred, 
    log_weights = res_sleep$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
})

testthat::test_that("measure_rps() with continuous data and scaled version works as expected", {
  res <- measure_rps(res_sleep$y, res_sleep$ypred, scaled = TRUE)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_sleep$y))
  expect_true(all(res$pointwise < 0))
})


# measure_brier() ---------------------------------------

testthat::test_that("measure_brier() works as expected", {
  res_brier <- measure_brier(y = res_binary$y, ypred = res_binary$ypred, log_weights = NULL)

  expect_equal(names(res_brier), c("estimates", "pointwise"))
  expect_equal(length(res_brier$estimates), 2)
  expect_equal(length(res_brier$pointwise), length(res_binary$y))
  expect_true(all(res_brier$pointwise >= 0 & res_brier$pointwise <= 1))
  
  expect_snapshot_output(measure_brier(y = res_binary$y, ypred = res_binary$ypred))
})

testthat::test_that("measure_brier() rejects out-of-range ypred", {
  bad_ypred <- res_binary$ypred
  bad_ypred[1, 1] <- 1.5
  expect_error(
    measure_brier(y = res_binary$y, ypred = bad_ypred),
    regexp = "`ypred` must contain values in \\[0, 1\\]"
  )
})

testthat::test_that("measure_brier() with log-weights works as expected", {
  res_brier <- measure_brier(
    y = res_binary$y,
    ypred = res_binary$ypred,
    log_weights = res_binary$log_weights
  )

  expect_equal(names(res_brier), c("estimates", "pointwise"))
  expect_equal(length(res_brier$estimates), 2)
  expect_equal(length(res_brier$pointwise), length(res_binary$y))
  expect_true(all(res_brier$pointwise >= 0 & res_brier$pointwise <= 1))

  res_brier2 <- measure_brier(
    y = res_binary$y,
    ypred = res_binary$ypred,
    log_weights = .normalize_log_weights(res_binary$log_weights)
  )
  expect_equal(res_brier2$pointwise, res_brier$pointwise)
})

testthat::test_that("measure_brier() errors when y not binary", {
  expect_error(
    measure_brier(
      y = res_binom$y,
      ypred = res_binary$ypred,
      log_weights = res_binary$log_weights
    ),
    regexp = "The brier score expects binary data 'y'."
  )
})

# measure_mae() ------------------------------------------------
testthat::test_that("measure_mae() works as expected", {
  res <- measure_mae(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  
  expect_snapshot_output(measure_mae(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("measure_mae() with log_weights works as expected", {
  res <- measure_mae(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = res_roaches$loo1$psis_object$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
})

# measure_rmse() / measure_mse() -----------------------------------------
testthat::test_that("measure_mse() and measure_rmse() work as expected", {
  res_mse <- measure_mse(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  res_rmse <- measure_rmse(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  expect_equal(names(res_mse), c("estimates", "pointwise"))
  expect_equal(length(res_mse$estimates), 2)
  expect_equal(length(res_mse$pointwise), length(res_roaches$y))

  expect_equal(sqrt(abs(res_mse$estimates[1])), res_rmse$estimates[1])
  expect_equal(res_mse$estimates[2]/(2*sqrt(abs(res_mse$estimates[1]))), 
  res_rmse$estimates[2])
  
  expect_snapshot_output(measure_mse(y = res_roaches$y, mupred = res_roaches$mupred))
  expect_snapshot_output(measure_rmse(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("measure_rmse() works with se=0", {
  mupred0 <- t(replicate(4000, res_roaches$y))

  res_rmse <- measure_rmse(y = res_roaches$y, mupred = mupred0,
    log_weights = NULL)

  expect_equal(names(res_rmse), c("estimates", "pointwise"))
  expect_equal(length(res_rmse$estimates), 2)
  expect_equal(length(res_rmse$pointwise), length(res_roaches$y))
  expect_equal(unname(res_rmse$estimates[2]), 0)
})

# measure_r2() ---------------------------------------------------
testthat::test_that("measure_r2() works as expected", {
  res <- measure_r2(y = res_roaches$y, mupred = res_roaches$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  expect_true(all(res$estimates[1] >= 0 & res$estimates[1] <= 1))
  
  expect_snapshot_output(measure_r2(y = res_roaches$y, mupred = res_roaches$mupred))
})

testthat::test_that("measure_r2() with log_weights works as expected", {
  res <- measure_r2(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = res_roaches$loo1$psis_object$log_weights)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  expect_true(all(res$estimates[1] >= 0 & res$estimates[1] <= 1))
})

# measure_acc() / measure_bacc() -----------------------------------------------------
testthat::skip_if_not_installed("brms")

testthat::test_that("measure_acc() works as expected", {
  res <- measure_acc(y = as.integer(res_cat$y), mupred = res_cat$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(all(res$pointwise >= 0 & res$pointwise <= 1))
  
  expect_snapshot_output(measure_acc(y = as.integer(res_cat$y), mupred = res_cat$mupred))
})

testthat::test_that("measure_acc() rejects out-of-range mupred", {
  bad_mupred <- res_cat$mupred
  bad_mupred[1, 1, 1] <- -0.1
  expect_error(
    measure_acc(y = as.integer(res_cat$y), mupred = bad_mupred),
    regexp = "`mupred` must contain values in \\[0, 1\\]"
  )
})

testthat::test_that("measure_acc() with log-weights works as expected", {
  res <- measure_acc(
    y = as.integer(res_cat$y),
    mupred = res_cat$mupred,
    log_weights = res_cat$loo$psis_object$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})

testthat::test_that("measure_bacc() pointwise contributions sum to estimate", {
  y <- c(1L, 1L, 2L, 2L)
  mupred <- array(
    c(0.8, 0.2, 0.7, 0.3, 0.3, 0.7, 0.2, 0.8),
    dim = c(1, 4, 2)
  )
  res <- measure_bacc(y, mupred)

  classes <- sort(unique(y))
  K <- length(classes)
  weights <- rep(1 / nrow(mupred), nrow(mupred))
  weighted_mupred <- apply(sweep(mupred, 1, weights, `*`), c(2, 3), sum)
  mupred_hat <- apply(weighted_mupred, 1, which.max)
  acc_i <- (mupred_hat == y) * 1L
  n_c <- tabulate(match(y, classes))
  expected_bacc_i <- acc_i / (K * n_c[match(y, classes)])

  expect_equal(as.numeric(res$pointwise), expected_bacc_i)
  expect_equal(sum(as.numeric(res$pointwise)), unname(res$estimates[1, "Estimate"]))
})

testthat::test_that("measure_bacc() works as expected", {
  res <- measure_bacc(y = as.integer(res_cat$y), mupred = res_cat$mupred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
  
  expect_snapshot_output(measure_bacc(y = as.integer(res_cat$y), mupred = res_cat$mupred))
})

testthat::test_that("measure_bacc() rejects out-of-range mupred", {
  bad_mupred <- res_cat$mupred
  bad_mupred[1, 1, 1] <- 1.2
  expect_error(
    measure_bacc(y = as.integer(res_cat$y), mupred = bad_mupred),
    regexp = "`mupred` must contain values in \\[0, 1\\]"
  )
})

testthat::test_that("measure_bacc() with log-weights works as expected", {
  res <- measure_bacc(
    y = as.integer(res_cat$y),
    mupred = res_cat$mupred,
    log_weights = res_cat$loo$psis_object$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})
