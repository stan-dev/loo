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

testthat::test_that("measure_rps() validates the shape of log_weights", {
  bad_log_weights <- res_binom$log_weights[, -1, drop = FALSE]
  expect_error(
    measure_rps(
      y = res_binom$y,
      ypred = res_binom$ypred,
      log_weights = bad_log_weights
    ),
    regexp = "`log_weights` must have"
  )
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
    log_weights = res_roaches$log_weights)

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

testthat::test_that("higher_is_better reorients loss measures to utility scale", {
  res_mse <- measure_mse(y = res_roaches$y, mupred = res_roaches$mupred)
  res_mse_utility <- measure_mse(
    y = res_roaches$y,
    mupred = res_roaches$mupred,
    higher_is_better = TRUE
  )

  expect_equal(
    unname(res_mse_utility$estimates["Estimate"]),
    -unname(res_mse$estimates["Estimate"])
  )
  expect_equal(res_mse_utility$pointwise, -res_mse$pointwise)
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
    log_weights = res_roaches$log_weights)

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
    log_weights = res_cat$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})

testthat::test_that("measure_acc() broadcasts weights over categories without warning", {
  y <- as.integer(res_cat$y)
  lw <- res_cat$log_weights

  expect_no_warning(res <- measure_acc(
    y = y, mupred = res_cat$mupred, log_weights = lw
  ))
  expect_no_warning(measure_bacc(
    y = y, mupred = res_cat$mupred, log_weights = lw
  ))

  w <- exp(.normalize_log_weights(lw))
  n_cat <- dim(res_cat$mupred)[3]
  manual <- vapply(
    seq_len(n_cat),
    function(k) colSums(w * res_cat$mupred[, , k]),
    numeric(length(y))
  )
  expect_equal(as.numeric(res$pointwise), (apply(manual, 1, which.max) == y) * 1)
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

testthat::test_that("measure_bacc() accepts precomputed pointwise values", {
  y <- c(1L, 1L, 2L, 2L)
  acc_i <- c(1L, 0L, 1L, 1L)
  res <- measure_bacc(y = y, mupred = NULL, pointwise = acc_i)

  expect_equal(unname(res$estimates[1, "Estimate"]), 0.75)
  expect_error(
    measure_bacc(y = y, mupred = NULL, pointwise = acc_i[-1]),
    regexp = "must have the same length"
  )
})

testthat::test_that("mlpd and ic count draws after the 3-D conversion", {
  LLarr <- example_loglik_array()
  dims_elpd <- attr(measure_elpd(LLarr), "dims")

  expect_equal(attr(measure_mlpd(LLarr), "dims"), dims_elpd)
  expect_equal(attr(measure_ic(LLarr), "dims"), dims_elpd)
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
    log_weights = res_cat$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_cat$y))
  expect_true(!all(res$pointwise < 0 | res$pointwise > 1))
})

# mlpd / ic derived from elpd ------------------------
testthat::test_that("mlpd and ic keep the sign of a positive pointwise lppd", {
  # lppd_i is a log density and is positive for concentrated predictions.
  set.seed(1)
  S <- 200
  n <- 30
  y <- rnorm(n, 0, 0.05)
  mu <- matrix(rnorm(S * n, 0, 0.01), S, n)
  ylp <- t(sapply(1:S, function(s) dnorm(y, mu[s, ], 0.05, log = TRUE)))
  lppd <- matrixStats::colLogSumExps(ylp) - log(S)

  # the removed heuristic negated every entry as soon as one was positive
  expect_true(any(lppd > 0))
  expect_true(any(lppd < 0))

  res <- insample_pred_measure(ylp = ylp, measure = c("mlpd", "ic"))

  expect_equal(unname(res$estimates["mlpd", 1]), mean(lppd))
  expect_equal(unname(res$estimates["ic", 1]), sum(-2 * lppd))
  expect_equal(unname(res$pointwise[, "mlpd"]), unname(lppd))
  expect_equal(unname(res$pointwise[, "ic"]), unname(-2 * lppd))
  expect_equal(cor(res$pointwise[, "elpd"], res$pointwise[, "mlpd"]), 1)
})

testthat::test_that("mlpd on the test source uses ylp_test, not ylp", {
  set.seed(2)
  S <- 100
  ylp <- matrix(rnorm(S * 20, -1, 0.1), S, 20)
  ylp_test <- matrix(rnorm(S * 8, -5, 0.1), S, 8)
  lppd_test <- matrixStats::colLogSumExps(ylp_test) - log(S)

  res <- test_pred_measure(ylp = ylp, ylp_test = ylp_test, measure = "mlpd")

  expect_equal(unname(res$estimates["mlpd_test", 1]), mean(lppd_test))
  expect_equal(unname(res$pointwise[, "mlpd_test"]), unname(lppd_test))
})

testthat::test_that("mlpd and ic work when only a loo object is given", {
  LL <- example_loglik_matrix()
  lo <- suppressWarnings(
    loo(LL, save_psis = TRUE, r_eff = rep(1, ncol(LL)))
  )

  res <- loo_pred_measure(loo = lo, measure = c("mlpd", "ic"))

  elpd_loo_i <- lo$pointwise[, "elpd_loo"]
  expect_equal(unname(res$estimates["mlpd_loo", 1]), mean(elpd_loo_i))
  expect_equal(unname(res$estimates["ic_loo", 1]), sum(-2 * elpd_loo_i))
})

# .exx_pwm ------------------------
testthat::test_that(".exx_pwm() reduces to the classic unbiased PWM estimator", {
  set.seed(3)
  S <- 500
  x <- matrix(rnorm(S * 6), S, 6)
  ref <- colMeans(apply(x, 2, sort) * (2 * (2 * seq_len(S) - S - 1) / (S - 1)))

  expect_equal(.exx_pwm(x), ref)
  expect_equal(.exx_pwm(x, matrix(1 / S, S, 6)), .exx_pwm(x))
})

testthat::test_that(".exx_pwm() stays non-negative under concentrated weights", {
  # E|X - X'| is non-negative by definition. The previous inline estimator
  # divided by (1 - w) and drove the estimate negative when one draw dominated,
  # which made log(EXX) in srps NaN.
  set.seed(4)
  S <- 200
  x <- matrix(rgamma(S * 4, 2, 1), S, 4)
  lw <- matrix(0, S, 4)
  lw[1, ] <- 8 # one draw carries most, but not all, of the weight
  w <- exp(.normalize_log_weights(lw))
  expect_gt(max(w), 0.9)

  EXX <- .exx_pwm(x, w)
  expect_true(all(EXX > 0))
  expect_true(all(is.finite(EXX)))

  srps <- measure_rps(y = colMeans(x), ypred = x, log_weights = lw, scaled = TRUE)
  expect_true(all(is.finite(srps$pointwise)))
})

testthat::test_that("the scaled score aborts on a point-mass predictive", {
  # All the weight on one draw makes E|X - X'| exactly 0, so the scaling term
  # log(E|X - X'|) is undefined. The unscaled score stays well defined.
  set.seed(5)
  S <- 100
  x <- matrix(rgamma(S * 3, 2, 1), S, 3)
  lw <- matrix(0, S, 3)
  lw[1, ] <- 60
  y <- colMeans(x)

  expect_equal(.exx_pwm(x, exp(.normalize_log_weights(lw))), rep(0, 3))
  expect_error(
    measure_rps(y = y, ypred = x, log_weights = lw, scaled = TRUE),
    "point mass"
  )
  expect_true(all(is.finite(
    measure_rps(y = y, ypred = x, log_weights = lw)$pointwise
  )))
})

testthat::test_that(".exx_pwm() requires at least two draws", {
  expect_error(.exx_pwm(matrix(1, 1, 3)), "at least 2 draws")
})
