# load test data --------------------------------------
res_roaches <- readRDS("data-for-tests/res_roaches.Rds")
res_binom <- readRDS("data-for-tests/res_binomial.Rds")
res_binary <- readRDS("data-for-tests/res_binary.Rds")
res_cat <- readRDS("data-for-tests/res_penguins.Rds")

set.seed(123456789)
n <- 10; S <- 100

res2 <- list(
  y = rnorm(n),
  x1 = matrix(rnorm(n * S), nrow = S),
  x2 = matrix(rnorm(n * S), nrow = S),
  ll = matrix(rnorm(n * S) * 0.1 - 1, nrow = S)
)

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
})

# rps() -------------------------------------

testthat::test_that("rps() works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    size = 10,
    log_weights = NULL
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binary$y))
})

testthat::test_that("srps() (scaled version) works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    size = 10,
    log_weights = NULL,
    scale = TRUE
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binary$y))
})

testthat::test_that("rps() with log-weights works as expected", {
  res <- rps(
    y = res_binom$y,
    ypred = res_binom$ypred,
    size = 10,
    log_weights = res_binom$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_binary$y))
  expect_true(all(res$pointwise >= 0))
})

# brier() ---------------------------------------

testthat::test_that("brier() works as expected", {
  res_brier <- brier(y = res_binary$y, ypred = res_binary$ypred, log_weights = NULL)

  expect_equal(names(res_brier), c("estimates", "pointwise"))
  expect_equal(length(res_brier$estimates), 2)
  expect_equal(length(res_brier$pointwise), length(res_binary$y))
  expect_true(all(res_brier$pointwise >= 0 & res_brier$pointwise <= 1))
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

# crps() ---------------------------------------
with_seed <- function(seed, code) {
  code <- substitute(code)
  orig.seed <- .Random.seed
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed)
  eval.parent(code)
}

test_that("crps computation is correct", {
  expect_equal(.crps_fun(2.0, 1.0), 0.0)
  expect_equal(.crps_fun(1.0, 2.0), -1.5)
  expect_equal(.crps_fun(pi, pi^2), 0.5 * pi - pi^2)

  expect_equal(.crps_fun(1.0, 0.0, scale = TRUE), 0.0)
  expect_equal(.crps_fun(1.0, 2.0, scale = TRUE), -2.0)
  expect_equal(.crps_fun(pi, pi^2, scale = TRUE), -pi^2/pi - 0.5 * log(pi))
})

test_that("crps matches snapshots", {
  expect_snapshot_value(with_seed(1, crps(res2$x1, res2$x2, res2$y)), 
  style = "serialize")
  expect_snapshot_value(with_seed(1, scrps(res2$x1, res2$x2, res2$y)), 
  style = "serialize")
  expect_snapshot_value(
    with_seed(1, loo_crps(res2$x1, res2$x2, res2$y, res2$ll)), 
    style = "serialize")
  expect_snapshot_value(
    with_seed(1, loo_scrps(res2$x1, res2$x2, res2$y, res2$ll)), 
    style = "serialize")
})

test_that("input validation throws correct errors", {
  expect_error(validate_crps_input(as.character(res2$x1), res2$x2, res2$y),
               "is.numeric(x) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(res2$x1, as.character(res2$x2), res2$y),
               "is.numeric(x2) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(res2$x1, res2$x2, c('a', 'b')),
               "is.numeric(y) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(res2$x1, t(res2$x2), res2$y),
               "identical(dim(x), dim(x2)) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(res2$x1, res2$x2, c(1, 2)),
               "ncol(x) == length(y) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(res2$x1, res2$x2, res2$y, t(res2$ll)),
               "ifelse(is.null(log_lik), TRUE, identical(dim(log_lik), dim(x))) is not TRUE",
               fixed = TRUE)
})

test_that("methods for single data point don't error", {
  expect_silent(crps(res2$x1[,1], res2$x2[,1], res2$y[1]))
  expect_silent(scrps(res2$x1[,1], res2$x2[,1], res2$y[1]))
})

testthat::test_that("crps2() works as expected", {
  res <- crps2(res_roaches$y, res_roaches$ypred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
})

testthat::test_that("crps2() with log-weights works as expected", {
  res <- crps2(
    res_roaches$y, 
    res_roaches$ypred, 
    log_weights = res_roaches$predperf_loo$log_weights
  )

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))

  # unnormalized log-weights as input
  res2 <- crps2(
    res_roaches$y, 
    res_roaches$ypred, 
    log_weights = res_roaches$predperf_loo$psis_object$log_weights
  )
  expect_equal(res$pointwise, res2$pointwise)
})

testthat::test_that("scrps2() (scaled version) works as expected", {
  res <- crps2(res_roaches$y, res_roaches$ypred, log_weights = NULL, scale = TRUE)
  res2 <- scrps2(res_roaches$y, res_roaches$ypred, log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
  expect_equal(res$estimates, res2$estimates)
  expect_equal(res$pointwise, res2$pointwise)
})

# mae() ------------------------------------------------
testthat::test_that("mae() works as expected", {
  res <- mae(y = res_roaches$y, mupred = res_roaches$mupred,
    log_weights = NULL)

  expect_equal(names(res), c("estimates", "pointwise"))
  expect_equal(length(res$estimates), 2)
  expect_equal(length(res$pointwise), length(res_roaches$y))
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
  expect_equal(res_mse$estimates[2]/(2*sqrt(abs(res_mse$estimates[1]))), res_rmse$estimates[2])
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
