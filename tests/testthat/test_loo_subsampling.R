options(mc.cores = 1)

test_that("overall loo_subampling works as expected (compared with loo) for diff_est", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    # each time called internally within loo the arguments will be equal to:
    # data_i: ith row of fdata (fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 500,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")

  # Check consistency
  expect_equal(
    loo_ss$pointwise[, "elpd_loo_approx"],
    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]],
    ignore_attr = TRUE
  )

  # Expect values
  z <- 2
  expect_lte(
    loo_ss$estimates["elpd_loo", "Estimate"] -
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["elpd_loo", "Estimate"] +
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["p_loo", "Estimate"] -
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["p_loo", "Estimate"] +
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["looic", "Estimate"] -
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["looic", "Estimate"] +
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  # Test that observations works as expected
  expect_message(
    loo_ss2 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = obs_idx(loo_ss),
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(loo_ss2$estimates, loo_ss$estimates, tolerance = 0.00000001)
  expect_silent(
    loo_ss2 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = loo_ss,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(loo_ss2$estimates, loo_ss$estimates, tolerance = 0.00000001)

  # Test lpd
  expect_silent(
    loo_ss_lpd <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 500,
      loo_approximation = "lpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss_lpd, "psis_loo_ss")
  z <- 2
  expect_lte(
    loo_ss_lpd$estimates["elpd_loo", "Estimate"] -
      z * loo_ss_lpd$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd$estimates["elpd_loo", "Estimate"] +
      z * loo_ss_lpd$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss_lpd$estimates["p_loo", "Estimate"] -
      z * loo_ss_lpd$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd$estimates["p_loo", "Estimate"] +
      z * loo_ss_lpd$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss_lpd$estimates["looic", "Estimate"] -
      z * loo_ss_lpd$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd$estimates["looic", "Estimate"] +
      z * loo_ss_lpd$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss_lpd$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss_lpd$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss_lpd$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  expect_silent(
    loo_ss_lpd10 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 500,
      loo_approximation = "lpd",
      loo_approximation_draws = 10,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss_lpd10, "psis_loo_ss")

  z <- 2
  expect_lte(
    loo_ss_lpd10$estimates["elpd_loo", "Estimate"] -
      z * loo_ss_lpd10$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd10$estimates["elpd_loo", "Estimate"] +
      z * loo_ss_lpd10$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss_lpd10$estimates["p_loo", "Estimate"] -
      z * loo_ss_lpd10$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd10$estimates["p_loo", "Estimate"] +
      z * loo_ss_lpd10$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss_lpd10$estimates["looic", "Estimate"] -
      z * loo_ss_lpd10$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss_lpd10$estimates["looic", "Estimate"] +
      z * loo_ss_lpd10$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss_lpd10$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss_lpd10$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss_lpd10$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  # Test conversion of objects
  expect_silent(true_loo_2 <- loo:::as.psis_loo.psis_loo(true_loo))
  expect_silent(true_loo_ss <- loo:::as.psis_loo_ss.psis_loo(true_loo))
  expect_s3_class(true_loo_ss, "psis_loo_ss")
  expect_silent(true_loo_conv <- loo:::as.psis_loo.psis_loo_ss(true_loo_ss))
  expect_failure(expect_s3_class(true_loo_conv, "psis_loo_ss"))
  expect_equal(true_loo_conv, true_loo)
  expect_error(loo:::as.psis_loo.psis_loo_ss(loo_ss))
})

test_that("loo with subsampling of all observations works as ordinary loo.", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 1000,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_error(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 1001,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )

  expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  )
  expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  )
  expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  )

  expect_equal(dim(true_loo), dim(loo_ss))
  expect_equal(true_loo$diagnostics, loo_ss$diagnostics)
  expect_equal(max(loo_ss$pointwise[, "m_i"]), 1)
})

test_that("overall loo_subsample works with diff_srs as expected (compared with loo)", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 200,
      loo_approximation = "plpd",
      estimator = "diff_srs",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(
    true_loo$estimates[1, 1],
    loo_ss$estimates[1, 1],
    tolerance = 0.1
  )
})

test_that("Test the srs estimator with 'none' approximation", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 200,
      loo_approximation = "none",
      estimator = "srs",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_error(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 1100,
      loo_approximation = "none",
      estimator = "srs",
      r_eff = rep(1, nrow(fake_data))
    )
  )

  expect_equal(length(obs_idx(loo_ss)), nobs(loo_ss))

  # Check consistency
  expect_equal(
    loo_ss$pointwise[, "elpd_loo_approx"],
    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]],
    ignore_attr = TRUE
  )

  # Expect values
  z <- 2
  expect_lte(
    loo_ss$estimates["elpd_loo", "Estimate"] -
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["elpd_loo", "Estimate"] +
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["p_loo", "Estimate"] -
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["p_loo", "Estimate"] +
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["looic", "Estimate"] -
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["looic", "Estimate"] +
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))
})

test_that("Test the Hansen-Hurwitz estimator", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 300,
      loo_approximation = "plpd",
      estimator = "hh_pps",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_silent(
    loo_ss_max <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 1100,
      loo_approximation = "plpd",
      estimator = "hh_pps",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss_max, "psis_loo_ss")
  expect_silent(
    loo_ss_max2 <- update(
      loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 1100,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(nobs(loo_ss_max2), 1100)
  expect_gt(max(loo_ss_max2$pointwise[, "m_i"]), 1)
  expect_error(
    loo_ss_max2 <- update(
      loo_ss_max2,
      draws = fake_posterior,
      data = fake_data,
      observations = 300,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss_max3 <- update(
      loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 1500,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss2 <- update(
      loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = loo_ss,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_error(
    loo_ss2 <- update(
      loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = loo_ss,
      loo_approximation = "lpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(loo_ss$estimates, loo_ss2$estimates)
  expect_equal(length(obs_idx(loo_ss_max)), length(obs_idx(loo_ss_max2)))
  expect_equal(length(obs_idx(loo_ss_max)), nobs(loo_ss_max))

  # Check consistency
  expect_equal(
    loo_ss$pointwise[, "elpd_loo_approx"],
    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]],
    ignore_attr = TRUE
  )
  # Check consistency
  expect_equal(
    loo_ss_max$pointwise[, "elpd_loo_approx"],
    loo_ss_max$loo_subsampling$elpd_loo_approx[loo_ss_max$pointwise[, "idx"]],
    ignore_attr = TRUE
  )

  # Expect values
  z <- 2
  expect_lte(
    loo_ss$estimates["elpd_loo", "Estimate"] -
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["elpd_loo", "Estimate"] +
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["p_loo", "Estimate"] -
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["p_loo", "Estimate"] +
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    true_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["looic", "Estimate"] -
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["looic", "Estimate"] +
      z * loo_ss$estimates["looic", "subsampling SE"],
    true_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    true_loo$estimates["elpd_loo", "Estimate"],
    loo_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["p_loo", "Estimate"],
    loo_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    true_loo$estimates["looic", "Estimate"],
    loo_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  expect_lte(
    loo_ss_max$estimates["elpd_loo", "Estimate"] -
      z * loo_ss_max$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss_max$estimates["elpd_loo", "Estimate"] +
      z * loo_ss_max$estimates["elpd_loo", "subsampling SE"],
    true_loo$estimates["elpd_loo", "Estimate"]
  )
})


test_that("update.psis_loo_ss works as expected (compared with loo)", {
  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    true_loo <- loo(
      llfun_test,
      draws = fake_posterior,
      data = fake_data,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 500,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")

  # Check error when draws and data dimensions differ
  expect_error(
    loo_ss2 <- update(
      object = loo_ss,
      draws = cbind(fake_posterior, 1),
      data = fake_data,
      observations = 600,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_error(
    loo_ss2 <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data[-1, ],
      observations = 600,
      r_eff = rep(1, nrow(fake_data))
    )
  )

  # Add tests for adding observations
  expect_silent(
    loo_ss2 <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 600,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(dim(loo_ss2)[2] - dim(loo_ss)[2], expected = 100)
  expect_equal(dim(loo_ss2)[2], expected = dim(loo_ss2$pointwise)[1])
  expect_length(loo_ss2$diagnostics$pareto_k, 600)
  expect_length(loo_ss2$diagnostics$n_eff, 600)
  for (i in 1:nrow(loo_ss2$estimates)) {
    expect_lt(
      loo_ss2$estimates[i, "subsampling SE"],
      loo_ss$estimates[i, "subsampling SE"]
    )
  }

  expect_silent(
    loo_ss2b <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data
    )
  )
  expect_equal(loo_ss2b$estimates, loo_ss$estimates)
  expect_equal(loo_ss2b$pointwise, loo_ss$pointwise)
  expect_equal(loo_ss2b$diagnostics$pareto_k, loo_ss$diagnostics$pareto_k)
  expect_equal(loo_ss2b$diagnostics$n_eff, loo_ss$diagnostics$n_eff)

  expect_silent(
    loo_ss3 <- update(
      object = loo_ss2,
      draws = fake_posterior,
      data = fake_data,
      observations = loo_ss
    )
  )
  expect_equal(loo_ss3$estimates, loo_ss$estimates)
  expect_equal(loo_ss3$pointwise, loo_ss$pointwise)
  expect_equal(loo_ss3$diagnostics$pareto_k, loo_ss$diagnostics$pareto_k)
  expect_equal(loo_ss3$diagnostics$n_eff, loo_ss$diagnostics$n_eff)

  expect_silent(
    loo_ss4 <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 1000,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(loo_ss4$estimates[, 1], true_loo$estimates[, 1])
  expect_equal(
    loo_ss4$estimates[, 2],
    true_loo$estimates[, 2],
    tolerance = 0.001
  )

  expect_silent(
    loo_ss5 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 1000,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )

  ss4_order <- order(loo_ss4$pointwise[, "idx"])
  expect_equal(
    loo_ss4$pointwise[ss4_order, c(1, 3, 4)],
    loo_ss5$pointwise[, c(1, 3, 4)]
  )
  expect_equal(
    loo_ss4$diagnostics$pareto_k[ss4_order],
    loo_ss5$diagnostics$pareto_k
  )
  expect_equal(loo_ss4$diagnostics$n_eff[ss4_order], loo_ss5$diagnostics$n_eff)
  expect_equal(
    loo_ss4$pointwise[ss4_order, c(1, 3, 4)],
    true_loo$pointwise[, c(1, 3, 4)]
  )
  expect_equal(
    loo_ss4$diagnostics$pareto_k[ss4_order],
    true_loo$diagnostics$pareto_k
  )
  expect_equal(loo_ss4$diagnostics$n_eff[ss4_order], true_loo$diagnostics$n_eff)

  expect_error(
    loo_ss_min <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 50,
      r_eff = rep(1, nrow(fake_data))
    )
  )

  expect_silent(true_loo_ss <- loo:::as.psis_loo_ss.psis_loo(true_loo))
  expect_silent(
    loo_ss_subset0 <- update(
      true_loo_ss,
      observations = loo_ss,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_true(identical(obs_idx(loo_ss_subset0), obs_idx(loo_ss)))
  expect_silent(
    loo_ss_subset1 <- update(
      object = loo_ss,
      observations = loo_ss,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_message(
    loo_ss_subset2 <- update(
      object = loo_ss,
      observations = obs_idx(loo_ss)[1:10],
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_equal(nobs(loo_ss_subset2), 10)

  expect_silent(true_loo_ss <- loo:::as.psis_loo_ss.psis_loo(true_loo))
  set.seed(4711)
  expect_silent(
    loo_ss2 <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data,
      observations = 600,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss2_subset0 <- update(
      object = true_loo_ss,
      observations = loo_ss2,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_true(setequal(obs_idx(loo_ss2), obs_idx(loo_ss2_subset0)))
  expect_true(identical(obs_idx(loo_ss2), obs_idx(loo_ss2_subset0)))
  expect_true(identical(loo_ss2$diagnostic, loo_ss2_subset0$diagnostic))

  # Add tests for changing approx variable
  expect_silent(
    loo_ss_lpd <- update(
      object = loo_ss,
      draws = fake_posterior,
      data = fake_data,
      loo_approximation = "lpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_failure(expect_equal(
    loo_ss_lpd$loo_subsampling$elpd_loo_approx,
    loo_ss$loo_subsampling$elpd_loo_approx
  ))
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss)[2])
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss_lpd$pointwise)[1])
  expect_length(loo_ss_lpd$diagnostics$pareto_k, 500)
  expect_length(loo_ss_lpd$diagnostics$n_eff, 500)
  expect_failure(expect_equal(
    loo_ss_lpd$estimates[1, "subsampling SE"],
    loo_ss$estimates[1, "subsampling SE"]
  ))
  expect_failure(expect_equal(
    loo_ss_lpd$estimates[3, "subsampling SE"],
    loo_ss$estimates[3, "subsampling SE"]
  ))
})

test_that("loo_compare_subsample", {
  skip_on_cran() # to get under cran check time limit

  set.seed(123)
  N <- 1000
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  sigma <- 2
  y <- rnorm(N, 1 + 2 * x1 - 2 * x2 - 1 * x3, sd = sigma)
  X <- cbind("x0" = rep(1, N), x1, x2, x3)

  # Generate samples from posterior
  samples_blin <- function(X, y, sigma, draws = 1000) {
    XtX <- t(X) %*% X
    b_hat <- solve(XtX) %*% (t(X) %*% y)
    Lambda_n = XtX + diag(ncol(X))
    mu_n <- solve(Lambda_n) %*%
      (XtX %*% b_hat + diag(ncol(X)) %*% rep(0, ncol(X)))
    L <- t(chol(sigma^2 * solve(Lambda_n)))
    draws_mat <- matrix(0, ncol = ncol(X), nrow = draws)
    for (i in 1:draws) {
      z <- rnorm(length(mu_n))
      draws_mat[i, ] <- L %*% z + mu_n
    }
    draws_mat
  }

  fake_posterior1 <- samples_blin(X[, 1:2], y, sigma, draws = 1000)
  fake_posterior2 <- samples_blin(X[, 1:3], y, sigma, draws = 1000)
  fake_posterior3 <- samples_blin(X, y, sigma, draws = 1000)

  fake_data1 <- data.frame(y, X[, 1:2])
  fake_data2 <- data.frame(y, X[, 1:3])
  fake_data3 <- data.frame(y, X)

  llfun_test <- function(data_i, draws) {
    dnorm(
      x = data_i$y,
      mean = draws %*% t(data_i[, -1, drop = FALSE]),
      sd = sigma,
      log = TRUE
    )
  }

  expect_silent(
    l1 <- loo(
      llfun_test,
      data = fake_data1,
      draws = fake_posterior1,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    l2 <- loo(
      llfun_test,
      data = fake_data2,
      draws = fake_posterior2,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    l3 <- loo(
      llfun_test,
      data = fake_data3,
      draws = fake_posterior3,
      r_eff = rep(1, N)
    )
  )

  expect_silent(
    lss1 <- loo_subsample(
      llfun_test,
      data = fake_data1,
      draws = fake_posterior1,
      observations = 100,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    lss2 <- loo_subsample(
      llfun_test,
      data = fake_data2,
      draws = fake_posterior2,
      observations = 100,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    lss3 <- loo_subsample(
      llfun_test,
      data = fake_data3,
      draws = fake_posterior3,
      observations = 100,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    lss2o1 <- loo_subsample(
      llfun_test,
      data = fake_data2,
      draws = fake_posterior2,
      observations = lss1,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    lss3o1 <- loo_subsample(
      llfun_test,
      data = fake_data3,
      draws = fake_posterior3,
      observations = lss1,
      r_eff = rep(1, N)
    )
  )
  expect_silent(
    lss2hh <- loo_subsample(
      llfun_test,
      data = fake_data2,
      draws = fake_posterior2,
      observations = 100,
      estimator = "hh_pps",
      r_eff = rep(1, N)
    )
  )

  expect_snapshot(
    lcss <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2, lss3))
  )
  expect_warning(
    lcss2 <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2, lss3o1))
  )
  expect_silent(
    lcsso <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2o1, lss3o1))
  )
  expect_warning(
    lcssohh <- loo:::loo_compare.psis_loo_ss_list(
      x = list(lss1, lss2hh, lss3o1)
    )
  )
  expect_message(
    lcssf1 <- loo:::loo_compare.psis_loo_ss_list(
      x = list(loo:::as.psis_loo_ss.psis_loo(l1), lss2o1, lss3o1)
    )
  )
  expect_message(
    lcssf2 <- loo:::loo_compare.psis_loo_ss_list(
      x = list(
        loo:::as.psis_loo_ss.psis_loo(l1),
        lss2o1,
        loo:::as.psis_loo_ss.psis_loo(l3)
      )
    )
  )

  expect_equal(lcss[, 1], lcsso[, 1], tolerance = 1)
  expect_equal(lcss2[, 1], lcsso[, 1], tolerance = 1)
  expect_equal(lcssohh[, 1], lcsso[, 1], tolerance = 1)
  expect_equal(lcssf1[, 1], lcsso[, 1], tolerance = 1)
  expect_equal(lcssf2[, 1], lcsso[, 1], tolerance = 1)

  expect_gt(lcss[, 2][2], lcsso[, 2][2])
  expect_gt(lcss[, 2][3], lcsso[, 2][3])
  expect_gt(lcss2[, 2][2], lcsso[, 2][2])
  expect_equal(lcss2[, 2][3], lcsso[, 2][3])
  expect_gt(lcssohh[, 2][2], lcsso[, 2][2])
  expect_equal(lcssohh[, 2][3], lcsso[, 2][3])

  expect_silent(
    lcss2m <- loo:::loo_compare.psis_loo_ss_list(x = list(lss2o1, lss3o1))
  )
  expect_equal(unname(lcss2m[,]), unname(lcsso[1:2, ]))

  expect_snapshot(lcssapi <- loo_compare(lss1, lss2, lss3))
  expect_equal(lcssapi, lcss)
  expect_warning(lcssohhapi <- loo_compare(lss1, lss2hh, lss3o1))
  expect_equal(lcssohhapi, lcssohh)
  expect_silent(lcss2mapi <- loo_compare(lss2o1, lss3o1))
  expect_equal(lcss2mapi, lcss2m)
})

test_that("Test 'tis' and 'sis'", {
  skip_on_cran()

  set.seed(123)
  N <- 1000
  K <- 10
  S <- 1000
  a0 <- 3
  b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y)
  b <- b0 + N * K - sum(y)
  fake_posterior <- draws <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y, K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(
    loo_ss_full <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 1000,
        loo_approximation = "plpd",
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_plpd <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "plpd",
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_tis_S1000 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "tis",
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_tis_S100 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "tis",
        loo_approximation_draws = 100,
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_tis_S10 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "tis",
        loo_approximation_draws = 10,
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_sis_S1000 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "sis",
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_sis_S100 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "sis",
        loo_approximation_draws = 100,
        r_eff = rep(1, nrow(fake_data))
      )
  )
  expect_silent(
    loo_ss_sis_S10 <-
      loo_subsample(
        x = llfun_test,
        draws = fake_posterior,
        data = fake_data,
        observations = 100,
        loo_approximation = "sis",
        loo_approximation_draws = 10,
        r_eff = rep(1, nrow(fake_data))
      )
  )

  SEs <- 4
  expect_gt(
    loo_ss_tis_S1000$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_tis_S1000$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_tis_S1000$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_tis_S1000$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_gt(
    loo_ss_tis_S100$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_tis_S100$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_tis_S100$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_tis_S100$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_gt(
    loo_ss_tis_S10$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_tis_S10$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_tis_S10$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_tis_S10$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )

  expect_gt(
    loo_ss_sis_S1000$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_sis_S1000$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_sis_S1000$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_sis_S1000$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_gt(
    loo_ss_sis_S100$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_sis_S100$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_sis_S100$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_sis_S100$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_gt(
    loo_ss_sis_S10$estimates["elpd_loo", "Estimate"] +
      SEs * loo_ss_sis_S10$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lt(
    loo_ss_sis_S10$estimates["elpd_loo", "Estimate"] -
      SEs * loo_ss_sis_S10$estimates["elpd_loo", "subsampling SE"],
    loo_ss_full$estimates["elpd_loo", "Estimate"]
  )
})
