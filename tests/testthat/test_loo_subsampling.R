library(loo)
options(mc.cores = 1)

context("loo_subsampling")

test_that("overall loo_subampling works as expected (compared with loo) for diff_est", {
  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    # each time called internally within loo the arguments will be equal to:
    # data_i: ith row of fdata (fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 500, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss, "psis_loo_ss")

  # Check consistency
  expect_equivalent(loo_ss$pointwise[, "elpd_loo_approx"],
                    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]])

  # Expect values
  z <- 2
  expect_lte(loo_ss$estimates["elpd_loo", "Estimate"] - z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss$estimates["elpd_loo", "Estimate"] + z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss$estimates["p_loo", "Estimate"] - z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss$estimates["p_loo", "Estimate"] + z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss$estimates["looic", "Estimate"] - z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss$estimates["looic", "Estimate"] + z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss$estimates["looic", "Estimate"], tol = 0.00000001))

  # Test that observations works as expected
  expect_message(loo_ss2 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = obs_idx(loo_ss), loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_equal(loo_ss2$estimates, loo_ss$estimates, tol = 0.00000001)
  expect_silent(loo_ss2 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = loo_ss, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_equal(loo_ss2$estimates, loo_ss$estimates, tol = 0.00000001)

  # Test lpd
  expect_silent(loo_ss_lpd <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 500, loo_approximation = "lpd", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss_lpd, "psis_loo_ss")
  z <- 2
  expect_lte(loo_ss_lpd$estimates["elpd_loo", "Estimate"] - z * loo_ss_lpd$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss_lpd$estimates["elpd_loo", "Estimate"] + z * loo_ss_lpd$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss_lpd$estimates["p_loo", "Estimate"] - z * loo_ss_lpd$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss_lpd$estimates["p_loo", "Estimate"] + z * loo_ss_lpd$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss_lpd$estimates["looic", "Estimate"] - z * loo_ss_lpd$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss_lpd$estimates["looic", "Estimate"] + z * loo_ss_lpd$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss_lpd$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss_lpd$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss_lpd$estimates["looic", "Estimate"], tol = 0.00000001))

  expect_silent(loo_ss_lpd10 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 500, loo_approximation = "lpd", loo_approximation_draws = 10, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss_lpd10, "psis_loo_ss")

  z <- 2
  expect_lte(loo_ss_lpd10$estimates["elpd_loo", "Estimate"] - z * loo_ss_lpd10$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss_lpd10$estimates["elpd_loo", "Estimate"] + z * loo_ss_lpd10$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss_lpd10$estimates["p_loo", "Estimate"] - z * loo_ss_lpd10$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss_lpd10$estimates["p_loo", "Estimate"] + z * loo_ss_lpd10$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss_lpd10$estimates["looic", "Estimate"] - z * loo_ss_lpd10$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss_lpd10$estimates["looic", "Estimate"] + z * loo_ss_lpd10$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss_lpd10$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss_lpd10$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss_lpd10$estimates["looic", "Estimate"], tol = 0.00000001))

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
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 1000, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_error(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 1001, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))

  expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001)
  expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss$estimates["p_loo", "Estimate"], tol = 0.00000001)
  expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss$estimates["looic", "Estimate"], tol = 0.00000001)

  expect_equal(dim(true_loo),dim(loo_ss))
  expect_equal(true_loo$diagnostics, loo_ss$diagnostics)
  expect_equal(max(loo_ss$pointwise[, "m_i"]), 1)

})

test_that("overall loo_subsample works with diff_srs as expected (compared with loo)", {
  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(x = llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 200, loo_approximation = "plpd", estimator = "diff_srs", r_eff  = rep(1, nrow(fake_data))))
  expect_equal(true_loo$estimates[1,1], loo_ss$estimates[1,1], tol = 0.1)

})

test_that("Test the srs estimator with 'none' approximation", {
  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 200, loo_approximation = "none", estimator = "srs", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_error(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 1100, loo_approximation = "none", estimator = "srs", r_eff = rep(1, nrow(fake_data))))

  expect_equal(length(obs_idx(loo_ss)), nobs(loo_ss))


  # Check consistency
  expect_equivalent(loo_ss$pointwise[, "elpd_loo_approx"],
                    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]])

  # Expect values
  z <- 2
  expect_lte(loo_ss$estimates["elpd_loo", "Estimate"] - z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss$estimates["elpd_loo", "Estimate"] + z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss$estimates["p_loo", "Estimate"] - z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss$estimates["p_loo", "Estimate"] + z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss$estimates["looic", "Estimate"] - z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss$estimates["looic", "Estimate"] + z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss$estimates["looic", "Estimate"], tol = 0.00000001))

})

test_that("Test the Hansen-Hurwitz estimator", {
  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 300, loo_approximation = "plpd", estimator = "hh_pps", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss, "psis_loo_ss")
  expect_silent(loo_ss_max <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 1100, loo_approximation = "plpd", estimator = "hh_pps", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss_max, "psis_loo_ss")
  expect_silent(loo_ss_max2 <- update(loo_ss, draws = fake_posterior, data = fake_data, observations = 1100, r_eff = rep(1, nrow(fake_data))))
  expect_equal(nobs(loo_ss_max2), 1100)
  expect_gt(max(loo_ss_max2$pointwise[, "m_i"]), 1)
  expect_error(loo_ss_max2 <- update(loo_ss_max2, draws = fake_posterior, data = fake_data, observations = 300, r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_max3 <- update(loo_ss, draws = fake_posterior, data = fake_data, observations = 1500, r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss2 <- update(loo_ss, draws = fake_posterior, data = fake_data, observations = loo_ss, r_eff = rep(1, nrow(fake_data))))
  expect_error(loo_ss2 <- update(loo_ss, draws = fake_posterior, data = fake_data, observations = loo_ss, loo_approximation = "lpd", r_eff = rep(1, nrow(fake_data))))
  expect_equal(loo_ss$estimates, loo_ss2$estimates)
  expect_equal(length(obs_idx(loo_ss_max)), length(obs_idx(loo_ss_max2)))
  expect_equal(length(obs_idx(loo_ss_max)), nobs(loo_ss_max))


  # Check consistency
  expect_equivalent(loo_ss$pointwise[, "elpd_loo_approx"],
                    loo_ss$loo_subsampling$elpd_loo_approx[loo_ss$pointwise[, "idx"]])
  # Check consistency
  expect_equivalent(loo_ss_max$pointwise[, "elpd_loo_approx"],
                    loo_ss_max$loo_subsampling$elpd_loo_approx[loo_ss_max$pointwise[, "idx"]])

  # Expect values
  z <- 2
  expect_lte(loo_ss$estimates["elpd_loo", "Estimate"] - z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss$estimates["elpd_loo", "Estimate"] + z * loo_ss$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss$estimates["p_loo", "Estimate"] - z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss$estimates["p_loo", "Estimate"] + z * loo_ss$estimates["p_loo", "subsampling SE"], true_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss$estimates["looic", "Estimate"] - z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss$estimates["looic", "Estimate"] + z * loo_ss$estimates["looic", "subsampling SE"], true_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(true_loo$estimates["elpd_loo", "Estimate"], loo_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["p_loo", "Estimate"], loo_ss$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(true_loo$estimates["looic", "Estimate"], loo_ss$estimates["looic", "Estimate"], tol = 0.00000001))

  expect_lte(loo_ss_max$estimates["elpd_loo", "Estimate"] - z * loo_ss_max$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss_max$estimates["elpd_loo", "Estimate"] + z * loo_ss_max$estimates["elpd_loo", "subsampling SE"], true_loo$estimates["elpd_loo", "Estimate"])

})


test_that("update.psis_loo_ss works as expected (compared with loo)", {


  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(true_loo <- loo(llfun_test, draws = fake_posterior, data = fake_data, r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(true_loo, "psis_loo")
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 500, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_s3_class(loo_ss, "psis_loo_ss")

  # Check error when draws and data dimensions differ
  expect_error(loo_ss2 <- update(object = loo_ss, draws = cbind(fake_posterior, 1), data = fake_data, observations = 600, r_eff = rep(1, nrow(fake_data))))
  expect_error(loo_ss2 <- update(object = loo_ss, draws = fake_posterior, data = fake_data[-1,], observations = 600, r_eff = rep(1, nrow(fake_data))))

  # Add tests for adding observations
  expect_silent(loo_ss2 <- update(object = loo_ss, draws = fake_posterior, data = fake_data, observations = 600, r_eff = rep(1, nrow(fake_data))))
  expect_equal(dim(loo_ss2)[2] - dim(loo_ss)[2], expected = 100)
  expect_equal(dim(loo_ss2)[2], expected = dim(loo_ss2$pointwise)[1])
  expect_length(loo_ss2$diagnostics$pareto_k, 600)
  expect_length(loo_ss2$diagnostics$n_eff, 600)
  for(i in 1:nrow(loo_ss2$estimates)) {
    expect_lt(loo_ss2$estimates[i, "subsampling SE"],
              loo_ss$estimates[i, "subsampling SE"])
  }

  expect_silent(loo_ss2b <- update(object = loo_ss, draws = fake_posterior, data = fake_data))
  expect_equal(loo_ss2b$estimates, loo_ss$estimates)
  expect_equal(loo_ss2b$pointwise, loo_ss$pointwise)
  expect_equal(loo_ss2b$diagnostics$pareto_k, loo_ss$diagnostics$pareto_k)
  expect_equal(loo_ss2b$diagnostics$n_eff, loo_ss$diagnostics$n_eff)

  expect_silent(loo_ss3 <- update(object = loo_ss2, draws = fake_posterior, data = fake_data, observations = loo_ss))
  expect_equal(loo_ss3$estimates, loo_ss$estimates)
  expect_equal(loo_ss3$pointwise, loo_ss$pointwise)
  expect_equal(loo_ss3$diagnostics$pareto_k, loo_ss$diagnostics$pareto_k)
  expect_equal(loo_ss3$diagnostics$n_eff, loo_ss$diagnostics$n_eff)

  expect_silent(loo_ss4 <- update(object = loo_ss, draws = fake_posterior, data = fake_data, observations = 1000, r_eff = rep(1, nrow(fake_data))))
  expect_equal(loo_ss4$estimates[,1], true_loo$estimates[,1])
  expect_equal(loo_ss4$estimates[,2], true_loo$estimates[,2], tol = 0.001)

  expect_silent(loo_ss5 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 1000, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))

  ss4_order <- order(loo_ss4$pointwise[, "idx"])
  expect_equal(loo_ss4$pointwise[ss4_order,c(1,3,4)], loo_ss5$pointwise[,c(1,3,4)])
  expect_equal(loo_ss4$diagnostics$pareto_k[ss4_order], loo_ss5$diagnostics$pareto_k)
  expect_equal(loo_ss4$diagnostics$n_eff[ss4_order], loo_ss5$diagnostics$n_eff)
  expect_equal(loo_ss4$pointwise[ss4_order,c(1,3,4)], true_loo$pointwise[,c(1,3,4)])
  expect_equal(loo_ss4$diagnostics$pareto_k[ss4_order], true_loo$diagnostics$pareto_k)
  expect_equal(loo_ss4$diagnostics$n_eff[ss4_order], true_loo$diagnostics$n_eff)

  expect_error(loo_ss_min <- update(object = loo_ss, draws = fake_posterior, data = fake_data, observations = 50, r_eff = rep(1, nrow(fake_data))))

  expect_silent(true_loo_ss <- loo:::as.psis_loo_ss.psis_loo(true_loo))
  expect_silent(loo_ss_subset0 <- update(true_loo_ss, observations = loo_ss, r_eff = rep(1, nrow(fake_data))))
  expect_true(identical(obs_idx(loo_ss_subset0), obs_idx(loo_ss)))
  expect_silent(loo_ss_subset1 <- update(object = loo_ss, observations = loo_ss, r_eff = rep(1, nrow(fake_data))))
  expect_message(loo_ss_subset2 <- update(object = loo_ss, observations = obs_idx(loo_ss)[1:10], r_eff = rep(1, nrow(fake_data))))
  expect_equal(nobs(loo_ss_subset2), 10)


  expect_silent(true_loo_ss <- loo:::as.psis_loo_ss.psis_loo(true_loo))
  set.seed(4711)
  expect_silent(loo_ss2 <- update(object = loo_ss, draws = fake_posterior, data = fake_data, observations = 600, r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss2_subset0 <- update(object = true_loo_ss, observations = loo_ss2, r_eff = rep(1, nrow(fake_data))))
  expect_true(setequal(obs_idx(loo_ss2), obs_idx(loo_ss2_subset0)))
  expect_true(identical(obs_idx(loo_ss2), obs_idx(loo_ss2_subset0)))
  expect_true(identical(loo_ss2$diagnostic, loo_ss2_subset0$diagnostic))

  # Add tests for changing approx variable
  expect_silent(loo_ss_lpd <- update(object = loo_ss, draws = fake_posterior, data = fake_data, loo_approximation = "lpd", r_eff = rep(1, nrow(fake_data))))
  expect_failure(expect_equal(loo_ss_lpd$loo_subsampling$elpd_loo_approx, loo_ss$loo_subsampling$elpd_loo_approx))
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss)[2])
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss_lpd$pointwise)[1])
  expect_length(loo_ss_lpd$diagnostics$pareto_k, 500)
  expect_length(loo_ss_lpd$diagnostics$n_eff, 500)
  expect_failure(expect_equal(loo_ss_lpd$estimates[1, "subsampling SE"], loo_ss$estimates[1, "subsampling SE"]))
  expect_failure(expect_equal(loo_ss_lpd$estimates[3, "subsampling SE"], loo_ss$estimates[3, "subsampling SE"]))

})







context("loo_subsampling_approximations")

geterate_test_elpd_dataset <- function() {
  N <- 10; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- draws <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)

  list(fake_posterior = fake_posterior, fake_data = fake_data)
}

test_elpd_loo_approximation <- function(cores) {
  set.seed(123)
  test_data <- geterate_test_elpd_dataset()
  fake_posterior <- test_data$fake_posterior
  fake_data <- test_data$fake_data

  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  # Compute plpd approximation
  expect_silent(pi_vals <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "plpd", cores = cores))
  # Compute it manually
  point <- mean(fake_posterior)
  llik <- dbinom(fake_data$y, size = fake_data$K, prob = point, log = TRUE)
  abs_lliks <- abs(llik)
  man_elpd_loo_approximation <- abs_lliks/sum(abs_lliks)
  expect_equal(abs(pi_vals)/sum(abs(pi_vals)), man_elpd_loo_approximation, tol = 0.00001)

  # Compute lpd approximation
  expect_silent(pi_vals <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "lpd", cores = cores))
  # Compute it manually
  llik <- numeric(10)
  for(i in seq_along(fake_data$y)){
    llik[i] <- loo:::logMeanExp(dbinom(fake_data$y[i], size = fake_data$K, prob = fake_posterior, log = TRUE))
  }
  abs_lliks <- abs(llik)
  man_approx_loo_variable <- abs_lliks/sum(abs_lliks)
  expect_equal(abs(pi_vals)/sum(abs(pi_vals)), man_approx_loo_variable, tol = 0.00001)

  # Compute waic approximation
  expect_silent(pi_vals_waic <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "waic", cores = cores))
  expect_true(all(pi_vals > pi_vals_waic))
  expect_true(sum(pi_vals) - sum(pi_vals_waic) < 1)

  # Compute tis approximation
  expect_silent(pi_vals_tis <- loo:::elpd_loo_approximation(.llfun = llfun_test,
                                                            data = fake_data,
                                                            draws = fake_posterior,
                                                            loo_approximation = "tis",
                                                            loo_approximation_draws = 100,
                                                            cores = cores))
  expect_true(all(pi_vals > pi_vals_tis))
  expect_true(sum(pi_vals) - sum(pi_vals_tis) < 1)
}

test_that("elpd_loo_approximation works as expected", {
  test_elpd_loo_approximation(1)
})

test_that("elpd_loo_approximation with multiple cores", {
  test_elpd_loo_approximation(2)
})

test_that("Test loo_approximation_draws", {


  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- draws <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(res1 <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "waic", loo_approximation_draws = NULL, cores = 1))
  expect_silent(res2 <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "waic", loo_approximation_draws = 10, cores = 1))
  expect_silent(res3 <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior[1:10*100,], loo_approximation = "waic", loo_approximation_draws = NULL, cores = 1))
  expect_silent(res4 <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior[1:10*100,, drop = FALSE], loo_approximation = "waic", loo_approximation_draws = NULL, cores = 1))
  expect_failure(expect_equal(res1, res3))
  expect_equal(res2, res3)

  expect_silent(loo_ss1 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 100, loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss2 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 100, loo_approximation = "plpd", loo_approximation_draws = 10, r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss3 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 100, loo_approximation = "plpd", loo_approximation_draws = 31, r_eff = rep(1, nrow(fake_data))))
  expect_error(loo_ss4 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = 100, loo_approximation = "plpd", loo_approximation_draws = 3100, r_eff = rep(1, nrow(fake_data))))

  expect_equal(names(loo_ss1$loo_subsampling), c("elpd_loo_approx", "loo_approximation", "loo_approximation_draws", "estimator", ".llfun", ".llgrad", ".llhess", "data_dim", "ndraws"))
  expect_null(loo_ss1$loo_subsampling$loo_approximation_draws)
  expect_equal(loo_ss2$loo_subsampling$loo_approximation_draws, 10L)
  expect_equal(loo_ss3$loo_subsampling$loo_approximation_draws, 31L)

})



test_that("waic using delta method and gradient", {


  if (FALSE){
    # Code to generate testdata - saved and loaded to avoid dependency of mvtnorm
    set.seed(123)
    N <- 400; beta <- c(1,2); X_full <- matrix(rep(1,N), ncol = 1); X_full <- cbind(X_full, runif(N)); S <- 1000
    y_full <- rnorm(n = N, mean = X_full%*%beta, sd = 1)
    X <- X_full; y <- y_full
    Lambda_0 <- diag(length(beta)); mu_0 <- c(0,0)
    b_hat <- solve(t(X)%*%X)%*%t(X)%*%y
    mu_n <- solve(t(X)%*%X)%*%(t(X)%*%X%*%b_hat + Lambda_0%*%mu_0)
    Lambda_n <- t(X)%*%X + Lambda_0
    # Uncomment row below when running. Commented out to remove CHECK warnings
    # fake_posterior <- mvtnorm::rmvnorm(n = S, mean = mu_n, sigma = solve(Lambda_n))
    colnames(fake_posterior) <- c("a", "b")
    fake_data <- data.frame(y, X)
    save(fake_posterior, fake_data, file = test_path("data-for-tests/normal_reg_waic_test_example.rda"))
  } else {
    load(file = test_path("data-for-tests/normal_reg_waic_test_example.rda"))
  }

  .llfun <- function(data_i, draws) {
    # data_i: ith row of fdata (fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dnorm(data_i$y, mean = draws[, c("a", "b")] %*% t(as.matrix(data_i[, c("X1", "X2")])), sd = 1, log = TRUE)
  }

  .llgrad <- function(data_i, draws) {
    x_i <- data_i[, "X2"]
    gr <- cbind(data_i$y - draws[,"a"] - draws[,"b"]*x_i,
                (data_i$y - draws[,"a"] - draws[,"b"]*x_i) * x_i)
    colnames(gr) <- c("a", "b")
    gr
  }

  fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))

  expect_silent(approx_loo_waic <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic"))
  expect_silent(approx_loo_waic_delta <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic_grad", .llgrad = .llgrad))
  expect_silent(approx_loo_waic_delta_diag <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic_grad_marginal", .llgrad = .llgrad))

  # Test that the approaches should not deviate too much
  diff_waic_delta <- mean(approx_loo_waic - approx_loo_waic_delta)
  diff_waic_delta_diag <- mean(approx_loo_waic - approx_loo_waic_delta_diag)
  expect_equal(approx_loo_waic,approx_loo_waic_delta_diag, tol = 0.1)
  expect_equal(approx_loo_waic,approx_loo_waic_delta, tol = 0.01)

  # Test usage in subsampling_loo
  expect_silent(loo_ss_waic <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic", observations = 50, llgrad = .llgrad))
  expect_silent(loo_ss_waic_delta <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_grad", observations = 50, llgrad = .llgrad))
  expect_silent(loo_ss_waic_delta_marginal <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_grad_marginal", observations = 50, llgrad = .llgrad))
  expect_silent(loo_ss_plpd <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "plpd", observations = 50, llgrad = .llgrad))
  expect_error(loo_ss_waic_delta <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_grad", observations = 50))
})

test_that("waic using delta 2nd order method", {


  if (FALSE){
    # Code to generate testdata - saved and loaded to avoid dependency of MCMCPack
    set.seed(123)
    N <- 100; beta <- c(1,2); X_full <- matrix(rep(1,N), ncol = 1); X_full <- cbind(X_full, runif(N)); S <- 1000
    y_full <- rnorm(n = N, mean = X_full%*%beta, sd = 0.5)
    X <- X_full; y <- y_full
    # Uncomment row below when running. Commented out to remove CHECK warnings
    # fake_posterior <- MCMCpack::MCMCregress(y~x, data = data.frame(y = y,x=X[,2]), thin = 10, mcmc = 10000) # Because Im lazy
    fake_posterior <- as.matrix(fake_posterior)
    fake_posterior[,"sigma2"] <- sqrt(fake_posterior[,"sigma2"])
    colnames(fake_posterior) <- c("a", "b", "sigma")
    fake_data <- data.frame(y, X)
    save(fake_posterior, fake_data, file = test_path("data-for-tests/normal_reg_waic_test_example2.rda"), compression_level = 9)
  } else {
    load(file = test_path("data-for-tests/normal_reg_waic_test_example2.rda"))
  }

  .llfun <- function(data_i, draws) {
    # data_i: ith row of fdata (data_i <- fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dnorm(data_i$y, mean = draws[, c("a", "b")] %*% t(as.matrix(data_i[, c("X1", "X2")])), sd = draws[, c("sigma")], log = TRUE)
  }

  .llgrad <- function(data_i, draws) {
    sigma <- draws[,"sigma"]
    sigma2 <- sigma^2
    b <- draws[,"b"]
    a <- draws[,"a"]
    x_i <- unlist(data_i[, c("X1", "X2")])
    e <- (data_i$y - draws[,"a"] * x_i[1] - draws[,"b"] * x_i[2])

    gr <- cbind(e * x_i[1] / sigma2,
                e * x_i[2] / sigma2,
                - 1 / sigma + e^2 / (sigma2 * sigma))
    colnames(gr) <- c("a", "b", "sigma")
    gr
  }

  .llhess <- function(data_i, draws) {
    hess_array <- array(0, dim = c(ncol(draws), ncol(draws), nrow(draws)), dimnames = list(colnames(draws),colnames(draws),NULL))
    sigma <- draws[,"sigma"]
    sigma2 <- sigma^2
    sigma3 <- sigma2*sigma
    b <- draws[,"b"]
    a <- draws[,"a"]
    x_i <- unlist(data_i[, c("X1", "X2")])
    e <- (data_i$y - draws[,"a"] * x_i[1] - draws[,"b"] * x_i[2])

    hess_array[1,1,] <- - x_i[1]^2 / sigma2
    hess_array[1,2,] <- hess_array[2,1,] <- - x_i[1] * x_i[2] / sigma2
    hess_array[2,2,] <- - x_i[2]^2 / sigma2
    hess_array[3,1,] <- hess_array[1,3,] <- -2 * x_i[1] * e / sigma3
    hess_array[3,2,] <- hess_array[2,3,] <- -2 * x_i[2] * e / sigma3
    hess_array[3,3,] <- 1 / sigma2 - 3 * e^2 / (sigma2^2)
    hess_array
  }

  #data <- fake_data
  fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))
  #draws <- fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))

  expect_silent(approx_loo_waic <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic"))
  expect_silent(approx_loo_waic_delta <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic_grad", .llgrad = .llgrad))
  expect_silent(approx_loo_waic_delta2 <- loo:::elpd_loo_approximation(.llfun, data = fake_data, draws = fake_posterior, cores = 1, loo_approximation = "waic_hess", .llgrad = .llgrad, .llhess = .llhess))

  # Test that the approaches should not deviate too much
  expect_equal(approx_loo_waic,approx_loo_waic_delta2, tol = 0.01)
  expect_equal(approx_loo_waic,approx_loo_waic_delta, tol = 0.01)

  expect_silent(test_loo_ss_waic <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic", observations = 50, llgrad = .llgrad))
  expect_error(test_loo_ss_delta2 <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_hess", observations = 50, llgrad = .llgrad))
  expect_silent(test_loo_ss_delta2 <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_hess", observations = 50, llgrad = .llgrad, llhess = .llhess))
  expect_silent(test_loo_ss_delta <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "waic_grad", observations = 50, llgrad = .llgrad))
  expect_silent(test_loo_ss_point <- loo_subsample(x = .llfun, data = fake_data, draws = fake_posterior, cores = 1, r_eff = rep(1, nrow(fake_data)), loo_approximation = "plpd", observations = 50, llgrad = .llgrad))
})







context("loo_subsampling_estimation")

test_that("whhest works as expected", {


  N <- 100
  m <- 10
  z <- rep(1/N, m)
  y <- 1:10
  m_i <- rep(1,m)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(whe$y_hat_ppz, 550)
  man_var <- (sum((whe$y_hat_ppz - y/z)^2)/(m-1))/m
  expect_equal(whe$v_hat_y_ppz, man_var)
  z <- 1:10/(sum(1:10)*10)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(whe$y_hat_ppz, 550)
  expect_equal(whe$v_hat_y_ppz, 0)

  # School book example
  # https://newonlinecourses.science.psu.edu/stat506/node/15/
  z <- c(650/15650, 2840/15650, 3200/15650)
  y <- c(420, 1785, 2198)
  m_i <- c(1,1,1)
  N <- 10
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(round(whe$y_hat_ppz, 2), 10232.75, tol = 0)
  expect_equal(whe$v_hat_y_ppz, 73125.74, tol = 0.01)
  # Double check that it is rounding error
  man_var_round <- (sum((round(y/z,2) - 10232.75)^2)) * (1/2) * (1/3)
  expect_equal(man_var_round, 73125.74, tol = 0.001)
  man_var_exact <- (sum((y/z - 10232.75)^2)) * (1/2) * (1/3)
  expect_equal(whe$v_hat_y_ppz, man_var_exact, tol = 0.001)

  # Add test for variance estimation
  N <- 100
  m <- 10
  y <- rep(1:10, 1)
  true_var <- var(rep(y, 10)) * (99)
  z <- rep(1/N, m)
  m_i <- rep(100000, m)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(true_var, whe$hat_v_y_ppz, tol = 0.01)

  # Add tests for m_i
  N <- 100
  y <- rep(1:10, 2)
  m <- length(y)
  z <- rep(1/N, m)
  m_i <- rep(1,m)
  expect_silent(whe1 <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  y <- rep(1:10)
  m <- length(y)
  z <- rep(1/N, m)
  m_i <- rep(2,m)
  expect_silent(whe2 <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(whe1$y_hat_ppz, whe2$y_hat_ppz)
  expect_equal(whe1$v_hat_y_ppz, whe2$v_hat_y_ppz)
  expect_equal(whe1$hat_v_y_ppz, whe1$hat_v_y_ppz)

})


test_that("srs_diff_est works as expected", {


  set.seed(1234)
  N <- 1000
  y_true <- 1:N
  sigma_hat_true <- sqrt(N * sum((y_true - mean(y_true))^2) / length(y_true))
  y_approx <- rnorm(N, y_true, 0.1)
  m <- 100
  sigma_hat <- y_hat <- se_y_hat <- numeric(10000)
  for(i in 1:10000){
    y_idx <- sample(1:N, size = m)
    y <- y_true[y_idx]
    res <- loo:::srs_diff_est(y_approx, y, y_idx)
    y_hat[i] <- res$y_hat
    se_y_hat[i] <- sqrt(res$v_y_hat)
    sigma_hat[i] <- sqrt(res$hat_v_y)
  }
  expect_equal(mean(y_hat), sum(y_true), tol = 0.1)

  in_ki <- y_hat + 2 * se_y_hat > sum(y_true) & y_hat - 2*se_y_hat < sum(y_true)
  expect_equal(mean(in_ki), 0.95, tol = 0.01)

  # Should be  unbiased
  expect_equal(mean(sigma_hat), sigma_hat_true, tol = 0.1)

  m <- N
  y_idx <- sample(1:N, size = m)
  y <- y_true[y_idx]
  res <- loo:::srs_diff_est(y_approx, y, y_idx)
  expect_equal(res$y_hat, 500500, tol = 0.0001)
  expect_equal(res$v_y_hat, 0, tol = 0.0001)
  expect_equal(sqrt(res$hat_v_y), sigma_hat_true, tol = 0.1)

})

test_that("srs_est works as expected", {


  set.seed(1234)
  # Cochran 1976 example Table 2.2

  y <- c(rep(42,23),rep(41,4), 36, 32, 29, 27, 27, 23, 19, 16, 16, 15, 15, 14, 11, 10, 9, 7, 6, 6, 6, 5, 5, 4, 3)
  expect_equal(sum(y), 1471)
  approx_loo <- rep(0L, 676)
  expect_equal(sum(y^2), 54497)
  res <- loo:::srs_est(y = y, approx_loo)
  expect_equal(res$y_hat, 19888, tol = 0.0001)
  expect_equal(res$v_y_hat, 676^2*229*(1-0.074)/50, tol = 0.0001)
  expect_equal(res$hat_v_y, 676 * var(y), tol = 0.0001)

  # Simulation example
  set.seed(1234)
  N <- 1000
  y_true <- 1:N
  sigma_hat_true <- sqrt(N * sum((y_true - mean(y_true))^2) / length(y_true))

  m <- 100
  y_hat <- se_y_hat <- sigma_hat <- numeric(10000)
  for(i in 1:10000){
    y_idx <- sample(1:N, size = m)
    y <- y_true[y_idx]
    res <- loo:::srs_est(y = y, y_approx = y_true)
    y_hat[i] <- res$y_hat
    se_y_hat[i] <- sqrt(res$v_y_hat)
    sigma_hat[i] <- sqrt(res$hat_v_y)
  }
  expect_equal(mean(y_hat), sum(y_true), tol = 0.1)

  in_ki <- y_hat + 2*se_y_hat > sum(y_true) & y_hat - 2*se_y_hat < sum(y_true)
  expect_equal(mean(in_ki), 0.95, tol = 0.01)

  # Should be  unbiased
  expect_equal(mean(sigma_hat), sigma_hat_true, tol = 0.1)

  m <- N
  y_idx <- sample(1:N, size = m)
  y <- y_true[y_idx]
  res <- loo:::srs_est(y, y_true)
  expect_equal(res$y_hat, 500500, tol = 0.0001)
  expect_equal(res$v_y_hat, 0, tol = 0.0001)

})




context("loo_subsampling cases")

test_that("Test loo_subsampling and loo_approx with radon data", {
  skip_on_cran() # avoid going over time limit for tests

  load(test_path("data-for-tests/test_radon_laplace_loo.rda"))
  # Rename to spot variable leaking errors
  llfun_test <- llfun
  log_p_test <- log_p
  log_g_test <- log_q
  draws_test <- draws
  data_test <- data
  rm(llfun, log_p,log_q, draws, data)

  set.seed(134)
  expect_silent(full_loo <- loo(llfun_test, draws = draws_test, data = data_test, r_eff = rep(1, nrow(data_test))))
  expect_s3_class(full_loo, "psis_loo")

  set.seed(134)
  expect_silent(loo_ss <- loo_subsample(x = llfun_test, draws = draws_test, data = data_test, observations = 200, loo_approximation = "plpd", r_eff = rep(1, nrow(data_test))))
  expect_s3_class(loo_ss, "psis_loo_ss")

  set.seed(134)
  expect_silent(loo_ap_ss <- loo_subsample(x = llfun_test, draws = draws_test, data = data_test, log_p = log_p_test, log_g = log_g_test, observations = 200, loo_approximation = "plpd", r_eff = rep(1, nrow(data_test))))
  expect_s3_class(loo_ap_ss, "psis_loo_ss")
  expect_s3_class(loo_ap_ss, "psis_loo_ap")

  expect_silent(loo_ap_ss_full <- loo_subsample(x = llfun_test, log_p = log_p_test, log_g = log_g_test, draws = draws_test, data = data_test, observations = NULL, loo_approximation = "plpd", r_eff = rep(1, nrow(data_test))))
  expect_failure(expect_s3_class(loo_ap_ss_full, "psis_loo_ss"))
  expect_s3_class(loo_ap_ss_full, "psis_loo_ap")

  # Expect similar results
  z <- 2
  expect_lte(loo_ss$estimates["elpd_loo", "Estimate"] - z * loo_ss$estimates["elpd_loo", "subsampling SE"], full_loo$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ss$estimates["elpd_loo", "Estimate"] + z * loo_ss$estimates["elpd_loo", "subsampling SE"], full_loo$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ss$estimates["p_loo", "Estimate"] - z * loo_ss$estimates["p_loo", "subsampling SE"], full_loo$estimates["p_loo", "Estimate"])
  expect_gte(loo_ss$estimates["p_loo", "Estimate"] + z * loo_ss$estimates["p_loo", "subsampling SE"], full_loo$estimates["p_loo", "Estimate"])
  expect_lte(loo_ss$estimates["looic", "Estimate"] - z * loo_ss$estimates["looic", "subsampling SE"], full_loo$estimates["looic", "Estimate"])
  expect_gte(loo_ss$estimates["looic", "Estimate"] + z * loo_ss$estimates["looic", "subsampling SE"], full_loo$estimates["looic", "Estimate"])

  expect_failure(expect_equal(full_loo$estimates["elpd_loo", "Estimate"], loo_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(full_loo$estimates["p_loo", "Estimate"], loo_ss$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(full_loo$estimates["looic", "Estimate"], loo_ss$estimates["looic", "Estimate"], tol = 0.00000001))

  z <- 2
  expect_lte(loo_ap_ss$estimates["elpd_loo", "Estimate"] - z * loo_ap_ss$estimates["elpd_loo", "subsampling SE"], loo_ap_ss_full$estimates["elpd_loo", "Estimate"])
  expect_gte(loo_ap_ss$estimates["elpd_loo", "Estimate"] + z * loo_ap_ss$estimates["elpd_loo", "subsampling SE"], loo_ap_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lte(loo_ap_ss$estimates["p_loo", "Estimate"] - z * loo_ap_ss$estimates["p_loo", "subsampling SE"], loo_ap_ss_full$estimates["p_loo", "Estimate"])
  expect_gte(loo_ap_ss$estimates["p_loo", "Estimate"] + z * loo_ap_ss$estimates["p_loo", "subsampling SE"], loo_ap_ss_full$estimates["p_loo", "Estimate"])
  expect_lte(loo_ap_ss$estimates["looic", "Estimate"] - z * loo_ap_ss$estimates["looic", "subsampling SE"], loo_ap_ss_full$estimates["looic", "Estimate"])
  expect_gte(loo_ap_ss$estimates["looic", "Estimate"] + z * loo_ap_ss$estimates["looic", "subsampling SE"], loo_ap_ss_full$estimates["looic", "Estimate"])

  expect_failure(expect_equal(loo_ap_ss_full$estimates["elpd_loo", "Estimate"], loo_ap_ss$estimates["elpd_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(loo_ap_ss_full$estimates["p_loo", "Estimate"], loo_ap_ss$estimates["p_loo", "Estimate"], tol = 0.00000001))
  expect_failure(expect_equal(loo_ap_ss_full$estimates["looic", "Estimate"], loo_ap_ss$estimates["looic", "Estimate"], tol = 0.00000001))

  # Correct printout
  expect_failure(expect_output(print(full_loo), "Posterior approximation correction used\\."))
  expect_failure(expect_output(print(full_loo), "subsampled log-likelihood\nvalues"))

  expect_failure(expect_output(print(loo_ss), "Posterior approximation correction used\\."))
  expect_output(print(loo_ss), "subsampled log-likelihood\nvalues")

  expect_output(print(loo_ap_ss), "Posterior approximation correction used\\.")
  expect_output(print(loo_ap_ss), "subsampled log-likelihood\nvalues")

  expect_output(print(loo_ap_ss_full), "Posterior approximation correction used\\.")
  expect_failure(expect_output(print(loo_ap_ss_full), "subsampled log-likelihood\nvalues"))

  # Test conversion of objects
  expect_silent(loo_ap_full <- loo:::as.psis_loo.psis_loo(loo_ap_ss_full))
  expect_s3_class(loo_ap_full, "psis_loo_ap")
  expect_silent(loo_ap_full_ss <- loo:::as.psis_loo_ss.psis_loo(loo_ap_full))
  expect_s3_class(loo_ap_full_ss, "psis_loo_ss")
  expect_s3_class(loo_ap_full_ss, "psis_loo_ap")
  expect_silent(loo_ap_full2 <- loo:::as.psis_loo.psis_loo_ss(loo_ap_full_ss))
  expect_s3_class(loo_ap_full2, "psis_loo_ap")
  expect_failure(expect_s3_class(loo_ap_full2, "psis_loo_ss"))
  expect_equal(loo_ap_full2,loo_ap_full)

  # Test update
  set.seed(4712)
  expect_silent(loo_ss2 <- update(loo_ss, draws = draws_test, data = data_test, observations = 1000, r_eff = rep(1, nrow(data_test))))
  expect_gt(dim(loo_ss2)[2], dim(loo_ss)[2])
  expect_gt(dim(loo_ss2$pointwise)[1], dim(loo_ss$pointwise)[1])
  expect_equal(nobs(loo_ss), 200)
  expect_equal(nobs(loo_ss2), 1000)
  for(i in 1:nrow(loo_ss2$estimates)) {
    expect_lt(loo_ss2$estimates[i, "subsampling SE"],
              loo_ss$estimates[i, "subsampling SE"])
  }

  set.seed(4712)
  expect_silent(loo_ap_ss2 <- update(object = loo_ap_ss, draws = draws_test, data = data_test, observations = 2000))
  expect_gt(dim(loo_ap_ss2)[2], dim(loo_ap_ss)[2])
  expect_gt(dim(loo_ap_ss2$pointwise)[1], dim(loo_ap_ss$pointwise)[1])
  expect_equal(nobs(loo_ap_ss), 200)
  expect_equal(nobs(loo_ap_ss2), 2000)
  for(i in 1:nrow(loo_ap_ss2$estimates)) {
    expect_lt(loo_ap_ss2$estimates[i, "subsampling SE"],
              loo_ap_ss$estimates[i, "subsampling SE"])
  }

  expect_equal(round(full_loo$estimates), round(loo_ap_ss_full$estimates))
  expect_failure(expect_equal(full_loo$estimates, loo_ap_ss_full$estimates))
  expect_equal(dim(full_loo), dim(loo_ap_ss_full))
  expect_s3_class(loo_ap_ss_full, "psis_loo_ap")

})


test_that("Test the vignette", {
  skip_on_cran()


  # NOTE
  # If any of these test fails, the vignette probably needs to be updated

  if (FALSE) {
    # Generate vignette test case
    library("rstan")
    stan_code <- "
    data {
      int<lower=0> N;             // number of data points
      int<lower=0> P;             // number of predictors (including intercept)
      matrix[N,P] X;              // predictors (including 1s for intercept)
      int<lower=0,upper=1> y[N];  // binary outcome
    }
    parameters {
      vector[P] beta;
    }
    model {
      beta ~ normal(0, 1);
      y ~ bernoulli_logit(X * beta);
    }
    "
    # logistic <- function(x) {1  / (1 + exp(-x))}
    # logit <- function(x) {log(x) - log(1-x)}
    llfun_logistic <- function(data_i, draws) {
      x_i <- as.matrix(data_i[, which(grepl(colnames(data_i), pattern = "X")), drop=FALSE])
      y_pred <- draws %*% t(x_i)
      dbinom(x = data_i$y, size = 1, prob = 1  / (1 + exp(-y_pred)), log = TRUE)
    }

    # Prepare data
    url <- "http://stat.columbia.edu/~gelman/arm/examples/arsenic/wells.dat"
    wells <- read.table(url)
    wells$dist100 <- with(wells, dist / 100)
    X <- model.matrix(~ dist100 + arsenic, wells)
    standata <- list(y = wells$switch, X = X, N = nrow(X), P = ncol(X))

    # Fit model
    set.seed(4711)
    fit_1 <- stan(model_code = stan_code, data = standata, seed = 4711)
    print(fit_1, pars = "beta")

    parameter_draws <- extract(fit_1)$beta
    stan_df <- as.data.frame(standata)
    loo_i(1, llfun_logistic, data = stan_df, draws = parameter_draws)

    sm <- stan_model(model_code = stan_code)
    set.seed(4711)
    fit_laplace <- optimizing(sm, data = standata, draws = 2000, seed = 42)
    parameter_draws_laplace <- fit_laplace$theta_tilde
    log_p <- fit_laplace$log_p # The log density of the posterior
    log_g <- fit_laplace$log_g # The log density of the approximation

    # For comparisons
    standata$X[, "arsenic"] <- log(standata$X[, "arsenic"])
    stan_df2 <- as.data.frame(standata)
    set.seed(4711)
    fit_2 <- stan(fit = fit_1, data = standata, seed = 4711)
    parameter_draws_2 <- extract(fit_2)$beta

    save(llfun_logistic,
         stan_df, stan_df2,
         parameter_draws, parameter_draws_laplace, parameter_draws_2,
         log_p, log_g,
         file = test_path("data-for-tests/loo_subsample_vignette.rda"), compression_level = 9)

  } else {
    load(test_path("data-for-tests/loo_subsample_vignette.rda"))
  }

  set.seed(4711)
  expect_no_warning(looss_1 <- loo_subsample(llfun_logistic, draws = parameter_draws, data = stan_df, observations = 100))
  expect_output(print(looss_1), "Computed from 4000 by 100 subsampled log-likelihood")
  expect_output(print(looss_1), "values from 3020 total observations.")
  expect_output(print(looss_1), "MCSE and ESS estimates assume independent draws")
  expect_output(print(looss_1), "elpd_loo  -1968.5 15.6            0.3")
  expect_output(print(looss_1), "p_loo         3.1  0.1            0.4")
  expect_s3_class(looss_1, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(looss_1b <- update(looss_1, draws = parameter_draws, data = stan_df, observations = 200))
  expect_output(print(looss_1b), "Computed from 4000 by 200 subsampled log-likelihood")
  expect_output(print(looss_1b), "values from 3020 total observations.")
  expect_output(print(looss_1b), "MCSE and ESS estimates assume independent draws")
  expect_output(print(looss_1b), "elpd_loo  -1968.3 15.6            0.2")
  expect_output(print(looss_1b), "p_loo         3.2  0.1            0.4")
  expect_s3_class(looss_1b, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(looss_2 <- loo_subsample(x = llfun_logistic, draws = parameter_draws, data = stan_df, observations = 100, estimator = "hh_pps", loo_approximation = "lpd", loo_approximation_draws = 100))
  expect_output(print(looss_2), "Computed from 4000 by 100 subsampled log-likelihood")
  expect_output(print(looss_2), "values from 3020 total observations.")
  expect_output(print(looss_2), "MCSE and ESS estimates assume independent draws")
  # Currently failing
  # expect_output(print(looss_2), "elpd_loo  -1968.9 15.4            0.5")
  # expect_output(print(looss_2), "p_loo         3.5  0.2            0.5")
  expect_s3_class(looss_2, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(aploo_1 <- loo_approximate_posterior(llfun_logistic, draws = parameter_draws_laplace, data = stan_df, log_p = log_p, log_g = log_g))
  expect_output(print(aploo_1), "Computed from 2000 by 3020 log-likelihood matrix")
  expect_output(print(aploo_1), "MCSE and ESS estimates assume independent draws")
  expect_output(print(aploo_1), "elpd_loo  -1968.4 15.6")
  expect_output(print(aploo_1), "p_loo         3.2  0.2")
  expect_output(print(aploo_1), "Posterior approximation correction used.")
  expect_output(print(aploo_1), "All Pareto k estimates are good")
  expect_equal(length(pareto_k_ids(aploo_1,threshold=0.5)), 31)
  expect_s3_class(aploo_1, c("psis_loo_ap", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(looapss_1 <- loo_subsample(llfun_logistic, draws = parameter_draws_laplace, data = stan_df, log_p = log_p, log_g = log_g, observations = 100))
  expect_output(print(looapss_1), "Computed from 2000 by 100 subsampled log-likelihood")
  expect_output(print(looapss_1), "MCSE and ESS estimates assume independent draws")
  expect_output(print(looapss_1), "values from 3020 total observations.")
  expect_output(print(looapss_1), "elpd_loo  -1968.2 15.6            0.4")
  expect_output(print(looapss_1), "p_loo         2.9  0.1            0.5")
  expect_output(print(looapss_1), "All Pareto k estimates are good")
  expect_equal(length(pareto_k_ids(looapss_1,threshold=0.5)), 3)

  # Loo compare
  set.seed(4711)
  expect_no_warning(looss_1 <- loo_subsample(llfun_logistic, draws = parameter_draws, data = stan_df, observations = 100))
  set.seed(4712)
  expect_no_warning(looss_2 <- loo_subsample(x = llfun_logistic, draws = parameter_draws_2, data = stan_df2, observations = 100))
  expect_output(print(looss_2), "Computed from 4000 by 100 subsampled log-likelihood")
  expect_output(print(looss_2), "MCSE and ESS estimates assume independent draws")
  expect_output(print(looss_2), "values from 3020 total observations.")
  expect_output(print(looss_2), "elpd_loo  -1952.0 16.2            0.2")
  expect_output(print(looss_2), "p_loo         2.6  0.1            0.3")

  expect_warning(comp <- loo_compare(looss_1, looss_2), "Different subsamples in 'model2' and 'model1'. Naive diff SE is used.")
  expect_output(print(comp), "model1 16.5      22.5     0.4")

  set.seed(4712)
  expect_no_warning(looss_2_m <- loo_subsample(x = llfun_logistic, draws = parameter_draws_2, data = stan_df2, observations = looss_1))
  expect_message(looss_2_m <- suppressWarnings(loo_subsample(x = llfun_logistic, draws = parameter_draws_2, data = stan_df2, observations = obs_idx(looss_1))),
                 "Simple random sampling with replacement assumed.")

  expect_silent(comp <- loo_compare(looss_1, looss_2_m))
  expect_output(print(comp), "model1 16.1       4.4     0.1")

  set.seed(4712)
  expect_no_warning(looss_1 <- update(looss_1, draws = parameter_draws, data = stan_df, observations = 200))
  expect_no_warning(looss_2_m <- update(looss_2_m, draws = parameter_draws_2, data = stan_df2, observations = looss_1))
  expect_silent(comp2 <- loo_compare(looss_1, looss_2_m))
  expect_output(print(comp2), "model1 16.3       4.4     0.1")

  expect_no_warning(looss_2_full <- loo(x = llfun_logistic, draws = parameter_draws_2, data = stan_df2))
  expect_message(comp3 <- loo_compare(x = list(looss_1, looss_2_full)),
                 "Estimated elpd_diff using observations included in loo calculations for all models.")
  expect_output(print(comp3), "model1 16.5       4.4     0.3")

})


context("loo_compare_subsample")

test_that("loo_compare_subsample", {
  skip_on_cran() # to get under cran check time limit

  set.seed(123)
  N <- 1000
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  x3 <- rnorm(N)
  sigma <- 2
  y <- rnorm(N, 1 + 2*x1 - 2*x2 - 1*x3, sd = sigma)
  X <- cbind("x0" = rep(1, N), x1, x2, x3)

  # Generate samples from posterior
  samples_blin <- function(X, y, sigma, draws = 1000){
    XtX <- t(X)%*%X
    b_hat <- solve(XtX) %*% (t(X) %*% y)
    Lambda_n = XtX + diag(ncol(X))
    mu_n <- solve(Lambda_n) %*% (XtX %*% b_hat + diag(ncol(X)) %*% rep(0,ncol(X)))
    L <- t(chol(sigma^2 * solve(Lambda_n)))
    draws_mat <- matrix(0, ncol=ncol(X), nrow = draws)
    for(i in 1:draws){
      z <- rnorm(length(mu_n))
      draws_mat[i,] <- L %*% z + mu_n
    }
    draws_mat
  }

  fake_posterior1 <- samples_blin(X[,1:2], y, sigma, draws = 1000)
  fake_posterior2 <- samples_blin(X[,1:3], y, sigma, draws = 1000)
  fake_posterior3 <- samples_blin(X, y, sigma, draws = 1000)

  fake_data1 <- data.frame(y, X[,1:2])
  fake_data2 <- data.frame(y, X[,1:3])
  fake_data3 <- data.frame(y, X)

  llfun_test <- function(data_i, draws) {
    dnorm(x = data_i$y, mean = draws %*% t(data_i[,-1, drop=FALSE]), sd = sigma, log = TRUE)
  }

  expect_silent(l1 <- loo(llfun_test, data = fake_data1, draws = fake_posterior1, r_eff = rep(1, N)))
  expect_silent(l2 <- loo(llfun_test, data = fake_data2, draws = fake_posterior2, r_eff = rep(1, N)))
  expect_silent(l3 <- loo(llfun_test, data = fake_data3, draws = fake_posterior3, r_eff = rep(1, N)))

  expect_silent(lss1 <- loo_subsample(llfun_test, data = fake_data1, draws = fake_posterior1, observations = 100, r_eff = rep(1, N)))
  expect_silent(lss2 <- loo_subsample(llfun_test, data = fake_data2, draws = fake_posterior2, observations = 100, r_eff = rep(1, N)))
  expect_silent(lss3 <- loo_subsample(llfun_test, data = fake_data3, draws = fake_posterior3, observations = 100, r_eff = rep(1, N)))
  expect_silent(lss2o1 <- loo_subsample(llfun_test, data = fake_data2, draws = fake_posterior2, observations = lss1, r_eff = rep(1, N)))
  expect_silent(lss3o1 <- loo_subsample(llfun_test, data = fake_data3, draws = fake_posterior3, observations = lss1, r_eff = rep(1, N)))
  expect_silent(lss2hh <- loo_subsample(llfun_test, data = fake_data2, draws = fake_posterior2, observations = 100, estimator = "hh_pps", r_eff = rep(1, N)))

  expect_warning(lcss <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2, lss3)))
  expect_warning(lcss2 <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2, lss3o1)))
  expect_silent(lcsso <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2o1, lss3o1)))
  expect_warning(lcssohh <- loo:::loo_compare.psis_loo_ss_list(x = list(lss1, lss2hh, lss3o1)))
  expect_message(lcssf1 <- loo:::loo_compare.psis_loo_ss_list(x = list(loo:::as.psis_loo_ss.psis_loo(l1), lss2o1, lss3o1)))
  expect_message(lcssf2 <- loo:::loo_compare.psis_loo_ss_list(x = list(loo:::as.psis_loo_ss.psis_loo(l1), lss2o1, loo:::as.psis_loo_ss.psis_loo(l3))))

  expect_equal(lcss[,1], lcsso[,1], tolerance = 1)
  expect_equal(lcss2[,1], lcsso[,1], tolerance = 1)
  expect_equal(lcssohh[,1], lcsso[,1], tolerance = 1)
  expect_equal(lcssf1[,1], lcsso[,1], tolerance = 1)
  expect_equal(lcssf2[,1], lcsso[,1], tolerance = 1)

  expect_gt(lcss[,2][2], lcsso[,2][2])
  expect_gt(lcss[,2][3], lcsso[,2][3])
  expect_gt(lcss2[,2][2], lcsso[,2][2])
  expect_equal(lcss2[,2][3], lcsso[,2][3])
  expect_gt(lcssohh[,2][2], lcsso[,2][2])
  expect_equal(lcssohh[,2][3], lcsso[,2][3])

  expect_silent(lcss2m <- loo:::loo_compare.psis_loo_ss_list(x = list(lss2o1, lss3o1)))
  expect_equal(unname(lcss2m[,]), unname(lcsso[1:2,]))

  expect_warning(lcssapi <- loo_compare(lss1, lss2, lss3))
  expect_equal(lcssapi, lcss)
  expect_warning(lcssohhapi <- loo_compare(lss1, lss2hh, lss3o1))
  expect_equal(lcssohhapi, lcssohh)
  expect_silent(lcss2mapi <- loo_compare(lss2o1, lss3o1))
  expect_equal(lcss2mapi, lcss2m)

})

context("subsample with tis, sis")

test_that("Test 'tis' and 'sis'", {
  skip_on_cran()

  set.seed(123)
  N <- 1000; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- draws <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  expect_silent(loo_ss_full <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 1000,
                                loo_approximation = "plpd",
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_plpd <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "plpd",
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_tis_S1000 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "tis",
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_tis_S100 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "tis",
                                loo_approximation_draws = 100,
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_tis_S10 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "tis",
                                loo_approximation_draws = 10,
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_sis_S1000 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "sis",
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_sis_S100 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "sis",
                                loo_approximation_draws = 100,
                                r_eff = rep(1, nrow(fake_data))))
  expect_silent(loo_ss_sis_S10 <-
                  loo_subsample(x = llfun_test,
                                draws = fake_posterior,
                                data = fake_data,
                                observations = 100,
                                loo_approximation = "sis",
                                loo_approximation_draws = 10,
                                r_eff = rep(1, nrow(fake_data))))


  SEs <- 4
  expect_gt(loo_ss_tis_S1000$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_tis_S1000$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_tis_S1000$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_tis_S1000$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_gt(loo_ss_tis_S100$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_tis_S100$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_tis_S100$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_tis_S100$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_gt(loo_ss_tis_S10$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_tis_S10$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_tis_S10$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_tis_S10$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])

  expect_gt(loo_ss_sis_S1000$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_sis_S1000$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_sis_S1000$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_sis_S1000$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_gt(loo_ss_sis_S100$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_sis_S100$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_sis_S100$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_sis_S100$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_gt(loo_ss_sis_S10$estimates["elpd_loo", "Estimate"] + SEs*loo_ss_sis_S10$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
  expect_lt(loo_ss_sis_S10$estimates["elpd_loo", "Estimate"] - SEs*loo_ss_sis_S10$estimates["elpd_loo", "subsampling SE"], loo_ss_full$estimates["elpd_loo", "Estimate"])
})
