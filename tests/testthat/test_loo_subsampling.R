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
                    loo_ss$subsamling_loo$elpd_loo_approx[loo_ss$pointwise[, "idx"]])

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
  expect_silent(loo_ss2 <- loo_subsample(x = llfun_test, draws = fake_posterior, data = fake_data, observations = obs_idx(loo_ss), loo_approximation = "plpd", r_eff = rep(1, nrow(fake_data))))
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

  skip("Add tests on supplying observation vector")
  skip("implement loo_subsampling.matrix and loo_subsampling.array.") # TODO
  skip("test all approximations and sample approaches.") # TODO

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

  # Check consistency
  expect_equivalent(loo_ss$pointwise[, "elpd_loo_approx"],
                    loo_ss$subsamling_loo$elpd_loo_approx[loo_ss$pointwise[, "idx"]])

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

  # Add tests for changing approx variable
  expect_silent(loo_ss_lpd <- update(object = loo_ss, draws = fake_posterior, data = fake_data, loo_approximation = "lpd", r_eff = rep(1, nrow(fake_data))))
  expect_failure(expect_equal(loo_ss_lpd$subsamling_loo$elpd_loo_approx, loo_ss$subsamling_loo$elpd_loo_approx))
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss)[2])
  expect_equal(dim(loo_ss_lpd)[2], dim(loo_ss_lpd$pointwise)[1])
  expect_length(loo_ss_lpd$diagnostics$pareto_k, 500)
  expect_length(loo_ss_lpd$diagnostics$n_eff, 500)
  expect_failure(expect_equal(loo_ss_lpd$estimates[1, "subsampling SE"], loo_ss$estimates[1, "subsampling SE"]))
  expect_failure(expect_equal(loo_ss_lpd$estimates[3, "subsampling SE"], loo_ss$estimates[3, "subsampling SE"]))

})




context("loo_subsampling_approximations")

test_that("elpd_loo_approximation works as expected", {
  set.seed(123)
  N <- 10; K <- 10; S <- 1000; a0 <- 3; b0 <- 2
  p <- 0.7
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  fake_posterior <- draws <- as.matrix(rbeta(S, a, b))
  fake_data <- data.frame(y,K)
  rm(N, K, S, a0, b0, p, y, a, b)
  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  # Compute plpd approximation
  expect_silent(pi_vals <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "plpd", cores = 1))
  # Compute it manually
  point <- mean(fake_posterior)
  llik <- dbinom(fake_data$y, size = fake_data$K, prob = point, log = TRUE)
  abs_lliks <- abs(llik)
  man_elpd_loo_approximation <- abs_lliks/sum(abs_lliks)
  expect_equal(abs(pi_vals)/sum(abs(pi_vals)), man_elpd_loo_approximation, tol = 0.00001)

  # Compute lpd approximation
  expect_silent(pi_vals <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "lpd", cores = 1))
  # Compute it manually
  llik <- numeric(10)
  for(i in seq_along(fake_data$y)){
    llik[i] <- loo:::logMeanExp(dbinom(fake_data$y[i], size = fake_data$K, prob = fake_posterior, log = TRUE))
  }
  abs_lliks <- abs(llik)
  man_approx_loo_variable <- abs_lliks/sum(abs_lliks)
  expect_equal(abs(pi_vals)/sum(abs(pi_vals)), man_approx_loo_variable, tol = 0.00001)

  # Compute waic approximation
  expect_silent(pi_vals_waic <- loo:::elpd_loo_approximation(.llfun = llfun_test, data = fake_data, draws = fake_posterior, loo_approximation = "waic", cores = 1))
  expect_true(all(pi_vals > pi_vals_waic))
  expect_true(sum(pi_vals) - sum(pi_vals_waic) < 1)

  skip("Compute delta methods waic") # TODO
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




