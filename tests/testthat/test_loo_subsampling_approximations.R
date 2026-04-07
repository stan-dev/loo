options(mc.cores = 1)

generate_test_elpd_dataset <- function() {
  N <- 10
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

  list(fake_posterior = fake_posterior, fake_data = fake_data)
}

test_elpd_loo_approximation <- function(cores) {
  set.seed(123)
  test_data <- generate_test_elpd_dataset()
  fake_posterior <- test_data$fake_posterior
  fake_data <- test_data$fake_data

  llfun_test <- function(data_i, draws) {
    dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
  }

  # Compute plpd approximation
  expect_silent(
    pi_vals <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "plpd",
      cores = cores
    )
  )
  # Compute it manually
  point <- mean(fake_posterior)
  llik <- dbinom(fake_data$y, size = fake_data$K, prob = point, log = TRUE)
  abs_lliks <- abs(llik)
  man_elpd_loo_approximation <- abs_lliks / sum(abs_lliks)
  expect_equal(
    abs(pi_vals) / sum(abs(pi_vals)),
    man_elpd_loo_approximation,
    tolerance = 0.00001
  )

  # Compute lpd approximation
  expect_silent(
    pi_vals <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "lpd",
      cores = cores
    )
  )
  # Compute it manually
  llik <- numeric(10)
  for (i in seq_along(fake_data$y)) {
    llik[i] <- loo:::logMeanExp(dbinom(
      fake_data$y[i],
      size = fake_data$K,
      prob = fake_posterior,
      log = TRUE
    ))
  }
  abs_lliks <- abs(llik)
  man_approx_loo_variable <- abs_lliks / sum(abs_lliks)
  expect_equal(
    abs(pi_vals) / sum(abs(pi_vals)),
    man_approx_loo_variable,
    tolerance = 0.00001
  )

  # Compute waic approximation
  expect_silent(
    pi_vals_waic <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "waic",
      cores = cores
    )
  )
  expect_true(all(pi_vals > pi_vals_waic))
  expect_true(sum(pi_vals) - sum(pi_vals_waic) < 1)

  # Compute tis approximation
  expect_silent(
    pi_vals_tis <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "tis",
      loo_approximation_draws = 100,
      cores = cores
    )
  )
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
    res1 <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "waic",
      loo_approximation_draws = NULL,
      cores = 1
    )
  )
  expect_silent(
    res2 <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior,
      loo_approximation = "waic",
      loo_approximation_draws = 10,
      cores = 1
    )
  )
  expect_silent(
    res3 <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior[1:10 * 100, ],
      loo_approximation = "waic",
      loo_approximation_draws = NULL,
      cores = 1
    )
  )
  expect_silent(
    res4 <- loo:::elpd_loo_approximation(
      .llfun = llfun_test,
      data = fake_data,
      draws = fake_posterior[1:10 * 100, , drop = FALSE],
      loo_approximation = "waic",
      loo_approximation_draws = NULL,
      cores = 1
    )
  )
  expect_failure(expect_equal(res1, res3))
  expect_equal(res2, res3)

  expect_silent(
    loo_ss1 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 100,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss2 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 100,
      loo_approximation = "plpd",
      loo_approximation_draws = 10,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_silent(
    loo_ss3 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 100,
      loo_approximation = "plpd",
      loo_approximation_draws = 31,
      r_eff = rep(1, nrow(fake_data))
    )
  )
  expect_error(
    loo_ss4 <- loo_subsample(
      x = llfun_test,
      draws = fake_posterior,
      data = fake_data,
      observations = 100,
      loo_approximation = "plpd",
      loo_approximation_draws = 3100,
      r_eff = rep(1, nrow(fake_data))
    )
  )

  expect_equal(
    names(loo_ss1$loo_subsampling),
    c(
      "elpd_loo_approx",
      "loo_approximation",
      "loo_approximation_draws",
      "estimator",
      ".llfun",
      ".llgrad",
      ".llhess",
      "data_dim",
      "ndraws"
    )
  )
  expect_null(loo_ss1$loo_subsampling$loo_approximation_draws)
  expect_equal(loo_ss2$loo_subsampling$loo_approximation_draws, 10L)
  expect_equal(loo_ss3$loo_subsampling$loo_approximation_draws, 31L)
})


test_that("waic using delta method and gradient", {
  if (FALSE) {
    # Code to generate testdata - saved and loaded to avoid dependency of mvtnorm
    set.seed(123)
    N <- 400
    beta <- c(1, 2)
    X_full <- matrix(rep(1, N), ncol = 1)
    X_full <- cbind(X_full, runif(N))
    S <- 1000
    y_full <- rnorm(n = N, mean = X_full %*% beta, sd = 1)
    X <- X_full
    y <- y_full
    Lambda_0 <- diag(length(beta))
    mu_0 <- c(0, 0)
    b_hat <- solve(t(X) %*% X) %*% t(X) %*% y
    mu_n <- solve(t(X) %*% X) %*% (t(X) %*% X %*% b_hat + Lambda_0 %*% mu_0)
    Lambda_n <- t(X) %*% X + Lambda_0
    # Uncomment row below when running. Commented out to remove CHECK warnings
    # fake_posterior <- mvtnorm::rmvnorm(n = S, mean = mu_n, sigma = solve(Lambda_n))
    colnames(fake_posterior) <- c("a", "b")
    fake_data <- data.frame(y, X)
    save(
      fake_posterior,
      fake_data,
      file = test_path("data-for-tests/normal_reg_waic_test_example.rda")
    )
  } else {
    load(file = test_path("data-for-tests/normal_reg_waic_test_example.rda"))
  }

  .llfun <- function(data_i, draws) {
    # data_i: ith row of fdata (fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dnorm(
      data_i$y,
      mean = draws[, c("a", "b")] %*% t(as.matrix(data_i[, c("X1", "X2")])),
      sd = 1,
      log = TRUE
    )
  }

  .llgrad <- function(data_i, draws) {
    x_i <- data_i[, "X2"]
    gr <- cbind(
      data_i$y - draws[, "a"] - draws[, "b"] * x_i,
      (data_i$y - draws[, "a"] - draws[, "b"] * x_i) * x_i
    )
    colnames(gr) <- c("a", "b")
    gr
  }

  fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))

  expect_silent(
    approx_loo_waic <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic"
    )
  )
  expect_silent(
    approx_loo_waic_delta <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic_grad",
      .llgrad = .llgrad
    )
  )
  expect_silent(
    approx_loo_waic_delta_diag <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic_grad_marginal",
      .llgrad = .llgrad
    )
  )

  # Test that the approaches should not deviate too much
  diff_waic_delta <- mean(approx_loo_waic - approx_loo_waic_delta)
  diff_waic_delta_diag <- mean(approx_loo_waic - approx_loo_waic_delta_diag)
  expect_equal(approx_loo_waic, approx_loo_waic_delta_diag, tolerance = 0.1)
  expect_equal(approx_loo_waic, approx_loo_waic_delta, tolerance = 0.01)

  # Test usage in subsampling_loo
  expect_silent(
    loo_ss_waic <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_silent(
    loo_ss_waic_delta <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_grad",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_silent(
    loo_ss_waic_delta_marginal <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_grad_marginal",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_silent(
    loo_ss_plpd <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "plpd",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_error(
    loo_ss_waic_delta <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_grad",
      observations = 50
    )
  )
})

test_that("waic using delta 2nd order method", {
  if (FALSE) {
    # Code to generate testdata - saved and loaded to avoid dependency of MCMCPack
    set.seed(123)
    N <- 100
    beta <- c(1, 2)
    X_full <- matrix(rep(1, N), ncol = 1)
    X_full <- cbind(X_full, runif(N))
    S <- 1000
    y_full <- rnorm(n = N, mean = X_full %*% beta, sd = 0.5)
    X <- X_full
    y <- y_full
    # Uncomment row below when running. Commented out to remove CHECK warnings
    # fake_posterior <- MCMCpack::MCMCregress(y~x, data = data.frame(y = y,x=X[,2]), thin = 10, mcmc = 10000) # Because Im lazy
    fake_posterior <- as.matrix(fake_posterior)
    fake_posterior[, "sigma2"] <- sqrt(fake_posterior[, "sigma2"])
    colnames(fake_posterior) <- c("a", "b", "sigma")
    fake_data <- data.frame(y, X)
    save(
      fake_posterior,
      fake_data,
      file = test_path("data-for-tests/normal_reg_waic_test_example2.rda"),
      compression_level = 9
    )
  } else {
    load(file = test_path("data-for-tests/normal_reg_waic_test_example2.rda"))
  }

  .llfun <- function(data_i, draws) {
    # data_i: ith row of fdata (data_i <- fake_data[i,, drop=FALSE])
    # draws: entire fake_posterior matrix
    dnorm(
      data_i$y,
      mean = draws[, c("a", "b")] %*% t(as.matrix(data_i[, c("X1", "X2")])),
      sd = draws[, c("sigma")],
      log = TRUE
    )
  }

  .llgrad <- function(data_i, draws) {
    sigma <- draws[, "sigma"]
    sigma2 <- sigma^2
    b <- draws[, "b"]
    a <- draws[, "a"]
    x_i <- unlist(data_i[, c("X1", "X2")])
    e <- (data_i$y - draws[, "a"] * x_i[1] - draws[, "b"] * x_i[2])

    gr <- cbind(
      e * x_i[1] / sigma2,
      e * x_i[2] / sigma2,
      -1 / sigma + e^2 / (sigma2 * sigma)
    )
    colnames(gr) <- c("a", "b", "sigma")
    gr
  }

  .llhess <- function(data_i, draws) {
    hess_array <- array(
      0,
      dim = c(ncol(draws), ncol(draws), nrow(draws)),
      dimnames = list(colnames(draws), colnames(draws), NULL)
    )
    sigma <- draws[, "sigma"]
    sigma2 <- sigma^2
    sigma3 <- sigma2 * sigma
    b <- draws[, "b"]
    a <- draws[, "a"]
    x_i <- unlist(data_i[, c("X1", "X2")])
    e <- (data_i$y - draws[, "a"] * x_i[1] - draws[, "b"] * x_i[2])

    hess_array[1, 1, ] <- -x_i[1]^2 / sigma2
    hess_array[1, 2, ] <- hess_array[2, 1, ] <- -x_i[1] * x_i[2] / sigma2
    hess_array[2, 2, ] <- -x_i[2]^2 / sigma2
    hess_array[3, 1, ] <- hess_array[1, 3, ] <- -2 * x_i[1] * e / sigma3
    hess_array[3, 2, ] <- hess_array[2, 3, ] <- -2 * x_i[2] * e / sigma3
    hess_array[3, 3, ] <- 1 / sigma2 - 3 * e^2 / (sigma2^2)
    hess_array
  }

  #data <- fake_data
  fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))
  #draws <- fake_posterior <- cbind(fake_posterior, runif(nrow(fake_posterior)))

  expect_silent(
    approx_loo_waic <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic"
    )
  )
  expect_silent(
    approx_loo_waic_delta <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic_grad",
      .llgrad = .llgrad
    )
  )
  expect_silent(
    approx_loo_waic_delta2 <- loo:::elpd_loo_approximation(
      .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      loo_approximation = "waic_hess",
      .llgrad = .llgrad,
      .llhess = .llhess
    )
  )

  # Test that the approaches should not deviate too much
  expect_equal(approx_loo_waic, approx_loo_waic_delta2, tolerance = 0.01)
  expect_equal(approx_loo_waic, approx_loo_waic_delta, tolerance = 0.01)

  expect_silent(
    test_loo_ss_waic <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_error(
    test_loo_ss_delta2 <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_hess",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_silent(
    test_loo_ss_delta2 <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_hess",
      observations = 50,
      llgrad = .llgrad,
      llhess = .llhess
    )
  )
  expect_silent(
    test_loo_ss_delta <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "waic_grad",
      observations = 50,
      llgrad = .llgrad
    )
  )
  expect_silent(
    test_loo_ss_point <- loo_subsample(
      x = .llfun,
      data = fake_data,
      draws = fake_posterior,
      cores = 1,
      r_eff = rep(1, nrow(fake_data)),
      loo_approximation = "plpd",
      observations = 50,
      llgrad = .llgrad
    )
  )
})


test_that("whhest works as expected", {
  N <- 100
  m <- 10
  z <- rep(1 / N, m)
  y <- 1:10
  m_i <- rep(1, m)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(whe$y_hat_ppz, 550)
  man_var <- (sum((whe$y_hat_ppz - y / z)^2) / (m - 1)) / m
  expect_equal(whe$v_hat_y_ppz, man_var)
  z <- 1:10 / (sum(1:10) * 10)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(whe$y_hat_ppz, 550)
  expect_equal(whe$v_hat_y_ppz, 0)

  # School book example
  # https://newonlinecourses.science.psu.edu/stat506/node/15/
  z <- c(650 / 15650, 2840 / 15650, 3200 / 15650)
  y <- c(420, 1785, 2198)
  m_i <- c(1, 1, 1)
  N <- 10
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(round(whe$y_hat_ppz, 2), 10232.75, tolerance = 0)
  expect_equal(whe$v_hat_y_ppz, 73125.74, tolerance = 0.01)
  # Double check that it is rounding error
  man_var_round <- (sum((round(y / z, 2) - 10232.75)^2)) * (1 / 2) * (1 / 3)
  expect_equal(man_var_round, 73125.74, tolerance = 0.001)
  man_var_exact <- (sum((y / z - 10232.75)^2)) * (1 / 2) * (1 / 3)
  expect_equal(whe$v_hat_y_ppz, man_var_exact, tolerance = 0.001)

  # Add test for variance estimation
  N <- 100
  m <- 10
  y <- rep(1:10, 1)
  true_var <- var(rep(y, 10)) * (99)
  z <- rep(1 / N, m)
  m_i <- rep(100000, m)
  expect_silent(whe <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  expect_equal(true_var, whe$hat_v_y_ppz, tolerance = 0.01)

  # Add tests for m_i
  N <- 100
  y <- rep(1:10, 2)
  m <- length(y)
  z <- rep(1 / N, m)
  m_i <- rep(1, m)
  expect_silent(whe1 <- loo:::whhest(z = z, m_i = m_i, y = y, N = N))
  y <- rep(1:10)
  m <- length(y)
  z <- rep(1 / N, m)
  m_i <- rep(2, m)
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
  for (i in 1:10000) {
    y_idx <- sample(1:N, size = m)
    y <- y_true[y_idx]
    res <- loo:::srs_diff_est(y_approx, y, y_idx)
    y_hat[i] <- res$y_hat
    se_y_hat[i] <- sqrt(res$v_y_hat)
    sigma_hat[i] <- sqrt(res$hat_v_y)
  }
  expect_equal(mean(y_hat), sum(y_true), tolerance = 0.1)

  in_ki <- y_hat + 2 * se_y_hat > sum(y_true) &
    y_hat - 2 * se_y_hat < sum(y_true)
  expect_equal(mean(in_ki), 0.95, tolerance = 0.01)

  # Should be  unbiased
  expect_equal(mean(sigma_hat), sigma_hat_true, tolerance = 0.1)

  m <- N
  y_idx <- sample(1:N, size = m)
  y <- y_true[y_idx]
  res <- loo:::srs_diff_est(y_approx, y, y_idx)
  expect_equal(res$y_hat, 500500, tolerance = 0.0001)
  expect_equal(res$v_y_hat, 0, tolerance = 0.0001)
  expect_equal(sqrt(res$hat_v_y), sigma_hat_true, tolerance = 0.1)
})

test_that("srs_est works as expected", {
  set.seed(1234)
  # Cochran 1976 example Table 2.2

  y <- c(
    rep(42, 23),
    rep(41, 4),
    36,
    32,
    29,
    27,
    27,
    23,
    19,
    16,
    16,
    15,
    15,
    14,
    11,
    10,
    9,
    7,
    6,
    6,
    6,
    5,
    5,
    4,
    3
  )
  expect_equal(sum(y), 1471)
  approx_loo <- rep(0L, 676)
  expect_equal(sum(y^2), 54497)
  res <- loo:::srs_est(y = y, approx_loo)
  expect_equal(res$y_hat, 19888, tolerance = 0.0001)
  expect_equal(res$v_y_hat, 676^2 * 229 * (1 - 0.074) / 50, tolerance = 0.0001)
  expect_equal(res$hat_v_y, 676 * var(y), tolerance = 0.0001)

  # Simulation example
  set.seed(1234)
  N <- 1000
  y_true <- 1:N
  sigma_hat_true <- sqrt(N * sum((y_true - mean(y_true))^2) / length(y_true))

  m <- 100
  y_hat <- se_y_hat <- sigma_hat <- numeric(10000)
  for (i in 1:10000) {
    y_idx <- sample(1:N, size = m)
    y <- y_true[y_idx]
    res <- loo:::srs_est(y = y, y_approx = y_true)
    y_hat[i] <- res$y_hat
    se_y_hat[i] <- sqrt(res$v_y_hat)
    sigma_hat[i] <- sqrt(res$hat_v_y)
  }
  expect_equal(mean(y_hat), sum(y_true), tolerance = 0.1)

  in_ki <- y_hat + 2 * se_y_hat > sum(y_true) &
    y_hat - 2 * se_y_hat < sum(y_true)
  expect_equal(mean(in_ki), 0.95, tolerance = 0.01)

  # Should be  unbiased
  expect_equal(mean(sigma_hat), sigma_hat_true, tolerance = 0.1)

  m <- N
  y_idx <- sample(1:N, size = m)
  y <- y_true[y_idx]
  res <- loo:::srs_est(y, y_true)
  expect_equal(res$y_hat, 500500, tolerance = 0.0001)
  expect_equal(res$v_y_hat, 0, tolerance = 0.0001)
})
