# Create test data
# Checked by Mans M and Paul B 24th of June 2019
set.seed(123)
N <- 50
K <- 10
S <- 1000
a0 <- 1
b0 <- 1
p <- 0.5
y <- rbinom(N, size = K, prob = p)
fake_data <- data.frame(y, K)

# The log posterior
log_post <- function(p, y, a0, b0, K) {
  log_lik <- sum(dbinom(x = y, size = K, prob = p, log = TRUE)) # the log likelihood
  log_post <- log_lik + dbeta(x = p, shape1 = a0, shape2 = b0, log = TRUE) # the log prior
  log_post
}
it <- optim(
  par = 0.5,
  fn = log_post,
  control = list(fnscale = -1),
  hessian = TRUE,
  y = y,
  a0 = a0,
  b0 = b0,
  K = K,
  lower = 0.01,
  upper = 0.99,
  method = "Brent"
)
lap_params <- c(mu = it$par, sd = sqrt(solve(-it$hessian)))

a <- a0 + sum(y)
b <- b0 + N * K - sum(y)
fake_true_posterior <- as.matrix(rbeta(S, a, b))
fake_laplace_posterior <- as.matrix(rnorm(
  n = S,
  mean = lap_params["mu"],
  sd = lap_params["sd"]
))
# mean(fake_laplace_posterior); sd(fake_laplace_posterior)

p_draws <- as.vector(fake_laplace_posterior)
log_p <- numeric(S)
for (s in 1:S) {
  log_p[s] <- log_post(p_draws[s], y = y, a0 = a0, b0 = b0, K = K)
}
log_g <- as.vector(dnorm(
  as.vector(fake_laplace_posterior),
  mean = lap_params["mu"],
  sd = lap_params["sd"],
  log = TRUE
))

llfun <- function(data_i, draws) {
  dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
}

ll <- matrix(0, nrow = S, ncol = N)
for (j in 1:N) {
  ll[, j] <- llfun(
    data_i = fake_data[j, , drop = FALSE],
    draws = fake_laplace_posterior
  )
}


test_that("loo_approximate_posterior.array works as loo_approximate_posterior.matrix", {
  # Create array with two "chains"
  log_p_mat <- matrix(log_p, nrow = (S / 2), ncol = 2)
  log_g_mat <- matrix(log_g, nrow = (S / 2), ncol = 2)
  ll_array <- array(0, dim = c((S / 2), 2, ncol(ll)))
  ll_array[, 1, ] <- ll[1:(S / 2), ]
  ll_array[, 2, ] <- ll[(S / 2 + 1):S, ]

  # Assert that they are ok
  expect_equal(ll_array[1:2, 1, 1:2], ll[1:2, 1:2], ignore_attr = TRUE)
  expect_equal(
    ll_array[1:2, 2, 1:2],
    ll[(S / 2 + 1):((S / 2) + 2), 1:2],
    ignore_attr = TRUE
  )

  # Compute aploo
  expect_silent(
    aploo1 <- loo_approximate_posterior.matrix(
      x = ll,
      log_p = log_p,
      log_g = log_g
    )
  )
  expect_silent(
    aploo2 <- loo_approximate_posterior.array(
      x = ll_array,
      log_p = log_p_mat,
      log_g = log_g_mat
    )
  )
  expect_silent(aploo1b <- loo.matrix(x = ll, r_eff = rep(1, N)))

  # Check equivalence
  expect_equal(aploo1$estimates, aploo2$estimates)
  expect_equal(class(aploo1), class(aploo2))
  expect_failure(expect_equal(aploo1b$estimates, aploo2$estimates))
  expect_failure(expect_equal(class(aploo1), class(aploo1b)))

  # Should fail with matrix
  expect_error(
    aploo2 <- loo_approximate_posterior.matrix(
      x = ll,
      log_p = as.matrix(log_p),
      log_g = log_g
    )
  )
  expect_error(
    aploo2 <- loo_approximate_posterior.matrix(
      x = ll,
      log_p = as.matrix(log_p),
      log_g = as.matrix(log_g)
    )
  )

  # Expect log_p and log_g be stored in the approximate_posterior in the same way
  expect_length(aploo1$approximate_posterior$log_p, nrow(ll))
  expect_length(aploo1$approximate_posterior$log_g, nrow(ll))
  expect_equal(
    aploo1$approximate_posterior$log_p,
    aploo2$approximate_posterior$log_p
  )
  expect_equal(
    aploo1$approximate_posterior$log_g,
    aploo2$approximate_posterior$log_g
  )
})

test_that("loo_approximate_posterior.function works as loo_approximate_posterior.matrix", {
  # Compute aploo
  expect_silent(
    aploo1 <- loo_approximate_posterior.matrix(
      x = ll,
      log_p = log_p,
      log_g = log_g
    )
  )
  expect_silent(aploo1b <- loo.matrix(x = ll, r_eff = rep(1, N)))
  expect_silent(
    aploo2 <- loo_approximate_posterior.function(
      x = llfun,
      log_p = log_p,
      log_g = log_g,
      data = fake_data,
      draws = fake_laplace_posterior
    )
  )

  # Check equivalence
  expect_equal(aploo1$estimates, aploo2$estimates)
  expect_equal(class(aploo1), class(aploo2))
  expect_failure(expect_equal(aploo1b$estimates, aploo2$estimates))

  # Check equivalence
  # Expect log_p and log_g be stored in the approximate_posterior in the same way
  expect_length(
    aploo2$approximate_posterior$log_p,
    nrow(fake_laplace_posterior)
  )
  expect_length(
    aploo2$approximate_posterior$log_g,
    nrow(fake_laplace_posterior)
  )
  expect_equal(
    aploo1$approximate_posterior$log_p,
    aploo2$approximate_posterior$log_p
  )
  expect_equal(
    aploo1$approximate_posterior$log_g,
    aploo2$approximate_posterior$log_g
  )
})
