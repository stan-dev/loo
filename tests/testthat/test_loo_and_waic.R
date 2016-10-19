library(loo)
options(loo.cores = 1)
set.seed(123)
x <- matrix(rnorm(5000), 100, 50)

context("loo and waic")
test_that("loo ok with > 1 core", {
  expect_warning(loo(x, cores = 2), "Some Pareto k diagnostic values are too high")
})

test_that("loo and waic return expected results", {
  expect_warning(ww <- waic(x), "p_waic")
  wnms <- names(ww)

  expect_warning(ll <- loo(x), "Some Pareto k diagnostic values are too high")
  lnms <- names(ll)

  waic_val <- unlist(ww[grep("pointwise", wnms, invert = TRUE)])
  loo_val <- unlist(ll[grep("pointwise|pareto_k", lnms, invert = TRUE)])

  waic_ans <- structure(c(-25.1291601624985, 49.5591373599165, 50.258320324997,
                          0.742909828603039, 1.05244649202819, 1.48581965720608),
                        .Names = c("elpd_waic", "p_waic", "waic", "se_elpd_waic",
                                   "se_p_waic", "se_waic"))
  loo_ans <- structure(c(-24.2339828660691, 48.6639600634871, 48.4679657321382,
                         0.706772789390649, 0.990389359906155, 1.4135455787813),
                       .Names = c("elpd_loo", "p_loo", "looic", "se_elpd_loo",
                                  "se_p_loo", "se_looic"))
  expect_equal(loo_val, loo_ans)
  expect_equal(waic_val, waic_ans)

  pareto_k <- ll$pareto_k
  pareto_k_val <- list(mean = mean(pareto_k), range = range(pareto_k))
  pareto_k_ans <- list(mean = 0.273775211941035,
                       range = c(-0.245572962408798, 0.859988648019406))
  expect_equal(pareto_k_val, pareto_k_ans)
})

test_that("loo and waic throw appropriate errors", {
  x[1,1] <- NaN
  expect_error(loo(x), regexp = "NA log-likelihood")
  expect_error(waic(x), regexp = "NA log-likelihood")

  vec <- 1:10
  arr <- array(1:100, dim = c(2,5,10))
  expect_error(loo(vec), regexp = "no applicable method")
  expect_error(waic(arr), regexp = "no applicable method")
})

test_that("function and matrix methods return same result", {
  set.seed(024)

  # fake data and posterior draws
  N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
  p <- rbeta(1, a0, b0)
  y <- rbinom(N, size = K, prob = p)
  a <- a0 + sum(y); b <- b0 + N * K - sum(y)
  draws <- rbeta(S, a, b)
  data <- data.frame(y,K)
  llfun <- function(i, data, draws) {
    dbinom(data$y, size = data$K, prob = draws, log = TRUE)
  }
  loo_with_fn <- loo(llfun, args = nlist(data, draws, N, S))
  waic_with_fn <- waic(llfun, args = nlist(data, draws, N, S))

  # Check that we get same answer if using log-likelihood matrix
  log_lik_mat <- sapply(1:N, function(i) llfun(i, data[i,, drop=FALSE], draws))
  loo_with_mat <- loo(log_lik_mat)
  waic_with_mat <- waic(log_lik_mat)
  expect_equal(loo_with_mat, loo_with_fn)
  expect_equal(waic_with_mat, waic_with_fn)
})


# pareto_k_ids ------------------------------------------------------------
test_that("pareto_k_ids identifies correct observations", {
  ll <- suppressWarnings(loo(x))
  expect_identical(pareto_k_ids(ll, threshold = 0.5),
                   which(ll$pareto_k > 0.5))
  expect_identical(pareto_k_ids(ll, threshold = 1),
                   which(ll$pareto_k > 1))
})
