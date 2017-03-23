library(loo)

context("E_loo")

test_that("E_loo return types are correct", {
  set.seed(123)
  x <- matrix(rnorm(200), 20, 10)
  lw <- matrix(rnorm(200), 20, 10)

  expect_null(dim(E_loo(x[,1], lw[,1], type="mean")))
  expect_equal(length(E_loo(x[,1], lw[,1], type="var")), 1)
  expect_equal(
    length(E_loo(x[,1], lw[,1], type="quantile", probs = c(.25, .75))),
    2
  )

  expect_null(dim(E_loo(x, lw, type="mean")))
  expect_null(dim(E_loo(x, lw, type="var")))
  expect_equal(length(E_loo(x, lw, type="mean")), 10)
  expect_equal(length(E_loo(x, lw, type="var")), 10)
  expect_equal(
    dim(E_loo(x, lw, type="quantile", probs = c(.25, .75))),
    c(2, 10)
  )
})

test_that("E_loo.default equal to reference", {
  set.seed(124)
  x <- rnorm(100)
  lw <- rnorm(100)

  expect_equal_to_reference(
    E_loo(x, lw, type = "mean"),
    "E_loo_default_mean.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "var"),
    "E_loo_default_var.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "quantile", probs = 0.5),
    "E_loo_default_quantile_50.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "quantile", probs = c(0.1, 0.9)),
    "E_loo_default_quantile_10_90.rds"
  )
})

test_that("E_loo.matrix equal to reference", {
  set.seed(125)
  x <- matrix(rnorm(200), 20, 10)
  lw <- matrix(rnorm(200), 20, 10)

  expect_equal_to_reference(
    E_loo(x, lw, type = "mean"),
    "E_loo_matrix_mean.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "var"),
    "E_loo_matrix_var.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "quantile", probs = 0.5),
    "E_loo_matrix_quantile_50.rds"
  )
  expect_equal_to_reference(
    E_loo(x, lw, type = "quantile", probs = c(0.1, 0.9)),
    "E_loo_matrix_quantile_10_90.rds"
  )
})

test_that("E_loo throws correct errors", {
  expect_error(E_loo(1, 1, type = "purple"))
  expect_error(
    E_loo(1, 1, type = "quantile", probs = 2),
    "all(probs > 0 & probs < 1) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo("a", "b"),
    "is.numeric(x) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo(1:10, 1:5),
    "length(lw) == length(x) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo(cbind(1:10, 1:10), cbind(1:5, 1:5)),
    "identical(dim(x), dim(lw)) is not TRUE",
    fixed = TRUE
  )
})




.wquant_rapprox <- function(x, w, probs) {
  stopifnot(all(probs > 0 & probs < 1))
  ord <- order(x)
  d <- x[ord]
  ww <- w[ord]
  p <- cumsum(ww) / sum(ww)
  stats::approx(p, d, probs, rule = 2)$y
}
.wquant_sim <- function(x, w, probs, n_sims) {
  xx <- sample(x, size = n_sims, replace = TRUE, prob = w / sum(w))
  quantile(xx, probs, names = FALSE)
}


test_that("weighted quantiles work", {
  set.seed(123)
  pr <- seq(0.025, 0.975, 0.025)

  x1 <- rnorm(100)
  w1 <- rlnorm(100)
  expect_equal(
    .wquant(x1, w1, pr),
    .wquant_rapprox(x1, w1, pr)
  )

  x1 <- rnorm(1e4)
  w1 <- rlnorm(1e4)
  expect_equal(
    .wquant(x1, w1, pr),
    .wquant_sim(x1, w1, pr, n_sim = 2e6),
    tol = 0.001
  )

  expect_equal(
    .wquant(x1, rep(1, length(x1)), pr),
    quantile(x1, probs = pr, names = FALSE)
  )
})

