library(loo)

context("E_loo")

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_mat <- relative_eff(exp(LLmat), chain_id)
r_eff_vec <- relative_eff(exp(LLvec), chain_id = chain_id)
psis_mat <- psis(-LLmat, r_eff = r_eff_mat, cores = 2)
psis_vec <- psis(-LLvec, r_eff = r_eff_vec)

set.seed(123)
x <- matrix(rnorm(length(LLmat)), nrow = nrow(LLmat), ncol = ncol(LLmat))
log_rats <- -LLmat

# matrix method
E_test_mean <- E_loo(x, psis_mat, type = "mean", log_ratios = log_rats)
E_test_var <- E_loo(x, psis_mat, type = "var", log_ratios = log_rats)
E_test_sd <- E_loo(x, psis_mat, type = "sd", log_ratios = log_rats)
E_test_quant <- E_loo(x, psis_mat, type = "quantile", probs = 0.5, log_ratios = log_rats)
E_test_quant2 <- E_loo(x, psis_mat, type = "quantile", probs = c(0.1, 0.9), log_ratios = log_rats)

# vector method
E_test_mean_vec <- E_loo(x[, 1], psis_vec, type = "mean", log_ratios = log_rats[,1])
E_test_var_vec <- E_loo(x[, 1], psis_vec, type = "var", log_ratios = log_rats[,1])
E_test_sd_vec <- E_loo(x[, 1], psis_vec, type = "sd", log_ratios = log_rats[,1])
E_test_quant_vec <- E_loo(x[, 1], psis_vec, type = "quant", probs = 0.5, log_ratios = log_rats[,1])
E_test_quant_vec2 <- E_loo(x[, 1], psis_vec, type = "quant", probs = c(0.1, 0.5, 0.9), log_ratios = log_rats[,1])

# E_loo_khat
khat <- loo:::E_loo_khat.matrix(x, psis_mat, log_rats)

test_that("E_loo return types correct for matrix method", {
  expect_type(E_test_mean, "list")
  expect_named(E_test_mean, c("value", "pareto_k"))
  expect_length(E_test_mean, 2)
  expect_length(E_test_mean$value, ncol(x))
  expect_length(E_test_mean$pareto_k, ncol(x))

  expect_type(E_test_var, "list")
  expect_named(E_test_var, c("value", "pareto_k"))
  expect_length(E_test_var, 2)
  expect_length(E_test_var$value, ncol(x))
  expect_length(E_test_var$pareto_k, ncol(x))

  expect_type(E_test_sd, "list")
  expect_named(E_test_sd, c("value", "pareto_k"))
  expect_length(E_test_sd, 2)
  expect_length(E_test_sd$value, ncol(x))
  expect_length(E_test_sd$pareto_k, ncol(x))

  expect_type(E_test_quant, "list")
  expect_named(E_test_quant, c("value", "pareto_k"))
  expect_length(E_test_quant, 2)
  expect_length(E_test_quant$value, ncol(x))
  expect_length(E_test_quant$pareto_k, ncol(x))

  expect_type(E_test_quant2, "list")
  expect_named(E_test_quant2, c("value", "pareto_k"))
  expect_length(E_test_quant2, 2)
  expect_equal(dim(E_test_quant2$value), c(2, ncol(x)))
  expect_length(E_test_quant2$pareto_k, ncol(x))
})

test_that("E_loo return types correct for default/vector method", {
  expect_type(E_test_mean_vec, "list")
  expect_named(E_test_mean_vec, c("value", "pareto_k"))
  expect_length(E_test_mean_vec, 2)
  expect_length(E_test_mean_vec$value, 1)
  expect_length(E_test_mean_vec$pareto_k, 1)

  expect_type(E_test_var_vec, "list")
  expect_named(E_test_var_vec, c("value", "pareto_k"))
  expect_length(E_test_var_vec, 2)
  expect_length(E_test_var_vec$value, 1)
  expect_length(E_test_var_vec$pareto_k, 1)

  expect_type(E_test_sd_vec, "list")
  expect_named(E_test_sd_vec, c("value", "pareto_k"))
  expect_length(E_test_sd_vec, 2)
  expect_length(E_test_sd_vec$value, 1)
  expect_length(E_test_sd_vec$pareto_k, 1)

  expect_type(E_test_quant_vec, "list")
  expect_named(E_test_quant_vec, c("value", "pareto_k"))
  expect_length(E_test_quant_vec, 2)
  expect_length(E_test_quant_vec$value, 1)
  expect_length(E_test_quant_vec$pareto_k, 1)

  expect_type(E_test_quant_vec2, "list")
  expect_named(E_test_quant_vec2, c("value", "pareto_k"))
  expect_length(E_test_quant_vec2, 2)
  expect_length(E_test_quant_vec2$value, 3)
  expect_length(E_test_quant_vec2$pareto_k, 1)
})

test_that("E_loo.default equal to reference", {
  expect_equal_to_reference(E_test_mean_vec, test_path("reference-results/E_loo_default_mean.rds"))
  expect_equal_to_reference(E_test_var_vec, test_path("reference-results/E_loo_default_var.rds"))
  expect_equal_to_reference(E_test_sd_vec, test_path("reference-results/E_loo_default_sd.rds"))
  expect_equal_to_reference(E_test_quant_vec, test_path("reference-results/E_loo_default_quantile_50.rds"))
  expect_equal_to_reference(E_test_quant_vec2, test_path("reference-results/E_loo_default_quantile_10_50_90.rds"))
})

test_that("E_loo.matrix equal to reference", {
  expect_equal_to_reference(E_test_mean, test_path("reference-results/E_loo_matrix_mean.rds"))
  expect_equal_to_reference(E_test_var, test_path("reference-results/E_loo_matrix_var.rds"))
  expect_equal_to_reference(E_test_sd, test_path("reference-results/E_loo_matrix_sd.rds"))
  expect_equal_to_reference(E_test_quant, test_path("reference-results/E_loo_matrix_quantile_50.rds"))
  expect_equal_to_reference(E_test_quant2, test_path("reference-results/E_loo_matrix_quantile_10_90.rds"))
})

test_that("E_loo throws correct errors and warnings", {
  # warnings
  expect_no_warning(E_loo.matrix(x, psis_mat))
  # no warnings if x is constant, binary, NA, NaN, Inf
  expect_no_warning(E_loo.matrix(x*0, psis_mat))
  expect_no_warning(E_loo.matrix(0+(x>0), psis_mat))
  expect_no_warning(E_loo.matrix(x+NA, psis_mat))   
  expect_no_warning(E_loo.matrix(x*NaN, psis_mat))
  expect_no_warning(E_loo.matrix(x*Inf, psis_mat))
  expect_no_warning(E_test <- E_loo.default(x[, 1], psis_vec))
  expect_length(E_test$pareto_k, 1)

  # errors
  expect_error(E_loo(x, 1), "is.psis")
  expect_error(
    E_loo(x, psis_mat, type = "quantile", probs = 2),
    "all(probs > 0 & probs < 1) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo(rep("a", nrow(x)), psis_vec),
    "is.numeric(x) is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo(1:10, psis_vec),
    "length(x) == dim(psis_object)[1] is not TRUE",
    fixed = TRUE
  )
  expect_error(
    E_loo(cbind(1:10, 1:10), psis_mat),
    "identical(dim(x), dim(psis_object)) is not TRUE",
    fixed = TRUE
  )
})


test_that("weighted quantiles work", {
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
  # expect_equal(
  #   .wquant(x1, w1, pr),
  #   .wquant_sim(x1, w1, pr, n_sim = 5e6),
  #   tol = 0.005
  # )

  expect_equal(
    .wquant(x1, rep(1, length(x1)), pr),
    quantile(x1, probs = pr, names = FALSE)
  )
})

test_that("weighted variance works", {
  x <- rnorm(100)
  w <- rep(0.01, 100)
  expect_equal(.wvar(x, w), var(x))
  expect_equal(.wsd(x, w), sqrt(.wvar(x, w)))

  w <- c(rep(0.1, 10), rep(0, 90))
  expect_equal(.wvar(x, w), var(x[w > 0]))
})
