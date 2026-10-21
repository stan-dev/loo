# generate fake data
set.seed(123)
y <- rnorm(50, 0, 1)
sd_sim1 <- abs(rnorm(500, 1.5, 0.1))
sd_sim2 <- abs(rnorm(500, 1.2, 0.1))
sd_sim3 <- abs(rnorm(500, 1, 0.05))
log_lik1 <- log_lik2 <- log_lik3 <- matrix(NA, 500, 50)
for (s in 1:500) {
  log_lik1[s, ] <- dnorm(y, -1, sd_sim1[s], log = T)
  log_lik2[s, ] <- dnorm(y, 0.7, sd_sim2[s], log = T)
  log_lik3[s, ] <- dnorm(y, 1, sd_sim3[s], log = T)
}

ll_list <- list(log_lik1, log_lik2, log_lik3)
r_eff_list <- list(rep(0.9, 50), rep(0.9, 50), rep(0.9, 50))

loo_list <- lapply(1:length(ll_list), function(j) {
  loo(ll_list[[j]], r_eff = r_eff_list[[j]])
})

tol <- 0.01 # absolute tolerance of weights

test_that("stacking handles equal infinite log predictive densities", {
  lpd <- matrix(
    c(-Inf, 0, -Inf, -1, -1, -1),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(as.numeric(stacking_weights(lpd)), c(0, 1, 0), tolerance = 1e-6)
})

test_that("loo_model_weights throws correct errors and warnings", {
  expect_error(
    loo_model_weights(log_lik1),
    "list of matrices or a list of 'psis_loo' objects"
  )
  expect_error(loo_model_weights(list(log_lik1)), "At least two models")
  expect_error(loo_model_weights(list(loo_list[[1]])), "At least two models")
  expect_error(
    loo_model_weights(list(log_lik1), method = "pseudobma"),
    "At least two models"
  )

  expect_error(
    loo_model_weights(list(log_lik1, log_lik2[-1, ])),
    "same dimensions"
  )
  expect_error(
    loo_model_weights(list(log_lik1, log_lik2, log_lik3[, -1])),
    "same dimensions"
  )

  loo_list2 <- loo_list
  attr(loo_list2[[3]], "dims") <- c(10, 10)
  expect_error(loo_model_weights(loo_list2), "same dimensions")

  expect_error(
    loo_model_weights(ll_list, r_eff_list = r_eff_list[-1]),
    "one component for each model"
  )

  r_eff_list[[3]] <- rep(0.9, 51)
  expect_error(
    loo_model_weights(ll_list, r_eff_list = r_eff_list),
    "same length as the number of columns"
  )

  expect_error(
    loo_model_weights(list(loo_list[[1]], 2)),
    "List elements must all be 'psis_loo' objects or log-likelihood matrices",
    fixed = TRUE
  )

  expect_no_warning(loo_model_weights(ll_list))
})


test_that("loo_model_weights (stacking and pseudo-BMA) gives expected result", {
  w1 <- loo_model_weights(ll_list, method = "stacking", r_eff_list = r_eff_list)
  expect_type(w1, "double")
  expect_s3_class(w1, "stacking_weights")
  expect_length(w1, 3)
  expect_named(w1, paste0("model", c(1:3)))
  expect_snapshot_value(as.numeric(w1), style = "serialize")
  expect_snapshot(print(w1))

  w1_b <- loo_model_weights(loo_list)
  expect_identical(w1, w1_b)

  w2 <- loo_model_weights(
    ll_list,
    r_eff_list = r_eff_list,
    method = "pseudobma",
    BB = TRUE
  )
  expect_type(w2, "double")
  expect_s3_class(w2, "pseudobma_bb_weights")
  expect_length(w2, 3)
  expect_named(w2, paste0("model", c(1:3)))
  expect_snapshot_value(as.numeric(w2), style = "serialize")
  expect_snapshot(print(w2))

  w3 <- loo_model_weights(
    ll_list,
    r_eff_list = r_eff_list,
    method = "pseudobma",
    BB = FALSE
  )
  expect_type(w3, "double")
  expect_length(w3, 3)
  expect_named(w3, paste0("model", c(1:3)))
  expect_equal(
    as.numeric(w3),
    c(5.365279e-05, 9.999436e-01, 2.707028e-06),
    tolerance = tol
  )
  expect_snapshot(print(w3))

  w3_b <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
  expect_identical(w3, w3_b)
})

# Upper bound of how far the log score of the weights `w` is below its
# maximum: N * log(max_k g_k), with g_k = mean_i p_ik / sum_j w_j p_ij
stacking_gap_bound <- function(w, lpd) {
  log_p <- matrixStats::rowLogSumExps(sweep(lpd, 2, log(w), "+"))
  nrow(lpd) * max(matrixStats::colLogSumExps(lpd - log_p) - log(nrow(lpd)))
}

# EM for mixture weights with fixed components: a slow independent solver
stacking_em <- function(lpd, iter = 20000) {
  w <- rep(1 / ncol(lpd), ncol(lpd))
  for (i in seq_len(iter)) {
    a <- sweep(lpd, 2, log(w), "+")
    w <- colMeans(exp(a - matrixStats::rowLogSumExps(a)))
  }
  w
}

test_that("stacking_weights gives expected result", {
  set.seed(10)
  lpd_point <- matrix(rnorm(2000, -1), 400, 5)

  actual <- stacking_weights(lpd_point)

  expect_s3_class(actual, "stacking_weights")
  expect_named(actual, paste0("model", 1:5))
  expect_equal(sum(actual), 1)
  expect_lt(stacking_gap_bound(as.numeric(actual), lpd_point), 400 * 1e-8)
  expect_equal(as.numeric(actual), stacking_em(lpd_point), tolerance = 1e-5)
})

test_that("stacking_weights finds the optimum for many models", {
  set.seed(119)
  N <- 100
  K <- 500
  y <- rnorm(N)
  mu <- rnorm(K, 0, 0.5)
  s <- exp(rnorm(K, 0, 0.3))
  lpd <- sapply(seq_len(K), function(k) dnorm(y, mu[k], s[k], log = TRUE))

  w <- as.numeric(stacking_weights(lpd))
  expect_equal(sum(w), 1)
  expect_true(all(w >= 0))
  expect_lt(stacking_gap_bound(w, lpd), N * 1e-8)
  # the optimal weights are sparse
  expect_lt(sum(w > 0), 20)
})

test_that("stacking_weights handles -Inf, duplicate and identical models", {
  set.seed(2)
  N <- 60
  lpd <- matrix(rnorm(3 * N, -1), N, 3)
  # each model is impossible for a third of the observations, so the
  # combination needs all three
  lpd[1:20, 1] <- -Inf
  lpd[21:40, 2] <- -Inf
  lpd[41:60, 3] <- -Inf
  w <- as.numeric(stacking_weights(lpd))
  expect_true(all(w > 0))
  expect_lt(stacking_gap_bound(w, lpd), N * 1e-8)

  # duplicate models share the weight that one of them would get
  lpd2 <- matrix(rnorm(4 * N, -1), N, 4)
  w2 <- as.numeric(stacking_weights(lpd2))
  w_dup <- as.numeric(stacking_weights(cbind(lpd2, lpd2[, 1])))
  expect_equal(c(w_dup[1] + w_dup[5], w_dup[2:4]), w2, tolerance = 1e-6)

  # identical models: all weight to one of them
  lpd3 <- matrix(rep(rnorm(N, -1), 3), N, 3)
  expect_equal(sort(as.numeric(stacking_weights(lpd3))), c(0, 0, 1))
})

test_that("optim_method and optim_control are deprecated", {
  set.seed(3)
  lpd <- matrix(rnorm(40, -1), 20, 2)
  w <- stacking_weights(lpd)
  expect_warning(
    w2 <- stacking_weights(lpd, optim_control = list(reltol = 1e-10)),
    "deprecated and ignored"
  )
  expect_identical(w2, w)
  expect_warning(
    loo_model_weights(ll_list, optim_method = "Nelder-Mead"),
    "deprecated and ignored"
  )
  expect_no_warning(loo_model_weights(ll_list))
})

test_that("Bayesian bootstrap gives expected result", {
  lpd_point <- rbind(
    c(-0.2, -0.8, -1.1),
    c(-1.4, -0.3, -0.5),
    c(-0.6, -0.7, -0.1),
    c(-1.0, -1.2, -0.4),
    c(-0.9, -0.5, -0.8)
  )
  BB_n <- 25
  alpha <- 0.7

  set.seed(0)
  weights <- pseudobma_weights(lpd_point, BB = TRUE, BB_n = BB_n, alpha = alpha)

  expect_s3_class(weights, "pseudobma_bb_weights")
  expect_named(weights, paste0("model", 1:3))
  expect_equal(sum(weights), 1)
  expect_snapshot_value(
    unname(as.numeric(weights)),
    style = "deparse"
  )

  # Compare to a reference loop over the same Bayesian bootstrap draws
  set.seed(0)
  BB_weighting <- dirichlet_rng(BB_n, rep(alpha, nrow(lpd_point)))
  expected <- matrix(NA, BB_n, ncol(lpd_point))
  for (bb in 1:BB_n) {
    z_bb <- BB_weighting[bb, ] %*% lpd_point * nrow(lpd_point)
    uwts <- exp(z_bb - max(z_bb))
    expected[bb, ] <- uwts / sum(uwts)
  }
  expect_equal(as.numeric(weights), colMeans(expected))
})

test_that("pseudo-BMA gives zero weight to an impossible model", {
  lpd <- cbind(rep(-Inf, 3), c(-2, -1, 0))

  expect_equal(as.numeric(pseudobma_weights(lpd, BB = FALSE)), c(0, 1))
})

lpd_bad_values <- list(`NA` = NA_real_, `NaN` = NaN, `Inf` = Inf)

# One test_that() per value, so a failure names which one broke
for (lpd_bad_name in names(lpd_bad_values)) {
  test_that(paste("model weighting rejects", lpd_bad_name, "predictive densities"), {
    bad <- lpd_bad_values[[lpd_bad_name]]
    set.seed(1)
    lpd <- cbind(rnorm(20, -2), rnorm(20, -2.2))
    error <- "All values in 'lpd_point' must be finite or -Inf."

    single <- lpd
    single[3, 2] <- bad
    expect_error(stacking_weights(single), error, fixed = TRUE)
    expect_error(pseudobma_weights(single, BB = FALSE), error, fixed = TRUE)
    expect_error(pseudobma_weights(single, BB = TRUE), error, fixed = TRUE)

    whole_column <- lpd
    whole_column[, 2] <- bad
    expect_error(stacking_weights(whole_column), error, fixed = TRUE)
    expect_error(pseudobma_weights(whole_column, BB = FALSE), error, fixed = TRUE)
  })
}

test_that("model weighting still accepts -Inf predictive densities", {
  set.seed(1)
  lpd <- cbind(rnorm(20, -2), rnorm(20, -2.2))
  lpd[3, 2] <- -Inf

  expect_no_error(stacking_weights(lpd))
  expect_no_error(pseudobma_weights(lpd, BB = FALSE))
  expect_no_error(pseudobma_weights(lpd, BB = TRUE))
})

test_that("model weighting rejects inputs with no finite predictive density", {
  stacking_lpd <- rbind(c(-Inf, -Inf), c(-1, -1))
  expect_error(
    stacking_weights(stacking_lpd),
    "Each observation must have a finite predictive density for at least one model.",
    fixed = TRUE
  )

  pseudobma_lpd <- matrix(-Inf, nrow = 2, ncol = 2)
  error <- "At least one model must have a finite total predictive density."
  expect_error(pseudobma_weights(pseudobma_lpd, BB = FALSE), error, fixed = TRUE)
  expect_error(pseudobma_weights(pseudobma_lpd, BB = TRUE), error, fixed = TRUE)
})

test_that("stacking_weights and pseudobma_weights throw correct errors", {
  xx <- cbind(rnorm(10))
  expect_error(stacking_weights(xx), "two models are required")
  expect_error(pseudobma_weights(xx), "two models are required")
})

test_that("loo_model_weights uses correct names for list of loo objects", {
  loo1 <- loo_list[[1]]
  loo2 <- loo_list[[2]]
  loo3 <- loo_list[[3]]

  expect_named(
    loo_model_weights(list(loo1, loo2, loo3)),
    c("model1", "model2", "model3")
  )
  expect_named(
    loo_model_weights(list("a" = loo1, loo2, "c" = loo3)),
    c("a", "model2", "c")
  )
  expect_named(
    loo_model_weights(list(`a` = loo1, `b` = loo2, `c` = loo3)),
    c("a", "b", "c")
  )
})
