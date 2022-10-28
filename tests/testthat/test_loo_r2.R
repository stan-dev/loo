options(mc.cores = 1)
set.seed(123)
context("r2")

LL <- example_loglik_matrix()
chain_id <- rep(1:2, each = dim(LL)[1] / 2)
r_eff <- relative_eff(exp(LL), chain_id)
psisd_obj <- psis(-LL, r_eff = r_eff, cores = 2)

yrep <- matrix(rnorm(length(LL)), nrow = nrow(LL), ncol = ncol(LL))
y <- rnorm(ncol(LL))


r2 <- loo_r2(y = y,
             yrep = yrep ,
             log_lik = LL,
             r_eff = r_eff)

test_that("loo_r2 stops with incorrect inputs", {
  expect_error(loo_r2(as.character(y), yrep, LL, r_eff = r_eff),
               "is.numeric(y) is not TRUE",
               fixed = TRUE)

  expect_error(loo_r2(y, as.character(yrep), LL, r_eff = r_eff),
               "is.numeric(yrep) is not TRUE",
               fixed = TRUE)

  yrep_invalid <- matrix(rnorm(9), nrow = 3)
  expect_error(loo_r2(y, yrep_invalid, LL, r_eff = r_eff),
               "identical(ncol(yrep), length(y)) is not TRUE",
               fixed = TRUE)

  yrep_invalid <- matrix(rnorm(64), nrow = 2)
  expect_error(loo_r2(y, yrep_invalid, LL, r_eff = r_eff),
               "identical(dim(yrep), dim(log_lik)) is not TRUE",
               fixed = TRUE)
})


test_that("loo_r2 return types are correct", {
  expect_type(r2, "list")
  expect_named(r2, c("estimate", "pointwise", "se", "diagnostics"))
})

test_that("loo_r2 results haven't changed", {
  expect_equal_to_reference(r2, "reference-results/loo_r2.rds")
})
