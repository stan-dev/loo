options(mc.cores = 1)
set.seed(123)

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
LLvec <- LLmat[, 1]
chain_id <- rep(1:2, each = nrow(LLarr))
r_eff_arr <- relative_eff(exp(LLarr))
r_eff_mat <- relative_eff(exp(LLmat), chain_id = chain_id)

loo1 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))
waic1 <- suppressWarnings(waic(LLarr))
elpd1 <- suppressWarnings(elpd(LLarr))

test_that("using loo.cores is deprecated", {
  options(mc.cores = NULL)
  options(loo.cores = 1)
  expect_warning(loo(LLarr, r_eff = r_eff_arr, cores = 2), "loo.cores")
  options(loo.cores = NULL)
  options(mc.cores = 1)
})

test_that("loo, waic and elpd results haven't changed", {
  expect_snapshot_value(loo1, style = "serialize")
  expect_snapshot_value(waic1, style = "serialize")
  expect_snapshot_value(elpd1, style = "serialize")
})

test_that("loo with cores=1 and cores=2 gives same results", {
  loo2 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr, cores = 2))
  expect_equal(loo1$estimates, loo2$estimates)
})

test_that("waic returns object with correct structure", {
  expect_true(is.waic(waic1))
  expect_true(is.loo(waic1))
  expect_false(is.psis_loo(waic1))
  expect_named(
    waic1,
    c(
      "estimates",
      "pointwise",

      # deprecated but still there
      "elpd_waic",
      "p_waic",
      "waic",
      "se_elpd_waic",
      "se_p_waic",
      "se_waic"
    )
  )
  est_names <- dimnames(waic1$estimates)
  expect_equal(est_names[[1]], c("elpd_waic", "p_waic", "waic"))
  expect_equal(est_names[[2]], c("Estimate", "SE"))
  expect_equal(colnames(waic1$pointwise), est_names[[1]])
  expect_equal(dim(waic1), dim(LLmat))
})

test_that("loo returns object with correct structure", {
  expect_false(is.waic(loo1))
  expect_true(is.loo(loo1))
  expect_true(is.psis_loo(loo1))
  expect_named(
    loo1,
    c(
      "estimates",
      "pointwise",
      "diagnostics",
      "psis_object",

      # deprecated but still there
      "elpd_loo",
      "p_loo",
      "looic",
      "se_elpd_loo",
      "se_p_loo",
      "se_looic"
    )
  )
  expect_named(loo1$diagnostics, c("pareto_k", "n_eff", "r_eff"))
  expect_equal(dimnames(loo1$estimates)[[1]], c("elpd_loo", "p_loo", "looic"))
  expect_equal(dimnames(loo1$estimates)[[2]], c("Estimate", "SE"))
  expect_equal(
    colnames(loo1$pointwise),
    c("elpd_loo", "mcse_elpd_loo", "p_loo", "looic", "influence_pareto_k")
  )
  expect_equal(dim(loo1), dim(LLmat))
})


test_that("elpd returns object with correct structure", {
  expect_true(is.loo(elpd1))
  expect_named(
    elpd1,
    c(
      "estimates",
      "pointwise"
    )
  )
  est_names <- dimnames(elpd1$estimates)
  expect_equal(est_names[[1]], c("elpd", "ic"))
  expect_equal(est_names[[2]], c("Estimate", "SE"))
  expect_equal(colnames(elpd1$pointwise), est_names[[1]])
  expect_equal(dim(elpd1), dim(LLmat))
})


test_that("two pareto k values are equal", {
  expect_identical(
    loo1$pointwise[, "influence_pareto_k"],
    loo1$diagnostics$pareto_k
  )
})

test_that("loo.array and loo.matrix give same result", {
  l2 <- suppressWarnings(loo(LLmat, r_eff = r_eff_mat))
  expect_identical(loo1$estimates, l2$estimates)
  expect_identical(loo1$diagnostics, l2$diagnostics)

  # the mcse_elpd_loo columns won't be identical because we use sampling
  expect_identical(loo1$pointwise[, -2], l2$pointwise[, -2])
  expect_equal(loo1$pointwise[, 2], l2$pointwise[, 2], tolerance = 0.005)
})

test_that("loo.array runs with multiple cores", {
  loo_with_arr1 <- loo(LLarr, cores = 1, r_eff = NA)
  loo_with_arr2 <- loo(LLarr, cores = 2, r_eff = NA)
  expect_identical(loo_with_arr1$estimates, loo_with_arr2$estimates)
})

test_that("waic.array and waic.matrix give same result", {
  waic2 <- suppressWarnings(waic(LLmat))
  expect_identical(waic1, waic2)
})

test_that("elpd.array and elpd.matrix give same result", {
  elpd2 <- suppressWarnings(elpd(LLmat))
  expect_identical(elpd1, elpd2)
})

test_that("loo, waic, and elpd error with vector input", {
  expect_error(loo(LLvec), regexp = "no applicable method")
  expect_error(waic(LLvec), regexp = "no applicable method")
  expect_error(elpd(LLvec), regexp = "no applicable method")
})


# testing function methods ------------------------------------------------
source(test_path("data-for-tests/function_method_stuff.R"))

waic_with_fn <- waic(llfun, data = data, draws = draws)
waic_with_mat <- waic(llmat_from_fn)

loo_with_fn <- loo(
  llfun,
  data = data,
  draws = draws,
  r_eff = rep(1, nrow(data))
)
loo_with_mat <- loo(
  llmat_from_fn,
  r_eff = rep(1, ncol(llmat_from_fn)),
  save_psis = TRUE
)

test_that("loo.cores deprecation warning works with function method", {
  options(loo.cores = 1)
  expect_warning(
    loo(
      llfun,
      cores = 2,
      data = data,
      draws = draws,
      r_eff = rep(1, nrow(data))
    ),
    "loo.cores"
  )
  options(loo.cores = NULL)
})

test_that("loo_i results match loo results for ith data point", {
  expect_no_warning(
    loo_i_val <- loo_i(i = 2, llfun = llfun, data = data, draws = draws),
  )
  expect_equal(
    loo_i_val$pointwise[, "elpd_loo"],
    loo_with_fn$pointwise[2, "elpd_loo"]
  )
  expect_equal(
    loo_i_val$pointwise[, "p_loo"],
    loo_with_fn$pointwise[2, "p_loo"]
  )
  expect_equal(
    loo_i_val$diagnostics$pareto_k,
    loo_with_fn$diagnostics$pareto_k[2]
  )
  expect_equal(loo_i_val$diagnostics$n_eff, loo_with_fn$diagnostics$n_eff[2])
})

test_that("function and matrix methods return same result", {
  expect_equal(waic_with_mat, waic_with_fn)
  expect_identical(loo_with_mat$estimates, loo_with_fn$estimates)
  expect_identical(loo_with_mat$diagnostics, loo_with_fn$diagnostics)
  expect_identical(dim(loo_with_mat), dim(loo_with_fn))
})

test_that("loo.function runs with multiple cores", {
  loo_with_fn1 <- loo(
    llfun,
    data = data,
    draws = draws,
    r_eff = rep(1, nrow(data)),
    cores = 1
  )
  loo_with_fn2 <- loo(
    llfun,
    data = data,
    draws = draws,
    r_eff = rep(1, nrow(data)),
    cores = 2
  )
  expect_identical(loo_with_fn2$estimates, loo_with_fn1$estimates)
})

test_that("save_psis option to loo.function makes correct psis object", {
  loo_with_fn2 <- loo.function(
    llfun,
    data = data,
    draws = draws,
    r_eff = rep(1, nrow(data)),
    save_psis = TRUE
  )
  expect_identical(loo_with_fn2$psis_object, loo_with_mat$psis_object)
})

test_that("loo doesn't throw r_eff warnings", {
  expect_no_warning(loo(-LLarr))
  expect_no_warning(loo(-LLmat))
  expect_no_warning(loo(llfun, data = data, draws = draws))
})
