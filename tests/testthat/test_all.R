library(loo)


# test helpers ------------------------------------------------------------
context("logColMeansExp")
set.seed(123)
x <- matrix(rnorm(100), 20, 5)
val <- logColMeansExp(x)
ans <- log(colMeans(exp(x)))
test_that("logColMeansExp(x) = log(colMeans(exp(x))", {
  expect_equivalent(val, ans)
})



# test loo_and_waic -------------------------------------------------------
context("loo_and_waic")
set.seed(123)
x <- matrix(rnorm(100), 20, 5)
loo <- loo_and_waic(x, cores = 1)
val <- loo[-grep("pointwise|pareto_k", names(loo))]
ans <- structure(list(elpd_loo = -1.5894575381457, p_loo = 4.02790091225809,
                      elpd_waic = -1.74710271515681, p_waic = 4.1855460892692,
                      looic = 3.1789150762914, waic = 3.49420543031363, se_elpd_loo = 0.559871005470356,
                      se_p_loo = 0.432836117275547, se_elpd_waic = 0.526738752057132,
                      se_p_waic = 0.30556903588444, se_looic = 1.11974201094071,
                      se_waic = 1.05347750411426), .Names = c("elpd_loo", "p_loo",
                                                              "elpd_waic", "p_waic", "looic", "waic", "se_elpd_loo", "se_p_loo",
                                                              "se_elpd_waic", "se_p_waic", "se_looic", "se_waic"))

test_that("loo_and_waic returns correct result", {
  expect_equal(val, ans)
})

context("loo_and_waic_diff")
loo_diff <- unlist(loo_and_waic_diff(loo, loo))
test_that("loo_and_waic_diff returns correct result", {
  expect_true(all(loo_diff == 0))
})

