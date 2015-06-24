library(loo)

# logColMeansExp_dat and loo_and_waic data are in R/sysdata.rda

context("logColMeansExp")
x <- logColMeansExp_dat$x
ans <- logColMeansExp_dat$ans
logColMeansExp_x <- logColMeansExp(x)
lcme_x <- logColMeansExp_x
test_that("logColMeansExp(x) = log(colMeans(exp(x))", {
  expect_equivalent(lcme_x, ans)
})

context("loo_and_waic")
seed <- loo_and_waic_dat$seed
x <- loo_and_waic_dat$x
ans <- loo_and_waic_dat$ans
set.seed(seed)
loo <- loo_and_waic(x, cores = 1)
loo_diff <- unlist(loo_and_waic_diff(loo, loo))
test_that("loo_and_waic returns correct result", {
  expect_identical(loo, ans)
})

context("loo_and_waic_diff")
test_that("loo_and_waic_diff returns correct result", {
  expect_true(all(loo_diff == 0))
})



