library(loo)

# logColMeansExp_dat in R/sysdata.rda

context("logColMeansExp")
x <- logColMeansExp_dat$x
ans <- logColMeansExp_dat$ans
logColMeansExp_x <- logColMeansExp(x)
lcme_x <- logColMeansExp_x
test_that("logColMeansExp(x) = log(colMeans(exp(x))", {
  expect_equivalent(lcme_x, ans)
})
