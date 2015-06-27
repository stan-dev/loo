library(loo)

# test loo and waic -------------------------------------------------------
context("loo and waic")
test_that("loo and waic return expected results", {

  set.seed(123)
  x <- matrix(rnorm(1000), 20, 50)
  ww <- waic(x)
  wnms <- names(ww)
  ll <- loo(x, cores = 1)
  lnms <- names(ll)
  waic_val <- unlist(ww[-grep("pointwise", wnms)])
  loo_val <- unlist(ll[-grep("pointwise|pareto_k", lnms)])

  loo_ans <- structure(c(-21.8049963976997, 45.7894711322857, 43.6099927953993,
                         1.58654243300237, 1.63674929831605, 3.17308486600474),
                       .Names = c("elpd_loo", "p_loo", "looic", "se_elpd_loo",
                                  "se_p_loo", "se_looic"))
  waic_ans <- structure(c(-25.6339619037423, 49.6184366383284, 51.2679238074847,
                          1.65934536754977, 1.81876567815575, 3.31869073509954),
                        .Names = c("elpd_waic", "p_waic", "waic", "se_elpd_waic",
                                   "se_p_waic", "se_waic"))


  expect_equivalent(loo_val, loo_ans)
  expect_equivalent(waic_val, waic_ans)
  vec <- 1:10
  arr <- array(1:100, dim = c(2,5,10))
  expect_error(loo(vec), 'log_lik should be a matrix')
  expect_error(waic(arr), 'log_lik should be a matrix')
})
