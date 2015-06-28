library(loo)

# test loo and waic -------------------------------------------------------
context("loo and waic")
test_that("loo and waic return expected results", {

  set.seed(123)
  x <- matrix(rnorm(5000), 100, 50)
  ww <- waic(x)
  wnms <- names(ww)
  ll <- loo(x, cores = 1)
  lnms <- names(ll)
  waic_val <- unlist(ww[-grep("pointwise", wnms)])
  loo_val <- unlist(ll[-grep("pointwise|pareto_k", lnms)])

  waic_ans <- structure(c(-25.1291601624985, 49.5591373599165, 50.258320324997,
                          0.742909828603039, 1.05244649202819, 1.48581965720608),
                        .Names = c("elpd_waic", "p_waic", "waic", "se_elpd_waic",
                                   "se_p_waic", "se_waic"))
  loo_ans <- structure(c(-24.2339828660691, 48.6639600634871, 48.4679657321382,
                         0.706772789390649, 0.990389359906155, 1.4135455787813),
                       .Names = c("elpd_loo", "p_loo", "looic", "se_elpd_loo",
                                  "se_p_loo", "se_looic"))


  pareto_k <- ll$pareto_k
  pareto_k_val <- list(mean = mean(pareto_k), range = range(pareto_k))
  pareto_k_ans <- list(mean = 0.273775211941035,
                       range = c(-0.245572962408798, 0.859988648019406))

  expect_equivalent(loo_val, loo_ans)
  expect_equivalent(waic_val, waic_ans)
  vec <- 1:10
  arr <- array(1:100, dim = c(2,5,10))
  expect_error(loo(vec), 'log_lik should be a matrix')
  expect_error(waic(arr), 'log_lik should be a matrix')

  expect_equivalent(pareto_k_val, pareto_k_ans)
})
